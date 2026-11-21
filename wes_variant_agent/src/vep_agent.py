
import pysam
import requests
import json
from pathlib import Path


class VEPAgent:

    def __init__(self, species="homo_sapiens"):
        self.base_url = "https://rest.ensembl.org"
        self.endpoint = f"/vep/{species}/region"

        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

        self.params = {
            "canonical": 1,
            "hgvs": 1,
            "protein": 1,
            "symbol": 1,
            "mane": 1
        }

    def read_vcf(self, vcf_path):
        variants = []

        with pysam.VariantFile(vcf_path) as vcf:
            for record in vcf:
                for alt in record.alts or []:
                    variant = (
                        f"{record.chrom} "
                        f"{record.pos} "
                        f". "
                        f"{record.ref} "
                        f"{alt} "
                        f". . ."
                    )
                    variants.append(variant)

        return variants

    def annotate_batch(self, variants):
        url = self.base_url + self.endpoint

        response = requests.post(
            url,
            headers=self.headers,
            params=self.params,
            json={"variants": variants},
            timeout=120
        )

        response.raise_for_status()
        return response.json()

    def flatten_annotations(self, annotations):

        flattened = []

        for annotation in annotations:

            input_fields = annotation.get("input", "").split()

            if len(input_fields) >= 5:
                chrom = input_fields[0]
                pos = int(input_fields[1])
                ref = input_fields[3]
                alt = input_fields[4]
            else:
                chrom = annotation.get("seq_region_name")
                pos = annotation.get("start")
                ref = None
                alt = None

            variant_summary = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "most_severe_consequence":
                    annotation.get("most_severe_consequence"),
                "assembly": annotation.get("assembly_name")
            }

            for transcript in annotation.get(
                "transcript_consequences", []
            ):

                row = {
                    **variant_summary,
                    "gene_id": transcript.get("gene_id"),
                    "gene_symbol": transcript.get("gene_symbol"),
                    "transcript_id": transcript.get("transcript_id"),
                    "consequence": transcript.get("consequence_terms"),
                    "impact": transcript.get("impact"),
                    "biotype": transcript.get("biotype"),
                    "hgvsc": transcript.get("hgvsc"),
                    "hgvsp": transcript.get("hgvsp"),
                    "protein_start": transcript.get("protein_start"),
                    "amino_acids": transcript.get("amino_acids"),
                    "canonical": transcript.get("canonical"),
                    "mane": transcript.get("mane_select")
                }

                flattened.append(row)

        return flattened

    def annotate_vcf(
        self,
        vcf_path,
        output_path="data/output/vep_annotation.json"
    ):

        variants = self.read_vcf(vcf_path)
        annotations = []

        batch_size = 100

        for i in range(0, len(variants), batch_size):

            batch = variants[i:i + batch_size]

            print(f"Annotating batch {i // batch_size + 1}...")

            results = self.annotate_batch(batch)
            annotations.extend(results)

        flattened = self.flatten_annotations(annotations)

        output = {
            "source": "Ensembl VEP REST API",
            "species": "homo_sapiens",
            "assembly": "GRCh38",
            "total_input_variants": len(variants),
            "total_annotated": len(annotations),
            "total_transcript_annotations": len(flattened),
            "annotations": annotations,
            "transcript_annotations": flattened
        }

        Path(output_path).parent.mkdir(
            parents=True,
            exist_ok=True
        )

        with open(output_path, "w") as f:
            json.dump(output, f, indent=4)

        return output


if __name__ == "__main__":

    agent = VEPAgent()

    report = agent.annotate_vcf(
        "data/output/sample.normalized.vcf.gz"
    )

    print("\nVEP Annotation Summary")
    print("Input variants:", report["total_input_variants"])
    print("Annotated variants:", report["total_annotated"])
    print(
        "Transcript annotations:",
        report["total_transcript_annotations"]
    )
    print("Output: data/output/vep_annotation.json")
