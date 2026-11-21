
import json
import requests
import pysam

from pathlib import Path


class AlphaMissenseAgent:

    def __init__(self, species="homo_sapiens"):

        self.base_url = "https://rest.ensembl.org"
        self.endpoint = f"/vep/{species}/region"

        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

        self.params = {
            "AlphaMissense": 1,
            "canonical": 1,
            "hgvs": 1,
            "protein": 1,
            "symbol": 1,
            "mane": 1
        }

        self.session = requests.Session()

    # ----------------------------------
    # 1. Read variants from VCF
    # ----------------------------------

    def read_vcf(self, vcf_path):

        variants = []

        with pysam.VariantFile(vcf_path) as vcf:

            for record in vcf:

                for alt in record.alts or []:

                    variants.append({
                        "chrom": str(record.chrom).replace("chr", ""),
                        "pos": int(record.pos),
                        "ref": record.ref,
                        "alt": alt
                    })

        return variants

    # ----------------------------------
    # 2. Convert variants to VEP format
    # ----------------------------------

    def format_variants(self, variants):

        formatted = []

        for variant in variants:

            formatted.append(
                f"{variant['chrom']} "
                f"{variant['pos']} "
                f". "
                f"{variant['ref']} "
                f"{variant['alt']} "
                f". . ."
            )

        return formatted

    # ----------------------------------
    # 3. Query Ensembl VEP + AlphaMissense
    # ----------------------------------

    def annotate_batch(self, variants):

        url = self.base_url + self.endpoint

        response = self.session.post(
            url,
            headers=self.headers,
            params=self.params,
            json={"variants": variants},
            timeout=120
        )

        response.raise_for_status()

        return response.json()

    # ----------------------------------
    # 4. Extract AlphaMissense predictions
    # ----------------------------------

    def extract_predictions(self, annotations):

        results = []

        for annotation in annotations:

            input_fields = annotation.get("input", "").split()

            if len(input_fields) < 5:
                continue

            chrom = input_fields[0]
            pos = int(input_fields[1])
            ref = input_fields[3]
            alt = input_fields[4]

            transcript_results = []

            for transcript in annotation.get(
                "transcript_consequences", []
            ):

                consequences = transcript.get(
                    "consequence_terms", []
                )

                if "missense_variant" not in consequences:
                    continue

                transcript_results.append({

                    "gene_id": transcript.get("gene_id"),
                    "gene_symbol": transcript.get("gene_symbol"),
                    "transcript_id": transcript.get("transcript_id"),
                    "protein_id": transcript.get("protein_id"),
                    "hgvsp": transcript.get("hgvsp"),
                    "amino_acids": transcript.get("amino_acids"),
                    "protein_position": transcript.get("protein_start"),

                    "am_pathogenicity": transcript.get(
                        "am_pathogenicity"
                    ),

                    "am_class": transcript.get("am_class")

                })

            if transcript_results:

                has_prediction = any(
                    item["am_pathogenicity"] is not None
                    for item in transcript_results
                )

                status = (
                    "prediction_found"
                    if has_prediction
                    else "score_unavailable"
                )

            else:

                status = "no_missense_consequence"

            results.append({

                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,

                "most_severe_consequence": annotation.get(
                    "most_severe_consequence"
                ),

                "status": status,

                "transcript_predictions": transcript_results

            })

        return results

    # ----------------------------------
    # 5. Run complete analysis
    # ----------------------------------

    def run(
        self,
        vcf_path,
        output_path="data/output/alphamissense_annotation.json"
    ):

        variants = self.read_vcf(vcf_path)

        formatted_variants = self.format_variants(variants)

        print("\nAlphaMissense Agent")
        print("------------------------------")
        print("Input variants:", len(variants))

        annotations = []

        batch_size = 100

        for i in range(0, len(formatted_variants), batch_size):

            batch = formatted_variants[i:i + batch_size]

            print(
                f"Annotating batch "
                f"{i // batch_size + 1}..."
            )

            raw_results = self.annotate_batch(batch)

            annotations.extend(
                self.extract_predictions(raw_results)
            )

        found = sum(
            item["status"] == "prediction_found"
            for item in annotations
        )

        missense = sum(
            item["status"] in [
                "prediction_found",
                "score_unavailable"
            ]
            for item in annotations
        )

        report = {

            "source": "Ensembl VEP with AlphaMissense",
            "species": "homo_sapiens",
            "assembly": "GRCh38",

            "total_input_variants": len(variants),
            "missense_variants": missense,
            "predictions_found": found,

            "annotations": annotations

        }

        Path(output_path).parent.mkdir(
            parents=True,
            exist_ok=True
        )

        with open(
            output_path,
            "w",
            encoding="utf-8"
        ) as file:

            json.dump(
                report,
                file,
                indent=4,
                ensure_ascii=False
            )

        print("\nAlphaMissense Analysis Completed")
        print("------------------------------")
        print("Input variants:", len(variants))
        print("Missense variants:", missense)
        print("Predictions found:", found)
        print("Output:", output_path)

        return report


if __name__ == "__main__":

    agent = AlphaMissenseAgent()

    agent.run(
        vcf_path="data/output/sample.normalized.vcf.gz"
    )
