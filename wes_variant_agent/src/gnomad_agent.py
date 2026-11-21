
import json
import time
import requests
import pysam

from pathlib import Path


class GnomADAgent:

    def __init__(self, dataset="gnomad_r4"):

        self.api_url = "https://gnomad.broadinstitute.org/api"
        self.dataset = dataset

        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
            "User-Agent": "WES-Variant-Agent/1.0"
        })

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
    # 2. Query gnomAD GraphQL API
    # ----------------------------------

    def query_variant(self, variant):

        variant_id = (
            f"{variant['chrom']}-"
            f"{variant['pos']}-"
            f"{variant['ref']}-"
            f"{variant['alt']}"
        )

        query = """
        query GetVariant($variantId: String!, $dataset: DatasetId!) {
            variant(variantId: $variantId, dataset: $dataset) {
                variant_id
                rsids
                chrom
                pos
                ref
                alt
                exome {
                    ac
                    an
                    af
                }
                genome {
                    ac
                    an
                    af
                }
            }
        }
        """

        payload = {
            "query": query,
            "variables": {
                "variantId": variant_id,
                "dataset": self.dataset
            }
        }

        response = self.session.post(
            self.api_url,
            json=payload,
            timeout=60
        )

        response.raise_for_status()

        data = response.json()

        if data.get("errors"):
            raise RuntimeError(
                f"GraphQL error: {data['errors']}"
            )

        return data.get("data", {}).get("variant")

    # ----------------------------------
    # 3. Extract population information
    # ----------------------------------

    def extract_frequency(self, data):

        if not data:
            return None

        return {
            "allele_count": data.get("ac"),
            "allele_number": data.get("an"),
            "allele_frequency": data.get("af")
        }

    # ----------------------------------
    # 4. Annotate variants
    # ----------------------------------

    def annotate_variants(self, variants):

        results = []

        for index, variant in enumerate(variants, start=1):

            variant_id = (
                f"{variant['chrom']}-"
                f"{variant['pos']}-"
                f"{variant['ref']}-"
                f"{variant['alt']}"
            )

            print(
                f"Querying gnomAD: {index}/{len(variants)} "
                f"- {variant_id}"
            )

            try:

                data = self.query_variant(variant)

                time.sleep(0.3)

                if data is None:

                    results.append({
                        **variant,
                        "gnomad_status": "not_found",
                        "gnomad_records": None
                    })

                    continue

                # Verify exact REF/ALT match
                exact_match = (
                    str(data.get("chrom")) == variant["chrom"]
                    and int(data.get("pos")) == variant["pos"]
                    and data.get("ref") == variant["ref"]
                    and data.get("alt") == variant["alt"]
                )

                if not exact_match:

                    results.append({
                        **variant,
                        "gnomad_status": "allele_mismatch",
                        "gnomad_records": None
                    })

                    continue

                results.append({
                    **variant,
                    "gnomad_status": "found",
                    "variant_id": data.get("variant_id"),
                    "rsids": data.get("rsids"),
                    "exome": self.extract_frequency(
                        data.get("exome")
                    ),
                    "genome": self.extract_frequency(
                        data.get("genome")
                    )
                })

            except requests.RequestException as error:

                print("API connection error:", error)

                results.append({
                    **variant,
                    "gnomad_status": "api_error",
                    "error": str(error)
                })

            except (RuntimeError, ValueError, KeyError) as error:

                print("Annotation error:", error)

                results.append({
                    **variant,
                    "gnomad_status": "annotation_error",
                    "error": str(error)
                })

        return results

    # ----------------------------------
    # 5. Run complete analysis
    # ----------------------------------

    def run(
        self,
        vcf_path,
        output_path="data/output/gnomad_annotation.json"
    ):

        variants = self.read_vcf(vcf_path)

        print("\nGnomAD Agent")
        print("------------------------------")
        print("Dataset:", self.dataset)
        print("Input variants:", len(variants))

        annotations = self.annotate_variants(variants)

        found = sum(
            item["gnomad_status"] == "found"
            for item in annotations
        )

        report = {
            "source": "gnomAD GraphQL API",
            "dataset": self.dataset,
            "assembly": "GRCh38",
            "total_input_variants": len(variants),
            "variants_found": found,
            "variants_not_found": sum(
                item["gnomad_status"] == "not_found"
                for item in annotations
            ),
            "annotations": annotations
        }

        Path(output_path).parent.mkdir(
            parents=True,
            exist_ok=True
        )

        with open(output_path, "w", encoding="utf-8") as file:

            json.dump(
                report,
                file,
                indent=4,
                ensure_ascii=False
            )

        print("\ngnomAD Analysis Completed")
        print("------------------------------")
        print("Input variants:", len(variants))
        print("Variants found:", found)
        print("Output:", output_path)

        return report


if __name__ == "__main__":

    agent = GnomADAgent(dataset="gnomad_r4")

    report = agent.run(
        vcf_path="data/output/sample.normalized.vcf.gz"
    )
