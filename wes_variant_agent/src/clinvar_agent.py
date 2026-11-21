
import json
import time
import requests
import pysam

from pathlib import Path
from urllib.parse import quote


class ClinVarAgent:

    def __init__(self, email=None, api_key=None):

        self.base_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        )

        self.email = email
        self.api_key = api_key

        self.session = requests.Session()

        self.session.headers.update({
            "User-Agent": "WES-Variant-Agent/1.0"
        })

        self.common_params = {}

        if email:
            self.common_params["email"] = email

        if api_key:
            self.common_params["api_key"] = api_key

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
    # 2. Search ClinVar
    # ----------------------------------

    def search_variant(self, variant):

        chromosome = variant["chrom"]
        position = variant["pos"]

        # Search using GRCh38 genomic coordinates
        term = f"{chromosome}:{position}[chrpos38]"

        url = f"{self.base_url}/esearch.fcgi"

        params = {
            **self.common_params,
            "db": "clinvar",
            "term": term,
            "retmode": "json",
            "retmax": 20
        }

        response = self.session.get(
            url,
            params=params,
            timeout=60
        )

        response.raise_for_status()

        data = response.json()

        result = data.get("esearchresult", {})

        return result.get("idlist", [])

    # ----------------------------------
    # 3. Retrieve ClinVar summaries
    # ----------------------------------

    def get_summaries(self, clinvar_ids):

        if not clinvar_ids:
            return {}

        url = f"{self.base_url}/esummary.fcgi"

        params = {
            **self.common_params,
            "db": "clinvar",
            "id": ",".join(clinvar_ids),
            "retmode": "json"
        }

        response = self.session.get(
            url,
            params=params,
            timeout=60
        )

        response.raise_for_status()

        data = response.json()

        return data.get("result", {})

    # ----------------------------------
    # 4. Extract clinical information
    # ----------------------------------

    def extract_clinical_data(self, summary):

        if not isinstance(summary, dict):
            return {}

        clinical = summary.get(
            "clinical_significance", {}
        )

        if not isinstance(clinical, dict):
            clinical = {}

        genes = summary.get("genes", [])

        gene_symbols = []

        if isinstance(genes, list):

            for gene in genes:

                if isinstance(gene, dict):
                    symbol = gene.get("symbol")

                    if symbol:
                        gene_symbols.append(symbol)

        return {
            "variation_id": summary.get("uid"),
            "accession": summary.get("accession"),
            "title": summary.get("title"),
            "variation_type": summary.get("variation_type"),
            "clinical_significance": clinical.get("description"),
            "review_status": clinical.get("review_status"),
            "last_evaluated": clinical.get("last_evaluated"),
            "genes": gene_symbols,
            "trait_set": summary.get("trait_set"),
            "assembly_set": summary.get("assembly_set")
        }

    # ----------------------------------
    # 5. Annotate variants
    # ----------------------------------

    def annotate_variants(self, variants):

        results = []

        for index, variant in enumerate(variants, start=1):

            print(
                f"Searching ClinVar: "
                f"{index}/{len(variants)} "
                f"- {variant['chrom']}:{variant['pos']}"
            )

            try:

                clinvar_ids = self.search_variant(variant)

                time.sleep(0.35)

                summaries = self.get_summaries(clinvar_ids)

                records = []

                for clinvar_id in clinvar_ids:

                    summary = summaries.get(clinvar_id, {})

                    clinical_data = self.extract_clinical_data(
                        summary
                    )

                    records.append(clinical_data)

                results.append({
                    **variant,
                    "search_status": (
                        "records_found"
                        if clinvar_ids
                        else "no_records_found"
                    ),
                    "clinvar_record_count": len(records),
                    "clinvar_records": records
                })

            except requests.RequestException as error:

                print(f"API error: {error}")

                results.append({
                    **variant,
                    "search_status": "api_error",
                    "error": str(error),
                    "clinvar_record_count": 0,
                    "clinvar_records": []
                })

            except (ValueError, KeyError) as error:

                print(f"Data parsing error: {error}")

                results.append({
                    **variant,
                    "search_status": "parsing_error",
                    "error": str(error),
                    "clinvar_record_count": 0,
                    "clinvar_records": []
                })

        return results

    # ----------------------------------
    # 6. Run complete ClinVar analysis
    # ----------------------------------

    def run(
        self,
        vcf_path,
        output_path="data/output/clinvar_annotation.json"
    ):

        variants = self.read_vcf(vcf_path)

        print("\nClinVar Agent")
        print("------------------------------")
        print("Input variants:", len(variants))

        annotations = self.annotate_variants(variants)

        matched = sum(
            1 for item in annotations
            if item["search_status"] == "records_found"
        )

        report = {
            "source": "NCBI ClinVar E-utilities",
            "assembly": "GRCh38",
            "total_input_variants": len(variants),
            "variants_with_records": matched,
            "variants_without_records": (
                len(variants) - matched
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

        print("\nClinVar Analysis Completed")
        print("------------------------------")
        print("Input variants:", len(variants))
        print("Variants with records:", matched)
        print("Output:", output_path)

        return report


if __name__ == "__main__":

    agent = ClinVarAgent()

    report = agent.run(
        vcf_path="data/output/sample.normalized.vcf.gz"
    )
