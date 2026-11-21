
import os
import json
import numpy as np
import pysam

from pathlib import Path
from alphagenome.data import genome
from alphagenome.atlas import atlas


class AlphaGenomeAgent:

    def __init__(self):

        api_key = os.environ.get("ALPHAGENOME_API_KEY")

        if not api_key:
            raise ValueError(
                "ALPHAGENOME_API_KEY environment variable is missing."
            )

        self.client = atlas.create(api_key=api_key)

        # Retrieve available scorer names from official API
        self.scorer_metadata = self.client.scorer_metadata()
        self.scorer_names = list(self.scorer_metadata.keys())

        print("AlphaGenome Atlas initialized.")
        print("Available scorers:", len(self.scorer_names))

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
    # 2. Query AlphaGenome Atlas
    # ----------------------------------

    def query_variant(self, variant):

        ag_variant = genome.Variant(
            chromosome=f"chr{variant['chrom']}",
            position=variant["pos"],
            reference_bases=variant["ref"],
            alternate_bases=variant["alt"]
        )

        return self.client.query_variant(
            ag_variant,
            requested_scorers=self.scorer_names
        )

    # ----------------------------------
    # 3. Summarize scorer results
    # ----------------------------------

    def summarize_scores(self, scores):

        summaries = {}

        for scorer_name, adata in scores.items():

            matrix = adata.X

            if hasattr(matrix, "toarray"):
                matrix = matrix.toarray()

            matrix = np.asarray(matrix, dtype=float)

            if matrix.size == 0:
                continue

            finite_mask = np.isfinite(matrix)

            if not finite_mask.any():
                continue

            valid_values = matrix[finite_mask]

            # Identify the largest absolute scores
            flat_indices = np.flatnonzero(finite_mask)
            top_indices = flat_indices[
                np.argsort(
                    np.abs(valid_values)
                )[::-1][:10]
            ]

            top_scores = []

            n_columns = matrix.shape[1]

            for flat_index in top_indices:

                row_index = flat_index // n_columns
                col_index = flat_index % n_columns

                row_info = (
                    adata.obs.iloc[row_index].to_dict()
                    if row_index < len(adata.obs)
                    else {}
                )

                track_info = (
                    adata.var.iloc[col_index].to_dict()
                    if col_index < len(adata.var)
                    else {}
                )

                top_scores.append({
                    "raw_score": float(matrix[row_index, col_index]),
                    "gene": row_info.get("gene_name"),
                    "gene_id": row_info.get("gene_id"),
                    "track": track_info.get("name"),
                    "biosample": track_info.get("biosample_name"),
                    "biosample_curie": track_info.get("biosample_curie")
                })

            summaries[scorer_name] = {
                "matrix_shape": list(matrix.shape),
                "max_absolute_score": float(
                    np.max(np.abs(valid_values))
                ),
                "mean_score": float(np.mean(valid_values)),
                "top_scores": top_scores
            }

        return summaries

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
                f"Querying AlphaGenome Atlas: "
                f"{index}/{len(variants)} - {variant_id}"
            )

            try:

                scores = self.query_variant(variant)

                if not scores:

                    results.append({
                        **variant,
                        "status": "no_atlas_scores",
                        "scorer_results": {}
                    })

                    continue

                summaries = self.summarize_scores(scores)

                results.append({
                    **variant,
                    "status": (
                        "scores_found"
                        if summaries
                        else "empty_scores"
                    ),
                    "scorer_results": summaries
                })

            except Exception as error:

                print("AlphaGenome error:", error)

                results.append({
                    **variant,
                    "status": "api_error",
                    "error": str(error),
                    "scorer_results": {}
                })

        return results

    # ----------------------------------
    # 5. Run complete analysis
    # ----------------------------------

    def run(
        self,
        vcf_path,
        output_path="data/output/alphagenome_annotation.json"
    ):

        variants = self.read_vcf(vcf_path)

        print("\nAlphaGenome Agent")
        print("------------------------------")
        print("Input variants:", len(variants))

        annotations = self.annotate_variants(variants)

        found = sum(
            item["status"] == "scores_found"
            for item in annotations
        )

        report = {
            "source": "Google DeepMind AlphaGenome Atlas API",
            "assembly": "GRCh38",
            "analysis_type": "Precomputed variant effect scores",
            "total_input_variants": len(variants),
            "variants_with_scores": found,
            "variants_without_scores": (
                len(variants) - found
            ),
            "available_scorers": self.scorer_names,
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
                ensure_ascii=False,
                default=str
            )

        print("\nAlphaGenome Analysis Completed")
        print("------------------------------")
        print("Input variants:", len(variants))
        print("Variants with scores:", found)
        print("Output:", output_path)

        return report


if __name__ == "__main__":

    agent = AlphaGenomeAgent()

    agent.run(
        vcf_path="data/output/sample.normalized.vcf.gz"
    )
