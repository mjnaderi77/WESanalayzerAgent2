
import pysam
import json
from pathlib import Path


class VariantQC:

    def __init__(self, min_qual=30, min_dp=10):
        self.min_qual = min_qual
        self.min_dp = min_dp

    def evaluate_variant(self, record):

        qual = record.qual
        filters = list(record.filter.keys())

        checks = {
            "qual_pass": (
                qual is not None
                and qual >= self.min_qual
            ),
            "filter_pass": filters == ["PASS"]
        }

        return {
            "chrom": record.chrom,
            "pos": record.pos,
            "ref": record.ref,
            "alt": record.alts,
            "qual": qual,
            "filters": filters,
            "checks": checks,
            "variant_qc_pass": all(checks.values())
        }

    def evaluate_sample(self, record, sample_name):

        sample = record.samples[sample_name]

        dp = sample.get("DP")
        gt = sample.get("GT")
        ad = sample.get("AD")

        allele_balance = None

        if (
            ad is not None
            and all(a is not None for a in ad)
            and sum(ad) > 0
        ):
            total_ad = sum(ad)

            allele_balance = [
                round(depth / total_ad, 4)
                for depth in ad
            ]

        checks = {
            "dp_pass": (
                dp is not None
                and dp >= self.min_dp
            ),
            "genotype_available": (
                gt is not None
                and all(a is not None for a in gt)
            )
        }

        return {
            "sample": sample_name,
            "dp": dp,
            "genotype": gt,
            "ad": ad,
            "allele_balance": allele_balance,
            "ad_available": (
                ad is not None
                and len(ad) >= 2
                and all(a is not None for a in ad)
            ),
            "checks": checks,
            "sample_qc_pass": all(checks.values())
        }


def run_qc(
    vcf_path,
    output_path="data/output/qc_report.json"
):

    qc = VariantQC()

    results = []

    with pysam.VariantFile(vcf_path) as vcf:

        samples = list(vcf.header.samples)

        for record in vcf:

            variant_result = qc.evaluate_variant(record)

            sample_results = []

            for sample_name in samples:

                sample_result = qc.evaluate_sample(
                    record,
                    sample_name
                )

                sample_results.append(sample_result)

            final_result = {
                **variant_result,
                "samples": sample_results,
                "qc_pass": (
                    variant_result["variant_qc_pass"]
                    and all(
                        s["sample_qc_pass"]
                        for s in sample_results
                    )
                )
            }

            results.append(final_result)

    output = {
        "total_variants": len(results),
        "passed": sum(
            r["qc_pass"] for r in results
        ),
        "failed": sum(
            not r["qc_pass"] for r in results
        ),
        "thresholds": {
            "min_qual": qc.min_qual,
            "min_dp": qc.min_dp
        },
        "variants": results
    }

    Path(output_path).parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(output_path, "w") as f:
        json.dump(
            output,
            f,
            indent=4,
            default=str
        )

    return output


if __name__ == "__main__":

    report = run_qc(
        "data/output/sample.normalized.vcf.gz"
    )

    print(json.dumps(report, indent=4, default=str))
