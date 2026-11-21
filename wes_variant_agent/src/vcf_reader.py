
from pathlib import Path
import pysam


class VCFReader:

    def __init__(self, vcf_path: str):

        self.vcf_path = Path(vcf_path)

        if not self.vcf_path.exists():
            raise FileNotFoundError(
                f"VCF file not found: {self.vcf_path}"
            )

        self.vcf = pysam.VariantFile(str(self.vcf_path))
        self.samples = list(self.vcf.header.samples)

    def get_metadata(self):

        return {
            "file": self.vcf_path.name,
            "samples": self.samples,
            "sample_count": len(self.samples),
            "contigs": list(self.vcf.header.contigs),
        }

    def read_variants(self):

        variants = []

        for record in self.vcf:

            for alt in record.alts or ():

                variants.append({
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "ref": record.ref,
                    "alt": alt,
                    "qual": record.qual,
                    "filter": list(record.filter.keys()),
                    "info": dict(record.info),
                })

        return variants

    def close(self):
        self.vcf.close()

