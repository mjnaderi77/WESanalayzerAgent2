
from src.vcf_reader import VCFReader

vcf = VCFReader("data/input/sample.vcf.gz")

print("METADATA")
print(vcf.get_metadata())

print("\nVARIANTS")

for variant in vcf.read_variants():
    print(variant)

vcf.close()
