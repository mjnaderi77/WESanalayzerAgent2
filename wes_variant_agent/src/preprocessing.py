
from pathlib import Path
import subprocess


def normalize_vcf(input_vcf: str, output_vcf: str):

    input_path = Path(input_vcf)
    output_path = Path(output_vcf)

    if not input_path.exists():
        raise FileNotFoundError(input_path)

    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True)

    if output_path.exists():
        raise FileExistsError(
            f"Output already exists: {output_path}"
        )

    command = [
        "bcftools",
        "norm",
        "-m", "-any",
        "-Oz",
        "-o", str(output_path),
        str(input_path)
    ]

    subprocess.run(command, check=True)

    subprocess.run(
        ["bcftools", "index", "-f", "-t", str(output_path)],
        check=True
    )

    return {
        "input": str(input_path),
        "output": str(output_path),
        "index": str(output_path) + ".tbi",
        "status": "success"
    }


if __name__ == "__main__":

    result = normalize_vcf(
        "data/input/sample.vcf.gz",
        "data/output/sample.normalized.vcf.gz"
    )

    print(result)
