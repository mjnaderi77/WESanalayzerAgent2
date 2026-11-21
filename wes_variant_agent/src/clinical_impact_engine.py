
import json
import os
from pathlib import Path

from google import genai
from google.genai import types


PROJECT_DIR = Path(__file__).resolve().parent.parent

DEFAULT_INPUT = PROJECT_DIR / "data/output/master_wes_report.json"
DEFAULT_OUTPUT = PROJECT_DIR / "data/output/clinical_impact_report.json"


OUTPUT_SCHEMA = {
    "type": "OBJECT",
    "properties": {
        "clinical_relevance": {
            "type": "STRING",
            "enum": [
                "potentially_relevant",
                "uncertain",
                "likely_low_relevance",
                "insufficient_evidence"
            ]
        },
        "molecular_interpretation": {"type": "STRING"},
        "disease_association": {"type": "STRING"},
        "supporting_evidence": {
            "type": "ARRAY",
            "items": {"type": "STRING"}
        },
        "conflicting_or_missing_evidence": {
            "type": "ARRAY",
            "items": {"type": "STRING"}
        },
        "confidence": {
            "type": "STRING",
            "enum": ["low", "moderate", "high"]
        },
        "recommended_follow_up": {
            "type": "ARRAY",
            "items": {"type": "STRING"}
        },
        "clinical_summary": {"type": "STRING"},
        "limitations": {
            "type": "ARRAY",
            "items": {"type": "STRING"}
        }
    },
    "required": [
        "clinical_relevance",
        "molecular_interpretation",
        "disease_association",
        "supporting_evidence",
        "conflicting_or_missing_evidence",
        "confidence",
        "recommended_follow_up",
        "clinical_summary",
        "limitations"
    ]
}


SYSTEM_PROMPT = """
You are a clinical genomics evidence interpretation assistant.

Analyze the supplied variant annotations and evidence.

STRICT RULES:
1. Use only the evidence supplied.
2. Never invent gene-disease associations or phenotypes.
3. Distinguish computational predictions from clinical evidence.
4. AlphaMissense and AlphaGenome predictions are not proof of pathogenicity.
5. Interpret gnomAD frequency in context.
6. Do not assume ClinVar records are exact allele matches unless
   chromosome, position, reference and alternate alleles are verified.
7. Do not assign ACMG/AMP classification.
8. Do not diagnose a patient.
9. Explicitly identify missing or conflicting evidence.
10. Recommend expert clinical genetics review when appropriate.
11. Do not infer patient phenotype from variant annotations.
12. If the input is insufficient, say so clearly.

Return a concise evidence-based interpretation.
"""


class ClinicalImpactEngine:

    def __init__(self, model=None):

        api_key = os.getenv("GEMINI_API_KEY")

        if not api_key:
            raise RuntimeError(
                "GEMINI_API_KEY environment variable is not set."
            )

        self.client = genai.Client(api_key=api_key)

        self.model = model or os.getenv(
            "GEMINI_MODEL",
            "gemini-2.5-flash"
        )

    @staticmethod
    def load_report(path):

        with open(path, "r", encoding="utf-8") as file:
            return json.load(file)

    @staticmethod
    def compact_variant(variant):

        payload = {
            "variant": {
                "chrom": variant.get("chrom"),
                "pos": variant.get("pos"),
                "ref": variant.get("ref"),
                "alt": variant.get("alt"),
                "qual": variant.get("qual"),
                "filter": variant.get("filter"),
                "genotype": variant.get("genotype")
            },
            "evidence": variant.get("evidence", {})
        }

        serialized = json.dumps(
            payload,
            ensure_ascii=False,
            default=str
        )

        if len(serialized) > 20000:
            payload["evidence"] = {
                "notice": "Evidence payload was too large and was omitted.",
                "available_sources": list(
                    variant.get("evidence", {}).keys()
                )
            }

        return payload

    def interpret_variant(self, variant):

        payload = self.compact_variant(variant)

        response = self.client.models.generate_content(
            model=self.model,
            contents=[
                SYSTEM_PROMPT,
                "Interpret this variant evidence:",
                json.dumps(payload, ensure_ascii=False, default=str)
            ],
            config=types.GenerateContentConfig(
                response_mime_type="application/json",
                response_schema=OUTPUT_SCHEMA,
                temperature=0.2
            )
        )

        if not response.text:
            raise RuntimeError("Gemini returned an empty response.")

        return json.loads(response.text)

    def run(
        self,
        input_path=DEFAULT_INPUT,
        output_path=DEFAULT_OUTPUT
    ):

        report = self.load_report(input_path)

        variants = report.get("variants", [])

        results = []

        for index, variant in enumerate(variants, start=1):

            variant_id = (
                f"{variant.get('chrom')}:"
                f"{variant.get('pos')}:"
                f"{variant.get('ref')}>"
                f"{variant.get('alt')}"
            )

            print(
                f"[Clinical Impact] "
                f"Processing {index}/{len(variants)}: {variant_id}"
            )

            try:

                interpretation = self.interpret_variant(variant)

                results.append({
                    "variant_id": variant_id,
                    "status": "success",
                    "interpretation": interpretation
                })

            except Exception as error:

                results.append({
                    "variant_id": variant_id,
                    "status": "failed",
                    "error": str(error)
                })

                print(f"[ERROR] {variant_id}: {error}")

        output = {
            "engine": "Clinical Impact Engine",
            "provider": "Google Gemini",
            "model": self.model,
            "input_file": str(input_path),
            "total_variants": len(variants),
            "successful_interpretations": sum(
                item["status"] == "success"
                for item in results
            ),
            "failed_interpretations": sum(
                item["status"] == "failed"
                for item in results
            ),
            "clinical_disclaimer": (
                "AI-generated evidence interpretation only. "
                "Not a diagnosis and not an ACMG/AMP classification. "
                "Requires review by a qualified clinical genetics professional."
            ),
            "results": results
        }

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as file:
            json.dump(
                output,
                file,
                indent=4,
                ensure_ascii=False
            )

        print("\nClinical Impact Engine completed.")
        print(
            f"Successful: {output['successful_interpretations']}"
        )
        print(
            f"Failed: {output['failed_interpretations']}"
        )
        print(f"Output: {output_path}")

        return output


if __name__ == "__main__":

    engine = ClinicalImpactEngine()
    engine.run()
