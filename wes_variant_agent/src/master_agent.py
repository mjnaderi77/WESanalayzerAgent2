
import os
import json
import argparse
import inspect
from pathlib import Path
from datetime import datetime, timezone

from dotenv import load_dotenv

from .quality_control import VariantQC
from .vep_agent import VEPAgent
from .clinvar_agent import ClinVarAgent
from .gnomad_agent import GnomADAgent
from .alphamissense_agent import AlphaMissenseAgent
from .alphagenome_agent import AlphaGenomeAgent
from .clinical_impact_engine import ClinicalImpactEngine


load_dotenv()


class WESVariantAgent:

    def __init__(self, tools=None, engine=None):

        self.output_dir = Path("data/output")
        self.output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        self.trace = []

        self.tools = tools or self.default_tools()

        self.engine = engine or self.initialize_engine()

    # -------------------------------------------------
    # Initialize tools
    # -------------------------------------------------

    def default_tools(self):

        tools = {}

        # Quality control
        tools["qc"] = VariantQC(
            min_qual=30,
            min_dp=10
        )

        # VEP
        tools["vep"] = VEPAgent()

        # ClinVar
        try:
            tools["clinvar"] = ClinVarAgent()
        except Exception as e:
            tools["clinvar"] = None
            self.log(
                "INIT",
                "ClinVar disabled",
                str(e)
            )

        # gnomAD
        try:
            tools["gnomad"] = GnomADAgent()
        except Exception as e:
            tools["gnomad"] = None
            self.log(
                "INIT",
                "gnomAD disabled",
                str(e)
            )

        # AlphaMissense
        try:
            tools["alphamissense"] = AlphaMissenseAgent()
        except Exception as e:
            tools["alphamissense"] = None
            self.log(
                "INIT",
                "AlphaMissense disabled",
                str(e)
            )

        # AlphaGenome
        if os.getenv("ALPHAGENOME_API_KEY"):

            try:
                tools["alphagenome"] = AlphaGenomeAgent()

            except Exception as e:
                tools["alphagenome"] = None
                self.log(
                    "INIT",
                    "AlphaGenome disabled",
                    str(e)
                )

        else:
            tools["alphagenome"] = None
            self.log(
                "INIT",
                "AlphaGenome disabled",
                "ALPHAGENOME_API_KEY is missing"
            )

        return tools

    # -------------------------------------------------
    # Initialize Gemini clinical engine
    # -------------------------------------------------

    def initialize_engine(self):

        if not os.getenv("GEMINI_API_KEY"):

            self.log(
                "INIT",
                "ClinicalImpactEngine disabled",
                "GEMINI_API_KEY is not set"
            )

            return None

        try:

            return ClinicalImpactEngine()

        except Exception as e:

            self.log(
                "INIT",
                "ClinicalImpactEngine disabled",
                str(e)
            )

            return None

    # -------------------------------------------------
    # Logging
    # -------------------------------------------------

    def log(
        self,
        step,
        decision,
        reason=""
    ):

        entry = {
            "timestamp": datetime.now(
                timezone.utc
            ).isoformat(),

            "step": step,

            "decision": decision,

            "reason": reason
        }

        self.trace.append(entry)

        print(
            f"[{step}] {decision}"
        )

        if reason:
            print(f"       {reason}")

    # -------------------------------------------------
    # Save JSON
    # -------------------------------------------------

    def save_json(
        self,
        filename,
        data
    ):

        path = self.output_dir / filename

        with open(
            path,
            "w",
            encoding="utf-8"
        ) as file:

            json.dump(
                data,
                file,
                indent=4,
                ensure_ascii=False,
                default=str
            )

        return path

    # -------------------------------------------------
    # Read existing JSON
    # -------------------------------------------------

    @staticmethod
    def read_json(path):

        path = Path(path)

        if not path.exists():
            return None

        try:

            with open(
                path,
                encoding="utf-8"
            ) as file:

                return json.load(file)

        except Exception:
            return None

    # -------------------------------------------------
    # Execute file-based agents
    # -------------------------------------------------

    def call_tool(
        self,
        name,
        input_path,
        output_path,
        method_name="run",
        top_n=None
    ):

        tool = self.tools.get(name)

        if tool is None:

            self.log(
                "ACT",
                f"{name} skipped",
                "Tool is unavailable"
            )

            return None

        method = getattr(
            tool,
            method_name,
            None
        )

        if method is None:

            self.log(
                "ACT",
                f"{name} failed",
                f"Method {method_name} does not exist"
            )

            return None

        try:

            parameters = inspect.signature(
                method
            ).parameters

            kwargs = {}

            if "output_path" in parameters:
                kwargs["output_path"] = str(
                    output_path
                )

            if (
                top_n is not None
                and "top_n" in parameters
            ):
                kwargs["top_n"] = top_n

            self.log(
                "ACT",
                f"Running {name}"
            )

            result = method(
                str(input_path),
                **kwargs
            )

            if result is not None:

                self.save_json(
                    output_path.name,
                    result
                )

            else:

                result = self.read_json(
                    output_path
                )

            self.log(
                "OBSERVE",
                f"{name} completed",
                f"Output: {output_path}"
            )

            return result

        except Exception as e:

            self.log(
                "ERROR",
                f"{name} failed",
                str(e)
            )

            return {
                "status": "error",
                "tool": name,
                "error": str(e)
            }

    # -------------------------------------------------
    # QC execution
    # -------------------------------------------------

    def run_qc(self, vcf_path):

        qc_tool = self.tools["qc"]

        output_path = (
            self.output_dir / "qc_report.json"
        )

        self.log(
            "QC",
            "Running quality control",
            str(vcf_path)
        )

        try:

            result = qc_tool.run_qc(
                str(vcf_path),
                output_path=str(output_path)
            )

            if result is not None:

                self.save_json(
                    output_path.name,
                    result
                )

            else:

                result = self.read_json(
                    output_path
                )

            return result

        except Exception as e:

            self.log(
                "ERROR",
                "QC failed",
                str(e)
            )

            return {
                "status": "error",
                "error": str(e)
            }

    # -------------------------------------------------
    # Extract VEP consequences
    # -------------------------------------------------

    @staticmethod
    def extract_consequences(data):

        consequences = set()

        def walk(obj):

            if isinstance(obj, dict):

                for key, value in obj.items():

                    if (
                        "consequence" in key.lower()
                        and isinstance(value, str)
                    ):
                        consequences.add(value)

                    elif (
                        "consequence" in key.lower()
                        and isinstance(value, list)
                    ):
                        for item in value:
                            if isinstance(item, str):
                                consequences.add(item)

                    walk(value)

            elif isinstance(obj, list):

                for item in obj:
                    walk(item)

        walk(data)

        return consequences

    # -------------------------------------------------
    # PLAN
    # -------------------------------------------------

    def plan(self, vep_results):

        plan = [
            "gnomad",
            "clinvar"
        ]

        consequences = self.extract_consequences(
            vep_results
        )

        missense_terms = {
            "missense_variant"
        }

        alphagenome_terms = {

            "intron_variant",
            "intergenic_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "regulatory_region_variant",
            "TF_binding_site_variant",
            "non_coding_transcript_variant",
            "non_coding_transcript_exon_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "splice_region_variant",
            "splice_donor_variant",
            "splice_acceptor_variant"
        }

        if consequences & missense_terms:

            plan.append("alphamissense")

        if consequences & alphagenome_terms:

            plan.append("alphagenome")

        self.log(
            "PLAN",
            "Tool selection completed",
            f"Consequences: {sorted(consequences)}"
        )

        self.log(
            "PLAN",
            "Selected tools",
            str(plan)
        )

        return plan

    # -------------------------------------------------
    # REFLECT
    # -------------------------------------------------

    def reflect(
        self,
        qc_result,
        tool_results
    ):

        observations = []

        if qc_result:

            total = qc_result.get(
                "total_records"
            )

            passed = qc_result.get(
                "passed"
            )

            failed = qc_result.get(
                "failed"
            )

            observations.append({

                "type": "quality_control",

                "total_variants": total,

                "passed": passed,

                "failed": failed

            })

        for name, result in tool_results.items():

            if result is None:

                observations.append({
                    "type": "tool_unavailable",
                    "tool": name
                })

            elif isinstance(result, dict):

                observations.append({

                    "type": "tool_execution",

                    "tool": name,

                    "status": result.get(
                        "status",
                        "completed"
                    )

                })

        self.log(
            "REFLECT",
            "Evidence review completed",
            f"{len(observations)} observations"
        )

        return observations

    # -------------------------------------------------
    # MASTER RUN
    # -------------------------------------------------

    def run(
        self,
        vcf_path,
        top_n=20
    ):

        vcf_path = Path(vcf_path)

        if not vcf_path.exists():

            raise FileNotFoundError(
                f"VCF file not found: {vcf_path}"
            )

        self.log(
            "START",
            "WES Variant Agent started",
            str(vcf_path)
        )

        # STEP 1: QC

        qc_result = self.run_qc(
            vcf_path
        )

        # STEP 2: VEP

        vep_result = self.call_tool(
            name="vep",
            input_path=vcf_path,
            output_path=(
                self.output_dir / "vep_annotation.json"
            ),
            method_name="annotate_vcf",
            top_n=top_n
        )

        # STEP 3: PLAN

        if vep_result:

            selected_tools = self.plan(
                vep_result
            )

        else:

            selected_tools = [
                "gnomad",
                "clinvar"
            ]

            self.log(
                "PLAN",
                "VEP unavailable",
                "Skipping prediction tools because consequence annotation is unavailable"
            )

        # STEP 4: ACT

        tool_results = {}

        for name in selected_tools:

            output_path = (
                self.output_dir /
                f"{name}_results.json"
            )

            tool_results[name] = self.call_tool(
                name=name,
                input_path=vcf_path,
                output_path=output_path
            )

        # STEP 5: REFLECT

        observations = self.reflect(
            qc_result,
            tool_results
        )

        # STEP 6: Clinical impact

        clinical_result = None

        if self.engine is not None and vep_result is not None:

            try:

                vep_path = (
                    self.output_dir /
                    "vep_annotation.json"
                )

                clinical_output = (
                    self.output_dir /
                    "clinical_impact.json"
                )

                self.log(
                    "SCORE",
                    "Running Clinical Impact Engine"
                )

                clinical_result = self.engine.run(
                    input_path=str(vep_path),
                    output_path=str(clinical_output)
                )

                if clinical_result is not None:

                    self.save_json(
                        clinical_output.name,
                        clinical_result
                    )

                else:

                    clinical_result = self.read_json(
                        clinical_output
                    )

            except Exception as e:

                self.log(
                    "ERROR",
                    "Clinical interpretation failed",
                    str(e)
                )

                clinical_result = {
                    "status": "error",
                    "error": str(e)
                }

        else:

            self.log(
                "SCORE",
                "Clinical interpretation skipped",
                "Gemini engine or VEP result unavailable"
            )

        # STEP 7: FINAL REPORT

        report = {

            "project": "WES Variant Interpretation Agent",

            "input_vcf": str(vcf_path),

            "timestamp": datetime.now(
                timezone.utc
            ).isoformat(),

            "workflow": [
                "QC",
                "VEP",
                "PLAN",
                "ACT",
                "OBSERVE",
                "REFLECT",
                "SCORE"
            ],

            "qc": qc_result,

            "vep": vep_result,

            "selected_tools": selected_tools,

            "tool_results": tool_results,

            "observations": observations,

            "clinical_interpretation": clinical_result

        }

        self.save_json(
            "master_wes_report.json",
            report
        )

        self.save_json(
            "master_trace.json",
            self.trace
        )

        self.log(
            "FINISH",
            "WES analysis completed",
            "master_wes_report.json"
        )

        return report


# -------------------------------------------------
# CLI
# -------------------------------------------------

def main():

    parser = argparse.ArgumentParser(
        description="WES Master Variant Interpretation Agent"
    )

    parser.add_argument(
        "--vcf",
        default="data/input/sample.vcf",
        help="Input VCF file"
    )

    parser.add_argument(
        "--top",
        type=int,
        default=20,
        help="Maximum variants requested by compatible tools"
    )

    args = parser.parse_args()

    agent = WESVariantAgent()

    agent.run(
        args.vcf,
        top_n=args.top
    )


if __name__ == "__main__":

    main()
