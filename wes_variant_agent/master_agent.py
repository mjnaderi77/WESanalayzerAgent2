"""Master agent: plans, acts, observes, reflects, prioritizes, reports.

Loop per variant:  PLAN (which tools?) -> ACT (call tools, retry/fallback)
-> OBSERVE (store evidence) -> REFLECT (flags, skip decisions) -> SCORE.
Every decision is written to `trace`, so the run can be audited.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

from .alphagenome_agent import AlphaGenomeAgent
from .alphamissense_agent import AlphaMissenseAgent
from .clinical_impact_engine import ClinicalImpactEngine
from .clinvar_agent import ClinVarAgent
from .gnomad_agent import GnomADAgent
from .quality_control import VariantQC
from .vep_agent import VEPAgent

OUTPUT_DIR = Path(__file__).resolve().parent.parent / "data" / "output"

COMMON_AF = 0.05      # above this the variant is treated as benign-like
LOW_FREQ_AF = 0.01
IMPACT_POINTS = {"HIGH": 3.0, "MODERATE": 2.0, "LOW": 0.5, "MODIFIER": 0.0}


def default_tools() -> dict:
    tools = {
        "qc": VariantQC(),
        "vep": VEPAgent(),
        "clinvar": ClinVarAgent(),
        "gnomad": GnomADAgent(),
        "alphamissense": AlphaMissenseAgent(),
    }

    # AlphaGenome is optional – if API key is missing, register a safe dummy
    try:
        tools["alphagenome"] = AlphaGenomeAgent()
    except (ValueError, KeyError, OSError) as e:
        class _DummyAlphaGenome:
            def available(self) -> bool:
                return False

            def run(self, *args, **kwargs):
                return {"ok": False, "error": "ALPHAGENOME_API_KEY is missing"}

        tools["alphagenome"] = _DummyAlphaGenome()
        print(f"[warning] AlphaGenome disabled: {e}")

    return tools


class WESVariantAgent:
    def __init__(self, tools: dict = None, engine: ClinicalImpactEngine = None):
        self.tools = tools or default_tools()
        self.engine = engine or ClinicalImpactEngine()
        self.trace = []

    # ---------- helpers ----------
    def log(self, step, decision, reason=""):
        self.trace.append({
            "time": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "step": step,
            "decision": decision,
            "reason": reason,
        })
        print(f"[{step}] {decision}" + (f" - {reason}" if reason else ""))

    def call_tool(self, name, variant) -> dict:
        """Run one tool; never let a tool failure crash the whole agent."""
        try:
            return self.tools[name].run(variant)
        except Exception as exc:
            self.log("act", f"{name} failed for {variant['key']}", str(exc)[:150])
            return {"ok": False, "error": str(exc)[:200]}

    # ---------- planning ----------
    def plan(self, variant, evidence) -> list:
        """Decide which tools are needed for THIS variant, given what we know."""
        steps = ["gnomad", "clinvar"]
        vep = evidence.get("vep", {})
        if "missense_variant" in vep.get("consequence_terms", []):
            if self.tools["alphamissense"].available():
                steps.append("alphamissense")
            else:
                self.log("plan", "skip alphamissense", "file not configured")
        elif vep.get("impact") == "MODIFIER" and self.tools["alphagenome"].available():
            steps.append("alphagenome")
        return steps

    # ---------- reflection ----------
    def reflect(self, rec) -> None:
        e, flags = rec["evidence"], rec["flags"]
        af = e.get("gnomad", {}).get("af")
        cls = (e.get("clinvar", {}).get("classification") or "").lower()
        if af is not None and af >= COMMON_AF:
            flags.append("common_variant")
        if "pathogenic" in cls and "benign" not in cls and af is not None and af > LOW_FREQ_AF:
            flags.append("conflict_clinvar_vs_gnomad")
        if "conflicting" in cls:
            flags.append("clinvar_conflicting_interpretations")
        if not e.get("vep", {}).get("gene_symbol"):
            flags.append("no_gene_annotation")
        failed = [n for n, v in e.items() if isinstance(v, dict) and v.get("ok") is False]
        if failed:
            flags.append("incomplete_evidence:" + "+".join(failed))
        if rec.get("allele_balance") and rec.get("genotype") == [0, 1]:
            ab = rec["allele_balance"][1]
            if ab < 0.2 or ab > 0.8:
                flags.append("skewed_allele_balance")

    # ---------- scoring ----------
    @staticmethod
    def score(rec) -> float:
        e, s = rec["evidence"], 0.0
        s += IMPACT_POINTS.get(e.get("vep", {}).get("impact"), 0.0)
        cls = (e.get("clinvar", {}).get("classification") or "").lower()
        if "conflicting" in cls:
            s += 2
        elif "pathogenic" in cls and "benign" not in cls:
            s += 5
        elif "uncertain" in cls:
            s += 1
        elif "benign" in cls:
            s -= 3
        am = e.get("alphamissense", {}).get("class") or ""
        if am == "likely_pathogenic":
            s += 2
        elif am == "likely_benign":
            s -= 1
        af = e.get("gnomad", {}).get("af")
        if af is None:
            s += 1.5 if e.get("gnomad", {}).get("ok") else 0.0
        elif af < 1e-4:
            s += 1.0
        elif af > COMMON_AF:
            s -= 5
        elif af > LOW_FREQ_AF:
            s -= 3
        return round(s, 2)

    # ---------- main run ----------
    def run(self, vcf_path, top_n=10, output_dir=OUTPUT_DIR) -> dict:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        self.log("qc", f"running QC on {vcf_path}")
        qc = self.tools["qc"].run(vcf_path, output_path=out / "qc_report.json")
        passed = [v for v in qc["variants"] if v["qc_pass"]]
        self.log("qc", f"{len(passed)}/{qc['total_records']} variants passed")
        if not passed:
            return self._finish(out, qc, [], [], "No variant passed QC")

        self.log("annotate", "calling VEP for all passed variants")
        try:
            vep = self.tools["vep"].run(passed)["by_variant"]
        except Exception as exp:
            vep = {}
            self.log("annotate", "VEP failed - continuing with limited evidence", str(exp)[:150])

        records = []
        for v in passed:
            rec = {**v, "evidence": {"vep": vep.get(v["key"], {})}, "flags": []}
            if v["key"] not in vep:
                rec["evidence"]["vep"] = {"ok": False}
            for tool_name in self.plan(v, rec["evidence"]):
                rec["evidence"][tool_name] = self.call_tool(tool_name, v)
                if tool_name == "gnomad":
                    af = rec["evidence"]["gnomad"].get("af")
                    if af is not None and af >= COMMON_AF:
                        self.log("plan", f"stop early for {v['key']}", f"AF={af:.3f} is common")
                        break
            self.reflect(rec)
            rec["priority_score"] = self.score(rec)
            records.append(rec)

        records.sort(key=lambda r: r["priority_score"], reverse=True)
        self.log("interpret", f"clinical impact for top {min(top_n, len(records))} variants")
        for rec in records[:top_n]:
            if "common_variant" in rec["flags"]:
                continue  # no need to spend LLM calls on common variants
            rec["clinical_impact"] = self.engine.interpret(rec)
        return self._finish(out, qc, records, records[:top_n], None)

    def _finish(self, out, qc, records, top, note) -> dict:
        report = {
            "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "note": note,
            "qc_summary": {k: qc[k] for k in ("total_records", "passed", "failed", "thresholds")},
            "ranked_variants": records,
            "agent_trace": self.trace,
        }
        (out / "clinical_impact_report.json").write_text(json.dumps(report, indent=2))
        (out / "clinical_impact_report.md").write_text(self.to_markdown(report, top))
        self.log("report", f"saved to {out}")
        return report

    @staticmethod
    def to_markdown(report, top) -> str:
        lines = [
            "# WES Variant Interpretation Report (research use)",
            "",
            f"Generated: {report['generated']}",
            "",
            f"QC: {report['qc_summary']['passed']}/{report['qc_summary']['total_records']} passed",
            "",
        ]
        if report.get("note"):
            lines += [f"**Note:** {report['note']}", ""]
        for i, r in enumerate(top, 1):
            vep = r["evidence"].get("vep", {})
            lines += [
                f"## {i}. {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} "
                f"({vep.get('gene_symbol') or 'no gene'})",
                f"- Priority score: {r['priority_score']}",
                f"- Consequence: {', '.join(vep.get('consequence_terms') or []) or 'n/a'} "
                f"({vep.get('impact')})",
                f"- Flags: {', '.join(r['flags']) or 'none'}",
            ]
            ci = r.get("clinical_impact")
            if ci:
                lines += [
                    f"- Interpretation ({ci['source']}): {ci['clinical_significance']}",
                    f"- Confidence: {ci['confidence']}",
                    "- Supporting: " + "; ".join(ci["evidence_summary"]["supporting"]),
                    "- Against: " + "; ".join(ci["evidence_summary"]["against"]),
                    "- Next steps: " + "; ".join(ci["further_investigation"]),
                ]
            lines.append("")
        lines.append("> Not a clinical diagnosis. Review by a clinical geneticist is required.")
        return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="WES variant interpretation agent")
    parser.add_argument("--vcf", default="data/input/sample.vcf")
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()
    WESVariantAgent().run(args.vcf, top_n=args.top)


if __name__ == "__main__":
    main()
