"""Conflict-aware final decisions for validated targets."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from docking.utils import write_json


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        import json

        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except Exception:
        return {}


def _downgrade(decision: str) -> str:
    return {
        "GO": "CONDITIONAL_GO",
        "CONDITIONAL_GO": "REVIEW",
        "REVIEW": "REVIEW",
        "NO_GO": "NO_GO",
    }.get(str(decision), "REVIEW")


def _action(decision: str, flags: list[str]) -> str:
    if decision == "NO_GO":
        return "NO_GO"
    if flags:
        return "REVIEW_CONFLICT"
    if decision == "GO":
        return "PROCEED_VALIDATION"
    if decision == "CONDITIONAL_GO":
        return "CONDITIONAL_PROCEED"
    return "REVIEW"


def build_decision_report(out_dir: Path) -> tuple[pd.DataFrame, dict]:
    """Combine target validation with global scientific-readiness conflicts."""
    validation = pd.read_csv(out_dir / "target_validation_scores.csv")
    evidence = _read_json(out_dir / "evidence_hub" / "evidence_hub_summary.json")
    external = _read_json(out_dir / "external_validation_summary.json")
    omics = _read_json(out_dir / "omics_qc_summary.json")
    structural = _read_json(out_dir / "structural_quality_summary.json")
    global_flags: list[str] = []
    if str(evidence.get("status", "")) != "completed":
        global_flags.append("evidence_hub_incomplete")
    if str(external.get("status", "")) != "completed":
        global_flags.append("external_validation_missing")
    elif float(external.get("auroc") or 0.0) < 0.70:
        global_flags.append("external_validation_weak")
    if str(omics.get("status", "")) != "completed" or not omics.get(
        "gate_passed",
        False,
    ):
        global_flags.append("omics_qc_not_passed")
    structural_gates = structural.get("gates") or {}
    for gate in ("positive_control", "replicate_consensus", "md_rmsd_stability"):
        if not structural_gates.get(gate, False):
            global_flags.append(f"structural_{gate}_missing")

    rows: list[dict] = []
    for _, row in validation.iterrows():
        flags = list(global_flags)
        safety = float(row.get("safety_risk", 0.0) or 0.0)
        if safety >= 0.5:
            flags.append("high_safety_risk")
        decision = str(row.get("decision") or "REVIEW")
        if safety >= 0.5:
            decision = _downgrade(decision)
        rows.append(
            {
                "gene": str(row.get("gene") or ""),
                "base_decision": row.get("decision"),
                "final_decision": decision,
                "action": _action(decision, flags),
                "adjusted_score": row.get("adjusted_score"),
                "safety_risk": safety,
                "conflict_flags": ";".join(flags),
            }
        )
    frame = pd.DataFrame(rows).sort_values(
        ["final_decision", "adjusted_score", "gene"],
        ascending=[True, False, True],
    )
    frame.to_csv(out_dir / "target_decision_report.csv", index=False)
    summary = {
        "targets": int(len(frame)),
        "global_conflicts": global_flags,
        "actions": {
            str(key): int(value)
            for key, value in frame["action"].value_counts().items()
        },
        "final_decisions": {
            str(key): int(value)
            for key, value in frame["final_decision"].value_counts().items()
        },
        "top_targets": frame.head(20).to_dict(orient="records"),
    }
    write_json(out_dir / "target_decision_report.json", summary)
    return frame, summary
