"""Conflict-aware final decisions for validated targets."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from docking.utils import write_json


ACTION_SEVERITY = {
    "NO_GO": 0,
    "REVIEW_CONFLICT": 1,
    "REVIEW": 2,
    "CONDITIONAL_PROCEED": 3,
    "PROCEED_VALIDATION": 4,
}


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
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


def _target_action(decision: str, flags: list[str]) -> str:
    if decision == "NO_GO":
        return "NO_GO"
    if decision == "REVIEW":
        return "REVIEW_CONFLICT" if flags else "REVIEW"
    if decision == "CONDITIONAL_GO":
        return "REVIEW_CONFLICT" if flags else "CONDITIONAL_PROCEED"
    if decision == "GO":
        return "REVIEW_CONFLICT" if flags else "PROCEED_VALIDATION"
    return "REVIEW"


def _platform_action(flags: list[str]) -> str:
    return "PLATFORM_REVIEW" if flags else "PLATFORM_READY"


def _composite_action(
    target_action: str,
    platform_flags: list[str],
) -> str:
    if ACTION_SEVERITY.get(target_action, 2) <= ACTION_SEVERITY["REVIEW"]:
        return target_action
    if platform_flags:
        return "REVIEW_CONFLICT"
    return target_action


def _composite_decision(decision: str, platform_flags: list[str]) -> str:
    if not platform_flags:
        return decision
    if decision in {"GO", "CONDITIONAL_GO"}:
        return "REVIEW"
    return decision


def build_decision_report(out_dir: Path) -> tuple[pd.DataFrame, dict]:
    """Separate target-level evidence from platform-level readiness conflicts."""
    validation = pd.read_csv(out_dir / "target_validation_scores.csv")
    evidence = _read_json(out_dir / "evidence_hub" / "evidence_hub_summary.json")
    external = _read_json(out_dir / "external_validation_summary.json")
    omics = _read_json(out_dir / "omics_qc_summary.json")
    structural = _read_json(out_dir / "structural_quality_summary.json")

    platform_flags: list[str] = []
    if str(evidence.get("status", "")) != "completed":
        platform_flags.append("evidence_hub_incomplete")
    if str(external.get("status", "")) != "completed":
        platform_flags.append("external_validation_missing")
    elif float(external.get("auroc") or 0.0) < 0.70:
        platform_flags.append("external_validation_weak")
    if str(omics.get("status", "")) != "completed" or not omics.get(
        "gate_passed",
        False,
    ):
        platform_flags.append("omics_qc_not_passed")
    structural_gates = structural.get("gates") or {}
    for gate in (
        "positive_control",
        "replicate_consensus",
        "md_rmsd_stability",
    ):
        if not structural_gates.get(gate, False):
            platform_flags.append(f"structural_{gate}_missing")

    rows: list[dict] = []
    for _, row in validation.iterrows():
        target_flags: list[str] = []
        safety = float(row.get("safety_risk", 0.0) or 0.0)
        if safety >= 0.5:
            target_flags.append("high_safety_risk")
        target_decision = str(row.get("decision") or "REVIEW")
        if safety >= 0.5:
            target_decision = _downgrade(target_decision)
        target_action = _target_action(target_decision, target_flags)
        platform_action = _platform_action(platform_flags)
        composite_action = _composite_action(
            target_action,
            platform_flags,
        )
        rows.append(
            {
                "gene": str(row.get("gene") or ""),
                "base_decision": row.get("decision"),
                "target_decision": target_decision,
                "final_decision": target_decision,
                "composite_decision": _composite_decision(
                    target_decision,
                    platform_flags,
                ),
                "target_action": target_action,
                "platform_action": platform_action,
                "composite_action": composite_action,
                "action": composite_action,
                "adjusted_score": row.get("adjusted_score"),
                "safety_risk": safety,
                "target_flags": ";".join(target_flags),
                "platform_conflicts": ";".join(platform_flags),
                "conflict_flags": ";".join(
                    [*target_flags, *platform_flags]
                ),
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        frame = pd.DataFrame(
            columns=[
                "gene",
                "base_decision",
                "target_decision",
                "final_decision",
                "composite_decision",
                "target_action",
                "platform_action",
                "composite_action",
                "action",
                "adjusted_score",
                "safety_risk",
                "target_flags",
                "platform_conflicts",
                "conflict_flags",
            ]
        )
    else:
        frame["_severity"] = frame["composite_action"].map(
            ACTION_SEVERITY
        ).fillna(2)
        frame = frame.sort_values(
            ["_severity", "adjusted_score", "gene"],
            ascending=[True, False, True],
        ).drop(columns=["_severity"])
    frame.to_csv(out_dir / "target_decision_report.csv", index=False)
    summary = {
        "targets": int(len(frame)),
        "platform_conflicts": platform_flags,
        "platform_action": _platform_action(platform_flags),
        "actions": {
            str(key): int(value)
            for key, value in frame["composite_action"].value_counts().items()
        }
        if not frame.empty
        else {},
        "target_actions": {
            str(key): int(value)
            for key, value in frame["target_action"].value_counts().items()
        }
        if not frame.empty
        else {},
        "final_decisions": {
            str(key): int(value)
            for key, value in frame["final_decision"].value_counts().items()
        }
        if not frame.empty
        else {},
        "top_targets": frame.head(20).to_dict(orient="records")
        if not frame.empty
        else [],
    }
    write_json(out_dir / "target_decision_report.json", summary)
    return frame, summary
