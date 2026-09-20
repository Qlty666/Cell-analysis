"""Regression coverage for the September correctness review."""
import io
import json
import logging
from pathlib import Path
from subprocess import CompletedProcess
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from common.fingerprints import atomic_json, file_hash, fingerprint
from common.http import HttpError, http_download
from docking.config import load_config
from docking import docking, pipeline, md_simulation
from pipeline.integrated_report import _valid_perturbation
from analysis.model_validation import nested_predictions, splits

ROOT = Path(__file__).resolve().parents[1]
LOG = logging.getLogger(__name__)


def test_bounded_submission_and_invalid_limit():
    from concurrent.futures import ThreadPoolExecutor
    from common.concurrency import bounded_futures
    consumed = []
    def items():
        for value in range(20):
            consumed.append(value)
            yield value
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = bounded_futures(pool, lambda value: value * 2, items(), 3)
        future, item = next(results)
        assert len(consumed) == 3
        values = [future.result()] + [f.result() for f, _ in results]
        assert sorted(values) == list(range(0, 40, 2))
        with pytest.raises(ValueError):
            list(bounded_futures(pool, lambda value: value, [], 0))


def make_cfg(root):
    cfg = load_config(ROOT / "config/docking_config.json", {"workdir": str(root)})
    cfg.receptor_output().parent.mkdir(parents=True, exist_ok=True)
    cfg.receptor_output().write_text("receptor")
    ligand = root / "ligand.pdbqt"
    ligand.write_text("ligand")
    cfg.manifest_path().parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"id": "L1", "pdbqt": str(ligand), "status": "ok"}]).to_csv(cfg.manifest_path(), index=False)
    return cfg, ligand


def fake_dock(cfg, tool, row, folder, log):
    pose = folder / "L1_rep1.pdbqt"
    pose.write_text("pose")
    return {"id": "L1", "replicate": 1, "seed": row["seed"], "affinity": -8, "status": "ok", "pose_file": str(pose)}


def test_resume_changes_and_force(tmp_path):
    cfg, ligand = make_cfg(tmp_path)
    with patch.object(docking, "find_tool", return_value="fake"), patch.object(docking, "_dock_one", side_effect=fake_dock) as run:
        docking.run_docking(cfg, LOG)
        docking.run_docking(cfg, LOG)
        assert run.call_count == 1
        cfg.data["docking"]["seed"] = 999
        docking.run_docking(cfg, LOG)
        assert run.call_count == 2
        ligand.write_text("changed ligand")
        docking.run_docking(cfg, LOG)
        assert run.call_count == 3
        cfg.receptor_output().write_text("changed receptor")
        docking.run_docking(cfg, LOG)
        assert run.call_count == 4
        (cfg.docked_dir() / "L1_rep1.pdbqt").unlink()
        docking.run_docking(cfg, LOG)
        assert run.call_count == 5
        pipeline.run_pipeline(cfg, force=True, stages=[("03", "dock", docking.run_docking)])
        assert run.call_count == 6
        assert len(pd.read_csv(cfg.results_path())) == 1
        assert cfg.get("docking", "resume") is True


def test_failed_jobs_are_retried_without_duplicate_rows(tmp_path):
    cfg, _ = make_cfg(tmp_path)
    with patch.object(docking, "find_tool", return_value="fake"), patch.object(docking, "_dock_one", side_effect=[RuntimeError("interrupted"), fake_dock(cfg, "", {"seed": 42}, tmp_path, LOG)]) as run:
        docking.run_docking(cfg, LOG)
        docking.run_docking(cfg, LOG)
        assert run.call_count == 2
        assert list(pd.read_csv(cfg.results_path()).status) == ["ok"]


class Response(io.BytesIO):
    def __init__(self, data, status=200, **headers):
        super().__init__(data)
        self.status = status
        self.headers = headers


def test_short_download_never_publishes(tmp_path):
    out = tmp_path / "file"
    with patch("urllib.request.urlopen", return_value=Response(b"abc", **{"Content-Length": "100"})):
        with pytest.raises(HttpError, match="incomplete"):
            http_download("https://example.test/file", out, retries=1)
    assert not out.exists()


def partial(out):
    out.with_name(out.name + ".part").write_bytes(b"abc")
    atomic_json(out.with_name(out.name + ".part.json"), {"url": "https://example.test/file", "etag": '"v1"'})


@pytest.mark.parametrize("offset,etag", [(0, '"v1"'), (3, '"v2"')])
def test_invalid_resume_rejected(tmp_path, offset, etag):
    out = tmp_path / "file"
    partial(out)
    with patch("urllib.request.urlopen", return_value=Response(b"def", 206, **{"Content-Length": "3", "Content-Range": f"bytes {offset}-{offset+2}/6", "ETag": etag})):
        with pytest.raises(HttpError):
            http_download("https://example.test/file", out, retries=1)
    assert not out.exists()


def test_valid_resume_and_ignored_range(tmp_path):
    out = tmp_path / "file"
    partial(out)
    with patch("urllib.request.urlopen", return_value=Response(b"def", 206, **{"Content-Length": "3", "Content-Range": "bytes 3-5/6", "ETag": '"v1"'})) as request:
        http_download("https://example.test/file", out, retries=1)
        assert request.call_args.args[0].get_header("If-range") == '"v1"'
    assert out.read_bytes() == b"abcdef"
    partial(out)
    with patch("urllib.request.urlopen", return_value=Response(b"new", **{"Content-Length": "3"})):
        http_download("https://example.test/file", out, retries=1)
    assert out.read_bytes() == b"new"


def test_tpr_only_is_not_completed(tmp_path):
    cfg, _ = make_cfg(tmp_path)
    folder = cfg.md_dir() / "L1"
    folder.mkdir(parents=True)
    (folder / "md.tpr").write_text("tpr")
    with patch.object(md_simulation, "_pick_hits", return_value=pd.DataFrame([{"id": "L1"}])), patch.object(md_simulation, "_prepare_hit_dir", return_value={}), patch.object(md_simulation, "_find_gmx", return_value="fake"), patch.object(md_simulation, "_write_report"):
        with pytest.raises(Exception, match="no MD simulation completed"):
            md_simulation.run_md_simulation(cfg, LOG, mode="auto")
    summary = json.loads((cfg.md_dir() / "md_simulation_summary.json").read_text())
    assert summary["completed"] == 0 and summary["failed"] == 1


@pytest.mark.parametrize("time_ps,exit_code,passed", [(100000, 0, True), (50000, 0, False), (100000, 1, False)])
def test_trajectory_completion(tmp_path, time_ps, exit_code, passed):
    cfg, _ = make_cfg(tmp_path)
    for name in ("md.tpr", "md.xtc", "md.gro", "md.log"):
        (tmp_path / name).write_text("Finished mdrun")
    with patch.object(md_simulation, "run_command", return_value=CompletedProcess([], exit_code, "", f"Last frame 100 time {time_ps}")):
        if passed:
            md_simulation._verify_md_completion("fake", cfg, tmp_path)
        else:
            with pytest.raises(Exception):
                md_simulation._verify_md_completion("fake", cfg, tmp_path)


def test_maxwarn_zero_and_invalid(tmp_path):
    cfg, _ = make_cfg(tmp_path)
    cfg.data["md_simulation"]["maxwarn"] = 0
    commands = []
    with patch.object(md_simulation, "run_command", side_effect=lambda cmd, **kw: (commands.append(cmd) or CompletedProcess(cmd, 0))), patch.object(md_simulation, "_make_whole_trajectory"):
        md_simulation._run_stages("fake", cfg, tmp_path)
    assert all(c[c.index("-maxwarn") + 1] == "0" for c in commands if "grompp" in c)
    for value in (-1, 0.5, float("nan")):
        cfg.data["md_simulation"]["maxwarn"] = value
        with pytest.raises(Exception):
            md_simulation._maxwarn(cfg)


def test_empty_and_valid_perturbation_provenance(tmp_path):
    config_path = ROOT / "config/docking_config.json"
    cfg = load_config(config_path)
    summary_path = tmp_path / "insilico_summary.json"
    for text in ("", "{}", "invalid", '{"status":"failed"}'):
        summary_path.write_text(text)
        assert not _valid_perturbation(summary_path, config_path)
    expression = tmp_path / "expression.csv"
    expression.write_text("gene,sample\nA,1\n")
    output = tmp_path / "targets.csv"
    output.write_text("gene,score\nA,1\n")
    cfg.data["knockout"]["expression_csv"] = str(expression)
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(cfg.data))
    summary = {
        "status": "completed", "schema_version": 2, "ko_gene": "A",
        "cells": 10, "genes_modeled": 5, "top_targets": ["A"],
        "config_signature": fingerprint({"knockout": cfg.data["knockout"], "insilico_knockout": cfg.data.get("insilico_knockout", {})}),
        "inputs": {str(expression.resolve()): file_hash(expression)},
        "output_hashes": {str(output.resolve()): file_hash(output)},
    }
    atomic_json(summary_path, summary)
    assert _valid_perturbation(summary_path, config_path)
    expression.write_text("gene,sample\nA,2\n")
    assert not _valid_perturbation(summary_path, config_path)


def test_integrated_expression_signature_tracks_input_and_environment(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from pipeline.integration import _stage_signature
    raw = tmp_path / "data"
    raw.mkdir()
    source = raw / "counts.csv"
    source.write_text("A,1")
    args = SimpleNamespace()
    context = {"single_cell_root": tmp_path}
    first = _stage_signature("01", args, tmp_path, context)
    source.write_text("A,2")
    second = _stage_signature("01", args, tmp_path, context)
    assert first != second
    monkeypatch.setenv("LIVER_DE_COVARIATES", "batch")
    assert second != _stage_signature("01", args, tmp_path, context)


def test_nested_validation_groups_and_noise():
    from sklearn.pipeline import Pipeline
    from sklearn.impute import SimpleImputer
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(42)
    X = pd.DataFrame(rng.normal(size=(160, 300)))
    y = np.tile([0, 0, 1, 1], 40)
    groups = np.repeat(np.arange(80), 2)
    for train, test in splits(X, y, 4, 42, groups):
        assert not set(groups[train]) & set(groups[test])
    model = Pipeline([("impute", SimpleImputer()), ("select", SelectKBest(f_classif, k=5)), ("clf", LogisticRegression())])
    p, records = nested_predictions({"logistic": model}, X, y, 1, 4, 42, groups)
    assert 0.30 < roc_auc_score(y, p) < 0.70
    assert all(len(r["features"]) == 5 and not set(r["train_samples"]) & set(r["test_samples"]) for r in records)
