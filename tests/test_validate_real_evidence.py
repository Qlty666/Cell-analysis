#!/usr/bin/env python3
"""Unit tests for the real evidence validation script."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import validate_real_evidence  # noqa: E402


class TestSkillNames(unittest.TestCase):
    def test_rcsb_skill_directory_is_correct(self):
        self.assertEqual(validate_real_evidence.SKILL_NAMES["rcsb"], "rcsb-pdb-skill")

    def test_all_skill_names_resolve(self):
        skills = Path.home() / ".codex" / "skills"
        scripts = [
            skills / name / "scripts" / "rest_request.py"
            for name in validate_real_evidence.SKILL_NAMES.values()
        ]
        if not any(script.exists() for script in scripts):
            self.skipTest("Codex skill scripts are not installed")
        for script in scripts:
            self.assertTrue(script.exists())

    def test_call_skill_returns_failure_when_script_missing(self):
        with patch.object(validate_real_evidence.Path, "exists", return_value=False):
            result = validate_real_evidence.call_skill("rcsb", {})
        self.assertFalse(result["ok"])
        self.assertIn("missing skill script", result["error"]["message"])


class TestThresholds(unittest.TestCase):
    def test_defaults_require_most_targets(self):
        self.assertEqual(validate_real_evidence.DEFAULT_MIN_OK_TARGETS, 10)
        self.assertEqual(validate_real_evidence.DEFAULT_MIN_LIGANDS, 10)
        self.assertEqual(
            validate_real_evidence.DEFAULT_MIN_OK_TARGETS,
            len(validate_real_evidence.TARGETS),
        )


if __name__ == "__main__":
    unittest.main()
