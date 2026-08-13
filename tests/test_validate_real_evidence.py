#!/usr/bin/env python3
"""Unit tests for the real evidence validation script."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import validate_real_evidence  # noqa: E402


class TestSkillNames(unittest.TestCase):
    def test_rcsb_skill_directory_is_correct(self):
        self.assertEqual(validate_real_evidence.SKILL_NAMES["rcsb"], "rcsb-pdb-skill")

    def test_all_skill_names_resolve(self):
        skills = Path.home() / ".codex" / "skills"
        for name in validate_real_evidence.SKILL_NAMES.values():
            self.assertTrue((skills / name / "scripts" / "rest_request.py").exists())


if __name__ == "__main__":
    unittest.main()
