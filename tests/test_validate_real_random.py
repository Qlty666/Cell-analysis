#!/usr/bin/env python3
"""Unit tests for the random real-data validation script."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import validate_real_random  # noqa: E402


class _Response:
    def __init__(self, payload: bytes):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        return False

    def read(self) -> bytes:
        return self.payload


class TestDownloadPdb(unittest.TestCase):
    def test_retries_transient_download_failure(self):
        calls = []

        def fake_urlopen(request, timeout):
            calls.append((request.full_url, timeout))
            if len(calls) == 1:
                raise TimeoutError("temporary timeout")
            return _Response(b"HEADER\nATOM\n")

        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "1ABC.pdb"
            with patch.object(
                validate_real_random.urllib.request,
                "urlopen",
                side_effect=fake_urlopen,
            ), patch.object(validate_real_random.time, "sleep"):
                validate_real_random.download_pdb("1ABC", dest)
            self.assertEqual(dest.read_bytes(), b"HEADER\nATOM\n")
        self.assertEqual(len(calls), 2)

    def test_removes_partial_file_after_all_attempts_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp) / "1ABC.pdb"
            with patch.object(
                validate_real_random.urllib.request,
                "urlopen",
                side_effect=TimeoutError("timeout"),
            ), patch.object(validate_real_random.time, "sleep"):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "failed after 2 attempts",
                ):
                    validate_real_random.download_pdb(
                        "1ABC",
                        dest,
                        attempts=2,
                    )
            self.assertFalse(dest.exists())


if __name__ == "__main__":
    unittest.main()
