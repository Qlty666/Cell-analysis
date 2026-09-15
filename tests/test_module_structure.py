"""Compatibility checks for the refactored pipeline and web module layout."""

from __future__ import annotations

import unittest

import pipeline.differential as differential
import pipeline.integration as integration
import pipeline.key_targets as key_targets
import pipeline.qc as qc
import web_analysis
import web_files
import web_ui
import web_validation


class TestPipelineModuleStructure(unittest.TestCase):
    def test_integration_reexports_domain_helpers(self):
        self.assertIs(
            integration.run_differential_abundance,
            differential.run_differential_abundance,
        )
        self.assertIs(
            integration._bh_adjust,
            differential._bh_adjust,
        )
        self.assertIs(
            integration._chi2_contingency,
            differential._chi2_contingency,
        )
        self.assertIs(
            integration.extract_key_genes,
            key_targets.extract_key_genes,
        )
        self.assertIs(
            integration.DEFAULT_GENE_BLACKLIST,
            key_targets.DEFAULT_GENE_BLACKLIST,
        )
        self.assertIs(
            integration.collect_qc_metrics,
            qc.collect_qc_metrics,
        )
        self.assertIs(
            integration.evaluate_qc_gate,
            qc.evaluate_qc_gate,
        )
        self.assertIs(
            integration.write_qc_metrics,
            qc.write_qc_metrics,
        )


class TestWebModuleStructure(unittest.TestCase):
    def test_web_ui_reexports_domain_helpers(self):
        self.assertIs(
            web_ui.start_analysis_job,
            web_analysis.start_analysis_job,
        )
        self.assertIs(
            web_ui.analysis_results,
            web_analysis.analysis_results,
        )
        self.assertIs(
            web_ui.start_validation_job,
            web_validation.start_validation_job,
        )
        self.assertIs(
            web_ui.validation_report_text,
            web_validation.validation_report_text,
        )
        self.assertIs(
            web_ui._analysis_files,
            web_files.analysis_files,
        )


if __name__ == "__main__":
    unittest.main()
