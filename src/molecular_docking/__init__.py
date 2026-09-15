"""Standalone molecular docking workflow.

This board is kept independent from the virtual screening pipeline. It
reuses shared AutoDock Vina preparation/docking helpers but owns its own
configuration, CLI, result report and web page.
"""
