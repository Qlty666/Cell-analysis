"""Shared exceptions for the integrated pipeline."""


class IntegrationError(RuntimeError):
    """Raised when the integrated pipeline cannot continue."""


class PauseRequested(Exception):
    """Raised when the user requested a clean pipeline pause."""
