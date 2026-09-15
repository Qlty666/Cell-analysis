#!/usr/bin/env python3
"""Small HTTP GET helper with retries and explicit failure reporting."""

from __future__ import annotations

import time
import urllib.request

USER_AGENT = "Mozilla/5.0 (compatible; liver-cancer-pipeline)"


def http_get(
    url: str,
    timeout: int = 90,
    retries: int = 3,
    backoff: float = 1.0,
    log=None,
    headers: dict | None = None,
) -> str:
    """Fetch ``url`` as text, retrying transient network failures.

    Raises the last exception when every attempt fails so callers can decide
    how to surface it (the evidence module returns ``{"ok": False, ...}``).
    """
    attempts = max(1, int(retries))
    last_exc: Exception | None = None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(
            url,
            headers={"User-Agent": USER_AGENT, **(headers or {})},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as resp:
                return resp.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - retried below
            last_exc = exc
            if log is not None:
                log.warning(
                    "http_get attempt %s/%s failed for %s: %s",
                    attempt,
                    attempts,
                    url,
                    exc,
                )
            if attempt < attempts:
                time.sleep(max(0.0, float(backoff)) * attempt)
    assert last_exc is not None
    raise last_exc
