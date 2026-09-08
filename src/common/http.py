"""Shared HTTP helpers with retries, backoff, size limits and clear errors."""

from __future__ import annotations

import logging
import time
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

# Replace with a real maintainer address (or set the caller's user_agent) when
# deploying so remote services can contact the operator about traffic.
CONTACT_PLACEHOLDER = "your-email@example.com"

DEFAULT_USER_AGENT = (
    "liver-cancer-pipeline/1.0 "
    f"(+https://github.com/Qlty666/Cell-analysis; mailto:{CONTACT_PLACEHOLDER})"
)


class HttpError(RuntimeError):
    """Raised when an HTTP request fails after all retries."""


def _build_request(
    url: str,
    *,
    user_agent: str | None,
) -> urllib.request.Request:
    return urllib.request.Request(
        url,
        headers={"User-Agent": user_agent or DEFAULT_USER_AGENT},
    )


def _content_length(response) -> int | None:
    declared = response.headers.get("Content-Length")
    if declared is not None and str(declared).isdigit():
        return int(declared)
    return None


def _check_declared_size(url: str, declared: int | None, max_bytes: int | None) -> None:
    if max_bytes is not None and declared is not None and declared > max_bytes:
        raise HttpError(
            f"response for {url} declares {declared} bytes, "
            f"above the {max_bytes} byte limit"
        )


def http_get(
    url: str,
    *,
    timeout: float = 60,
    retries: int = 3,
    backoff: float = 1.5,
    max_bytes: int | None = None,
    user_agent: str | None = None,
) -> bytes:
    """GET ``url`` and return the raw response body.

    Transient failures are retried ``retries`` times with exponential backoff
    (``backoff ** attempt`` seconds). ``max_bytes`` caps the accepted response
    size so a misbehaving server cannot exhaust memory. Every failure ends as
    an :class:`HttpError` that names the URL and the last underlying error.
    """
    attempts = max(1, int(retries))
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            request = _build_request(url, user_agent=user_agent)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                _check_declared_size(url, _content_length(response), max_bytes)
                if max_bytes is None:
                    return response.read()
                body = response.read(max_bytes + 1)
                if len(body) > max_bytes:
                    raise HttpError(
                        f"response for {url} exceeds the {max_bytes} byte limit"
                    )
                return body
        except HttpError:
            raise
        except Exception as exc:  # noqa: BLE001 - normalized into HttpError below
            last_error = exc
            if attempt < attempts:
                delay = backoff ** attempt
                logger.warning(
                    "GET %s failed (attempt %d/%d): %s; retrying in %.1fs",
                    url,
                    attempt,
                    attempts,
                    exc,
                    delay,
                )
                time.sleep(delay)
    raise HttpError(
        f"GET {url} failed after {attempts} attempt(s): {last_error}"
    ) from last_error


def http_download(
    url: str,
    out: Path,
    *,
    timeout: float = 600,
    retries: int = 3,
    backoff: float = 3.0,
    max_bytes: int | None = None,
    user_agent: str | None = None,
    log=None,
) -> Path:
    """Stream ``url`` into ``out``, resuming a partial file when supported.

    The response is written in chunks so large archives never have to fit in
    memory. A partial ``out`` from an earlier attempt is resumed with a Range
    request; servers that ignore the range cause a clean restart. Failures
    raise :class:`HttpError` naming the URL.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    attempts = max(1, int(retries))
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        existing = out.stat().st_size if out.exists() else 0
        try:
            request = _build_request(url, user_agent=user_agent)
            if existing:
                request.add_header("Range", f"bytes={existing}-")
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", 200)
                resuming = bool(existing) and status == 206
                if existing and not resuming:
                    existing = 0
                declared = _content_length(response)
                expected = None if declared is None else declared + existing
                _check_declared_size(url, expected, max_bytes)
                written = existing
                with out.open("ab" if resuming else "wb") as fh:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        written += len(chunk)
                        if max_bytes is not None and written > max_bytes:
                            raise HttpError(
                                f"download of {url} exceeds the "
                                f"{max_bytes} byte limit"
                            )
                        fh.write(chunk)
            if not out.exists() or out.stat().st_size == 0:
                raise RuntimeError("downloaded file is empty")
            return out
        except HttpError:
            raise
        except Exception as exc:  # noqa: BLE001 - normalized into HttpError below
            last_error = exc
            if attempt < attempts:
                delay = backoff ** attempt
                if log is not None:
                    log(
                        f"download attempt {attempt}/{attempts} failed ({exc}); "
                        f"retrying in {delay:.0f}s"
                    )
                time.sleep(delay)
    raise HttpError(
        f"GET {url} failed after {attempts} attempt(s): {last_error}"
    ) from last_error
