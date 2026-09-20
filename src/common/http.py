"""Shared HTTP helpers with retries, backoff, size limits and clear errors."""

from __future__ import annotations

import logging
import os
import re
import urllib.error
import time
import urllib.request
from pathlib import Path
from .fingerprints import atomic_json, read_state

logger = logging.getLogger(__name__)

# Replace with a real maintainer address (or set the caller's user_agent) when
# deploying so remote services can contact the operator about traffic.
CONTACT_EMAIL = (
    os.environ.get("LIVER_CONTACT_EMAIL", "").strip()
    or "not-configured@example.invalid"
)

DEFAULT_USER_AGENT = (
    "liver-cancer-pipeline/1.0 "
    f"(+https://github.com/Qlty666/Cell-analysis; mailto:{CONTACT_EMAIL})"
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
    memory. A ``.part`` file with a matching resource validator is resumed
    with Range/If-Range; servers ignoring the range cause a clean restart. Failures
    raise :class:`HttpError` naming the URL.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    partial = out.with_name(out.name + ".part")
    metadata = partial.with_name(partial.name + ".json")
    attempts = max(1, int(retries))
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        saved = read_state(metadata)
        validator = saved.get("etag") or saved.get("last_modified")
        existing = partial.stat().st_size if partial.exists() and saved.get("url") == url and validator else 0
        try:
            request = _build_request(url, user_agent=user_agent)
            if existing:
                request.add_header("Range", f"bytes={existing}-")
                request.add_header("If-Range", validator)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", 200)
                resuming = bool(existing) and status == 206
                total = None
                if status == 206:
                    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", ""))
                    if not match:
                        raise HttpError("missing or invalid Content-Range")
                    start, end, total = map(int, match.groups())
                    if start != existing or end < start or end >= total:
                        raise HttpError("incorrect Content-Range for requested offset")
                    if existing and ((saved.get("etag") and response.headers.get("ETag") != saved["etag"]) or (not saved.get("etag") and response.headers.get("Last-Modified") != saved.get("last_modified"))):
                        raise HttpError("resource changed during resumed download")
                if existing and not resuming:
                    existing = 0
                declared = _content_length(response)
                expected = None if declared is None else declared + existing
                if total is not None:
                    if declared is not None and declared != end - start + 1:
                        raise HttpError("Content-Length disagrees with Content-Range")
                    expected = total
                _check_declared_size(url, expected, max_bytes)
                etag = response.headers.get("ETag")
                if etag and etag.startswith("W/"):
                    etag = None
                atomic_json(metadata, {"url": url, "etag": etag, "last_modified": response.headers.get("Last-Modified")})
                written = existing
                with partial.open("ab" if resuming else "wb") as fh:
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
            if not partial.exists() or partial.stat().st_size == 0:
                raise RuntimeError("downloaded file is empty")
            if expected is not None and written != expected:
                raise OSError(f"incomplete download: expected {expected} bytes, received {written}")
            os.replace(partial, out)
            metadata.unlink(missing_ok=True)
            return out
        except HttpError:
            partial.unlink(missing_ok=True)
            metadata.unlink(missing_ok=True)
            raise
        except urllib.error.HTTPError as exc:
            # An unsatisfied range does not prove that a local file is complete.
            if exc.code == 416:
                partial.unlink(missing_ok=True)
                metadata.unlink(missing_ok=True)
            elif exc.code not in (408, 429, 500, 502, 503, 504):
                raise HttpError(f"GET {url} failed: HTTP {exc.code}") from exc
            last_error = exc
            if attempt < attempts:
                time.sleep(backoff ** attempt)
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
