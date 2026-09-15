"""Shared parsing and path helpers for the web console."""

from __future__ import annotations

import re
from pathlib import Path


def first_value(data: dict, key: str, default: str = "") -> str:
    values = data.get(key)
    if not values:
        return default
    return str(values[0])


def float3(data: dict, prefix: str):
    values = [
        first_value(data, f"{prefix}_x"),
        first_value(data, f"{prefix}_y"),
        first_value(data, f"{prefix}_z"),
    ]
    if any(str(value).strip() == "" for value in values):
        return None
    return [float(value) for value in values]


def integer_field(data: dict, key: str):
    value = first_value(data, key, "")
    if str(value).strip() == "":
        return None
    return int(float(value))


def float_field(data: dict, key: str):
    value = first_value(data, key, "")
    if str(value).strip() == "":
        return None
    return float(value)


def integer_list_field(data: dict, key: str) -> list[int] | None:
    value = first_value(data, key, "")
    if str(value).strip() == "":
        return None
    tokens = [part for part in re.split(r"[\s,;]+", str(value)) if part]
    return [int(float(part)) for part in tokens]


def raw_count_flag(value: str) -> bool | None:
    text = str(value or "").strip().lower()
    if text in ("", "auto"):
        return None
    return text in ("1", "true", "yes", "on")


def truthy(value: object) -> bool:
    return str(value or "").strip().lower() in ("1", "true", "on", "yes")


def cli_path(value: str, base: Path) -> str:
    """Return an absolute path relative to the project root."""
    value = str(value or "").strip()
    if not value:
        return ""
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return str(path.resolve())


def analysis_flags(data: dict, names: tuple[str, ...]) -> list[str]:
    return [
        "--" + name.replace("_", "-")
        for name in names
        if truthy(first_value(data, name, ""))
    ]


def require_existing_file(value: str, label: str, base: Path) -> str:
    resolved = cli_path(value, base)
    if not Path(resolved).is_file():
        raise ValueError(f"{label}不存在：{value}")
    return resolved
