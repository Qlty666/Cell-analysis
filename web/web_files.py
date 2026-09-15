"""Result-file discovery helpers shared by web result pages."""

from __future__ import annotations

from pathlib import Path

from web_data import (
    RESULT_FILE_SUFFIXES,
    RESULT_IMAGE_SUFFIXES,
)


def is_result_file(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in RESULT_FILE_SUFFIXES


def list_result_files(root: Path) -> list[str]:
    if not root.exists() or not root.is_dir():
        return []
    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if is_result_file(path)
    )


def list_result_images(root: Path) -> list[str]:
    if not root.exists() or not root.is_dir():
        return []
    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file()
        and path.suffix.lower() in RESULT_IMAGE_SUFFIXES
    )


def analysis_files(root: Path, images_only: bool = False) -> list[str]:
    if not root.exists() or not root.is_dir():
        return []
    suffixes = (
        RESULT_IMAGE_SUFFIXES
        if images_only
        else RESULT_IMAGE_SUFFIXES
        | {".csv", ".html", ".json", ".md", ".xlsx"}
    )
    return sorted(
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes
    )
