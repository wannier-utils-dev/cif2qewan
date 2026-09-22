"""The output boundary: write rendered inputs into a directory."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Mapping, Union

PathLike = Union[str, Path]


def write_rendered(rendered: Mapping[str, str], directory: PathLike) -> List[Path]:
    """Write every ``relative path -> text`` pair below ``directory``.

    Sub-directories (``band/``, ``check_wannier/``) are created as needed.
    Existing files are overwritten. Returns the written paths in order.
    """
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    written = []
    for relative, text in rendered.items():
        target = root / relative
        if not _is_below(target, root):
            raise ValueError(f"refusing to write outside {root}: {relative}")
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
        written.append(target)
    return written


def _is_below(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def listed(paths: Iterable[Path], root: PathLike) -> List[str]:
    """The paths relative to ``root`` as strings, for messages."""
    root = Path(root)
    return [str(p.relative_to(root)) for p in paths]
