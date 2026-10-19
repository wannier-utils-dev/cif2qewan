#!/usr/bin/env python
"""Generate a cif2qewan pseudopotential table for a PseudoDojo set.

The bundled PseudoDojo tables follow one rule set, which this script makes
reproducible:

- ``pp_file_name`` is ``<element>_sr`` (the ``.UPF`` extension is added by
  cif2qewan, and ``--so`` turns ``_sr`` into ``_fr``), so the files of the
  scalar-relativistic set have to be renamed ``X.upf`` -> ``X_sr.UPF`` and
  those of the fully relativistic set ``X.upf`` -> ``X_fr.UPF``;
- ``ecutwfc`` is twice the "high" hint of the ``.djrepo`` file (hartree ->
  rydberg) and ``ecutrho`` is four times ``ecutwfc``; elements whose
  ``.djrepo`` has no hints take them from a fallback set (``--fallback-djrepo``);
- ``zval`` is ``z_valence`` of the UPF file;
- ``nexclude`` and ``orbitals`` are the choices of a reference table (one of
  the bundled tables) for the same element. When the valence configuration
  of the new pseudopotential differs from the reference's (compared through
  ``--reference-upf``), the semicore shells that were added or removed are
  counted into ``nexclude`` and the element is reported for a manual check.

Example (the bundled v0.4 PBEsol table)::

    python tools/pseudodojo_table.py \
        --upf nc-sr-04_pbesol_standard_upf --djrepo nc-sr-04_pbesol_standard_djrepo \
        --reference cif2qewan/nc-sr-05_pbe_standard_upf.csv \
        --reference-upf nc-sr-05_pbe_standard_upf \
        --output cif2qewan/nc-sr-04_pbesol_standard_upf.csv

``--check`` compares the generated rows with ``--output`` instead of
writing it (used to confirm that the bundled tables follow the rules).
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from cif2qewan.structure.model import ELEMENTS

COLUMNS = ("atom", "pp_file_name", "nexclude", "orbitals", "ecutwfc", "ecutrho", "zval")
ANGULAR = {0: "s", 1: "p", 2: "d", 3: "f"}
HARTREE_TO_RYDBERG = 2.0
ECUTRHO_FACTOR = 4.0

Shell = Tuple[int, str]  # (n, l letter)


def valence_configuration(upf: Path) -> Tuple[Shell, ...]:
    """The valence shells listed in the ONCVPSP input echoed in PP_INFO."""
    text = upf.read_text(errors="replace")
    header = re.search(
        r"# atsym\s+z\s+nc\s+nv.*?\n(\S+)\s+([\d.]+)\s+(\d+)\s+(\d+)", text
    )
    if header is None:
        raise ValueError(f"{upf}: no ONCVPSP input block (atsym/nc/nv) found")
    n_core, n_valence = int(header.group(3)), int(header.group(4))
    rows = re.findall(r"^\s*(\d+)\s+(\d+)\s+([\d.]+)\s*$", text[header.end() :], re.M)[
        : n_core + n_valence
    ]
    if len(rows) != n_core + n_valence:
        raise ValueError(f"{upf}: expected {n_core + n_valence} configuration rows")
    return tuple((int(n), ANGULAR[int(l)]) for n, l, _ in rows[n_core:])


def valence_charge(upf: Path) -> Optional[float]:
    """``z_valence`` of the UPF header, or None if it is not found."""
    match = re.search(
        r'z_valence\s*=\s*"?\s*([\d.Ee+-]+)', upf.read_text(errors="replace")
    )
    return float(match.group(1)) if match else None


def bands(shell: Shell) -> int:
    return {"s": 1, "p": 3, "d": 5, "f": 7}[shell[1]]


def high_hint(djrepo: Path) -> Optional[float]:
    """The "high" cutoff hint in hartree, or None if the file has no hints."""
    data = json.loads(djrepo.read_text())
    hints = data.get("hints")
    if not hints:
        return None
    return float(hints["high"]["ecut"])


def cutoffs(element: str, djrepo_dirs: List[Path]) -> Tuple[Optional[float], str]:
    """(ecutwfc in Ry, where it came from) from the first set that has hints."""
    for directory in djrepo_dirs:
        path = directory / f"{element}.djrepo"
        if path.is_file():
            hint = high_hint(path)
            if hint is not None:
                return HARTREE_TO_RYDBERG * hint, directory.name
    return None, ""


def read_table(path: Path) -> Dict[str, Dict[str, str]]:
    with path.open(newline="") as handle:
        return {row["atom"]: row for row in csv.DictReader(handle)}


def format_number(value: float) -> str:
    return f"{value:.1f}"


def build_rows(args: argparse.Namespace, notes: List[str]) -> List[Dict[str, str]]:
    reference = read_table(args.reference) if args.reference else {}
    djrepo_dirs = [args.djrepo] + list(args.fallback_djrepo)
    rows = []
    for element in ELEMENTS:
        row = {column: "" for column in COLUMNS}
        row["atom"] = element
        row["pp_file_name"] = f"{element}{args.suffix}"
        rows.append(row)
        upf = args.upf / f"{element}.upf"
        if not upf.is_file():
            continue

        zval = valence_charge(upf)
        if zval is None:
            notes.append(
                f"WARNING {element}: no z_valence in {upf.name}; zval left empty"
            )
        else:
            row["zval"] = f"{zval:g}"

        ecutwfc, source = cutoffs(element, djrepo_dirs)
        if ecutwfc is None:
            notes.append(f"WARNING {element}: no cutoff hints in any djrepo set")
        else:
            row["ecutwfc"] = format_number(ecutwfc)
            row["ecutrho"] = format_number(ECUTRHO_FACTOR * ecutwfc)
            if source != args.djrepo.name:
                notes.append(f"NOTE {element}: cutoff hints taken from {source}")

        ref = reference.get(element)
        if ref is None or ref["nexclude"] == "":
            notes.append(
                f"WARNING {element}: no reference entry; nexclude/orbitals left empty"
            )
            continue
        row["orbitals"] = ref["orbitals"]
        nexclude = int(float(ref["nexclude"]))
        if args.reference_upf is not None:
            ref_upf = args.reference_upf / f"{element}.upf"
            if ref_upf.is_file():
                new = valence_configuration(upf)
                old = valence_configuration(ref_upf)
                if set(new) != set(old):
                    added = sorted(set(new) - set(old))
                    removed = sorted(set(old) - set(new))
                    nexclude += sum(bands(s) for s in added)
                    nexclude -= sum(bands(s) for s in removed)
                    notes.append(
                        f"CHECK {element}: valence {' '.join(f'{n}{l}' for n, l in old)}"
                        f" -> {' '.join(f'{n}{l}' for n, l in new)}; nexclude "
                        f"{ref['nexclude']} -> {nexclude} (added/removed shells "
                        "counted as semicore)"
                    )
            else:
                notes.append(
                    f"NOTE {element}: not in the reference UPF set; "
                    "nexclude copied unchanged"
                )
        row["nexclude"] = str(nexclude)
    return rows


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--upf", type=Path, required=True, help="UPF directory of the set"
    )
    parser.add_argument(
        "--djrepo", type=Path, required=True, help="djrepo directory of the set"
    )
    parser.add_argument(
        "--fallback-djrepo",
        type=Path,
        action="append",
        default=[],
        help="djrepo directory used for elements without hints (repeatable)",
    )
    parser.add_argument(
        "--reference", type=Path, help="table providing nexclude and orbitals"
    )
    parser.add_argument(
        "--reference-upf", type=Path, help="UPF directory the reference table describes"
    )
    parser.add_argument(
        "--suffix", default="_sr", help="file-name suffix (default _sr)"
    )
    parser.add_argument("--output", type=Path, required=True, help="CSV to write")
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare with --output instead of writing it; exit 1 on differences",
    )
    args = parser.parse_args(argv)

    notes: List[str] = []
    rows = build_rows(args, notes)
    for note in notes:
        print(note, file=sys.stderr)

    if args.check:
        existing = read_table(args.output)
        differences = []
        for row in rows:
            old = existing.get(row["atom"])
            if old is None:
                differences.append(f"{row['atom']}: missing in {args.output}")
            elif {k: old[k] for k in COLUMNS} != row:
                differences.append(f"{row['atom']}: {dict(old)} != {row}")
        for line in differences:
            print(line)
        print(f"{len(differences)} differences for {args.output}")
        return 1 if differences else 0

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"wrote {args.output} ({sum(1 for r in rows if r['ecutwfc'])} elements with data)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
