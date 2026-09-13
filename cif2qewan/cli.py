"""Command-line interface of cif2qewan.

::

    cif2qewan <cif_file> <toml_file> [--so] [--mag]
              [--reader {cif2cell,pymatgen}] [--cif2cell-output FILE]
              [--output-dir DIR]

The CLI parses the arguments, runs the reader -> builder -> renderer ->
writer pipeline and reports errors; every decision about the calculation
lives in :mod:`cif2qewan.workflow.builder`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from cif2qewan import __version__
from cif2qewan.config import Config
from cif2qewan.exceptions import Cif2qewanError
from cif2qewan.structure.readers import Cif2cellReader, PymatgenReader
from cif2qewan.workflow.builder import StructureHints, build_plan
from cif2qewan.workflow.output import write_rendered
from cif2qewan.workflow.render import render_plan

CIF2CELL_OUTPUT_NAME = Cif2cellReader.output_name  # cif_scf.in


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cif2qewan",
        description="Generate Quantum ESPRESSO and Wannier90 input files from a CIF file.",
    )
    parser.add_argument(
        "cif_file", help="crystal structure (CIF; other formats with --reader pymatgen)"
    )
    parser.add_argument("toml_file", help="cif2qewan.toml with paths and parameters")
    parser.add_argument("--so", action="store_true", help="include spin-orbit coupling")
    parser.add_argument("--mag", action="store_true", help="ferromagnetic calculation")
    parser.add_argument(
        "--reader",
        choices=("cif2cell", "pymatgen"),
        default="cif2cell",
        help="how to read the structure (default: run cif2cell as in earlier versions)",
    )
    parser.add_argument(
        "--cif2cell-output",
        metavar="FILE",
        help=f"use this existing cif2cell output ({CIF2CELL_OUTPUT_NAME}) instead of running cif2cell; "
        "the CIF file is then not read",
    )
    parser.add_argument(
        "--output-dir",
        "-o",
        metavar="DIR",
        default=".",
        help="where to write the inputs (default: .)",
    )
    parser.add_argument(
        "--version", action="version", version=f"cif2qewan {__version__}"
    )
    return parser


def run(args: argparse.Namespace) -> List[str]:
    """Generate the inputs; returns the written paths relative to the output directory."""
    config = Config.from_toml(args.toml_file, so=args.so, mag=args.mag)
    output_dir = Path(args.output_dir)
    rendered = {}

    if args.cif2cell_output is not None:
        output = Cif2cellReader.read_output(args.cif2cell_output)
        structure, hints = output.structure, StructureHints.from_cif2cell(output)
        if (
            Path(args.cif2cell_output).resolve()
            != (output_dir / CIF2CELL_OUTPUT_NAME).resolve()
        ):
            rendered[CIF2CELL_OUTPUT_NAME] = output.text
    elif args.reader == "pymatgen":
        structure, hints = PymatgenReader().read(args.cif_file), StructureHints()
    else:
        reader = Cif2cellReader(
            config.cif2cell_path, k_resolution=config.scf_k_resolution
        )
        output = reader.run(args.cif_file)
        structure, hints = output.structure, StructureHints.from_cif2cell(output)
        rendered[CIF2CELL_OUTPUT_NAME] = (
            output.text
        )  # kept for reproducibility, as before

    plan = build_plan(structure, config, hints)
    rendered.update(render_plan(plan))
    written = write_rendered(rendered, output_dir)
    return [str(path.relative_to(output_dir)) for path in written]


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point of the ``cif2qewan`` console script."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        written = run(args)
    except Cif2qewanError as exc:
        print(f"cif2qewan: error: {exc}", file=sys.stderr)
        return 1
    print("cif2qewan: wrote " + ", ".join(written))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
