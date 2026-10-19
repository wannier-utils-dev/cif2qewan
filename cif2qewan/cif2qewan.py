"""Deprecated entry point kept for ``python -m cif2qewan.cif2qewan``.

The 0.2.x implementation (``qe_wannier_in``, ``pseudo_list``) was removed in
0.3.0; the command line lives in :mod:`cif2qewan.cli` and the input
generation in :mod:`cif2qewan.workflow`. This module only forwards to
:func:`cif2qewan.cli.main` and will be removed in 0.4.0.
"""

from __future__ import annotations

import sys
import warnings
from typing import Optional, Sequence


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the ``cif2qewan`` command (deprecated alias of :func:`cif2qewan.cli.main`)."""
    warnings.warn(
        "cif2qewan.cif2qewan.main is deprecated since 0.3.0 and will be removed in 0.4.0; "
        "use cif2qewan.cli.main or the cif2qewan command",
        DeprecationWarning,
        stacklevel=2,
    )
    from cif2qewan.cli import main as cli_main

    return cli_main(argv)


if __name__ == "__main__":
    from cif2qewan.cli import main as cli_main

    sys.exit(cli_main())
