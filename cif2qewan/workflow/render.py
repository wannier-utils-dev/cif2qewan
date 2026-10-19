"""Render every input of a CalculationPlan to text, keyed by relative path."""

from __future__ import annotations

from typing import Dict

from cif2qewan.exceptions import InputModelError
from cif2qewan.qe.model import NamelistInput, PwInput
from cif2qewan.qe.writer import render as render_qe
from cif2qewan.wannier90.model import Wannier90Input
from cif2qewan.wannier90.writer import render_win
from cif2qewan.workflow.model import CalculationPlan


def render_plan(plan: CalculationPlan) -> Dict[str, str]:
    """Relative path -> file content for every input of the plan."""
    rendered = {}
    for path, model in plan.files().items():
        if isinstance(model, (PwInput, NamelistInput)):
            rendered[path] = render_qe(model)
        elif isinstance(model, Wannier90Input):
            rendered[path] = render_win(model)
        else:
            raise InputModelError(f"cannot render {path}: {model!r}")
    return rendered
