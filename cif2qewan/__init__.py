"""
cif2qewan: Quantum ESPRESSO and Wannier90 input generation from CIF files.

Modules
-------
cli
    The ``cif2qewan`` command: reader -> workflow builder -> renderers -> files.
structure, qe, wannier90, workflow
    Structure readers, input models, renderers and the scientific policy.
cif2qewan
    Deprecated alias of ``cli`` for ``python -m cif2qewan.cif2qewan``.
band_comp
    Compare DFT and Wannier90 band structures.
wannier_conv
    Check the accuracy of the Wannier90 interpolation.
"""

__version__ = "0.3.0rc1"
__all__ = ["__version__"]
