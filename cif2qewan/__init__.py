"""
cif2qewan: Quantum ESPRESSO and Wannier90 input generation from CIF files.

Modules
-------
cli
    The ``cif2qewan`` command: reader -> workflow builder -> renderers -> files.
structure, qe, wannier90, workflow
    Structure readers, input models, renderers and the scientific policy.
cif2qewan
    Deprecated 0.2.x implementation, kept for existing scripts.
band_comp
    Compare DFT and Wannier90 band structures.
wannier_conv
    Check the accuracy of the Wannier90 interpolation.
"""

__version__ = "0.2.0"
__all__ = ["__version__"]
