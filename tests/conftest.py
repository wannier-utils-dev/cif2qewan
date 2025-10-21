# -*- coding: utf-8 -*-
"""
Pytest configuration and fixtures for cif2qewan tests.

This module provides common fixtures and configuration for all tests.
"""

import os
import tempfile
import pytest
import shutil
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for testing."""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path)


@pytest.fixture
def example_cif_file(temp_dir):
    """Create a sample CIF file for testing."""
    cif_file = os.path.join(temp_dir, "test.cif")
    cif_content = """data_test
_cell_length_a    2.863
_cell_length_b    2.863
_cell_length_c    2.863
_cell_angle_alpha 90.0
_cell_angle_beta  90.0
_cell_angle_gamma 90.0
_space_group_name_H-M 'I m -3 m'
_space_group_IT_number 229
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Fe1 Fe 0.0 0.0 0.0
"""
    with open(cif_file, 'w') as f:
        f.write(cif_content)
    return cif_file


@pytest.fixture
def example_toml_file(temp_dir):
    """Create a sample TOML configuration file for testing."""
    toml_file = os.path.join(temp_dir, "test.toml")
    toml_content = """cif2cell_path = "/usr/bin/cif2cell"
pseudo_dir = "/tmp/pseudos"
pp_list_path = "/tmp/pp_list.csv"
scf_k_resolution = 0.15
degauss = 0.01

[pw2wan]
write_unk = ".true."
"""
    with open(toml_file, 'w') as f:
        f.write(toml_content)
    return toml_file


@pytest.fixture
def example_csv_file(temp_dir):
    """Create a sample pseudopotential CSV file for testing."""
    csv_file = os.path.join(temp_dir, "test_pp.csv")
    csv_content = """atom,pp_file_name,nexclude,orbitals,ecutwfc,ecutrho
Fe,Fe.pbe-n-rrkjus_psl.1.0.0.UPF,0,dsp,40.0,200.0
O,O.pbe-n-rrkjus_psl.1.0.0.UPF,0,sp,40.0,200.0
"""
    with open(csv_file, 'w') as f:
        f.write(csv_content)
    return csv_file


@pytest.fixture
def project_root():
    """Get the project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def example_dir(project_root):
    """Get the examples directory."""
    return project_root / "examples"
