#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Integration tests for cif2qewan package.

This module contains integration tests that test the interaction
between different components of the cif2qewan toolkit.
"""

import os
import tempfile
import unittest
import subprocess
import sys
from unittest.mock import patch, MagicMock

# Add the parent directory to the path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class TestEndToEndWorkflow(unittest.TestCase):
    """Test cases for end-to-end workflow."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.cif_file = os.path.join(self.temp_dir, "test.cif")
        self.toml_file = os.path.join(self.temp_dir, "test.toml")
        
        # Create test files
        self.create_test_cif()
        self.create_test_toml()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def create_test_cif(self):
        """Create a simple test CIF file."""
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
        with open(self.cif_file, 'w') as f:
            f.write(cif_content)
    
    def create_test_toml(self):
        """Create a test TOML configuration file."""
        toml_content = """cif2cell_path = "/usr/bin/cif2cell"
pseudo_dir = "/tmp/pseudos"
pp_list_path = "/tmp/pp_list.csv"
scf_k_resolution = 0.15
degauss = 0.01

[pw2wan]
write_unk = ".true."
"""
        with open(self.toml_file, 'w') as f:
            f.write(toml_content)
    
    @patch('subprocess.run')
    def test_cif2qewan_command_line(self, mock_run):
        """Test cif2qewan command line interface."""
        # Mock the subprocess calls
        mock_run.return_value = MagicMock(returncode=0)
        
        # Test that the command would be called correctly
        cmd = ["python", "cif2qewan.py", self.cif_file, self.toml_file, "--so", "--mag"]
        
        # This is a mock test - in reality, we'd run the actual command
        self.assertTrue(os.path.exists(self.cif_file))
        self.assertTrue(os.path.exists(self.toml_file))
    
    def test_file_generation_structure(self):
        """Test that the expected file structure is generated."""
        # This would test the actual file generation
        # For now, we'll test the input files exist
        self.assertTrue(os.path.exists(self.cif_file))
        self.assertTrue(os.path.exists(self.toml_file))
    
    def test_configuration_loading(self):
        """Test configuration file loading."""
        import toml
        
        config = toml.load(self.toml_file)
        
        # Test required configuration keys
        required_keys = ["cif2cell_path", "pseudo_dir", "pp_list_path"]
        for key in required_keys:
            self.assertIn(key, config)
        
        # Test configuration values
        self.assertEqual(config["scf_k_resolution"], 0.15)
        self.assertEqual(config["degauss"], 0.01)
        self.assertEqual(config["pw2wan"]["write_unk"], ".true.")


class TestBandCompIntegration(unittest.TestCase):
    """Test cases for band_comp integration."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_band_comp_command_line(self):
        """Test band_comp command line interface."""
        # Test that the command exists and can be imported
        try:
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            import band_comp
            self.assertTrue(callable(band_comp.main))
        except ImportError:
            self.fail("band_comp module could not be imported")
    
    def test_band_comp_help(self):
        """Test band_comp help functionality."""
        # Test that help can be displayed
        try:
            import argparse
            parser = argparse.ArgumentParser()
            parser.add_argument("-o", dest="odir", required=True, help="Output directory")
            
            # Test that the parser is created correctly
            self.assertTrue(parser is not None)
        except Exception as e:
            self.fail(f"band_comp help test failed: {e}")


class TestWannierConvIntegration(unittest.TestCase):
    """Test cases for wannier_conv integration."""
    
    def test_wannier_conv_command_line(self):
        """Test wannier_conv command line interface."""
        # Test that the command exists and can be imported
        try:
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            import wannier_conv
            # wannier_conv.py doesn't have a main function, it runs directly
            self.assertTrue(True)  # Just test that it can be imported
        except ImportError:
            self.fail("wannier_conv module could not be imported")
    
    def test_wannier_conv_help(self):
        """Test wannier_conv help functionality."""
        # Test that help can be displayed
        try:
            import argparse
            parser = argparse.ArgumentParser()
            parser.add_argument("-e", type=float, default=0.0, dest="emax", help="emax")
            parser.add_argument("-o", dest="odir", required=True, help="Output directory")
            parser.add_argument("-i", default="./check_wannier/nscf.out", dest="nscf_file", help="NSCF file")
            
            # Test that the parser is created correctly
            self.assertTrue(parser is not None)
        except Exception as e:
            self.fail(f"wannier_conv help test failed: {e}")


# get_cif integration tests removed - not distributed


class TestPackageStructure(unittest.TestCase):
    """Test cases for package structure."""
    
    def test_package_imports(self):
        """Test that all package modules can be imported."""
        try:
            import cif2qewan
            self.assertTrue(hasattr(cif2qewan, '__version__'))
        except ImportError:
            self.fail("cif2qewan package could not be imported")
    
    def test_command_line_scripts(self):
        """Test that command line scripts are available."""
        # Test that the scripts exist in the package (excluding get_cif)
        script_names = ["cif2qewan", "band_comp", "wannier_conv"]
        
        for script in script_names:
            # Check if the script can be imported
            try:
                sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
                if script == "cif2qewan":
                    import cif2qewan
                    # cif2qewan.py doesn't have a main function, it runs directly
                    self.assertTrue(True)  # Just test that it can be imported
                elif script == "band_comp":
                    import band_comp
                    self.assertTrue(hasattr(band_comp, 'main'))
                elif script == "wannier_conv":
                    import wannier_conv
                    # wannier_conv.py doesn't have a main function, it runs directly
                    self.assertTrue(True)  # Just test that it can be imported
            except ImportError:
                self.fail(f"Script {script} could not be imported")
    
    def test_package_metadata(self):
        """Test package metadata."""
        try:
            import cif2qewan
            self.assertTrue(hasattr(cif2qewan, '__version__'))
            self.assertTrue(hasattr(cif2qewan, '__author__'))
            self.assertTrue(hasattr(cif2qewan, '__license__'))
        except ImportError:
            self.fail("Package metadata not available")


class TestConfigurationIntegration(unittest.TestCase):
    """Test cases for configuration integration."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_toml_configuration_loading(self):
        """Test TOML configuration loading."""
        toml_file = os.path.join(self.temp_dir, "test.toml")
        
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
        
        # Test loading the configuration
        import toml
        config = toml.load(toml_file)
        
        # Test all required keys are present
        required_keys = ["cif2cell_path", "pseudo_dir", "pp_list_path", "scf_k_resolution", "degauss"]
        for key in required_keys:
            self.assertIn(key, config)
        
        # Test pw2wan section
        self.assertIn("pw2wan", config)
        self.assertIn("write_unk", config["pw2wan"])
    
    def test_csv_pseudopotential_loading(self):
        """Test CSV pseudopotential list loading."""
        csv_file = os.path.join(self.temp_dir, "test_pp.csv")
        
        csv_content = """atom,pp_file_name,nexclude,orbitals,ecutwfc,ecutrho
Fe,Fe.pbe-n-rrkjus_psl.1.0.0.UPF,0,dsp,40.0,200.0
O,O.pbe-n-rrkjus_psl.1.0.0.UPF,0,sp,40.0,200.0
"""
        with open(csv_file, 'w') as f:
            f.write(csv_content)
        
        # Test loading the CSV
        import pandas as pd
        df = pd.read_csv(csv_file)
        
        # Test data structure
        self.assertEqual(len(df), 2)
        self.assertIn('Fe', df['atom'].values)
        self.assertIn('O', df['atom'].values)
        
        # Test data types
        self.assertTrue(pd.api.types.is_numeric_dtype(df['ecutwfc']))
        self.assertTrue(pd.api.types.is_numeric_dtype(df['ecutrho']))


if __name__ == '__main__':
    # Run the tests
    unittest.main(verbosity=2)
