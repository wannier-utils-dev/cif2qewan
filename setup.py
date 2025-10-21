#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Setup script for cif2qewan package.

This script handles the installation and distribution of the cif2qewan package.
It uses setuptools for building and installing the package.
"""

from setuptools import setup, find_packages
import os

# Read the README file for long description
def read_readme():
    """Read the README file for long description."""
    readme_path = os.path.join(os.path.dirname(__file__), "README.md")
    if os.path.exists(readme_path):
        with open(readme_path, "r", encoding="utf-8") as f:
            return f.read()
    return "A comprehensive Python toolkit for generating Quantum ESPRESSO and Wannier90 input files from CIF structures"

# Read requirements from pyproject.toml or create minimal requirements
def read_requirements():
    """Read requirements from pyproject.toml or create minimal requirements."""
    requirements = [
        "numpy>=1.19.0",
        "pandas>=1.3.0", 
        "pymatgen>=2022.0.0",
        "toml>=0.10.0",
        "docopt>=0.6.0",
        "matplotlib>=3.3.0",
    ]
    return requirements

if __name__ == "__main__":
    setup(
        name="cif2qewan",
        version="1.0.0",
        description="A comprehensive Python toolkit for generating Quantum ESPRESSO and Wannier90 input files from CIF structures",
        long_description=read_readme(),
        long_description_content_type="text/markdown",
        author="cif2qewan developers",
        author_email="cif2qewan@example.com",
        url="https://github.com/wannier-utils-dev/cif2qewan",
        project_urls={
            "Homepage": "https://github.com/wannier-utils-dev/cif2qewan",
            "Documentation": "https://github.com/wannier-utils-dev/cif2qewan#readme",
            "Repository": "https://github.com/wannier-utils-dev/cif2qewan.git",
            "Issues": "https://github.com/wannier-utils-dev/cif2qewan/issues",
        },
        packages=find_packages(),
        package_data={
            "cif2qewan": [
                "*.toml",
                "*.csv", 
                "*.sh",
            ],
        },
        include_package_data=True,
        install_requires=read_requirements(),
        extras_require={
            "dev": [
                "pytest>=6.0",
                "pytest-cov>=2.0",
                "flake8>=3.8",
                "black>=21.0",
                "isort>=5.0",
                "mypy>=0.800",
            ],
            "test": [
                "pytest>=6.0",
                "pytest-cov>=2.0",
            ],
            "docs": [
                "sphinx>=4.0",
                "sphinx-rtd-theme>=1.0",
                "sphinx-autodoc-typehints>=1.0",
            ],
        },
        entry_points={
            "console_scripts": [
                "cif2qewan=cif2qewan.cif2qewan:main",
                "band_comp=cif2qewan.band_comp:main", 
                "wannier_conv=cif2qewan.wannier_conv:main",
            ],
        },
        python_requires=">=3.7",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: 3.12",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Chemistry",
        ],
        keywords=[
            "quantum-espresso",
            "wannier90",
            "dft", 
            "band-structure",
            "cif",
            "materials-science",
            "first-principles",
        ],
        zip_safe=False,
    )
