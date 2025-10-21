#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test runner script for cif2qewan.

This script provides a convenient way to run all tests with different
configurations and options.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def run_tests(test_type="all", verbose=True, coverage=False, parallel=False):
    """
    Run tests with specified options.

    Parameters
    ----------
    test_type : str
        Type of tests to run ('all', 'unit', 'integration', 'functional')
    verbose : bool
        Whether to run tests in verbose mode
    coverage : bool
        Whether to generate coverage report
    parallel : bool
        Whether to run tests in parallel
    """
    # Get project root directory
    project_root = Path(__file__).parent

    # Build pytest command
    cmd = ["python", "-m", "pytest"]

    # Add test path
    cmd.append("tests/")

    # Add verbosity
    if verbose:
        cmd.extend(["-v", "--tb=short"])

    # Add coverage if requested
    if coverage:
        cmd.extend(
            ["--cov=cif2qewan", "--cov-report=html", "--cov-report=term-missing"]
        )

    # Add parallel execution if requested
    if parallel:
        cmd.extend(["-n", "auto"])

    # Add test type filtering
    if test_type == "unit":
        cmd.extend(["-m", "unit"])
    elif test_type == "integration":
        cmd.extend(["-m", "integration"])
    elif test_type == "functional":
        cmd.extend(["-m", "functional"])
    elif test_type == "fast":
        cmd.extend(["-m", "not slow"])

    # Add markers for external dependencies
    cmd.extend(
        [
            "-m",
            "not requires_qe and not requires_wannier90 and not requires_cif2cell and not requires_api",
        ]
    )

    print(f"Running tests with command: {' '.join(cmd)}")
    print(f"Project root: {project_root}")

    # Change to project root directory
    os.chdir(project_root)

    # Run tests
    try:
        result = subprocess.run(cmd, check=True)
        print("\n[SUCCESS] All tests passed!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] Tests failed with exit code {e.returncode}")
        return False
    except Exception as e:
        print(f"\n[ERROR] Error running tests: {e}")
        return False


def main():
    """Main function for test runner."""
    parser = argparse.ArgumentParser(description="Run cif2qewan tests")
    parser.add_argument(
        "--type",
        choices=["all", "unit", "integration", "functional", "fast"],
        default="all",
        help="Type of tests to run",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=True,
        help="Run tests in verbose mode",
    )
    parser.add_argument(
        "--coverage", "-c", action="store_true", help="Generate coverage report"
    )
    parser.add_argument(
        "--parallel", "-p", action="store_true", help="Run tests in parallel"
    )
    parser.add_argument(
        "--install-deps", action="store_true", help="Install test dependencies"
    )

    args = parser.parse_args()

    # Install dependencies if requested
    if args.install_deps:
        print("Installing test dependencies...")
        try:
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "-r", "requirements-test.txt"],
                check=True,
            )
            print("✅ Test dependencies installed successfully!")
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to install test dependencies: {e}")
            return False

    # Run tests
    success = run_tests(
        test_type=args.type,
        verbose=args.verbose,
        coverage=args.coverage,
        parallel=args.parallel,
    )

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
