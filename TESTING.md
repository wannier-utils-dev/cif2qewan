# Testing Guide for cif2qewan

This document provides a comprehensive overview of the testing framework for the cif2qewan project, including test categories, execution methods, and coverage information.

## Table of Contents

- [Test Structure](#test-structure)
- [Test Categories](#test-categories)
- [Running Tests](#running-tests)
- [GitHub Actions Workflows](#github-actions-workflows)
- [Test Coverage](#test-coverage)
- [Troubleshooting](#troubleshooting)

## Test Structure

The testing framework is organized into three main categories:

```
tests/
├── __init__.py              # Package initialization
├── conftest.py              # Shared fixtures and configuration
├── test_cif2qewan.py        # Unit tests
├── test_integration.py      # Integration tests
└── test_functional.py       # Functional tests
```

## Test Categories

### 1. Unit Tests (`test_cif2qewan.py`)

**Purpose**: Test individual components in isolation

**Test Classes**:
- `TestQEWannierIn`: Tests for the main QE Wannier input generator
  - `test_initialization`: Class initialization
  - `test_config_parsing`: TOML configuration parsing
- `TestPseudoList`: Tests for pseudopotential list handling
  - `test_csv_parsing`: CSV file parsing
  - `test_pseudopotential_info`: Pseudopotential information extraction
- `TestBandComp`: Tests for band structure comparison
  - `test_get_ef_from_scfout`: Fermi energy extraction
  - `test_get_band_data`: Band data extraction
- `TestWannierConv`: Tests for Wannier90 convergence checking
  - `test_hamiltonian_initialization`: Hamiltonian class initialization
  - `test_nscfout_initialization`: NSCF output parsing
- `TestConfiguration`: Tests for configuration handling
  - `test_toml_config_validation`: TOML file validation
  - `test_invalid_config_handling`: Error handling for invalid configs
- `TestFileOperations`: Tests for file and directory operations
  - `test_directory_creation`: Directory creation
  - `test_file_creation`: File creation

**Total Unit Tests**: 12 tests

### 2. Integration Tests (`test_integration.py`)

**Purpose**: Test interaction between different components

**Test Classes**:
- `TestEndToEndWorkflow`: End-to-end workflow testing
  - `test_cif2qewan_command_line`: Command-line interface testing
  - `test_configuration_loading`: Configuration loading
  - `test_file_generation_structure`: File generation structure
- `TestBandCompIntegration`: Band comparison integration
  - `test_band_comp_command_line`: Command-line interface
  - `test_band_comp_help`: Help functionality
- `TestWannierConvIntegration`: Wannier90 convergence integration
  - `test_wannier_conv_command_line`: Command-line interface
  - `test_wannier_conv_help`: Help functionality
- `TestPackageStructure`: Package structure validation
  - `test_package_imports`: Module imports
  - `test_package_metadata`: Package metadata
  - `test_command_line_scripts`: Command-line scripts
- `TestConfigurationIntegration`: Configuration integration
  - `test_toml_configuration_loading`: TOML configuration
  - `test_csv_pseudopotential_loading`: CSV pseudopotential loading

**Total Integration Tests**: 12 tests

### 3. Functional Tests (`test_functional.py`)

**Purpose**: Test complete workflows and user scenarios

**Test Classes**:
- `TestFeExample`: Iron example workflow testing
  - `test_cif2qewan_execution`: cif2qewan execution
  - `test_generated_files_structure`: Generated file structure
  - `test_output_directory_structure`: Output directory structure
  - `test_scf_input_content`: SCF input file content
  - `test_wannier90_input_content`: Wannier90 input file content
- `TestBandCompFunctional`: Band comparison functionality
  - `test_band_comp_help`: Help functionality
  - `test_band_comp_output_directory_creation`: Output directory creation
- `TestWannierConvFunctional`: Wannier90 convergence functionality
  - `test_wannier_conv_help`: Help functionality
  - `test_wannier_conv_output_creation`: Output creation
- `TestPackageInstallation`: Package installation testing
  - `test_package_import`: Package import
  - `test_command_line_scripts`: Command-line scripts
  - `test_package_metadata`: Package metadata

**Total Functional Tests**: 12 tests

**Total Tests**: 36 tests

## Running Tests

### Quick Test Execution

```bash
# Run all tests
python run_tests.py

# Run specific test types
python run_tests.py --type unit
python run_tests.py --type integration
python run_tests.py --type functional
python run_tests.py --type fast

# Run with coverage
python run_tests.py --coverage

# Run in parallel
python run_tests.py --parallel
```

### Direct pytest Execution

```bash
# Run all tests
python -m pytest tests/

# Run specific test files
python -m pytest tests/test_cif2qewan.py
python -m pytest tests/test_integration.py
python -m pytest tests/test_functional.py

# Run with markers
python -m pytest tests/ -m unit
python -m pytest tests/ -m integration
python -m pytest tests/ -m functional
python -m pytest tests/ -m "not slow"

# Run with coverage
python -m pytest tests/ --cov=cif2qewan --cov-report=html
```

### Test Markers

The following markers are available for test categorization:

- `unit`: Unit tests
- `integration`: Integration tests
- `functional`: Functional tests
- `slow`: Slow-running tests
- `requires_qe`: Tests requiring Quantum ESPRESSO
- `requires_wannier90`: Tests requiring Wannier90
- `requires_cif2cell`: Tests requiring cif2cell
- `requires_api`: Tests requiring external API access

## GitHub Actions Workflows

### CI Workflow (`ci.yml`)
- **Python Versions**: 3.9, 3.10, 3.11, 3.12
- **Platforms**: Ubuntu Latest
- **Tests**: Unit, Integration, Functional
- **Coverage**: Code coverage reporting
- **Linting**: flake8, black, isort, mypy

### Test Matrix (`test-matrix.yml`)
- **Python Versions**: 3.9, 3.10, 3.11, 3.12
- **Platforms**: Ubuntu, macOS, Windows
- **Tests**: Cross-platform compatibility
- **Schedule**: Weekly execution

### Compatibility (`compatibility.yml`)
- **Python Versions**: 3.9, 3.10, 3.11, 3.12
- **Platforms**: Ubuntu 22.04/24.04, macOS 12/13, Windows 2022
- **Tests**: Platform-specific compatibility

### Lint (`lint.yml`)
- **Tools**: flake8, black, isort, mypy
- **Checks**: Code style, formatting, type checking
- **Standards**: PEP 8, line length 88, complexity limits

### Validate (`validate.yml`)
- **Package Structure**: Required files, directories
- **Examples**: PseudoDojo/PSLibrary examples
- **Configuration**: TOML file validation
- **Python Code**: Syntax checking

### Security (`security.yml`)
- **Tools**: safety, bandit
- **Checks**: Vulnerability scanning, security linting
- **Schedule**: Daily execution

### Quality (`quality.yml`)
- **Tools**: radon, xenon, vulture
- **Metrics**: Complexity, maintainability, dead code
- **Reports**: Quality metrics and analysis

### Performance (`performance.yml`)
- **Tools**: memory-profiler, psutil
- **Metrics**: Memory usage, execution time
- **Benchmarks**: Performance regression detection

### Coverage (`coverage.yml`)
- **Tool**: pytest-cov
- **Reports**: HTML, XML, terminal
- **Upload**: Codecov integration

## Test Coverage

### Current Coverage
- **Package Coverage**: 23%
- **Test Coverage**: 100% (all tests pass)
- **Success Rate**: 36/36 tests passing

### Coverage Goals
- **Target**: 80% package coverage
- **Critical Paths**: Core functionality
- **Edge Cases**: Error handling and boundary conditions

### Coverage Reports
- **HTML Report**: `htmlcov/index.html`
- **XML Report**: `coverage.xml`
- **Terminal Report**: Console output with missing lines

## Troubleshooting

### Common Issues

#### 1. Import Errors
```bash
# Solution: Install package in development mode
pip install -e .
```

#### 2. Missing Dependencies
```bash
# Solution: Install test dependencies
pip install -r requirements-test.txt
```

#### 3. Permission Errors
```bash
# Solution: Check file permissions
chmod +x *.py *.sh
```

#### 4. Unicode Errors (Windows)
```bash
# Solution: Use UTF-8 encoding
set PYTHONIOENCODING=utf-8
```

### Test Debugging

#### Verbose Output
```bash
python -m pytest tests/ -v -s
```

#### Specific Test Execution
```bash
python -m pytest tests/test_cif2qewan.py::TestQEWannierIn::test_initialization -v
```

#### Coverage Analysis
```bash
python -m pytest tests/ --cov=cif2qewan --cov-report=term-missing
```

### Continuous Integration

#### Local CI Simulation
```bash
# Run all CI checks locally
python run_tests.py --type fast
flake8 .
black --check .
isort --check-only .
mypy .
```

#### GitHub Actions Debugging
- Check Actions tab in GitHub repository
- Review workflow logs for specific errors
- Test locally with same Python versions

## Test Development

### Adding New Tests

1. **Unit Tests**: Add to `test_cif2qewan.py`
2. **Integration Tests**: Add to `test_integration.py`
3. **Functional Tests**: Add to `test_functional.py`

### Test Naming Convention

```python
def test_<functionality>_<condition>():
    """Test <functionality> when <condition>."""
    pass
```

### Test Documentation

- Use descriptive test names
- Include docstrings for complex tests
- Add comments for test setup and assertions
- Document expected behavior and edge cases

## Performance Testing

### Benchmark Tests
- **Import Time**: Module loading performance
- **Memory Usage**: RAM consumption monitoring
- **Execution Time**: Function call timing
- **File I/O**: Disk operation performance

### Load Testing
- **Large Files**: Processing large CIF files
- **Batch Operations**: Multiple file processing
- **Memory Limits**: Resource usage under stress

## Security Testing

### Vulnerability Scanning
- **Dependencies**: Known security issues
- **Code Analysis**: Static security analysis
- **Input Validation**: Malicious input handling
- **File Operations**: Path traversal protection

### Best Practices
- Regular dependency updates
- Security-focused code reviews
- Automated vulnerability scanning
- Secure coding guidelines

---

For more information about testing, see the [Contributing Guide](CONTRIBUTING.md) or [GitHub Issues](https://github.com/wannier-utils-dev/cif2qewan/issues).
