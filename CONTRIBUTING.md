# Contributing to cif2qewan

Thank you for your interest in contributing to cif2qewan! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Contributing Guidelines](#contributing-guidelines)
- [Pull Request Process](#pull-request-process)
- [Issue Reporting](#issue-reporting)
- [Coding Standards](#coding-standards)
- [Testing](#testing)
- [Documentation](#documentation)

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/). By participating, you agree to uphold this code.

## Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/your-username/cif2qewan.git
   cd cif2qewan
   ```
3. **Add the upstream repository**:
   ```bash
   git remote add upstream https://github.com/wannier-utils-dev/cif2qewan.git
   ```

## Development Setup

### Prerequisites

- Python 3.7+
- Git
- Virtual environment (recommended)

### Setup Development Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e .

# Install development dependencies
pip install -e ".[dev,test,docs]"

# Install pre-commit hooks (optional)
pip install pre-commit
pre-commit install
```

### Development Dependencies

```bash
# Install all development dependencies
pip install -r requirements.txt
pip install pytest pytest-cov flake8 black isort mypy sphinx
```

## Contributing Guidelines

### Types of Contributions

- **Bug fixes**: Fix issues in the code
- **New features**: Add new functionality
- **Documentation**: Improve documentation
- **Tests**: Add or improve test coverage
- **Performance**: Optimize code performance
- **Refactoring**: Improve code structure

### Workflow

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**:
   - Write code following the coding standards
   - Add tests for new functionality
   - Update documentation as needed

3. **Test your changes**:
   ```bash
   # Run tests
   pytest

   # Run linting
   flake8 cif2qewan/
   black --check cif2qewan/
   isort --check-only cif2qewan/

   # Type checking
   mypy cif2qewan/
   ```

4. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Add: brief description of changes"
   ```

5. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Create a Pull Request** on GitHub

## Pull Request Process

### Before Submitting

- [ ] Code follows the project's coding standards
- [ ] Tests pass locally
- [ ] Documentation is updated
- [ ] Commit messages are clear and descriptive
- [ ] Branch is up to date with main/develop

### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Testing
- [ ] Tests pass locally
- [ ] New tests added for new functionality
- [ ] Manual testing completed

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No breaking changes (or clearly documented)
```

## Issue Reporting

### Before Creating an Issue

1. **Search existing issues** to avoid duplicates
2. **Check the documentation** for solutions
3. **Try the latest version** to see if the issue is already fixed

### Issue Template

```markdown
## Bug Report / Feature Request

### Description
Clear and concise description of the issue or feature request.

### Steps to Reproduce (for bugs)
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

### Expected Behavior
What you expected to happen.

### Actual Behavior
What actually happened.

### Environment
- OS: [e.g., Ubuntu 20.04]
- Python version: [e.g., 3.8.5]
- cif2qewan version: [e.g., 1.0.0]

### Additional Context
Any other context about the problem.
```

## Coding Standards

### Python Code Style

- **PEP 8**: Follow Python PEP 8 style guidelines
- **Line length**: Maximum 88 characters (Black formatter)
- **Imports**: Organize imports (standard library, third-party, local)
- **Type hints**: Use type hints for all functions
- **Docstrings**: Use NumPy-style docstrings

### Code Formatting

```bash
# Format code with Black
black cif2qewan/

# Sort imports with isort
isort cif2qewan/

# Check style with flake8
flake8 cif2qewan/
```

### Documentation Standards

- **Docstrings**: All functions and classes must have docstrings
- **Comments**: Complex code should be commented
- **README**: Keep README.md updated
- **Type hints**: Use type hints for better documentation

### Example Code Style

```python
def calculate_energy(
    k_point: np.ndarray, 
    hamiltonian: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate band energies at a given k-point.
    
    Parameters
    ----------
    k_point : np.ndarray
        K-point coordinates (shape: [3]).
    hamiltonian : np.ndarray
        Hamiltonian matrix (shape: [n_bands, n_bands]).
        
    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Eigenvalues and eigenvectors.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    return eigenvalues, eigenvectors
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=cif2qewan --cov-report=html

# Run specific test file
pytest tests/test_cif2qewan.py

# Run with verbose output
pytest -v
```

### Writing Tests

- **Test coverage**: Aim for >80% coverage
- **Test naming**: Use descriptive test names
- **Test organization**: Group related tests in classes
- **Fixtures**: Use pytest fixtures for common setup

### Example Test

```python
import pytest
import numpy as np
from cif2qewan.cif2qewan import qe_wannier_in

class TestQEWannierIn:
    """Test cases for QEWannierIn class."""
    
    def test_initialization(self):
        """Test class initialization."""
        # Test initialization
        pass
    
    def test_band_calculation(self):
        """Test band structure calculation."""
        # Test band calculation
        pass
```

## Documentation

### Building Documentation

```bash
# Install documentation dependencies
pip install sphinx sphinx-rtd-theme

# Build documentation
cd docs/
make html
```

### Documentation Standards

- **API documentation**: Document all public functions and classes
- **Examples**: Include usage examples
- **Tutorials**: Provide step-by-step tutorials
- **Changelog**: Keep CHANGELOG.md updated

## Release Process

### Version Numbering

We use [Semantic Versioning](https://semver.org/):
- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality (backward compatible)
- **PATCH**: Bug fixes (backward compatible)

### Release Checklist

- [ ] All tests pass
- [ ] Documentation updated
- [ ] Changelog updated
- [ ] Version number updated
- [ ] Tag created
- [ ] Release notes written

## Getting Help

- **GitHub Issues**: For bug reports and feature requests
- **Discussions**: For questions and general discussion
- **Email**: For security issues (use private communication)

## Recognition

Contributors will be recognized in:
- **CONTRIBUTORS.md**: List of all contributors
- **Release notes**: Major contributors mentioned
- **Documentation**: Contributors acknowledged

Thank you for contributing to cif2qewan!
