# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Type hints improvements with `py.typed` marker
- Comprehensive test suite with pytest fixtures
- CLI tests using Click's CliRunner
- Core functionality tests with mock calculators
- Type definitions module (`types.py`) with TypedDict and Protocol classes
- `calculate_elastic_constants()` function for full elastic tensor calculation
- `find_adsorption_sites()` function for automatic void space detection
- Batch processing API: `batch_single_point()` and `batch_optimize()` methods
- Progress callback support for long-running tasks (MD simulations, batch operations)
- New example: D3 dispersion correction usage (`06_d3_correction.py`)
- New example: Batch processing and progress callbacks (`07_batch_processing.py`)
- Examples README with overview of all examples

### Changed
- Version management now uses `importlib.metadata` for single source of truth
- CLI version option dynamically reads package version
- Improved `utils/__init__.py` exports with complete API surface
- Python version requirement upgraded from 3.8 to 3.9+
- Enhanced type annotations across all modules
- Mypy configuration added for static type checking

### Fixed
- Missing `py.typed` marker file for type checker support
- Missing `tests/__init__.py` for proper test discovery
- Incomplete utility function exports in `utils/__init__.py`

## [0.1.0] - 2026-01-08

### Added
- Initial release of MACE Inference library
- `MACEInference` class with unified API for MACE force field inference
- Single-point energy calculation (`single_point()`)
- Structure optimization (`optimize()`) with LBFGS, BFGS, FIRE optimizers
- Molecular dynamics simulations (`run_md()`) with NVT and NPT ensembles
- Phonon calculations (`phonon()`) using Phonopy integration
- Thermal properties calculation from phonon data
- Bulk modulus calculation (`bulk_modulus()`) via equation of state fitting
- Gas adsorption energy calculation (`adsorption_energy()`)
- Coordination environment analysis (`coordination()`)
- DFT-D3 dispersion correction support via torch-dftd
- Automatic device detection (CPU/CUDA)
- Command-line interface (`mace-infer`) with subcommands:
  - `energy` - Single-point energy calculation
  - `optimize` - Structure optimization
  - `md` - Molecular dynamics simulation
  - `phonon` - Phonon calculation
  - `bulk-modulus` - Bulk modulus calculation
  - `adsorption` - Adsorption energy calculation
  - `info` - Environment information
- Example scripts for common use cases
- Comprehensive documentation

### Supported MACE Models
- MACE-MP-0 (small, medium, large)
- Custom fine-tuned models from file

### Dependencies
- mace-torch >= 0.3.10
- ase >= 3.22.0
- phonopy >= 2.20.0
- numpy >= 1.20.0
- click >= 8.0.0

[Unreleased]: https://github.com/lichman0405/mace-inference/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/lichman0405/mace-inference/releases/tag/v0.1.0
