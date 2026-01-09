# Changelog

All notable changes to MACE Inference will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **CI/CD**: GitHub Actions workflow for automated testing (Python 3.9, 3.10, 3.11)
- **Documentation**: MkDocs Material documentation with GitHub Pages deployment
- **API docs**: Auto-generated API reference using mkdocstrings
- Type hints improvements with `py.typed` marker
- Comprehensive test suite with pytest fixtures
- Batch processing API: `batch_single_point()` and `batch_optimize()` methods
- Progress callback support for long-running tasks
- Example structure files (CIF) for all examples

### Changed
- Updated README with GitHub Pages documentation links
- Enhanced type annotations across all modules
- Python version requirement: 3.9+

### Fixed
- API parameter mismatch issues in core functions
- All 7 examples updated with correct API usage

## [0.1.0] - 2024-01-XX

### Added
- Initial release of MACE Inference
- `MACEInference` main class with high-level API
- Single-point energy, force, and stress calculations
- Geometry optimization with BFGS and FIRE
- NVT molecular dynamics with Langevin thermostat
- NPT molecular dynamics with Berendsen barostat
- Phonon calculations via Phonopy
- Thermal property calculations (free energy, entropy, heat capacity)
- Bulk modulus via equation of state fitting
- Gas adsorption energy calculations
- Coordination number analysis
- DFT-D3 dispersion correction support via torch-dftd
- Command-line interface (CLI)
- Example scripts for all features
- Comprehensive test suite

### Supported Models
- Materials Project pretrained models: small, medium, large
- Custom MACE models (.model files)

### Supported Devices
- CPU
- CUDA (NVIDIA GPU)
- MPS (Apple Silicon)

---

## Version History

| Version | Date | Description |
|---------|------|-------------|
| 0.1.0 | 2024-01 | Initial release |
