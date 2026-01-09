# MACE Inference

<div align="center">

[![CI](https://github.com/lichman0405/mace-inference/workflows/CI/badge.svg)](https://github.com/lichman0405/mace-inference/actions)
[![PyPI version](https://img.shields.io/pypi/v/mace-inference.svg)](https://pypi.org/project/mace-inference/)
[![Downloads](https://img.shields.io/pypi/dm/mace-inference.svg)](https://pypi.org/project/mace-inference/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/mace-inference.svg)](https://pypi.org/project/mace-inference/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**High-level Python library for MACE machine learning force field inference tasks**

[Documentation](https://lichman0405.github.io/mace-inference/) • [Installation](#-installation) • [Quick Start](#-quick-start) • [Examples](#-examples)

</div>

---

## 🌟 Features

- **🚀 Simple API**: High-level interface for common inference tasks
- **🔌 Flexible**: Works with all MACE models (MACE-MP-0, MACE-MPA-0, MACE-OMAT-0, etc.)
- **⚡ GPU Support**: Automatic device detection (CPU/CUDA)
- **🧪 Complete Workflow**: Single-point energy, optimization, MD, phonon, adsorption, etc.
- **🔧 CLI Tools**: Command-line interface for non-programmers
- **📦 Pure Python**: Built on ASE and Phonopy - no extra software required

---

## 📦 Installation

### Basic Installation (CPU)

```bash
pip install mace-inference
```

### GPU Support

```bash
pip install mace-inference[gpu]
```

### With D3 Dispersion Correction

```bash
pip install mace-inference[d3]
```

### Development Installation

```bash
git clone https://github.com/lichman0405/mace-inference.git
cd mace-inference
pip install -e ".[all]"
```

---

## 🚀 Quick Start

### Python API

```python
from ase.io import read
from mace_inference import MACEInference

# Initialize MACE calculator
calc = MACEInference(model="medium", device="auto")

# Load structure
atoms = read("structure.cif")

# Single-point energy calculation
result = calc.single_point(atoms)
print(f"Energy: {result['energy']:.4f} eV")
print(f"Max Force: {result['max_force']:.4f} eV/Å")

# Structure optimization
optimized = calc.optimize(atoms, fmax=0.05)  # Returns Atoms object
# Save with: from ase.io import write; write("optimized.cif", optimized)

# Molecular dynamics (NVT)
trajectory = calc.run_md(
    atoms, 
    ensemble="nvt", 
    temperature_K=300, 
    steps=1000
)

# Phonon calculation
phonon_result = calc.phonon(atoms, supercell_matrix=[2, 2, 2])

# Adsorption energy
result = calc.adsorption_energy(
    framework=mof,
    adsorbate="CO2",
    site_position=[10.0, 10.0, 10.0]
)
print(f"Adsorption Energy: {result['E_ads']:.3f} eV")

# Batch processing with progress callback
def progress(current, total):
    print(f"Processing {current}/{total}")

results = calc.batch_single_point(
    structures,
    progress_callback=progress
)
```

### Command Line Interface

```bash
# Single-point energy
mace-infer energy structure.cif --model medium

# Structure optimization
mace-infer optimize structure.cif --fmax 0.05 --output optimized.cif

# Molecular dynamics
mace-infer md structure.cif --ensemble nvt --temp 300 --steps 10000

# Phonon calculation
mace-infer phonon structure.cif --supercell 2 2 2 --temp-range 0 1000 10

# Bulk modulus
mace-infer bulk-modulus structure.cif

# Adsorption energy
mace-infer adsorption mof.cif --gas CO2 --site 10.0 10.0 10.0
```

---

## 📚 Documentation

Full documentation is available at **[lichman0405.github.io/mace-inference](https://lichman0405.github.io/mace-inference/)**

- **[Installation Guide](https://lichman0405.github.io/mace-inference/getting-started/installation/)** - Detailed installation instructions
- **[Quick Start](https://lichman0405.github.io/mace-inference/getting-started/quickstart/)** - Getting started guide
- **[User Guide](https://lichman0405.github.io/mace-inference/user-guide/overview/)** - Detailed usage guides
- **[API Reference](https://lichman0405.github.io/mace-inference/api/core/)** - Complete API documentation
- **[Examples](examples/)** - Usage examples and tutorials
- **[Changelog](CHANGELOG.md)** - Version history and updates

### Available Tasks

| Task | Method | Description |
|------|--------|-------------|
| **Single-Point Energy** | `single_point()` | Calculate energy, forces, stress |
| **Structure Optimization** | `optimize()` | Minimize energy (atoms + cell) |
| **Molecular Dynamics** | `run_md()` | NVT/NPT simulations |
| **Phonon Calculation** | `phonon()` | Phonon dispersion & thermodynamics |
| **Bulk Modulus** | `bulk_modulus()` | Equation of state fitting |
| **Coordination Analysis** | `coordination()` | Metal-ligand bond lengths |
| **Adsorption Energy** | `adsorption_energy()` | Gas molecule binding energy |
| **Batch Single-Point** | `batch_single_point()` | Process multiple structures |
| **Batch Optimize** | `batch_optimize()` | Batch structure optimization |

### Supported MACE Models

- **MACE-MP-0a/0b3** (Materials Project)
- **MACE-MPA-0** (MPtraj + Alexandria)
- **MACE-OMAT-0** (OMat24)
- **MACE-MATPES-r2SCAN-0** (MatPES + MPtraj)
- Custom fine-tuned models

---

## 📖 Examples

See the [examples/](examples/) directory for detailed usage:

- [01_basic_usage.py](examples/01_basic_usage.py) - Basic single-point and optimization
- [02_molecular_dynamics.py](examples/02_molecular_dynamics.py) - NVT/NPT MD simulations
- [03_phonon_calculation.py](examples/03_phonon_calculation.py) - Phonon and thermal properties
- [04_adsorption_study.py](examples/04_adsorption_study.py) - Gas adsorption in MOFs
- [05_high_throughput.py](examples/05_high_throughput.py) - Batch processing
- [06_d3_correction.py](examples/06_d3_correction.py) - Using D3 dispersion correction
- [07_batch_processing.py](examples/07_batch_processing.py) - Batch API and progress callbacks

---

## 🧪 Testing

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=mace_inference --cov-report=html
```

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **MACE**: [ACEsuit/mace](https://github.com/ACEsuit/mace)
- **ASE**: [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
- **Phonopy**: [Phonon calculations](https://phonopy.github.io/phonopy/)

---

## 📮 Citation

If you use this library in your research, please cite:

```bibtex
@software{mace_inference,
  title = {MACE Inference: High-level Python library for MACE force field inference},
  author = {Shibo Li},
  year = {2026},
  url = {https://github.com/lichman0405/mace-inference}
}
```

And the original MACE paper:

```bibtex
@inproceedings{Batatia2022mace,
  title={{MACE}: Higher Order Equivariant Message Passing Neural Networks for Fast and Accurate Force Fields},
  author={Ilyes Batatia and David Peter Kovacs and Gregor N. C. Simm and Christoph Ortner and Gabor Csanyi},
  booktitle={Advances in Neural Information Processing Systems},
  year={2022}
}
```
