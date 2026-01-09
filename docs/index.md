# MACE Inference

<div align="center" markdown>

**High-level Python interface for MACE machine learning interatomic potentials**

[![CI](https://github.com/lichman0405/mace-inference/workflows/CI/badge.svg)](https://github.com/lichman0405/mace-inference/actions)
[![PyPI version](https://img.shields.io/pypi/v/mace-inference.svg)](https://pypi.org/project/mace-inference/)
[![Downloads](https://img.shields.io/pypi/dm/mace-inference.svg)](https://pypi.org/project/mace-inference/)
[![Python Version](https://img.shields.io/pypi/pyversions/mace-inference.svg)](https://pypi.org/project/mace-inference/)
[![License](https://img.shields.io/github/license/lichman0405/mace-inference.svg)](https://github.com/lichman0405/mace-inference/blob/main/LICENSE)

</div>

---

## What is MACE Inference?

**MACE Inference** is a user-friendly Python library that provides a high-level interface for performing atomistic simulations using [MACE](https://github.com/ACEsuit/mace) machine learning interatomic potentials.

Instead of writing complex setup code, you can perform common computational chemistry tasks with just a few lines:

```python
from mace_inference import MACEInference
from ase.io import read

# Initialize calculator
calc = MACEInference(model="medium", device="auto")

# Load structure and calculate
atoms = read("structure.cif")
result = calc.single_point(atoms)

print(f"Energy: {result['energy']:.4f} eV")
print(f"Max force: {result['max_force']:.4f} eV/Å")
```

## Key Features

<div class="grid cards" markdown>

- :material-lightning-bolt:{ .lg .middle } **Easy to Use**

    ---

    Simple, intuitive API for common atomistic simulation tasks

- :material-atom:{ .lg .middle } **Comprehensive Tasks**

    ---

    Single-point, optimization, MD, phonons, mechanics, adsorption

- :material-rocket-launch:{ .lg .middle } **GPU Accelerated**

    ---

    Automatic device detection with CUDA/MPS support

- :material-puzzle:{ .lg .middle } **D3 Correction**

    ---

    Optional DFT-D3 dispersion correction for van der Waals systems

</div>

## Supported Tasks

| Task | Method | Description |
|------|--------|-------------|
| Single-point | `single_point()` | Calculate energy, forces, stress |
| Optimization | `optimize()` | Geometry relaxation with BFGS/FIRE |
| NVT MD | `run_md()` | Constant temperature dynamics |
| NPT MD | `run_md()` | Constant pressure dynamics |
| Phonons | `phonon()` | Phonon spectrum & thermal properties |
| Bulk Modulus | `bulk_modulus()` | Elastic properties via EOS fitting |
| Adsorption | `adsorption_energy()` | Gas-framework binding energy |
| Coordination | `coordination()` | Coordination number analysis |

## Quick Links

<div class="grid cards" markdown>

- [:material-download: **Installation**](getting-started/installation.md)

    Get started by installing MACE Inference

- [:material-play: **Quick Start**](getting-started/quickstart.md)

    Learn the basics with a simple example

- [:material-book-open-variant: **User Guide**](user-guide/overview.md)

    Detailed guides for each feature

- [:material-api: **API Reference**](api/core.md)

    Complete API documentation

</div>

## Requirements

- Python 3.9+
- ASE (Atomic Simulation Environment)
- mace-torch
- phonopy (for phonon calculations)
- torch-dftd (optional, for D3 correction)

## License

MACE Inference is released under the MIT License. See [LICENSE](https://github.com/lichman0405/mace-inference/blob/main/LICENSE) for details.
