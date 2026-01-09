# Examples

This section contains practical examples demonstrating how to use MACE Inference for various computational tasks.

## Example Files

All examples are available in the [`examples/`](https://github.com/lichman0405/mace-inference/tree/main/examples) directory of the repository.

| Example | Description |
|---------|-------------|
| [Basic Usage](basic-usage.md) | Single-point, optimization, basic workflow |
| [Molecular Dynamics](molecular-dynamics.md) | NVT and NPT simulations |
| [Phonon Calculation](phonon-calculation.md) | Phonon DOS and thermal properties |
| [Adsorption Study](adsorption-study.md) | Gas adsorption in frameworks |
| [High Throughput](high-throughput.md) | Screening multiple structures |
| [D3 Correction](d3-correction.md) | Using DFT-D3 dispersion |
| [Batch Processing](batch-processing.md) | Processing multiple files |

## Running Examples

To run the examples, first clone the repository:

```bash
git clone https://github.com/lichman0405/mace-inference.git
cd mace-inference
```

Install the package:

```bash
pip install -e .
```

Run an example:

```bash
python examples/01_basic_usage.py
```

## Structure Files

Examples use CIF structure files from `examples/structures/`:

- `cu_fcc.cif` - FCC Copper
- `cu_paddlewheel.cif` - Cu paddlewheel cluster (MOF building block)
- `si_diamond.cif` - Diamond cubic Silicon

## Output

Examples that generate output save files to `examples/output/`.

## Quick Start Example

Here's a minimal working example:

```python
from mace_inference import MACEInference
from ase.io import read

# Initialize
calc = MACEInference(model="medium", device="auto")

# Load structure
atoms = read("structure.cif")

# Calculate
result = calc.single_point(atoms)
print(f"Energy: {result['energy']:.4f} eV")

# Optimize
opt = calc.optimize(atoms, fmax=0.01)
print(f"Optimized energy: {opt['energy']:.4f} eV")
```
