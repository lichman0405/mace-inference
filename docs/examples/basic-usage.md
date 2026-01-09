# Example 1: Basic Usage

This example demonstrates the fundamental operations with MACE Inference.

## Source Code

See [`examples/01_basic_usage.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/01_basic_usage.py)

## What You'll Learn

- Initializing the MACE calculator
- Loading structures from CIF files
- Single-point energy calculation
- Structure optimization
- Saving results

## Code Walkthrough

### Import and Initialize

```python
from pathlib import Path
from ase.io import read, write
from mace_inference import MACEInference

# Initialize calculator
calc = MACEInference(model="medium", device="auto")
```

### Load Structure

```python
STRUCTURES_DIR = Path(__file__).parent / "structures"
atoms = read(str(STRUCTURES_DIR / "cu_fcc.cif"))

print(f"Structure: {atoms.get_chemical_formula()}")
print(f"Number of atoms: {len(atoms)}")
```

### Single-Point Calculation

```python
result = calc.single_point(atoms)

print(f"Energy: {result['energy']:.6f} eV")
print(f"Energy/atom: {result['energy_per_atom']:.6f} eV")
print(f"Max force: {result['max_force']:.6f} eV/Å")
```

### Structure Optimization

```python
opt_result = calc.optimize(
    atoms,
    fmax=0.01,           # Force convergence
    steps=500,           # Max iterations
    optimizer="BFGS"     # Algorithm
)

print(f"Converged: {opt_result['converged']}")
print(f"Final energy: {opt_result['energy']:.6f} eV")
print(f"Steps taken: {opt_result['steps']}")
```

### Full Cell Relaxation

```python
opt_cell = calc.optimize(
    atoms,
    fmax=0.01,
    relax_cell=True      # Also optimize cell vectors
)

# Compare volumes
initial_vol = atoms.get_volume()
final_vol = opt_cell['atoms'].get_volume()
print(f"Volume change: {(final_vol/initial_vol - 1)*100:.2f}%")
```

## Expected Output

```
============================================================
Example 1: Basic Usage
============================================================

1. Loading structure...
   Structure: Cu4
   Number of atoms: 4

2. Single-point calculation...
   Energy: -16.3376 eV
   Energy/atom: -4.0844 eV
   Max force: 0.0001 eV/Å

3. Structure optimization...
   Converged: True
   Steps: 1
   Final energy: -16.3376 eV

4. Cell relaxation...
   Volume change: +0.00%

✅ Example 1 completed successfully!
```

## Key Points

1. **Model Selection**: Use `"small"`, `"medium"`, or `"large"` for pretrained models
2. **Device Selection**: `"auto"` picks the best available (GPU if available)
3. **Force Units**: Forces are in eV/Å
4. **Convergence**: Check `converged` flag to verify optimization success
