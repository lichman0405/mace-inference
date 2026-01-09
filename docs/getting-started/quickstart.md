# Quick Start

This guide will get you started with MACE Inference in minutes.

## Basic Workflow

The typical workflow with MACE Inference consists of three steps:

1. **Initialize** the calculator
2. **Load** your atomic structure
3. **Run** calculations

```python
from mace_inference import MACEInference
from ase.io import read

# 1. Initialize calculator
calc = MACEInference(model="medium", device="auto")

# 2. Load structure
atoms = read("your_structure.cif")

# 3. Run calculation
result = calc.single_point(atoms)
print(f"Energy: {result['energy']:.4f} eV")
```

## Available Models

MACE Inference supports pretrained models from the Materials Project:

| Model | Size | Accuracy | Speed |
|-------|------|----------|-------|
| `small` | ~5M params | Good | Fast |
| `medium` | ~10M params | Better | Medium |
| `large` | ~25M params | Best | Slower |

```python
# Use different models
calc_fast = MACEInference(model="small")
calc_accurate = MACEInference(model="large")

# Or use a custom model
calc_custom = MACEInference(model="/path/to/your_model.model")
```

## Device Selection

```python
# Automatic device selection (recommended)
calc = MACEInference(model="medium", device="auto")

# Force CPU
calc = MACEInference(model="medium", device="cpu")

# Use specific GPU
calc = MACEInference(model="medium", device="cuda:0")

# Use Apple Silicon GPU
calc = MACEInference(model="medium", device="mps")
```

## Common Tasks

### Single-Point Calculation

Calculate energy, forces, and stress:

```python
result = calc.single_point(atoms)

print(f"Energy: {result['energy']:.4f} eV")
print(f"Energy/atom: {result['energy_per_atom']:.4f} eV")
print(f"Max force: {result['max_force']:.4f} eV/Å")
print(f"Forces shape: {result['forces'].shape}")
print(f"Stress: {result['stress']}")  # Voigt notation
```

### Structure Optimization

Relax the geometry to minimize forces:

```python
optimized = calc.optimize(
    atoms,
    fmax=0.01,            # Force convergence (eV/Å)
    steps=500,            # Max iterations
    optimizer="LBFGS",    # LBFGS, BFGS or FIRE
    optimize_cell=False,  # Also optimize cell
    output="relaxed.cif"  # Save result to file
)

# optimized is an Atoms object
print(f"Final energy: {optimized.get_potential_energy():.4f} eV")
```

### Molecular Dynamics

Run NVT or NPT simulations:

```python
# NVT (constant temperature)
final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,   # Kelvin
    timestep=1.0,        # fs
    steps=1000,
    trajectory="nvt.traj"  # Save trajectory to file
)

# NPT (constant pressure)
final = calc.run_md(
    atoms,
    ensemble="npt",
    temperature_K=300,
    pressure_GPa=0.0001, # GPa (~1 atm)
    timestep=1.0,
    steps=1000
)

# final is the last Atoms object
print(f"Final energy: {final.get_potential_energy():.4f} eV")
```

### Phonon Calculation

Calculate phonon properties:

```python
result = calc.phonon(
    atoms,
    supercell_matrix=[2, 2, 2],
    displacement=0.01,              # Å
    temperature_range=(0, 500, 10)  # (min, max, step) K
)

# Access results
frequencies = result['frequencies']  # THz at Gamma
phonon_obj = result['phonon']        # Phonopy object

# Thermal properties (if temperature_range provided)
if 'thermal' in result:
    print(f"Free energy: {result['thermal']['free_energy']} eV")
```

## D3 Dispersion Correction

Enable DFT-D3 for van der Waals interactions:

```python
# Enable D3 correction
calc = MACEInference(
    model="medium",
    device="auto",
    enable_d3=True,
    d3_xc="pbe",         # Exchange-correlation functional
    d3_damping="bj"      # Damping scheme
)

result = calc.single_point(atoms)
# Energy now includes D3 dispersion
```

!!! warning "Check your MACE model"
    Make sure your MACE model wasn't already trained with D3 correction, 
    or you'll double-count dispersion interactions.

## Command Line Interface

MACE Inference also provides a CLI:

```bash
# Single-point calculation
mace-inference single-point structure.cif --model medium

# Optimization
mace-inference optimize structure.cif --fmax 0.01 --output optimized.xyz

# Molecular dynamics
mace-inference md structure.cif --ensemble nvt --temperature 300 --steps 1000

# Phonon calculation
mace-inference phonon structure.cif --supercell 2 2 2
```

## Next Steps

- [User Guide](../user-guide/overview.md) - Detailed tutorials for each task
- [API Reference](../api/core.md) - Complete API documentation
- [Examples](../examples/index.md) - Real-world usage examples
