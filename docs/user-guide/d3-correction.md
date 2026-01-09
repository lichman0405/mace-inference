# D3 Dispersion Correction

DFT-D3 is a widely-used method to add van der Waals (dispersion) interactions to DFT calculations. MACE Inference supports D3 correction through the `torch-dftd` library.

## When to Use D3

D3 correction is recommended for systems where dispersion interactions are important:

| System Type | D3 Recommended |
|-------------|----------------|
| Van der Waals crystals | ✅ Yes |
| Metal-organic frameworks (MOFs) | ✅ Yes |
| Layered materials (graphite, MoS₂) | ✅ Yes |
| Molecular crystals | ✅ Yes |
| Adsorption on surfaces | ✅ Yes |
| Polymer chains | ✅ Yes |
| Covalent/ionic crystals (Si, NaCl) | ❌ Usually not needed |
| Pure metals | ❌ Usually not needed |

!!! warning "Check your MACE model"
    If your MACE model was trained on DFT data that already includes D3 
    correction, adding D3 again will **double-count** dispersion interactions.
    The pretrained Materials Project models do NOT include D3.

## Installation

D3 correction requires the `torch-dftd` package:

```bash
pip install torch-dftd
```

## Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read

# Enable D3 correction
calc = MACEInference(
    model="medium",
    device="auto",
    enable_d3=True,
    d3_xc="pbe"  # Exchange-correlation functional
)

atoms = read("mof.cif")
result = calc.single_point(atoms)
# Energy now includes D3 dispersion contribution
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `enable_d3` | `bool` | `False` | Enable D3 correction |
| `d3_xc` | `str` | `"pbe"` | Exchange-correlation functional |
| `d3_damping` | `str` | `"bj"` | Damping function |

### Supported Functionals

Common functionals for `d3_xc`:

- `"pbe"` - PBE (most common)
- `"rpbe"` - Revised PBE
- `"blyp"` - BLYP
- `"b3lyp"` - B3LYP
- `"tpss"` - TPSS
- `"scan"` - SCAN

### Damping Functions

- `"bj"` - Becke-Johnson damping (recommended, default)
- `"zero"` - Zero damping (original D3)

## Examples

### Comparing With and Without D3

```python
# Without D3
calc_no_d3 = MACEInference(model="medium", enable_d3=False)

# With D3
calc_d3 = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")

atoms = read("mof.cif")

result_no_d3 = calc_no_d3.single_point(atoms)
result_d3 = calc_d3.single_point(atoms)

d3_contribution = result_d3['energy'] - result_no_d3['energy']

print(f"MACE energy:      {result_no_d3['energy']:.4f} eV")
print(f"MACE + D3 energy: {result_d3['energy']:.4f} eV")
print(f"D3 contribution:  {d3_contribution:.4f} eV")
print(f"D3 per atom:      {d3_contribution / len(atoms) * 1000:.1f} meV/atom")
```

### Volume Optimization with D3

D3 affects equilibrium volume, especially for vdW materials:

```python
from ase.io import read

atoms = read("layered_material.cif")

# Optimize without D3
calc1 = MACEInference(model="medium", enable_d3=False)
opt1 = calc1.optimize(atoms, fmax=0.01, relax_cell=True)

# Optimize with D3
calc2 = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")
opt2 = calc2.optimize(atoms, fmax=0.01, relax_cell=True)

V1 = opt1['atoms'].get_volume()
V2 = opt2['atoms'].get_volume()

print(f"Volume without D3: {V1:.2f} Å³")
print(f"Volume with D3:    {V2:.2f} Å³")
print(f"Difference:        {(V2/V1 - 1)*100:.1f}%")
```

Typically, D3 causes contraction because it adds attractive interactions.

### Adsorption with D3

D3 is important for accurate adsorption energies:

```python
from ase.build import molecule

framework = read("HKUST-1.cif")
co2 = molecule("CO2")

calc = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")

result = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=[10.0, 10.0, 10.0]
)

print(f"Adsorption energy (with D3): {result['E_ads']:.4f} eV")
```

### Different Damping Schemes

```python
# Becke-Johnson damping (default, recommended)
calc_bj = MACEInference(
    model="medium",
    enable_d3=True,
    d3_xc="pbe",
    d3_damping="bj"
)

# Zero damping (original)
calc_zero = MACEInference(
    model="medium",
    enable_d3=True,
    d3_xc="pbe",
    d3_damping="zero"
)

result_bj = calc_bj.single_point(atoms)
result_zero = calc_zero.single_point(atoms)

print(f"D3-BJ energy: {result_bj['energy']:.4f} eV")
print(f"D3-zero energy: {result_zero['energy']:.4f} eV")
```

## Force Contributions

D3 also contributes to forces:

```python
import numpy as np

calc_no_d3 = MACEInference(model="medium", enable_d3=False)
calc_d3 = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")

result_no_d3 = calc_no_d3.single_point(atoms)
result_d3 = calc_d3.single_point(atoms)

force_diff = result_d3['forces'] - result_no_d3['forces']
max_d3_force = np.abs(force_diff).max()

print(f"Max D3 force contribution: {max_d3_force:.4f} eV/Å")
```

## Performance

D3 calculation adds minimal overhead:

```python
import time

atoms = read("large_mof.cif")  # Large system

# Without D3
calc1 = MACEInference(model="medium", enable_d3=False)
t1 = time.time()
calc1.single_point(atoms)
t1 = time.time() - t1

# With D3
calc2 = MACEInference(model="medium", enable_d3=True)
t2 = time.time()
calc2.single_point(atoms)
t2 = time.time() - t2

print(f"Time without D3: {t1:.2f} s")
print(f"Time with D3:    {t2:.2f} s")
print(f"Overhead:        {(t2/t1 - 1)*100:.1f}%")
```

## Common Issues

### torch-dftd Not Installed

```
ImportError: D3 correction requires torch-dftd. Install with: pip install torch-dftd
```

**Solution**: Install the package:
```bash
pip install torch-dftd
```

### Unsupported Functional

```
ValueError: Unknown XC functional 'xyz'
```

**Solution**: Use a supported functional like `"pbe"`, `"rpbe"`, etc.

### Double Counting

If your results seem too negative (over-binding):

1. Check if your MACE model was trained with D3
2. Try without D3 and compare to reference data
3. Check the functional matches the training data
