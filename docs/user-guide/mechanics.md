# Mechanical Properties

MACE Inference can calculate elastic properties such as bulk modulus by fitting equations of state.

## Bulk Modulus

The bulk modulus measures a material's resistance to uniform compression.

### Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read

calc = MACEInference(model="medium", device="auto")
atoms = read("structure.cif")

result = calc.bulk_modulus(atoms)
print(f"Bulk modulus: {result['bulk_modulus']:.1f} GPa")
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `atoms` | `Atoms` | required | Structure (should be optimized) |
| `scale_range` | `tuple` | `(0.95, 1.05)` | Volume scaling range |
| `n_points` | `int` | `7` | Number of points to fit |
| `eos` | `str` | `"birchmurnaghan"` | Equation of state |

### Return Values

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `bulk_modulus` | `float` | GPa | Bulk modulus |
| `equilibrium_volume` | `float` | Å³ | Equilibrium volume |
| `equilibrium_energy` | `float` | eV | Energy at equilibrium |
| `volumes` | `list` | Å³ | Volumes used in fit |
| `energies` | `list` | eV | Energies at each volume |

### Available EOS

- `"birchmurnaghan"` - Birch-Murnaghan (default, most common)
- `"murnaghan"` - Murnaghan
- `"vinet"` - Vinet

## Examples

### Basic Bulk Modulus

```python
result = calc.bulk_modulus(atoms)

print(f"Bulk modulus: {result['bulk_modulus']:.1f} GPa")
print(f"Equilibrium volume: {result['equilibrium_volume']:.2f} Å³")
print(f"Equilibrium energy: {result['equilibrium_energy']:.6f} eV")
```

### Custom Range and Points

```python
# For softer materials, use wider range
result = calc.bulk_modulus(
    atoms,
    scale_range=(0.90, 1.10),  # ±10% volume
    n_points=11
)

# For harder materials, narrower range is fine
result = calc.bulk_modulus(
    atoms,
    scale_range=(0.97, 1.03),  # ±3% volume
    n_points=5
)
```

### Plotting E-V Curve

```python
import matplotlib.pyplot as plt
import numpy as np

result = calc.bulk_modulus(atoms, n_points=11)

plt.figure(figsize=(8, 6))
plt.plot(result['volumes'], result['energies'], 'o-', markersize=8)
plt.xlabel("Volume (Å³)")
plt.ylabel("Energy (eV)")
plt.title(f"E-V Curve (B = {result['bulk_modulus']:.1f} GPa)")
plt.axvline(result['equilibrium_volume'], color='r', linestyle='--', 
            label=f"V₀ = {result['equilibrium_volume']:.2f} Å³")
plt.legend()
plt.savefig("ev_curve.png", dpi=150)
```

### Comparing Materials

```python
materials = ["Cu.cif", "Al.cif", "Si.cif", "Fe.cif"]

print(f"{'Material':<15} {'B (GPa)':<12} {'V₀ (Å³)':<12}")
print("-" * 40)

for mat_file in materials:
    atoms = read(mat_file)
    # Optimize first
    opt = calc.optimize(atoms, fmax=0.01)
    # Calculate bulk modulus
    result = calc.bulk_modulus(opt['atoms'])
    print(f"{mat_file:<15} {result['bulk_modulus']:<12.1f} {result['equilibrium_volume']:<12.2f}")
```

### Different EOS

```python
eos_types = ["birchmurnaghan", "murnaghan", "vinet"]

for eos in eos_types:
    result = calc.bulk_modulus(atoms, eos=eos)
    print(f"{eos}: B = {result['bulk_modulus']:.1f} GPa")
```

## Workflow Tips

### Always Optimize First

```python
# Bulk modulus requires a well-relaxed structure
opt = calc.optimize(atoms, fmax=0.001, relax_cell=True)

# Then calculate bulk modulus
result = calc.bulk_modulus(opt['atoms'])
```

### Check Fit Quality

If the bulk modulus seems unreasonable:

1. Check if structure is properly relaxed
2. Try different volume range
3. Increase number of points
4. Plot E-V curve to visually inspect fit

```python
result = calc.bulk_modulus(atoms, n_points=15)

# Check if E-V curve is smooth
import numpy as np
energies = np.array(result['energies'])
if np.any(np.diff(energies[:len(energies)//2]) > 0):
    print("Warning: E-V curve may not be well-behaved at low volumes")
```

## Common Issues

### Negative or Zero Bulk Modulus

This indicates the structure is unstable:

1. Structure not at equilibrium - optimize first
2. Wrong crystal structure
3. Volume range too wide

### Very Large Bulk Modulus

1. Volume range too narrow
2. Not enough points for fit
3. Structure is under strain

### Inconsistent Results

Try different EOS and compare:

```python
for eos in ["birchmurnaghan", "murnaghan", "vinet"]:
    result = calc.bulk_modulus(atoms, eos=eos)
    print(f"{eos}: {result['bulk_modulus']:.1f} GPa")
```

If results differ significantly, the volume range may need adjustment.
