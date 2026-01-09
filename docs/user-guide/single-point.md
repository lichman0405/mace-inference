# Single Point Calculations

Single-point calculations compute the energy, forces, and stress of a structure at its current geometry without any optimization.

## Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read

calc = MACEInference(model="medium", device="auto")
atoms = read("structure.cif")

result = calc.single_point(atoms)
```

## Return Values

The `single_point()` method returns a dictionary:

| Key | Type | Description |
|-----|------|-------------|
| `energy` | `float` | Total potential energy (eV) |
| `energy_per_atom` | `float` | Energy divided by number of atoms (eV) |
| `forces` | `ndarray` | Forces on each atom, shape (N, 3) (eV/Å) |
| `max_force` | `float` | Maximum force magnitude (eV/Å) |
| `stress` | `ndarray` | Stress tensor in Voigt notation (eV/Å³) |

## Examples

### Basic Energy Calculation

```python
result = calc.single_point(atoms)

print(f"Total energy: {result['energy']:.6f} eV")
print(f"Energy per atom: {result['energy_per_atom']:.6f} eV")
print(f"Number of atoms: {len(atoms)}")
```

### Force Analysis

```python
import numpy as np

result = calc.single_point(atoms)
forces = result['forces']

print(f"Force shape: {forces.shape}")
print(f"Max force: {result['max_force']:.4f} eV/Å")
print(f"Mean force magnitude: {np.linalg.norm(forces, axis=1).mean():.4f} eV/Å")

# Check if structure is relaxed
if result['max_force'] < 0.05:
    print("Structure appears to be well-relaxed")
else:
    print("Structure may need optimization")
```

### Stress Analysis

```python
result = calc.single_point(atoms)
stress = result['stress']  # Voigt notation: xx, yy, zz, yz, xz, xy

print("Stress tensor (eV/Å³):")
print(f"  σxx = {stress[0]:.6f}")
print(f"  σyy = {stress[1]:.6f}")
print(f"  σzz = {stress[2]:.6f}")

# Calculate pressure (negative trace / 3)
pressure = -np.mean(stress[:3])
print(f"Pressure: {pressure:.6f} eV/Å³")
print(f"Pressure: {pressure * 160.2:.2f} GPa")  # Convert to GPa
```

### Comparing Structures

```python
from ase.io import read

structures = [read(f"structure_{i}.cif") for i in range(5)]

print("Structure comparison:")
print("-" * 50)
print(f"{'File':<20} {'E (eV)':<12} {'E/atom (eV)':<12}")
print("-" * 50)

for i, atoms in enumerate(structures):
    result = calc.single_point(atoms)
    print(f"structure_{i}.cif    {result['energy']:<12.4f} {result['energy_per_atom']:<12.4f}")
```

## With D3 Correction

For systems where dispersion interactions are important:

```python
# Without D3
calc_no_d3 = MACEInference(model="medium", enable_d3=False)
result_no_d3 = calc_no_d3.single_point(atoms)

# With D3
calc_d3 = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")
result_d3 = calc_d3.single_point(atoms)

d3_contribution = result_d3['energy'] - result_no_d3['energy']
print(f"MACE energy: {result_no_d3['energy']:.4f} eV")
print(f"MACE + D3 energy: {result_d3['energy']:.4f} eV")
print(f"D3 contribution: {d3_contribution:.4f} eV")
```

## Performance Tips

### Batch Calculations

For many structures, reuse the calculator:

```python
calc = MACEInference(model="medium", device="auto")

results = []
for atoms in structure_list:
    result = calc.single_point(atoms)
    results.append(result)
```

### GPU Acceleration

For large systems, GPU provides significant speedup:

```python
calc = MACEInference(model="medium", device="cuda")
```

### Lower Precision for Screening

For quick screening, use float32:

```python
calc = MACEInference(
    model="medium", 
    device="auto",
    default_dtype="float32"  # Faster but less accurate
)
```

## Common Use Cases

1. **Quick energy check** before expensive calculations
2. **Screening** many candidate structures
3. **Validation** of ML model predictions
4. **Force verification** to check if structure is relaxed
5. **Stress analysis** for pressure calculations
