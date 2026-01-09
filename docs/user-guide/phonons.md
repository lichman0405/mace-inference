# Phonon Calculations

Phonon calculations analyze the vibrational properties of crystals, providing information about dynamical stability, thermal properties, and heat capacity.

## Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read

calc = MACEInference(model="medium", device="auto")
atoms = read("structure.cif")

result = calc.phonon(atoms, supercell=(2, 2, 2))
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `atoms` | `Atoms` | required | Unit cell structure |
| `supercell` | `tuple` | `(2, 2, 2)` | Supercell dimensions |
| `delta` | `float` | `0.01` | Displacement magnitude (Å) |
| `mesh` | `tuple` | `(20, 20, 20)` | q-point mesh for DOS |
| `temperature` | `float` | `300` | Temperature for thermal properties (K) |

## Return Values

| Key | Type | Description |
|-----|------|-------------|
| `frequencies` | `ndarray` | Phonon frequencies at Γ (THz) |
| `band_structure` | `dict` | Band structure data |
| `dos` | `dict` | Phonon density of states |
| `thermal_properties` | `dict` | Free energy, entropy, Cv |

### Thermal Properties

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `free_energy` | `float` | eV | Helmholtz free energy |
| `entropy` | `float` | eV/K | Vibrational entropy |
| `heat_capacity` | `float` | eV/K | Heat capacity at constant V |

## Examples

### Basic Phonon Calculation

```python
result = calc.phonon(
    atoms,
    supercell=(2, 2, 2),
    delta=0.01
)

# Frequencies at Gamma point
freqs = result['frequencies']
print(f"Number of modes: {len(freqs)}")
print(f"Lowest frequency: {freqs.min():.4f} THz")
print(f"Highest frequency: {freqs.max():.4f} THz")
```

### Check Dynamical Stability

Imaginary frequencies indicate instability:

```python
result = calc.phonon(atoms, supercell=(2, 2, 2))

freqs = result['frequencies']
imaginary = freqs[freqs < 0]

if len(imaginary) > 0:
    print("⚠️ Structure is dynamically UNSTABLE")
    print(f"Imaginary modes: {imaginary}")
else:
    print("✓ Structure is dynamically stable")
```

!!! note "Imaginary frequencies"
    Small negative frequencies near zero (< 0.1 THz) are often numerical 
    artifacts from acoustic modes and can usually be ignored.

### Thermal Properties at Different Temperatures

```python
# Calculate at 300 K (default)
result = calc.phonon(atoms, supercell=(2, 2, 2), temperature=300)

thermal = result['thermal_properties']
print(f"At 300 K:")
print(f"  Free energy: {thermal['free_energy']:.4f} eV")
print(f"  Entropy: {thermal['entropy']:.6f} eV/K")
print(f"  Heat capacity: {thermal['heat_capacity']:.6f} eV/K")
```

### Phonon DOS Analysis

```python
result = calc.phonon(atoms, supercell=(2, 2, 2), mesh=(30, 30, 30))

dos = result['dos']
frequencies = dos['frequencies']
density = dos['total_dos']

# Find peak
peak_idx = density.argmax()
print(f"DOS peak at {frequencies[peak_idx]:.2f} THz")
```

### Plotting Phonon Band Structure

```python
import matplotlib.pyplot as plt

result = calc.phonon(atoms, supercell=(3, 3, 3))
band = result['band_structure']

plt.figure(figsize=(10, 6))
for i in range(len(band['frequencies'][0])):
    plt.plot(band['distances'], [f[i] for f in band['frequencies']], 'b-', lw=0.5)

plt.axhline(0, color='k', linestyle='--', lw=0.5)
plt.xlabel("Wave vector")
plt.ylabel("Frequency (THz)")
plt.title("Phonon Band Structure")
plt.savefig("phonon_bands.png", dpi=150)
```

### Plotting Phonon DOS

```python
import matplotlib.pyplot as plt

result = calc.phonon(atoms, supercell=(3, 3, 3))
dos = result['dos']

plt.figure(figsize=(8, 5))
plt.fill_between(dos['frequencies'], dos['total_dos'], alpha=0.5)
plt.plot(dos['frequencies'], dos['total_dos'])
plt.xlabel("Frequency (THz)")
plt.ylabel("DOS")
plt.title("Phonon Density of States")
plt.axvline(0, color='k', linestyle='--', lw=0.5)
plt.savefig("phonon_dos.png", dpi=150)
```

## Supercell Considerations

### Choosing Supercell Size

The supercell must be large enough to capture phonon interactions:

| Material Type | Recommended Supercell |
|--------------|----------------------|
| Simple metals (Cu, Al) | (2, 2, 2) - (3, 3, 3) |
| Semiconductors (Si, GaAs) | (3, 3, 3) - (4, 4, 4) |
| Complex oxides | (2, 2, 2) - (3, 3, 3) |
| Layered materials | Larger in stacking direction |

```python
# For a material with small unit cell
result = calc.phonon(atoms, supercell=(4, 4, 4))

# For a material with large unit cell
result = calc.phonon(atoms, supercell=(2, 2, 2))
```

### Anisotropic Supercells

For non-cubic or layered materials:

```python
# Layered material (e.g., graphite)
result = calc.phonon(atoms, supercell=(4, 4, 2))

# 1D chain
result = calc.phonon(atoms, supercell=(1, 1, 6))
```

## Advanced Options

### Custom q-point Mesh

```python
# Finer mesh for accurate DOS
result = calc.phonon(
    atoms,
    supercell=(2, 2, 2),
    mesh=(40, 40, 40)  # Finer q-mesh
)
```

### Displacement Magnitude

```python
# Larger displacement for better signal
result = calc.phonon(atoms, supercell=(2, 2, 2), delta=0.02)

# Smaller displacement for more linear regime
result = calc.phonon(atoms, supercell=(2, 2, 2), delta=0.005)
```

## Workflow Tips

### Always Optimize First

Phonon calculations require a well-relaxed structure:

```python
# Step 1: Optimize with tight tolerance
opt = calc.optimize(atoms, fmax=0.001)

# Step 2: Calculate phonons
phonon = calc.phonon(opt['atoms'], supercell=(3, 3, 3))
```

### Convergence Testing

Test convergence with supercell size:

```python
for size in [(2, 2, 2), (3, 3, 3), (4, 4, 4)]:
    result = calc.phonon(atoms, supercell=size)
    freqs = result['frequencies']
    print(f"Supercell {size}: max freq = {freqs.max():.2f} THz")
```

## Common Issues

### Imaginary Frequencies at Γ

1. **Structure not relaxed**: Optimize with tighter tolerance
2. **Wrong space group**: Check if symmetry is correct
3. **Numerical noise**: Try smaller `delta`

### Calculation Too Slow

1. Reduce supercell size
2. Use GPU acceleration
3. Use smaller q-mesh for DOS

### Memory Issues

Large supercells require significant memory:

```python
# Reduce mesh for initial testing
result = calc.phonon(atoms, supercell=(2, 2, 2), mesh=(10, 10, 10))
```
