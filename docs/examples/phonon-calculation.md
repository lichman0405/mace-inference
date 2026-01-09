# Example 3: Phonon Calculation

This example demonstrates phonon calculations and thermal property analysis.

## Source Code

See [`examples/03_phonon_calculation.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/03_phonon_calculation.py)

## What You'll Learn

- Calculating phonon frequencies
- Checking dynamical stability
- Computing thermal properties
- Plotting phonon DOS

## Code Walkthrough

### Setup

```python
from mace_inference import MACEInference
from ase.io import read

calc = MACEInference(model="medium", device="auto")
atoms = read("structures/si_diamond.cif")
```

### Optimize Unit Cell

Phonon calculations require well-relaxed structures:

```python
opt = calc.optimize(atoms, fmax=0.001)
atoms = opt['atoms']
```

### Calculate Phonons

```python
result = calc.phonon(
    atoms,
    supercell=(2, 2, 2),  # Supercell for force constants
    delta=0.01,           # Displacement in Å
    temperature=300       # For thermal properties
)
```

### Check Stability

```python
frequencies = result['frequencies']
imaginary = frequencies[frequencies < -0.1]  # Threshold for numerical noise

if len(imaginary) > 0:
    print("⚠️ Structure is dynamically unstable!")
    print(f"Imaginary modes: {imaginary}")
else:
    print("✓ Structure is dynamically stable")
```

### Thermal Properties

```python
thermal = result['thermal_properties']

print(f"At 300 K:")
print(f"  Free energy: {thermal['free_energy']:.4f} eV")
print(f"  Entropy: {thermal['entropy']*1000:.4f} meV/K")
print(f"  Heat capacity: {thermal['heat_capacity']*1000:.4f} meV/K")
```

### Plot DOS

```python
import matplotlib.pyplot as plt

dos = result['dos']
plt.figure(figsize=(8, 5))
plt.fill_between(dos['frequencies'], dos['total_dos'], alpha=0.5)
plt.plot(dos['frequencies'], dos['total_dos'])
plt.xlabel("Frequency (THz)")
plt.ylabel("DOS (states/THz)")
plt.title("Phonon Density of States")
plt.axvline(0, color='k', linestyle='--', lw=0.5)
plt.savefig("phonon_dos.png", dpi=150)
```

## Expected Output

```
============================================================
Example 3: Phonon Calculation
============================================================

1. Loading structure: Si8

2. Optimizing structure...
   Converged in 1 steps

3. Calculating phonons...
   Supercell: (2, 2, 2) = 64 atoms
   Displacements: 6
   ✓ Force constants calculated

4. Analyzing results...
   Gamma-point frequencies (THz):
     0.00, 0.00, 0.00 (acoustic)
     12.45, 12.45, 12.45
     15.23, 15.23, 15.23
     ...
   
   ✓ No imaginary frequencies - structure is stable

5. Thermal properties at 300 K:
   Free energy: 0.2341 eV
   Entropy: 0.4521 meV/K
   Heat capacity: 0.3892 meV/K

✅ Example 3 completed successfully!
```

## Key Points

1. **Supercell Size**: Larger supercells give more accurate results but are slower
2. **Optimization**: Tight force tolerance (< 0.001 eV/Å) is important
3. **Imaginary Modes**: Small negatives near zero are numerical noise; large negatives indicate instability
4. **Units**: Frequencies in THz, thermal properties per atom
