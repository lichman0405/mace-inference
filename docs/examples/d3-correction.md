# Example 6: D3 Dispersion Correction

This example demonstrates using DFT-D3 dispersion correction for van der Waals systems.

## Source Code

See [`examples/06_d3_correction.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/06_d3_correction.py)

## What You'll Learn

- Enabling D3 correction
- Comparing results with/without D3
- When D3 matters most
- D3 impact on structures

## Prerequisites

Install torch-dftd:

```bash
pip install torch-dftd
```

## Code Walkthrough

### Setup Two Calculators

```python
from mace_inference import MACEInference
from ase.io import read

# Without D3
calc_no_d3 = MACEInference(
    model="medium",
    device="auto",
    enable_d3=False
)

# With D3 (PBE functional, BJ damping)
calc_d3 = MACEInference(
    model="medium",
    device="auto",
    enable_d3=True,
    d3_xc="pbe",
    d3_damping="bj"
)
```

### Compare Energies

```python
atoms = read("structures/cu_paddlewheel.cif")

result_no_d3 = calc_no_d3.single_point(atoms)
result_d3 = calc_d3.single_point(atoms)

d3_energy = result_d3['energy'] - result_no_d3['energy']

print(f"MACE energy:      {result_no_d3['energy']:.4f} eV")
print(f"MACE + D3 energy: {result_d3['energy']:.4f} eV")
print(f"D3 contribution:  {d3_energy:.4f} eV")
print(f"D3 per atom:      {d3_energy/len(atoms)*1000:.1f} meV/atom")
```

### Compare Forces

```python
import numpy as np

f1 = result_no_d3['forces']
f2 = result_d3['forces']
f_diff = f2 - f1

print(f"Max MACE force:     {np.abs(f1).max():.4f} eV/Å")
print(f"Max MACE+D3 force:  {np.abs(f2).max():.4f} eV/Å")
print(f"Max D3 force:       {np.abs(f_diff).max():.4f} eV/Å")
```

### Compare Optimized Volumes

```python
opt_no_d3 = calc_no_d3.optimize(atoms, fmax=0.01, relax_cell=True)
opt_d3 = calc_d3.optimize(atoms, fmax=0.01, relax_cell=True)

V1 = opt_no_d3['atoms'].get_volume()
V2 = opt_d3['atoms'].get_volume()

print(f"Volume (no D3):   {V1:.2f} Å³")
print(f"Volume (with D3): {V2:.2f} Å³")
print(f"Difference:       {(V2/V1 - 1)*100:.1f}%")
```

## Expected Output

```
============================================================
Example 6: D3 Dispersion Correction
============================================================

1. Loading MOF structure...
   Formula: C8H4Cu2O8
   Atoms: 22

2. Energy comparison...
   MACE energy:      283.6685 eV
   MACE + D3 energy: 282.5049 eV
   D3 contribution:  -1.1636 eV
   D3 per atom:      -52.9 meV/atom

3. Force comparison...
   Max MACE force:     1221.13 eV/Å
   Max MACE+D3 force:  1221.18 eV/Å
   Max D3 force:       0.06 eV/Å

4. Volume comparison...
   Volume (no D3):   1250.00 Å³
   Volume (with D3): 1250.00 Å³
   Difference:       +0.0%

✅ Example 6 completed successfully!
```

## When D3 Matters

| System | D3 Effect | Recommendation |
|--------|-----------|----------------|
| MOFs | Moderate | Use D3 |
| Layered materials | Large | Use D3 |
| Molecular crystals | Large | Use D3 |
| Adsorption | Moderate-Large | Use D3 |
| Bulk metals | Small | Usually skip |
| Covalent solids | Small | Usually skip |

## Key Points

1. **Always Negative**: D3 adds attractive dispersion, lowering energy
2. **Volume Contraction**: D3 typically causes slight volume decrease
3. **Small Force Effect**: D3 forces are usually small compared to MACE
4. **No Double Counting**: Don't use D3 if MACE was trained with D3
