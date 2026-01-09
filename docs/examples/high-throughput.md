# Example 5: High-Throughput Screening

This example demonstrates screening multiple structures efficiently.

## Source Code

See [`examples/05_high_throughput.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/05_high_throughput.py)

## What You'll Learn

- Processing multiple structures
- Ranking by computed properties
- Parallel-like workflows
- Result aggregation and export

## Code Walkthrough

### Setup

```python
from pathlib import Path
from mace_inference import MACEInference
from ase.io import read
import json

calc = MACEInference(model="medium", device="auto")
```

### Load Multiple Structures

```python
structures_dir = Path("structures")
structures = []

for cif_path in structures_dir.glob("*.cif"):
    atoms = read(str(cif_path))
    atoms.info['source'] = cif_path.name
    structures.append(atoms)

print(f"Loaded {len(structures)} structures")
```

### Screen for Stability

```python
results = []

for atoms in structures:
    # Single-point for quick screening
    sp = calc.single_point(atoms)
    
    # Quick optimization to check stability
    opt = calc.optimize(atoms, fmax=0.1, steps=50)
    
    results.append({
        'name': atoms.info['source'],
        'formula': atoms.get_chemical_formula(),
        'E_per_atom': sp['energy_per_atom'],
        'max_force': sp['max_force'],
        'opt_converged': opt['converged'],
        'opt_steps': opt['steps']
    })
```

### Rank Results

```python
# Sort by energy per atom
ranked = sorted(results, key=lambda x: x['E_per_atom'])

print("\nRanking by stability (E/atom):")
print("-" * 50)
for i, r in enumerate(ranked, 1):
    status = "✓" if r['opt_converged'] else "○"
    print(f"{i}. {r['name']:<20} {r['E_per_atom']:.4f} eV/atom {status}")
```

### Export Results

```python
# Save to JSON
with open("screening_results.json", 'w') as f:
    json.dump(results, f, indent=2)

# Or CSV
with open("screening_results.csv", 'w') as f:
    f.write("name,formula,E_per_atom,converged\n")
    for r in results:
        f.write(f"{r['name']},{r['formula']},{r['E_per_atom']:.6f},{r['opt_converged']}\n")
```

## Expected Output

```
============================================================
Example 5: High-Throughput Screening
============================================================

1. Loading structures...
   ✓ cu_fcc.cif: Cu4
   ✓ cu_paddlewheel.cif: C8H4Cu2O8
   ✓ si_diamond.cif: Si8
   Total: 3 structures

2. Screening...
   [1/3] cu_fcc.cif... E/atom = -4.08 eV ✓
   [2/3] cu_paddlewheel.cif... E/atom = 12.89 eV ✓
   [3/3] si_diamond.cif... E/atom = -5.34 eV ✓

3. Ranking by stability:
   --------------------------------------------------
   1. si_diamond.cif        -5.3412 eV/atom ✓
   2. cu_fcc.cif            -4.0844 eV/atom ✓
   3. cu_paddlewheel.cif    12.8940 eV/atom ✓

4. Results saved:
   ✓ screening_results.json
   ✓ screening_results.csv

✅ Example 5 completed successfully!
```

## Key Points

1. **Reuse Calculator**: Initialize once, use for all structures
2. **Quick Screening**: Use loose tolerances first, then refine
3. **Error Handling**: Catch failures to continue processing
4. **Progress Tracking**: Print status for long-running jobs
