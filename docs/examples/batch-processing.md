# Example 7: Batch Processing

This example demonstrates processing multiple structures with comprehensive I/O handling.

## Source Code

See [`examples/07_batch_processing.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/07_batch_processing.py)

## What You'll Learn

- Loading multiple structures from a directory
- Batch single-point calculations
- Batch optimizations
- Exporting results in multiple formats
- Error handling for robustness

## Code Walkthrough

### Setup

```python
from pathlib import Path
from mace_inference import MACEInference
from ase.io import read, write
import json

calc = MACEInference(model="medium", device="auto")

STRUCTURES_DIR = Path("structures")
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)
```

### Load Structures

```python
structures = []
for cif_path in STRUCTURES_DIR.glob("*.cif"):
    atoms = read(str(cif_path))
    atoms.info['source_file'] = cif_path.name
    structures.append(atoms)
    print(f"   ✓ {cif_path.name}: {atoms.get_chemical_formula()}")
```

### Batch Single-Point

```python
results = []

for i, atoms in enumerate(structures, 1):
    source = atoms.info.get('source_file', f'structure_{i}')
    print(f"\n   [{i}/{len(structures)}] Processing {source}...")
    
    try:
        result = calc.single_point(atoms)
        atoms.info['mace_energy'] = result['energy']
        
        results.append({
            'source': source,
            'formula': atoms.get_chemical_formula(),
            'n_atoms': len(atoms),
            'energy': result['energy'],
            'energy_per_atom': result['energy'] / len(atoms),
            'max_force': float(result['max_force']),
            'status': 'success'
        })
        print(f"       Energy: {result['energy']:.4f} eV ✓")
        
    except Exception as e:
        results.append({
            'source': source,
            'status': 'failed',
            'error': str(e)
        })
        print(f"       ✗ Failed: {e}")
```

### Export Results

```python
# JSON
with open(OUTPUT_DIR / "batch_results.json", 'w') as f:
    json.dump(results, f, indent=2)

# CSV
with open(OUTPUT_DIR / "batch_results.csv", 'w') as f:
    f.write("source,formula,n_atoms,energy,energy_per_atom,max_force,status\n")
    for r in results:
        if r['status'] == 'success':
            f.write(f"{r['source']},{r['formula']},{r['n_atoms']},"
                    f"{r['energy']:.6f},{r['energy_per_atom']:.6f},"
                    f"{r['max_force']:.6f},{r['status']}\n")
```

### Save Annotated Structures

```python
for atoms in structures:
    source = atoms.info.get('source_file', 'unknown')
    if 'mace_energy' in atoms.info:
        base_name = Path(source).stem
        xyz_path = OUTPUT_DIR / f"{base_name}_annotated.xyz"
        atoms_copy = atoms.copy()
        atoms_copy.info = atoms.info.copy()
        write(str(xyz_path), atoms_copy, format='extxyz')
```

### Batch Optimization

```python
for i, atoms in enumerate(structures, 1):
    source = atoms.info.get('source_file', f'structure_{i}')
    print(f"\n   [{i}/{len(structures)}] Optimizing {source}...")
    
    initial_E = calc.single_point(atoms)['energy']
    opt = calc.optimize(atoms, fmax=0.05, steps=200)
    
    delta_E = opt['energy'] - initial_E
    print(f"       ΔE: {delta_E:.4f} eV")
    
    # Save optimized structure
    base_name = Path(source).stem
    write(str(OUTPUT_DIR / f"{base_name}_optimized.xyz"), opt['atoms'])
```

## Expected Output

```
============================================================
Example 7: Batch Processing and Utilities
============================================================

1. Initializing MACE calculator...

2. Loading structures...
   ✓ cu_fcc.cif: Cu4
   ✓ cu_paddlewheel.cif: C8H4Cu2O8
   ✓ si_diamond.cif: Si8
   Total: 3 structures

3. Batch Single-Point Calculations...
   [1/3] Processing cu_fcc.cif...
       Energy: -16.3376 eV ✓
   [2/3] Processing cu_paddlewheel.cif...
       Energy: 283.6685 eV ✓
   [3/3] Processing si_diamond.cif...
       Energy: -42.7293 eV ✓

4. Exporting Results...
   ✓ JSON results: batch_results.json
   ✓ CSV results: batch_results.csv

5. Saving Annotated Structures...
   ✓ cu_fcc_annotated.xyz
   ✓ cu_paddlewheel_annotated.xyz
   ✓ si_diamond_annotated.xyz

6. Batch Optimization...
   [1/3] Optimizing cu_fcc.cif...
       ΔE: 0.0000 eV
   [2/3] Optimizing cu_paddlewheel.cif...
       ΔE: -414.7552 eV
   [3/3] Optimizing si_diamond.cif...
       ΔE: 0.0000 eV

7. Summary Statistics...
   Total structures: 3
   Successful: 3
   Failed: 0

✅ Example 7 completed successfully!
```

## Output Files

After running, you'll have:

```
output/
├── batch_results.json      # Full results in JSON
├── batch_results.csv       # Summary table
├── cu_fcc_annotated.xyz    # With energy in info
├── cu_fcc_optimized.xyz    # Optimized structure
├── cu_paddlewheel_annotated.xyz
├── cu_paddlewheel_optimized.xyz
├── si_diamond_annotated.xyz
└── si_diamond_optimized.xyz
```

## Key Points

1. **Error Handling**: Use try/except to continue on failures
2. **Progress Tracking**: Print status for monitoring
3. **Multiple Formats**: Export both JSON (detailed) and CSV (analysis)
4. **Copy for Write**: Use `.copy()` to avoid issues with cached arrays
