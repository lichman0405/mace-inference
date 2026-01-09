#!/usr/bin/env python
"""
Example 7: Batch Processing and Utilities

This example demonstrates:
- Batch single-point calculations
- File I/O utilities
- Error handling and progress tracking
"""

from pathlib import Path
from ase.io import read, write
from ase import Atoms
from mace_inference import MACEInference
import json

EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"
OUTPUT_DIR = EXAMPLES_DIR / "output"


def main():
    print("=" * 60)
    print("Example 7: Batch Processing and Utilities")
    print("=" * 60)
    
    # Create output directory
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Initialize calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="medium", device="auto")
    
    # Load multiple structures
    print("\n2. Loading structures...")
    structures = []
    
    for cif_path in STRUCTURES_DIR.glob("*.cif"):
        atoms = read(str(cif_path))
        atoms.info['source_file'] = cif_path.name
        structures.append(atoms)
        print(f"   ✓ {cif_path.name}: {atoms.get_chemical_formula()}")
    
    print(f"\n   Total: {len(structures)} structures loaded")
    
    # Batch single-point calculation
    print("\n" + "=" * 60)
    print("3. Batch Single-Point Calculations")
    print("=" * 60)
    
    results = []
    
    for i, atoms in enumerate(structures, 1):
        source = atoms.info.get('source_file', f'structure_{i}')
        print(f"\n   [{i}/{len(structures)}] Processing {source}...")
        
        try:
            # Calculate
            result = calc.single_point(atoms)
            
            # Store results in atoms.info
            atoms.info['mace_energy'] = result['energy']
            atoms.info['mace_energy_per_atom'] = result['energy'] / len(atoms)
            atoms.info['mace_max_force'] = float(abs(result['forces']).max())
            
            results.append({
                'source': source,
                'formula': atoms.get_chemical_formula(),
                'n_atoms': len(atoms),
                'energy': result['energy'],
                'energy_per_atom': result['energy'] / len(atoms),
                'max_force': float(abs(result['forces']).max()),
                'status': 'success'
            })
            
            print(f"       Energy: {result['energy']:.4f} eV")
            print(f"       Status: ✓ Success")
            
        except Exception as e:
            results.append({
                'source': source,
                'formula': atoms.get_chemical_formula(),
                'status': 'failed',
                'error': str(e)
            })
            print(f"       Status: ✗ Failed - {e}")
    
    # Export results
    print("\n" + "=" * 60)
    print("4. Exporting Results")
    print("=" * 60)
    
    # Save as JSON
    json_path = OUTPUT_DIR / "batch_results.json"
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"   ✓ JSON results: {json_path.name}")
    
    # Save as CSV
    csv_path = OUTPUT_DIR / "batch_results.csv"
    with open(csv_path, 'w') as f:
        f.write("source,formula,n_atoms,energy,energy_per_atom,max_force,status\n")
        for r in results:
            if r['status'] == 'success':
                f.write(f"{r['source']},{r['formula']},{r['n_atoms']},{r['energy']:.6f},{r['energy_per_atom']:.6f},{r['max_force']:.6f},{r['status']}\n")
            else:
                f.write(f"{r['source']},{r['formula']},,,,,{r['status']}\n")
    print(f"   ✓ CSV results: {csv_path.name}")
    
    # Save structures with energies to XYZ
    print("\n" + "=" * 60)
    print("5. Saving Annotated Structures")
    print("=" * 60)
    
    for atoms in structures:
        source = atoms.info.get('source_file', 'unknown')
        if 'mace_energy' in atoms.info:
            base_name = Path(source).stem
            xyz_path = OUTPUT_DIR / f"{base_name}_annotated.xyz"
            # Copy atoms to avoid issues with cached calculator arrays
            atoms_copy = atoms.copy()
            atoms_copy.info = atoms.info.copy()
            write(str(xyz_path), atoms_copy, format='extxyz')
            print(f"   ✓ {xyz_path.name}")
    
    # Batch optimization
    print("\n" + "=" * 60)
    print("6. Batch Optimization")
    print("=" * 60)
    
    opt_results = []
    
    for i, atoms in enumerate(structures, 1):
        source = atoms.info.get('source_file', f'structure_{i}')
        print(f"\n   [{i}/{len(structures)}] Optimizing {source}...")
        
        try:
            # Copy to preserve original
            atoms_opt = atoms.copy()
            
            # Optimize
            optimized = calc.optimize(atoms_opt, fmax=0.05, steps=30)
            
            # Calculate energy change
            initial_e = atoms.info.get('mace_energy', 0)
            final_e = calc.single_point(optimized)['energy']
            delta_e = final_e - initial_e
            
            opt_results.append({
                'source': source,
                'initial_energy': initial_e,
                'final_energy': final_e,
                'delta_e': delta_e,
                'status': 'success'
            })
            
            # Save optimized structure
            base_name = Path(source).stem
            opt_path = OUTPUT_DIR / f"{base_name}_optimized.xyz"
            write(str(opt_path), optimized, format='extxyz')
            
            print(f"       ΔE: {delta_e:.4f} eV")
            print(f"       Saved: {opt_path.name}")
            
        except Exception as e:
            print(f"       Failed: {e}")
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("7. Summary Statistics")
    print("=" * 60)
    
    successful = [r for r in results if r['status'] == 'success']
    failed = [r for r in results if r['status'] == 'failed']
    
    print(f"\n   Total structures: {len(results)}")
    print(f"   Successful: {len(successful)}")
    print(f"   Failed: {len(failed)}")
    
    if successful:
        energies = [r['energy_per_atom'] for r in successful]
        print(f"\n   Energy per atom statistics:")
        print(f"     Min: {min(energies):.4f} eV")
        print(f"     Max: {max(energies):.4f} eV")
        print(f"     Mean: {sum(energies)/len(energies):.4f} eV")
    
    print(f"\n   Output files saved to: {OUTPUT_DIR}")
    print("\n✅ Example 7 completed successfully!")


if __name__ == "__main__":
    main()
