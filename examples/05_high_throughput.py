#!/usr/bin/env python
"""
Example 5: High-Throughput Screening

This example demonstrates:
- Processing multiple structures in batch
- Comparing properties across materials
- Simple screening workflow
"""

from pathlib import Path
from ase.io import read
from mace_inference import MACEInference
import json

EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"


def main():
    print("=" * 60)
    print("Example 5: High-Throughput Screening")
    print("=" * 60)
    
    # Initialize calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="medium", device="auto")
    
    # Find all CIF files
    print("\n2. Finding structure files...")
    cif_files = list(STRUCTURES_DIR.glob("*.cif"))
    print(f"   Found {len(cif_files)} structure files:")
    for f in cif_files:
        print(f"     - {f.name}")
    
    # Process each structure
    print("\n" + "=" * 60)
    print("3. Processing Structures")
    print("=" * 60)
    
    results = []
    
    for cif_path in cif_files:
        print(f"\n   Processing: {cif_path.name}")
        
        # Load structure
        atoms = read(str(cif_path))
        formula = atoms.get_chemical_formula()
        n_atoms = len(atoms)
        
        # Single-point calculation
        sp_result = calc.single_point(atoms)
        
        # Structure optimization
        opt_result = calc.optimize(atoms, fmax=0.05, steps=30)
        
        # Collect results
        result = {
            'filename': cif_path.name,
            'formula': formula,
            'n_atoms': n_atoms,
            'energy': sp_result['energy'],
            'energy_per_atom': sp_result['energy'] / n_atoms,
            'max_force': max(abs(sp_result['forces'].flatten())),
            'optimized_energy': opt_result.get_potential_energy() if hasattr(opt_result, 'get_potential_energy') else None,
            'volume': atoms.cell.volume,
        }
        results.append(result)
        
        print(f"     Formula: {formula}")
        print(f"     Energy: {result['energy']:.4f} eV")
        print(f"     Energy/atom: {result['energy_per_atom']:.4f} eV")
        print(f"     Max force: {result['max_force']:.4f} eV/Å")
    
    # Comparative analysis
    print("\n" + "=" * 60)
    print("4. Comparative Analysis")
    print("=" * 60)
    
    # Sort by energy per atom
    results_sorted = sorted(results, key=lambda x: x['energy_per_atom'])
    
    print("\n   Structures ranked by stability (energy/atom):")
    print(f"   {'Rank':<6} {'Name':<25} {'E/atom (eV)':<12} {'N atoms':<10}")
    print("   " + "-" * 55)
    
    for i, r in enumerate(results_sorted, 1):
        print(f"   {i:<6} {r['formula']:<25} {r['energy_per_atom']:<12.4f} {r['n_atoms']:<10}")
    
    # Screening criteria
    print("\n" + "=" * 60)
    print("5. Screening with Criteria")
    print("=" * 60)
    
    # Example: screen for well-converged structures
    converged_threshold = 0.1  # eV/Å
    
    converged = [r for r in results if r['max_force'] < converged_threshold]
    not_converged = [r for r in results if r['max_force'] >= converged_threshold]
    
    print(f"\n   Force convergence threshold: {converged_threshold} eV/Å")
    print(f"   Converged structures: {len(converged)}")
    print(f"   Need optimization: {len(not_converged)}")
    
    if not_converged:
        print("\n   Structures needing optimization:")
        for r in not_converged:
            print(f"     - {r['formula']}: max force = {r['max_force']:.4f} eV/Å")
    
    # Save results
    print("\n" + "=" * 60)
    print("6. Saving Results")
    print("=" * 60)
    
    output_file = EXAMPLES_DIR / "screening_results.json"
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"   Results saved to: {output_file.name}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"   Structures processed: {len(results)}")
    print(f"   Most stable: {results_sorted[0]['formula']} ({results_sorted[0]['energy_per_atom']:.4f} eV/atom)")
    if len(results_sorted) > 1:
        print(f"   Least stable: {results_sorted[-1]['formula']} ({results_sorted[-1]['energy_per_atom']:.4f} eV/atom)")
    print("\n✅ Example 5 completed successfully!")


if __name__ == "__main__":
    main()
