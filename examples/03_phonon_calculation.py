#!/usr/bin/env python
"""
Example 3: Phonon Calculation and Thermal Properties

This example demonstrates:
- Phonon dispersion calculation using finite displacement method
- Thermal properties (free energy, entropy, heat capacity)
"""

from pathlib import Path
from ase.io import read
from mace_inference import MACEInference
import numpy as np

EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"


def main():
    print("=" * 60)
    print("Example 3: Phonon Calculation")
    print("=" * 60)
    
    # Initialize calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="medium", device="auto")
    
    # Load structure
    print("\n2. Loading structure...")
    cif_path = STRUCTURES_DIR / "si_diamond.cif"
    atoms = read(str(cif_path))
    
    print(f"   Structure: {atoms.get_chemical_formula()}")
    print(f"   Number of atoms: {len(atoms)}")
    print(f"   Cell: {atoms.cell.lengths()}")
    
    # First optimize the structure
    print("\n3. Optimizing structure...")
    atoms_opt = calc.optimize(atoms, fmax=0.01, steps=50)
    print("   ✓ Structure optimized")
    
    # Phonon calculation
    print("\n" + "=" * 60)
    print("4. Phonon Calculation")
    print("=" * 60)
    print("   Supercell: 2x2x2")
    print("   Displacement: 0.01 Å")
    print("   This may take a minute...")
    
    result = calc.phonon(
        atoms_opt,
        supercell_matrix=[2, 2, 2],
        displacement=0.01,
        mesh=[10, 10, 10],
        temperature_range=(0, 800, 50)  # 0-800 K, 50 K steps
    )
    
    print(f"\n   Number of displacements: {len(result['phonon'].supercells_with_displacements)}")
    print("   ✓ Force constants calculated")
    
    # Thermal properties
    if 'thermal' in result:
        thermal = result['thermal']
        temps = thermal['temperatures']
        
        print("\n" + "=" * 60)
        print("5. Thermal Properties")
        print("=" * 60)
        print(f"\n   {'T (K)':<10} {'F (kJ/mol)':<15} {'S (J/mol·K)':<15} {'Cv (J/mol·K)':<15}")
        print("   " + "-" * 55)
        
        for T_target in [100, 200, 300, 500, 800]:
            if T_target <= temps[-1]:
                idx = np.argmin(np.abs(temps - T_target))
                T = temps[idx]
                F = thermal['free_energy'][idx]
                S = thermal['entropy'][idx]
                Cv = thermal['heat_capacity'][idx]
                print(f"   {T:<10.0f} {F:<15.3f} {S:<15.3f} {Cv:<15.3f}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"   Material: Si (diamond structure)")
    print(f"   Supercell size: 2×2×2 = 16 atoms")
    if 'thermal' in result:
        idx_300 = np.argmin(np.abs(temps - 300))
        print(f"   Heat capacity at 300 K: {thermal['heat_capacity'][idx_300]:.2f} J/(mol·K)")
    print("\n✅ Example 3 completed successfully!")


if __name__ == "__main__":
    main()
