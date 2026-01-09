#!/usr/bin/env python
"""
Example 6: D3 Dispersion Correction

This example demonstrates:
- Using DFT-D3 dispersion correction for van der Waals interactions
- Comparing energies with and without D3 correction
- When to use dispersion corrections
"""

from pathlib import Path
from ase.io import read
from mace_inference import MACEInference

EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"


def main():
    print("=" * 60)
    print("Example 6: D3 Dispersion Correction")
    print("=" * 60)
    
    # Load a structure where dispersion is important
    print("\n1. Loading structure...")
    cif_path = STRUCTURES_DIR / "cu_paddlewheel.cif"
    atoms = read(str(cif_path))
    
    print(f"   Structure: {atoms.get_chemical_formula()}")
    print(f"   Number of atoms: {len(atoms)}")
    
    # Initialize calculators with and without D3
    print("\n2. Initializing calculators...")
    
    calc_no_d3 = MACEInference(
        model="medium",
        device="auto",
        enable_d3=False
    )
    print("   ✓ MACE calculator (without D3)")
    
    calc_with_d3 = MACEInference(
        model="medium",
        device="auto",
        enable_d3=True,
        d3_xc="pbe",      # Exchange-correlation functional
        d3_damping="bj"   # Damping function
    )
    print("   ✓ MACE + D3 calculator (with D3)")
    
    # Compare single-point energies
    print("\n" + "=" * 60)
    print("3. Single-Point Energy Comparison")
    print("=" * 60)
    
    result_no_d3 = calc_no_d3.single_point(atoms)
    result_with_d3 = calc_with_d3.single_point(atoms)
    
    e_mace = result_no_d3['energy']
    e_mace_d3 = result_with_d3['energy']
    e_d3_correction = e_mace_d3 - e_mace
    
    print(f"\n   MACE energy:        {e_mace:.4f} eV")
    print(f"   MACE + D3 energy:   {e_mace_d3:.4f} eV")
    print(f"   D3 correction:      {e_d3_correction:.4f} eV")
    print(f"   D3 per atom:        {e_d3_correction / len(atoms):.4f} eV/atom")
    
    # Compare forces
    print("\n" + "=" * 60)
    print("4. Force Comparison")
    print("=" * 60)
    
    forces_no_d3 = result_no_d3['forces']
    forces_with_d3 = result_with_d3['forces']
    force_diff = forces_with_d3 - forces_no_d3
    
    max_force_no_d3 = abs(forces_no_d3).max()
    max_force_with_d3 = abs(forces_with_d3).max()
    max_force_diff = abs(force_diff).max()
    
    print(f"\n   Max force (MACE):      {max_force_no_d3:.4f} eV/Å")
    print(f"   Max force (MACE+D3):   {max_force_with_d3:.4f} eV/Å")
    print(f"   Max D3 force contrib:  {max_force_diff:.4f} eV/Å")
    
    # Structure optimization comparison
    print("\n" + "=" * 60)
    print("5. Optimization Comparison")
    print("=" * 60)
    
    print("\n   Optimizing without D3...")
    opt_no_d3 = calc_no_d3.optimize(atoms.copy(), fmax=0.05, steps=30)
    
    print("   Optimizing with D3...")
    opt_with_d3 = calc_with_d3.optimize(atoms.copy(), fmax=0.05, steps=30)
    
    # Compare optimized structures
    vol_no_d3 = opt_no_d3.cell.volume
    vol_with_d3 = opt_with_d3.cell.volume
    vol_diff = (vol_with_d3 - vol_no_d3) / vol_no_d3 * 100
    
    print(f"\n   Volume (no D3):     {vol_no_d3:.2f} Å³")
    print(f"   Volume (with D3):   {vol_with_d3:.2f} Å³")
    print(f"   Volume change:      {vol_diff:+.2f}%")
    
    # When to use D3
    print("\n" + "=" * 60)
    print("6. When to Use D3 Correction")
    print("=" * 60)
    
    print("""
   D3 dispersion correction is recommended for:
   
   ✓ Van der Waals crystals (molecular crystals, layered materials)
   ✓ Metal-organic frameworks (MOFs)
   ✓ Adsorption on surfaces
   ✓ Molecular dimers and complexes
   ✓ Graphene and 2D materials with stacking
   
   D3 may be less important for:
   
   ○ Covalent/ionic bulk crystals (Si, NaCl)
   ○ Metals with strong metallic bonding
   ○ Systems where MACE was trained with D3-corrected data
   
   Note: Check if your MACE model was trained with D3 correction!
   If so, adding D3 again will double-count dispersion.
""")
    
    # Summary
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"   D3 energy contribution: {e_d3_correction:.4f} eV ({e_d3_correction/len(atoms)*1000:.2f} meV/atom)")
    print(f"   D3 effect on structure: {vol_diff:+.2f}% volume change")
    print("\n✅ Example 6 completed successfully!")


if __name__ == "__main__":
    main()
