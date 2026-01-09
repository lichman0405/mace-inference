#!/usr/bin/env python
"""
Example 2: Molecular Dynamics - NVT and NPT simulations

This example demonstrates how to run MD simulations using MACE:
- NVT ensemble (constant temperature)
- NPT ensemble (constant temperature and pressure)
- Trajectory analysis
"""

from pathlib import Path
from ase.io import read
from mace_inference import MACEInference
import numpy as np

EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"


def main():
    print("=" * 60)
    print("Example 2: Molecular Dynamics Simulations")
    print("=" * 60)
    
    # Initialize calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="medium", device="auto")
    
    # Load and prepare structure
    print("\n2. Loading structure...")
    cif_path = STRUCTURES_DIR / "cu_fcc.cif"
    atoms = read(str(cif_path))
    atoms = atoms * (3, 3, 3)  # 3x3x3 supercell for MD
    
    print(f"   Structure: {atoms.get_chemical_formula()}")
    print(f"   Number of atoms: {len(atoms)}")
    print(f"   Initial volume: {atoms.get_volume():.2f} Å³")
    
    # NVT Molecular Dynamics
    print("\n" + "=" * 60)
    print("3. NVT Molecular Dynamics (300 K)")
    print("=" * 60)
    
    print("   Running NVT MD (100 steps, timestep=1 fs)...")
    final_nvt = calc.run_md(
        atoms,
        ensemble="nvt",
        temperature_K=300,
        steps=100,
        timestep=1.0,
        log_interval=20
    )
    
    result_nvt = calc.single_point(final_nvt)
    print(f"\n   Final energy: {result_nvt['energy']:.4f} eV")
    print(f"   Final volume: {final_nvt.get_volume():.2f} Å³ (unchanged in NVT)")
    print("   ✓ NVT simulation completed")
    
    # NPT Molecular Dynamics
    print("\n" + "=" * 60)
    print("4. NPT Molecular Dynamics (300 K, 0 GPa)")
    print("=" * 60)
    
    print("   Running NPT MD (100 steps, timestep=1 fs)...")
    atoms_npt = atoms.copy()  # Fresh copy
    
    final_npt = calc.run_md(
        atoms_npt,
        ensemble="npt",
        temperature_K=300,
        pressure_GPa=0.0,
        steps=100,
        timestep=1.0,
        log_interval=20
    )
    
    result_npt = calc.single_point(final_npt)
    volume_change = (final_npt.get_volume() - atoms.get_volume()) / atoms.get_volume() * 100
    
    print(f"\n   Final energy: {result_npt['energy']:.4f} eV")
    print(f"   Final volume: {final_npt.get_volume():.2f} Å³")
    print(f"   Volume change: {volume_change:+.2f}%")
    print("   ✓ NPT simulation completed")
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"   NVT final energy: {result_nvt['energy']:.4f} eV")
    print(f"   NPT final energy: {result_npt['energy']:.4f} eV")
    print(f"   NPT volume change: {volume_change:+.2f}%")
    print("\n✅ Example 2 completed successfully!")


if __name__ == "__main__":
    main()
