#!/usr/bin/env python
"""
Example 1: Basic Usage - Single-point energy and structure optimization

This example demonstrates the basic functionality of mace-inference:
- Loading a structure from CIF file
- Calculating single-point energy
- Structure optimization
"""

from pathlib import Path
from ase.io import read
from mace_inference import MACEInference
import numpy as np

# Get the examples directory
EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"


def main():
    print("=" * 60)
    print("Example 1: Basic Usage")
    print("=" * 60)
    
    # Initialize MACE calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="medium", device="auto")
    print(f"   Model: {calc.model_name}")
    print(f"   Device: {calc.device}")
    
    # Load structure from CIF file
    print("\n2. Loading structure from CIF file...")
    cif_path = STRUCTURES_DIR / "cu_fcc.cif"
    atoms = read(str(cif_path))
    atoms = atoms * (2, 2, 2)  # Create 2x2x2 supercell
    
    print(f"   Source: {cif_path.name}")
    print(f"   Formula: {atoms.get_chemical_formula()}")
    print(f"   Number of atoms: {len(atoms)}")
    print(f"   Cell: {atoms.cell.lengths()}")
    
    # Single-point energy calculation
    print("\n3. Single-Point Energy Calculation")
    print("-" * 40)
    result = calc.single_point(atoms)
    
    print(f"   Total Energy:    {result['energy']:.6f} eV")
    print(f"   Energy/atom:     {result['energy_per_atom']:.6f} eV")
    print(f"   Max Force:       {result['max_force']:.6f} eV/Å")
    print(f"   RMS Force:       {result['rms_force']:.6f} eV/Å")
    if result.get('pressure_GPa') is not None:
        print(f"   Pressure:        {result['pressure_GPa']:.4f} GPa")
    
    # Perturb structure slightly
    print("\n4. Perturbing Structure...")
    np.random.seed(42)
    atoms.rattle(stdev=0.1)  # Random displacement
    
    result_perturbed = calc.single_point(atoms)
    print(f"   Energy after perturbation: {result_perturbed['energy']:.6f} eV")
    print(f"   Max Force after perturbation: {result_perturbed['max_force']:.6f} eV/Å")
    
    # Optimize structure
    print("\n5. Structure Optimization")
    print("-" * 40)
    optimized = calc.optimize(
        atoms,
        fmax=0.01,
        steps=200,
        optimizer="LBFGS"
    )
    
    result_opt = calc.single_point(optimized)
    print(f"   Energy after optimization: {result_opt['energy']:.6f} eV")
    print(f"   Max Force after optimization: {result_opt['max_force']:.6f} eV/Å")
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"   Energy change: {result_opt['energy'] - result_perturbed['energy']:.6f} eV")
    print(f"   Force reduction: {result_perturbed['max_force']:.4f} → {result_opt['max_force']:.6f} eV/Å")
    print("\n✅ Example 1 completed successfully!")


if __name__ == "__main__":
    main()
