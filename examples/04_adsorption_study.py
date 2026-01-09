#!/usr/bin/env python
"""
Example 4: Adsorption Energy Study

This example demonstrates:
- Loading a metal-organic framework (MOF) structure
- Calculating adsorption energy for gas molecules
- Coordination number analysis
"""

from pathlib import Path
from ase.io import read
from ase import Atoms
from mace_inference import MACEInference

EXAMPLES_DIR = Path(__file__).parent
STRUCTURES_DIR = EXAMPLES_DIR / "structures"


def create_gas_molecule(molecule_type: str) -> Atoms:
    """Create common gas molecules."""
    if molecule_type == "CO2":
        return Atoms(
            "CO2",
            positions=[[0, 0, 0], [0, 0, 1.16], [0, 0, -1.16]]
        )
    elif molecule_type == "H2O":
        return Atoms(
            "H2O",
            positions=[[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]]
        )
    elif molecule_type == "N2":
        return Atoms(
            "N2",
            positions=[[0, 0, 0], [0, 0, 1.10]]
        )
    elif molecule_type == "CH4":
        return Atoms(
            "CH4",
            positions=[
                [0, 0, 0],
                [0.63, 0.63, 0.63],
                [-0.63, -0.63, 0.63],
                [0.63, -0.63, -0.63],
                [-0.63, 0.63, -0.63],
            ]
        )
    else:
        raise ValueError(f"Unknown molecule: {molecule_type}")


def main():
    print("=" * 60)
    print("Example 4: Adsorption Energy Study")
    print("=" * 60)
    
    # Initialize calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="medium", device="auto")
    
    # Load MOF structure
    print("\n2. Loading MOF framework...")
    cif_path = STRUCTURES_DIR / "cu_paddlewheel.cif"
    framework = read(str(cif_path))
    
    print(f"   Framework formula: {framework.get_chemical_formula()}")
    print(f"   Number of atoms: {len(framework)}")
    print(f"   Cell volume: {framework.cell.volume:.1f} Å³")
    
    # Single adsorption calculation
    print("\n" + "=" * 60)
    print("3. Single Molecule Adsorption (CO2)")
    print("=" * 60)
    
    co2 = create_gas_molecule("CO2")
    print(f"   Adsorbate: {co2.get_chemical_formula()}")
    
    result = calc.adsorption_energy(
        framework=framework,
        adsorbate=co2,
        site_position=[5.0, 5.0, 5.0],  # Near Cu site
        optimize=True,
        fix_framework=True  # Fix framework during optimization
    )
    
    print(f"\n   Framework energy: {result['E_mof']:.4f} eV")
    print(f"   Adsorbate energy: {result['E_gas']:.4f} eV")
    print(f"   Complex energy: {result['E_complex']:.4f} eV")
    print(f"   Adsorption energy: {result['E_ads']:.4f} eV")
    print(f"                    = {result['E_ads'] * 96.485:.2f} kJ/mol")
    
    if result['E_ads'] < 0:
        print("   → Favorable adsorption (exothermic)")
    else:
        print("   → Unfavorable adsorption (endothermic)")
    
    # Coordination analysis
    print("\n" + "=" * 60)
    print("4. Coordination Analysis")
    print("=" * 60)
    
    coord_result = calc.coordination(
        framework,
        cutoff_multiplier=1.3
    )
    
    # Print coordination info for each metal
    print("\n   Metal atom coordination:")
    for metal_idx, info in coord_result['coordination'].items():
        cn = info['coordination_number']
        print(f"     Atom {metal_idx}: CN = {cn}")
    
    # Multi-molecule comparison (quick overview)
    print("\n" + "=" * 60)
    print("5. Multi-Molecule Comparison")
    print("=" * 60)
    
    molecules = ["CO2", "H2O", "N2", "CH4"]
    results = {}
    
    for mol_name in molecules:
        mol = create_gas_molecule(mol_name)
        res = calc.adsorption_energy(
            framework=framework,
            adsorbate=mol,
            site_position=[5.0, 5.0, 5.0],
            optimize=False,  # Skip optimization for speed
            fix_framework=True
        )
        results[mol_name] = res['E_ads']
        print(f"   {mol_name:4s}: {res['E_ads']:8.4f} eV ({res['E_ads'] * 96.485:8.2f} kJ/mol)")
    
    # Selectivity
    print("\n   Adsorption selectivity (vs N2):")
    e_n2 = results["N2"]
    for mol_name, e_ads in results.items():
        if mol_name != "N2":
            # Simple selectivity based on energy difference
            selectivity = abs(e_ads - e_n2)
            print(f"     {mol_name}/N2: ΔE = {selectivity:.3f} eV")
    
    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"   Framework: Cu-paddlewheel cluster")
    print(f"   Strongest binding: {min(results.items(), key=lambda x: x[1])[0]}")
    print(f"   Weakest binding: {max(results.items(), key=lambda x: x[1])[0]}")
    print("\n✅ Example 4 completed successfully!")


if __name__ == "__main__":
    main()
