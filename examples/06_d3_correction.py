"""
Example 6: D3 Dispersion Correction

This example demonstrates the use of DFT-D3 dispersion correction with MACE.

D3 correction is important for:
- Van der Waals interactions
- Molecular crystals
- Layered materials (graphite, MoS2)
- Gas adsorption in porous materials
- Weakly bound complexes

Requirements:
    pip install mace-inference[d3]
    # or
    pip install torch-dftd
"""

from ase.build import bulk, molecule
from ase.io import write
from mace_inference import MACEInference
from mace_inference.utils import check_d3_available
import numpy as np

# Check D3 availability
print("=== D3 Dispersion Correction Example ===\n")

if check_d3_available():
    print("✓ torch-dftd is installed, D3 correction available")
else:
    print("✗ torch-dftd not installed")
    print("  Install with: pip install mace-inference[d3]")
    print("  Continuing without D3 for demonstration...\n")

# =============================================================================
# Part 1: Compare with and without D3 for layered materials
# =============================================================================
print("\n" + "="*60)
print("Part 1: Layered Material (Graphite-like)")
print("="*60)

# Create a simple layered structure (Graphene bilayer proxy)
# In real applications, use actual graphite or MoS2 structures
from ase import Atoms

# Simple bilayer structure
layer1 = Atoms('C2', positions=[[0, 0, 0], [1.42, 0, 0]], 
               cell=[2.46, 4.26, 10.0], pbc=True)
layer2 = layer1.copy()
layer2.translate([0, 0, 3.35])  # Interlayer distance

bilayer = layer1 + layer2
bilayer.center(vacuum=5.0, axis=2)

print(f"Bilayer structure: {len(bilayer)} atoms")
print(f"Cell: {bilayer.cell.lengths()}")

# Calculate without D3
print("\nCalculating WITHOUT D3 correction...")
calc_no_d3 = MACEInference(model="medium", device="auto", enable_d3=False)
result_no_d3 = calc_no_d3.single_point(bilayer)
print(f"  Energy (no D3):     {result_no_d3['energy']:.6f} eV")
print(f"  Energy/atom (no D3): {result_no_d3['energy_per_atom']:.6f} eV")

# Calculate with D3
if check_d3_available():
    print("\nCalculating WITH D3 correction...")
    calc_d3 = MACEInference(
        model="medium", 
        device="auto", 
        enable_d3=True,
        d3_damping="bj",  # Becke-Johnson damping (recommended)
        d3_xc="pbe"       # PBE functional parameters
    )
    result_d3 = calc_d3.single_point(bilayer)
    print(f"  Energy (with D3):     {result_d3['energy']:.6f} eV")
    print(f"  Energy/atom (with D3): {result_d3['energy_per_atom']:.6f} eV")
    
    # D3 contribution
    d3_contribution = result_d3['energy'] - result_no_d3['energy']
    print(f"\n  D3 contribution: {d3_contribution:.6f} eV")
    print(f"  D3 per atom:     {d3_contribution/len(bilayer):.6f} eV")

# =============================================================================
# Part 2: Different D3 damping functions
# =============================================================================
print("\n" + "="*60)
print("Part 2: D3 Damping Function Comparison")
print("="*60)

if check_d3_available():
    damping_types = ["bj", "zero", "zerom", "bjm"]
    
    # Use a molecular system for clearer comparison
    benzene = molecule('C6H6')
    benzene.center(vacuum=10.0)
    benzene.pbc = True
    
    print(f"\nTest molecule: Benzene ({len(benzene)} atoms)")
    print("-" * 40)
    
    for damping in damping_types:
        try:
            calc = MACEInference(
                model="medium",
                device="auto",
                enable_d3=True,
                d3_damping=damping,
                d3_xc="pbe"
            )
            result = calc.single_point(benzene)
            print(f"  {damping:>6} damping: E = {result['energy']:.6f} eV")
        except Exception as e:
            print(f"  {damping:>6} damping: Error - {e}")

# =============================================================================
# Part 3: D3 for different XC functionals
# =============================================================================
print("\n" + "="*60)
print("Part 3: D3 with Different XC Functional Parameters")
print("="*60)

if check_d3_available():
    xc_functionals = ["pbe", "b3lyp", "tpss", "revpbe"]
    
    print(f"\nTest structure: Benzene")
    print("-" * 40)
    
    for xc in xc_functionals:
        try:
            calc = MACEInference(
                model="medium",
                device="auto",
                enable_d3=True,
                d3_damping="bj",
                d3_xc=xc
            )
            result = calc.single_point(benzene)
            print(f"  {xc:>8} functional: E = {result['energy']:.6f} eV")
        except Exception as e:
            print(f"  {xc:>8} functional: Error - {e}")

# =============================================================================
# Part 4: Interlayer binding energy with D3
# =============================================================================
print("\n" + "="*60)
print("Part 4: Interlayer Binding Energy Curve")
print("="*60)

if check_d3_available():
    # Vary interlayer distance
    distances = np.linspace(2.5, 5.0, 11)  # Å
    energies_no_d3 = []
    energies_d3 = []
    
    print("\nScanning interlayer distance...")
    print(f"{'Distance (Å)':<15} {'E (no D3)':<15} {'E (D3)':<15}")
    print("-" * 45)
    
    for d in distances:
        # Create bilayer with specified distance
        test_bilayer = layer1.copy()
        layer2_copy = layer1.copy()
        layer2_copy.translate([0, 0, d])
        test_bilayer = test_bilayer + layer2_copy
        test_bilayer.center(vacuum=5.0, axis=2)
        
        # No D3
        result1 = calc_no_d3.single_point(test_bilayer)
        energies_no_d3.append(result1['energy'])
        
        # With D3
        result2 = calc_d3.single_point(test_bilayer)
        energies_d3.append(result2['energy'])
        
        print(f"{d:<15.2f} {result1['energy']:<15.6f} {result2['energy']:<15.6f}")
    
    # Find equilibrium distance
    idx_min_d3 = np.argmin(energies_d3)
    print(f"\nEquilibrium distance (with D3): {distances[idx_min_d3]:.2f} Å")
    print(f"Binding energy: {min(energies_d3) - max(energies_d3):.4f} eV")

# =============================================================================
# Part 5: Best practices
# =============================================================================
print("\n" + "="*60)
print("D3 Correction: Best Practices")
print("="*60)
print("""
1. When to use D3:
   - Molecular crystals and organic solids
   - Layered materials (graphite, MoS2, BN)
   - Gas adsorption in MOFs/zeolites
   - Molecular complexes and dimers
   - Any system with significant dispersion interactions

2. Recommended settings:
   - d3_damping="bj"  (Becke-Johnson, most widely used)
   - d3_xc="pbe"      (match your reference DFT functional)

3. Performance notes:
   - D3 adds ~10-20% computational overhead
   - Scales well with system size
   - GPU acceleration supported

4. Caveats:
   - D3 parameters are optimized for specific DFT functionals
   - May double-count some correlation already in MACE
   - Test on benchmark systems when accuracy is critical
""")

print("\n✓ D3 correction example completed!")
