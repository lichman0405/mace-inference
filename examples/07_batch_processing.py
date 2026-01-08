#!/usr/bin/env python
"""
Example 07: Batch Processing and Progress Callbacks

This example demonstrates:
1. Batch single-point calculations on multiple structures
2. Batch structure optimization
3. Using progress callbacks for long-running tasks
4. Error handling in batch operations

Batch processing is more efficient than individual calls as it reuses
the calculator and provides better progress tracking.
"""

from mace_inference import MACEInference
from ase.build import bulk, molecule, fcc111, add_adsorbate
from ase import Atoms
import numpy as np
import tempfile
from pathlib import Path


def simple_progress(current: int, total: int) -> None:
    """Simple progress callback that prints to console."""
    percent = 100.0 * current / total
    bar_length = 30
    filled = int(bar_length * current / total)
    bar = "=" * filled + "-" * (bar_length - filled)
    print(f"\rProgress: [{bar}] {percent:5.1f}% ({current}/{total})", end="", flush=True)
    if current == total:
        print()  # New line when complete


def main():
    print("=" * 60)
    print("MACE Inference: Batch Processing Example")
    print("=" * 60)
    
    # Initialize calculator
    print("\n1. Initializing MACE calculator...")
    calc = MACEInference(model="small", device="auto")
    print(f"   Calculator: {calc}")
    
    # =========================================================================
    # Part 1: Batch Single-Point Calculations
    # =========================================================================
    print("\n" + "=" * 60)
    print("Part 1: Batch Single-Point Calculations")
    print("=" * 60)
    
    # Create multiple test structures
    structures = [
        bulk("Cu", "fcc", a=3.6),           # FCC copper
        bulk("Al", "fcc", a=4.05),          # FCC aluminum
        bulk("Fe", "bcc", a=2.87),          # BCC iron
        bulk("Ni", "fcc", a=3.52),          # FCC nickel
        bulk("Pt", "fcc", a=3.92),          # FCC platinum
    ]
    
    print(f"\n   Created {len(structures)} bulk metal structures")
    
    # Perform batch single-point calculations with progress callback
    print("\n   Running batch single-point calculations...")
    results = calc.batch_single_point(
        structures,
        progress_callback=simple_progress
    )
    
    # Display results
    print("\n   Results:")
    print("   " + "-" * 50)
    elements = ["Cu", "Al", "Fe", "Ni", "Pt"]
    for elem, result in zip(elements, results):
        if result.get("success", True):
            energy = result["energy"]
            forces_max = np.max(np.abs(result["forces"]))
            print(f"   {elem}: E = {energy:10.4f} eV, max|F| = {forces_max:.4f} eV/A")
        else:
            print(f"   {elem}: ERROR - {result.get('error', 'Unknown')}")
    
    # =========================================================================
    # Part 2: Batch Structure Optimization
    # =========================================================================
    print("\n" + "=" * 60)
    print("Part 2: Batch Structure Optimization")
    print("=" * 60)
    
    # Create structures with small perturbations
    perturbed_structures = []
    for i, atoms in enumerate(structures[:3]):
        perturbed = atoms.copy()
        # Add random perturbations to positions
        perturbed.positions += np.random.randn(*perturbed.positions.shape) * 0.1
        perturbed_structures.append(perturbed)
    
    print(f"\n   Created {len(perturbed_structures)} perturbed structures")
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "optimized"
        
        print("\n   Running batch optimization...")
        opt_results = calc.batch_optimize(
            perturbed_structures,
            fmax=0.05,
            steps=100,
            output_dir=str(output_dir),
            progress_callback=simple_progress
        )
        
        # Display results
        print("\n   Optimization Results:")
        print("   " + "-" * 50)
        for elem, result in zip(["Cu", "Al", "Fe"], opt_results):
            if result.get("success", False):
                print(f"   {elem}: E = {result['final_energy']:10.4f} eV, "
                      f"max|F| = {result['final_max_force']:.4f} eV/A, "
                      f"converged = {result['converged']}")
            else:
                print(f"   {elem}: ERROR - {result.get('error', 'Unknown')}")
    
    # =========================================================================
    # Part 3: Progress Callbacks for MD Simulation
    # =========================================================================
    print("\n" + "=" * 60)
    print("Part 3: MD Simulation with Progress Callback")
    print("=" * 60)
    
    # Create a small copper cluster
    cu_bulk = bulk("Cu", "fcc", a=3.6) * (2, 2, 2)
    print(f"\n   Structure: Cu supercell with {len(cu_bulk)} atoms")
    
    print("\n   Running NVT MD (100 steps)...")
    
    # Custom progress callback with timing
    import time
    start_time = time.time()
    
    def md_progress(current: int, total: int) -> None:
        elapsed = time.time() - start_time
        if current > 0:
            eta = elapsed / current * (total - current)
            print(f"\r   Step {current:5d}/{total} | "
                  f"Elapsed: {elapsed:6.1f}s | ETA: {eta:6.1f}s", end="", flush=True)
        if current == total:
            print()
    
    with tempfile.NamedTemporaryFile(suffix=".traj", delete=False) as traj_file:
        final_atoms = calc.run_md(
            cu_bulk,
            ensemble="nvt",
            temperature_K=300,
            steps=100,
            timestep=1.0,
            log_interval=10,
            trajectory=traj_file.name,
            progress_callback=md_progress
        )
    
    print(f"\n   Final temperature: ~300 K (target)")
    print(f"   Total time: {time.time() - start_time:.1f}s")
    
    # =========================================================================
    # Part 4: Error Handling in Batch Operations
    # =========================================================================
    print("\n" + "=" * 60)
    print("Part 4: Error Handling in Batch Operations")
    print("=" * 60)
    
    # Create a mix of valid and "problematic" structures
    mixed_structures = [
        bulk("Cu", "fcc", a=3.6),           # Valid
        bulk("Ag", "fcc", a=4.09),          # Valid
        Atoms("X"),                          # Invalid element (may cause error)
    ]
    
    print(f"\n   Testing batch calculations with {len(mixed_structures)} structures")
    print("   (including one potentially problematic structure)")
    
    results = calc.batch_single_point(mixed_structures)
    
    print("\n   Results:")
    labels = ["Cu (valid)", "Ag (valid)", "X (invalid)"]
    for label, result in zip(labels, results):
        if result.get("success", True):
            print(f"   {label}: E = {result['energy']:.4f} eV")
        else:
            print(f"   {label}: FAILED - {result.get('error', 'Unknown error')}")
    
    # Count successes
    successes = sum(1 for r in results if r.get("success", True))
    print(f"\n   Summary: {successes}/{len(results)} successful calculations")
    
    print("\n" + "=" * 60)
    print("Batch processing example complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
