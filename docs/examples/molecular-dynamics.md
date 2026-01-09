# Example 2: Molecular Dynamics

This example demonstrates NVT and NPT molecular dynamics simulations.

## Source Code

See [`examples/02_molecular_dynamics.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/02_molecular_dynamics.py)

## What You'll Learn

- Running NVT (constant temperature) dynamics
- Running NPT (constant pressure) dynamics
- Temperature and pressure control
- Saving and analyzing trajectories

## Code Walkthrough

### Setup

```python
from mace_inference import MACEInference
from ase.io import read, write

calc = MACEInference(model="medium", device="auto")
atoms = read("structures/cu_fcc.cif")
```

### Optimize First

Always relax the structure before running MD:

```python
# optimize() returns Atoms object directly
atoms = calc.optimize(atoms, fmax=0.01)
```

### NVT Simulation

```python
# run_md() returns the final Atoms object
final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,    # Kelvin
    timestep=1.0,         # fs
    steps=1000,
    trajectory="nvt.traj", # Save trajectory to file
    log_interval=100
)

print(f"Final energy: {final.get_potential_energy():.4f} eV")
```

### NPT Simulation

```python
final = calc.run_md(
    atoms,
    ensemble="npt",
    temperature_K=300,
    pressure_GPa=0.0001,  # GPa (~1 atm)
    timestep=1.0,
    steps=1000
)

# Check volume equilibration
V_initial = atoms.get_volume()
V_final = final.get_volume()
print(f"Volume change: {(V_final/V_initial - 1)*100:.2f}%")
```

### Save and Analyze Trajectory

```python
from ase.io import Trajectory

# Read saved trajectory
traj = Trajectory("nvt.traj")
print(f"Trajectory frames: {len(traj)}")

# Analyze energies
energies = [frame.get_potential_energy() for frame in traj]
```

## Expected Output

```
============================================================
Example 2: Molecular Dynamics
============================================================

1. Loading and optimizing structure...
   Structure: Cu4
   Optimized in 1 steps

2. NVT Simulation (300 K)...
   Step 100: T = 287.3 K, E = -16.12 eV
   Step 200: T = 312.5 K, E = -16.08 eV
   ...
   Trajectory saved: 100 frames
   Average T: 298.5 ± 45.2 K

3. NPT Simulation (300 K, 1 atm)...
   Step 100: T = 295.1 K, V = 46.8 Å³
   Step 200: T = 308.2 K, V = 46.9 Å³
   ...
   Volume change: +0.2%

✅ Example 2 completed successfully!
```

## Tips

1. **Timestep**: 1 fs is safe for most systems; reduce for light atoms (H)
2. **Equilibration**: Discard initial frames for analysis
3. **Friction**: Higher friction = stronger temperature control but less realistic dynamics
4. **Pressure units**: MACE Inference uses GPa (1 atm ≈ 0.0001 GPa)
