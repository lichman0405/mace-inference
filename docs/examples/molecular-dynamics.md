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
opt = calc.optimize(atoms, fmax=0.01)
atoms = opt['atoms']
```

### NVT Simulation

```python
nvt_result = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature=300,      # K
    timestep=1.0,         # fs
    steps=1000,
    friction=0.01,        # Langevin friction
    save_interval=10,
    log_interval=100
)

print(f"Trajectory frames: {len(nvt_result['trajectory'])}")
print(f"Final T: {nvt_result['temperatures'][-1]:.1f} K")
```

### NPT Simulation

```python
npt_result = calc.run_md(
    atoms,
    ensemble="npt",
    temperature=300,
    pressure=0.0001,      # GPa (~1 atm)
    timestep=1.0,
    steps=1000
)

# Check volume equilibration
V_initial = atoms.get_volume()
V_final = npt_result['final_atoms'].get_volume()
print(f"Volume change: {(V_final/V_initial - 1)*100:.2f}%")
```

### Save Trajectory

```python
write("trajectory.xyz", nvt_result['trajectory'])
```

### Analyze Temperature

```python
import numpy as np

temps = np.array(nvt_result['temperatures'])
print(f"Temperature: {temps.mean():.1f} ± {temps.std():.1f} K")
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
