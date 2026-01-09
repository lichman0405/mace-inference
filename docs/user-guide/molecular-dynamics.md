# Molecular Dynamics

Molecular dynamics (MD) simulates atomic motion over time, sampling configurations at finite temperature.

## Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read

calc = MACEInference(model="medium", device="auto")
atoms = read("structure.cif")

# NVT simulation at 300 K
final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,
    timestep=1.0,
    steps=1000
)
```

## Ensembles

MACE Inference supports two common ensembles:

### NVT (Canonical)

Constant number of particles, volume, and temperature.

```python
final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,  # Kelvin
    timestep=1.0,       # fs
    steps=5000,
    taut=100            # Temperature coupling time (fs)
)
```

Uses Langevin thermostat for temperature control.

### NPT (Isothermal-Isobaric)

Constant number of particles, pressure, and temperature.

```python
final = calc.run_md(
    atoms,
    ensemble="npt",
    temperature_K=300,   # Kelvin
    pressure_GPa=0.0001, # GPa (0.1 MPa ≈ 1 atm)
    timestep=1.0,        # fs
    steps=5000
)
```

Uses Berendsen barostat for pressure control.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `atoms` | `Atoms` | required | Initial structure |
| `ensemble` | `str` | `"nvt"` | `"nvt"` or `"npt"` |
| `temperature_K` | `float` | `300` | Temperature in Kelvin |
| `pressure_GPa` | `float` | `None` | Pressure in GPa (NPT only) |
| `timestep` | `float` | `1.0` | Time step in fs |
| `steps` | `int` | `1000` | Number of MD steps |
| `taut` | `float` | `None` | Temperature coupling time (fs) |
| `taup` | `float` | `None` | Pressure coupling time (fs) |
| `trajectory` | `str` | `None` | Save trajectory to file |
| `log_interval` | `int` | `100` | Print every N steps |
| `progress_callback` | `Callable` | `None` | Callback function |

## Return Values

The `run_md` method returns the **final Atoms object** after the simulation.

## Examples

### Basic NVT Simulation

```python
# Run MD and get final structure
final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,
    timestep=1.0,
    steps=5000,
    trajectory="md.traj",  # Save trajectory to file
    log_interval=100
)

# Final structure is an Atoms object
print(f"Final energy: {final.get_potential_energy():.4f} eV")
```

### Using Trajectory File

```python
from ase.io import read, Trajectory

# Run MD with trajectory saved
final = calc.run_md(
    atoms, 
    ensemble="nvt", 
    temperature_K=300, 
    steps=5000,
    trajectory="md.traj"
)

# Read trajectory for analysis
traj = Trajectory("md.traj")
print(f"Trajectory frames: {len(traj)}")

# Analyze energies
energies = [frame.get_potential_energy() for frame in traj]
```

### NPT for Volume Equilibration

```python
# Equilibrate at 1 atm
final = calc.run_md(
    atoms,
    ensemble="npt",
    temperature_K=300,
    pressure_GPa=0.0001,  # ~1 atm
    timestep=1.0,
    steps=10000,
    trajectory="npt.traj"
)

# Check volume change
print(f"Initial volume: {atoms.get_volume():.2f} Å³")
print(f"Final volume: {final.get_volume():.2f} Å³")
```

### High-Temperature Dynamics

```python
# Simulate at 1000 K
final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=1000,
    timestep=0.5,     # Smaller timestep for stability
    steps=10000,
    taut=50           # Tighter coupling at high T
)
```

### With Progress Callback

```python
def progress(current, total):
    if current % 100 == 0:
        print(f"MD: {current}/{total} steps ({100*current/total:.1f}%)")

final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,
    steps=5000,
    progress_callback=progress
)
```

## Workflow Tips

### Always Optimize First

```python
# Step 1: Optimize (returns Atoms object)
optimized = calc.optimize(atoms, fmax=0.01)

# Step 2: Run MD on relaxed structure
final = calc.run_md(optimized, ensemble="nvt", temperature_K=300, steps=5000)
```

### Equilibration + Production

```python
# Equilibration run
equil = calc.run_md(atoms, ensemble="nvt", temperature_K=300, steps=2000)

# Production run from equilibrated structure
final = calc.run_md(equil, ensemble="nvt", temperature_K=300, steps=10000, trajectory="prod.traj")

# Analyze production trajectory from file
```

### Temperature Ramp

```python
temps = [100, 200, 300, 400, 500]
current_atoms = atoms.copy()

for T in temps:
    current_atoms = calc.run_md(
        current_atoms,
        ensemble="nvt",
        temperature_K=T,
        steps=1000
    )
    print(f"T={T} K: E={current_atoms.get_potential_energy():.4f} eV")
```

## Common Issues

### Temperature Instability

If temperature fluctuates wildly:

1. Reduce timestep
2. Increase friction
3. Check initial structure (optimize first)

### Atoms "Explode"

If atoms fly apart:

1. Structure may not be relaxed - optimize first
2. Timestep too large
3. Check for atoms too close together

### Slow Simulation

For faster simulations:

1. Use GPU: `device="cuda"`
2. Increase `save_interval`
3. Use float32: `default_dtype="float32"`
