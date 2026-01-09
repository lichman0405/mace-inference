# Structure Optimization

Structure optimization (geometry relaxation) minimizes the potential energy by adjusting atomic positions and optionally cell parameters.

## Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read

calc = MACEInference(model="medium", device="auto")
atoms = read("structure.cif")

result = calc.optimize(atoms, fmax=0.01)
optimized_atoms = result['atoms']
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `atoms` | `Atoms` | required | Structure to optimize |
| `fmax` | `float` | `0.05` | Force convergence criterion (eV/Å) |
| `steps` | `int` | `500` | Maximum optimization steps |
| `optimizer` | `str` | `"BFGS"` | Algorithm: `"BFGS"` or `"FIRE"` |
| `relax_cell` | `bool` | `False` | Also optimize cell parameters |
| `fix_symmetry` | `bool` | `False` | Preserve crystal symmetry |
| `trajectory` | `str` | `None` | Save trajectory to file |

## Return Values

| Key | Type | Description |
|-----|------|-------------|
| `atoms` | `Atoms` | Optimized structure |
| `energy` | `float` | Final energy (eV) |
| `forces` | `ndarray` | Final forces (eV/Å) |
| `max_force` | `float` | Final maximum force (eV/Å) |
| `converged` | `bool` | Whether optimization converged |
| `steps` | `int` | Number of steps taken |
| `trajectory` | `list` | List of atoms at each step |

## Examples

### Basic Optimization

```python
optimized = calc.optimize(atoms, fmax=0.01)

# Get properties from the optimized structure
print(f"Final energy: {optimized.get_potential_energy():.6f} eV")
max_force = np.max(np.linalg.norm(optimized.get_forces(), axis=1))
print(f"Final max force: {max_force:.6f} eV/Å")
```

### Full Cell Relaxation

Optimize both atomic positions and cell parameters:

```python
optimized = calc.optimize(
    atoms,
    fmax=0.01,
    optimize_cell=True  # Also optimize cell vectors
)

# Compare cell parameters
print("Original cell:")
print(atoms.cell[:])

print("\nOptimized cell:")
print(optimized.cell[:])
```

### Choosing Optimizer

```python
# LBFGS - Default, memory-efficient quasi-Newton method
optimized = calc.optimize(atoms, optimizer="LBFGS", fmax=0.01)

# BFGS - Traditional quasi-Newton method
optimized = calc.optimize(atoms, optimizer="BFGS", fmax=0.01)

# FIRE - Better for difficult cases, molecular dynamics-based
optimized = calc.optimize(atoms, optimizer="FIRE", fmax=0.01, steps=1000)
```

!!! tip "When to use FIRE"
    FIRE is often better for:
    
    - Highly distorted structures
    - Systems with soft modes
    - Cases where BFGS oscillates or fails

### Saving Output

Save the optimized structure directly:

```python
optimized = calc.optimize(
    atoms,
    fmax=0.01,
    output="optimized.cif",     # Save structure to file
    trajectory="opt.traj"       # Save trajectory to file
)

# Read trajectory for analysis
from ase.io import Trajectory
import matplotlib.pyplot as plt

traj = Trajectory("opt.traj")
print(f"Trajectory length: {len(traj)} frames")

# Analyze convergence
energies = [a.get_potential_energy() for a in traj]
plt.plot(energies)
plt.xlabel("Step")
plt.ylabel("Energy (eV)")
plt.savefig("convergence.png")
```

### Saving Trajectory

You can save the optimization trajectory for later analysis:

```python
optimized = calc.optimize(
    atoms,
    fmax=0.01,
    trajectory="optimization.traj"
)
```

!!! warning "Symmetry requirements"
    This requires the `spglib` package:
    ```bash
    pip install spglib
    ```

### Handling Non-Convergence

```python
import numpy as np

optimized = calc.optimize(atoms, fmax=0.01, steps=500)

# Check convergence by examining forces
max_force = np.max(np.linalg.norm(optimized.get_forces(), axis=1))
if max_force > 0.01:
    print(f"Warning: May not have converged, max force: {max_force:.4f} eV/Å")
    
    # Options:
    # 1. Continue optimization
    optimized = calc.optimize(optimized, fmax=0.01, steps=500)
    
    # 2. Try different optimizer
    optimized = calc.optimize(optimized, optimizer="FIRE", fmax=0.01)
```

## Optimization Strategies

### Two-Stage Optimization

For difficult structures:

```python
# Stage 1: Quick rough optimization with FIRE
rough = calc.optimize(atoms, optimizer="FIRE", fmax=0.1, steps=200)

# Stage 2: Fine optimization with LBFGS
final = calc.optimize(rough, optimizer="LBFGS", fmax=0.01)
```

### Variable Cell with Fixed Volume

Not directly supported, but you can use ASE filters:

```python
from ase.constraints import StrainFilter

# Fix volume, allow shape changes
sf = StrainFilter(atoms, mask=[1, 1, 1, 1, 1, 1])
# Then optimize...
```

## Comparing Initial and Final

```python
from ase.io import read, write

atoms = read("initial.cif")
result = calc.optimize(atoms, fmax=0.01, relax_cell=True)

# Energy change
initial_E = calc.single_point(atoms)['energy']
final_E = result['energy']
print(f"ΔE = {final_E - initial_E:.4f} eV")

# Volume change
initial_V = atoms.get_volume()
final_V = result['atoms'].get_volume()
print(f"ΔV = {final_V - initial_V:.2f} Å³ ({(final_V/initial_V - 1)*100:.1f}%)")

# Save optimized structure
write("optimized.cif", result['atoms'])
```

## Common Issues

### Optimization Oscillates

Try FIRE optimizer or reduce step size (requires ASE modification).

### Atoms Move Too Far

Your initial structure may be very far from a minimum. Try running a short MD at low temperature first to relax.

### Cell Optimization Fails

Cell optimization is more sensitive. Try:

1. First optimize positions only
2. Then optimize cell with fixed positions
3. Finally optimize both together
