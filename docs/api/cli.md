# Command Line Interface

MACE Inference provides a command-line interface for common tasks.

## Installation

The CLI is automatically installed with the package:

```bash
pip install mace-inference
```

## Usage

```bash
mace-inference <command> [options]
```

## Commands

### single-point

Calculate energy, forces, and stress for a structure.

```bash
mace-inference single-point structure.cif [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--model` | str | `medium` | Model name or path |
| `--device` | str | `auto` | Device (auto/cpu/cuda) |
| `--d3` | flag | False | Enable D3 correction |
| `--d3-xc` | str | `pbe` | D3 functional |

**Example:**

```bash
mace-inference single-point cu_fcc.cif --model medium --device auto
```

**Output:**

```
Energy: -16.3376 eV
Energy/atom: -4.0844 eV
Max force: 0.0001 eV/Å
```

---

### optimize

Optimize atomic positions and cell.

```bash
mace-inference optimize structure.cif [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--model` | str | `medium` | Model name or path |
| `--device` | str | `auto` | Device |
| `--fmax` | float | `0.05` | Force convergence (eV/Å) |
| `--steps` | int | `500` | Maximum steps |
| `--optimizer` | str | `BFGS` | BFGS or FIRE |
| `--relax-cell` | flag | False | Also optimize cell |
| `--output` | str | None | Output file |

**Example:**

```bash
mace-inference optimize structure.cif --fmax 0.01 --relax-cell --output optimized.cif
```

---

### md

Run molecular dynamics simulation.

```bash
mace-inference md structure.cif [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--model` | str | `medium` | Model name or path |
| `--device` | str | `auto` | Device |
| `--ensemble` | str | `nvt` | nvt or npt |
| `--temperature` | float | `300` | Temperature (K) |
| `--pressure` | float | None | Pressure (GPa, NPT only) |
| `--timestep` | float | `1.0` | Time step (fs) |
| `--steps` | int | `1000` | Number of steps |
| `--output` | str | None | Trajectory output file |

**Example:**

```bash
mace-inference md structure.cif --ensemble nvt --temperature 300 --steps 5000 --output traj.xyz
```

---

### phonon

Calculate phonon properties.

```bash
mace-inference phonon structure.cif [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--model` | str | `medium` | Model name or path |
| `--device` | str | `auto` | Device |
| `--supercell` | int int int | `2 2 2` | Supercell dimensions |
| `--delta` | float | `0.01` | Displacement (Å) |
| `--temperature` | float | `300` | Temperature for thermal props |

**Example:**

```bash
mace-inference phonon si_diamond.cif --supercell 3 3 3 --temperature 300
```

---

### bulk-modulus

Calculate bulk modulus.

```bash
mace-inference bulk-modulus structure.cif [OPTIONS]
```

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--model` | str | `medium` | Model name or path |
| `--device` | str | `auto` | Device |
| `--n-points` | int | `7` | Number of E-V points |
| `--eos` | str | `birchmurnaghan` | Equation of state |

**Example:**

```bash
mace-inference bulk-modulus cu_fcc.cif --n-points 11
```

---

### adsorption

Calculate adsorption energy.

```bash
mace-inference adsorption framework.cif adsorbate [OPTIONS]
```

**Arguments:**

- `framework.cif` - Path to framework structure
- `adsorbate` - Gas molecule name (e.g., CO2, H2O, CH4)

**Options:**

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--model` | str | `medium` | Model name or path |
| `--device` | str | `auto` | Device |
| `--site` | float float float | required | Adsorption site [x, y, z] |
| `--fmax` | float | `0.05` | Force convergence |
| `--fix-framework` | flag | True | Keep framework fixed |

**Example:**

```bash
mace-inference adsorption mof.cif CO2 --site 10.0 10.0 10.0 --fmax 0.02
```

---

## Global Options

These options are available for all commands:

| Option | Description |
|--------|-------------|
| `--help` | Show help message |
| `--version` | Show version |
| `--verbose` | Verbose output |
| `--quiet` | Minimal output |

## Examples

### Batch Processing with Shell

```bash
# Process multiple structures
for f in structures/*.cif; do
    mace-inference single-point "$f" --model medium >> results.txt
done
```

### Pipeline

```bash
# Optimize then calculate phonons
mace-inference optimize input.cif --output opt.cif
mace-inference phonon opt.cif --supercell 3 3 3
```
