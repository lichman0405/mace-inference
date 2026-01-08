# MACE Inference API Reference

> **Version**: 0.1.0  
> **Python**: >= 3.9  
> **Last Updated**: 2026-01-08

This document provides complete API reference for the `mace-inference` Python library.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [MACEInference Class](#maceinference-class)
   - [Constructor](#constructor)
   - [Single-Point Calculations](#single-point-calculations)
   - [Structure Optimization](#structure-optimization)
   - [Molecular Dynamics](#molecular-dynamics)
   - [Phonon Calculations](#phonon-calculations)
   - [Mechanical Properties](#mechanical-properties)
   - [Adsorption Studies](#adsorption-studies)
   - [Batch Processing](#batch-processing)
3. [Task Functions](#task-functions)
4. [Type Definitions](#type-definitions)
5. [Utility Functions](#utility-functions)
6. [CLI Reference](#cli-reference)

---

## Quick Start

```python
from mace_inference import MACEInference

# Initialize calculator
calc = MACEInference(model="medium", device="auto")

# Single-point energy
result = calc.single_point("structure.cif")
print(f"Energy: {result['energy']:.4f} eV")

# Structure optimization
optimized = calc.optimize("structure.cif", fmax=0.05)

# Batch processing with progress
def progress(current, total):
    print(f"Processing {current}/{total}")

results = calc.batch_single_point(structures, progress_callback=progress)
```

---

## MACEInference Class

The main interface for all MACE inference tasks.

```python
from mace_inference import MACEInference
```

### Constructor

```python
MACEInference(
    model: str = "medium",
    device: Literal["auto", "cpu", "cuda"] = "auto",
    enable_d3: bool = False,
    d3_damping: str = "bj",
    d3_xc: str = "pbe",
    default_dtype: str = "float64"
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model` | `str` | `"medium"` | Model name (`"small"`, `"medium"`, `"large"`) or path to custom model file |
| `device` | `str` | `"auto"` | Compute device: `"auto"` (auto-detect), `"cpu"`, or `"cuda"` |
| `enable_d3` | `bool` | `False` | Enable DFT-D3 dispersion correction |
| `d3_damping` | `str` | `"bj"` | D3 damping function: `"bj"`, `"zero"`, `"zerom"`, `"bjm"` |
| `d3_xc` | `str` | `"pbe"` | D3 exchange-correlation functional |
| `default_dtype` | `str` | `"float64"` | Default data type: `"float32"` or `"float64"` |

**Available Models:**

| Model | Description | Speed | Accuracy |
|-------|-------------|-------|----------|
| `"small"` | MACE-MP-0 small | Fast | Lower |
| `"medium"` | MACE-MP-0 medium | Balanced | **Recommended** |
| `"large"` | MACE-MP-0 large | Slow | Highest |
| `/path/to/model.pt` | Custom fine-tuned model | Varies | Varies |

**Example:**

```python
# Basic usage
calc = MACEInference(model="medium", device="auto")

# With D3 dispersion correction (for MOFs, layered materials)
calc = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")

# GPU with custom model
calc = MACEInference(model="/path/to/custom.pt", device="cuda")

# Fast inference with lower precision
calc = MACEInference(model="small", default_dtype="float32")
```

---

### Single-Point Calculations

#### `single_point()`

Calculate energy, forces, and stress for a structure.

```python
def single_point(
    self,
    structure: Union[str, Path, Atoms],
    properties: Optional[List[str]] = None
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | `str`, `Path`, `Atoms` | Required | Structure file path (CIF, XYZ, POSCAR) or ASE Atoms object |
| `properties` | `List[str]` | `None` | Properties to calculate. Default: `["energy", "forces", "stress"]` |

**Returns:** `Dict[str, Any]`

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `energy` | `float` | eV | Total potential energy |
| `forces` | `np.ndarray` | eV/Å | Forces on atoms, shape `(N, 3)` |
| `stress` | `np.ndarray` | eV/Å³ | Stress tensor (Voigt notation), shape `(6,)` |
| `energy_per_atom` | `float` | eV | Energy per atom |

**Example:**

```python
result = calc.single_point("structure.cif")
print(f"Energy: {result['energy']:.4f} eV")
print(f"Energy/atom: {result['energy_per_atom']:.4f} eV")
print(f"Max force: {np.max(np.abs(result['forces'])):.4f} eV/Å")
```

---

### Structure Optimization

#### `optimize()`

Minimize energy by relaxing atomic positions and optionally cell parameters.

```python
def optimize(
    self,
    structure: Union[str, Path, Atoms],
    fmax: float = 0.05,
    steps: int = 500,
    optimizer: str = "LBFGS",
    optimize_cell: bool = False,
    trajectory: Optional[str] = None,
    logfile: Optional[str] = None,
    output: Optional[str] = None
) -> Atoms
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | `str`, `Path`, `Atoms` | Required | Input structure |
| `fmax` | `float` | `0.05` | Force convergence criterion (eV/Å) |
| `steps` | `int` | `500` | Maximum optimization steps |
| `optimizer` | `str` | `"LBFGS"` | Algorithm: `"LBFGS"`, `"BFGS"`, `"FIRE"` |
| `optimize_cell` | `bool` | `False` | Also optimize cell parameters |
| `trajectory` | `str` | `None` | Path to save trajectory file |
| `logfile` | `str` | `None` | Path to save optimization log |
| `output` | `str` | `None` | Path to save optimized structure |

**Returns:** `Atoms` - Optimized ASE Atoms object

**Example:**

```python
# Basic optimization
optimized = calc.optimize("structure.cif", fmax=0.05)

# Full cell relaxation with output
optimized = calc.optimize(
    "structure.cif",
    fmax=0.01,
    optimize_cell=True,
    trajectory="opt.traj",
    output="optimized.cif"
)
```

---

### Molecular Dynamics

#### `run_md()`

Run molecular dynamics simulation with NVT or NPT ensemble.

```python
def run_md(
    self,
    structure: Union[str, Path, Atoms],
    ensemble: Literal["nvt", "npt"] = "nvt",
    temperature_K: float = 300,
    steps: int = 1000,
    timestep: float = 1.0,
    pressure_GPa: Optional[float] = None,
    taut: Optional[float] = None,
    taup: Optional[float] = None,
    trajectory: Optional[str] = None,
    logfile: Optional[str] = None,
    log_interval: int = 100,
    progress_callback: Optional[Callable[[int, int], None]] = None
) -> Atoms
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | `str`, `Path`, `Atoms` | Required | Input structure |
| `ensemble` | `str` | `"nvt"` | Ensemble type: `"nvt"` or `"npt"` |
| `temperature_K` | `float` | `300` | Target temperature (K) |
| `steps` | `int` | `1000` | Number of MD steps |
| `timestep` | `float` | `1.0` | Time step (fs) |
| `pressure_GPa` | `float` | `None` | Target pressure for NPT (GPa) |
| `taut` | `float` | `None` | Temperature coupling time (fs) |
| `taup` | `float` | `None` | Pressure coupling time (fs) |
| `trajectory` | `str` | `None` | Path to save trajectory |
| `logfile` | `str` | `None` | Path to save MD log |
| `log_interval` | `int` | `100` | Logging interval (steps) |
| `progress_callback` | `Callable` | `None` | Progress callback `(current_step, total_steps)` |

**Returns:** `Atoms` - Final structure after MD

**Example:**

```python
# NVT simulation
final = calc.run_md(
    "structure.cif",
    ensemble="nvt",
    temperature_K=300,
    steps=10000,
    trajectory="nvt.traj"
)

# NPT simulation with progress callback
def progress(step, total):
    print(f"\rMD: {step}/{total} ({100*step/total:.1f}%)", end="")

final = calc.run_md(
    "structure.cif",
    ensemble="npt",
    temperature_K=300,
    pressure_GPa=0.1,
    steps=5000,
    progress_callback=progress
)
```

---

### Phonon Calculations

#### `phonon()`

Calculate phonon properties and optionally thermal properties.

```python
def phonon(
    self,
    structure: Union[str, Path, Atoms],
    supercell_matrix: Union[List[int], int] = 2,
    displacement: float = 0.01,
    mesh: List[int] = None,
    temperature_range: Optional[Tuple[float, float, float]] = None,
    output_dir: Optional[str] = None
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | `str`, `Path`, `Atoms` | Required | Input structure (primitive cell) |
| `supercell_matrix` | `int` or `List[int]` | `2` | Supercell size (e.g., `2` or `[2,2,2]`) |
| `displacement` | `float` | `0.01` | Atomic displacement distance (Å) |
| `mesh` | `List[int]` | `[20,20,20]` | k-point mesh for phonon DOS |
| `temperature_range` | `Tuple` | `None` | `(T_min, T_max, T_step)` in K for thermal properties |
| `output_dir` | `str` | `None` | Directory for output files |

**Returns:** `Dict[str, Any]`

| Key | Type | Description |
|-----|------|-------------|
| `phonon` | `Phonopy` | Phonopy object with force constants |
| `frequencies` | `np.ndarray` | Phonon frequencies at Gamma point (THz) |
| `thermal` | `Dict` | Thermal properties (if `temperature_range` specified) |

**Example:**

```python
# Basic phonon calculation
result = calc.phonon("structure.cif", supercell_matrix=[2, 2, 2])
print(f"Frequencies at Gamma: {result['frequencies']}")

# With thermal properties
result = calc.phonon(
    "structure.cif",
    supercell_matrix=2,
    temperature_range=(0, 1000, 10)
)
thermal = result['thermal']
print(f"Free energy at 300K: {thermal['free_energy'][30]:.2f} kJ/mol")
```

---

### Mechanical Properties

#### `bulk_modulus()`

Calculate bulk modulus via equation of state fitting.

```python
def bulk_modulus(
    self,
    structure: Union[str, Path, Atoms],
    strain_range: float = 0.05,
    n_points: int = 7,
    optimize_first: bool = True,
    fmax: float = 0.01
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | `str`, `Path`, `Atoms` | Required | Input structure |
| `strain_range` | `float` | `0.05` | Maximum strain for volume sweep (±5%) |
| `n_points` | `int` | `7` | Number of strain points |
| `optimize_first` | `bool` | `True` | Optimize structure before calculation |
| `fmax` | `float` | `0.01` | Force criterion for optimization |

**Returns:** `Dict[str, Any]`

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `v0` | `float` | Å³ | Equilibrium volume |
| `e0` | `float` | eV | Equilibrium energy |
| `B` | `float` | eV/Å³ | Bulk modulus (raw) |
| `B_GPa` | `float` | GPa | Bulk modulus |
| `volumes` | `np.ndarray` | Å³ | Sampled volumes |
| `energies` | `np.ndarray` | eV | Corresponding energies |

**Example:**

```python
result = calc.bulk_modulus("structure.cif")
print(f"Bulk modulus: {result['B_GPa']:.1f} GPa")
print(f"Equilibrium volume: {result['v0']:.2f} Å³")
```

---

### Adsorption Studies

#### `adsorption_energy()`

Calculate gas molecule adsorption energy in porous materials.

```python
def adsorption_energy(
    self,
    framework: Union[str, Path, Atoms],
    adsorbate: Union[str, Atoms],
    site_position: List[float],
    optimize: bool = True,
    fmax: float = 0.05,
    fix_framework: bool = True
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `framework` | `str`, `Path`, `Atoms` | Required | Framework structure (MOF, zeolite) |
| `adsorbate` | `str` or `Atoms` | Required | Adsorbate molecule (`"CO2"`, `"H2O"`, etc.) or Atoms |
| `site_position` | `List[float]` | Required | Adsorption site `[x, y, z]` in Å |
| `optimize` | `bool` | `True` | Optimize combined structure |
| `fmax` | `float` | `0.05` | Force criterion for optimization |
| `fix_framework` | `bool` | `True` | Fix framework atoms during optimization |

**Supported Adsorbates:** `"H2"`, `"N2"`, `"O2"`, `"CO"`, `"CO2"`, `"H2O"`, `"CH4"`, `"NH3"`, etc.

**Returns:** `Dict[str, Any]`

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `adsorption_energy` | `float` | eV | Adsorption energy (negative = favorable) |
| `E_framework` | `float` | eV | Framework energy |
| `E_adsorbate` | `float` | eV | Isolated adsorbate energy |
| `E_complex` | `float` | eV | Combined system energy |

**Example:**

```python
result = calc.adsorption_energy(
    framework="mof.cif",
    adsorbate="CO2",
    site_position=[10.0, 10.0, 10.0]
)
print(f"Adsorption energy: {result['adsorption_energy']:.3f} eV")
# Negative value indicates favorable adsorption
```

#### `coordination()`

Analyze metal coordination environment.

```python
def coordination(
    self,
    structure: Union[str, Path, Atoms],
    metal_indices: Optional[List[int]] = None,
    cutoff_multiplier: float = 1.3
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structure` | `str`, `Path`, `Atoms` | Required | Input structure |
| `metal_indices` | `List[int]` | `None` | Metal atom indices (auto-detect if None) |
| `cutoff_multiplier` | `float` | `1.3` | Multiplier for covalent radii cutoff |

**Returns:** `Dict[str, Any]`

| Key | Type | Description |
|-----|------|-------------|
| `coordination` | `Dict` | Coordination info per metal atom |
| `metal_indices` | `List[int]` | Detected/specified metal indices |

**Example:**

```python
result = calc.coordination("mof.cif")
for idx, info in result["coordination"].items():
    print(f"Metal {idx}: CN={info['coordination_number']}, "
          f"avg distance={info['average_distance']:.3f} Å")
```

---

### Batch Processing

#### `batch_single_point()`

Perform single-point calculations on multiple structures efficiently.

```python
def batch_single_point(
    self,
    structures: List[Union[str, Path, Atoms]],
    properties: Optional[List[str]] = None,
    progress_callback: Optional[Callable[[int, int], None]] = None
) -> List[Dict[str, Any]]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structures` | `List` | Required | List of structure files or Atoms objects |
| `properties` | `List[str]` | `None` | Properties to calculate |
| `progress_callback` | `Callable` | `None` | Callback `(current, total)` for progress |

**Returns:** `List[Dict[str, Any]]` - Results for each structure

Each result dict includes:
- All standard single-point properties
- `success`: `bool` - Whether calculation succeeded
- `structure_index`: `int` - Index in input list
- `error`: `str` - Error message (if failed)

**Example:**

```python
structures = ["s1.cif", "s2.cif", "s3.cif"]

def progress(current, total):
    print(f"Processing {current}/{total}")

results = calc.batch_single_point(structures, progress_callback=progress)

for r in results:
    if r["success"]:
        print(f"Structure {r['structure_index']}: E = {r['energy']:.4f} eV")
    else:
        print(f"Structure {r['structure_index']}: FAILED - {r['error']}")
```

#### `batch_optimize()`

Optimize multiple structures in batch.

```python
def batch_optimize(
    self,
    structures: List[Union[str, Path, Atoms]],
    fmax: float = 0.05,
    steps: int = 500,
    optimizer: str = "LBFGS",
    optimize_cell: bool = False,
    output_dir: Optional[Union[str, Path]] = None,
    progress_callback: Optional[Callable[[int, int], None]] = None
) -> List[Dict[str, Any]]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `structures` | `List` | Required | List of input structures |
| `fmax` | `float` | `0.05` | Force convergence criterion |
| `steps` | `int` | `500` | Max steps per structure |
| `optimizer` | `str` | `"LBFGS"` | Optimization algorithm |
| `optimize_cell` | `bool` | `False` | Optimize cell parameters |
| `output_dir` | `str`, `Path` | `None` | Directory to save optimized structures |
| `progress_callback` | `Callable` | `None` | Progress callback |

**Returns:** `List[Dict[str, Any]]`

Each result dict includes:
- `success`: `bool` - Whether optimization succeeded
- `structure_index`: `int` - Index in input list
- `atoms`: `Atoms` - Optimized structure (or None if failed)
- `final_energy`: `float` - Final energy
- `final_max_force`: `float` - Final maximum force
- `converged`: `bool` - Whether optimization converged
- `output_file`: `str` - Path to saved file (if `output_dir` specified)

**Example:**

```python
results = calc.batch_optimize(
    structures,
    fmax=0.05,
    output_dir="optimized/",
    progress_callback=lambda c, t: print(f"{c}/{t}")
)

for r in results:
    if r["converged"]:
        print(f"Structure {r['structure_index']}: converged, E = {r['final_energy']:.4f} eV")
```

---

## Task Functions

Low-level task functions are available for advanced usage:

```python
from mace_inference.tasks import (
    single_point_energy,
    optimize_structure,
    run_nvt_md,
    run_npt_md,
    calculate_phonon,
    calculate_thermal_properties,
    calculate_bulk_modulus,
    calculate_elastic_constants,
    calculate_adsorption_energy,
    analyze_coordination,
    find_adsorption_sites,
)
```

### `calculate_elastic_constants()`

Calculate full elastic tensor using stress-strain analysis.

```python
def calculate_elastic_constants(
    atoms: Atoms,
    calculator: Calculator,
    strain_magnitude: float = 0.01,
    optimize_first: bool = True,
    fmax: float = 0.01
) -> Dict[str, Any]
```

**Returns:**

| Key | Type | Description |
|-----|------|-------------|
| `C` | `np.ndarray` | 6×6 elastic tensor (GPa) |
| `K_V`, `K_R`, `K_H` | `float` | Voigt/Reuss/Hill bulk modulus (GPa) |
| `G_V`, `G_R`, `G_H` | `float` | Voigt/Reuss/Hill shear modulus (GPa) |
| `E_H` | `float` | Hill Young's modulus (GPa) |
| `nu_H` | `float` | Hill Poisson's ratio |
| `C11`, `C12`, `C44` | `float` | Cubic elastic constants (GPa, if applicable) |
| `A` | `float` | Zener anisotropy ratio (if cubic) |

### `find_adsorption_sites()`

Find potential adsorption sites in porous structures.

```python
def find_adsorption_sites(
    atoms: Atoms,
    grid_spacing: float = 0.5,
    probe_radius: float = 1.2,
    min_distance: float = 2.0,
    energy_cutoff: float = None
) -> List[np.ndarray]
```

**Returns:** List of `[x, y, z]` positions for potential adsorption sites.

---

## Type Definitions

Type hints are provided via `TypedDict` classes:

```python
from mace_inference import (
    SinglePointResult,
    BulkModulusResult,
    PhononResult,
    AdsorptionResult,
    CoordinationResult,
    ThermalPropertiesResult,
    Calculator,  # Protocol for calculator interface
)
```

---

## Utility Functions

```python
from mace_inference import get_device
from mace_inference.utils import (
    parse_structure_input,
    save_structure,
    create_supercell,
    atoms_to_dict,
    get_device_info,
    check_d3_available,
)
```

### `get_device()`

```python
def get_device(device: str = "auto") -> str
```

Detect and return appropriate compute device.

**Example:**

```python
from mace_inference import get_device

device = get_device("auto")  # Returns "cuda" if available, else "cpu"
print(f"Using device: {device}")
```

---

## CLI Reference

The `mace-infer` command provides a command-line interface.

```bash
mace-infer --help
mace-infer --version
```

### Commands

| Command | Description |
|---------|-------------|
| `energy` | Single-point energy calculation |
| `optimize` | Structure optimization |
| `md` | Molecular dynamics simulation |
| `phonon` | Phonon calculation |
| `bulk-modulus` | Bulk modulus calculation |
| `adsorption` | Adsorption energy calculation |
| `info` | Display environment information |

### Examples

```bash
# Single-point energy
mace-infer energy structure.cif --model medium

# Optimization
mace-infer optimize structure.cif --fmax 0.05 --output optimized.cif

# MD simulation
mace-infer md structure.cif --ensemble nvt --temp 300 --steps 10000

# Phonon calculation
mace-infer phonon structure.cif --supercell 2 2 2

# With D3 correction
mace-infer adsorption mof.cif --gas CO2 --site 10 10 10 --d3

# System info
mace-infer info --verbose
```

---

## See Also

- [README.md](README.md) - Project overview and quick start
- [QUICKSTART.md](QUICKSTART.md) - Quick start guide
- [INSTALL_GUIDE.md](INSTALL_GUIDE.md) - Detailed installation instructions
- [CHANGELOG.md](CHANGELOG.md) - Version history
- [examples/](examples/) - Example scripts

---

*Documentation generated: 2026-01-08*
