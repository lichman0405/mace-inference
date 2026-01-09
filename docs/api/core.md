# MACEInference

::: mace_inference.MACEInference
    options:
      show_root_heading: true
      show_source: true
      heading_level: 2
      members:
        - __init__
        - single_point
        - optimize
        - run_md
        - phonon
        - bulk_modulus
        - adsorption_energy
        - coordination
        - get_calculator

## Constructor

```python
MACEInference(
    model: str = "medium",
    device: str = "auto",
    enable_d3: bool = False,
    d3_xc: str = "pbe",
    d3_damping: str = "bj",
    default_dtype: str = "float64"
)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `model` | `str` | `"medium"` | Model name (`"small"`, `"medium"`, `"large"`) or path to `.model` file |
| `device` | `str` | `"auto"` | Device: `"auto"`, `"cpu"`, `"cuda"`, `"cuda:0"`, `"mps"` |
| `enable_d3` | `bool` | `False` | Enable DFT-D3 dispersion correction |
| `d3_xc` | `str` | `"pbe"` | Exchange-correlation functional for D3 |
| `d3_damping` | `str` | `"bj"` | D3 damping function (`"bj"` or `"zero"`) |
| `default_dtype` | `str` | `"float64"` | Default precision (`"float32"` or `"float64"`) |

### Example

```python
from mace_inference import MACEInference

# Basic initialization
calc = MACEInference(model="medium", device="auto")

# With D3 correction
calc = MACEInference(
    model="medium",
    device="cuda",
    enable_d3=True,
    d3_xc="pbe"
)

# Custom model
calc = MACEInference(model="/path/to/model.model")
```

## Methods

### single_point

```python
def single_point(self, atoms: Atoms) -> dict
```

Calculate energy, forces, and stress for a structure.

**Parameters:**

- `atoms` (`Atoms`): ASE Atoms object

**Returns:** `dict` with keys:

- `energy` (`float`): Total energy in eV
- `energy_per_atom` (`float`): Energy per atom in eV
- `forces` (`ndarray`): Forces array (N, 3) in eV/Å
- `max_force` (`float`): Maximum force magnitude in eV/Å
- `stress` (`ndarray`): Stress tensor in Voigt notation (eV/Å³)

---

### optimize

```python
def optimize(
    self,
    atoms: Atoms,
    fmax: float = 0.05,
    steps: int = 500,
    optimizer: str = "LBFGS",
    optimize_cell: bool = False,
    trajectory: Optional[str] = None,
    logfile: Optional[str] = None,
    output: Optional[str] = None
) -> Atoms
```

Optimize atomic positions and optionally cell parameters.

**Parameters:**

- `atoms` (`Atoms`): Structure to optimize
- `fmax` (`float`): Force convergence criterion in eV/Å
- `steps` (`int`): Maximum optimization steps
- `optimizer` (`str`): `"LBFGS"`, `"BFGS"` or `"FIRE"`
- `optimize_cell` (`bool`): Also optimize cell parameters
- `trajectory` (`str`): Save trajectory to this file
- `logfile` (`str`): Optimization log file
- `output` (`str`): Save optimized structure to this file

**Returns:** `Atoms` - Optimized structure

---

### run_md

```python
def run_md(
    self,
    atoms: Atoms,
    ensemble: str = "nvt",
    temperature_K: float = 300,
    steps: int = 1000,
    timestep: float = 1.0,
    pressure_GPa: Optional[float] = None,
    taut: Optional[float] = None,
    taup: Optional[float] = None,
    trajectory: Optional[str] = None,
    logfile: Optional[str] = None,
    log_interval: int = 100,
    progress_callback: Optional[Callable] = None
) -> Atoms
```

Run molecular dynamics simulation.

**Parameters:**

- `atoms` (`Atoms`): Initial structure
- `ensemble` (`str`): `"nvt"` or `"npt"`
- `temperature_K` (`float`): Temperature in Kelvin
- `steps` (`int`): Number of MD steps
- `timestep` (`float`): Time step in fs
- `pressure_GPa` (`float`): Pressure in GPa (NPT only)
- `taut` (`float`): Temperature coupling time in fs
- `taup` (`float`): Pressure coupling time in fs
- `trajectory` (`str`): Save trajectory to this file
- `logfile` (`str`): MD log file
- `log_interval` (`int`): Logging interval in steps
- `progress_callback` (`Callable`): Callback function(current_step, total_steps)

**Returns:** `Atoms` - Final structure after MD simulation

---

### phonon

```python
def phonon(
    self,
    atoms: Atoms,
    supercell_matrix: Union[List[int], int] = 2,
    displacement: float = 0.01,
    mesh: List[int] = [20, 20, 20],
    temperature_range: Optional[tuple] = None,
    output_dir: Optional[str] = None
) -> dict
```

Calculate phonon properties.

**Parameters:**

- `atoms` (`Atoms`): Unit cell structure
- `supercell_matrix` (`int` or `list`): Supercell dimensions (e.g., 2 or [2, 2, 2])
- `displacement` (`float`): Displacement magnitude in Å
- `mesh` (`list`): q-point mesh for DOS
- `temperature_range` (`tuple`): Temperature range for thermal properties (min, max, step)
- `output_dir` (`str`): Directory for output files

**Returns:** `dict` with keys:

- `frequencies` (`ndarray`): Frequencies at Γ in THz
- `phonon` (`Phonopy`): Phonopy object
- `thermal` (`dict`): Thermal properties (if temperature_range provided)

---

### bulk_modulus

```python
def bulk_modulus(
    self,
    atoms: Atoms,
    scale_range: Tuple[float, float] = (0.95, 1.05),
    n_points: int = 7,
    eos: str = "birchmurnaghan"
) -> dict
```

Calculate bulk modulus by fitting equation of state.

**Parameters:**

- `atoms` (`Atoms`): Structure (should be optimized)
- `scale_range` (`tuple`): Volume scaling range
- `n_points` (`int`): Number of points for fit
- `eos` (`str`): Equation of state type

**Returns:** `dict` with keys:

- `bulk_modulus` (`float`): Bulk modulus in GPa
- `equilibrium_volume` (`float`): Equilibrium volume in ų
- `equilibrium_energy` (`float`): Energy at equilibrium in eV
- `volumes` (`list`): Volumes used in fit
- `energies` (`list`): Energies at each volume

---

### adsorption_energy

```python
def adsorption_energy(
    self,
    framework: Atoms,
    adsorbate: Atoms,
    site_position: List[float],
    optimize: bool = True,
    fmax: float = 0.05,
    fix_framework: bool = True
) -> dict
```

Calculate gas adsorption energy in a framework.

**Parameters:**

- `framework` (`Atoms`): Framework structure (MOF, zeolite, etc.)
- `adsorbate` (`Atoms`): Gas molecule
- `site_position` (`list`): Position [x, y, z] to place adsorbate
- `optimize` (`bool`): Optimize the complex
- `fmax` (`float`): Force tolerance in eV/Å
- `fix_framework` (`bool`): Keep framework fixed

**Returns:** `dict` with keys:

- `E_ads` (`float`): Adsorption energy in eV
- `E_mof` (`float`): Framework energy in eV
- `E_gas` (`float`): Gas molecule energy in eV
- `E_complex` (`float`): Complex energy in eV
- `complex_structure` (`Atoms`): Final complex structure
- `optimized` (`bool`): Whether structure was optimized

---

### coordination

```python
def coordination(
    self,
    atoms: Atoms,
    metal_indices: Optional[List[int]] = None,
    cutoff_multiplier: float = 1.3
) -> dict
```

Analyze coordination environment of metal atoms.

**Parameters:**

- `atoms` (`Atoms`): Structure to analyze
- `metal_indices` (`list`): Indices of metal atoms (auto-detect if None)
- `cutoff_multiplier` (`float`): Multiplier for covalent radii cutoff

**Returns:** `dict` with keys:

- `coordination` (`dict`): Coordination info for each metal atom
- `metal_indices` (`list`): Indices of metal atoms analyzed

---

### get_calculator

```python
def get_calculator(self) -> Calculator
```

Get the underlying ASE calculator.

**Returns:** ASE `Calculator` object (MACE or MACE+D3)
