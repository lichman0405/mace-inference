# Tasks Module

The tasks module contains the implementation of various simulation tasks.

## Static Calculations

::: mace_inference.tasks.static
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - calculate_single_point

### calculate_single_point

```python
def calculate_single_point(atoms: Atoms, calculator: Calculator) -> dict
```

Calculate energy, forces, and stress for a structure.

**Parameters:**

- `atoms` (`Atoms`): ASE Atoms object
- `calculator` (`Calculator`): ASE calculator

**Returns:** `dict` with energy, forces, stress, etc.

---

## Optimization

::: mace_inference.tasks.static
    options:
      show_root_heading: false
      heading_level: 3
      members:
        - optimize_structure

### optimize_structure

```python
def optimize_structure(
    atoms: Atoms,
    calculator: Calculator,
    fmax: float = 0.05,
    steps: int = 500,
    optimizer: str = "BFGS",
    relax_cell: bool = False,
    fix_symmetry: bool = False,
    trajectory: Optional[str] = None
) -> dict
```

Optimize atomic positions and cell.

---

## Molecular Dynamics

::: mace_inference.tasks.dynamics
    options:
      show_root_heading: false
      heading_level: 3

### run_nvt

```python
def run_nvt(
    atoms: Atoms,
    calculator: Calculator,
    temperature: float,
    timestep: float = 1.0,
    steps: int = 1000,
    friction: float = 0.01,
    save_interval: int = 10,
    log_interval: int = 100
) -> dict
```

Run NVT molecular dynamics with Langevin thermostat.

### run_npt

```python
def run_npt(
    atoms: Atoms,
    calculator: Calculator,
    temperature: float,
    pressure: float,
    timestep: float = 1.0,
    steps: int = 1000,
    save_interval: int = 10,
    log_interval: int = 100
) -> dict
```

Run NPT molecular dynamics with Berendsen barostat.

---

## Phonon Calculations

::: mace_inference.tasks.phonon
    options:
      show_root_heading: false
      heading_level: 3

### calculate_phonon

```python
def calculate_phonon(
    atoms: Atoms,
    calculator: Calculator,
    supercell: Tuple[int, int, int] = (2, 2, 2),
    delta: float = 0.01,
    mesh: Tuple[int, int, int] = (20, 20, 20)
) -> dict
```

Calculate phonon frequencies and density of states.

### calculate_thermal_properties

```python
def calculate_thermal_properties(
    phonon: Phonopy,
    t_min: float = 0,
    t_max: float = 1000,
    t_step: float = 10
) -> dict
```

Calculate thermal properties from phonon data.

---

## Mechanical Properties

::: mace_inference.tasks.mechanics
    options:
      show_root_heading: false
      heading_level: 3

### calculate_bulk_modulus

```python
def calculate_bulk_modulus(
    atoms: Atoms,
    calculator: Calculator,
    scale_range: Tuple[float, float] = (0.95, 1.05),
    n_points: int = 7,
    eos: str = "birchmurnaghan"
) -> dict
```

Calculate bulk modulus by fitting equation of state.

---

## Adsorption

::: mace_inference.tasks.adsorption
    options:
      show_root_heading: false
      heading_level: 3

### calculate_adsorption_energy

```python
def calculate_adsorption_energy(
    framework: Atoms,
    adsorbate: Atoms,
    site_position: List[float],
    calculator: Calculator,
    optimize: bool = True,
    fmax: float = 0.05,
    fix_framework: bool = True
) -> dict
```

Calculate adsorption energy of a gas molecule in a framework.

**The adsorption energy is calculated as:**

$$E_{ads} = E_{complex} - E_{framework} - E_{gas}$$
