# MACE Inference Examples

This directory contains example scripts demonstrating the various capabilities of the `mace-inference` library.

## Prerequisites

Before running the examples, ensure you have installed the library:

```bash
# Basic installation
pip install mace-inference

# With D3 dispersion correction support
pip install mace-inference[d3]

# With GPU support
pip install mace-inference[gpu]
```

## Examples Overview

| Example | Description | Key Features |
|---------|-------------|--------------|
| [01_basic_usage.py](01_basic_usage.py) | Getting started | Single-point energy, optimization |
| [02_molecular_dynamics.py](02_molecular_dynamics.py) | MD simulations | NVT, NPT ensembles |
| [03_phonon_calculation.py](03_phonon_calculation.py) | Phonon properties | Dispersion, thermal properties |
| [04_adsorption_study.py](04_adsorption_study.py) | Gas adsorption | MOFs, binding energies |
| [05_high_throughput.py](05_high_throughput.py) | Batch processing | Screening workflows |
| [06_d3_correction.py](06_d3_correction.py) | D3 correction | Dispersion interactions |
| [07_batch_processing.py](07_batch_processing.py) | Batch API & Progress | Callbacks, error handling |

---

## 01. Basic Usage

**File:** `01_basic_usage.py`

Learn the fundamentals of using MACE Inference:
- Creating a `MACEInference` calculator
- Loading and building structures with ASE
- Single-point energy calculations
- Structure optimization

```python
from mace_inference import MACEInference

calc = MACEInference(model="medium", device="auto")
result = calc.single_point(atoms)
optimized = calc.optimize(atoms, fmax=0.05)
```

---

## 02. Molecular Dynamics

**File:** `02_molecular_dynamics.py`

Run molecular dynamics simulations:
- NVT ensemble (constant temperature)
- NPT ensemble (constant temperature and pressure)
- Trajectory output and analysis

```python
final_atoms = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,
    steps=10000,
    trajectory="md.traj"
)
```

---

## 03. Phonon Calculation

**File:** `03_phonon_calculation.py`

Calculate phonon properties and thermodynamics:
- Phonon dispersion curves
- Density of states (DOS)
- Thermal properties (free energy, entropy, heat capacity)

```python
result = calc.phonon(
    atoms,
    supercell_matrix=[2, 2, 2],
    temperature_range=(0, 1000, 10)
)
```

---

## 04. Adsorption Study

**File:** `04_adsorption_study.py`

Study gas adsorption in porous materials:
- Adsorption energy calculation
- Multiple gas molecules (CO2, H2O, CH4, N2)
- D3 correction for van der Waals interactions

```python
calc = MACEInference(model="medium", enable_d3=True)
result = calc.adsorption_energy(
    framework=mof,
    adsorbate="CO2",
    site_position=[10.0, 10.0, 10.0]
)
```

---

## 05. High-Throughput Screening

**File:** `05_high_throughput.py`

Batch processing multiple structures:
- Efficient calculator reuse
- Parallel-friendly workflow
- Results aggregation and export

```python
for name, atoms in structures.items():
    result = calc.single_point(atoms)
    optimized = calc.optimize(atoms)
    bm = calc.bulk_modulus(optimized)
```

---

## 06. D3 Dispersion Correction

**File:** `06_d3_correction.py`

Using DFT-D3 dispersion correction:
- When to use D3 (layered materials, molecular crystals)
- Comparing results with/without D3
- Different damping functions (bj, zero, zerom, bjm)
- Different XC functional parameters

```python
calc = MACEInference(
    model="medium",
    enable_d3=True,
    d3_damping="bj",
    d3_xc="pbe"
)
```

---

## 07. Batch Processing and Progress Callbacks

**File:** `07_batch_processing.py`

Efficient batch processing with progress tracking:
- Batch single-point calculations on multiple structures
- Batch structure optimization
- Progress callbacks for long-running tasks
- Error handling in batch operations

```python
# Batch calculations with progress callback
def progress(current, total):
    print(f"Processing {current}/{total}")

results = calc.batch_single_point(
    structures,
    progress_callback=progress
)

# Batch optimization with output
opt_results = calc.batch_optimize(
    structures,
    fmax=0.05,
    output_dir="optimized/"
)

# MD with progress callback
final = calc.run_md(
    atoms,
    steps=10000,
    progress_callback=progress
)
```

---

## Running Examples

```bash
# Run a specific example
python examples/01_basic_usage.py

# Run with GPU acceleration
CUDA_VISIBLE_DEVICES=0 python examples/02_molecular_dynamics.py
```

## Output Files

Examples may generate the following output files:
- `*.cif`, `*.xyz` - Structure files
- `*.traj` - ASE trajectory files
- `*.log` - Optimization/MD log files
- `phonopy.yaml` - Phonopy output files

## Tips

1. **First run**: The first run downloads MACE models (~100-500 MB) and may take longer.

2. **Device selection**: Use `device="auto"` for automatic GPU detection, or specify `device="cuda"` or `device="cpu"`.

3. **Memory**: For large systems, consider using the "small" model or running on CPU with more RAM.

4. **D3 correction**: Enable for systems with significant dispersion interactions (layered materials, molecular crystals, gas adsorption).

## Need Help?

- [Documentation](https://github.com/lichman0405/mace-inference#readme)
- [Issue Tracker](https://github.com/lichman0405/mace-inference/issues)
