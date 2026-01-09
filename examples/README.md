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

## Structure Files

All examples use CIF structure files from the `structures/` directory:

| File | Description |
|------|-------------|
| `cu_fcc.cif` | Copper FCC crystal (4 atoms) |
| `si_diamond.cif` | Silicon diamond structure (8 atoms) |
| `cu_paddlewheel.cif` | Cu-paddlewheel cluster (MOF building block) |

You can add your own CIF files to this directory for testing.

## Examples Overview

| Example | Description | Key Features |
|---------|-------------|--------------|
| [01_basic_usage.py](01_basic_usage.py) | Getting started | Single-point energy, optimization |
| [02_molecular_dynamics.py](02_molecular_dynamics.py) | MD simulations | NVT, NPT ensembles |
| [03_phonon_calculation.py](03_phonon_calculation.py) | Phonon properties | Dispersion, thermal properties |
| [04_adsorption_study.py](04_adsorption_study.py) | Gas adsorption | MOFs, binding energies |
| [05_high_throughput.py](05_high_throughput.py) | Multi-structure screening | Batch comparison |
| [06_d3_correction.py](06_d3_correction.py) | D3 correction | Dispersion interactions |
| [07_batch_processing.py](07_batch_processing.py) | Batch API & I/O | Export results, file utilities |

---

## 01. Basic Usage

**File:** `01_basic_usage.py`

Learn the fundamentals of using MACE Inference:
- Creating a `MACEInference` calculator
- Loading structures from CIF files
- Single-point energy calculations
- Structure optimization

```python
from pathlib import Path
from ase.io import read
from mace_inference import MACEInference

# Load structure from CIF
atoms = read("structures/cu_fcc.cif")

# Initialize calculator
calc = MACEInference(model="medium", device="auto")

# Single-point calculation
result = calc.single_point(atoms)
print(f"Energy: {result['energy']:.4f} eV")

# Optimization
optimized = calc.optimize(atoms, fmax=0.05)
```

**Run:**
```bash
python 01_basic_usage.py
```

---

## 02. Molecular Dynamics

**File:** `02_molecular_dynamics.py`

Run molecular dynamics simulations:
- NVT ensemble (constant temperature)
- NPT ensemble (constant temperature and pressure)
- Temperature analysis

```python
# NVT simulation
final_atoms = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,
    steps=100,
    timestep=1.0
)

# NPT simulation
final_atoms = calc.run_md(
    atoms,
    ensemble="npt",
    temperature_K=300,
    pressure_GPa=0.0001,
    steps=100,
    timestep=1.0
)
```

**Run:**
```bash
python 02_molecular_dynamics.py
```

---

## 03. Phonon Calculation

**File:** `03_phonon_calculation.py`

Calculate phonon properties and thermodynamics:
- Force constant calculation via finite displacement
- Thermal properties (free energy, entropy, heat capacity)

```python
result = calc.phonon(
    atoms,
    supercell_matrix=[2, 2, 2],
    displacement=0.01,
    mesh=[10, 10, 10],
    temperature_range=(0, 800, 50)
)

# Access thermal properties
thermal = result['thermal']
print(f"Cv at 300K: {thermal['heat_capacity'][6]:.2f} J/(mol·K)")
```

**Run:**
```bash
python 03_phonon_calculation.py
```

---

## 04. Adsorption Study

**File:** `04_adsorption_study.py`

Study gas adsorption in porous materials:
- Adsorption energy calculation
- Multiple gas molecules (CO2, H2O, CH4, N2)
- Coordination number analysis

```python
# Calculate adsorption energy
result = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    adsorption_site=[5.0, 5.0, 5.0],
    optimize=True,
    fix_framework=True
)
print(f"Adsorption energy: {result['adsorption_energy']:.4f} eV")

# Coordination analysis
coord = calc.coordination_analysis(framework, cutoff=3.0)
print(f"Average CN: {coord['average_coordination']:.2f}")
```

**Run:**
```bash
python 04_adsorption_study.py
```

---

## 05. High-Throughput Screening

**File:** `05_high_throughput.py`

Process multiple structures efficiently:
- Automatic CIF file discovery
- Property comparison across materials
- Results export to JSON

```python
# Process all CIF files
for cif_path in structures_dir.glob("*.cif"):
    atoms = read(str(cif_path))
    result = calc.single_point(atoms)
    
    # Collect and compare
    results.append({
        'formula': atoms.get_chemical_formula(),
        'energy_per_atom': result['energy'] / len(atoms)
    })
```

**Run:**
```bash
python 05_high_throughput.py
```

---

## 06. D3 Dispersion Correction

**File:** `06_d3_correction.py`

Use DFT-D3 dispersion correction:
- Enable D3 for van der Waals interactions
- Compare with and without D3
- Guidelines for when to use D3

```python
# Without D3
calc_no_d3 = MACEInference(model="medium", use_d3=False)

# With D3
calc_with_d3 = MACEInference(
    model="medium",
    use_d3=True,
    d3_xc="pbe",
    d3_cutoff=40.0
)

# Compare energies
e_mace = calc_no_d3.single_point(atoms)['energy']
e_d3 = calc_with_d3.single_point(atoms)['energy']
print(f"D3 correction: {e_d3 - e_mace:.4f} eV")
```

**Run:**
```bash
python 06_d3_correction.py
```

---

## 07. Batch Processing

**File:** `07_batch_processing.py`

Efficient batch processing utilities:
- Batch calculations with error handling
- Export to JSON and CSV
- Save annotated structures (XYZ with energies)
- Batch optimization

```python
# Process with error handling
for atoms in structures:
    try:
        result = calc.single_point(atoms)
        atoms.info['mace_energy'] = result['energy']
    except Exception as e:
        print(f"Failed: {e}")

# Export results
with open("results.json", "w") as f:
    json.dump(results, f)
```

**Run:**
```bash
python 07_batch_processing.py
```

---

## Running All Examples

To run all examples sequentially:

```bash
cd examples
python 01_basic_usage.py
python 02_molecular_dynamics.py
python 03_phonon_calculation.py
python 04_adsorption_study.py
python 05_high_throughput.py
python 06_d3_correction.py
python 07_batch_processing.py
```

## Output Files

Some examples generate output files:
- `screening_results.json` - Results from example 05
- `output/` directory - Results from example 07

These are automatically created and can be safely deleted.

## Tips

1. **First run**: The first run will download the MACE model (~400 MB) to `~/.cache/mace/`
2. **Device selection**: Use `device="auto"` to automatically select GPU if available
3. **Custom structures**: Replace CIF files in `structures/` with your own
4. **Memory**: For large systems, reduce supercell size in phonon calculations
