# MACE Inference - Complete Installation and Usage Guide

> **Related Documents:**
> - [API Reference](API_REFERENCE.md) - Complete API documentation
> - [README](README.md) - Project overview
> - [Quick Start](QUICKSTART.md) - Get started in 5 minutes

---

## 📋 Prerequisites

Before installing `mace-inference`, ensure you have a proper Python environment with PyTorch.

### System Requirements

| Component | CPU Version | GPU Version |
|-----------|-------------|-------------|
| Python | 3.9 - 3.12 | 3.9 - 3.12 |
| PyTorch | >= 2.0 (recommended 2.3.x) | >= 2.0 with CUDA |
| RAM | >= 8 GB | >= 8 GB |
| GPU | Not required | NVIDIA with >= 8 GB VRAM |

### Step 0: Set Up Python Environment

**Option A: Conda (Recommended)**

```bash
# CPU version
conda create -n mace-inference python=3.10 -y
conda activate mace-inference
pip install torch==2.3.1 --index-url https://download.pytorch.org/whl/cpu

# GPU version (CUDA 12.1)
conda create -n mace-inference-gpu python=3.10 -y
conda activate mace-inference-gpu
pip install torch==2.3.1 --index-url https://download.pytorch.org/whl/cu121
```

**Option B: venv**

```bash
python -m venv mace-env
# Windows
.\mace-env\Scripts\activate
# Linux/Mac
source mace-env/bin/activate

pip install torch==2.3.1
```

---

## 📦 Project Structure

```
mace-inference/
├── src/mace_inference/          # Core source code
│   ├── __init__.py              # Package initialization
│   ├── core.py                  # MACEInference main class
│   ├── cli.py                   # Command-line interface
│   ├── types.py                 # Type definitions (TypedDict, Protocol)
│   ├── py.typed                 # PEP 561 type marker
│   ├── tasks/                   # Task modules
│   │   ├── __init__.py
│   │   ├── static.py            # Single-point energy, structure optimization
│   │   ├── dynamics.py          # Molecular dynamics
│   │   ├── phonon.py            # Phonon calculations
│   │   ├── mechanics.py         # Mechanical properties (bulk modulus, elastic constants)
│   │   └── adsorption.py        # Adsorption energy, coordination, site finding
│   └── utils/                   # Utility modules
│       ├── __init__.py
│       ├── device.py            # Device management
│       ├── d3_correction.py     # D3 correction
│       └── io.py                # Structure I/O
├── examples/                    # Usage examples
│   ├── README.md                # Examples overview
│   ├── 01_basic_usage.py
│   ├── 02_molecular_dynamics.py
│   ├── 03_phonon_calculation.py
│   ├── 04_adsorption_study.py
│   ├── 05_high_throughput.py
│   ├── 06_d3_correction.py
│   └── 07_batch_processing.py
├── tests/                       # Unit tests
│   ├── __init__.py
│   ├── conftest.py              # Pytest fixtures
│   ├── test_install.py
│   ├── test_utils.py
│   ├── test_core.py
│   └── test_cli.py
├── pyproject.toml               # Project configuration
├── README.md                    # Project description
├── QUICKSTART.md                # Quick start guide
├── API_REFERENCE.md             # API documentation
├── CHANGELOG.md                 # Version history
├── LICENSE                      # MIT License
└── INSTALL_GUIDE.md             # This file
```

## 🚀 Installation Steps

### Method 1: From PyPI (Coming Soon)

```bash
# Basic installation
pip install mace-inference

# With D3 dispersion correction
pip install mace-inference[d3]

# With GPU support
pip install mace-inference[gpu]

# Full installation
pip install mace-inference[all]
```

### Method 2: Local Development Installation (Recommended for now)

```bash
# 1. Clone or navigate to project directory
cd mace-inference

# 2. Basic installation (CPU version)
pip install -e .

# 3. Or install full version (including D3, GPU, dev tools)
pip install -e ".[all]"
```

### Method 3: Modular Installation

```bash
# CPU version only
pip install -e .

# CPU + D3 correction
pip install -e ".[d3]"

# GPU support
pip install -e ".[gpu]"

# Development version (includes testing and code checking tools)
pip install -e ".[dev]"
```

### Post-Installation Verification

```bash
# Run installation test script
python test_install.py

# Test command-line tools
mace-infer --version
mace-infer info
```

## 📚 Python API Usage

### 1. Basic Usage

```python
from mace_inference import MACEInference

# Initialize calculator (auto-detect device)
calc = MACEInference(model="medium", device="auto")

# Single-point energy calculation
result = calc.single_point("structure.cif")
print(f"Energy: {result['energy']:.6f} eV")
print(f"Max Force: {result['max_force']:.6f} eV/Å")
```

### 2. Structure Optimization

```python
# Optimize atomic coordinates only
optimized = calc.optimize(
    "structure.cif",
    fmax=0.05,
    output="optimized.cif"
)

# Optimize both atoms and cell
optimized = calc.optimize(
    "structure.cif",
    fmax=0.05,
    optimize_cell=True,
    output="optimized.cif"
)
```

### 3. Molecular Dynamics

```python
# NVT simulation
final = calc.run_md(
    "structure.cif",
    ensemble="nvt",
    temperature_K=300,
    steps=10000,
    timestep=1.0,
    trajectory="md.traj"
)

# NPT simulation
final = calc.run_md(
    "structure.cif",
    ensemble="npt",
    temperature_K=300,
    pressure_GPa=0.0,
    steps=10000,
    trajectory="npt.traj"
)
```

### 4. Phonon Calculations

```python
# Phonon + thermodynamic properties
result = calc.phonon(
    "structure.cif",
    supercell_matrix=[2, 2, 2],
    temperature_range=(0, 1000, 10)
)

# Extract heat capacity data
thermal = result['thermal_properties']
print(f"Heat capacity at 300K: {thermal['heat_capacity'][30]:.3f} J/(mol·K)")
```

### 5. Mechanical Properties

```python
# Bulk modulus
bm_result = calc.bulk_modulus("structure.cif")
print(f"Bulk Modulus: {bm_result['B_GPa']:.2f} GPa")
```

### 6. Adsorption Energy Calculation

```python
# Gas adsorption
result = calc.adsorption_energy(
    framework="mof.cif",
    adsorbate="CO2",
    site_position=[10.0, 10.0, 10.0],
    optimize=True
)

print(f"Adsorption Energy: {result['adsorption_energy']:.3f} eV")
```

### 7. Coordination Environment Analysis

```python
# Analyze metal coordination
coord_result = calc.coordination("mof.cif")

for metal_idx, info in coord_result["coordination"].items():
    print(f"Metal {metal_idx}: CN = {info['coordination_number']}")
    print(f"  Average distance: {info['average_distance']:.3f} Å")
```

## 🖥️ Command-Line Tools Usage

### Basic Commands

```bash
# View system information
mace-infer info --verbose

# Single-point energy
mace-infer energy structure.cif --model medium

# Structure optimization
mace-infer optimize structure.cif \
    --fmax 0.05 \
    --cell \
    --output optimized.cif

# Molecular dynamics
mace-infer md structure.cif \
    --ensemble nvt \
    --temp 300 \
    --steps 10000 \
    --trajectory md.traj

# Phonon calculation
mace-infer phonon structure.cif \
    --supercell 2 2 2 \
    --temp-range 0 1000 10

# Bulk modulus
mace-infer bulk-modulus structure.cif

# Adsorption energy
mace-infer adsorption mof.cif \
    --gas CO2 \
    --site 10.0 10.0 10.0 \
    --d3
```

## 🔧 Advanced Configuration

### Using D3 Correction

```python
# Python API
calc = MACEInference(
    model="medium",
    device="cuda",
    enable_d3=True,
    d3_damping="bj",
    d3_xc="pbe"
)

# CLI
mace-infer adsorption mof.cif --gas CO2 --site 10 10 10 --d3
```

### GPU Acceleration

```python
# Force GPU usage
calc = MACEInference(model="large", device="cuda")

# CLI auto-detects
mace-infer energy structure.cif --device cuda
```

### Custom Models

```python
# Use custom trained model
calc = MACEInference(model="/path/to/custom_model.pt", device="cuda")
```

## 📊 Batch Processing Example

```python
from mace_inference import MACEInference
from pathlib import Path

calc = MACEInference(model="medium", device="auto")

# Batch single-point calculations with progress callback
def progress(current, total):
    print(f"Processing {current}/{total}")

structures = list(Path("structures/").glob("*.cif"))
results = calc.batch_single_point(structures, progress_callback=progress)

# Batch optimization with output directory
opt_results = calc.batch_optimize(
    structures,
    fmax=0.05,
    output_dir="optimized/",
    progress_callback=progress
)

# Process results
for result in results:
    if result.get("success", True):
        print(f"Structure {result['structure_index']}: E = {result['energy']:.4f} eV")
    else:
        print(f"Structure {result['structure_index']}: FAILED - {result.get('error')}")
```

### MD with Progress Callback

```python
# Monitor MD progress
def md_progress(current_step, total_steps):
    percent = 100 * current_step / total_steps
    print(f"\rMD: {current_step}/{total_steps} ({percent:.1f}%)", end="")

final = calc.run_md(
    atoms,
    ensemble="nvt",
    temperature_K=300,
    steps=10000,
    progress_callback=md_progress
)
```

## 🧪 Running Examples

```bash
cd examples

# Basic usage
python 01_basic_usage.py

# Molecular dynamics
python 02_molecular_dynamics.py

# Phonon calculation (takes a few minutes)
python 03_phonon_calculation.py

# Adsorption study
python 04_adsorption_study.py

# High-throughput screening
python 05_high_throughput.py

# D3 dispersion correction
python 06_d3_correction.py

# Batch processing with progress callbacks
python 07_batch_processing.py
```

## 🧰 Development and Testing

```bash
# Run unit tests
pytest

# Code coverage
pytest --cov=mace_inference --cov-report=html

# Code formatting
black src/

# Code linting
flake8 src/

# Type checking
mypy src/
```

## 📖 API Documentation

For detailed API documentation, refer to:
- [API Reference](API_REFERENCE.md) - Complete API documentation
- [Quick Start](QUICKSTART.md) - Getting started guide
- [Examples](examples/README.md) - Usage examples

## ⚠️ Troubleshooting

### 1. mace-torch Installation Fails

```bash
# Install PyTorch first
pip install torch==2.3.1

# Then install mace-torch
pip install mace-torch
```

### 2. CUDA Not Available

Check PyTorch CUDA installation:
```python
import torch
print(torch.cuda.is_available())
print(torch.version.cuda)
```

### 3. Out of Memory (MD/Phonon Calculations)

- Reduce supercell size
- Use CPU instead of GPU (larger memory)
- Reduce MD steps or phonon mesh density

### 4. CLI Command Not Found

```bash
# Ensure installation in editable mode
pip install -e .

# Or run directly
python -m mace_inference.cli --version
```

## 🚀 Performance Optimization Tips

1. **GPU vs CPU**
   - Small systems (<100 atoms): CPU and GPU comparable
   - Large systems (>500 atoms): GPU 10-100x faster

2. **Batch Processing**
   - Reuse same `MACEInference` instance
   - Avoid repeated model loading

3. **MD Simulations**
   - Minimum supercell: 3×3×3
   - Timestep: 1-2 fs
   - Equilibration: at least 10000 steps

4. **Phonon Calculations**
   - Supercell: 2×2×2 usually sufficient
   - Mesh: [20, 20, 20] is a good starting point

## 📝 Next Steps

- Read [README.md](README.md) for project overview
- Check [QUICKSTART.md](QUICKSTART.md) for quick start
- Run examples in [examples/](examples/)
- Review complete API documentation

## 🤝 Contributing

Contributions, bug reports, and suggestions are welcome!

1. Fork the project
2. Create a feature branch
3. Commit your changes
4. Push to the branch
5. Open a Pull Request

## 📄 License

MIT License - see [LICENSE](LICENSE) for details
