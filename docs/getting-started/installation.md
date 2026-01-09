# Installation

This guide covers how to install MACE Inference on your system.

## Requirements

- **Python**: 3.9, 3.10, or 3.11
- **Operating System**: Linux, macOS, or Windows
- **Hardware**: CPU or NVIDIA GPU (for acceleration)

## Installation Methods

### Method 1: pip (Recommended)

=== "CPU Only"

    ```bash
    pip install mace-inference
    ```

=== "GPU (CUDA)"

    ```bash
    # Install PyTorch with CUDA first
    pip install torch --index-url https://download.pytorch.org/whl/cu118
    
    # Then install mace-inference
    pip install mace-inference
    ```

### Method 2: From Source

```bash
# Clone the repository
git clone https://github.com/lichman0405/mace-inference.git
cd mace-inference

# Install in development mode
pip install -e ".[dev]"
```

### Method 3: Conda Environment

We recommend using a dedicated conda environment:

=== "CPU"

    ```bash
    # Create environment
    conda create -n mace-inference python=3.10
    conda activate mace-inference
    
    # Install dependencies
    pip install -r requirements-cpu.txt
    pip install -e .
    ```

=== "GPU"

    ```bash
    # Create environment
    conda create -n mace-inference-gpu python=3.10
    conda activate mace-inference-gpu
    
    # Install PyTorch with CUDA
    conda install pytorch pytorch-cuda=11.8 -c pytorch -c nvidia
    
    # Install dependencies
    pip install -r requirements-gpu.txt
    pip install -e .
    ```

## Optional Dependencies

### D3 Dispersion Correction

For DFT-D3 dispersion correction support:

```bash
pip install torch-dftd
```

### Development Tools

For contributing to the project:

```bash
pip install -e ".[dev]"
```

This includes:

- `pytest` - Testing framework
- `ruff` - Linting
- `mypy` - Type checking

## Verify Installation

After installation, verify everything works:

```python
from mace_inference import MACEInference

# This will download the model on first run
calc = MACEInference(model="medium", device="cpu")
print("✓ MACE Inference installed successfully!")
```

Or use the CLI:

```bash
mace-inference --help
```

## Troubleshooting

### CUDA Not Detected

If GPU is not detected:

```python
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"CUDA version: {torch.version.cuda}")
```

Make sure your PyTorch installation matches your CUDA version.

### Model Download Issues

Models are cached in `~/.cache/mace/`. If download fails:

1. Check internet connection
2. Try manual download from [MACE Models](https://github.com/ACEsuit/mace)
3. Specify local model path:

```python
calc = MACEInference(model="/path/to/model.model")
```

### Memory Issues

For large systems, try:

```python
# Use lower precision
calc = MACEInference(model="medium", default_dtype="float32")

# Or reduce batch size in dynamics
result = calc.run_md(atoms, steps=1000, save_interval=100)
```

## Next Steps

- [Quick Start Guide](quickstart.md) - Learn the basics
- [User Guide](../user-guide/overview.md) - Detailed tutorials
- [API Reference](../api/core.md) - Complete documentation
