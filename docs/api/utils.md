# Utilities Module

Utility functions for device management, I/O, and D3 correction.

## Device Utilities

::: mace_inference.utils.device
    options:
      show_root_heading: false
      heading_level: 3

### get_device

```python
def get_device(device: str = "auto") -> str
```

Determine the best available device for computation.

**Parameters:**

- `device` (`str`): Device specification
  - `"auto"`: Automatically select best available
  - `"cpu"`: Force CPU
  - `"cuda"`: Use CUDA GPU
  - `"cuda:0"`, `"cuda:1"`: Specific GPU
  - `"mps"`: Apple Silicon GPU

**Returns:** `str` - The device string to use

**Example:**

```python
from mace_inference.utils import get_device

device = get_device("auto")
print(f"Using device: {device}")  # e.g., "cuda" or "cpu"
```

---

## I/O Utilities

::: mace_inference.utils.io
    options:
      show_root_heading: false
      heading_level: 3

### load_structure

```python
def load_structure(path: str) -> Atoms
```

Load atomic structure from file.

**Parameters:**

- `path` (`str`): Path to structure file (CIF, POSCAR, XYZ, etc.)

**Returns:** `Atoms` - ASE Atoms object

### save_structure

```python
def save_structure(atoms: Atoms, path: str, format: Optional[str] = None)
```

Save atomic structure to file.

**Parameters:**

- `atoms` (`Atoms`): Structure to save
- `path` (`str`): Output file path
- `format` (`str`): File format (auto-detected if None)

---

## D3 Correction

::: mace_inference.utils.d3_correction
    options:
      show_root_heading: false
      heading_level: 3

### create_d3_calculator

```python
def create_d3_calculator(
    base_calculator: Calculator,
    xc: str = "pbe",
    damping: str = "bj"
) -> Calculator
```

Create a combined MACE+D3 calculator.

**Parameters:**

- `base_calculator` (`Calculator`): The MACE calculator
- `xc` (`str`): Exchange-correlation functional
- `damping` (`str`): Damping function (`"bj"` or `"zero"`)

**Returns:** `Calculator` - Sum calculator with MACE + D3

**Example:**

```python
from mace_inference.utils.d3_correction import create_d3_calculator
from mace.calculators import MACECalculator

mace_calc = MACECalculator(model_path="medium", device="cpu")
combined_calc = create_d3_calculator(mace_calc, xc="pbe", damping="bj")
```

### D3_AVAILABLE

```python
D3_AVAILABLE: bool
```

Boolean indicating if torch-dftd is installed.

```python
from mace_inference.utils.d3_correction import D3_AVAILABLE

if D3_AVAILABLE:
    print("D3 correction is available")
else:
    print("Install torch-dftd for D3 support")
```
