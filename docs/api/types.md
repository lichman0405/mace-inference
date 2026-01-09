# Types Module

Type definitions and data classes used throughout MACE Inference.

::: mace_inference.types
    options:
      show_root_heading: false
      heading_level: 2

## Type Aliases

### DeviceType

```python
DeviceType = Literal["auto", "cpu", "cuda", "mps"]
```

Valid device specifications.

### OptimizerType

```python
OptimizerType = Literal["BFGS", "FIRE"]
```

Available optimizers for geometry relaxation.

### EnsembleType

```python
EnsembleType = Literal["nvt", "npt"]
```

Molecular dynamics ensemble types.

### EOSType

```python
EOSType = Literal["birchmurnaghan", "murnaghan", "vinet"]
```

Equation of state types for bulk modulus calculation.

## Usage

These types are used for type hints throughout the codebase:

```python
from mace_inference.types import DeviceType, OptimizerType

def my_function(device: DeviceType = "auto") -> None:
    ...
```

They help with:

- IDE autocompletion
- Type checking with mypy
- Documentation
