"""
Type definitions for MACE Inference.

This module provides TypedDict definitions for structured return types
and Protocol classes for interface definitions.
"""

from typing import (
    Any,
    Dict,
    List,
    Optional,
    Protocol,
    Sequence,
    Union,
    runtime_checkable,
)

import numpy as np
from numpy.typing import NDArray

# Try to import TypedDict from typing (Python 3.8+)
# For Python 3.8/3.9, we need typing_extensions for NotRequired
import sys

if sys.version_info >= (3, 11):
    from typing import TypedDict, NotRequired
else:
    from typing import TypedDict
    try:
        from typing_extensions import NotRequired
    except ImportError:
        # Fallback: NotRequired not available, use total=False pattern
        NotRequired = None  # type: ignore


# =============================================================================
# Array Types
# =============================================================================

ArrayLike = Union[List[float], np.ndarray, Sequence[float]]
ForceArray = NDArray[np.floating[Any]]  # Shape: (n_atoms, 3)
StressArray = NDArray[np.floating[Any]]  # Shape: (6,) Voigt notation
PositionArray = NDArray[np.floating[Any]]  # Shape: (n_atoms, 3)
CellArray = NDArray[np.floating[Any]]  # Shape: (3, 3)


# =============================================================================
# Result TypedDicts
# =============================================================================

class SinglePointResult(TypedDict, total=False):
    """Result from single-point energy calculation."""
    energy: float
    """Total potential energy in eV."""
    
    energy_per_atom: float
    """Energy per atom in eV/atom."""
    
    forces: ForceArray
    """Atomic forces in eV/Å, shape (n_atoms, 3)."""
    
    max_force: float
    """Maximum force component magnitude in eV/Å."""
    
    rms_force: float
    """Root mean square of forces in eV/Å."""
    
    stress: Optional[StressArray]
    """Stress tensor in Voigt notation (eV/Å³), shape (6,)."""
    
    pressure_GPa: Optional[float]
    """Hydrostatic pressure in GPa."""


class BulkModulusResult(TypedDict):
    """Result from bulk modulus calculation."""
    v0: float
    """Equilibrium volume in Å³."""
    
    e0: float
    """Equilibrium energy in eV."""
    
    B: float
    """Bulk modulus in eV/Å³."""
    
    B_GPa: float
    """Bulk modulus in GPa."""
    
    volumes: List[float]
    """List of volumes used for fitting in Å³."""
    
    energies: List[float]
    """List of energies used for fitting in eV."""
    
    eos_type: str
    """Equation of state type used for fitting."""


class ThermalPropertiesResult(TypedDict):
    """Result from thermal properties calculation."""
    temperatures: NDArray[np.floating[Any]]
    """Temperature array in K."""
    
    free_energy: NDArray[np.floating[Any]]
    """Helmholtz free energy in kJ/mol."""
    
    entropy: NDArray[np.floating[Any]]
    """Entropy in J/(mol·K)."""
    
    heat_capacity: NDArray[np.floating[Any]]
    """Heat capacity at constant volume in J/(mol·K)."""


class PhononResult(TypedDict, total=False):
    """Result from phonon calculation."""
    phonon: Any  # Phonopy object
    """Phonopy object with calculated force constants."""
    
    supercell_matrix: List[List[int]]
    """Supercell transformation matrix."""
    
    displacement: float
    """Atomic displacement distance in Å."""
    
    mesh: List[int]
    """k-point mesh for phonon DOS."""
    
    thermal_properties: ThermalPropertiesResult
    """Thermal properties if temperature_range was specified."""


class AdsorptionResult(TypedDict):
    """Result from adsorption energy calculation."""
    E_ads: float
    """Adsorption energy in eV."""
    
    E_mof: float
    """MOF (substrate) energy in eV."""
    
    E_gas: float
    """Gas molecule energy in eV."""
    
    E_complex: float
    """Adsorption complex energy in eV."""
    
    complex_structure: Any  # Atoms object
    """Final structure of the adsorption complex."""
    
    optimized: bool
    """Whether the complex was optimized."""


class NeighborInfo(TypedDict):
    """Information about a neighboring atom."""
    index: int
    """Atom index in the structure."""
    
    symbol: str
    """Chemical symbol of the atom."""
    
    distance: float
    """Distance from the central atom in Å."""


class MetalCoordinationInfo(TypedDict):
    """Coordination information for a single metal center."""
    metal_symbol: str
    """Chemical symbol of the metal atom."""
    
    coordination_number: int
    """Number of coordinating atoms."""
    
    neighbors: List[NeighborInfo]
    """List of neighboring atoms with distances."""
    
    average_distance: float
    """Average distance to neighbors in Å."""


class CoordinationResult(TypedDict):
    """Result from coordination analysis."""
    coordination: Dict[int, MetalCoordinationInfo]
    """Coordination info for each metal center, keyed by atom index."""
    
    n_metal_centers: int
    """Total number of metal centers found."""
    
    metal_indices: List[int]
    """List of metal atom indices."""


class DeviceInfo(TypedDict, total=False):
    """Information about compute devices."""
    cuda_available: bool
    """Whether CUDA is available."""
    
    cuda_count: int
    """Number of CUDA devices."""
    
    pytorch_version: str
    """PyTorch version string."""
    
    cuda_version: str
    """CUDA version string (if available)."""
    
    devices: List[Dict[str, Any]]
    """List of GPU device info dictionaries."""


# =============================================================================
# Protocol Classes
# =============================================================================

@runtime_checkable
class Calculator(Protocol):
    """Protocol for ASE-compatible calculators."""
    
    def get_potential_energy(self, atoms: Any = None, force_consistent: bool = False) -> float:
        """Calculate and return the potential energy."""
        ...
    
    def get_forces(self, atoms: Any = None) -> ForceArray:
        """Calculate and return the atomic forces."""
        ...
    
    def get_stress(self, atoms: Any = None, voigt: bool = True) -> StressArray:
        """Calculate and return the stress tensor."""
        ...


@runtime_checkable
class Optimizable(Protocol):
    """Protocol for objects that can be optimized."""
    
    def run(self, fmax: float = 0.05, steps: int = 500) -> bool:
        """Run the optimization."""
        ...


# =============================================================================
# Type Aliases for Common Patterns
# =============================================================================

# Structure input types
StructureInput = Union[str, "PathLike[str]", Any]  # str, Path, or Atoms

# Supercell specification
SupercellMatrix = Union[int, List[int], NDArray[np.integer[Any]]]

# Device specification
DeviceType = Union[str, None]  # "auto", "cpu", "cuda", or None

# Model specification
ModelType = str  # "small", "medium", "large", or path to model file

# Try to import PathLike for proper typing
try:
    from os import PathLike
except ImportError:
    PathLike = Any  # type: ignore


# =============================================================================
# Exports
# =============================================================================

__all__ = [
    # Array types
    "ArrayLike",
    "ForceArray",
    "StressArray",
    "PositionArray",
    "CellArray",
    # Result types
    "SinglePointResult",
    "BulkModulusResult",
    "ThermalPropertiesResult",
    "PhononResult",
    "AdsorptionResult",
    "NeighborInfo",
    "MetalCoordinationInfo",
    "CoordinationResult",
    "DeviceInfo",
    # Protocols
    "Calculator",
    "Optimizable",
    # Type aliases
    "StructureInput",
    "SupercellMatrix",
    "DeviceType",
    "ModelType",
]
