"""
MACE Inference: High-level Python library for MACE machine learning force field inference tasks
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("mace-inference")
except PackageNotFoundError:
    # Package is not installed (development mode without editable install)
    __version__ = "0.1.0.dev0"

__author__ = "Shibo Li"
__email__ = "shadow.li981@gmail.com"

from mace_inference.core import MACEInference
from mace_inference.utils.device import get_device
from mace_inference.types import (
    SinglePointResult,
    BulkModulusResult,
    PhononResult,
    AdsorptionResult,
    CoordinationResult,
    ThermalPropertiesResult,
    Calculator,
)

__all__ = [
    "MACEInference",
    "get_device",
    "__version__",
    # Type exports
    "SinglePointResult",
    "BulkModulusResult",
    "PhononResult",
    "AdsorptionResult",
    "CoordinationResult",
    "ThermalPropertiesResult",
    "Calculator",
]
