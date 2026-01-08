"""
Pytest configuration and shared fixtures for MACE inference tests.
"""

import pytest
import numpy as np
from unittest.mock import Mock, MagicMock
from ase import Atoms
from ase.build import bulk, molecule


# =============================================================================
# Structure Fixtures
# =============================================================================

@pytest.fixture
def cu_bulk():
    """Create a simple Cu FCC bulk structure."""
    return bulk('Cu', 'fcc', a=3.6)


@pytest.fixture
def cu_supercell(cu_bulk):
    """Create a 2x2x2 Cu supercell."""
    from mace_inference.utils import create_supercell
    return create_supercell(cu_bulk, 2)


@pytest.fixture
def water_molecule():
    """Create a water molecule."""
    return molecule('H2O')


@pytest.fixture
def co2_molecule():
    """Create a CO2 molecule."""
    return molecule('CO2')


@pytest.fixture
def simple_mof():
    """Create a simple MOF-like structure for testing."""
    # Simple Cu-based structure (Cu paddlewheel mimic)
    atoms = Atoms(
        symbols=['Cu', 'Cu', 'O', 'O', 'O', 'O', 'C', 'C'],
        positions=[
            [0.0, 0.0, 0.0],      # Cu1
            [2.5, 0.0, 0.0],      # Cu2
            [1.25, 1.0, 0.0],     # O1
            [1.25, -1.0, 0.0],    # O2
            [1.25, 0.0, 1.0],     # O3
            [1.25, 0.0, -1.0],    # O4
            [1.25, 1.5, 0.5],     # C1
            [1.25, -1.5, -0.5],   # C2
        ],
        cell=[10.0, 10.0, 10.0],
        pbc=True
    )
    return atoms


# =============================================================================
# Mock Calculator Fixtures
# =============================================================================

@pytest.fixture
def mock_calculator():
    """Create a mock ASE calculator for testing."""
    calc = Mock()
    
    # Default return values
    calc.get_potential_energy = Mock(return_value=-10.0)
    calc.get_forces = Mock(return_value=np.zeros((1, 3)))
    calc.get_stress = Mock(return_value=np.zeros(6))
    
    return calc


@pytest.fixture
def mock_calculator_for_atoms(cu_bulk):
    """Create a mock calculator configured for Cu bulk structure."""
    n_atoms = len(cu_bulk)
    
    calc = Mock()
    calc.get_potential_energy = Mock(return_value=-3.5 * n_atoms)
    calc.get_forces = Mock(return_value=np.random.randn(n_atoms, 3) * 0.01)
    calc.get_stress = Mock(return_value=np.array([0.01, 0.01, 0.01, 0.0, 0.0, 0.0]))
    
    return calc


@pytest.fixture
def mock_mace_inference():
    """Create a mock MACEInference instance for testing."""
    from unittest.mock import patch
    
    mock_calc = Mock()
    mock_calc.get_potential_energy = Mock(return_value=-10.0)
    mock_calc.get_forces = Mock(return_value=np.zeros((4, 3)))
    mock_calc.get_stress = Mock(return_value=np.zeros(6))
    
    with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
        from mace_inference import MACEInference
        
        # Create instance with patched calculator
        instance = MagicMock(spec=MACEInference)
        instance.device = "cpu"
        instance.model_name = "medium"
        instance.enable_d3 = False
        instance.calculator = mock_calc
        
        yield instance


# =============================================================================
# Temporary File Fixtures
# =============================================================================

@pytest.fixture
def temp_structure_file(tmp_path, cu_bulk):
    """Create a temporary structure file for testing."""
    from ase.io import write
    
    filepath = tmp_path / "test_structure.xyz"
    write(str(filepath), cu_bulk)
    return filepath


@pytest.fixture
def temp_cif_file(tmp_path, cu_bulk):
    """Create a temporary CIF file for testing."""
    from ase.io import write
    
    filepath = tmp_path / "test_structure.cif"
    write(str(filepath), cu_bulk)
    return filepath


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create a temporary output directory."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir


# =============================================================================
# Skip Markers
# =============================================================================

def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "requires_mace: mark test as requiring mace-torch installation"
    )
    config.addinivalue_line(
        "markers", "requires_gpu: mark test as requiring CUDA GPU"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )


@pytest.fixture
def skip_without_mace():
    """Skip test if mace-torch is not installed."""
    try:
        import mace
        return False
    except ImportError:
        pytest.skip("mace-torch not installed")
        return True


@pytest.fixture
def skip_without_gpu():
    """Skip test if CUDA GPU is not available."""
    import torch
    if not torch.cuda.is_available():
        pytest.skip("CUDA GPU not available")
        return True
    return False


# =============================================================================
# Helper Functions
# =============================================================================

def assert_atoms_equal(atoms1: Atoms, atoms2: Atoms, rtol: float = 1e-5):
    """Assert that two Atoms objects are approximately equal."""
    assert len(atoms1) == len(atoms2), "Number of atoms differs"
    assert atoms1.get_chemical_formula() == atoms2.get_chemical_formula(), "Formula differs"
    assert np.allclose(atoms1.positions, atoms2.positions, rtol=rtol), "Positions differ"
    assert np.allclose(atoms1.cell[:], atoms2.cell[:], rtol=rtol), "Cell differs"


def assert_result_has_keys(result: dict, required_keys: list):
    """Assert that a result dictionary has all required keys."""
    for key in required_keys:
        assert key in result, f"Missing required key: {key}"
