"""
Pytest configuration and shared fixtures for MACE inference tests.

Uses real MACE models - no mocking.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
from ase import Atoms
from ase.build import bulk, molecule


# =============================================================================
# Structure Fixtures
# =============================================================================

@pytest.fixture
def cu_bulk():
    """Create a simple Cu FCC bulk structure."""
    atoms = bulk('Cu', 'fcc', a=3.6)
    return atoms


@pytest.fixture
def cu_supercell():
    """Create a 2x2x2 Cu supercell."""
    atoms = bulk('Cu', 'fcc', a=3.6) * (2, 2, 2)
    return atoms


@pytest.fixture
def water_molecule():
    """Create a water molecule in a box."""
    atoms = molecule('H2O')
    atoms.center(vacuum=5.0)
    atoms.pbc = False
    return atoms


@pytest.fixture
def co2_molecule():
    """Create a CO2 molecule in a box."""
    atoms = molecule('CO2')
    atoms.center(vacuum=5.0)
    atoms.pbc = False
    return atoms


@pytest.fixture
def simple_mof():
    """Create a simple MOF-like structure for testing."""
    atoms = Atoms(
        symbols=['Cu', 'Cu', 'O', 'O', 'O', 'O', 'C', 'C'],
        positions=[
            [0.0, 0.0, 0.0],
            [2.5, 0.0, 0.0],
            [1.25, 1.0, 0.0],
            [1.25, -1.0, 0.0],
            [1.25, 0.0, 1.0],
            [1.25, 0.0, -1.0],
            [1.25, 1.5, 0.5],
            [1.25, -1.5, -0.5],
        ],
        cell=[10.0, 10.0, 10.0],
        pbc=True
    )
    return atoms


@pytest.fixture
def nacl_bulk():
    """Create a NaCl rocksalt structure."""
    atoms = bulk('NaCl', 'rocksalt', a=5.64)
    return atoms


# =============================================================================
# Calculator Fixture
# =============================================================================

@pytest.fixture(scope="module")
def mace_calc():
    """Create a real MACEInference calculator (module-scoped for efficiency)."""
    from mace_inference import MACEInference
    calc = MACEInference(model="small", device="cpu")
    return calc


# =============================================================================
# File Fixtures
# =============================================================================

@pytest.fixture
def temp_structure_file(cu_bulk):
    """Create a temporary structure file."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        from ase.io import write
        write(f.name, cu_bulk, format='extxyz')
        filepath = f.name
    yield filepath
    if os.path.exists(filepath):
        os.remove(filepath)


@pytest.fixture
def temp_cif_file():
    """Create a temporary CIF file with Cu bulk."""
    cif_content = """data_Cu
_symmetry_space_group_name_H-M   'F m -3 m'
_cell_length_a   3.6
_cell_length_b   3.6
_cell_length_c   3.6
_cell_angle_alpha   90.0
_cell_angle_beta    90.0
_cell_angle_gamma   90.0

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Cu1 Cu 0.0 0.0 0.0
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
        f.write(cif_content)
        filepath = f.name
    yield filepath
    if os.path.exists(filepath):
        os.remove(filepath)


@pytest.fixture
def temp_output_dir():
    """Create a temporary output directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir
