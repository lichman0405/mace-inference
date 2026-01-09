"""
Tests for utility functions.

Uses real implementations - no mocking.
"""

import pytest
import numpy as np
from ase.build import bulk


class TestDeviceUtils:
    """Test device utility functions"""

    def test_get_device_auto(self):
        """Test automatic device detection"""
        from mace_inference.utils import get_device
        
        device = get_device("auto")
        assert device in ["cpu", "cuda"]

    def test_get_device_cpu(self):
        """Test explicit CPU device"""
        from mace_inference.utils import get_device
        
        device = get_device("cpu")
        assert device == "cpu"

    def test_invalid_device(self):
        """Test invalid device raises error"""
        from mace_inference.utils import validate_device
        
        with pytest.raises(ValueError):
            validate_device("invalid_device")

    def test_get_device_info(self):
        """Test device info retrieval"""
        from mace_inference.utils import get_device_info
        
        info = get_device_info()
        
        assert "pytorch_version" in info
        assert "cuda_available" in info
        assert "cuda_count" in info


class TestIOUtils:
    """Test I/O utility functions"""

    def test_parse_structure_input_atoms(self, cu_bulk):
        """Test parsing Atoms object"""
        from mace_inference.utils import parse_structure_input
        
        result = parse_structure_input(cu_bulk)
        
        assert len(result) == len(cu_bulk)
        assert result.get_chemical_symbols() == cu_bulk.get_chemical_symbols()

    def test_parse_structure_input_file(self, temp_structure_file):
        """Test parsing structure from file"""
        from mace_inference.utils import parse_structure_input
        
        result = parse_structure_input(temp_structure_file)
        
        assert len(result) > 0

    def test_parse_structure_input_invalid(self):
        """Test parsing invalid input raises error"""
        from mace_inference.utils import parse_structure_input
        
        with pytest.raises((FileNotFoundError, ValueError, TypeError)):
            parse_structure_input("nonexistent_file.xyz")

    def test_create_supercell(self, cu_bulk):
        """Test supercell creation"""
        from mace_inference.utils import create_supercell
        
        supercell = create_supercell(cu_bulk, [2, 2, 2])
        
        assert len(supercell) == len(cu_bulk) * 8

    def test_create_supercell_isotropic(self, cu_bulk):
        """Test isotropic supercell creation"""
        from mace_inference.utils import create_supercell
        
        supercell = create_supercell(cu_bulk, 2)
        
        assert len(supercell) == len(cu_bulk) * 8

    def test_save_structure(self, cu_bulk, temp_output_dir):
        """Test saving structure to file"""
        import os
        from mace_inference.utils import save_structure
        
        output_file = os.path.join(temp_output_dir, "test.xyz")
        save_structure(cu_bulk, output_file)
        
        assert os.path.exists(output_file)

    def test_atoms_to_dict(self, cu_bulk):
        """Test converting Atoms to dictionary"""
        from mace_inference.utils import atoms_to_dict
        
        result = atoms_to_dict(cu_bulk)
        
        assert "symbols" in result
        assert "positions" in result
        assert "cell" in result


class TestD3Correction:
    """Test D3 dispersion correction utilities"""

    def test_check_d3_available(self):
        """Test D3 availability check"""
        from mace_inference.utils import check_d3_available
        
        # Should return boolean
        result = check_d3_available()
        assert isinstance(result, bool)


class TestStaticTasks:
    """Test static calculation tasks"""

    def test_single_point_energy_structure(self, mace_calc, cu_bulk):
        """Test single-point energy calculation"""
        from mace_inference.tasks import single_point_energy
        
        result = single_point_energy(cu_bulk, mace_calc.calculator)
        
        assert "energy" in result
        assert "forces" in result
        assert isinstance(result["energy"], float)

    def test_optimize_structure(self, mace_calc, cu_bulk):
        """Test structure optimization"""
        from mace_inference.tasks import optimize_structure
        from ase import Atoms
        
        result = optimize_structure(cu_bulk, mace_calc.calculator, fmax=0.1, steps=5)
        
        assert isinstance(result, Atoms)


class TestDynamicsTasks:
    """Test molecular dynamics tasks"""

    def test_nvt_md(self, mace_calc, cu_supercell, temp_output_dir):
        """Test NVT molecular dynamics"""
        import os
        from mace_inference.tasks import run_nvt_md
        from ase import Atoms
        
        traj_file = os.path.join(temp_output_dir, "nvt.traj")
        
        result = run_nvt_md(
            cu_supercell,
            mace_calc.calculator,
            temperature_K=300,
            steps=3,
            timestep=1.0,
            trajectory=traj_file
        )
        
        assert isinstance(result, Atoms)
        assert os.path.exists(traj_file)
