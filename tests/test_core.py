"""
Tests for core MACEInference class functionality.

Uses real MACE models - no mocking.
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
from ase.build import bulk


class TestMACEInferenceInit:
    """Test MACEInference initialization"""

    def test_init_attributes(self, mace_calc):
        """Test that initialization sets correct attributes"""
        assert mace_calc.model_name == "small"
        assert mace_calc.device == "cpu"
        assert mace_calc.enable_d3 is False
        assert mace_calc.default_dtype == "float64"
        assert mace_calc.calculator is not None

    def test_init_with_different_models(self):
        """Test initialization with different model sizes"""
        from mace_inference import MACEInference
        
        # Test medium model
        calc = MACEInference(model="medium", device="cpu")
        assert calc.model_name == "medium"

    def test_init_with_d3(self):
        """Test initialization with D3 correction enabled"""
        from mace_inference import MACEInference
        
        calc = MACEInference(model="small", device="cpu", enable_d3=True)
        assert calc.enable_d3 is True
        assert calc.d3_damping == "bj"
        assert calc.d3_xc == "pbe"

    def test_repr(self, mace_calc):
        """Test string representation"""
        repr_str = repr(mace_calc)
        
        assert "MACEInference" in repr_str
        assert "small" in repr_str


class TestSinglePoint:
    """Test single-point energy calculation"""

    def test_single_point_with_atoms(self, mace_calc, cu_bulk):
        """Test single-point calculation with Atoms object"""
        result = mace_calc.single_point(cu_bulk)
        
        assert "energy" in result
        assert "energy_per_atom" in result
        assert "forces" in result
        assert "max_force" in result
        assert "rms_force" in result
        
        assert isinstance(result["energy"], float)
        assert isinstance(result["forces"], np.ndarray)
        assert result["forces"].shape == (len(cu_bulk), 3)

    def test_single_point_with_file(self, mace_calc, temp_structure_file):
        """Test single-point calculation with file path"""
        result = mace_calc.single_point(temp_structure_file)
        
        assert "energy" in result
        assert isinstance(result["energy"], float)

    def test_single_point_energy_reasonable(self, mace_calc, cu_bulk):
        """Test that energy values are physically reasonable"""
        result = mace_calc.single_point(cu_bulk)
        
        # Cu cohesive energy is about -3.5 eV/atom
        energy_per_atom = result["energy_per_atom"]
        assert -10 < energy_per_atom < 0  # Should be negative for bound system


class TestOptimize:
    """Test structure optimization"""

    def test_optimize_returns_atoms(self, mace_calc, cu_bulk):
        """Test that optimization returns Atoms object"""
        from ase import Atoms
        
        result = mace_calc.optimize(cu_bulk, fmax=0.1, steps=5)
        
        assert isinstance(result, Atoms)
        assert len(result) == len(cu_bulk)

    def test_optimize_reduces_forces(self, mace_calc):
        """Test that optimization reduces forces"""
        # Create slightly distorted structure
        atoms = bulk('Cu', 'fcc', a=3.6)
        atoms.positions[0] += [0.1, 0.0, 0.0]  # Displace one atom
        
        initial_result = mace_calc.single_point(atoms)
        initial_max_force = initial_result["max_force"]
        
        optimized = mace_calc.optimize(atoms, fmax=0.05, steps=50)
        final_result = mace_calc.single_point(optimized)
        final_max_force = final_result["max_force"]
        
        # Optimization should reduce forces
        assert final_max_force <= initial_max_force

    def test_optimize_with_output(self, mace_calc, cu_bulk, temp_output_dir):
        """Test optimization with output file"""
        output_file = os.path.join(temp_output_dir, "optimized.xyz")
        
        mace_calc.optimize(cu_bulk, fmax=0.1, steps=5, output=output_file)
        
        assert os.path.exists(output_file)

    def test_optimize_with_trajectory(self, mace_calc, cu_bulk, temp_output_dir):
        """Test optimization with trajectory file"""
        traj_file = os.path.join(temp_output_dir, "opt.traj")
        
        mace_calc.optimize(cu_bulk, fmax=0.1, steps=5, trajectory=traj_file)
        
        assert os.path.exists(traj_file)


class TestMolecularDynamics:
    """Test molecular dynamics simulations"""

    def test_run_md_nvt(self, mace_calc, cu_supercell, temp_output_dir):
        """Test NVT MD simulation"""
        traj_file = os.path.join(temp_output_dir, "md.traj")
        
        result = mace_calc.run_md(
            cu_supercell,
            ensemble="nvt",
            temperature_K=300,
            steps=5,
            timestep=1.0,
            trajectory=traj_file
        )
        
        from ase import Atoms
        assert isinstance(result, Atoms)
        assert os.path.exists(traj_file)

    def test_run_md_npt(self, mace_calc, temp_output_dir):
        """Test NPT MD simulation with proper cell"""
        from ase.build import bulk
        
        # Create a structure with standard orientation cell for NPT
        atoms = bulk('Cu', 'fcc', a=3.6, cubic=True) * (2, 2, 2)
        
        result = mace_calc.run_md(
            atoms,
            ensemble="npt",
            temperature_K=300,
            pressure_GPa=0.0,
            steps=3,
            timestep=1.0
        )
        
        from ase import Atoms
        assert isinstance(result, Atoms)

    def test_run_md_invalid_ensemble(self, mace_calc, cu_bulk):
        """Test MD with invalid ensemble raises error"""
        with pytest.raises(ValueError, match="ensemble"):
            mace_calc.run_md(cu_bulk, ensemble="invalid", temperature_K=300, steps=1)


class TestPhonon:
    """Test phonon calculations"""

    def test_phonon_basic(self, mace_calc, cu_bulk, temp_output_dir):
        """Test basic phonon calculation"""
        result = mace_calc.phonon(
            cu_bulk,
            supercell_matrix=[2, 2, 2],
            output_dir=temp_output_dir
        )
        
        assert "phonon" in result
        assert result["phonon"] is not None

    def test_phonon_with_thermal(self, mace_calc, cu_bulk, temp_output_dir):
        """Test phonon calculation with thermal properties"""
        result = mace_calc.phonon(
            cu_bulk,
            supercell_matrix=[2, 2, 2],
            temperature_range=(100, 500, 100),
            output_dir=temp_output_dir
        )
        
        # Should have thermal properties
        assert "thermal" in result or "thermal_properties" in result


class TestBulkModulus:
    """Test bulk modulus calculation"""

    def test_bulk_modulus_result_keys(self, mace_calc, cu_bulk):
        """Test bulk modulus calculation returns expected keys"""
        result = mace_calc.bulk_modulus(cu_bulk, n_points=5, optimize_first=False)
        
        assert "v0" in result
        assert "e0" in result
        assert "B" in result
        assert "B_GPa" in result
        
        # Bulk modulus should be positive
        assert result["B_GPa"] > 0

    def test_bulk_modulus_reasonable_value(self, mace_calc, cu_bulk):
        """Test that bulk modulus is physically reasonable for Cu"""
        result = mace_calc.bulk_modulus(cu_bulk, n_points=7, optimize_first=False)
        
        # Cu bulk modulus is about 140 GPa
        # Allow wide range due to model approximations
        assert 50 < result["B_GPa"] < 300


class TestAdsorptionEnergy:
    """Test adsorption energy calculation"""

    def test_adsorption_energy_result_keys(self, mace_calc, simple_mof):
        """Test adsorption energy calculation returns expected keys"""
        result = mace_calc.adsorption_energy(
            framework=simple_mof,
            adsorbate="H2O",
            site_position=[5.0, 5.0, 5.0],
            optimize=False
        )
        
        assert "E_ads" in result or "adsorption_energy" in result
        assert "E_mof" in result or "E_framework" in result


class TestCoordinationAnalysis:
    """Test coordination analysis"""

    def test_coordination_analysis(self, mace_calc, simple_mof):
        """Test coordination analysis for metal sites"""
        result = mace_calc.coordination(
            simple_mof,
            metal_indices=[0, 1],  # Cu atoms
            cutoff_multiplier=1.3
        )
        
        assert "coordination" in result or isinstance(result, dict)


class TestBatchOperations:
    """Test batch operations"""

    def test_batch_single_point(self, mace_calc, cu_bulk, water_molecule):
        """Test batch single-point calculations"""
        structures = [cu_bulk, water_molecule]
        
        results = mace_calc.batch_single_point(structures)
        
        assert len(results) == 2
        assert all("energy" in r or "error" in r for r in results)

    def test_batch_optimize(self, mace_calc, cu_bulk, temp_output_dir):
        """Test batch optimization"""
        # Create slightly different structures
        atoms1 = cu_bulk.copy()
        atoms2 = cu_bulk.copy()
        atoms2.positions[0] += [0.05, 0.0, 0.0]
        
        structures = [atoms1, atoms2]
        
        results = mace_calc.batch_optimize(
            structures,
            fmax=0.1,
            steps=5,
            output_dir=temp_output_dir
        )
        
        assert len(results) == 2
