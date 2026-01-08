"""
Tests for core MACEInference class functionality.

These tests use mock calculators to avoid requiring mace-torch installation.
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from ase.build import bulk


class TestMACEInferenceInit:
    """Test MACEInference initialization"""

    def test_init_attributes(self):
        """Test that initialization sets correct attributes"""
        with patch('mace_inference.core.mace_mp') as mock_mace_mp:
            mock_calc = Mock()
            mock_mace_mp.return_value = mock_calc
            
            from mace_inference import MACEInference
            
            with patch.object(MACEInference, '_create_calculator', return_value=mock_calc):
                calc = MACEInference(model="medium", device="cpu", enable_d3=False)
                
                assert calc.model_name == "medium"
                assert calc.device == "cpu"
                assert calc.enable_d3 is False
                assert calc.default_dtype == "float64"

    def test_init_with_d3(self):
        """Test initialization with D3 correction enabled"""
        mock_calc = Mock()
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu", enable_d3=True)
            assert calc.enable_d3 is True
            assert calc.d3_damping == "bj"
            assert calc.d3_xc == "pbe"

    def test_repr(self):
        """Test string representation"""
        mock_calc = Mock()
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="large", device="cpu", enable_d3=True)
            repr_str = repr(calc)
            
            assert "MACEInference" in repr_str
            assert "large" in repr_str
            assert "enable_d3=True" in repr_str


class TestSinglePoint:
    """Test single-point energy calculation"""

    def test_single_point_with_atoms(self, cu_bulk):
        """Test single-point calculation with Atoms object"""
        n_atoms = len(cu_bulk)
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(return_value=-3.5 * n_atoms)
        mock_calc.get_forces = Mock(return_value=np.zeros((n_atoms, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.single_point(cu_bulk)
            
            assert "energy" in result
            assert "energy_per_atom" in result
            assert "forces" in result
            assert "max_force" in result
            assert "rms_force" in result
            assert result["energy"] == pytest.approx(-3.5 * n_atoms)

    def test_single_point_with_file(self, temp_structure_file):
        """Test single-point calculation with file path"""
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(return_value=-10.0)
        mock_calc.get_forces = Mock(return_value=np.zeros((4, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.single_point(str(temp_structure_file))
            
            assert "energy" in result
            mock_calc.get_potential_energy.assert_called()


class TestOptimize:
    """Test structure optimization"""

    def test_optimize_returns_atoms(self, cu_bulk):
        """Test that optimization returns Atoms object"""
        n_atoms = len(cu_bulk)
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(return_value=-10.0)
        mock_calc.get_forces = Mock(return_value=np.zeros((n_atoms, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            from ase import Atoms
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.optimize(cu_bulk, fmax=0.05, steps=5)
            
            assert isinstance(result, Atoms)
            assert len(result) == n_atoms

    def test_optimize_with_output(self, cu_bulk, tmp_path):
        """Test optimization with output file"""
        n_atoms = len(cu_bulk)
        output_file = tmp_path / "optimized.xyz"
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(return_value=-10.0)
        mock_calc.get_forces = Mock(return_value=np.zeros((n_atoms, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu")
            calc.optimize(cu_bulk, fmax=0.05, steps=5, output=str(output_file))
            
            assert output_file.exists()


class TestMolecularDynamics:
    """Test molecular dynamics simulations"""

    def test_run_md_nvt(self, cu_bulk):
        """Test NVT molecular dynamics"""
        n_atoms = len(cu_bulk)
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(return_value=-10.0)
        mock_calc.get_forces = Mock(return_value=np.zeros((n_atoms, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            from ase import Atoms
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.run_md(
                cu_bulk,
                ensemble="nvt",
                temperature_K=300,
                steps=10,
                timestep=1.0
            )
            
            assert isinstance(result, Atoms)

    def test_run_md_npt(self, cu_bulk):
        """Test NPT molecular dynamics"""
        n_atoms = len(cu_bulk)
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(return_value=-10.0)
        mock_calc.get_forces = Mock(return_value=np.zeros((n_atoms, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            from ase import Atoms
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.run_md(
                cu_bulk,
                ensemble="npt",
                temperature_K=300,
                pressure_GPa=0.0,
                steps=10,
                timestep=1.0
            )
            
            assert isinstance(result, Atoms)

    def test_run_md_invalid_ensemble(self, cu_bulk):
        """Test invalid ensemble raises error"""
        mock_calc = Mock()
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu")
            
            with pytest.raises(ValueError, match="Invalid ensemble"):
                calc.run_md(cu_bulk, ensemble="invalid")


class TestBulkModulus:
    """Test bulk modulus calculation"""

    def test_bulk_modulus_result_keys(self, cu_bulk):
        """Test bulk modulus returns expected keys"""
        n_atoms = len(cu_bulk)
        
        # Create mock that returns different energies for different volumes
        call_count = [0]
        def mock_energy():
            call_count[0] += 1
            # Simulate E-V curve (parabolic)
            return -10.0 + 0.1 * (call_count[0] - 6) ** 2
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(side_effect=mock_energy)
        mock_calc.get_forces = Mock(return_value=np.zeros((n_atoms, 3)))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.bulk_modulus(cu_bulk, n_points=11)
            
            assert "v0" in result
            assert "e0" in result
            assert "B_GPa" in result
            assert "volumes" in result
            assert "energies" in result


class TestAdsorptionEnergy:
    """Test adsorption energy calculation"""

    def test_adsorption_energy_result_keys(self, simple_mof):
        """Test adsorption energy returns expected keys"""
        n_mof = len(simple_mof)
        
        # Mock calculator with different energies for different structures
        energy_values = [-100.0, -5.0, -108.0]  # MOF, gas, complex
        call_idx = [0]
        
        def mock_energy():
            idx = min(call_idx[0], len(energy_values) - 1)
            call_idx[0] += 1
            return energy_values[idx]
        
        mock_calc = Mock()
        mock_calc.get_potential_energy = Mock(side_effect=mock_energy)
        mock_calc.get_forces = Mock(return_value=np.zeros((20, 3)))
        mock_calc.get_stress = Mock(return_value=np.zeros(6))
        
        with patch('mace_inference.core.MACEInference._create_calculator', return_value=mock_calc):
            from mace_inference import MACEInference
            
            calc = MACEInference(model="medium", device="cpu")
            result = calc.adsorption_energy(
                mof_structure=simple_mof,
                gas_molecule="CO2",
                site_position=[5.0, 5.0, 5.0],
                optimize_complex=False  # Skip optimization for faster test
            )
            
            assert "E_ads" in result
            assert "E_mof" in result
            assert "E_gas" in result
            assert "E_complex" in result


class TestCoordination:
    """Test coordination analysis"""

    def test_coordination_analysis(self, simple_mof):
        """Test coordination analysis returns expected structure"""
        from mace_inference.tasks.adsorption import analyze_coordination
        
        result = analyze_coordination(simple_mof)
        
        assert "coordination" in result
        assert "n_metal_centers" in result
        assert "metal_indices" in result
        # Cu atoms should be detected
        assert result["n_metal_centers"] == 2
