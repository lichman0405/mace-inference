"""
Tests for CLI (Command Line Interface) functionality.

Uses Click's CliRunner for isolated testing without subprocess.
"""

import pytest
from click.testing import CliRunner
from unittest.mock import Mock, patch, MagicMock
import numpy as np


class TestCLIBasic:
    """Test basic CLI functionality"""

    @pytest.fixture
    def runner(self):
        """Create CLI runner"""
        return CliRunner()

    def test_main_help(self, runner):
        """Test main help message"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['--help'])
        assert result.exit_code == 0
        assert 'MACE Inference' in result.output
        assert 'energy' in result.output
        assert 'optimize' in result.output
        assert 'md' in result.output

    def test_main_version(self, runner):
        """Test version display"""
        from mace_inference.cli import main
        import mace_inference
        
        result = runner.invoke(main, ['--version'])
        assert result.exit_code == 0


class TestEnergyCommand:
    """Test energy command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_energy_help(self, runner):
        """Test energy command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['energy', '--help'])
        assert result.exit_code == 0
        assert 'structure' in result.output.lower()
        assert '--model' in result.output
        assert '--device' in result.output

    def test_energy_missing_structure(self, runner):
        """Test energy command with missing structure"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['energy'])
        assert result.exit_code != 0

    def test_energy_with_mock(self, runner, temp_structure_file):
        """Test energy command with mocked calculator"""
        from mace_inference.cli import main
        
        mock_result = {
            'energy': -10.0,
            'energy_per_atom': -2.5,
            'max_force': 0.01,
            'rms_force': 0.005,
            'pressure_GPa': 0.1
        }
        
        with patch('mace_inference.cli.MACEInference') as MockCalc:
            mock_instance = MagicMock()
            mock_instance.single_point.return_value = mock_result
            MockCalc.return_value = mock_instance
            
            result = runner.invoke(main, ['energy', str(temp_structure_file)])
            
            assert result.exit_code == 0
            assert 'Energy' in result.output


class TestOptimizeCommand:
    """Test optimize command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_optimize_help(self, runner):
        """Test optimize command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['optimize', '--help'])
        assert result.exit_code == 0
        assert '--fmax' in result.output
        assert '--steps' in result.output
        assert '--optimizer' in result.output

    def test_optimize_with_mock(self, runner, temp_structure_file, tmp_path):
        """Test optimize command with mocked calculator"""
        from mace_inference.cli import main
        from ase.build import bulk
        
        output_file = tmp_path / "optimized.xyz"
        mock_atoms = bulk('Cu', 'fcc', a=3.6)
        
        with patch('mace_inference.cli.MACEInference') as MockCalc:
            mock_instance = MagicMock()
            mock_instance.optimize.return_value = mock_atoms
            MockCalc.return_value = mock_instance
            
            result = runner.invoke(main, [
                'optimize', str(temp_structure_file),
                '--fmax', '0.05',
                '--steps', '10',
                '--output', str(output_file)
            ])
            
            assert result.exit_code == 0
            assert 'completed' in result.output.lower() or 'Optimization' in result.output


class TestMDCommand:
    """Test md command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_md_help(self, runner):
        """Test md command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['md', '--help'])
        assert result.exit_code == 0
        assert '--ensemble' in result.output
        assert '--temp' in result.output
        assert '--steps' in result.output
        assert '--timestep' in result.output

    def test_md_nvt_with_mock(self, runner, temp_structure_file):
        """Test NVT MD command with mocked calculator"""
        from mace_inference.cli import main
        from ase.build import bulk
        
        mock_atoms = bulk('Cu', 'fcc', a=3.6)
        
        with patch('mace_inference.cli.MACEInference') as MockCalc:
            mock_instance = MagicMock()
            mock_instance.run_md.return_value = mock_atoms
            MockCalc.return_value = mock_instance
            
            result = runner.invoke(main, [
                'md', str(temp_structure_file),
                '--ensemble', 'nvt',
                '--temp', '300',
                '--steps', '10'
            ])
            
            assert result.exit_code == 0


class TestPhononCommand:
    """Test phonon command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_phonon_help(self, runner):
        """Test phonon command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['phonon', '--help'])
        assert result.exit_code == 0
        assert '--supercell' in result.output
        assert '--mesh' in result.output


class TestBulkModulusCommand:
    """Test bulk-modulus command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_bulk_modulus_help(self, runner):
        """Test bulk-modulus command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['bulk-modulus', '--help'])
        assert result.exit_code == 0
        assert '--points' in result.output

    def test_bulk_modulus_with_mock(self, runner, temp_structure_file):
        """Test bulk-modulus command with mocked calculator"""
        from mace_inference.cli import main
        
        mock_result = {
            'v0': 100.0,
            'e0': -10.0,
            'B_GPa': 150.0,
            'volumes': [95, 100, 105],
            'energies': [-9.5, -10.0, -9.5]
        }
        
        with patch('mace_inference.cli.MACEInference') as MockCalc:
            mock_instance = MagicMock()
            mock_instance.bulk_modulus.return_value = mock_result
            MockCalc.return_value = mock_instance
            
            result = runner.invoke(main, ['bulk-modulus', str(temp_structure_file)])
            
            assert result.exit_code == 0
            assert 'Bulk Modulus' in result.output


class TestAdsorptionCommand:
    """Test adsorption command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_adsorption_help(self, runner):
        """Test adsorption command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['adsorption', '--help'])
        assert result.exit_code == 0
        assert '--gas' in result.output
        assert '--site' in result.output

    def test_adsorption_with_mock(self, runner, temp_structure_file):
        """Test adsorption command with mocked calculator"""
        from mace_inference.cli import main
        
        mock_result = {
            'E_ads': -0.5,
            'E_mof': -100.0,
            'E_gas': -5.0,
            'E_complex': -105.5
        }
        
        with patch('mace_inference.cli.MACEInference') as MockCalc:
            mock_instance = MagicMock()
            mock_instance.adsorption_energy.return_value = mock_result
            MockCalc.return_value = mock_instance
            
            result = runner.invoke(main, [
                'adsorption', str(temp_structure_file),
                '--gas', 'CO2',
                '--site', '5.0', '5.0', '5.0'
            ])
            
            assert result.exit_code == 0
            assert 'Adsorption' in result.output


class TestInfoCommand:
    """Test info command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_info_basic(self, runner):
        """Test info command basic output"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['info'])
        assert result.exit_code == 0
        assert 'Version' in result.output
        assert 'PyTorch' in result.output

    def test_info_verbose(self, runner):
        """Test info command with verbose flag"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['info', '--verbose'])
        assert result.exit_code == 0
        assert 'CUDA' in result.output
