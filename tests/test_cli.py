"""
Tests for CLI (Command Line Interface) functionality.

Uses Click's CliRunner with real MACE models.
"""

import pytest
import os
import tempfile
from click.testing import CliRunner
from ase.build import bulk


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

    def test_energy_real_calculation(self, runner, temp_structure_file):
        """Test energy command with real MACE model"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, [
            'energy', 
            temp_structure_file,
            '--model', 'small',
            '--device', 'cpu'
        ])
        
        assert result.exit_code == 0
        assert 'Energy' in result.output or 'energy' in result.output.lower()


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

    def test_optimize_real_calculation(self, runner, temp_structure_file, temp_output_dir):
        """Test optimize command with real MACE model"""
        from mace_inference.cli import main
        
        output_file = os.path.join(temp_output_dir, "optimized.xyz")
        
        result = runner.invoke(main, [
            'optimize',
            temp_structure_file,
            '--model', 'small',
            '--device', 'cpu',
            '--fmax', '0.1',
            '--steps', '5',
            '--output', output_file
        ])
        
        assert result.exit_code == 0


class TestMDCommand:
    """Test md command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_md_help(self, runner):
        """Test MD command help"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['md', '--help'])
        assert result.exit_code == 0
        assert '--ensemble' in result.output or 'ensemble' in result.output.lower()
        assert '--temp' in result.output or 'temperature' in result.output.lower()


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
        assert '--supercell' in result.output or 'supercell' in result.output.lower()


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

    def test_bulk_modulus_real_calculation(self, runner, temp_structure_file):
        """Test bulk-modulus command with real MACE model"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, [
            'bulk-modulus',
            temp_structure_file,
            '--model', 'small',
            '--device', 'cpu',
            '--points', '5'
        ])
        
        assert result.exit_code == 0
        assert 'Bulk' in result.output or 'modulus' in result.output.lower()


class TestInfoCommand:
    """Test info command"""

    @pytest.fixture
    def runner(self):
        return CliRunner()

    def test_info_basic(self, runner):
        """Test info command"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['info'])
        
        # Should succeed and show package info
        assert result.exit_code == 0

    def test_info_verbose(self, runner):
        """Test info command with verbose flag"""
        from mace_inference.cli import main
        
        result = runner.invoke(main, ['info', '--verbose'])
        
        assert result.exit_code == 0
