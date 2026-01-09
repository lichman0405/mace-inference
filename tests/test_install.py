"""
Tests for package installation and imports.

Tests basic package functionality without heavy computation.
"""

import pytest


class TestPackageInstall:
    """Test package installation and imports"""

    def test_import_package(self):
        """Test that package can be imported"""
        import mace_inference
        
        assert mace_inference is not None

    def test_import_mace_inference_class(self):
        """Test MACEInference class import"""
        from mace_inference import MACEInference
        
        assert MACEInference is not None

    def test_import_get_device(self):
        """Test get_device function import"""
        from mace_inference import get_device
        
        assert callable(get_device)

    def test_version_format(self):
        """Test version string format"""
        from mace_inference import __version__
        
        assert isinstance(__version__, str)
        # Version should be in semver format
        parts = __version__.split('.')
        assert len(parts) >= 2


class TestDeviceUtils:
    """Test device utilities"""

    def test_get_device_auto(self):
        """Test automatic device detection"""
        from mace_inference import get_device
        
        device = get_device("auto")
        assert device in ["cpu", "cuda"]

    def test_get_device_cpu(self):
        """Test CPU device selection"""
        from mace_inference import get_device
        
        device = get_device("cpu")
        assert device == "cpu"

    def test_get_device_info(self):
        """Test device info function"""
        from mace_inference.utils import get_device_info
        
        info = get_device_info()
        
        assert isinstance(info, dict)
        assert "pytorch_version" in info
        assert "cuda_available" in info


class TestIOUtils:
    """Test I/O utilities"""

    def test_create_supercell(self, cu_bulk):
        """Test supercell creation"""
        from mace_inference.utils import create_supercell
        
        supercell = create_supercell(cu_bulk, 2)
        
        assert len(supercell) == len(cu_bulk) * 8

    def test_parse_structure_input_atoms(self, cu_bulk):
        """Test parsing Atoms object"""
        from mace_inference.utils import parse_structure_input
        
        result = parse_structure_input(cu_bulk)
        
        assert len(result) == len(cu_bulk)

    def test_parse_structure_input_invalid(self):
        """Test parsing invalid path"""
        from mace_inference.utils import parse_structure_input
        
        with pytest.raises((FileNotFoundError, ValueError, TypeError)):
            parse_structure_input("/nonexistent/path/structure.xyz")


class TestCoreInit:
    """Test MACEInference initialization"""

    def test_init_small_model(self):
        """Test initialization with small model"""
        from mace_inference import MACEInference
        
        calc = MACEInference(model="small", device="cpu")
        
        assert calc.model_name == "small"
        assert calc.device == "cpu"
        assert calc.calculator is not None

    def test_init_with_d3(self):
        """Test initialization with D3 correction"""
        from mace_inference import MACEInference
        
        calc = MACEInference(model="small", device="cpu", enable_d3=True)
        
        assert calc.enable_d3 is True

    def test_invalid_model_path(self):
        """Test invalid model path raises error"""
        from mace_inference import MACEInference
        
        with pytest.raises((FileNotFoundError, ValueError, RuntimeError)):
            MACEInference(model="/nonexistent/model.pt", device="cpu")


class TestCLI:
    """Test CLI functionality"""

    def test_cli_import(self):
        """Test CLI module import"""
        from mace_inference.cli import main
        
        assert main is not None

    def test_cli_help(self):
        """Test CLI help command"""
        from click.testing import CliRunner
        from mace_inference.cli import main
        
        runner = CliRunner()
        result = runner.invoke(main, ['--help'])
        
        assert result.exit_code == 0
        assert 'MACE Inference' in result.output

    def test_cli_version(self):
        """Test CLI version command"""
        from click.testing import CliRunner
        from mace_inference.cli import main
        
        runner = CliRunner()
        result = runner.invoke(main, ['--version'])
        
        assert result.exit_code == 0
