"""Test installation and basic functionality"""

import pytest


class TestPackageInstall:
    """Test package installation and imports"""

    def test_import_package(self):
        """Test that main package can be imported"""
        import mace_inference
        assert hasattr(mace_inference, "__version__")
        assert isinstance(mace_inference.__version__, str)

    def test_import_mace_inference_class(self):
        """Test that MACEInference class can be imported"""
        from mace_inference import MACEInference
        assert MACEInference is not None

    def test_import_get_device(self):
        """Test that get_device can be imported"""
        from mace_inference import get_device
        assert callable(get_device)

    def test_version_format(self):
        """Test that version string has valid format"""
        import mace_inference
        version = mace_inference.__version__
        # Version should contain at least one dot (e.g., "0.1.0")
        assert "." in version or "dev" in version


class TestDeviceUtils:
    """Test device utilities from package level"""

    def test_get_device_auto(self):
        """Test auto device detection"""
        from mace_inference.utils import get_device
        device = get_device("auto")
        assert device in ["cpu", "cuda"]

    def test_get_device_cpu(self):
        """Test CPU device selection"""
        from mace_inference.utils import get_device
        device = get_device("cpu")
        assert device == "cpu"

    def test_get_device_info(self):
        """Test device info retrieval"""
        from mace_inference.utils import get_device_info
        info = get_device_info()
        assert "cuda_available" in info
        assert "pytorch_version" in info
        assert isinstance(info["cuda_available"], bool)


class TestIOUtils:
    """Test I/O utilities"""

    def test_create_supercell(self):
        """Test supercell creation"""
        from ase.build import bulk
        from mace_inference.utils import create_supercell

        atoms = bulk('Cu', 'fcc', a=3.6)
        supercell = create_supercell(atoms, 2)
        assert len(supercell) == len(atoms) * 8

    def test_parse_structure_input_atoms(self):
        """Test parsing ASE Atoms input"""
        from ase.build import bulk
        from mace_inference.utils import parse_structure_input

        atoms = bulk('Cu', 'fcc', a=3.6)
        result = parse_structure_input(atoms)
        assert len(result) == len(atoms)

    def test_parse_structure_input_invalid(self):
        """Test parsing invalid input raises error"""
        from mace_inference.utils import parse_structure_input

        with pytest.raises(TypeError):
            parse_structure_input(12345)


class TestCoreInit:
    """Test core MACEInference initialization"""

    def test_init_requires_mace_torch(self):
        """Test that initialization imports mace-torch"""
        from mace_inference import MACEInference

        try:
            calc = MACEInference(model="medium", device="cpu")
            assert calc.device == "cpu"
            assert calc.model_name == "medium"
        except ImportError as e:
            # Expected if mace-torch is not installed
            assert "mace" in str(e).lower()
            pytest.skip("mace-torch not installed")

    def test_invalid_model_path(self):
        """Test that invalid model path raises error"""
        from mace_inference import MACEInference

        try:
            with pytest.raises(ValueError, match="Invalid model"):
                MACEInference(model="nonexistent_model_xyz", device="cpu")
        except ImportError:
            pytest.skip("mace-torch not installed")


class TestCLI:
    """Test CLI availability"""

    def test_cli_import(self):
        """Test that CLI module can be imported"""
        from mace_inference import cli
        assert hasattr(cli, "main")

    def test_cli_help(self):
        """Test CLI help option"""
        from click.testing import CliRunner
        from mace_inference.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "MACE Inference" in result.output

    def test_cli_version(self):
        """Test CLI version option"""
        from click.testing import CliRunner
        from mace_inference.cli import main
        import mace_inference

        runner = CliRunner()
        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0
        # Version should be in output
        assert mace_inference.__version__ in result.output or "version" in result.output.lower()
