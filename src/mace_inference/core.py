"""Core MACE Inference class"""

from typing import Optional, Union, List, Literal, Dict, Any, Callable
from pathlib import Path
import numpy as np
from ase import Atoms

from mace_inference.utils.device import get_device, validate_device
from mace_inference.utils.d3_correction import create_combined_calculator
from mace_inference.utils.io import parse_structure_input, save_structure, create_supercell
from mace_inference.tasks import (
    single_point_energy,
    optimize_structure,
    run_nvt_md,
    run_npt_md,
    calculate_phonon,
    calculate_thermal_properties,
    calculate_bulk_modulus,
    calculate_adsorption_energy,
    analyze_coordination,
)


class MACEInference:
    """
    High-level interface for MACE machine learning force field inference.
    
    This class provides a unified API for common inference tasks including:
    - Single-point energy calculations
    - Structure optimization
    - Molecular dynamics (NVT/NPT)
    - Phonon calculations
    - Mechanical properties
    - Adsorption energies
    
    Args:
        model: MACE model name ("small", "medium", "large") or path to custom model
        device: Compute device ("auto", "cpu", or "cuda")
        enable_d3: Enable DFT-D3 dispersion correction
        d3_damping: D3 damping function ("bj", "zero", "zerom", "bjm")
        d3_xc: D3 exchange-correlation functional (e.g., "pbe", "b3lyp")
        default_dtype: Default data type ("float32" or "float64")
        
    Examples:
        >>> # Basic usage
        >>> calc = MACEInference(model="medium", device="auto")
        >>> result = calc.single_point("structure.cif")
        
        >>> # With D3 correction
        >>> calc = MACEInference(model="medium", enable_d3=True)
        
        >>> # GPU acceleration
        >>> calc = MACEInference(model="large", device="cuda")
    """
    
    def __init__(
        self,
        model: str = "medium",
        device: Literal["auto", "cpu", "cuda"] = "auto",
        enable_d3: bool = False,
        d3_damping: str = "bj",
        d3_xc: str = "pbe",
        default_dtype: str = "float64",
    ):
        # Device setup
        self.device = get_device(device)
        validate_device(self.device)
        
        # Model configuration
        self.model_name = model
        self.enable_d3 = enable_d3
        self.d3_damping = d3_damping
        self.d3_xc = d3_xc
        self.default_dtype = default_dtype
        
        # Initialize calculator
        self.calculator = self._create_calculator()
        
    def _create_calculator(self):
        """Create MACE calculator with optional D3 correction."""
        try:
            from mace.calculators import mace_mp, mace_off, MACECalculator
        except ImportError:
            raise ImportError(
                "mace-torch is required. Install with: pip install mace-torch"
            )
        
        # Create base MACE calculator
        if self.model_name in ["small", "medium", "large"]:
            # Use pre-trained MACE-MP models
            mace_calc = mace_mp(
                model=self.model_name,
                device=self.device,
                default_dtype=self.default_dtype
            )
        elif Path(self.model_name).exists():
            # Load custom model from file
            mace_calc = MACECalculator(
                model_paths=self.model_name,
                device=self.device,
                default_dtype=self.default_dtype
            )
        else:
            raise ValueError(
                f"Invalid model: {self.model_name}. "
                "Use 'small', 'medium', 'large', or path to custom model."
            )
        
        # Add D3 correction if enabled
        return create_combined_calculator(
            mace_calc,
            enable_d3=self.enable_d3,
            d3_device=self.device,
            d3_damping=self.d3_damping,
            d3_xc=self.d3_xc
        )
    
    def single_point(
        self,
        structure: Union[str, Path, Atoms],
        properties: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Perform single-point energy calculation.
        
        Args:
            structure: Structure file path or Atoms object
            properties: Properties to calculate (default: ["energy", "forces", "stress"])
            
        Returns:
            Dictionary with calculated properties
            
        Examples:
            >>> result = calc.single_point("structure.cif")
            >>> print(f"Energy: {result['energy']:.4f} eV")
            >>> print(f"Max force: {result['forces'].max():.4f} eV/A")
        """
        atoms = parse_structure_input(structure)
        return single_point_energy(atoms, self.calculator, properties)
    
    def optimize(
        self,
        structure: Union[str, Path, Atoms],
        fmax: float = 0.05,
        steps: int = 500,
        optimizer: str = "LBFGS",
        optimize_cell: bool = False,
        trajectory: Optional[str] = None,
        logfile: Optional[str] = None,
        output: Optional[str] = None
    ) -> Atoms:
        """
        Optimize atomic structure.
        
        Args:
            structure: Input structure
            fmax: Force convergence criterion (eV/A)
            steps: Maximum optimization steps
            optimizer: Optimization algorithm ("LBFGS", "BFGS", "FIRE")
            optimize_cell: Whether to optimize cell parameters
            trajectory: Trajectory file path
            logfile: Optimization log file path
            output: Output structure file path
            
        Returns:
            Optimized Atoms object
            
        Examples:
            >>> optimized = calc.optimize("structure.cif", fmax=0.05)
            >>> optimized = calc.optimize("structure.cif", optimize_cell=True, output="opt.cif")
        """
        atoms = parse_structure_input(structure)
        optimized_atoms = optimize_structure(
            atoms=atoms,
            calculator=self.calculator,
            fmax=fmax,
            steps=steps,
            optimizer=optimizer,
            optimize_cell=optimize_cell,
            trajectory=trajectory,
            logfile=logfile
        )
        
        if output:
            save_structure(optimized_atoms, output)
        
        return optimized_atoms
    
    def run_md(
        self,
        structure: Union[str, Path, Atoms],
        ensemble: Literal["nvt", "npt"] = "nvt",
        temperature_K: float = 300,
        steps: int = 1000,
        timestep: float = 1.0,
        pressure_GPa: Optional[float] = None,
        taut: Optional[float] = None,
        taup: Optional[float] = None,
        trajectory: Optional[str] = None,
        logfile: Optional[str] = None,
        log_interval: int = 100,
        progress_callback: Optional[Callable[[int, int], None]] = None
    ) -> Atoms:
        """
        Run molecular dynamics simulation.
        
        Args:
            structure: Input structure
            ensemble: MD ensemble ("nvt" or "npt")
            temperature_K: Target temperature (K)
            steps: Number of MD steps
            timestep: Time step (fs)
            pressure_GPa: Target pressure for NPT (GPa)
            taut: Temperature coupling time (fs)
            taup: Pressure coupling time (fs)
            trajectory: Trajectory file path
            logfile: MD log file path
            log_interval: Logging interval (steps)
            progress_callback: Optional callback function(current_step, total_steps)
            
        Returns:
            Final Atoms object
            
        Examples:
            >>> # NVT simulation
            >>> final = calc.run_md("structure.cif", ensemble="nvt", temperature_K=300, steps=10000)
            
            >>> # NPT simulation
            >>> final = calc.run_md("structure.cif", ensemble="npt", temperature_K=300, 
            ...                     pressure_GPa=1.0, steps=10000)
            
            >>> # With progress callback
            >>> def progress(current, total):
            ...     print(f"MD: {current}/{total} steps")
            >>> final = calc.run_md("structure.cif", steps=1000, progress_callback=progress)
        """
        atoms = parse_structure_input(structure)
        
        if ensemble == "nvt":
            return run_nvt_md(
                atoms=atoms,
                calculator=self.calculator,
                temperature_K=temperature_K,
                timestep=timestep,
                steps=steps,
                trajectory=trajectory,
                logfile=logfile,
                log_interval=log_interval,
                taut=taut,
                progress_callback=progress_callback
            )
        elif ensemble == "npt":
            if pressure_GPa is None:
                pressure_GPa = 0.0  # Default to 0 GPa (1 atm = 0 GPa approx)
            
            return run_npt_md(
                atoms=atoms,
                calculator=self.calculator,
                temperature_K=temperature_K,
                pressure_GPa=pressure_GPa,
                timestep=timestep,
                steps=steps,
                trajectory=trajectory,
                logfile=logfile,
                log_interval=log_interval,
                taut=taut,
                taup=taup,
                progress_callback=progress_callback
            )
        else:
            raise ValueError(f"Invalid ensemble: {ensemble}. Use 'nvt' or 'npt'")
    
    def phonon(
        self,
        structure: Union[str, Path, Atoms],
        supercell_matrix: Union[List[int], int] = 2,
        displacement: float = 0.01,
        mesh: List[int] = None,
        temperature_range: Optional[tuple] = None,
        output_dir: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Calculate phonon properties.
        
        Args:
            structure: Input structure
            supercell_matrix: Supercell size for phonon calculation
            displacement: Atomic displacement distance (A)
            mesh: k-point mesh for phonon DOS (default: [20, 20, 20])
            temperature_range: Temperature range for thermal properties (min, max, step)
            output_dir: Directory for output files
            
        Returns:
            Dictionary with phonon properties
            
        Examples:
            >>> result = calc.phonon("structure.cif", supercell_matrix=[2, 2, 2])
            >>> print(f"Frequencies: {result['frequencies']}")
        """
        if mesh is None:
            mesh = [20, 20, 20]
            
        atoms = parse_structure_input(structure)
        
        # Create supercell matrix
        if isinstance(supercell_matrix, int):
            sc_matrix = [[supercell_matrix, 0, 0],
                        [0, supercell_matrix, 0],
                        [0, 0, supercell_matrix]]
        else:
            sc_matrix = [[supercell_matrix[0], 0, 0],
                        [0, supercell_matrix[1], 0],
                        [0, 0, supercell_matrix[2]]]
        
        result = calculate_phonon(
            atoms=atoms,
            calculator=self.calculator,
            supercell_matrix=sc_matrix,
            displacement=displacement,
            mesh=mesh
        )
        
        # Calculate thermal properties if requested
        if temperature_range is not None:
            thermal = calculate_thermal_properties(
                result["phonon"],
                temperature_range=temperature_range
            )
            result["thermal"] = thermal
        
        return result
    
    def bulk_modulus(
        self,
        structure: Union[str, Path, Atoms],
        strain_range: float = 0.05,
        n_points: int = 7,
        optimize_first: bool = True,
        fmax: float = 0.01
    ) -> Dict[str, Any]:
        """
        Calculate bulk modulus.
        
        Args:
            structure: Input structure
            strain_range: Maximum strain for volume sweep
            n_points: Number of strain points
            optimize_first: Whether to optimize structure first
            fmax: Force criterion for optimization
            
        Returns:
            Dictionary with bulk modulus results
            
        Examples:
            >>> result = calc.bulk_modulus("structure.cif")
            >>> print(f"Bulk modulus: {result['bulk_modulus_GPa']:.1f} GPa")
        """
        atoms = parse_structure_input(structure)
        
        # Optimize first if requested
        if optimize_first:
            atoms = optimize_structure(
                atoms=atoms,
                calculator=self.calculator,
                fmax=fmax,
                optimize_cell=True
            )
        
        return calculate_bulk_modulus(
            atoms=atoms,
            calculator=self.calculator,
            strain_range=strain_range,
            n_points=n_points
        )
    
    def adsorption_energy(
        self,
        framework: Union[str, Path, Atoms],
        adsorbate: Union[str, Atoms],
        site_position: List[float],
        optimize: bool = True,
        fmax: float = 0.05,
        fix_framework: bool = True
    ) -> Dict[str, Any]:
        """
        Calculate adsorption energy.
        
        Args:
            framework: Framework structure (e.g., MOF, zeolite)
            adsorbate: Adsorbate molecule (formula or Atoms)
            site_position: Adsorption site position [x, y, z]
            optimize: Whether to optimize the combined structure
            fmax: Force criterion for optimization
            fix_framework: Whether to fix framework atoms during optimization
            
        Returns:
            Dictionary with adsorption energy results
            
        Examples:
            >>> result = calc.adsorption_energy(
            ...     framework="mof.cif",
            ...     adsorbate="CO2",
            ...     site_position=[5.0, 5.0, 5.0]
            ... )
            >>> print(f"Adsorption energy: {result['adsorption_energy']:.3f} eV")
        """
        framework_atoms = parse_structure_input(framework)
        
        return calculate_adsorption_energy(
            framework=framework_atoms,
            adsorbate=adsorbate,
            calculator=self.calculator,
            site_position=site_position,
            optimize=optimize,
            fmax=fmax,
            fix_framework=fix_framework
        )
    
    def coordination(
        self,
        structure: Union[str, Path, Atoms],
        metal_indices: Optional[List[int]] = None,
        cutoff_multiplier: float = 1.3
    ) -> Dict[str, Any]:
        """
        Analyze coordination environment of metal atoms.
        
        Args:
            structure: Input structure
            metal_indices: Indices of metal atoms to analyze (auto-detect if None)
            cutoff_multiplier: Multiplier for covalent radii cutoff
            
        Returns:
            Dictionary with coordination analysis results
            
        Examples:
            >>> result = calc.coordination("mof.cif")
            >>> for metal_idx, info in result["coordination"].items():
            ...     print(f"Metal {metal_idx}: CN = {info['coordination_number']}")
        """
        atoms = parse_structure_input(structure)
        
        return analyze_coordination(
            atoms=atoms,
            metal_indices=metal_indices,
            cutoff_multiplier=cutoff_multiplier
        )
    
    # =====================================================================
    # Batch Processing Methods
    # =====================================================================
    
    def batch_single_point(
        self,
        structures: List[Union[str, Path, Atoms]],
        properties: Optional[List[str]] = None,
        progress_callback: Optional[Callable[[int, int], None]] = None
    ) -> List[Dict[str, Any]]:
        """
        Perform single-point calculations on multiple structures.
        
        More efficient than calling single_point() in a loop as it
        reuses the calculator and can potentially batch operations.
        
        Args:
            structures: List of structure files or Atoms objects
            properties: Properties to calculate (default: ["energy", "forces", "stress"])
            progress_callback: Optional callback function(current, total) for progress updates
            
        Returns:
            List of result dictionaries
            
        Examples:
            >>> structures = ["struct1.cif", "struct2.cif", "struct3.cif"]
            >>> results = calc.batch_single_point(structures)
            >>> for i, result in enumerate(results):
            ...     print(f"Structure {i}: E = {result['energy']:.4f} eV")
            
            >>> # With progress callback
            >>> def progress(current, total):
            ...     print(f"Processing {current}/{total}")
            >>> results = calc.batch_single_point(structures, progress_callback=progress)
        """
        results = []
        total = len(structures)
        
        for i, structure in enumerate(structures):
            if progress_callback:
                progress_callback(i + 1, total)
            
            try:
                atoms = parse_structure_input(structure)
                result = single_point_energy(atoms, self.calculator, properties)
                result["success"] = True
                result["structure_index"] = i
            except Exception as e:
                result = {
                    "success": False,
                    "error": str(e),
                    "structure_index": i
                }
            
            results.append(result)
        
        return results
    
    def batch_optimize(
        self,
        structures: List[Union[str, Path, Atoms]],
        fmax: float = 0.05,
        steps: int = 500,
        optimizer: str = "LBFGS",
        optimize_cell: bool = False,
        output_dir: Optional[Union[str, Path]] = None,
        progress_callback: Optional[Callable[[int, int], None]] = None
    ) -> List[Dict[str, Any]]:
        """
        Optimize multiple structures in batch.
        
        Args:
            structures: List of structure files or Atoms objects
            fmax: Force convergence criterion (eV/A)
            steps: Maximum optimization steps per structure
            optimizer: Optimization algorithm
            optimize_cell: Whether to optimize cell parameters
            output_dir: Directory to save optimized structures (optional)
            progress_callback: Optional callback function(current, total) for progress updates
            
        Returns:
            List of dictionaries containing optimized atoms and metadata
            
        Examples:
            >>> structures = ["struct1.cif", "struct2.cif"]
            >>> results = calc.batch_optimize(structures, output_dir="optimized/")
            >>> for result in results:
            ...     if result["success"]:
            ...         print(f"Structure {result['structure_index']}: converged={result['converged']}")
        """
        results = []
        total = len(structures)
        
        if output_dir:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)
        
        for i, structure in enumerate(structures):
            if progress_callback:
                progress_callback(i + 1, total)
            
            try:
                atoms = parse_structure_input(structure)
                
                # Determine output file if output_dir specified
                if output_dir:
                    if isinstance(structure, (str, Path)):
                        output_file = output_path / f"opt_{Path(structure).stem}.xyz"
                    else:
                        output_file = output_path / f"opt_structure_{i}.xyz"
                else:
                    output_file = None
                
                optimized = optimize_structure(
                    atoms=atoms,
                    calculator=self.calculator,
                    fmax=fmax,
                    steps=steps,
                    optimizer=optimizer,
                    optimize_cell=optimize_cell
                )
                
                if output_file:
                    save_structure(optimized, str(output_file))
                
                # Get final properties
                final_props = single_point_energy(optimized, self.calculator)
                
                result = {
                    "success": True,
                    "structure_index": i,
                    "atoms": optimized,
                    "final_energy": final_props["energy"],
                    "final_max_force": np.max(np.abs(final_props["forces"])),
                    "converged": np.max(np.abs(final_props["forces"])) < fmax,
                    "output_file": str(output_file) if output_file else None
                }
                
            except Exception as e:
                result = {
                    "success": False,
                    "error": str(e),
                    "structure_index": i,
                    "atoms": None,
                    "converged": False
                }
            
            results.append(result)
        
        return results
    
    def __repr__(self) -> str:
        return (
            f"MACEInference(model='{self.model_name}', device='{self.device}', "
            f"enable_d3={self.enable_d3})"
        )
