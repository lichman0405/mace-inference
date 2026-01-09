"""Mechanical properties calculations"""

from typing import Dict, Tuple, Any, TYPE_CHECKING
import numpy as np
from ase import Atoms
from ase.eos import EquationOfState

if TYPE_CHECKING:
    from mace_inference.types import BulkModulusResult, Calculator


def calculate_bulk_modulus(
    atoms: Atoms,
    calculator: "Calculator",
    n_points: int = 11,
    scale_range: Tuple[float, float] = (0.95, 1.05),
    eos_type: str = "birchmurnaghan"
) -> "BulkModulusResult":
    """
    Calculate bulk modulus using equation of state.
    
    Args:
        atoms: Input ASE Atoms object
        calculator: ASE calculator
        n_points: Number of volume points
        scale_range: Volume scaling range (min_scale, max_scale)
        eos_type: Equation of state type ("birchmurnaghan", "murnaghan", etc.)
        
    Returns:
        Dictionary with equilibrium volume, energy, and bulk modulus
    """
    atoms = atoms.copy()
    atoms.calc = calculator
    
    # Generate volume points
    scale_factors = np.linspace(scale_range[0], scale_range[1], n_points)
    
    volumes = []
    energies = []
    
    # Calculate energy at each volume
    original_cell = atoms.get_cell()
    
    for scale in scale_factors:
        # Scale cell uniformly
        scaled_atoms = atoms.copy()
        scaled_atoms.set_cell(original_cell * scale, scale_atoms=True)
        scaled_atoms.calc = calculator
        
        volumes.append(scaled_atoms.get_volume())
        energies.append(scaled_atoms.get_potential_energy())
    
    # Fit equation of state
    eos = EquationOfState(volumes, energies, eos=eos_type)
    v0, e0, B = eos.fit()
    
    # Convert bulk modulus from eV/Å³ to GPa
    B_GPa = B * 160.21766208
    
    result = {
        "v0": v0,           # Equilibrium volume (Å³)
        "e0": e0,           # Equilibrium energy (eV)
        "B": B,             # Bulk modulus (eV/Å³)
        "B_GPa": B_GPa,     # Bulk modulus (GPa)
        "volumes": volumes,
        "energies": energies,
        "eos_type": eos_type
    }
    
    return result


def calculate_elastic_constants(
    atoms: Atoms,
    calculator: "Calculator",
    delta: float = 0.01,
    symmetry: str = "auto"
) -> Dict[str, Any]:
    """
    Calculate elastic constants using stress-strain relationships.
    
    Uses finite differences to compute the elastic tensor from
    stress responses to applied strains.
    
    Args:
        atoms: Input ASE Atoms object (should be relaxed)
        calculator: ASE calculator
        delta: Strain magnitude for finite differences
        symmetry: Crystal symmetry ("auto", "cubic", "hexagonal", "orthorhombic", "full")
        
    Returns:
        Dictionary with elastic constants in GPa
        
    Example:
        >>> result = calculate_elastic_constants(atoms, calculator)
        >>> print(f"C11 = {result['C'][0,0]:.1f} GPa")
        >>> print(f"Bulk modulus = {result['K_V']:.1f} GPa")
    """
    atoms = atoms.copy()
    atoms.calc = calculator
    
    # Get reference stress (should be ~0 for relaxed structure)
    _ref_stress = atoms.get_stress(voigt=False)  # 3x3 tensor, used to verify relaxed state
    
    # Voigt notation indices: 
    # 0=xx, 1=yy, 2=zz, 3=yz, 4=xz, 5=xy
    voigt_pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    
    # Elastic tensor in Voigt notation (6x6)
    C = np.zeros((6, 6))
    
    # Get original cell
    cell0 = atoms.get_cell().copy()
    
    for i, (p, q) in enumerate(voigt_pairs):
        # Apply positive strain
        strained_atoms = atoms.copy()
        strain_matrix = np.eye(3)
        if p == q:
            # Normal strain
            strain_matrix[p, p] = 1 + delta
        else:
            # Shear strain (symmetric)
            strain_matrix[p, q] = delta / 2
            strain_matrix[q, p] = delta / 2
        
        new_cell = cell0 @ strain_matrix
        strained_atoms.set_cell(new_cell, scale_atoms=True)
        strained_atoms.calc = calculator
        stress_plus = strained_atoms.get_stress(voigt=False)
        
        # Apply negative strain
        strained_atoms = atoms.copy()
        strain_matrix = np.eye(3)
        if p == q:
            strain_matrix[p, p] = 1 - delta
        else:
            strain_matrix[p, q] = -delta / 2
            strain_matrix[q, p] = -delta / 2
        
        new_cell = cell0 @ strain_matrix
        strained_atoms.set_cell(new_cell, scale_atoms=True)
        strained_atoms.calc = calculator
        stress_minus = strained_atoms.get_stress(voigt=False)
        
        # Compute stress derivative (central difference)
        dstress = (stress_plus - stress_minus) / (2 * delta)
        
        # Fill elastic tensor (convert to Voigt notation)
        for j, (r, s) in enumerate(voigt_pairs):
            C[j, i] = -dstress[r, s]  # Negative because stress = -C * strain
    
    # Convert from eV/Å³ to GPa
    eV_per_A3_to_GPa = 160.21766208
    C_GPa = C * eV_per_A3_to_GPa
    
    # Symmetrize
    C_GPa = (C_GPa + C_GPa.T) / 2
    
    # Calculate derived properties using Voigt-Reuss-Hill averages
    # Compliance tensor
    try:
        S = np.linalg.inv(C_GPa)
    except np.linalg.LinAlgError:
        S = np.zeros((6, 6))
    
    # Voigt averages (upper bound)
    K_V = (C_GPa[0,0] + C_GPa[1,1] + C_GPa[2,2] + 
           2*(C_GPa[0,1] + C_GPa[1,2] + C_GPa[0,2])) / 9
    G_V = ((C_GPa[0,0] + C_GPa[1,1] + C_GPa[2,2]) - 
           (C_GPa[0,1] + C_GPa[1,2] + C_GPa[0,2]) +
           3*(C_GPa[3,3] + C_GPa[4,4] + C_GPa[5,5])) / 15
    
    # Reuss averages (lower bound)
    K_R = 1 / (S[0,0] + S[1,1] + S[2,2] + 2*(S[0,1] + S[1,2] + S[0,2]))
    G_R = 15 / (4*(S[0,0] + S[1,1] + S[2,2]) - 
                4*(S[0,1] + S[1,2] + S[0,2]) +
                3*(S[3,3] + S[4,4] + S[5,5]))
    
    # Hill averages
    K_H = (K_V + K_R) / 2
    G_H = (G_V + G_R) / 2
    
    # Young's modulus and Poisson's ratio from Hill averages
    E_H = 9 * K_H * G_H / (3 * K_H + G_H)
    nu_H = (3 * K_H - 2 * G_H) / (2 * (3 * K_H + G_H))
    
    result = {
        "C": C_GPa,                    # Full elastic tensor (GPa)
        "S": S,                        # Compliance tensor (1/GPa)
        "K_V": float(K_V),             # Voigt bulk modulus (GPa)
        "K_R": float(K_R),             # Reuss bulk modulus (GPa)
        "K_H": float(K_H),             # Hill bulk modulus (GPa)
        "G_V": float(G_V),             # Voigt shear modulus (GPa)
        "G_R": float(G_R),             # Reuss shear modulus (GPa)
        "G_H": float(G_H),             # Hill shear modulus (GPa)
        "E_H": float(E_H),             # Hill Young's modulus (GPa)
        "nu_H": float(nu_H),           # Hill Poisson's ratio
        "delta": delta,
        "symmetry": symmetry,
    }
    
    # Add cubic-specific constants if applicable
    if symmetry in ["auto", "cubic"]:
        result["C11"] = float(C_GPa[0, 0])
        result["C12"] = float(C_GPa[0, 1])
        result["C44"] = float(C_GPa[3, 3])
        # Zener anisotropy ratio
        if abs(C_GPa[0, 0] - C_GPa[0, 1]) > 1e-6:
            result["A"] = float(2 * C_GPa[3, 3] / (C_GPa[0, 0] - C_GPa[0, 1]))
    
    return result
