"""Adsorption energy and coordination analysis"""

from typing import Union, List, Optional, Dict, Any, TYPE_CHECKING
import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.optimize import LBFGS

if TYPE_CHECKING:
    from mace_inference.types import AdsorptionResult, CoordinationResult, Calculator


def calculate_adsorption_energy(
    framework: Atoms,
    adsorbate: Union[str, Atoms],
    site_position: List[float],
    calculator: "Calculator",
    optimize: bool = True,
    fmax: float = 0.05,
    fix_framework: bool = True
) -> "AdsorptionResult":
    """
    Calculate gas adsorption energy in framework (MOF, zeolite, etc.).
    
    E_ads = E(framework+adsorbate) - E(framework) - E(adsorbate)
    
    Args:
        framework: Framework structure (MOF, zeolite, etc.) as ASE Atoms
        adsorbate: Adsorbate molecule name (e.g., "CO2") or Atoms object
        site_position: Adsorption site position [x, y, z]
        calculator: ASE calculator
        optimize: Whether to optimize the adsorption complex
        fmax: Force convergence for optimization
        fix_framework: Whether to fix framework atoms during optimization
        
    Returns:
        Dictionary with adsorption energy and structures
    """
    # 1. Calculate framework energy
    framework_copy = framework.copy()
    framework_copy.calc = calculator
    E_framework = framework_copy.get_potential_energy()
    
    # 2. Get adsorbate molecule
    if isinstance(adsorbate, str):
        ads = molecule(adsorbate)
    elif isinstance(adsorbate, Atoms):
        ads = adsorbate.copy()
    else:
        raise TypeError("adsorbate must be str or Atoms object")
    
    # Calculate adsorbate energy (in vacuum)
    ads.calc = calculator
    ads.center(vacuum=10.0)  # Add vacuum around molecule
    E_ads_mol = ads.get_potential_energy()
    
    # 3. Create adsorption complex
    complex_atoms = framework_copy.copy()
    
    # Position adsorbate molecule at adsorption site
    ads_centered = ads.copy()
    ads_com = ads_centered.get_center_of_mass()
    translation = np.array(site_position) - ads_com
    ads_centered.translate(translation)
    
    # Combine structures
    complex_atoms += ads_centered
    complex_atoms.calc = calculator
    
    # 4. Optimize complex if requested
    if optimize:
        # Fix framework atoms if requested, only optimize adsorbate
        if fix_framework:
            from ase.constraints import FixAtoms
            n_framework_atoms = len(framework)
            constraint = FixAtoms(indices=range(n_framework_atoms))
            complex_atoms.set_constraint(constraint)
        
        opt = LBFGS(complex_atoms)
        opt.run(fmax=fmax)
    
    E_complex = complex_atoms.get_potential_energy()
    
    # 5. Calculate adsorption energy
    E_ads = E_complex - E_framework - E_ads_mol
    
    result = {
        "E_ads": E_ads,              # Adsorption energy (eV)
        "E_mof": E_framework,        # Framework energy (eV) - kept as E_mof for backward compatibility
        "E_gas": E_ads_mol,          # Adsorbate energy (eV) - kept as E_gas for backward compatibility
        "E_complex": E_complex,      # Complex energy (eV)
        "complex_structure": complex_atoms,
        "optimized": optimize
    }
    
    return result


def analyze_coordination(
    atoms: Atoms,
    metal_indices: Optional[List[int]] = None,
    cutoff_multiplier: float = 1.2
) -> "CoordinationResult":
    """
    Analyze coordination environment around metal centers.
    
    Args:
        atoms: ASE Atoms object
        metal_indices: Indices of metal atoms (auto-detect if None)
        cutoff_multiplier: Multiplier for natural cutoff radii
        
    Returns:
        Dictionary with coordination analysis
    """
    # Auto-detect metal atoms if not provided
    if metal_indices is None:
        # Common metal elements in MOFs
        metal_symbols = ['Cu', 'Zn', 'Zr', 'Fe', 'Ni', 'Co', 'Mn', 'Cr', 'Ti', 'V',
                        'Al', 'Mg', 'Ca', 'Sr', 'Ba', 'Cd', 'Hg', 'Pd', 'Pt']
        metal_indices = [i for i, sym in enumerate(atoms.get_chemical_symbols()) 
                        if sym in metal_symbols]
    
    if not metal_indices:
        raise ValueError("No metal atoms found. Specify metal_indices manually.")
    
    # Create neighbor list with natural cutoffs
    cutoffs = natural_cutoffs(atoms, mult=cutoff_multiplier)
    nl = NeighborList(cutoffs, skin=0.3, self_interaction=False, bothways=False)
    nl.update(atoms)
    
    coordination_data = {}
    
    for metal_idx in metal_indices:
        indices, offsets = nl.get_neighbors(metal_idx)
        
        # Calculate distances
        metal_pos = atoms.positions[metal_idx]
        neighbor_data = []
        
        for neighbor_idx, offset in zip(indices, offsets):
            neighbor_pos = atoms.positions[neighbor_idx] + offset @ atoms.cell
            distance = np.linalg.norm(neighbor_pos - metal_pos)
            
            neighbor_data.append({
                "index": int(neighbor_idx),
                "symbol": atoms[neighbor_idx].symbol,
                "distance": float(distance)
            })
        
        # Sort by distance
        neighbor_data.sort(key=lambda x: x["distance"])
        
        coordination_data[int(metal_idx)] = {
            "metal_symbol": atoms[metal_idx].symbol,
            "coordination_number": len(neighbor_data),
            "neighbors": neighbor_data,
            "average_distance": float(np.mean([n["distance"] for n in neighbor_data])) if neighbor_data else 0.0
        }
    
    result = {
        "coordination": coordination_data,
        "n_metal_centers": len(metal_indices),
        "metal_indices": metal_indices
    }
    
    return result


def find_adsorption_sites(
    atoms: Atoms,
    grid_spacing: float = 0.5,
    probe_radius: float = 1.2,
    min_distance: float = 2.0,
    energy_cutoff: float = None
) -> List[np.ndarray]:
    """
    Find potential adsorption sites in porous structure using grid-based search.
    
    This function identifies void spaces in the structure that could
    accommodate adsorbate molecules by searching for grid points that
    are sufficiently far from all atoms.
    
    Args:
        atoms: ASE Atoms object (porous structure like MOF, zeolite)
        grid_spacing: Grid spacing for site search (Å)
        probe_radius: Minimum distance from any atom to be considered a void (Å)
        min_distance: Minimum distance between identified sites (Å)
        energy_cutoff: Optional energy cutoff for filtering sites (not implemented)
        
    Returns:
        List of potential adsorption site positions [x, y, z]
        
    Example:
        >>> sites = find_adsorption_sites(mof, probe_radius=1.5)
        >>> print(f"Found {len(sites)} potential sites")
        >>> for site in sites[:5]:
        ...     print(f"  Site at {site}")
    """
    # Get cell parameters
    cell = atoms.get_cell()
    
    # Create grid points within the unit cell
    n_grid = [max(1, int(np.linalg.norm(cell[i]) / grid_spacing)) for i in range(3)]
    
    # Generate fractional coordinates
    grid_points = []
    for i in range(n_grid[0]):
        for j in range(n_grid[1]):
            for k in range(n_grid[2]):
                frac = np.array([i / n_grid[0], j / n_grid[1], k / n_grid[2]])
                grid_points.append(frac)
    
    grid_points = np.array(grid_points)
    
    # Convert to Cartesian coordinates
    cart_points = grid_points @ cell
    
    # Get atomic positions
    positions = atoms.get_positions()
    
    # Find points that are far enough from all atoms
    void_sites = []
    
    for point in cart_points:
        # Calculate minimum distance to any atom (considering PBC)
        min_dist = float('inf')
        
        for pos in positions:
            # Simple distance (could be improved with proper PBC handling)
            diff = point - pos
            # Apply minimum image convention
            frac_diff = np.linalg.solve(cell.T, diff)
            frac_diff = frac_diff - np.round(frac_diff)
            diff_pbc = frac_diff @ cell
            dist = np.linalg.norm(diff_pbc)
            min_dist = min(min_dist, dist)
        
        if min_dist >= probe_radius:
            void_sites.append(point)
    
    if not void_sites:
        return []
    
    void_sites = np.array(void_sites)
    
    # Cluster nearby sites to avoid redundancy
    selected_sites = []
    used = np.zeros(len(void_sites), dtype=bool)
    
    for i, site in enumerate(void_sites):
        if used[i]:
            continue
        
        # Mark nearby sites as used
        for j in range(i + 1, len(void_sites)):
            if used[j]:
                continue
            diff = site - void_sites[j]
            # Apply PBC
            frac_diff = np.linalg.solve(cell.T, diff)
            frac_diff = frac_diff - np.round(frac_diff)
            diff_pbc = frac_diff @ cell
            if np.linalg.norm(diff_pbc) < min_distance:
                used[j] = True
        
        selected_sites.append(site)
        used[i] = True
    
    return selected_sites
