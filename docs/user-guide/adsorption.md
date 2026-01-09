# Adsorption Studies

Calculate gas adsorption energies in frameworks like MOFs or zeolites.

## Basic Usage

```python
from mace_inference import MACEInference
from ase.io import read
from ase.build import molecule

calc = MACEInference(model="medium", device="auto")

# Load framework and create adsorbate
framework = read("mof.cif")
adsorbate = molecule("CO2")

# Calculate adsorption energy
result = calc.adsorption_energy(
    framework=framework,
    adsorbate=adsorbate,
    site_position=[5.0, 5.0, 5.0]
)

print(f"Adsorption energy: {result['E_ads']:.4f} eV")
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `framework` | `Atoms` | required | Framework structure (MOF, zeolite, etc.) |
| `adsorbate` | `Atoms` | required | Gas molecule to adsorb |
| `site_position` | `list` | required | Position to place adsorbate center [x, y, z] |
| `optimize` | `bool` | `True` | Optimize the complex after placing |
| `fmax` | `float` | `0.05` | Force tolerance for optimization (eV/Å) |
| `fix_framework` | `bool` | `True` | Keep framework atoms fixed during optimization |

## Return Values

| Key | Type | Unit | Description |
|-----|------|------|-------------|
| `E_ads` | `float` | eV | Adsorption energy (E_complex - E_mof - E_gas) |
| `E_mof` | `float` | eV | Energy of bare framework |
| `E_gas` | `float` | eV | Energy of isolated gas molecule |
| `E_complex` | `float` | eV | Energy of framework + adsorbate |
| `complex_structure` | `Atoms` | - | Final complex structure |
| `optimized` | `bool` | - | Whether structure was optimized |

## The Adsorption Energy

The adsorption energy is calculated as:

$$E_{ads} = E_{complex} - E_{framework} - E_{gas}$$

Where:

- $E_{complex}$ = Energy of framework with adsorbed molecule
- $E_{framework}$ = Energy of empty framework
- $E_{gas}$ = Energy of isolated gas molecule

**Negative values** indicate favorable (exothermic) adsorption.

## Examples

### CO₂ Adsorption in MOF

```python
from ase.build import molecule

framework = read("HKUST-1.cif")
co2 = molecule("CO2")

result = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=[10.0, 10.0, 10.0],
    optimize=True,
    fmax=0.05
)

print(f"E(MOF): {result['E_mof']:.4f} eV")
print(f"E(CO2): {result['E_gas']:.4f} eV")
print(f"E(complex): {result['E_complex']:.4f} eV")
print(f"E_ads: {result['E_ads']:.4f} eV ({result['E_ads']*1000:.1f} meV)")
print(f"E_ads: {result['E_ads']*96.485:.1f} kJ/mol")
```

### Different Adsorbates

```python
from ase.build import molecule

framework = read("zeolite.cif")
gases = ["H2O", "CO2", "N2", "CH4"]
site = [5.0, 5.0, 5.0]

print(f"{'Gas':<10} {'E_ads (eV)':<15} {'E_ads (kJ/mol)':<15}")
print("-" * 40)

for gas_name in gases:
    gas = molecule(gas_name)
    result = calc.adsorption_energy(
        framework=framework,
        adsorbate=gas,
        site_position=site
    )
    E_kJmol = result['E_ads'] * 96.485
    print(f"{gas_name:<10} {result['E_ads']:<15.4f} {E_kJmol:<15.1f}")
```

### Finding Optimal Adsorption Site

Screen multiple sites to find the strongest binding:

```python
import numpy as np

framework = read("mof.cif")
co2 = molecule("CO2")

# Define candidate sites
sites = [
    [5.0, 5.0, 5.0],
    [10.0, 5.0, 5.0],
    [5.0, 10.0, 5.0],
    [7.5, 7.5, 7.5],
]

best_E = 0
best_site = None

for site in sites:
    result = calc.adsorption_energy(
        framework=framework,
        adsorbate=co2,
        site_position=site
    )
    print(f"Site {site}: E_ads = {result['E_ads']:.4f} eV")
    
    if result['E_ads'] < best_E:
        best_E = result['E_ads']
        best_site = site

print(f"\nBest site: {best_site} with E_ads = {best_E:.4f} eV")
```

### Fixed vs Relaxed Framework

```python
# Fixed framework (faster, for screening)
result_fixed = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=site,
    fix_framework=True
)

# Relaxed framework (more accurate)
result_relaxed = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=site,
    fix_framework=False,
    fmax=0.02
)

print(f"Fixed framework: E_ads = {result_fixed['E_ads']:.4f} eV")
print(f"Relaxed framework: E_ads = {result_relaxed['E_ads']:.4f} eV")
```

### Using D3 Correction

Dispersion is important for physisorption:

```python
# With D3 correction for better van der Waals description
calc_d3 = MACEInference(model="medium", enable_d3=True, d3_xc="pbe")

result = calc_d3.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=site
)

print(f"E_ads (with D3): {result['E_ads']:.4f} eV")
```

## Workflow Tips

### Pre-optimize Framework

```python
# Optimize framework first
opt = calc.optimize(framework, fmax=0.01)
framework_opt = opt['atoms']

# Use optimized framework for adsorption
result = calc.adsorption_energy(
    framework=framework_opt,
    adsorbate=co2,
    site_position=site
)
```

### Save Complex Structure

```python
from ase.io import write

result = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=site
)

# Save the optimized complex
write("complex_optimized.cif", result['complex'])
```

### Unit Conversions

```python
# eV to kJ/mol
E_kJmol = result['E_ads'] * 96.485

# eV to kcal/mol
E_kcalmol = result['E_ads'] * 23.061

# eV to meV
E_meV = result['E_ads'] * 1000
```

## Common Issues

### Adsorbate Too Close to Framework

If atoms overlap, optimization may fail:

1. Choose a site further from framework atoms
2. Visualize the structure first
3. Reduce `fmax` and increase optimization steps

### Slow Optimization

For large frameworks:

1. Use `fix_framework=True`
2. Use GPU acceleration
3. Reduce framework to relevant region

### Unrealistic Energies

1. Check that framework structure is reasonable
2. Ensure adsorbate is correctly oriented
3. Consider using D3 correction for physisorption
