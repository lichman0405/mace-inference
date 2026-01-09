# Example 4: Adsorption Study

This example demonstrates gas adsorption energy calculations in metal-organic frameworks.

## Source Code

See [`examples/04_adsorption_study.py`](https://github.com/lichman0405/mace-inference/blob/main/examples/04_adsorption_study.py)

## What You'll Learn

- Setting up adsorption calculations
- Placing adsorbates at specific sites
- Calculating binding energies
- Comparing different gases

## Code Walkthrough

### Setup

```python
from mace_inference import MACEInference
from ase.io import read
from ase.build import molecule

calc = MACEInference(model="medium", device="auto")
framework = read("structures/cu_paddlewheel.cif")
```

### Create Adsorbate

```python
co2 = molecule("CO2")
print(f"Adsorbate: {co2.get_chemical_formula()}")
```

### Define Adsorption Site

```python
# Choose a site in the framework
# This could be near metal centers, in pores, etc.
site = [5.0, 5.0, 5.0]  # Fractional or Cartesian coordinates
```

### Calculate Adsorption Energy

```python
result = calc.adsorption_energy(
    framework=framework,
    adsorbate=co2,
    site_position=site,
    optimize=True,
    fmax=0.05,
    fix_framework=True
)

print(f"E(framework): {result['E_mof']:.4f} eV")
print(f"E(CO2): {result['E_gas']:.4f} eV")
print(f"E(complex): {result['E_complex']:.4f} eV")
print(f"E_ads: {result['E_ads']:.4f} eV")
```

### Convert Units

```python
E_kJ = result['E_ads'] * 96.485  # kJ/mol
print(f"Adsorption energy: {E_kJ:.1f} kJ/mol")
```

### Compare Multiple Gases

```python
gases = ["CO2", "H2O", "N2", "CH4"]

print(f"{'Gas':<8} {'E_ads (eV)':<12} {'E_ads (kJ/mol)':<12}")
print("-" * 35)

for gas_name in gases:
    gas = molecule(gas_name)
    result = calc.adsorption_energy(
        framework=framework,
        adsorbate=gas,
        site_position=site
    )
    E_kJ = result['E_ads'] * 96.485
    print(f"{gas_name:<8} {result['E_ads']:<12.4f} {E_kJ:<12.1f}")
```

## Expected Output

```
============================================================
Example 4: Adsorption Study
============================================================

1. Loading framework...
   Formula: C8H4Cu2O8
   Atoms: 22

2. Single CO2 adsorption...
   E(framework): 283.6685 eV
   E(CO2): -22.8234 eV
   E(complex): 260.3245 eV
   E_ads: -0.5206 eV (-50.2 kJ/mol)

3. Multiple gas comparison...
   Gas      E_ads (eV)   E_ads (kJ/mol)
   -----------------------------------
   CO2      -0.5206      -50.2
   H2O      -0.6821      -65.8
   N2       -0.2134      -20.6
   CH4      -0.3567      -34.4

   Binding strength: H2O > CO2 > CH4 > N2

✅ Example 4 completed successfully!
```

## Key Points

1. **Negative E_ads**: Favorable (exothermic) adsorption
2. **Site Selection**: Results depend strongly on placement
3. **Fixed Framework**: Faster but less accurate; relax for production
4. **D3 Correction**: Important for physisorption - consider enabling
