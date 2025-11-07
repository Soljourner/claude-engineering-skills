---
name: coolprop-db
description: "Query thermodynamic properties for 100+ fluids from CoolProp database"
category: databases
domain: fluids
complexity: basic
dependencies:
  - CoolProp
---

# CoolProp Database Skill

Query thermodynamic and transport properties for over 100 pure and pseudo-pure fluids using the CoolProp open-source thermophysical property library.

## Overview

CoolProp is a comprehensive thermophysical property database that provides:

- **Pure Fluids**: Water, air, ammonia, carbon dioxide, and 100+ others
- **Refrigerants**: R134a, R410A, R32, R404A, R407C, R507A, and many more
- **Industrial Fluids**: Methane, ethane, propane, butane, nitrogen, oxygen, argon
- **Transport Properties**: Viscosity, thermal conductivity, surface tension
- **Thermodynamic Properties**: Enthalpy, entropy, density, specific heat
- **Phase Information**: Saturation properties, two-phase behavior, critical points

CoolProp uses high-accuracy equations of state (Helmholtz energy formulations) and is validated against NIST REFPROP data.

## Installation

### Python
```bash
pip install CoolProp
```

### Verify Installation
```python
import CoolProp
print(CoolProp.__version__)
print(CoolProp.get_global_param_string("version"))
```

## Core API Functions

### PropsSI - Primary Property Query Function
```python
from CoolProp.CoolProp import PropsSI

# Syntax: PropsSI(output, input1_name, input1_value, input2_name, input2_value, fluid)
value = PropsSI('D', 'T', 298.15, 'P', 101325, 'Water')
```

### Props1SI - Single Input Properties
```python
from CoolProp.CoolProp import Props1SI

# For properties requiring only fluid name
T_crit = Props1SI('Tcrit', 'Water')  # Critical temperature
P_crit = Props1SI('Pcrit', 'Water')  # Critical pressure
```

## Property Codes (Input/Output Parameters)

### Thermodynamic Properties
| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `T` | Temperature | K | Absolute temperature |
| `P` | Pressure | Pa | Absolute pressure |
| `D` | Density | kg/m³ | Mass density |
| `H` | Enthalpy | J/kg | Specific enthalpy |
| `S` | Entropy | J/kg/K | Specific entropy |
| `U` | Internal Energy | J/kg | Specific internal energy |
| `Q` | Quality | - | Vapor mass fraction (0-1) |
| `Dmolar` | Molar Density | mol/m³ | Molar density |
| `Hmolar` | Molar Enthalpy | J/mol | Molar enthalpy |
| `Smolar` | Molar Entropy | J/mol/K | Molar entropy |

### Transport Properties
| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `V` | Viscosity | Pa·s | Dynamic viscosity |
| `L` | Thermal Conductivity | W/m/K | Thermal conductivity |
| `C` | Specific Heat (const P) | J/kg/K | Cp at constant pressure |
| `O` | Specific Heat (const V) | J/kg/K | Cv at constant volume |
| `PRANDTL` | Prandtl Number | - | Pr = μ·Cp/k |
| `I` | Surface Tension | N/m | Liquid-vapor interface |

### Phase Properties
| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `Phase` | Phase Index | - | 0=Liquid, 3=Supercritical, 5=Gas, 6=Two-phase |
| `Q` | Quality | - | 0=Saturated liquid, 1=Saturated vapor |

### Critical/Triple Point Properties (Use with Props1SI)
| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `Tcrit` | Critical Temperature | K | Critical point temperature |
| `Pcrit` | Critical Pressure | Pa | Critical point pressure |
| `Ttriple` | Triple Point Temp | K | Triple point temperature |
| `Ptriple` | Triple Point Press | Pa | Triple point pressure |
| `M` | Molar Mass | kg/mol | Molecular weight |
| `ACENTRIC` | Acentric Factor | - | Pitzer acentric factor |

## Common Fluids

### Water and Air
- `Water` - Pure water (H₂O)
- `Air` - Dry air (pseudo-pure mixture)

### Common Refrigerants
- `R134a` - HFC, common in automotive AC
- `R410A` - HFC blend, residential AC/heat pumps
- `R32` - HFC, lower GWP alternative
- `R404A` - HFC blend, commercial refrigeration
- `R407C` - HFC blend, AC systems
- `R507A` - HFC blend, low-temperature refrigeration
- `R22` - HCFC (being phased out)
- `R717` - Ammonia (NH₃)
- `R744` - Carbon dioxide (CO₂)

### Hydrocarbons
- `Methane`, `Ethane`, `Propane`, `n-Butane`, `IsoButane`
- `n-Pentane`, `Isopentane`, `n-Hexane`, `n-Heptane`, `n-Octane`
- `n-Nonane`, `n-Decane`

### Cryogenic Gases
- `Nitrogen`, `Oxygen`, `Argon`, `Helium`, `Neon`, `Hydrogen`

### Industrial Gases
- `CO2` - Carbon dioxide
- `CO` - Carbon monoxide
- `H2S` - Hydrogen sulfide
- `SO2` - Sulfur dioxide
- `Ammonia` - NH₃

## Query Examples

### Example 1: Water Properties at Standard Conditions
```python
from CoolProp.CoolProp import PropsSI

# Water at 25°C (298.15 K) and 1 atm (101325 Pa)
T = 298.15  # K
P = 101325  # Pa

density = PropsSI('D', 'T', T, 'P', P, 'Water')  # kg/m³
enthalpy = PropsSI('H', 'T', T, 'P', P, 'Water')  # J/kg
entropy = PropsSI('S', 'T', T, 'P', P, 'Water')  # J/kg/K
viscosity = PropsSI('V', 'T', T, 'P', P, 'Water')  # Pa·s
cp = PropsSI('C', 'T', T, 'P', P, 'Water')  # J/kg/K

print(f"Water at {T-273.15}°C and {P/1000:.1f} kPa:")
print(f"  Density: {density:.2f} kg/m³")
print(f"  Enthalpy: {enthalpy/1000:.2f} kJ/kg")
print(f"  Entropy: {entropy/1000:.4f} kJ/kg·K")
print(f"  Viscosity: {viscosity*1000:.4f} mPa·s")
print(f"  Cp: {cp/1000:.4f} kJ/kg·K")
```

### Example 2: Refrigerant Saturation Properties
```python
from CoolProp.CoolProp import PropsSI

# R134a saturation properties at 25°C
T_sat = 298.15  # K
fluid = 'R134a'

# Get saturation pressure at this temperature
P_sat = PropsSI('P', 'T', T_sat, 'Q', 0, fluid)  # Pa

# Saturated liquid properties (Q=0)
rho_liquid = PropsSI('D', 'T', T_sat, 'Q', 0, fluid)
h_liquid = PropsSI('H', 'T', T_sat, 'Q', 0, fluid)
s_liquid = PropsSI('S', 'T', T_sat, 'Q', 0, fluid)

# Saturated vapor properties (Q=1)
rho_vapor = PropsSI('D', 'T', T_sat, 'Q', 1, fluid)
h_vapor = PropsSI('H', 'T', T_sat, 'Q', 1, fluid)
s_vapor = PropsSI('S', 'T', T_sat, 'Q', 1, fluid)

# Latent heat
h_fg = h_vapor - h_liquid

print(f"{fluid} at {T_sat-273.15}°C:")
print(f"  Saturation pressure: {P_sat/1000:.2f} kPa")
print(f"  Liquid density: {rho_liquid:.2f} kg/m³")
print(f"  Vapor density: {rho_vapor:.2f} kg/m³")
print(f"  Latent heat: {h_fg/1000:.2f} kJ/kg")
```

### Example 3: Pressure-Enthalpy (P-H) State Point
```python
from CoolProp.CoolProp import PropsSI

# Find temperature at known pressure and enthalpy
P = 500000  # 5 bar = 500 kPa
h = 250000  # 250 kJ/kg

T = PropsSI('T', 'P', P, 'H', h, 'R134a')
Q = PropsSI('Q', 'P', P, 'H', h, 'R134a')

print(f"R134a at {P/1000:.0f} kPa and {h/1000:.0f} kJ/kg:")
print(f"  Temperature: {T-273.15:.2f}°C")
print(f"  Quality: {Q:.4f} (0=liquid, 1=vapor)")
```

### Example 4: Critical and Triple Point Data
```python
from CoolProp.CoolProp import Props1SI

fluids = ['Water', 'CO2', 'Nitrogen', 'R134a']

for fluid in fluids:
    T_crit = Props1SI('Tcrit', fluid)
    P_crit = Props1SI('Pcrit', fluid)
    T_triple = Props1SI('Ttriple', fluid)
    M = Props1SI('M', fluid)

    print(f"\n{fluid}:")
    print(f"  Critical point: {T_crit-273.15:.2f}°C, {P_crit/1e6:.2f} MPa")
    print(f"  Triple point: {T_triple-273.15:.2f}°C")
    print(f"  Molar mass: {M*1000:.2f} g/mol")
```

### Example 5: Viscosity Temperature Dependence
```python
from CoolProp.CoolProp import PropsSI
import numpy as np

# Calculate water viscosity from 0°C to 100°C at atmospheric pressure
P = 101325  # Pa
temperatures = np.linspace(273.15, 373.15, 11)  # 0 to 100°C

print("Water viscosity vs temperature:")
print("T(°C)    μ(mPa·s)")
for T in temperatures:
    mu = PropsSI('V', 'T', T, 'P', P, 'Water') * 1000  # Convert to mPa·s
    print(f"{T-273.15:5.0f}    {mu:.4f}")
```

### Example 6: Two-Phase Properties
```python
from CoolProp.CoolProp import PropsSI

# R134a at 10 bar with 50% quality
P = 1000000  # 10 bar = 1 MPa
Q = 0.5  # 50% vapor

T = PropsSI('T', 'P', P, 'Q', Q, 'R134a')
h = PropsSI('H', 'P', P, 'Q', Q, 'R134a')
s = PropsSI('S', 'P', P, 'Q', Q, 'R134a')
rho = PropsSI('D', 'P', P, 'Q', Q, 'R134a')

print(f"R134a two-phase at {P/1e6:.1f} MPa, quality = {Q}:")
print(f"  Temperature: {T-273.15:.2f}°C")
print(f"  Enthalpy: {h/1000:.2f} kJ/kg")
print(f"  Entropy: {s/1000:.4f} kJ/kg·K")
print(f"  Density: {rho:.2f} kg/m³")
```

## Temperature and Pressure Effects

### Temperature Ranges
- Each fluid has valid temperature ranges between triple point and maximum temperature
- Typical range: Triple point temperature < T < 2000 K (varies by fluid)
- Check limits using `Props1SI('Tmin', fluid)` and `Props1SI('Tmax', fluid)`

### Pressure Ranges
- Each fluid has valid pressure ranges
- Typical range: Triple point pressure < P < 1000 MPa (varies by fluid)
- Check limits using `Props1SI('pmin', fluid)` and `Props1SI('pmax', fluid)`

### Phase Regions
1. **Subcooled Liquid**: T < T_sat at given P, or P > P_sat at given T
2. **Two-Phase**: T = T_sat and 0 < Q < 1
3. **Superheated Vapor**: T > T_sat at given P, or P < P_sat at given T
4. **Supercritical**: T > T_crit and P > P_crit

### Input Pair Restrictions
Not all input pairs are valid in all regions:
- **(T, P)**: Valid in single-phase regions only (not in two-phase)
- **(P, Q)**: Valid for two-phase and saturation (0 ≤ Q ≤ 1)
- **(T, Q)**: Valid for two-phase and saturation (0 ≤ Q ≤ 1)
- **(P, H)**: Valid in all regions
- **(P, S)**: Valid in all regions
- **(H, S)**: Valid in all regions (useful for isentropic processes)

## Error Handling

### Common Errors and Solutions

**Error: "CoolProp error: [PropsSI] Two saturation inputs are not valid"**
- Problem: Trying to use (T, P) in two-phase region
- Solution: Use (T, Q) or (P, Q) for two-phase states

**Error: "CoolProp error: Value is outside range"**
- Problem: Temperature or pressure outside valid range
- Solution: Check fluid limits with Props1SI('Tmin', fluid), etc.

**Error: "CoolProp error: Fluid not found"**
- Problem: Incorrect fluid name or spelling
- Solution: Use exact fluid names (case-sensitive), check documentation

**Error: "Unable to match the inputs"**
- Problem: Invalid input combination or iteration failed
- Solution: Check that input values are physically reasonable

### Safe Query Pattern
```python
from CoolProp.CoolProp import PropsSI

def safe_props(output, input1, value1, input2, value2, fluid):
    """Query CoolProp with error handling"""
    try:
        result = PropsSI(output, input1, value1, input2, value2, fluid)
        return result
    except ValueError as e:
        print(f"Error querying {fluid}: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

# Usage
density = safe_props('D', 'T', 300, 'P', 101325, 'Water')
if density is not None:
    print(f"Density: {density:.2f} kg/m³")
```

## Engineering Applications

### HVAC System Design
- Refrigerant cycle analysis (evaporator, condenser, compressor, expansion)
- Heat exchanger calculations
- Coefficient of performance (COP) calculations

### Power Cycles
- Rankine cycle (steam turbines)
- Brayton cycle (gas turbines)
- Organic Rankine Cycle (ORC)

### Process Engineering
- Heat transfer calculations
- Pipe flow and pressure drop
- Phase equilibrium
- Storage tank design

### Cryogenic Systems
- Liquefaction processes
- Storage and handling
- Heat leak calculations

## References

### Official Documentation
- **CoolProp Website**: http://www.coolprop.org/
- **Python Documentation**: http://www.coolprop.org/coolprop/HighLevelAPI.html
- **Fluid Properties**: http://www.coolprop.org/fluid_properties/PurePseudoPure.html
- **Validation Data**: http://www.coolprop.org/validation/index.html

### Key Papers
- Bell, I. H., Wronski, J., Quoilin, S., & Lemort, V. (2014). "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp." *Industrial & Engineering Chemistry Research*, 53(6), 2498-2508.

### Equation of State References
- Lemmon, E. W., Huber, M. L., & McLinden, M. O. (2013). NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP), Version 9.1. National Institute of Standards and Technology.

### Source Code
- **GitHub**: https://github.com/CoolProp/CoolProp
- **Issue Tracker**: https://github.com/CoolProp/CoolProp/issues

## Best Practices

1. **Always use SI units** for inputs and outputs (K, Pa, J/kg, etc.)
2. **Check phase regions** before choosing input pairs
3. **Handle errors** gracefully with try-except blocks
4. **Validate ranges** before querying (especially temperature and pressure)
5. **Use quality (Q)** for two-phase calculations
6. **Cache critical properties** if querying multiple times
7. **Use Props1SI** for fluid-only properties (critical points, molar mass)
8. **Prefer (P, H) or (P, S)** for robust queries across all phase regions

## Quick Reference Table

| Task | Function | Example |
|------|----------|---------|
| Two-input property | `PropsSI(output, in1, val1, in2, val2, fluid)` | `PropsSI('D', 'T', 300, 'P', 101325, 'Water')` |
| Single-input property | `Props1SI(param, fluid)` | `Props1SI('Tcrit', 'Water')` |
| Saturation liquid | Use Q=0 | `PropsSI('H', 'T', 300, 'Q', 0, 'R134a')` |
| Saturation vapor | Use Q=1 | `PropsSI('H', 'T', 300, 'Q', 1, 'R134a')` |
| Two-phase | Use 0<Q<1 | `PropsSI('D', 'P', 500000, 'Q', 0.5, 'R134a')` |
| List all fluids | `CoolProp.__fluids__` | `import CoolProp; print(CoolProp.__fluids__)` |

---

*This skill provides access to one of the most comprehensive open-source thermophysical property databases available, suitable for research, engineering design, and educational applications.*
