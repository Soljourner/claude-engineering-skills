---
name: fluid-property-calculator
description: "Quick fluid property calculations using empirical formulas without database queries"
category: helpers
domain: fluids
complexity: basic
dependencies: []
---

# Fluid Property Calculator

Quick fluid property calculations using empirical correlations and analytical formulas. Provides instant property estimates without requiring external databases or data files.

## Overview

This helper provides fast, practical calculations for common fluids using well-established empirical correlations. All formulas include validity ranges and are verified against reference data.

## Available Calculations

### Water Properties (0-100°C)
Calculate temperature-dependent properties of liquid water:
- **Density** (kg/m³) - Polynomial correlation
- **Dynamic viscosity** (Pa·s) - Vogel equation
- **Kinematic viscosity** (m²/s) - Derived from dynamic viscosity
- **Thermal conductivity** (W/m·K) - Polynomial correlation
- **Specific heat capacity** (J/kg·K) - Polynomial correlation
- **Vapor pressure** (Pa) - Antoine equation
- **Prandtl number** - Dimensionless heat transfer parameter

**Validity Range**: 0-100°C at atmospheric pressure
**Typical Accuracy**: ±1-2% for most properties

### Air Properties (Standard Atmosphere)
Calculate temperature-dependent properties of air at atmospheric pressure:
- **Density** (kg/m³) - Ideal gas law
- **Dynamic viscosity** (Pa·s) - Sutherland's formula
- **Kinematic viscosity** (m²/s) - Derived from dynamic viscosity
- **Thermal conductivity** (W/m·K) - Polynomial correlation
- **Specific heat capacity** (J/kg·K) - Temperature-dependent correlation
- **Prandtl number** - Dimensionless heat transfer parameter

**Validity Range**: -50 to 200°C at 101.325 kPa
**Typical Accuracy**: ±1-3% for most properties

### Viscosity Correlations

#### Sutherland's Formula (Gases)
Temperature-dependent viscosity for gases:
```
μ = μ₀ × (T/T₀)^(3/2) × (T₀ + S)/(T + S)
```

Available for:
- Air (S = 110.4 K)
- Nitrogen (S = 111 K)
- Oxygen (S = 127 K)
- Carbon dioxide (S = 240 K)

#### Andrade Equation (Liquids)
Temperature-dependent viscosity for liquids:
```
μ = A × exp(B/T)
```

Provides empirical correlation for various liquids with custom parameters.

### Vapor Pressure (Antoine Equation)
Calculate saturation vapor pressure:
```
log₁₀(P) = A - B/(C + T)
```

Available for:
- Water
- Ethanol
- Methanol
- Acetone
- Benzene
- Toluene

**Units**: Temperature in °C, Pressure in mmHg or Pa (depending on constants)

### Dimensionless Numbers

#### Reynolds Number
Calculate flow regime indicator:
```
Re = ρ × V × L / μ = V × L / ν
```

Where:
- ρ = density (kg/m³)
- V = velocity (m/s)
- L = characteristic length (m)
- μ = dynamic viscosity (Pa·s)
- ν = kinematic viscosity (m²/s)

**Interpretation**:
- Re < 2300: Laminar flow (pipe)
- 2300 < Re < 4000: Transition
- Re > 4000: Turbulent flow (pipe)

#### Friction Factor Calculator

**Laminar Flow (Re < 2300)**:
```
f = 64 / Re
```

**Turbulent Flow - Smooth Pipes (Blasius)**:
Valid for Re < 100,000:
```
f = 0.316 / Re^0.25
```

**Turbulent Flow - Rough Pipes (Colebrook-White)**:
Iterative solution for:
```
1/√f = -2 log₁₀(ε/3.7D + 2.51/(Re√f))
```

Where:
- ε = absolute roughness (m)
- D = pipe diameter (m)

**Swamee-Jain Approximation** (non-iterative):
```
f = 0.25 / [log₁₀(ε/3.7D + 5.74/Re^0.9)]²
```

### Ideal Gas Properties

Calculate properties using ideal gas law and kinetic theory:
- **Density**: ρ = P/(R×T)
- **Specific heat ratio**: γ (for common gases)
- **Speed of sound**: a = √(γ×R×T)
- **Molar mass**: M (for common gases)

## Usage Examples

### Example 1: Water Properties at 20°C
```python
from calc import water_properties

props = water_properties(20)
print(f"Density: {props['density']:.2f} kg/m³")
print(f"Viscosity: {props['dynamic_viscosity']:.6f} Pa·s")
print(f"Thermal conductivity: {props['thermal_conductivity']:.4f} W/m·K")
```

### Example 2: Reynolds Number for Pipe Flow
```python
from calc import reynolds_number, water_properties

T = 25  # °C
V = 1.5  # m/s
D = 0.05  # m

props = water_properties(T)
Re = reynolds_number(V, D, props['kinematic_viscosity'])
print(f"Reynolds number: {Re:.0f}")
```

### Example 3: Friction Factor Calculation
```python
from calc import friction_factor

Re = 50000
roughness = 0.045e-3  # 0.045 mm for commercial steel
diameter = 0.1  # m

f = friction_factor(Re, roughness, diameter)
print(f"Friction factor: {f:.5f}")
```

### Example 4: Vapor Pressure of Water
```python
from calc import antoine_vapor_pressure

T = 80  # °C
P_vap = antoine_vapor_pressure('water', T)
print(f"Vapor pressure at {T}°C: {P_vap/1000:.2f} kPa")
```

### Example 5: Air Viscosity using Sutherland's Formula
```python
from calc import sutherland_viscosity

T = 100  # °C
mu = sutherland_viscosity('air', T + 273.15)  # Convert to Kelvin
print(f"Air viscosity at {T}°C: {mu:.6f} Pa·s")
```

## Quick Reference

### Common Water Properties
| T (°C) | ρ (kg/m³) | μ (mPa·s) | ν (mm²/s) | k (W/m·K) | Pr |
|--------|-----------|-----------|-----------|-----------|-----|
| 0      | 999.8     | 1.787     | 1.787     | 0.561     | 13.5 |
| 20     | 998.2     | 1.002     | 1.004     | 0.598     | 7.0 |
| 40     | 992.2     | 0.653     | 0.658     | 0.631     | 4.3 |
| 60     | 983.2     | 0.467     | 0.475     | 0.654     | 3.0 |
| 80     | 971.8     | 0.355     | 0.365     | 0.670     | 2.2 |
| 100    | 958.4     | 0.282     | 0.294     | 0.680     | 1.8 |

### Common Air Properties (at 101.325 kPa)
| T (°C) | ρ (kg/m³) | μ (μPa·s) | ν (mm²/s) | k (W/m·K) | Pr |
|--------|-----------|-----------|-----------|-----------|-----|
| 0      | 1.293     | 17.16     | 13.27     | 0.0243    | 0.71 |
| 20     | 1.205     | 18.24     | 15.14     | 0.0257    | 0.71 |
| 50     | 1.093     | 19.57     | 17.90     | 0.0279    | 0.71 |
| 100    | 0.946     | 21.67     | 22.90     | 0.0314    | 0.71 |

## Limitations

1. **Temperature Ranges**: Correlations are only valid within specified ranges
2. **Pressure Effects**: Most correlations assume atmospheric pressure
3. **Pure Substances**: Mixtures require different approaches
4. **Accuracy**: Empirical formulas provide estimates (±1-5% typical)
5. **Phase Changes**: Properties near phase transitions may be less accurate

## When to Use This Helper

**Good for**:
- Quick engineering calculations
- Preliminary design work
- Educational purposes
- When databases are unavailable
- Rapid prototyping

**Not suitable for**:
- High-precision scientific work
- Properties outside validity ranges
- Non-standard conditions (high pressure, etc.)
- Complex mixtures
- When accuracy better than ±1% is required

## Best Practices

1. **Check validity ranges** before using any correlation
2. **Verify units** - most functions use SI units (K for temperature in Sutherland, °C elsewhere)
3. **Compare results** with reference data when possible
4. **Use appropriate significant figures** based on correlation accuracy
5. **Document assumptions** in your calculations

## Additional Resources

See `reference.md` for:
- Complete correlation equations
- Literature sources
- Validation data
- Accuracy comparisons
- Alternative formulations
