---
name: material-properties-db
description: "Query fluid viscosities, densities, and material properties vs temperature"
category: databases
domain: materials
complexity: basic
dependencies: []
---

# Material Properties Database Skill

Query temperature-dependent fluid and material properties essential for pump design, heat transfer, and fluid mechanics calculations. This skill provides verified correlations and empirical data for common engineering fluids.

## Overview

Material property databases provide critical data for engineering calculations:

- **Fluid Properties**: Viscosity, density, surface tension, vapor pressure
- **Temperature Dependence**: Polynomial fits, Sutherland's law, Andrade equation
- **Phase Data**: Saturation properties, freezing/boiling points
- **Transport Properties**: Thermal conductivity, specific heat
- **Dimensionless Numbers**: Reynolds, Prandtl, kinematic viscosity

This skill focuses on practical correlations for fluids commonly encountered in pumping applications, chemical processing, and HVAC systems.

## Common Fluids for Pumps

### Water (H₂O)
The most common pumping fluid with well-established properties:
- **Temperature Range**: 0°C to 100°C (273.15 K to 373.15 K)
- **Density**: ~1000 kg/m³ (decreases slightly with temperature)
- **Viscosity**: Highly temperature-dependent (1.79 mPa·s at 0°C to 0.28 mPa·s at 100°C)
- **Applications**: HVAC, cooling systems, water supply, municipal systems
- **Standards**: IAPWS-95 formulation (International Association for Properties of Water and Steam)

### Hydraulic Oils
Mineral-based and synthetic oils used in hydraulic systems:
- **ISO VG Grades**: VG 32, VG 46, VG 68, VG 100 (viscosity at 40°C)
- **Temperature Range**: -20°C to 100°C typical
- **Density**: 850-900 kg/m³ (relatively constant)
- **Viscosity**: Strong temperature dependence (follows Walther equation)
- **Applications**: Hydraulic pumps, power transmission, control systems
- **Viscosity Index (VI)**: Measure of viscosity-temperature relationship (higher = less change)

### Lubricating Oils
Engine oils and industrial lubricants:
- **SAE Grades**: SAE 10W, 20W, 30, 40, 50
- **Multigrade**: SAE 10W-30, 15W-40, 20W-50
- **Temperature Range**: -40°C to 150°C
- **Density**: 870-920 kg/m³
- **Viscosity**: Engineered for specific temperature ranges
- **Applications**: Bearings, gearboxes, engines, turbines

### Refrigerants
HFC and natural refrigerants for cooling cycles:
- **Common**: R134a, R410A, R32, R717 (ammonia), R744 (CO₂)
- **Temperature Range**: -50°C to 70°C typical
- **Two-Phase Properties**: Critical for evaporators and condensers
- **Pressure Dependent**: Properties vary significantly with pressure
- **Applications**: Chillers, air conditioning, heat pumps, industrial refrigeration
- **Note**: Use CoolProp database for accurate refrigerant properties

### Chemicals and Process Fluids
Common industrial chemicals:
- **Ethylene Glycol**: Antifreeze, heat transfer fluid (-40°C to 100°C)
- **Propylene Glycol**: Food-grade antifreeze, pharmaceuticals
- **Acids/Bases**: Sulfuric acid, caustic soda (corrosive, density ~1.2-1.8 kg/L)
- **Solvents**: Acetone, toluene, methanol, ethanol
- **Hydrocarbons**: Gasoline, diesel, kerosene, crude oil
- **Brines**: Sodium chloride, calcium chloride solutions

### Gases (Compressed)
For gas handling and pipeline calculations:
- **Air**: Standard reference fluid (ideal gas at low pressure)
- **Natural Gas**: Primarily methane, compressible flow
- **Nitrogen**: Inert atmosphere, purging
- **Oxygen**: Medical, combustion applications
- **Note**: Compressibility effects significant at high pressure

## Temperature-Dependent Correlations

### Viscosity Models

#### Andrade Equation (Liquids)
Simple exponential model for liquid viscosity:

```
μ(T) = A · exp(B/T)
```

Where:
- μ = dynamic viscosity (Pa·s or mPa·s)
- T = absolute temperature (K)
- A, B = fluid-specific constants

**Good for**: Quick estimates, limited temperature ranges
**Accuracy**: ±5-10% for moderate temperature ranges

#### Vogel-Fulcher-Tammann Equation (Better for Oils)
More accurate for oils and high-viscosity fluids:

```
μ(T) = A · exp(B/(T - C))
```

Where:
- C = typically 95-140 K for oils
- Better fit over wide temperature ranges

#### Walther Equation (Petroleum Products)
ASTM D341 standard for petroleum oils:

```
log₁₀(log₁₀(ν + 0.7)) = A - B·log₁₀(T)
```

Where:
- ν = kinematic viscosity (cSt = mm²/s)
- T = absolute temperature (K)
- A, B = constants from two-point calibration

**Used for**: ISO VG oils, SAE grades, ASTM viscosity indices
**Accuracy**: Excellent for petroleum products

#### Sutherland's Law (Gases)
For gas viscosity temperature dependence:

```
μ(T) = μ₀ · (T/T₀)^(3/2) · (T₀ + S)/(T + S)
```

Where:
- μ₀ = reference viscosity at T₀
- T₀ = reference temperature (often 273.15 K)
- S = Sutherland constant (K)
  - Air: S = 110.4 K
  - Nitrogen: S = 111 K
  - Oxygen: S = 127 K

**Good for**: Ideal gases at moderate pressures
**Range**: Valid from ~100 K to 2000 K

### Density Models

#### Linear Approximation (Liquids)
For incompressible liquids over moderate temperature ranges:

```
ρ(T) = ρ₀ · [1 - β(T - T₀)]
```

Where:
- ρ₀ = density at reference temperature T₀ (kg/m³)
- β = volumetric thermal expansion coefficient (1/K)
  - Water: β ≈ 0.0002 K⁻¹ near 20°C
  - Oils: β ≈ 0.0007 K⁻¹

#### Polynomial Fit (Water)
IAPWS-IF97 simplified for atmospheric pressure:

```
ρ(T) = a₀ + a₁·T + a₂·T² + a₃·T³
```

**For water (0-100°C at 1 atm)**:
- High accuracy (±0.01%)
- Coefficients from NIST or steam tables

#### Ideal Gas Law (Gases)
For gases at low to moderate pressure:

```
ρ = P·M / (R·T)
```

Where:
- P = absolute pressure (Pa)
- M = molar mass (kg/mol)
- R = universal gas constant = 8.314 J/(mol·K)
- T = absolute temperature (K)

### Vapor Pressure Models

#### Antoine Equation
Most common correlation for vapor pressure:

```
log₁₀(P_vap) = A - B/(T + C)
```

Where:
- P_vap = vapor pressure (mmHg, kPa, or bar depending on constants)
- T = temperature (°C or K, depending on constants)
- A, B, C = fluid-specific constants

**Common fluids (T in °C, P in mmHg)**:
- **Water**: A=8.07131, B=1730.63, C=233.426 (1-100°C)
- **Ethanol**: A=8.04494, B=1554.3, C=222.65 (20-93°C)
- **Methanol**: A=7.89750, B=1474.08, C=229.13

**Applications**:
- NPSH calculations (Net Positive Suction Head)
- Cavitation prediction
- Flash point estimation
- Boiling point at altitude

#### Clausius-Clapeyron Equation
Thermodynamic basis for vapor pressure:

```
ln(P₂/P₁) = -ΔH_vap/R · (1/T₂ - 1/T₁)
```

Where:
- ΔH_vap = heat of vaporization (J/mol)
- R = gas constant = 8.314 J/(mol·K)

**Good for**: Extrapolation from known point, theoretical calculations

### Kinematic Viscosity
Relationship between dynamic and kinematic viscosity:

```
ν = μ / ρ
```

Where:
- ν = kinematic viscosity (m²/s or cSt)
- μ = dynamic viscosity (Pa·s)
- ρ = density (kg/m³)
- Conversion: 1 cSt = 1 mm²/s = 10⁻⁶ m²/s

**Important for**:
- Reynolds number calculations
- ISO VG oil ratings (viscosity at 40°C in cSt)
- Viscometer measurements

## Data Sources and Standards

### Primary Sources

#### NIST (National Institute of Standards and Technology)
- **NIST Chemistry WebBook**: https://webbook.nist.gov/chemistry/
- **Properties**: Thermophysical data for thousands of compounds
- **Accuracy**: Research-grade, high reliability
- **Coverage**: Density, viscosity, vapor pressure, thermal properties

#### IAPWS (International Association for Properties of Water and Steam)
- **IAPWS-95**: Water and steam properties formulation
- **IAPWS-IF97**: Industrial formulation (simpler, faster)
- **Coverage**: 0-1000°C, 0-1000 MPa
- **Accuracy**: Best available for water/steam

#### Perry's Chemical Engineers' Handbook
- **Publisher**: McGraw-Hill
- **Content**: Comprehensive physical property data
- **Correlations**: Empirical equations for thousands of fluids
- **Industry Standard**: Widely used in chemical engineering

#### ASHRAE Handbooks
- **Coverage**: HVAC fluids, refrigerants, psychrometrics
- **Updates**: Annual updates for refrigerants
- **Applications**: Building systems, refrigeration

### Standards Organizations

#### ASTM International
- **ASTM D341**: Viscosity-temperature charts for petroleum products
- **ASTM D445**: Kinematic viscosity measurement
- **ASTM D2270**: Viscosity index calculation
- **ASTM D6751**: Biodiesel specifications

#### ISO (International Organization for Standardization)
- **ISO 3448**: Industrial liquid lubricant viscosity grades (VG system)
- **ISO 12185**: Crude petroleum and petroleum products density
- **ISO 2909**: Petroleum measurement tables

#### API (American Petroleum Institute)
- **API gravity**: Oil density scale (°API)
- **Technical Data Book**: Petroleum refining properties

### Software and Databases

#### CoolProp
- Open-source thermophysical property library
- 100+ pure and pseudo-pure fluids
- High-accuracy equations of state
- See `coolprop-db` skill for details

#### REFPROP (NIST)
- Reference fluid thermodynamic properties
- Gold standard for accuracy
- Commercial license required
- Based on peer-reviewed equations of state

#### Engineering Equation Solver (EES)
- Built-in property database
- Automatic unit conversion
- Educational and professional versions

## Practical Usage Guidelines

### Property Selection for Pump Design

1. **Viscosity**: Critical for Reynolds number, friction losses
   - Use kinematic viscosity (ν) for Re calculations
   - Dynamic viscosity (μ) for wall shear stress

2. **Density**: Affects head-pressure conversion, power requirements
   - Use average density for approximate calculations
   - Temperature-corrected for accurate NPSH

3. **Vapor Pressure**: Essential for NPSH available calculations
   - Must be evaluated at pumping temperature
   - Critical for hot fluids or low suction pressure

4. **Specific Gravity**: Ratio to water density (dimensionless)
   - SG = ρ_fluid / ρ_water @ 4°C
   - Simplifies pump curve scaling

### Temperature Considerations

- **Design Point**: Select properties at maximum/minimum operating temperature
- **Startup**: Consider cold start conditions (high viscosity)
- **Seasonal Variation**: Account for ambient temperature effects
- **Heat Generation**: Pump inefficiency adds heat to fluid

### Uncertainty and Safety Factors

- **Property Uncertainty**: ±5% typical for correlations
- **Viscosity Range**: Design for ±20% variation if uncertain
- **NPSH Margin**: Add 0.5-1.0 m safety margin above required
- **Verification**: Always verify critical properties against multiple sources

## Query Methods

### Manual Calculation
Use empirical equations with fluid-specific constants:
```python
import math

def water_viscosity(T_celsius):
    """Vogel equation for water viscosity"""
    A = 0.02414  # mPa·s
    B = 247.8    # K
    C = 140      # K
    T_kelvin = T_celsius + 273.15
    mu = A * 10**(B / (T_kelvin - C))
    return mu  # mPa·s
```

### Tabular Interpolation
Linear or polynomial interpolation from standard tables:
```python
import numpy as np

# Example: Water density table
T_data = np.array([0, 20, 40, 60, 80, 100])  # °C
rho_data = np.array([999.8, 998.2, 992.2, 983.2, 971.8, 958.4])  # kg/m³

def interpolate_density(T):
    return np.interp(T, T_data, rho_data)
```

### Database Lookup
Use libraries like CoolProp for high-accuracy data:
```python
from CoolProp.CoolProp import PropsSI

# Water viscosity at 25°C, 1 atm
mu = PropsSI('V', 'T', 298.15, 'P', 101325, 'Water')
```

## Engineering Applications

### Reynolds Number Calculation
```python
Re = ρ · v · D / μ = v · D / ν
```
- Determines flow regime (laminar vs turbulent)
- Critical for friction factor selection
- Typical pump range: Re = 10⁵ to 10⁷

### NPSH Available
```python
NPSH_a = (P_atm - P_vap) / (ρ·g) + h_static - h_friction
```
- Requires vapor pressure at pumping temperature
- Prevents cavitation
- Must exceed NPSH_required by margin

### Pressure-Head Conversion
```python
H = ΔP / (ρ·g)
```
- H = head (m)
- ΔP = pressure rise (Pa)
- ρ = fluid density (kg/m³)
- g = 9.81 m/s²

### Power Calculation
```python
P_hydraulic = ρ · g · Q · H
P_shaft = P_hydraulic / η_pump
```
- Density affects power requirements directly
- Higher specific gravity = higher power

## Best Practices

1. **Always use absolute temperature** (Kelvin) for correlations
2. **Verify units** - many correlations use mixed units (°C, mmHg, cSt)
3. **Check validity range** - don't extrapolate beyond calibrated range
4. **Use multiple sources** for critical applications
5. **Document assumptions** - property source, temperature, pressure
6. **Consider impurities** - real fluids differ from pure substance data
7. **Account for aging** - oil degradation changes viscosity over time
8. **Validate with measurements** when possible (viscometer, hydrometer)

## Quick Reference Data

### Water at Atmospheric Pressure

| T (°C) | ρ (kg/m³) | μ (mPa·s) | ν (mm²/s) | P_vap (kPa) |
|--------|-----------|-----------|-----------|-------------|
| 0      | 999.8     | 1.787     | 1.787     | 0.611       |
| 10     | 999.7     | 1.307     | 1.307     | 1.228       |
| 20     | 998.2     | 1.002     | 1.004     | 2.339       |
| 25     | 997.0     | 0.890     | 0.893     | 3.169       |
| 30     | 995.7     | 0.798     | 0.801     | 4.246       |
| 40     | 992.2     | 0.653     | 0.658     | 7.384       |
| 50     | 988.0     | 0.547     | 0.554     | 12.35       |
| 60     | 983.2     | 0.467     | 0.475     | 19.94       |
| 70     | 977.8     | 0.404     | 0.413     | 31.19       |
| 80     | 971.8     | 0.355     | 0.365     | 47.39       |
| 90     | 965.3     | 0.315     | 0.326     | 70.14       |
| 100    | 958.4     | 0.282     | 0.294     | 101.3       |

### Common Oil Viscosities at 40°C

| ISO VG | ν @ 40°C (cSt) | ρ (kg/m³) | μ @ 40°C (mPa·s) |
|--------|----------------|-----------|------------------|
| VG 32  | 32             | 865       | 27.7             |
| VG 46  | 46             | 870       | 40.0             |
| VG 68  | 68             | 875       | 59.5             |
| VG 100 | 100            | 880       | 88.0             |
| VG 150 | 150            | 885       | 132.8            |

### Sutherland Constants for Common Gases

| Gas      | μ₀ @ 273K (μPa·s) | S (K)  | Valid Range |
|----------|-------------------|--------|-------------|
| Air      | 17.16             | 110.4  | 100-1900 K  |
| N₂       | 16.66             | 111    | 100-1900 K  |
| O₂       | 19.20             | 127    | 100-1900 K  |
| CO₂      | 13.73             | 240    | 200-1900 K  |
| H₂       | 8.41              | 72     | 100-1900 K  |

---

*This skill provides practical correlations and data sources for material properties essential to pump design, fluid mechanics, and thermal engineering applications.*
