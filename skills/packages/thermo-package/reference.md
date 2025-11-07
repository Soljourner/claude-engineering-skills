# Thermo Package Reference

## Quick Reference Guide

### Installation

```bash
pip install thermo chemicals
```

### Basic Usage

```python
from thermo.chemical import Chemical
from thermo.mixture import Mixture

# Single component
water = Chemical('water', T=298.15, P=101325)  # 25°C, 1 atm
print(f"Density: {water.rho} kg/m³")
print(f"Viscosity: {water.mu} Pa·s")
print(f"Vapor pressure: {water.Psat} Pa")

# Mixture
mix = Mixture(['ethanol', 'water'], ws=[0.3, 0.7], T=298.15, P=101325)
print(f"Mixture density: {mix.rho} kg/m³")
```

---

## Chemical Database

### Common Pump Fluids

The thermo library provides access to 20,000+ chemicals via the Chemicals database. Below are common fluids encountered in pump applications.

#### Water and Aqueous Solutions

| Chemical Name | CAS Number | Formula | Notes |
|--------------|------------|---------|-------|
| water | 7732-18-5 | H₂O | Universal solvent, most common pump fluid |
| heavy water | 7789-20-0 | D₂O | Neutron moderator in nuclear applications |
| seawater | N/A | Mixed | Use mixture with appropriate salts |

#### Alcohols and Glycols

| Chemical Name | CAS Number | Formula | Applications |
|--------------|------------|---------|--------------|
| methanol | 67-56-1 | CH₃OH | Fuel, solvent, antifreeze |
| ethanol | 64-17-5 | C₂H₅OH | Fuel, beverage, solvent |
| isopropanol | 67-63-0 | C₃H₇OH | Cleaning, disinfectant |
| glycerol | 56-81-5 | C₃H₈O₃ | Pharmaceutical, food processing |
| ethylene glycol | 107-21-1 | C₂H₆O₂ | Antifreeze, coolant |
| propylene glycol | 57-55-6 | C₃H₈O₂ | Food-grade antifreeze |

#### Hydrocarbons

| Chemical Name | CAS Number | Formula | Applications |
|--------------|------------|---------|--------------|
| methane | 74-82-8 | CH₄ | Natural gas, LNG |
| ethane | 74-84-0 | C₂H₆ | Natural gas component |
| propane | 74-98-6 | C₃H₈ | LPG, fuel |
| butane | 106-97-8 | C₄H₁₀ | LPG, lighter fuel |
| pentane | 109-66-0 | C₅H₁₂ | Solvent, blowing agent |
| hexane | 110-54-3 | C₆H₁₄ | Solvent, extraction |
| heptane | 142-82-5 | C₇H₁₆ | Reference fuel, solvent |
| octane | 111-65-9 | C₈H₁₈ | Gasoline component |
| benzene | 71-43-2 | C₆H₆ | Chemical intermediate |
| toluene | 108-88-3 | C₇H₈ | Solvent, fuel additive |
| xylene | 1330-20-7 | C₈H₁₀ | Solvent, chemical intermediate |

#### Petroleum Products

| Chemical Name | Mixture Type | Applications |
|--------------|-------------|--------------|
| gasoline | Complex mixture | Automotive fuel |
| diesel | Complex mixture | Diesel engines, heating |
| jet fuel | Complex mixture | Aviation fuel |
| kerosene | Complex mixture | Heating, aviation |
| crude oil | Complex mixture | Petroleum feedstock |

**Note:** For petroleum products, use representative single components or create custom mixtures.

#### Refrigerants

| Chemical Name | CAS Number | Type | Applications |
|--------------|------------|------|--------------|
| ammonia | 7664-41-7 | NH₃ | Industrial refrigeration |
| R-134a | 811-97-2 | HFC | Automotive A/C, refrigeration |
| R-410A | Mixed | HFC blend | Air conditioning |
| R-22 | 75-45-6 | HCFC | Legacy refrigerant (phasing out) |
| propane | 74-98-6 | HC | Natural refrigerant |
| carbon dioxide | 124-38-9 | CO₂ | Transcritical refrigeration |

#### Acids and Bases

| Chemical Name | CAS Number | Formula | Concentration Notes |
|--------------|------------|---------|-------------------|
| sulfuric acid | 7664-93-9 | H₂SO₄ | Specify concentration (e.g., 98%) |
| hydrochloric acid | 7647-01-0 | HCl | Specify concentration (e.g., 37%) |
| nitric acid | 7697-37-2 | HNO₃ | Specify concentration |
| acetic acid | 64-19-7 | CH₃COOH | Vinegar (5%), glacial (99.5%) |
| sodium hydroxide | 1310-73-2 | NaOH | Caustic soda (50% solution) |
| ammonia solution | 1336-21-6 | NH₃·H₂O | Specify concentration (28%) |

#### Solvents

| Chemical Name | CAS Number | Formula | Applications |
|--------------|------------|---------|--------------|
| acetone | 67-64-1 | C₃H₆O | Solvent, cleaning |
| MEK | 78-93-3 | C₄H₈O | Solvent, adhesives |
| ethyl acetate | 141-78-6 | C₄H₈O₂ | Solvent, coatings |
| dichloromethane | 75-09-2 | CH₂Cl₂ | Paint stripper, solvent |
| chloroform | 67-66-3 | CHCl₃ | Solvent, reagent |
| DMF | 68-12-2 | C₃H₇NO | Polar aprotic solvent |
| DMSO | 67-68-5 | C₂H₆OS | Polar aprotic solvent |

#### Thermal Fluids

| Chemical Name | CAS Number | Type | Max Temp (°C) |
|--------------|------------|------|---------------|
| Dowtherm A | 8004-13-5 | Diphenyl/diphenyl oxide | 400 |
| Therminol VP-1 | Mixed | Diphenyl/diphenyl oxide | 400 |
| ethylene glycol | 107-21-1 | Glycol | 160 |
| propylene glycol | 57-55-6 | Glycol | 160 |

**Note:** Use generic diphenyl oxide or create mixture for proprietary thermal oils.

#### Gases (at standard conditions)

| Chemical Name | CAS Number | Formula | Pump Type Required |
|--------------|------------|---------|-------------------|
| nitrogen | 7727-37-9 | N₂ | Gas booster, compressor |
| oxygen | 7782-44-7 | O₂ | Cryogenic pump |
| hydrogen | 1333-74-0 | H₂ | Diaphragm pump |
| carbon dioxide | 124-38-9 | CO₂ | Cryogenic or gas booster |
| argon | 7440-37-1 | Ar | Cryogenic pump |
| helium | 7440-59-7 | He | Cryogenic pump |

---

## Equation of State Options

### Available EOS Models

The thermo library implements multiple equations of state (EOS) for different applications:

#### 1. Ideal Gas Law
```python
from thermo.eos import IG
eos = IG(T=298.15, P=101325)
```
- **Use when:** Low pressure (<5 bar), high temperature
- **Accuracy:** Poor for liquids, good for gases at low pressure
- **Applications:** Approximate gas calculations

#### 2. Peng-Robinson (PR)
```python
from thermo.eos import PR
eos = PR(Tc=647.14, Pc=22048320, omega=0.344, T=298.15, P=101325)
```
- **Use when:** General purpose, VLE calculations
- **Accuracy:** Good for hydrocarbons and non-polar mixtures
- **Applications:** Oil & gas, refrigeration, most process fluids
- **Limitations:** Less accurate for highly polar systems

#### 3. Soave-Redlich-Kwong (SRK)
```python
from thermo.eos import SRK
eos = SRK(Tc=647.14, Pc=22048320, omega=0.344, T=298.15, P=101325)
```
- **Use when:** Similar to PR, older correlation
- **Accuracy:** Slightly less accurate than PR for vapor pressure
- **Applications:** Legacy systems, when PR not available
- **Advantages:** More stable numerically in some cases

#### 4. Van der Waals
```python
from thermo.eos import VDW
eos = VDW(Tc=647.14, Pc=22048320, T=298.15, P=101325)
```
- **Use when:** Educational purposes, simple estimates
- **Accuracy:** Poor compared to modern EOS
- **Applications:** Teaching, rough estimates

#### 5. IAPWS-95 (Water)
```python
from thermo.chemical import Chemical
water = Chemical('water', T=298.15, P=101325)
# Automatically uses IAPWS-95 for water
```
- **Use when:** Water or steam applications
- **Accuracy:** Highest accuracy for water (±0.01%)
- **Applications:** Steam tables, boiler feed, power generation
- **Range:** 0-1000°C, 0-1000 MPa

#### 6. CoolProp Integration
```python
# Thermo can use CoolProp as backend for enhanced accuracy
# Install: pip install CoolProp
```
- **Use when:** Highest accuracy needed
- **Accuracy:** Reference-quality for supported fluids
- **Applications:** Refrigeration, cryogenics, precise calculations

### EOS Selection Guide for Pump Applications

| Application | Recommended EOS | Reason |
|-------------|----------------|--------|
| Water (any T, P) | IAPWS-95 | Highest accuracy for water |
| Hydrocarbons | Peng-Robinson | Excellent for non-polar fluids |
| Natural gas | Peng-Robinson | Industry standard for gas mixtures |
| Refrigerants | CoolProp or PR | CoolProp more accurate if available |
| Thermal fluids | Peng-Robinson | Good for organic liquids |
| Steam | IAPWS-95 | Steam table standard |
| Air (low P) | Ideal gas | Sufficient accuracy |
| Cryogenic fluids | CoolProp | Better low-temperature accuracy |

### EOS Implementation Example

```python
from thermo import ChemicalConstantsPackage, PRMIX, CEOSGas, CEOSLiquid
from thermo.interaction_parameters import IPDB

# Define mixture
IDs = ['methane', 'ethane', 'propane']
zs = [0.85, 0.10, 0.05]  # Mole fractions
T = 298.15  # K
P = 5e6  # Pa (50 bar)

# Get chemical constants
constants = ChemicalConstantsPackage.constants_from_IDs(IDs)

# Binary interaction parameters
kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')

# EOS parameters
eos_kwargs = {
    'Tcs': constants.Tcs,
    'Pcs': constants.Pcs,
    'omegas': constants.omegas,
    'kijs': kijs
}

# Create EOS objects
gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=constants.HeatCapacityGases,
              T=T, P=P, zs=zs)
liquid = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=constants.HeatCapacityGases,
                    T=T, P=P, zs=zs)

print(f"Gas density: {gas.rho_mass()} kg/m³")
print(f"Liquid density: {liquid.rho_mass()} kg/m³")
```

---

## Property Methods and Correlations

### Temperature-Dependent Viscosity

Multiple correlations available (automatically selected based on data):

- **Yaws correlation:** Wide temperature range
- **DIPPR 101:** Industry standard
- **Viswanath-Natarajan:** High-temperature liquids
- **Vogel equation:** Simple, accurate near ambient

```python
water = Chemical('water', T=298.15)
print(f"Viscosity method: {water.mu_methods}")
# Shows available methods, first is default
```

### Vapor Pressure

Correlations (in order of typical preference):

- **Antoine equation:** Simple, accurate in limited range
- **Wagner equation:** Wide range, high accuracy
- **DIPPR 101:** Industry standard
- **Ambrose-Walton:** Corresponding states

```python
water = Chemical('water', T=298.15)
print(f"Vapor pressure methods: {water.Psat_methods}")
```

### Heat Capacity

- **DIPPR 107:** Gases
- **DIPPR 100:** Liquids
- **Poling correlation:** Estimation method
- **TRC tables:** High accuracy

---

## Property Symbols and Units

### All properties use SI units

| Property | Symbol | Units | Description |
|----------|--------|-------|-------------|
| Temperature | T | K | Absolute temperature |
| Pressure | P | Pa | Absolute pressure |
| Density | rho, ρ | kg/m³ | Mass density |
| Viscosity (dynamic) | mu, μ | Pa·s | Dynamic viscosity |
| Viscosity (kinematic) | nu, ν | m²/s | Kinematic viscosity (ν = μ/ρ) |
| Heat capacity (Cp) | Cp | J/kg·K | Constant pressure heat capacity |
| Heat capacity (Cv) | Cv | J/kg·K | Constant volume heat capacity |
| Enthalpy | H | J/kg | Specific enthalpy |
| Entropy | S | J/kg·K | Specific entropy |
| Vapor pressure | Psat | Pa | Saturation pressure |
| Thermal conductivity | k | W/m·K | Thermal conductivity |
| Surface tension | sigma, σ | N/m | Surface tension |
| Molecular weight | MW | g/mol | Molecular weight |

### Unit Conversions

```python
# Pressure
Pa = 1
bar = 100000  # Pa
psi = 6894.76  # Pa
atm = 101325  # Pa
kPa = 1000  # Pa
MPa = 1e6  # Pa

# Viscosity
Pa_s = 1
cP = 0.001  # Pa·s (centipoise)
cSt_to_m2s = 1e-6  # m²/s (centistokes to m²/s)

# Temperature
def C_to_K(T_C):
    return T_C + 273.15

def K_to_C(T_K):
    return T_K - 273.15

# Density
kg_m3 = 1
g_cm3 = 1000  # kg/m³
lb_ft3 = 16.0185  # kg/m³
```

---

## Validation and Data Sources

### Primary Data Sources

1. **DIPPR 801 Database**
   - AIChE Design Institute for Physical Properties
   - Industry-standard correlations
   - Temperature-dependent properties
   - Wide chemical coverage

2. **NIST Chemistry WebBook**
   - Reference-quality data
   - Experimental values
   - Validated correlations

3. **Yaws' Handbook**
   - Comprehensive property data
   - Temperature-dependent correlations
   - Critical properties

4. **IAPWS-95** (Water)
   - International standard for water/steam
   - Accuracy: ±0.01% for density
   - Range: 0-1000°C, 0-1000 MPa

5. **RefProp** (NIST)
   - Reference fluid properties
   - Used for validation
   - Highest accuracy available

### Validation Procedure

When using thermo for pump calculations, validate against:

1. **Steam tables** (for water/steam)
2. **Perry's Chemical Engineers' Handbook**
3. **Crane TP-410** (for common pump fluids)
4. **Manufacturer data sheets** (for proprietary fluids)
5. **CoolProp** (for cross-validation)

Example validation:
```python
from thermo.chemical import Chemical

# Validate water density at 20°C
water = Chemical('water', T=293.15, P=101325)
rho_calculated = water.rho  # Should be ~998.2 kg/m³
rho_reference = 998.2  # From steam tables

error = abs(rho_calculated - rho_reference) / rho_reference * 100
print(f"Error: {error:.3f}%")  # Should be < 0.1%
```

---

## Performance and Optimization

### Caching Chemical Objects

Create Chemical objects once and reuse:

```python
# BAD: Creates new object every iteration
for T in temperature_range:
    water = Chemical('water', T=T, P=101325)
    rho = water.rho

# GOOD: Reuse object
water = Chemical('water', T=298.15, P=101325)
densities = []
for T in temperature_range:
    water.T = T  # Update temperature
    densities.append(water.rho)
```

### Vectorized Calculations

For large datasets, consider using property methods directly:

```python
import numpy as np
from chemicals.viscosity import mu_IAPWS

# Vectorized viscosity calculation
T_array = np.linspace(273.15, 373.15, 100)
mu_array = [mu_IAPWS(T, 101325) for T in T_array]
```

### Property Calculation Cost

Relative computational cost (1 = fastest):

| Property | Cost | Notes |
|----------|------|-------|
| MW | 1 | Constant |
| Tc, Pc | 1 | Constant |
| rho | 10 | Correlation evaluation |
| mu | 10 | Correlation evaluation |
| Psat | 15 | May require iteration |
| H, S | 20 | Integration of Cp |
| Flash | 100+ | Iterative VLE solve |

---

## Troubleshooting

### Common Issues and Solutions

#### 1. Property Returns None

**Symptom:**
```python
chemical.rho  # Returns None
```

**Causes:**
- Temperature outside valid range
- Property not available for this chemical
- Phase not properly defined

**Solutions:**
```python
# Check valid temperature range
print(f"T min: {chemical.Tmin}, T max: {chemical.Tmax}")

# Check available methods
print(f"Density methods: {chemical.rho_methods}")

# Verify phase
print(f"Phase: {chemical.phase}")
```

#### 2. Flash Calculation Fails

**Symptom:**
```python
flasher.flash(T=T, P=P, zs=zs)  # Raises error or doesn't converge
```

**Solutions:**
- Check if T and P are in reasonable range
- Verify composition sums to 1.0
- Try different initial guess
- Use simpler EOS (PR instead of complex mixing rules)

#### 3. Mixture Properties Incorrect

**Symptom:**
- Properties don't match expected behavior
- Negative values or unrealistic results

**Solutions:**
```python
# Verify mass fractions sum to 1
mix = Mixture(['water', 'ethanol'], ws=[0.6, 0.4], T=298.15, P=101325)
print(f"Sum of ws: {sum([0.6, 0.4])}")  # Should be 1.0

# Check for phase separation
print(f"Phase: {mix.phase}")

# Validate against pure components
water = Chemical('water', T=298.15, P=101325)
ethanol = Chemical('ethanol', T=298.15, P=101325)
```

#### 4. Temperature/Pressure Out of Range

**Symptom:**
```python
Chemical('water', T=1000, P=101325)  # May fail or extrapolate
```

**Solutions:**
- Check correlation validity range
- Use appropriate EOS for conditions
- For water above 374°C (critical point), use IAPWS-95

---

## Links and Resources

### Official Documentation
- **Thermo Documentation:** https://thermo.readthedocs.io/
- **Chemicals Documentation:** https://chemicals.readthedocs.io/
- **GitHub Repository:** https://github.com/CalebBell/thermo

### Tutorials and Examples
- **Tutorial Notebooks:** https://github.com/CalebBell/thermo/tree/master/docs
- **API Reference:** https://thermo.readthedocs.io/thermo.html
- **Examples Gallery:** https://thermo.readthedocs.io/tutorial.html

### Related Libraries
- **CoolProp:** http://www.coolprop.org/ (Reference fluid properties)
- **Cantera:** https://cantera.org/ (Combustion and thermodynamics)
- **PyFluent:** Ansys Fluent Python interface
- **OpenFOAM:** Open-source CFD with Python bindings

### Reference Books
1. **Poling, B.E., Prausnitz, J.M., O'Connell, J.P.**
   "The Properties of Gases and Liquids" (5th Edition, 2001)
   - Standard reference for property estimation

2. **Perry's Chemical Engineers' Handbook** (9th Edition)
   - Comprehensive chemical engineering reference
   - Property data tables

3. **Yaws, C.L.**
   "Chemical Properties Handbook"
   - Extensive property correlations
   - Temperature-dependent data

4. **Crane Technical Paper 410 (TP-410)**
   "Flow of Fluids Through Valves, Fittings, and Pipe"
   - Fluid properties for pipe flow
   - Pump system design

5. **GPSA Engineering Data Book** (14th Edition)
   - Gas processing applications
   - Property data for hydrocarbons

### Standards Organizations
- **NIST:** National Institute of Standards and Technology
- **AIChE DIPPR:** Design Institute for Physical Properties
- **IAPWS:** International Association for Properties of Water and Steam
- **ASHRAE:** HVAC and refrigeration properties

### Community and Support
- **GitHub Issues:** Report bugs, request features
- **Stack Overflow:** Tag questions with `thermo` and `python`
- **Chemical Engineering Forums:** Discussion of applications

---

## Quick Start Template

```python
"""
Template for thermo package usage in pump applications
"""

from thermo.chemical import Chemical
from thermo.mixture import Mixture

# ===== SINGLE COMPONENT =====
# Define operating conditions
T = 60 + 273.15  # K (convert from °C)
P = 3e5  # Pa (3 bar absolute)

# Create chemical object
fluid = Chemical('water', T=T, P=P)

# Get properties
rho = fluid.rho  # kg/m³
mu = fluid.mu  # Pa·s
Psat = fluid.Psat  # Pa
Cp = fluid.Cp  # J/kg·K

print(f"Density: {rho:.2f} kg/m³")
print(f"Viscosity: {mu*1000:.4f} cP")
print(f"Vapor pressure: {Psat/1000:.2f} kPa")

# ===== MIXTURE =====
# Define composition (mass fractions)
components = ['ethanol', 'water']
mass_fractions = [0.3, 0.7]

# Create mixture
mixture = Mixture(components, ws=mass_fractions, T=T, P=P)

# Get mixture properties
rho_mix = mixture.rho
mu_mix = mixture.mu

print(f"\nMixture density: {rho_mix:.2f} kg/m³")
print(f"Mixture viscosity: {mu_mix*1000:.4f} cP")

# ===== TEMPERATURE SWEEP =====
temperatures = range(20, 101, 10)  # 20°C to 100°C

for T_c in temperatures:
    fluid.T = T_c + 273.15
    print(f"T = {T_c}°C: μ = {fluid.mu*1000:.4f} cP")

# ===== NPSH CALCULATION =====
def calculate_NPSHa(fluid, P_suction, v_suction, z_suction):
    """Calculate Net Positive Suction Head Available"""
    g = 9.81
    pressure_head = (P_suction - fluid.Psat) / (fluid.rho * g)
    velocity_head = v_suction**2 / (2*g)
    return pressure_head + velocity_head - z_suction

NPSHa = calculate_NPSHa(fluid, P_suction=150000, v_suction=2.0, z_suction=1.5)
print(f"\nNPSHa = {NPSHa:.2f} m")
```

---

*Last updated: 2025-11-07*
*For the latest documentation, visit: https://thermo.readthedocs.io/*
