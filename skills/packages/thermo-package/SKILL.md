---
name: thermo-package
description: "Calculate thermodynamic properties and mixture behavior for pump applications"
category: packages
domain: thermodynamics
complexity: intermediate
dependencies:
  - thermo
  - chemicals
---

# Thermo Package Skill

## Overview

The `thermo` library is a comprehensive Python package for chemical and mechanical engineers working with thermodynamic properties and phase equilibria. It provides validated correlations for:

- Pure component thermodynamic properties
- Mixture property calculations
- Vapor-liquid equilibrium (VLE)
- Enthalpy, entropy, and heat capacity calculations
- Temperature-dependent fluid properties
- Equation of state (EOS) implementations
- Phase equilibria for mixtures

The library includes extensive chemical databases and implements industry-standard correlations from DIPPR, NIST, and literature sources. Essential for pump applications involving temperature-dependent properties, mixtures, and multiphase flow.

## Installation

```bash
pip install thermo chemicals
```

The `chemicals` package is a required dependency that provides the chemical database and pure component property correlations.

## Key Modules

### thermo.chemical
Interface for pure component properties including:
- Temperature-dependent viscosity, density, thermal conductivity
- Vapor pressure and phase behavior
- Heat capacity (Cp, Cv)
- Enthalpy and entropy

### thermo.mixture
Mixture property calculations:
- Composition-weighted properties
- Mixing rules for EOS
- Partial molar properties
- Excess properties

### thermo.stream
Process stream calculations:
- Mass and energy balances
- Flash calculations (VLE)
- Phase fraction determination
- Enthalpy of streams

### thermo.eos
Equation of state implementations:
- Peng-Robinson (PR)
- Soave-Redlich-Kwong (SRK)
- Ideal gas law
- Cubic equations of state

## Core Capabilities for Pump Engineering

### 1. Temperature-Dependent Viscosity

Critical for pump performance across operating temperature ranges.

```python
from thermo.chemical import Chemical

# Water viscosity variation with temperature
water = Chemical('water', T=298.15, P=101325)  # 25°C, 1 atm

# At different temperatures
temps = [20, 40, 60, 80, 100]  # °C
for T_c in temps:
    water.T = T_c + 273.15
    print(f"T = {T_c}°C: μ = {water.mu*1000:.4f} cP, ρ = {water.rho:.1f} kg/m³")

# Output demonstrates viscosity decreasing with temperature:
# T = 20°C: μ = 1.0016 cP, ρ = 998.2 kg/m³
# T = 40°C: μ = 0.6529 cP, ρ = 992.2 kg/m³
# T = 60°C: μ = 0.4665 cP, ρ = 983.2 kg/m³
```

### 2. Pump Power Correction for Temperature

Account for viscosity and density changes when pumping at different temperatures.

```python
from thermo.chemical import Chemical

def pump_power_correction(fluid_name, T_design, T_actual, Q, H, eta=0.75):
    """
    Calculate pump power at actual temperature vs design temperature.

    Parameters:
    -----------
    fluid_name : str, chemical name
    T_design : float, design temperature (K)
    T_actual : float, actual temperature (K)
    Q : float, flow rate (m³/s)
    H : float, head (m)
    eta : float, pump efficiency

    Returns:
    --------
    Power correction factor and actual power required
    """
    # Properties at design temperature
    fluid_design = Chemical(fluid_name, T=T_design, P=101325)
    rho_design = fluid_design.rho
    mu_design = fluid_design.mu

    # Properties at actual temperature
    fluid_actual = Chemical(fluid_name, T=T_actual, P=101325)
    rho_actual = fluid_actual.rho
    mu_actual = fluid_actual.mu

    # Power calculation: P = ρ*g*Q*H/η
    g = 9.81  # m/s²
    P_design = rho_design * g * Q * H / eta
    P_actual = rho_actual * g * Q * H / eta

    # Viscosity affects efficiency (simplified)
    # Higher viscosity reduces efficiency
    eta_correction = (mu_design / mu_actual) ** 0.1  # Empirical
    eta_actual = eta * eta_correction
    P_actual_corrected = rho_actual * g * Q * H / eta_actual

    return {
        'P_design': P_design / 1000,  # kW
        'P_actual': P_actual / 1000,  # kW
        'P_corrected': P_actual_corrected / 1000,  # kW
        'rho_ratio': rho_actual / rho_design,
        'mu_ratio': mu_actual / mu_design,
        'eta_correction': eta_correction
    }

# Example: Pumping water at elevated temperature
result = pump_power_correction('water', T_design=298.15, T_actual=353.15,
                                Q=0.05, H=50, eta=0.75)

print(f"Design power: {result['P_design']:.2f} kW")
print(f"Actual power (density only): {result['P_actual']:.2f} kW")
print(f"Actual power (with efficiency): {result['P_corrected']:.2f} kW")
```

### 3. Mixture Properties

Essential for pumps handling mixed fluids or process streams.

```python
from thermo.mixture import Mixture

# Ethanol-water mixture (common in chemical processing)
mixture = Mixture(['ethanol', 'water'], ws=[0.3, 0.7], T=298.15, P=101325)

print(f"Mixture density: {mixture.rho:.1f} kg/m³")
print(f"Mixture viscosity: {mixture.mu*1000:.3f} cP")
print(f"Mixture heat capacity: {mixture.Cp:.1f} J/kg·K")

# Compare to pure components
from thermo.chemical import Chemical
ethanol = Chemical('ethanol', T=298.15, P=101325)
water = Chemical('water', T=298.15, P=101325)

print(f"\nPure ethanol: ρ = {ethanol.rho:.1f} kg/m³, μ = {ethanol.mu*1000:.3f} cP")
print(f"Pure water: ρ = {water.rho:.1f} kg/m³, μ = {water.mu*1000:.3f} cP")
```

### 4. Vapor Pressure and NPSH Calculations

Determine available NPSH and cavitation risk.

```python
from thermo.chemical import Chemical

def calculate_NPSHa(fluid_name, T, P_suction, z_suction, v_suction):
    """
    Calculate Net Positive Suction Head Available.

    Parameters:
    -----------
    fluid_name : str, chemical name
    T : float, fluid temperature (K)
    P_suction : float, suction pressure (Pa absolute)
    z_suction : float, suction elevation relative to pump (m)
    v_suction : float, velocity at suction (m/s)

    Returns:
    --------
    NPSHa in meters
    """
    fluid = Chemical(fluid_name, T=T, P=P_suction)

    rho = fluid.rho  # kg/m³
    Psat = fluid.Psat  # Vapor pressure, Pa

    g = 9.81  # m/s²

    # NPSHa = (P_suction - Psat)/(ρ*g) + v²/(2g) - z_suction
    pressure_head = (P_suction - Psat) / (rho * g)
    velocity_head = v_suction**2 / (2 * g)

    NPSHa = pressure_head + velocity_head - z_suction

    return {
        'NPSHa': NPSHa,
        'P_suction': P_suction / 1000,  # kPa
        'Psat': Psat / 1000,  # kPa
        'rho': rho,
        'margin': (P_suction - Psat) / 1000  # kPa
    }

# Example: Water at 80°C (elevated temperature)
result = calculate_NPSHa('water', T=353.15, P_suction=150000,
                         z_suction=2.0, v_suction=2.5)

print(f"NPSHa = {result['NPSHa']:.2f} m")
print(f"Suction pressure = {result['P_suction']:.1f} kPa")
print(f"Vapor pressure = {result['Psat']:.1f} kPa")
print(f"Pressure margin = {result['margin']:.1f} kPa")

# Warning if NPSHa is low
if result['NPSHa'] < 3.0:
    print("⚠ WARNING: Low NPSHa - high cavitation risk!")
```

### 5. Enthalpy Calculations for Pump Staging

Calculate enthalpy rise and temperature increase in multistage pumps.

```python
from thermo.chemical import Chemical

def pump_temperature_rise(fluid_name, T_inlet, P_inlet, P_outlet, eta):
    """
    Calculate temperature rise due to pump inefficiency.

    Energy not converted to useful work appears as heat in the fluid:
    ΔT = (1-η) * ΔP / (ρ * Cp)

    Parameters:
    -----------
    fluid_name : str
    T_inlet : float, inlet temperature (K)
    P_inlet : float, inlet pressure (Pa)
    P_outlet : float, outlet pressure (Pa)
    eta : float, pump efficiency
    """
    fluid = Chemical(fluid_name, T=T_inlet, P=P_inlet)

    rho = fluid.rho  # kg/m³
    Cp = fluid.Cp  # J/kg·K

    delta_P = P_outlet - P_inlet  # Pa

    # Temperature rise from inefficiency
    delta_T = (1 - eta) * delta_P / (rho * Cp)

    # Recalculate properties at outlet
    T_outlet = T_inlet + delta_T
    fluid_out = Chemical(fluid_name, T=T_outlet, P=P_outlet)

    # Enthalpy change
    H_inlet = fluid.H  # J/kg
    H_outlet = fluid_out.H  # J/kg
    delta_H = H_outlet - H_inlet

    return {
        'T_inlet': T_inlet - 273.15,  # °C
        'T_outlet': T_outlet - 273.15,  # °C
        'delta_T': delta_T,  # K
        'delta_H': delta_H / 1000,  # kJ/kg
        'delta_P': delta_P / 1e6  # MPa
    }

# Example: Multistage boiler feed pump
result = pump_temperature_rise('water', T_inlet=333.15, P_inlet=500000,
                               P_outlet=15000000, eta=0.82)

print(f"Inlet temperature: {result['T_inlet']:.1f} °C")
print(f"Outlet temperature: {result['T_outlet']:.1f} °C")
print(f"Temperature rise: {result['delta_T']:.2f} K")
print(f"Enthalpy rise: {result['delta_H']:.1f} kJ/kg")
print(f"Pressure rise: {result['delta_P']:.1f} MPa")
```

### 6. Phase Equilibrium for Multiphase Pumps

Determine if two-phase flow will occur in pump.

```python
from thermo.chemical import Chemical

def check_phase_state(fluid_name, T, P):
    """
    Check if fluid is single-phase or will flash to vapor.
    """
    fluid = Chemical(fluid_name, T=T, P=P)

    Psat = fluid.Psat  # Vapor pressure
    Tsat = fluid.Tb  # Boiling point at 1 atm

    # Determine phase
    if P > Psat:
        phase = "Liquid (subcooled)"
        margin = (P - Psat) / Psat * 100  # % margin
    elif P < Psat:
        phase = "Two-phase or Vapor (FLASHING!)"
        margin = (Psat - P) / P * 100  # % above operating pressure
    else:
        phase = "Saturated liquid"
        margin = 0

    return {
        'phase': phase,
        'T': T - 273.15,
        'P': P / 1000,  # kPa
        'Psat': Psat / 1000,  # kPa
        'margin': margin
    }

# Example: Check if water will flash at pump suction
temps = [40, 60, 80, 100]  # °C
P_suction = 120000  # Pa (1.2 bar absolute)

for T_c in temps:
    result = check_phase_state('water', T=T_c+273.15, P=P_suction)
    print(f"T = {result['T']:.0f}°C, P = {result['P']:.1f} kPa: {result['phase']}")
    if 'subcooled' in result['phase']:
        print(f"  Margin = {result['margin']:.1f}%")
    print()

# Output will show flashing occurs above ~105°C at 1.2 bar
```

## Equation of State Applications

### Peng-Robinson EOS for High-Pressure Pumps

```python
from thermo import ChemicalConstantsPackage, PRMIX, CEOSGas, CEOSLiquid, FlashVL
from thermo.interaction_parameters import IPDB

# High-pressure natural gas mixture
constants = ChemicalConstantsPackage.constants_from_IDs(['methane', 'ethane', 'propane'])

# Operating conditions
T = 298.15  # K
P = 10e6  # Pa (100 bar)
zs = [0.85, 0.10, 0.05]  # Mole fractions

# Peng-Robinson EOS
kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')
eos_kwargs = {'Tcs': constants.Tcs, 'Pcs': constants.Pcs, 'omegas': constants.omegas, 'kijs': kijs}

gas = CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=constants.HeatCapacityGases, T=T, P=P, zs=zs)
liquid = CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=constants.HeatCapacityGases, T=T, P=P, zs=zs)

# Flash calculation
flasher = FlashVL(constants, HeatCapacityGases=constants.HeatCapacityGases, gas=gas, liquids=[liquid])
res = flasher.flash(T=T, P=P, zs=zs)

print(f"Temperature: {T-273.15:.1f} °C")
print(f"Pressure: {P/1e6:.1f} MPa")
print(f"Phase: {res.phase}")
if res.phase == 'L':
    print(f"Liquid density: {res.rho_mass():.1f} kg/m³")
    print(f"Liquid viscosity: {res.mu()*1e6:.2f} μPa·s")
```

## Complete Engineering Examples

See `examples.py` for verified, production-ready examples including:

1. **Temperature-viscosity correction for pump curves**
2. **NPSH calculations with temperature variation**
3. **Mixture property calculations for blended fuels**
4. **Enthalpy rise in multistage pumps**
5. **Flash calculation for pump discharge conditions**
6. **Chemical database queries**

## Pump Application Summary

### When to Use Thermo Package

1. **Temperature-dependent properties**: When fluid temperature varies significantly
2. **Mixture handling**: Pumping blended fluids or process streams
3. **High-temperature applications**: Boiler feed pumps, thermal oil pumps
4. **Cavitation analysis**: Accurate vapor pressure for NPSH calculations
5. **Multiphase flow**: Determining if flashing will occur
6. **Energy calculations**: Heat balances for pump temperature rise

### Integration with Fluids Package

The `thermo` package complements `fluids` for complete pump analysis:

- **thermo**: Provides temperature-dependent properties (ρ, μ, Psat)
- **fluids**: Uses these properties for hydraulic calculations (friction, head loss)

```python
from thermo.chemical import Chemical
from fluids.core import Reynolds
from fluids.friction import friction_factor

# Get properties from thermo
water = Chemical('water', T=333.15, P=101325)  # 60°C
rho = water.rho
mu = water.mu

# Use in fluids calculations
Re = Reynolds(V=2.0, D=0.1, rho=rho, mu=mu)
f = friction_factor(Re=Re, eD=0.0001)

print(f"At 60°C: Re = {Re:.0f}, f = {f:.5f}")
```

## Common Chemicals Database

Access to 20,000+ chemicals via CAS number or name:
- Industrial fluids (water, oils, glycols)
- Hydrocarbons (methane, propane, gasoline)
- Refrigerants (R-134a, ammonia)
- Acids and bases (HCl, NaOH)
- Organic solvents (ethanol, acetone)

See `reference.md` for database details and common pump fluids.

## Best Practices

1. **Always specify both T and P** when creating Chemical objects
2. **Check phase state** before using liquid properties
3. **Validate against known data** (steam tables, RefProp, NIST)
4. **Use appropriate EOS** for high-pressure applications (>20 bar)
5. **Consider temperature range** of correlations (some limited to <150°C)
6. **Cache Chemical objects** for performance in iterative calculations

## Units

All properties in SI units:
- Temperature: K (convert from °C: T_K = T_C + 273.15)
- Pressure: Pa (1 bar = 100,000 Pa)
- Density: kg/m³
- Viscosity: Pa·s (1 cP = 0.001 Pa·s)
- Enthalpy: J/kg
- Heat capacity: J/kg·K

## Troubleshooting

### Issue: Property returns None
- **Cause**: Property not available for this chemical or outside valid range
- **Solution**: Check Chemical.Tmin, Chemical.Tmax for valid temperature range

### Issue: Mixture properties unrealistic
- **Cause**: Missing interaction parameters or invalid composition
- **Solution**: Verify sum of mass/mole fractions equals 1.0

### Issue: Flash calculation fails
- **Cause**: Initial guess outside valid region or bad EOS parameters
- **Solution**: Use simpler EOS (ideal gas) or adjust initial T/P estimate

## References

- DIPPR 801 Database (AIChE Design Institute for Physical Properties)
- NIST Chemistry WebBook
- Poling, B.E., Prausnitz, J.M., O'Connell, J.P., "The Properties of Gases and Liquids" (5th Edition)
- Perry's Chemical Engineers' Handbook (9th Edition)
- Peng, D.Y., Robinson, D.B., "A New Two-Constant Equation of State" (1976)

## Further Reading

- Official documentation: https://thermo.readthedocs.io/
- Chemicals documentation: https://chemicals.readthedocs.io/
- Source code: https://github.com/CalebBell/thermo
- Tutorial notebooks: https://github.com/CalebBell/thermo/tree/master/docs
