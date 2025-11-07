---
name: pint-units
description: "Handle engineering units with automatic conversion and dimensional analysis"
category: packages
domain: general
complexity: basic
dependencies:
  - pint
---

# Pint Units Skill

## Overview

`Pint` is a Python package for handling physical quantities with units. It provides:

- Automatic unit conversions
- Dimensional analysis and consistency checking
- Unit-aware arithmetic operations
- Integration with NumPy arrays
- Support for custom unit definitions
- Temperature conversions (including offset units)
- Scientific notation and formatting

Pint eliminates unit conversion errors and ensures dimensional consistency in engineering calculations, making it essential for any technical work involving physical quantities.

## Installation

```bash
pip install pint
```

For integration with NumPy and Pandas:

```bash
pip install pint[numpy]
```

## Basic Usage

### Creating Quantities

A quantity in Pint consists of a magnitude (number) and a unit.

```python
from pint import UnitRegistry

# Create unit registry (do this once per module)
ureg = UnitRegistry()

# Define quantities with units
distance = 100 * ureg.meter
time = 5 * ureg.second
pressure = 50 * ureg.psi

# Alternative syntax using string parsing
flow_rate = ureg('250 gallons/minute')
temperature = ureg('75 degF')
viscosity = ureg('1.5 centipoise')

# Access magnitude and units
print(f"Distance magnitude: {distance.magnitude}")  # 100
print(f"Distance units: {distance.units}")  # meter
```

### Unit Conversions

Convert between compatible units automatically:

```python
# Length conversions
distance = 100 * ureg.meter
print(distance.to('feet'))  # 328.084 foot
print(distance.to('kilometer'))  # 0.1 kilometer

# Pressure conversions
pressure = 100 * ureg.psi
print(pressure.to('bar'))  # 6.895 bar
print(pressure.to('pascal'))  # 689475.7 pascal

# Flow rate conversions
flow = 500 * ureg.gallons / ureg.minute
print(flow.to('liter/second'))  # 31.546 liter/second
print(flow.to('m**3/hour'))  # 113.562 meter**3/hour

# Temperature conversions (handles offset units)
temp = 75 * ureg.degF
print(temp.to('degC'))  # 23.889 degree_Celsius
print(temp.to('kelvin'))  # 297.039 kelvin
```

## Dimensional Analysis

Pint ensures dimensional consistency in calculations:

```python
# Velocity = Distance / Time
distance = 100 * ureg.meter
time = 5 * ureg.second
velocity = distance / time
print(velocity)  # 20.0 meter/second

# Force = Mass × Acceleration
mass = 10 * ureg.kg
acceleration = 9.81 * ureg.meter / ureg.second**2
force = mass * acceleration
print(force)  # 98.1 kilogram·meter/second²
print(force.to('newton'))  # 98.1 newton

# Dimensional consistency check
try:
    invalid = (100 * ureg.meter) + (50 * ureg.second)  # ❌ ERROR
except Exception as e:
    print(f"Error: Cannot add length and time")

# Correct operation
total_distance = (100 * ureg.meter) + (50 * ureg.feet)  # ✓ OK
print(total_distance)  # 115.24 meter
```

## Unit-Aware Calculations

### Fluid Flow Calculations

```python
# Calculate volumetric flow rate from velocity and area
velocity = 3.5 * ureg.meter / ureg.second
diameter = 150 * ureg.millimeter
area = 3.14159 * (diameter/2)**2

flow_rate = velocity * area
print(flow_rate.to('liter/second'))  # 61.86 liter/second
print(flow_rate.to('gpm'))  # 1638.1 gallon/minute

# Reynolds number calculation
density = 998 * ureg.kg / ureg.meter**3
viscosity = 1.0 * ureg.centipoise  # Common unit in industry
Re = (density * velocity * diameter) / viscosity
print(Re.to_base_units())  # Dimensionless: 523170.0
```

### Pump Power Calculations

```python
# Hydraulic power: P = ρ × g × Q × H
flow_rate = 100 * ureg.meter**3 / ureg.hour
head = 50 * ureg.meter
density = 1000 * ureg.kg / ureg.meter**3
gravity = 9.81 * ureg.meter / ureg.second**2

hydraulic_power = density * gravity * flow_rate * head
print(hydraulic_power.to('kilowatt'))  # 13.625 kilowatt
print(hydraulic_power.to('horsepower'))  # 18.27 horsepower

# With pump efficiency
efficiency = 0.75 * ureg.dimensionless
shaft_power = hydraulic_power / efficiency
print(shaft_power.to('kilowatt'))  # 18.17 kilowatt
```

### Pressure Drop Calculations

```python
# Darcy-Weisbach equation: ΔP = f × (L/D) × (ρV²/2)
friction_factor = 0.018 * ureg.dimensionless
length = 100 * ureg.meter
diameter = 0.15 * ureg.meter
velocity = 2.5 * ureg.meter / ureg.second
density = 1000 * ureg.kg / ureg.meter**3

pressure_drop = friction_factor * (length/diameter) * (density * velocity**2 / 2)
print(pressure_drop.to('pascal'))  # 37500.0 pascal
print(pressure_drop.to('psi'))  # 5.44 pound_force_per_square_inch
print(pressure_drop.to('bar'))  # 0.375 bar

# Convert to head loss
head_loss = pressure_drop / (density * gravity)
print(head_loss.to('meter'))  # 3.822 meter
print(head_loss.to('feet'))  # 12.54 foot
```

## Working with Arrays

Pint integrates seamlessly with NumPy for array operations:

```python
import numpy as np
from pint import UnitRegistry

ureg = UnitRegistry()

# Create array with units
flow_rates = np.array([50, 75, 100, 125, 150]) * ureg.gpm
heads = np.array([80, 75, 65, 50, 30]) * ureg.meter

# Array operations preserve units
powers = (flow_rates * heads * ureg.kg/ureg.meter**3 * 9.81*ureg.meter/ureg.second**2)
print(powers.to('kilowatt'))
# [2.48 3.72 4.95 6.19 7.43] kilowatt

# Statistical operations
mean_flow = np.mean(flow_rates)
std_flow = np.std(flow_rates)
print(f"Mean flow: {mean_flow.to('liter/minute'):.1f}")
print(f"Std dev: {std_flow.to('liter/minute'):.1f}")

# Unit-aware interpolation
target_head = 70 * ureg.meter
target_flow = np.interp(target_head.magnitude, heads[::-1].magnitude,
                        flow_rates[::-1].magnitude) * ureg.gpm
print(f"Flow at {target_head}: {target_flow.to('m**3/hour'):.1f}")
```

## Custom Unit Definitions

Define domain-specific units:

```python
# Add custom units to registry
ureg.define('barrel_oil = 42 * gallon = bbl')
ureg.define('standard_cubic_foot = foot**3 = scf')
ureg.define('darcy = centipoise * centimeter**2 / (second * atmosphere) = D')

# Use custom units
oil_volume = 1000 * ureg.barrel_oil
print(oil_volume.to('gallon'))  # 42000.0 gallon
print(oil_volume.to('liter'))  # 158987.3 liter

gas_flow = 5000 * ureg.standard_cubic_foot / ureg.day
print(gas_flow.to('m**3/hour'))  # 5.90 meter**3/hour

permeability = 100 * ureg.darcy
print(permeability.to_base_units())  # Base SI units
```

## Context Managers for Unit Systems

Switch between unit systems easily:

```python
from pint import UnitRegistry

ureg = UnitRegistry()

length = 100 * ureg.meter
mass = 50 * ureg.kg

# Default SI units
print(f"Length: {length}")  # 100 meter
print(f"Mass: {mass}")  # 50 kilogram

# Use US customary units
with ureg.context('US'):
    print(f"Length: {length.to('feet')}")  # 328.084 foot
    print(f"Mass: {mass.to('pound')}")  # 110.231 pound

# Use imperial units
with ureg.context('imperial'):
    print(f"Length: {length.to('yard')}")  # 109.361 yard
```

## Temperature Conversions

Handle absolute and relative temperature correctly:

```python
# Absolute temperatures
temp1 = 25 * ureg.degC
temp2 = temp1.to('degF')
print(temp2)  # 77.0 degree_Fahrenheit

temp3 = 300 * ureg.kelvin
print(temp3.to('degC'))  # 26.85 degree_Celsius

# Temperature differences (delta)
temp_rise = ureg.Quantity(20, ureg.delta_degC)
print(temp_rise.to(ureg.delta_degF))  # 36.0 delta_degree_Fahrenheit

# Heat transfer calculation
mass = 10 * ureg.kg
specific_heat = 4.18 * ureg.kJ / (ureg.kg * ureg.kelvin)
temp_change = 50 * ureg.delta_degC

heat = mass * specific_heat * temp_change
print(heat.to('kJ'))  # 2090.0 kilojoule
print(heat.to('BTU'))  # 1981.9 BTU
```

## Formatting and Display

Control how quantities are displayed:

```python
pressure = 150000 * ureg.pascal

# Default format
print(pressure)  # 150000 pascal

# Compact notation
print(f"{pressure:~}")  # 150000 Pa

# Specific precision
print(f"{pressure:.2f}")  # 150000.00 pascal

# Scientific notation
print(f"{pressure:~.2e}")  # 1.50e+05 Pa

# Pretty format
print(f"{pressure:~P}")  # 150000 Pa

# Custom format
print(f"{pressure.to('bar'):.3f~P}")  # 1.500 bar
```

## Common Engineering Unit Conversions

### Pressure

```python
pressure = 100 * ureg.psi
# Common conversions
print(pressure.to('bar'))  # 6.895 bar
print(pressure.to('kPa'))  # 689.476 kilopascal
print(pressure.to('MPa'))  # 0.689 megapascal
print(pressure.to('atm'))  # 6.805 atmosphere
print(pressure.to('mmHg'))  # 5171.5 millimeter_Hg
print(pressure.to('inch_H2O'))  # 2767.7 inch_H2O
```

### Flow Rate

```python
flow = 100 * ureg.gpm  # gallons per minute
# Common conversions
print(flow.to('liter/minute'))  # 378.541 liter/minute
print(flow.to('m**3/hour'))  # 22.712 meter**3/hour
print(flow.to('ft**3/second'))  # 0.223 foot**3/second
print(flow.to('barrel_oil/day'))  # 3428.57 barrel_oil/day (if defined)
```

### Viscosity

```python
# Dynamic viscosity
mu = 1.0 * ureg.centipoise
print(mu.to('pascal*second'))  # 0.001 pascal·second
print(mu.to('lbf*second/ft**2'))  # 2.089e-05 pound_force·second/foot²

# Kinematic viscosity
nu = 1.0 * ureg.centistokes
print(nu.to('m**2/second'))  # 1e-06 meter²/second
print(nu.to('ft**2/second'))  # 1.076e-05 foot²/second
```

### Energy and Power

```python
energy = 100 * ureg.kWh
print(energy.to('MJ'))  # 360.0 megajoule
print(energy.to('BTU'))  # 341214.2 BTU
print(energy.to('therm'))  # 3.412 therm

power = 50 * ureg.horsepower
print(power.to('kilowatt'))  # 37.285 kilowatt
print(power.to('BTU/hour'))  # 127259.0 BTU/hour
```

### Mass Flow Rate

```python
mass_flow = 1000 * ureg.kg / ureg.hour
print(mass_flow.to('lb/minute'))  # 36.74 pound/minute
print(mass_flow.to('ton/day'))  # 0.024 metric_ton/day
print(mass_flow.to('g/second'))  # 277.78 gram/second
```

## Best Practices

1. **Create UnitRegistry Once**: Define `ureg = UnitRegistry()` at module level, not inside functions
2. **Use Base Units for Storage**: Store values in consistent base units (SI) in databases
3. **Validate Input Units**: Always check that input quantities have expected dimensionality
4. **Handle Dimensionless Quantities**: Use `ureg.dimensionless` for unitless ratios
5. **Convert at Boundaries**: Convert to display units only when presenting results
6. **Check Compatibility**: Use `quantity.check(dimension)` to verify dimensional consistency
7. **Be Explicit**: Use explicit unit definitions rather than assuming defaults
8. **Temperature Care**: Use `delta_` prefix for temperature differences vs absolute temperatures

## Dimensional Consistency Checking

```python
def calculate_reynolds_number(velocity, diameter, density, viscosity):
    """
    Calculate Reynolds number with automatic dimensional checking.

    Re = ρVD/μ (dimensionless)
    """
    # Pint automatically checks dimensions
    Re = (density * velocity * diameter) / viscosity

    # Verify result is dimensionless
    assert Re.dimensionality == ureg.dimensionless.dimensionality

    # Return magnitude (pure number)
    return Re.to_base_units().magnitude

# Correct usage
rho = 1000 * ureg.kg / ureg.meter**3
V = 2.5 * ureg.meter / ureg.second
D = 0.15 * ureg.meter
mu = 1e-3 * ureg.pascal * ureg.second

Re = calculate_reynolds_number(V, D, rho, mu)
print(f"Reynolds number: {Re:.0f}")  # 375000

# Incorrect usage will raise error
try:
    Re_wrong = calculate_reynolds_number(V, D, rho, rho)  # Wrong dimension
except Exception as e:
    print("Error: Dimensional inconsistency detected")
```

## Integration with Engineering Workflows

### Example: Pump Performance Curve

```python
import numpy as np
from pint import UnitRegistry

ureg = UnitRegistry()

def pump_curve(flow_rates, coefficients):
    """
    Calculate pump head from flow rate using curve fit.
    H = H0 - A*Q - B*Q²

    Parameters
    ----------
    flow_rates : Quantity array
        Flow rates with units
    coefficients : dict
        'H0', 'A', 'B' with appropriate units

    Returns
    -------
    heads : Quantity array
        Pump heads with units
    """
    H0 = coefficients['H0']
    A = coefficients['A']
    B = coefficients['B']

    heads = H0 - A * flow_rates - B * flow_rates**2
    return heads

# Define pump curve coefficients
coeffs = {
    'H0': 80 * ureg.meter,
    'A': 200 * ureg.meter / (ureg.meter**3/ureg.second),
    'B': 3000 * ureg.meter / (ureg.meter**3/ureg.second)**2
}

# Calculate performance at different flow rates
Q = np.linspace(0, 0.1, 11) * ureg.meter**3 / ureg.second
H = pump_curve(Q, coeffs)

# Display in different units
print("Flow (GPM)  Head (ft)  Head (m)")
print("-" * 40)
for q, h in zip(Q, H):
    print(f"{q.to('gpm'):8.0f~P}  {h.to('feet'):8.1f~P}  {h.to('meter'):7.1f~P}")
```

## Troubleshooting

### Issue: DimensionalityError
**Cause**: Attempting to add/compare quantities with incompatible units
**Solution**: Check that all terms have the same dimensionality, or convert explicitly

### Issue: UndefinedUnitError
**Cause**: Using a unit that's not defined in the registry
**Solution**: Define custom units using `ureg.define()` or check spelling

### Issue: Offset units (temperature) errors
**Cause**: Mixing absolute and relative temperature units
**Solution**: Use `delta_degC` for temperature differences, plain `degC` for absolute

### Issue: Lost units in calculations
**Cause**: Using `.magnitude` too early in calculations
**Solution**: Keep quantities as Pint objects until final output

## Quick Reference Card

### Common Units

| Quantity | Units |
|----------|-------|
| Length | meter, foot, inch, mile, kilometer |
| Area | meter**2, foot**2, acre, hectare |
| Volume | liter, gallon, barrel_oil, ft**3, m**3 |
| Mass | kilogram, pound, ton, tonne |
| Time | second, minute, hour, day |
| Velocity | meter/second, ft/second, mph, kph |
| Flow (Vol) | m**3/hour, gpm, liter/minute, ft**3/second |
| Flow (Mass) | kg/second, lb/minute, ton/hour |
| Pressure | pascal, bar, psi, atm, mmHg |
| Force | newton, lbf, kgf |
| Power | watt, horsepower, BTU/hour |
| Energy | joule, kWh, BTU, calorie |
| Temperature | degC, degF, kelvin, rankine |
| Viscosity (dyn) | pascal*second, poise, centipoise |
| Viscosity (kin) | m**2/second, stokes, centistokes |

### Quick Conversions

```python
# Create quantity
Q = 100 * ureg.gpm

# Convert
Q.to('liter/minute')  # To specific unit
Q.to_base_units()     # To SI base units
Q.to_compact()        # To compact SI (kilo, mega, etc)

# Check dimensionality
Q.dimensionality      # [length]³/[time]
Q.check('[volume]/[time]')  # Verify dimension

# Access components
Q.magnitude          # Numeric value
Q.units              # Unit object
```

## References

- Official documentation: https://pint.readthedocs.io/
- Source code: https://github.com/hgrecco/pint
- Unit definitions: https://github.com/hgrecco/pint/blob/master/pint/default_en.txt
- NIST Guide to SI: https://www.nist.gov/pml/special-publication-811

## Further Reading

- NIST Special Publication 811: Guide for the Use of the International System of Units (SI)
- ISO 80000: Quantities and units
- API Standards for petroleum industry units
- ASME standards for engineering calculations
