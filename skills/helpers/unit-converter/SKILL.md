---
name: unit-converter
description: "Convert flow rates, pressures, viscosities, and other engineering units"
category: helpers
domain: general
complexity: basic
dependencies: [pint]
---

# Unit Converter Skill

A comprehensive unit conversion tool for engineering work using the Pint library. This skill handles common engineering unit conversions with dimensional analysis and consistency checking.

## Overview of Pint Library

Pint is a Python library for working with physical quantities. It provides:

- **Unit Registry**: Central registry of units and their relationships
- **Quantity Objects**: Combines numerical values with units
- **Automatic Conversion**: Converts between compatible units
- **Dimensional Analysis**: Verifies dimensional consistency
- **Custom Units**: Allows definition of domain-specific units
- **Context Management**: Handles different unit systems (SI, Imperial, etc.)

### Installation
```bash
pip install pint
```

### Basic Usage
```python
import pint

ureg = pint.UnitRegistry()

# Create a quantity
flow = 100 * ureg.liter / ureg.minute

# Convert to another unit
flow_gpm = flow.to('gallon/minute')
print(f"{flow_gpm:.2f}")  # Output: 26.42 gallon / minute
```

## Common Engineering Unit Conversions

### Flow Rate Conversions

Flow rates are fundamental in fluid mechanics and process engineering.

**Common Units:**
- **m³/s** (cubic meters per second): SI standard
- **m³/h** (cubic meters per hour): Common in HVAC
- **L/min** (liters per minute): Small flow applications
- **gpm** (gallons per minute): US standard for liquid flow
- **cfm** (cubic feet per minute): US standard for gas flow

**Conversion Examples:**
```python
# Water flow in a pipe
flow = 500 * ureg.liter / ureg.minute
print(f"Flow: {flow.to('m^3/s'):.4f}")        # 0.0083 m³/s
print(f"Flow: {flow.to('m^3/hour'):.2f}")     # 30.00 m³/h
print(f"Flow: {flow.to('gallon/minute'):.2f}") # 132.09 gpm

# Air flow in ventilation
air_flow = 1000 * ureg.cfm
print(f"Air flow: {air_flow.to('m^3/s'):.3f}") # 0.472 m³/s
print(f"Air flow: {air_flow.to('L/s'):.1f}")   # 471.9 L/s
```

### Pressure Conversions

Pressure measurements vary widely across industries and regions.

**Common Units:**
- **Pa** (Pascal): SI base unit (N/m²)
- **kPa** (kiloPascal): Common engineering unit
- **bar**: Metric unit (100 kPa)
- **psi** (pounds per square inch): US standard
- **mmHg** (millimeters of mercury): Medical, vacuum
- **atm** (atmosphere): Standard atmospheric pressure

**Conversion Examples:**
```python
# Hydraulic system pressure
pressure = 150 * ureg.bar
print(f"Pressure: {pressure.to('Pa'):.0f}")     # 15000000 Pa
print(f"Pressure: {pressure.to('kPa'):.0f}")    # 15000 kPa
print(f"Pressure: {pressure.to('psi'):.1f}")    # 2176.0 psi
print(f"Pressure: {pressure.to('atm'):.1f}")    # 147.9 atm

# Vacuum measurement
vacuum = 500 * ureg.mmHg
print(f"Vacuum: {vacuum.to('Pa'):.0f}")         # 66661 Pa
print(f"Vacuum: {vacuum.to('kPa'):.2f}")        # 66.66 kPa
print(f"Vacuum: {vacuum.to('psi'):.2f}")        # 9.67 psi

# Tire pressure
tire = 32 * ureg.psi
print(f"Tire: {tire.to('kPa'):.1f}")            # 220.6 kPa
print(f"Tire: {tire.to('bar'):.2f}")            # 2.21 bar
```

### Viscosity Conversions

Viscosity can be expressed as dynamic or kinematic viscosity.

**Dynamic Viscosity:**
- **Pa·s** (Pascal-second): SI unit
- **cP** (centipoise): Common in industry (1 cP = 1 mPa·s)
- **P** (Poise): CGS unit (1 P = 0.1 Pa·s)

**Kinematic Viscosity:**
- **m²/s**: SI unit
- **cSt** (centistokes): Common in industry
- **St** (Stokes): CGS unit

**Conversion Examples:**
```python
# Dynamic viscosity - Water at 20°C
viscosity_dynamic = 1.002 * ureg.centipoise
print(f"Dynamic: {viscosity_dynamic.to('Pa*s'):.6f}")  # 0.001002 Pa·s
print(f"Dynamic: {viscosity_dynamic.to('poise'):.5f}") # 0.01002 P

# Dynamic viscosity - Motor oil
oil_dynamic = 250 * ureg.centipoise
print(f"Oil: {oil_dynamic.to('Pa*s'):.3f}")            # 0.250 Pa·s

# Kinematic viscosity - Water at 20°C
viscosity_kinematic = 1.004 * ureg.centistokes
print(f"Kinematic: {viscosity_kinematic.to('m^2/s'):.9f}") # 0.000001004 m²/s
print(f"Kinematic: {viscosity_kinematic.to('stokes'):.5f}") # 0.01004 St

# Kinematic viscosity - Motor oil SAE 30
oil_kinematic = 100 * ureg.centistokes
print(f"Oil SAE 30: {oil_kinematic.to('m^2/s'):.7f}")  # 0.0000100 m²/s
```

**Note on Viscosity Conversion:**
To convert between dynamic and kinematic viscosity:
```
Kinematic Viscosity = Dynamic Viscosity / Density
ν (m²/s) = μ (Pa·s) / ρ (kg/m³)
```

### Length Conversions

**Common Units:**
- **m** (meter): SI base unit
- **mm** (millimeter): Engineering drawings
- **cm** (centimeter): General measurements
- **ft** (foot): US standard
- **in** (inch): US manufacturing

**Conversion Examples:**
```python
# Pipe diameter
diameter = 4 * ureg.inch
print(f"Diameter: {diameter.to('mm'):.1f}")     # 101.6 mm
print(f"Diameter: {diameter.to('m'):.4f}")      # 0.1016 m
print(f"Diameter: {diameter.to('ft'):.3f}")     # 0.333 ft

# Building dimensions
length = 15.5 * ureg.meter
print(f"Length: {length.to('ft'):.2f}")         # 50.85 ft
print(f"Length: {length.to('inch'):.1f}")       # 610.2 in

# Sheet thickness
thickness = 0.5 * ureg.mm
print(f"Thickness: {thickness.to('inch'):.4f}") # 0.0197 in
print(f"Thickness: {thickness.to('m'):.6f}")    # 0.000500 m
```

### Power Conversions

**Common Units:**
- **W** (Watt): SI unit
- **kW** (kilowatt): Common engineering unit
- **HP** (horsepower): Mechanical power (1 HP = 745.7 W)
- **MW** (megawatt): Large power systems

**Conversion Examples:**
```python
# Motor power
motor = 50 * ureg.horsepower
print(f"Motor: {motor.to('kW'):.2f}")           # 37.28 kW
print(f"Motor: {motor.to('W'):.0f}")            # 37285 W

# Pump power
pump = 15 * ureg.kW
print(f"Pump: {pump.to('hp'):.2f}")             # 20.11 hp
print(f"Pump: {pump.to('W'):.0f}")              # 15000 W

# Generator output
generator = 2.5 * ureg.MW
print(f"Generator: {generator.to('kW'):.0f}")   # 2500 kW
print(f"Generator: {generator.to('hp'):.1f}")   # 3352.0 hp
```

### Temperature Conversions

**Common Units:**
- **K** (Kelvin): SI base unit, absolute scale
- **°C** (Celsius): Metric, water freezes at 0°C
- **°F** (Fahrenheit): US standard, water freezes at 32°F

**Conversion Examples:**
```python
# Process temperature
temp_c = 85 * ureg.degC
print(f"Temperature: {temp_c.to('degF'):.1f}")  # 185.0 °F
print(f"Temperature: {temp_c.to('kelvin'):.2f}") # 358.15 K

# Cryogenic temperature
temp_k = 77 * ureg.kelvin  # Liquid nitrogen
print(f"LN2: {temp_k.to('degC'):.2f}")          # -196.15 °C
print(f"LN2: {temp_k.to('degF'):.2f}")          # -321.07 °F

# Ambient temperature
temp_f = 72 * ureg.degF
print(f"Room: {temp_f.to('degC'):.1f}")         # 22.2 °C
print(f"Room: {temp_f.to('kelvin'):.2f}")       # 295.37 K
```

### Additional Engineering Units

**Area:**
```python
area = 500 * ureg.meter**2
print(f"Area: {area.to('ft^2'):.1f}")           # 5381.9 ft²
print(f"Area: {area.to('hectare'):.2f}")        # 0.05 ha

# Pipe cross-section
pipe_area = 3.14 * ureg.inch**2
print(f"Pipe area: {pipe_area.to('cm^2'):.2f}") # 20.26 cm²
```

**Volume:**
```python
tank = 5000 * ureg.liter
print(f"Tank: {tank.to('m^3'):.1f}")            # 5.0 m³
print(f"Tank: {tank.to('gallon'):.1f}")         # 1320.9 gal
print(f"Tank: {tank.to('ft^3'):.2f}")           # 176.57 ft³
```

**Velocity:**
```python
velocity = 2.5 * ureg.meter / ureg.second
print(f"Velocity: {velocity.to('ft/s'):.2f}")   # 8.20 ft/s
print(f"Velocity: {velocity.to('km/h'):.1f}")   # 9.0 km/h
print(f"Velocity: {velocity.to('mph'):.2f}")    # 5.59 mph
```

**Mass Flow:**
```python
mass_flow = 100 * ureg.kg / ureg.hour
print(f"Mass flow: {mass_flow.to('kg/s'):.5f}") # 0.02778 kg/s
print(f"Mass flow: {mass_flow.to('lb/min'):.3f}") # 3.667 lb/min
```

**Density:**
```python
density = 1000 * ureg.kg / ureg.meter**3
print(f"Density: {density.to('g/cm^3'):.1f}")   # 1.0 g/cm³
print(f"Density: {density.to('lb/ft^3'):.2f}")  # 62.43 lb/ft³
```

**Energy:**
```python
energy = 1000 * ureg.joule
print(f"Energy: {energy.to('kJ'):.1f}")         # 1.0 kJ
print(f"Energy: {energy.to('BTU'):.3f}")        # 0.948 BTU
print(f"Energy: {energy.to('kWh'):.6f}")        # 0.000278 kWh
print(f"Energy: {energy.to('cal'):.1f}")        # 239.0 cal
```

## Dimensional Analysis

Dimensional analysis ensures that equations are dimensionally consistent.

### Checking Dimensional Consistency

```python
# Reynolds number (dimensionless)
rho = 1000 * ureg.kg / ureg.m**3      # Density
v = 2 * ureg.m / ureg.s               # Velocity
D = 0.1 * ureg.m                       # Diameter
mu = 0.001 * ureg.Pa * ureg.s         # Dynamic viscosity

Re = (rho * v * D / mu).to_base_units()
print(f"Reynolds number: {Re:.0f}")    # Dimensionless: 200000

# Verify it's dimensionless
print(f"Dimensionality: {Re.dimensionality}")  # dimensionless
```

### Common Dimensional Relationships

```python
# Verify force = mass × acceleration
mass = 10 * ureg.kg
accel = 9.81 * ureg.m / ureg.s**2
force = mass * accel
print(f"Force: {force.to('N'):.1f}")   # 98.1 N
print(f"Dimensionality: {force.dimensionality}")  # [length] * [mass] / [time]^2

# Verify pressure = force / area
force = 1000 * ureg.newton
area = 0.01 * ureg.m**2
pressure = force / area
print(f"Pressure: {pressure.to('kPa'):.0f}")  # 100 kPa
print(f"Dimensionality: {pressure.dimensionality}")  # [mass] / ([length] * [time]^2)

# Verify power = force × velocity
force = 500 * ureg.N
velocity = 10 * ureg.m / ureg.s
power = force * velocity
print(f"Power: {power.to('kW'):.1f}")  # 5.0 kW
print(f"Dimensionality: {power.dimensionality}")  # [length]^2 * [mass] / [time]^3
```

## Unit Consistency Checking

Pint automatically checks unit consistency and raises errors for invalid operations.

### Valid Operations

```python
# Adding compatible units
length1 = 10 * ureg.meter
length2 = 5 * ureg.feet
total_length = length1 + length2.to('meter')
print(f"Total: {total_length:.3f}")    # 11.524 m

# Multiplying different units
area = (5 * ureg.meter) * (3 * ureg.meter)
print(f"Area: {area:.0f}")             # 15 m²

# Dividing units
velocity = (100 * ureg.meter) / (10 * ureg.second)
print(f"Velocity: {velocity:.0f}")     # 10 m/s
```

### Invalid Operations (Will Raise Errors)

```python
# Cannot add incompatible units
try:
    result = (10 * ureg.meter) + (5 * ureg.second)
except pint.errors.DimensionalityError as e:
    print(f"Error: Cannot add length and time")

# Cannot compare incompatible units
try:
    result = (10 * ureg.meter) > (5 * ureg.kg)
except pint.errors.DimensionalityError as e:
    print(f"Error: Cannot compare length and mass")
```

### Best Practices

1. **Always specify units explicitly:**
   ```python
   # Good
   pressure = 150 * ureg.psi

   # Bad - units unclear
   pressure = 150
   ```

2. **Use base units for calculations:**
   ```python
   # Convert to base units for complex calculations
   result = calculation().to_base_units()
   ```

3. **Check dimensionality:**
   ```python
   # Verify result is what you expect
   if result.dimensionality == ureg.pascal.dimensionality:
       print("Result is a pressure")
   ```

4. **Handle unit strings carefully:**
   ```python
   # Parse unit strings safely
   try:
       quantity = ureg.Quantity(value, unit_string)
   except pint.errors.UndefinedUnitError:
       print(f"Unknown unit: {unit_string}")
   ```

## Using the converter.py Module

The included `converter.py` module provides convenient functions for common conversions:

```python
from converter import (
    convert_flow,
    convert_pressure,
    convert_viscosity_dynamic,
    convert_length,
    convert_power,
    convert_temperature,
    check_dimensional_consistency,
    ureg
)

# Simple conversions
flow_gpm = convert_flow(100, 'L/min', 'gal/min')
print(f"Flow: {flow_gpm:.2f} gpm")

pressure_bar = convert_pressure(150, 'psi', 'bar')
print(f"Pressure: {pressure_bar:.2f} bar")

# Check if operation is dimensionally valid
is_valid = check_dimensional_consistency('m/s', 'ft/min')
print(f"Can convert velocity: {is_valid}")  # True

is_valid = check_dimensional_consistency('m', 'kg')
print(f"Can convert length to mass: {is_valid}")  # False
```

## Quick Reference

See the included quick reference table for common conversion factors and formulas.

## Error Handling

Always wrap conversions in try-except blocks for production code:

```python
try:
    result = convert_pressure(150, 'psi', 'bar')
    print(f"Pressure: {result:.2f} bar")
except pint.errors.UndefinedUnitError:
    print("Unknown unit specified")
except pint.errors.DimensionalityError:
    print("Units are not compatible")
except Exception as e:
    print(f"Conversion error: {e}")
```

## Tips for Engineering Applications

1. **Always document units in comments and variable names**
2. **Use consistent unit systems within calculations**
3. **Verify dimensional analysis for derived quantities**
4. **Store values with units throughout calculations**
5. **Only convert to display units at the end**
6. **Use Pint's formatting for clear output**

## Resources

- [Pint Documentation](https://pint.readthedocs.io/)
- [Unit Registry](https://github.com/hgrecco/pint/blob/master/pint/default_en.txt)
- [NIST Guide to SI Units](https://www.nist.gov/pml/weights-and-measures/metric-si)
