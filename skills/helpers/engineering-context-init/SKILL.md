---
name: engineering-context-init
description: "Initialize engineering session with standard constants, units, and imports"
category: helpers
domain: general
complexity: basic
dependencies:
  - numpy
  - pint
  - scipy
  - matplotlib
---

# Engineering Context Initialization Helper

## Purpose

Standardize engineering sessions with a consistent environment that includes:
- Standard physical constants
- Unit registry for dimensional analysis
- Common scientific computing libraries
- Default plotting configurations
- Precision settings for numerical calculations

This helper eliminates the need to repeatedly set up the same environment at the start of each engineering session, ensuring consistency and saving time.

## What It Sets Up

### Standard Constants
- **g**: Standard gravity (9.80665 m/s²)
- **R**: Universal gas constant (8.314462618 J/(mol·K))
- **R_air**: Specific gas constant for air (287.05 J/(kg·K))
- **gamma_air**: Specific heat ratio for air (1.4)
- **P_atm**: Standard atmospheric pressure (101325 Pa)
- **T_std**: Standard temperature (288.15 K / 15°C)
- **rho_air_std**: Standard air density (1.225 kg/m³)
- **mu_air_std**: Dynamic viscosity of air at standard conditions (1.789e-5 Pa·s)
- **nu_air_std**: Kinematic viscosity of air at standard conditions (1.460e-5 m²/s)
- **c_light**: Speed of light (299792458 m/s)
- **sigma_sb**: Stefan-Boltzmann constant (5.670374419e-8 W/(m²·K⁴))

### Unit Registry (Pint)
- Configured pint UnitRegistry for dimensional analysis
- Support for unit conversions
- Unit-aware calculations

### Common Imports
- **NumPy**: Array operations and mathematical functions
- **SciPy**: Scientific computing (optimize, integrate, interpolate, special functions)
- **Matplotlib**: Plotting and visualization
- **Pandas**: Data manipulation (optional, if available)

### Default Plotting Style
- Figure size: 10x6 inches
- Grid enabled
- LaTeX-style font rendering (if available)
- Professional color schemes

### Precision Settings
- NumPy print precision: 6 decimal places
- Print suppression for small values
- Consistent formatting across sessions

## Usage

### Quick Start

Simply run the initialization script at the start of your session:

```python
exec(open('init-script.py').read())
```

Or copy-paste the contents of `init-script.py` directly into your Python session.

### Example 1: Basic Calculations with Units

```python
# After running init-script.py
from pint import UnitRegistry
ureg = UnitRegistry()

# Define a pressure with units
pressure = 500 * ureg.kPa
print(f"Pressure: {pressure}")
print(f"Pressure in psi: {pressure.to('psi')}")

# Calculate ideal gas properties
T = 300 * ureg.kelvin
V = 2 * ureg.m**3
n = (pressure * V / (R * T)).to('mol')
print(f"Moles of gas: {n:.2f}")
```

### Example 2: Atmospheric Calculations

```python
# After running init-script.py

# Calculate density at altitude using standard constants
altitude = 5000  # meters
T = T_std - 0.0065 * altitude  # Temperature lapse rate
P = P_atm * (T / T_std) ** (g / (0.0065 * R_air))
rho = P / (R_air * T)

print(f"At {altitude}m altitude:")
print(f"Temperature: {T:.2f} K ({T-273.15:.2f} °C)")
print(f"Pressure: {P:.2f} Pa ({P/1000:.2f} kPa)")
print(f"Density: {rho:.4f} kg/m³")
```

### Example 3: Flow Calculations

```python
# After running init-script.py

# Reynolds number calculation
velocity = 50  # m/s
length = 2     # m characteristic length
Re = velocity * length / nu_air_std

print(f"Reynolds number: {Re:.2e}")
if Re > 5e5:
    print("Flow is turbulent")
else:
    print("Flow is laminar")
```

### Example 4: Quick Plotting

```python
# After running init-script.py
import numpy as np
import matplotlib.pyplot as plt

# The plotting style is already configured
x = np.linspace(0, 2*np.pi, 100)
y = np.sin(x)

plt.figure()
plt.plot(x, y, linewidth=2)
plt.xlabel('x [rad]')
plt.ylabel('sin(x)')
plt.title('Example Plot with Pre-configured Style')
plt.show()
```

### Example 5: Unit Conversions

```python
# After running init-script.py

# Power calculation with unit conversions
force = 1000 * ureg.N
velocity = 30 * ureg.m / ureg.s
power = (force * velocity).to('kW')

print(f"Power: {power}")
print(f"Power in hp: {power.to('hp'):.2f}")
print(f"Power in BTU/hr: {power.to('BTU/hr'):.2f}")
```

## Files Included

- **SKILL.md**: This documentation file
- **init-script.py**: Complete initialization script ready to run
- **constants.md**: Reference list of all constants and their values

## Notes

- The script checks for library availability and imports only what's available
- If matplotlib is not available, plotting configuration is skipped
- If pandas is not available, it's simply not imported (optional)
- All constants use SI units unless otherwise specified
- The pint UnitRegistry is created as `ureg` for unit-aware calculations

## Customization

You can modify `init-script.py` to:
- Add domain-specific constants
- Change default plotting styles
- Adjust precision settings
- Include additional library imports
- Set up custom unit definitions in pint

## Best Practices

1. Run the initialization at the start of each new engineering session
2. Keep the script updated with your most commonly used constants
3. Document any custom constants you add
4. Use the unit registry for all calculations involving physical quantities
5. Store the script in your project directory for easy access
