# Physical Constants Reference

This document lists all physical constants defined in the engineering context initialization script.

## Fundamental Constants

| Constant | Symbol | Value | Units | Description |
|----------|--------|-------|-------|-------------|
| Speed of light | `c_light` | 299,792,458 | m/s | Speed of light in vacuum |
| Planck constant | `h_planck` | 6.62607015 × 10⁻³⁴ | J·s | Planck constant |
| Boltzmann constant | `k_B` | 1.380649 × 10⁻²³ | J/K | Boltzmann constant |
| Avogadro's number | `N_A` | 6.02214076 × 10²³ | 1/mol | Avogadro constant |
| Pi | `pi` | 3.141592653589793 | - | Mathematical constant π |

## Thermodynamic Constants

| Constant | Symbol | Value | Units | Description |
|----------|--------|-------|-------|-------------|
| Universal gas constant | `R` | 8.314462618 | J/(mol·K) | Universal gas constant |
| Stefan-Boltzmann constant | `sigma_sb` | 5.670374419 × 10⁻⁸ | W/(m²·K⁴) | Stefan-Boltzmann constant |

## Gravitational Constants

| Constant | Symbol | Value | Units | Description |
|----------|--------|-------|-------|-------------|
| Standard gravity | `g` | 9.80665 | m/s² | Standard gravitational acceleration |

## Air Properties at Standard Conditions

### Gas Constants for Air

| Constant | Symbol | Value | Units | Description |
|----------|--------|-------|-------|-------------|
| Specific gas constant (air) | `R_air` | 287.05 | J/(kg·K) | Specific gas constant for dry air |
| Specific heat ratio (air) | `gamma_air` | 1.4 | - | Ratio of specific heats (cp/cv) for air |
| Specific heat at constant pressure | `cp_air` | 1005.0 | J/(kg·K) | cp for air at standard conditions |
| Specific heat at constant volume | `cv_air` | 718.0 | J/(kg·K) | cv for air at standard conditions |

### Standard Atmospheric Conditions

| Constant | Symbol | Value | Units | Description |
|----------|--------|-------|-------|-------------|
| Standard atmospheric pressure | `P_atm` | 101,325 | Pa | Sea level standard atmospheric pressure |
| Standard temperature | `T_std` | 288.15 | K | Standard temperature (15°C) |
| Standard air density | `rho_air_std` | 1.225 | kg/m³ | Air density at standard conditions |

### Air Viscosity at Standard Conditions

| Constant | Symbol | Value | Units | Description |
|----------|--------|-------|-------|-------------|
| Dynamic viscosity (air) | `mu_air_std` | 1.789 × 10⁻⁵ | Pa·s | Dynamic viscosity at standard conditions |
| Kinematic viscosity (air) | `nu_air_std` | 1.460 × 10⁻⁵ | m²/s | Kinematic viscosity at standard conditions |

## Quick Reference Table (Copy-Paste Ready)

```python
# Fundamental Constants
c_light = 299792458.0          # m/s - Speed of light
h_planck = 6.62607015e-34      # J·s - Planck constant
k_B = 1.380649e-23             # J/K - Boltzmann constant
N_A = 6.02214076e23            # 1/mol - Avogadro's number
pi = 3.141592653589793         # Mathematical constant

# Gravity
g = 9.80665                    # m/s² - Standard gravity

# Thermodynamic
R = 8.314462618                # J/(mol·K) - Universal gas constant
sigma_sb = 5.670374419e-8      # W/(m²·K⁴) - Stefan-Boltzmann constant

# Air Properties
R_air = 287.05                 # J/(kg·K) - Specific gas constant for air
gamma_air = 1.4                # Specific heat ratio for air
cp_air = 1005.0                # J/(kg·K) - Specific heat at constant pressure
cv_air = 718.0                 # J/(kg·K) - Specific heat at constant volume

# Standard Atmosphere
P_atm = 101325.0               # Pa - Standard atmospheric pressure
T_std = 288.15                 # K - Standard temperature (15°C)
rho_air_std = 1.225            # kg/m³ - Standard air density

# Air Viscosity (Standard Conditions)
mu_air_std = 1.789e-5          # Pa·s - Dynamic viscosity
nu_air_std = 1.460e-5          # m²/s - Kinematic viscosity
```

## Usage Examples

### Converting Units

```python
# Using pint unit registry
pressure = 500 * ureg.kPa
print(pressure.to('psi'))

temperature = 300 * ureg.kelvin
print(temperature.to('degC'))
```

### Atmospheric Calculations

```python
# Calculate properties at 5000m altitude
props = atm_properties(5000)
print(f"Temperature: {props['Temperature_C']:.2f} °C")
print(f"Pressure: {props['Pressure_kPa']:.2f} kPa")
print(f"Density: {props['Density_kg_m3']:.4f} kg/m³")
```

### Reynolds Number

```python
# Flow over a wing
velocity = 50  # m/s
chord = 2      # m
Re = reynolds_number(velocity, chord)
print(f"Reynolds number: {Re:.2e}")
```

### Mach Number

```python
# Aircraft speed
velocity = 250  # m/s
altitude = 10000  # m
props = atm_properties(altitude)
M = mach_number(velocity, props['Temperature_K'])
print(f"Mach number: {M:.3f}")
```

## Notes

- All constants use SI units unless otherwise specified
- Standard atmospheric conditions are defined at sea level and 15°C (288.15 K)
- Air properties assume dry air composition
- Constants are based on CODATA 2018 and ISO 2533:1975 standards

## References

1. CODATA Recommended Values of the Fundamental Physical Constants: 2018
2. ISO 2533:1975 - Standard Atmosphere
3. NIST Reference on Constants, Units, and Uncertainty
4. Engineering Toolbox - Air Properties

## Related Constants Not Included

If you need additional constants, consider adding:

- **Water properties**: Density, viscosity, specific heat
- **Material properties**: Young's modulus, thermal conductivity, etc.
- **Electromagnetic constants**: Permittivity, permeability
- **Astronomical constants**: Earth radius, solar constant, etc.

These can be easily added to the `init-script.py` file as needed for your specific engineering domain.
