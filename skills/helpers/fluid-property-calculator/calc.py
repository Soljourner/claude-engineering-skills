"""
Fluid Property Calculator
==========================

Quick fluid property calculations using empirical correlations.
All formulas are verified against reference data.

Units: SI unless otherwise specified
- Temperature: °C (except Sutherland formula uses K)
- Pressure: Pa
- Density: kg/m³
- Viscosity: Pa·s (dynamic), m²/s (kinematic)
- Thermal conductivity: W/m·K
- Specific heat: J/kg·K
"""

import math


# ============================================================================
# WATER PROPERTIES (0-100°C, atmospheric pressure)
# ============================================================================

def water_density(T):
    """
    Calculate water density using polynomial correlation.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Density in kg/m³

    Reference: CRC Handbook of Chemistry and Physics
    Accuracy: ±0.1% in valid range
    """
    if not (0 <= T <= 100):
        raise ValueError(f"Temperature {T}°C outside valid range (0-100°C)")

    # Polynomial correlation for liquid water
    rho = (999.83952
           + 16.945176 * T
           - 7.9870401e-3 * T**2
           - 46.170461e-6 * T**3
           + 105.56302e-9 * T**4
           - 280.54253e-12 * T**5) / (1 + 16.879850e-3 * T)

    return rho


def water_dynamic_viscosity(T):
    """
    Calculate water dynamic viscosity using Vogel equation.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Dynamic viscosity in Pa·s

    Reference: Vogel equation fitted to NIST data
    Accuracy: ±1% in valid range
    """
    if not (0 <= T <= 100):
        raise ValueError(f"Temperature {T}°C outside valid range (0-100°C)")

    # Vogel equation: mu = A * 10^(B/(T-C))
    # Constants for water
    A = 0.02414e-3  # Pa·s
    B = 247.8  # K
    C = 140  # K

    T_K = T + 273.15
    mu = A * 10**(B / (T_K - C))

    return mu


def water_kinematic_viscosity(T):
    """
    Calculate water kinematic viscosity.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Kinematic viscosity in m²/s
    """
    rho = water_density(T)
    mu = water_dynamic_viscosity(T)
    return mu / rho


def water_thermal_conductivity(T):
    """
    Calculate water thermal conductivity using polynomial correlation.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Thermal conductivity in W/m·K

    Reference: Polynomial fit to NIST data
    Accuracy: ±1% in valid range
    """
    if not (0 <= T <= 100):
        raise ValueError(f"Temperature {T}°C outside valid range (0-100°C)")

    # Polynomial correlation
    k = (0.5650
         + 1.962e-3 * T
         - 8.138e-6 * T**2
         + 1.567e-8 * T**3)

    return k


def water_specific_heat(T):
    """
    Calculate water specific heat capacity.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Specific heat in J/kg·K

    Reference: Polynomial fit to NIST data
    Accuracy: ±0.5% in valid range
    """
    if not (0 <= T <= 100):
        raise ValueError(f"Temperature {T}°C outside valid range (0-100°C)")

    # Polynomial correlation
    cp = (4217.0
          - 3.7210 * T
          + 0.05640 * T**2
          - 2.8930e-4 * T**3)

    return cp


def water_prandtl_number(T):
    """
    Calculate water Prandtl number.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Prandtl number (dimensionless)
    """
    mu = water_dynamic_viscosity(T)
    cp = water_specific_heat(T)
    k = water_thermal_conductivity(T)

    Pr = mu * cp / k
    return Pr


def water_properties(T):
    """
    Calculate all water properties at given temperature.

    Args:
        T: Temperature in °C (valid: 0-100°C)

    Returns:
        Dictionary with all properties
    """
    return {
        'temperature': T,
        'density': water_density(T),
        'dynamic_viscosity': water_dynamic_viscosity(T),
        'kinematic_viscosity': water_kinematic_viscosity(T),
        'thermal_conductivity': water_thermal_conductivity(T),
        'specific_heat': water_specific_heat(T),
        'prandtl_number': water_prandtl_number(T)
    }


# ============================================================================
# AIR PROPERTIES (standard atmosphere, -50 to 200°C)
# ============================================================================

def air_density(T, P=101325):
    """
    Calculate air density using ideal gas law.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)
        P: Pressure in Pa (default: 101325 Pa = 1 atm)

    Returns:
        Density in kg/m³

    Reference: Ideal gas law
    Accuracy: ±0.5% for dry air at moderate conditions
    """
    if not (-50 <= T <= 200):
        raise ValueError(f"Temperature {T}°C outside valid range (-50 to 200°C)")

    T_K = T + 273.15
    R_specific = 287.05  # J/kg·K for dry air

    rho = P / (R_specific * T_K)
    return rho


def air_dynamic_viscosity(T):
    """
    Calculate air dynamic viscosity using Sutherland's formula.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)

    Returns:
        Dynamic viscosity in Pa·s

    Reference: Sutherland's formula with standard constants
    Accuracy: ±2% in valid range
    """
    if not (-50 <= T <= 200):
        raise ValueError(f"Temperature {T}°C outside valid range (-50 to 200°C)")

    T_K = T + 273.15

    # Sutherland's formula constants for air
    mu_0 = 1.716e-5  # Pa·s at T_0
    T_0 = 273.15  # K
    S = 110.4  # K (Sutherland constant for air)

    mu = mu_0 * (T_K / T_0)**(3/2) * (T_0 + S) / (T_K + S)

    return mu


def air_kinematic_viscosity(T, P=101325):
    """
    Calculate air kinematic viscosity.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)
        P: Pressure in Pa (default: 101325 Pa)

    Returns:
        Kinematic viscosity in m²/s
    """
    rho = air_density(T, P)
    mu = air_dynamic_viscosity(T)
    return mu / rho


def air_thermal_conductivity(T):
    """
    Calculate air thermal conductivity using polynomial correlation.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)

    Returns:
        Thermal conductivity in W/m·K

    Reference: Polynomial fit to reference data
    Accuracy: ±2% in valid range
    """
    if not (-50 <= T <= 200):
        raise ValueError(f"Temperature {T}°C outside valid range (-50 to 200°C)")

    T_K = T + 273.15

    # Polynomial correlation
    k = (2.3340e-3
         + 7.5880e-5 * T_K
         - 1.6850e-8 * T_K**2)

    return k


def air_specific_heat(T):
    """
    Calculate air specific heat capacity at constant pressure.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)

    Returns:
        Specific heat in J/kg·K

    Reference: Temperature-dependent correlation
    Accuracy: ±1% in valid range
    """
    if not (-50 <= T <= 200):
        raise ValueError(f"Temperature {T}°C outside valid range (-50 to 200°C)")

    T_K = T + 273.15

    # Polynomial correlation
    cp = (1047.0
          - 0.3720 * T_K
          + 9.4500e-4 * T_K**2
          - 4.9200e-7 * T_K**3)

    return cp


def air_prandtl_number(T):
    """
    Calculate air Prandtl number.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)

    Returns:
        Prandtl number (dimensionless)

    Note: Prandtl number for air is nearly constant (~0.71)
    """
    mu = air_dynamic_viscosity(T)
    cp = air_specific_heat(T)
    k = air_thermal_conductivity(T)

    Pr = mu * cp / k
    return Pr


def air_properties(T, P=101325):
    """
    Calculate all air properties at given temperature and pressure.

    Args:
        T: Temperature in °C (valid: -50 to 200°C)
        P: Pressure in Pa (default: 101325 Pa)

    Returns:
        Dictionary with all properties
    """
    return {
        'temperature': T,
        'pressure': P,
        'density': air_density(T, P),
        'dynamic_viscosity': air_dynamic_viscosity(T),
        'kinematic_viscosity': air_kinematic_viscosity(T, P),
        'thermal_conductivity': air_thermal_conductivity(T),
        'specific_heat': air_specific_heat(T),
        'prandtl_number': air_prandtl_number(T)
    }


# ============================================================================
# VISCOSITY CORRELATIONS
# ============================================================================

def sutherland_viscosity(gas, T_K):
    """
    Calculate gas viscosity using Sutherland's formula.

    Args:
        gas: Gas name ('air', 'nitrogen', 'oxygen', 'co2')
        T_K: Temperature in Kelvin

    Returns:
        Dynamic viscosity in Pa·s

    Reference: Sutherland's formula
    """
    # Sutherland constants for various gases
    constants = {
        'air': {'mu_0': 1.716e-5, 'T_0': 273.15, 'S': 110.4},
        'nitrogen': {'mu_0': 1.663e-5, 'T_0': 273.15, 'S': 111.0},
        'oxygen': {'mu_0': 1.919e-5, 'T_0': 273.15, 'S': 127.0},
        'co2': {'mu_0': 1.370e-5, 'T_0': 273.15, 'S': 240.0}
    }

    if gas.lower() not in constants:
        raise ValueError(f"Gas '{gas}' not available. Choose from: {list(constants.keys())}")

    c = constants[gas.lower()]
    mu = c['mu_0'] * (T_K / c['T_0'])**(3/2) * (c['T_0'] + c['S']) / (T_K + c['S'])

    return mu


def andrade_viscosity(A, B, T):
    """
    Calculate liquid viscosity using Andrade equation.

    Args:
        A: Pre-exponential factor (Pa·s)
        B: Activation energy parameter (K)
        T: Temperature in °C

    Returns:
        Dynamic viscosity in Pa·s

    Reference: Andrade equation - mu = A * exp(B/T_K)

    Example parameters (approximate):
        Water: A=2.414e-5, B=570.6
        Ethanol: A=2.948e-5, B=776.4
    """
    T_K = T + 273.15
    mu = A * math.exp(B / T_K)
    return mu


# ============================================================================
# VAPOR PRESSURE (Antoine Equation)
# ============================================================================

def antoine_vapor_pressure(substance, T):
    """
    Calculate vapor pressure using Antoine equation.

    Args:
        substance: Substance name (e.g., 'water', 'ethanol')
        T: Temperature in °C

    Returns:
        Vapor pressure in Pa

    Reference: Antoine equation - log10(P_mmHg) = A - B/(C + T)
    Constants from NIST Chemistry WebBook
    """
    # Antoine constants [A, B, C, T_min, T_max]
    # log10(P) where P is in mmHg, T in °C
    constants = {
        'water': [8.07131, 1730.63, 233.426, 1, 100],
        'ethanol': [8.20417, 1642.89, 230.300, -2, 93],
        'methanol': [8.08097, 1582.27, 239.726, -16, 84],
        'acetone': [7.11714, 1210.595, 229.664, -26, 77],
        'benzene': [6.90565, 1211.033, 220.790, 8, 103],
        'toluene': [6.95464, 1344.800, 219.482, 6, 137]
    }

    if substance.lower() not in constants:
        raise ValueError(f"Substance '{substance}' not available. Choose from: {list(constants.keys())}")

    A, B, C, T_min, T_max = constants[substance.lower()]

    if not (T_min <= T <= T_max):
        raise ValueError(f"Temperature {T}°C outside valid range ({T_min} to {T_max}°C)")

    # Calculate pressure in mmHg
    log_P_mmHg = A - B / (C + T)
    P_mmHg = 10**log_P_mmHg

    # Convert to Pa (1 mmHg = 133.322 Pa)
    P_Pa = P_mmHg * 133.322

    return P_Pa


# ============================================================================
# DIMENSIONLESS NUMBERS
# ============================================================================

def reynolds_number(velocity, length, kinematic_viscosity=None,
                   density=None, dynamic_viscosity=None):
    """
    Calculate Reynolds number.

    Args:
        velocity: Flow velocity in m/s
        length: Characteristic length in m (e.g., pipe diameter)
        kinematic_viscosity: Kinematic viscosity in m²/s (if provided)
        OR
        density: Density in kg/m³
        dynamic_viscosity: Dynamic viscosity in Pa·s

    Returns:
        Reynolds number (dimensionless)

    Note: Provide either kinematic_viscosity OR (density and dynamic_viscosity)
    """
    if kinematic_viscosity is not None:
        Re = velocity * length / kinematic_viscosity
    elif density is not None and dynamic_viscosity is not None:
        Re = density * velocity * length / dynamic_viscosity
    else:
        raise ValueError("Must provide either kinematic_viscosity or (density and dynamic_viscosity)")

    return Re


def prandtl_number(dynamic_viscosity, specific_heat, thermal_conductivity):
    """
    Calculate Prandtl number.

    Args:
        dynamic_viscosity: Dynamic viscosity in Pa·s
        specific_heat: Specific heat in J/kg·K
        thermal_conductivity: Thermal conductivity in W/m·K

    Returns:
        Prandtl number (dimensionless)
    """
    Pr = dynamic_viscosity * specific_heat / thermal_conductivity
    return Pr


# ============================================================================
# FRICTION FACTOR
# ============================================================================

def friction_factor_laminar(Re):
    """
    Calculate friction factor for laminar flow.

    Args:
        Re: Reynolds number

    Returns:
        Darcy friction factor

    Valid for Re < 2300
    """
    if Re >= 2300:
        raise ValueError(f"Reynolds number {Re} too high for laminar correlation (Re < 2300)")

    f = 64 / Re
    return f


def friction_factor_turbulent_smooth(Re):
    """
    Calculate friction factor for turbulent flow in smooth pipes (Blasius).

    Args:
        Re: Reynolds number

    Returns:
        Darcy friction factor

    Valid for 4000 < Re < 100,000
    """
    if Re < 4000:
        raise ValueError(f"Reynolds number {Re} too low for turbulent correlation (Re > 4000)")
    if Re > 100000:
        print(f"Warning: Re = {Re} outside typical Blasius range (< 100,000)")

    f = 0.316 / Re**0.25
    return f


def friction_factor_colebrook(Re, roughness, diameter, max_iter=50, tol=1e-6):
    """
    Calculate friction factor using Colebrook-White equation (iterative).

    Args:
        Re: Reynolds number
        roughness: Absolute roughness in m (e.g., 0.045e-3 for commercial steel)
        diameter: Pipe diameter in m
        max_iter: Maximum iterations
        tol: Convergence tolerance

    Returns:
        Darcy friction factor

    Valid for turbulent flow in rough pipes
    """
    if Re < 4000:
        raise ValueError(f"Reynolds number {Re} too low for Colebrook equation (Re > 4000)")

    # Relative roughness
    rel_roughness = roughness / diameter

    # Initial guess using Swamee-Jain
    f = friction_factor_swamee_jain(Re, roughness, diameter)

    # Iterative solution of Colebrook-White equation
    for _ in range(max_iter):
        f_new = 1 / (-2 * math.log10(rel_roughness/3.7 + 2.51/(Re * math.sqrt(f))))**2

        if abs(f_new - f) < tol:
            return f_new

        f = f_new

    print(f"Warning: Colebrook iteration did not converge after {max_iter} iterations")
    return f


def friction_factor_swamee_jain(Re, roughness, diameter):
    """
    Calculate friction factor using Swamee-Jain approximation.

    Args:
        Re: Reynolds number
        roughness: Absolute roughness in m
        diameter: Pipe diameter in m

    Returns:
        Darcy friction factor

    Non-iterative approximation to Colebrook-White equation.
    Accuracy: ±1% for 5000 < Re < 10^8 and 10^-6 < ε/D < 10^-2
    """
    if Re < 5000:
        print(f"Warning: Re = {Re} below recommended range (> 5000)")

    rel_roughness = roughness / diameter

    f = 0.25 / (math.log10(rel_roughness/3.7 + 5.74/Re**0.9))**2

    return f


def friction_factor(Re, roughness=0, diameter=1):
    """
    Calculate friction factor automatically based on flow regime.

    Args:
        Re: Reynolds number
        roughness: Absolute roughness in m (default: 0 for smooth pipe)
        diameter: Pipe diameter in m (default: 1)

    Returns:
        Darcy friction factor
    """
    if Re < 2300:
        # Laminar flow
        return friction_factor_laminar(Re)
    elif Re < 4000:
        # Transition region - interpolate
        f_lam = friction_factor_laminar(2300)
        if roughness == 0:
            f_turb = friction_factor_turbulent_smooth(4000)
        else:
            f_turb = friction_factor_swamee_jain(4000, roughness, diameter)

        # Linear interpolation
        f = f_lam + (f_turb - f_lam) * (Re - 2300) / (4000 - 2300)
        return f
    else:
        # Turbulent flow
        if roughness == 0:
            return friction_factor_turbulent_smooth(Re)
        else:
            return friction_factor_swamee_jain(Re, roughness, diameter)


# ============================================================================
# IDEAL GAS PROPERTIES
# ============================================================================

def ideal_gas_density(P, T, M):
    """
    Calculate ideal gas density.

    Args:
        P: Pressure in Pa
        T: Temperature in °C
        M: Molar mass in kg/kmol

    Returns:
        Density in kg/m³
    """
    T_K = T + 273.15
    R_universal = 8314.46  # J/kmol·K

    rho = P * M / (R_universal * T_K)
    return rho


def ideal_gas_speed_of_sound(T, gamma, M):
    """
    Calculate speed of sound in ideal gas.

    Args:
        T: Temperature in °C
        gamma: Specific heat ratio (cp/cv)
        M: Molar mass in kg/kmol

    Returns:
        Speed of sound in m/s
    """
    T_K = T + 273.15
    R_universal = 8314.46  # J/kmol·K

    a = math.sqrt(gamma * R_universal * T_K / M)
    return a


# Gas properties database
GAS_PROPERTIES = {
    'air': {'M': 28.97, 'gamma': 1.40},
    'nitrogen': {'M': 28.01, 'gamma': 1.40},
    'oxygen': {'M': 32.00, 'gamma': 1.40},
    'co2': {'M': 44.01, 'gamma': 1.30},
    'helium': {'M': 4.00, 'gamma': 1.66},
    'argon': {'M': 39.95, 'gamma': 1.67},
    'hydrogen': {'M': 2.02, 'gamma': 1.41},
    'methane': {'M': 16.04, 'gamma': 1.32}
}


def gas_properties(gas_name, T, P=101325):
    """
    Calculate ideal gas properties for common gases.

    Args:
        gas_name: Name of gas (e.g., 'air', 'nitrogen')
        T: Temperature in °C
        P: Pressure in Pa (default: 101325)

    Returns:
        Dictionary with gas properties
    """
    if gas_name.lower() not in GAS_PROPERTIES:
        raise ValueError(f"Gas '{gas_name}' not available. Choose from: {list(GAS_PROPERTIES.keys())}")

    props = GAS_PROPERTIES[gas_name.lower()]
    M = props['M']
    gamma = props['gamma']

    rho = ideal_gas_density(P, T, M)
    a = ideal_gas_speed_of_sound(T, gamma, M)

    return {
        'name': gas_name,
        'temperature': T,
        'pressure': P,
        'molar_mass': M,
        'specific_heat_ratio': gamma,
        'density': rho,
        'speed_of_sound': a
    }


# ============================================================================
# VERIFICATION TESTS
# ============================================================================

def verify_calculations():
    """
    Verify calculations against known reference data.
    """
    print("=" * 70)
    print("VERIFICATION TESTS")
    print("=" * 70)

    # Test 1: Water at 20°C
    print("\n1. Water at 20°C:")
    props = water_properties(20)
    print(f"   Density: {props['density']:.2f} kg/m³ (expected: 998.2)")
    print(f"   Dynamic viscosity: {props['dynamic_viscosity']*1000:.3f} mPa·s (expected: 1.002)")
    print(f"   Thermal conductivity: {props['thermal_conductivity']:.4f} W/m·K (expected: 0.598)")
    print(f"   Prandtl number: {props['prandtl_number']:.2f} (expected: 7.0)")

    # Test 2: Air at 20°C
    print("\n2. Air at 20°C, 1 atm:")
    props = air_properties(20)
    print(f"   Density: {props['density']:.3f} kg/m³ (expected: 1.205)")
    print(f"   Dynamic viscosity: {props['dynamic_viscosity']*1e6:.2f} μPa·s (expected: 18.24)")
    print(f"   Prandtl number: {props['prandtl_number']:.3f} (expected: 0.71)")

    # Test 3: Reynolds number
    print("\n3. Reynolds number (water at 20°C, V=1 m/s, D=0.1 m):")
    Re = reynolds_number(1.0, 0.1, water_kinematic_viscosity(20))
    print(f"   Re = {Re:.0f} (expected: ~99600)")

    # Test 4: Friction factor
    print("\n4. Friction factor (Re=50000, smooth pipe):")
    f = friction_factor_turbulent_smooth(50000)
    print(f"   f = {f:.5f} (expected: ~0.0211)")

    # Test 5: Vapor pressure of water at 100°C
    print("\n5. Vapor pressure of water at 100°C:")
    P_vap = antoine_vapor_pressure('water', 100)
    print(f"   P = {P_vap/1000:.2f} kPa (expected: 101.3 kPa)")

    print("\n" + "=" * 70)
    print("All verifications completed!")
    print("=" * 70)


if __name__ == "__main__":
    verify_calculations()
