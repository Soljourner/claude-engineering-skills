"""
Unit Converter Module for Engineering Applications

This module provides convenient functions for converting between common
engineering units using the Pint library. It includes pre-configured
unit registry and validation functions.

Usage:
    from converter import convert_flow, convert_pressure, ureg

    flow_gpm = convert_flow(100, 'L/min', 'gal/min')
    pressure_bar = convert_pressure(150, 'psi', 'bar')
"""

import pint
from typing import Union, Tuple

# Initialize unit registry
ureg = pint.UnitRegistry()

# Enable automatic conversion in calculations
ureg.default_format = '.4g'


# ============================================================================
# FLOW RATE CONVERSIONS
# ============================================================================

def convert_flow(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert flow rates between different units.

    Supported units:
        - m^3/s, m^3/h, m^3/min
        - L/s, L/min, L/h
        - gal/s, gal/min, gal/h (US gallons)
        - ft^3/s, ft^3/min (cfm), ft^3/h

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'L/min', 'm^3/s')
        to_unit: Target unit (e.g., 'gal/min', 'cfm')

    Returns:
        Converted value as float

    Examples:
        >>> convert_flow(100, 'L/min', 'gal/min')
        26.417205235815

        >>> convert_flow(1000, 'cfm', 'm^3/s')
        0.47194744525

        >>> convert_flow(50, 'm^3/h', 'L/min')
        833.33333333
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_volumetric_flow(value: float, from_unit: str, to_unit: str) -> float:
    """Alias for convert_flow() for clarity."""
    return convert_flow(value, from_unit, to_unit)


# ============================================================================
# PRESSURE CONVERSIONS
# ============================================================================

def convert_pressure(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert pressure between different units.

    Supported units:
        - Pa, kPa, MPa
        - bar, mbar
        - psi, psf
        - atm
        - mmHg, inHg
        - torr

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'psi', 'bar')
        to_unit: Target unit (e.g., 'kPa', 'atm')

    Returns:
        Converted value as float

    Examples:
        >>> convert_pressure(150, 'psi', 'bar')
        10.342137959999

        >>> convert_pressure(2.5, 'bar', 'kPa')
        250.0

        >>> convert_pressure(760, 'mmHg', 'atm')
        0.9999999974
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_gauge_to_absolute(gauge_pressure: float, unit: str,
                               atmospheric_pressure: float = 101.325) -> float:
    """
    Convert gauge pressure to absolute pressure.

    Args:
        gauge_pressure: Gauge pressure value
        unit: Pressure unit (e.g., 'kPa', 'psi')
        atmospheric_pressure: Atmospheric pressure in kPa (default: 101.325 kPa)

    Returns:
        Absolute pressure in the same unit

    Example:
        >>> convert_gauge_to_absolute(100, 'psi')
        114.69594877551
    """
    # Convert atmospheric pressure to the same unit
    atm_in_unit = convert_pressure(atmospheric_pressure, 'kPa', unit)
    return gauge_pressure + atm_in_unit


def convert_absolute_to_gauge(absolute_pressure: float, unit: str,
                               atmospheric_pressure: float = 101.325) -> float:
    """
    Convert absolute pressure to gauge pressure.

    Args:
        absolute_pressure: Absolute pressure value
        unit: Pressure unit (e.g., 'kPa', 'psi')
        atmospheric_pressure: Atmospheric pressure in kPa (default: 101.325 kPa)

    Returns:
        Gauge pressure in the same unit

    Example:
        >>> convert_absolute_to_gauge(114.7, 'psi')
        0.00405122449
    """
    atm_in_unit = convert_pressure(atmospheric_pressure, 'kPa', unit)
    return absolute_pressure - atm_in_unit


# ============================================================================
# VISCOSITY CONVERSIONS
# ============================================================================

def convert_viscosity_dynamic(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert dynamic viscosity between different units.

    Supported units:
        - Pa*s (Pascal-second)
        - cP, centipoise
        - P, poise
        - mPa*s (millipascal-second, equivalent to cP)

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'cP', 'Pa*s')
        to_unit: Target unit (e.g., 'Pa*s', 'poise')

    Returns:
        Converted value as float

    Examples:
        >>> convert_viscosity_dynamic(1.002, 'cP', 'Pa*s')
        0.001002

        >>> convert_viscosity_dynamic(0.25, 'Pa*s', 'cP')
        250.0
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_viscosity_kinematic(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert kinematic viscosity between different units.

    Supported units:
        - m^2/s
        - cSt, centistokes
        - St, stokes
        - mm^2/s (equivalent to cSt)

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'cSt', 'm^2/s')
        to_unit: Target unit (e.g., 'm^2/s', 'stokes')

    Returns:
        Converted value as float

    Examples:
        >>> convert_viscosity_kinematic(1.004, 'cSt', 'm^2/s')
        1.004e-06

        >>> convert_viscosity_kinematic(100, 'cSt', 'stokes')
        1.0
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def dynamic_to_kinematic(dynamic_viscosity: float, density: float,
                         visc_unit: str = 'Pa*s', dens_unit: str = 'kg/m^3',
                         result_unit: str = 'm^2/s') -> float:
    """
    Convert dynamic viscosity to kinematic viscosity.

    Formula: ν = μ / ρ

    Args:
        dynamic_viscosity: Dynamic viscosity value
        density: Fluid density value
        visc_unit: Dynamic viscosity unit (default: 'Pa*s')
        dens_unit: Density unit (default: 'kg/m^3')
        result_unit: Desired kinematic viscosity unit (default: 'm^2/s')

    Returns:
        Kinematic viscosity in specified unit

    Example:
        >>> dynamic_to_kinematic(1.002e-3, 1000, 'Pa*s', 'kg/m^3', 'cSt')
        1.002
    """
    mu = dynamic_viscosity * ureg(visc_unit)
    rho = density * ureg(dens_unit)
    nu = (mu / rho).to(result_unit)
    return nu.magnitude


def kinematic_to_dynamic(kinematic_viscosity: float, density: float,
                         visc_unit: str = 'm^2/s', dens_unit: str = 'kg/m^3',
                         result_unit: str = 'Pa*s') -> float:
    """
    Convert kinematic viscosity to dynamic viscosity.

    Formula: μ = ν × ρ

    Args:
        kinematic_viscosity: Kinematic viscosity value
        density: Fluid density value
        visc_unit: Kinematic viscosity unit (default: 'm^2/s')
        dens_unit: Density unit (default: 'kg/m^3')
        result_unit: Desired dynamic viscosity unit (default: 'Pa*s')

    Returns:
        Dynamic viscosity in specified unit

    Example:
        >>> kinematic_to_dynamic(1.004, 1000, 'cSt', 'kg/m^3', 'cP')
        1004.0
    """
    nu = kinematic_viscosity * ureg(visc_unit)
    rho = density * ureg(dens_unit)
    mu = (nu * rho).to(result_unit)
    return mu.magnitude


# ============================================================================
# LENGTH CONVERSIONS
# ============================================================================

def convert_length(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert length between different units.

    Supported units:
        - m, cm, mm, km
        - ft, inch, yd, mile
        - nm (nanometer), um (micrometer)

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'inch', 'm')
        to_unit: Target unit (e.g., 'mm', 'ft')

    Returns:
        Converted value as float

    Examples:
        >>> convert_length(4, 'inch', 'mm')
        101.6

        >>> convert_length(15.5, 'm', 'ft')
        50.85301837

        >>> convert_length(0.5, 'mm', 'inch')
        0.019685039
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_area(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert area between different units.

    Supported units:
        - m^2, cm^2, mm^2, km^2
        - ft^2, inch^2, yd^2, mile^2
        - hectare, acre

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'ft^2', 'm^2')
        to_unit: Target unit (e.g., 'cm^2', 'inch^2')

    Returns:
        Converted value as float

    Examples:
        >>> convert_area(500, 'm^2', 'ft^2')
        5381.955208

        >>> convert_area(3.14, 'inch^2', 'cm^2')
        20.258064
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_volume(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert volume between different units.

    Supported units:
        - m^3, cm^3, mm^3
        - L, mL
        - ft^3, inch^3
        - gal, qt, pt (US liquid)

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'L', 'gal')
        to_unit: Target unit (e.g., 'm^3', 'ft^3')

    Returns:
        Converted value as float

    Examples:
        >>> convert_volume(5000, 'L', 'm^3')
        5.0

        >>> convert_volume(5000, 'L', 'gal')
        1320.860129
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


# ============================================================================
# POWER CONVERSIONS
# ============================================================================

def convert_power(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert power between different units.

    Supported units:
        - W, kW, MW
        - hp (horsepower - mechanical)
        - BTU/h, BTU/s

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'hp', 'kW')
        to_unit: Target unit (e.g., 'W', 'BTU/h')

    Returns:
        Converted value as float

    Examples:
        >>> convert_power(50, 'hp', 'kW')
        37.28499962

        >>> convert_power(15, 'kW', 'hp')
        20.11842799

        >>> convert_power(2.5, 'MW', 'hp')
        3352.140199
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_energy(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert energy between different units.

    Supported units:
        - J, kJ, MJ
        - kWh, Wh
        - BTU
        - cal, kcal
        - eV

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'kWh', 'BTU')
        to_unit: Target unit (e.g., 'MJ', 'J')

    Returns:
        Converted value as float

    Examples:
        >>> convert_energy(1000, 'J', 'BTU')
        0.9478171203

        >>> convert_energy(100, 'kWh', 'MJ')
        360.0
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


# ============================================================================
# TEMPERATURE CONVERSIONS
# ============================================================================

def convert_temperature(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert temperature between different units.

    Supported units:
        - degC (Celsius)
        - degF (Fahrenheit)
        - kelvin or K (Kelvin)
        - degR (Rankine)

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'degC', 'degF')
        to_unit: Target unit (e.g., 'kelvin', 'degF')

    Returns:
        Converted value as float

    Examples:
        >>> convert_temperature(85, 'degC', 'degF')
        185.0

        >>> convert_temperature(77, 'kelvin', 'degC')
        -196.15

        >>> convert_temperature(72, 'degF', 'degC')
        22.222222
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def temperature_difference(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert temperature difference (not absolute temperature).

    This is important because temperature differences convert differently
    than absolute temperatures (no offset needed).

    Args:
        value: Temperature difference value
        from_unit: Source unit (e.g., 'delta_degC', 'delta_degF')
        to_unit: Target unit (e.g., 'delta_kelvin', 'delta_degF')

    Returns:
        Converted temperature difference

    Examples:
        >>> temperature_difference(10, 'delta_degC', 'delta_degF')
        18.0

        >>> temperature_difference(50, 'delta_degF', 'delta_degC')
        27.777778
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


# ============================================================================
# VELOCITY CONVERSIONS
# ============================================================================

def convert_velocity(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert velocity between different units.

    Supported units:
        - m/s, km/h
        - ft/s, ft/min, mph
        - knot

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'm/s', 'mph')
        to_unit: Target unit (e.g., 'ft/s', 'km/h')

    Returns:
        Converted value as float

    Examples:
        >>> convert_velocity(2.5, 'm/s', 'ft/s')
        8.202099738

        >>> convert_velocity(100, 'km/h', 'mph')
        62.13711922
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


# ============================================================================
# MASS AND MASS FLOW CONVERSIONS
# ============================================================================

def convert_mass(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert mass between different units.

    Supported units:
        - kg, g, mg, ton (metric ton)
        - lb, oz, ton (US ton)

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'kg', 'lb')
        to_unit: Target unit (e.g., 'g', 'oz')

    Returns:
        Converted value as float

    Examples:
        >>> convert_mass(100, 'kg', 'lb')
        220.4622622

        >>> convert_mass(1, 'ton', 'kg')
        1000.0
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_mass_flow(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert mass flow rate between different units.

    Supported units:
        - kg/s, kg/min, kg/h
        - g/s, g/min
        - lb/s, lb/min, lb/h
        - ton/h

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'kg/h', 'lb/min')
        to_unit: Target unit (e.g., 'kg/s', 'lb/h')

    Returns:
        Converted value as float

    Examples:
        >>> convert_mass_flow(100, 'kg/h', 'kg/s')
        0.02777777778

        >>> convert_mass_flow(100, 'kg/h', 'lb/min')
        3.666666667
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


def convert_density(value: float, from_unit: str, to_unit: str) -> float:
    """
    Convert density between different units.

    Supported units:
        - kg/m^3, g/cm^3, g/L
        - lb/ft^3, lb/gal

    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'kg/m^3', 'lb/ft^3')
        to_unit: Target unit (e.g., 'g/cm^3', 'lb/gal')

    Returns:
        Converted value as float

    Examples:
        >>> convert_density(1000, 'kg/m^3', 'g/cm^3')
        1.0

        >>> convert_density(1000, 'kg/m^3', 'lb/ft^3')
        62.42796058
    """
    quantity = value * ureg(from_unit)
    result = quantity.to(to_unit)
    return result.magnitude


# ============================================================================
# DIMENSIONAL ANALYSIS AND VALIDATION
# ============================================================================

def check_dimensional_consistency(unit1: str, unit2: str) -> bool:
    """
    Check if two units have the same dimensionality.

    Args:
        unit1: First unit string (e.g., 'm/s')
        unit2: Second unit string (e.g., 'ft/min')

    Returns:
        True if units are dimensionally consistent, False otherwise

    Examples:
        >>> check_dimensional_consistency('m/s', 'ft/min')
        True

        >>> check_dimensional_consistency('m', 'kg')
        False

        >>> check_dimensional_consistency('Pa', 'psi')
        True
    """
    try:
        q1 = ureg(unit1)
        q2 = ureg(unit2)
        return q1.dimensionality == q2.dimensionality
    except (pint.errors.UndefinedUnitError, AttributeError):
        return False


def get_dimensionality(unit: str) -> str:
    """
    Get the dimensionality of a unit.

    Args:
        unit: Unit string (e.g., 'Pa', 'm/s')

    Returns:
        String representation of dimensionality

    Examples:
        >>> get_dimensionality('Pa')
        '[mass] / [length] / [time] ** 2'

        >>> get_dimensionality('m/s')
        '[length] / [time]'
    """
    try:
        q = ureg(unit)
        return str(q.dimensionality)
    except pint.errors.UndefinedUnitError:
        return "Unknown unit"


def get_base_units(value: float, unit: str) -> Tuple[float, str]:
    """
    Convert a quantity to SI base units.

    Args:
        value: Numerical value
        unit: Unit string

    Returns:
        Tuple of (value in base units, base unit string)

    Examples:
        >>> get_base_units(100, 'psi')
        (689475.7293, 'kilogram / meter / second ** 2')

        >>> get_base_units(10, 'km/h')
        (2.7777777778, 'meter / second')
    """
    quantity = value * ureg(unit)
    base_quantity = quantity.to_base_units()
    return base_quantity.magnitude, str(base_quantity.units)


def is_dimensionless(unit: str) -> bool:
    """
    Check if a unit is dimensionless.

    Args:
        unit: Unit string

    Returns:
        True if dimensionless, False otherwise

    Examples:
        >>> is_dimensionless('')
        True

        >>> is_dimensionless('m')
        False
    """
    try:
        q = ureg(unit)
        return q.dimensionality == ureg.dimensionless.dimensionality
    except pint.errors.UndefinedUnitError:
        return False


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def format_quantity(value: float, unit: str, precision: int = 4) -> str:
    """
    Format a quantity with its unit for display.

    Args:
        value: Numerical value
        unit: Unit string
        precision: Number of significant figures (default: 4)

    Returns:
        Formatted string

    Example:
        >>> format_quantity(123.456789, 'm/s', 3)
        '123 m / s'
    """
    quantity = value * ureg(unit)
    return f"{quantity:.{precision}g~P}"


def convert_with_context(value: float, from_unit: str, to_unit: str,
                         context: str = None) -> float:
    """
    Convert units with a specific context.

    Contexts help with special conversions like:
    - 'spectroscopy': wavelength <-> frequency conversions
    - 'chemistry': molar concentration conversions

    Args:
        value: Numerical value to convert
        from_unit: Source unit
        to_unit: Target unit
        context: Context name (optional)

    Returns:
        Converted value
    """
    quantity = value * ureg(from_unit)
    if context:
        with ureg.context(context):
            result = quantity.to(to_unit)
    else:
        result = quantity.to(to_unit)
    return result.magnitude


# ============================================================================
# ENGINEERING CALCULATIONS
# ============================================================================

def reynolds_number(density: float, velocity: float, diameter: float, viscosity: float,
                   density_unit: str = 'kg/m^3', velocity_unit: str = 'm/s',
                   diameter_unit: str = 'm', viscosity_unit: str = 'Pa*s') -> float:
    """
    Calculate Reynolds number (dimensionless).

    Formula: Re = ρ × v × D / μ

    Args:
        density: Fluid density
        velocity: Flow velocity
        diameter: Characteristic length (pipe diameter)
        viscosity: Dynamic viscosity
        density_unit: Density unit (default: 'kg/m^3')
        velocity_unit: Velocity unit (default: 'm/s')
        diameter_unit: Length unit (default: 'm')
        viscosity_unit: Viscosity unit (default: 'Pa*s')

    Returns:
        Reynolds number (dimensionless)

    Example:
        >>> reynolds_number(1000, 2, 0.1, 0.001)
        200000.0
    """
    rho = density * ureg(density_unit)
    v = velocity * ureg(velocity_unit)
    D = diameter * ureg(diameter_unit)
    mu = viscosity * ureg(viscosity_unit)

    Re = (rho * v * D / mu).to_base_units()
    return Re.magnitude


def pressure_drop_pipe(flow_rate: float, diameter: float, length: float, roughness: float,
                      density: float, viscosity: float) -> float:
    """
    Simplified pressure drop calculation (for demonstration).

    This is a simplified Darcy-Weisbach equation example.
    For real applications, use more comprehensive hydraulic libraries.

    Args:
        flow_rate: Volumetric flow rate (m^3/s)
        diameter: Pipe diameter (m)
        length: Pipe length (m)
        roughness: Pipe roughness (m)
        density: Fluid density (kg/m^3)
        viscosity: Dynamic viscosity (Pa*s)

    Returns:
        Pressure drop (Pa)
    """
    # This is a placeholder - real calculation would need friction factor iteration
    # Just demonstrating unit handling
    Q = flow_rate * ureg('m^3/s')
    D = diameter * ureg('m')
    L = length * ureg('m')
    rho = density * ureg('kg/m^3')

    # Area
    A = 3.14159 * (D/2)**2

    # Velocity
    v = Q / A

    # Simplified pressure drop (not accurate, just for demonstration)
    dp = 0.02 * (L/D) * (rho * v**2 / 2)

    return dp.to('Pa').magnitude


# ============================================================================
# EXAMPLE USAGE AND TESTS
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("UNIT CONVERTER - EXAMPLE USAGE")
    print("=" * 70)

    # Flow conversions
    print("\n--- FLOW RATE CONVERSIONS ---")
    flow_lpm = 500
    print(f"{flow_lpm} L/min = {convert_flow(flow_lpm, 'L/min', 'm^3/s'):.4f} m³/s")
    print(f"{flow_lpm} L/min = {convert_flow(flow_lpm, 'L/min', 'm^3/h'):.2f} m³/h")
    print(f"{flow_lpm} L/min = {convert_flow(flow_lpm, 'L/min', 'gal/min'):.2f} gpm")

    flow_cfm = 1000
    print(f"\n{flow_cfm} cfm = {convert_flow(flow_cfm, 'ft^3/min', 'm^3/s'):.3f} m³/s")
    print(f"{flow_cfm} cfm = {convert_flow(flow_cfm, 'ft^3/min', 'L/s'):.1f} L/s")

    # Pressure conversions
    print("\n--- PRESSURE CONVERSIONS ---")
    press_bar = 150
    print(f"{press_bar} bar = {convert_pressure(press_bar, 'bar', 'Pa'):.0f} Pa")
    print(f"{press_bar} bar = {convert_pressure(press_bar, 'bar', 'kPa'):.0f} kPa")
    print(f"{press_bar} bar = {convert_pressure(press_bar, 'bar', 'psi'):.1f} psi")
    print(f"{press_bar} bar = {convert_pressure(press_bar, 'bar', 'atm'):.1f} atm")

    press_psi = 32
    print(f"\n{press_psi} psi = {convert_pressure(press_psi, 'psi', 'kPa'):.1f} kPa")
    print(f"{press_psi} psi = {convert_pressure(press_psi, 'psi', 'bar'):.2f} bar")

    # Gauge to absolute
    print(f"\n100 psig = {convert_gauge_to_absolute(100, 'psi'):.1f} psia")

    # Viscosity conversions
    print("\n--- VISCOSITY CONVERSIONS ---")
    visc_cp = 1.002
    print(f"{visc_cp} cP = {convert_viscosity_dynamic(visc_cp, 'cP', 'Pa*s'):.6f} Pa·s")

    visc_cst = 100
    print(f"{visc_cst} cSt = {convert_viscosity_kinematic(visc_cst, 'cSt', 'm^2/s'):.7f} m²/s")

    # Dynamic to kinematic
    nu = dynamic_to_kinematic(1.002e-3, 1000, 'Pa*s', 'kg/m^3', 'cSt')
    print(f"μ = 1.002 mPa·s, ρ = 1000 kg/m³ → ν = {nu:.3f} cSt")

    # Length conversions
    print("\n--- LENGTH CONVERSIONS ---")
    diam_inch = 4
    print(f"{diam_inch} inch = {convert_length(diam_inch, 'inch', 'mm'):.1f} mm")
    print(f"{diam_inch} inch = {convert_length(diam_inch, 'inch', 'm'):.4f} m")

    length_m = 15.5
    print(f"{length_m} m = {convert_length(length_m, 'm', 'ft'):.2f} ft")
    print(f"{length_m} m = {convert_length(length_m, 'm', 'inch'):.1f} in")

    # Power conversions
    print("\n--- POWER CONVERSIONS ---")
    power_hp = 50
    print(f"{power_hp} HP = {convert_power(power_hp, 'hp', 'kW'):.2f} kW")
    print(f"{power_hp} HP = {convert_power(power_hp, 'hp', 'W'):.0f} W")

    power_kw = 15
    print(f"{power_kw} kW = {convert_power(power_kw, 'kW', 'hp'):.2f} HP")

    # Temperature conversions
    print("\n--- TEMPERATURE CONVERSIONS ---")
    temp_c = 85
    print(f"{temp_c} °C = {convert_temperature(temp_c, 'degC', 'degF'):.1f} °F")
    print(f"{temp_c} °C = {convert_temperature(temp_c, 'degC', 'kelvin'):.2f} K")

    temp_k = 77
    print(f"{temp_k} K = {convert_temperature(temp_k, 'kelvin', 'degC'):.2f} °C")
    print(f"{temp_k} K = {convert_temperature(temp_k, 'kelvin', 'degF'):.2f} °F")

    # Velocity conversions
    print("\n--- VELOCITY CONVERSIONS ---")
    vel_ms = 2.5
    print(f"{vel_ms} m/s = {convert_velocity(vel_ms, 'm/s', 'ft/s'):.2f} ft/s")
    print(f"{vel_ms} m/s = {convert_velocity(vel_ms, 'm/s', 'km/h'):.1f} km/h")

    # Mass flow conversions
    print("\n--- MASS FLOW CONVERSIONS ---")
    mflow_kgh = 100
    print(f"{mflow_kgh} kg/h = {convert_mass_flow(mflow_kgh, 'kg/h', 'kg/s'):.5f} kg/s")
    print(f"{mflow_kgh} kg/h = {convert_mass_flow(mflow_kgh, 'kg/h', 'lb/min'):.3f} lb/min")

    # Dimensional consistency checking
    print("\n--- DIMENSIONAL CONSISTENCY ---")
    print(f"m/s and ft/min compatible? {check_dimensional_consistency('m/s', 'ft/min')}")
    print(f"Pa and psi compatible? {check_dimensional_consistency('Pa', 'psi')}")
    print(f"m and kg compatible? {check_dimensional_consistency('m', 'kg')}")

    # Get dimensionality
    print(f"\nDimensionality of Pa: {get_dimensionality('Pa')}")
    print(f"Dimensionality of m/s: {get_dimensionality('m/s')}")

    # Base units
    val, unit = get_base_units(100, 'psi')
    print(f"\n100 psi in base units: {val:.2f} {unit}")

    # Engineering calculation
    print("\n--- REYNOLDS NUMBER ---")
    Re = reynolds_number(1000, 2, 0.1, 0.001)
    print(f"Re = {Re:.0f}")

    print("\n" + "=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
