# Engineering Helper Functions and Utilities

This document provides a comprehensive collection of utility functions, helpers, and quick calculators for engineering calculations in Python.

## Table of Contents
- [Unit Conversion Utilities](#unit-conversion-utilities)
- [Engineering Constants](#engineering-constants)
- [Quick Calculators](#quick-calculators)
- [Resource Listers](#resource-listers)
- [Context Initialization](#context-initialization)
- [Error Handling](#error-handling)

---

## 1. Unit Conversion Utilities

### 1.1 Using Pint Library

The `pint` library provides comprehensive unit conversion capabilities with dimensional analysis.

#### Installation
```bash
pip install pint
```

#### Basic Usage

```python
import pint

# Create unit registry
ureg = pint.UnitRegistry()

# Define quantities with units
pressure = 150 * ureg.psi
print(f"Pressure: {pressure}")
print(f"Pressure in bar: {pressure.to(ureg.bar):.2f}")
print(f"Pressure in Pa: {pressure.to(ureg.Pa):.2f}")

# Temperature conversions
temp_f = 68 * ureg.degF
temp_c = temp_f.to(ureg.degC)
temp_k = temp_f.to(ureg.kelvin)
print(f"{temp_f} = {temp_c:.2f} = {temp_k:.2f}")

# Flow rate conversions
flow = 100 * ureg.gallon / ureg.minute
print(f"Flow: {flow}")
print(f"Flow in m³/h: {flow.to(ureg.meter**3 / ureg.hour):.4f}")
print(f"Flow in L/s: {flow.to(ureg.liter / ureg.second):.4f}")
```

#### Common Engineering Conversions

```python
from pint import UnitRegistry

ureg = UnitRegistry()

def convert_pressure(value, from_unit, to_unit):
    """
    Convert pressure between different units.

    Parameters:
    -----------
    value : float
        Pressure value
    from_unit : str
        Source unit (psi, bar, Pa, kPa, MPa, atm, torr, mmHg)
    to_unit : str
        Target unit

    Returns:
    --------
    float : Converted value
    """
    quantity = value * ureg(from_unit)
    return quantity.to(ureg(to_unit)).magnitude


def convert_flow(value, from_unit, to_unit):
    """
    Convert flow rate between different units.

    Parameters:
    -----------
    value : float
        Flow rate value
    from_unit : str
        Source unit (gpm, m3/h, L/s, ft3/s, etc.)
    to_unit : str
        Target unit

    Returns:
    --------
    float : Converted value
    """
    quantity = value * ureg(from_unit)
    return quantity.to(ureg(to_unit)).magnitude


def convert_temperature(value, from_unit, to_unit):
    """
    Convert temperature between different scales.

    Parameters:
    -----------
    value : float
        Temperature value
    from_unit : str
        Source unit (degC, degF, kelvin, degR)
    to_unit : str
        Target unit

    Returns:
    --------
    float : Converted value
    """
    quantity = value * ureg(from_unit)
    return quantity.to(ureg(to_unit)).magnitude


# Example usage
print(f"150 psi = {convert_pressure(150, 'psi', 'bar'):.2f} bar")
print(f"100 gpm = {convert_flow(100, 'gallon/minute', 'liter/second'):.2f} L/s")
print(f"25°C = {convert_temperature(25, 'degC', 'degF'):.2f}°F")
```

### 1.2 Custom Unit Converters

For specific engineering applications, custom converters can be useful:

```python
class EngineeringConverter:
    """Collection of common engineering unit conversions."""

    # Pressure conversions (all to Pa)
    PRESSURE_TO_PA = {
        'Pa': 1.0,
        'kPa': 1e3,
        'MPa': 1e6,
        'bar': 1e5,
        'psi': 6894.76,
        'atm': 101325,
        'torr': 133.322,
        'mmHg': 133.322,
        'inHg': 3386.39,
        'kgf/cm2': 98066.5
    }

    # Length conversions (all to meters)
    LENGTH_TO_M = {
        'm': 1.0,
        'cm': 0.01,
        'mm': 0.001,
        'km': 1000.0,
        'in': 0.0254,
        'ft': 0.3048,
        'yd': 0.9144,
        'mile': 1609.34
    }

    # Volume conversions (all to m³)
    VOLUME_TO_M3 = {
        'm3': 1.0,
        'L': 0.001,
        'mL': 1e-6,
        'gallon': 0.00378541,
        'ft3': 0.0283168,
        'in3': 1.63871e-5,
        'barrel': 0.158987
    }

    # Mass conversions (all to kg)
    MASS_TO_KG = {
        'kg': 1.0,
        'g': 0.001,
        'mg': 1e-6,
        'ton': 1000.0,
        'lb': 0.453592,
        'oz': 0.0283495
    }

    # Power conversions (all to W)
    POWER_TO_W = {
        'W': 1.0,
        'kW': 1000.0,
        'MW': 1e6,
        'hp': 745.7,
        'BTU/h': 0.293071
    }

    @classmethod
    def convert(cls, value, from_unit, to_unit, category):
        """
        Convert between units within a category.

        Parameters:
        -----------
        value : float
            Value to convert
        from_unit : str
            Source unit
        to_unit : str
            Target unit
        category : str
            Unit category (pressure, length, volume, mass, power)

        Returns:
        --------
        float : Converted value
        """
        conversion_dict = {
            'pressure': cls.PRESSURE_TO_PA,
            'length': cls.LENGTH_TO_M,
            'volume': cls.VOLUME_TO_M3,
            'mass': cls.MASS_TO_KG,
            'power': cls.POWER_TO_W
        }

        if category not in conversion_dict:
            raise ValueError(f"Unknown category: {category}")

        conv = conversion_dict[category]

        if from_unit not in conv:
            raise ValueError(f"Unknown unit: {from_unit}")
        if to_unit not in conv:
            raise ValueError(f"Unknown unit: {to_unit}")

        # Convert to base unit, then to target unit
        base_value = value * conv[from_unit]
        result = base_value / conv[to_unit]

        return result


# Example usage
converter = EngineeringConverter()
print(f"150 psi = {converter.convert(150, 'psi', 'bar', 'pressure'):.2f} bar")
print(f"12 in = {converter.convert(12, 'in', 'mm', 'length'):.2f} mm")
print(f"100 hp = {converter.convert(100, 'hp', 'kW', 'power'):.2f} kW")
```

### 1.3 Dimensional Analysis

```python
from pint import UnitRegistry

ureg = UnitRegistry()

def check_dimensional_consistency(equation_terms):
    """
    Check if all terms in an equation have consistent dimensions.

    Parameters:
    -----------
    equation_terms : list
        List of pint Quantity objects

    Returns:
    --------
    bool : True if consistent, False otherwise
    """
    if len(equation_terms) < 2:
        return True

    base_dim = equation_terms[0].dimensionality

    for term in equation_terms[1:]:
        if term.dimensionality != base_dim:
            print(f"Inconsistent dimensions: {base_dim} vs {term.dimensionality}")
            return False

    return True


# Example: Check Bernoulli's equation terms
P1 = 100000 * ureg.Pa  # Pressure
rho = 1000 * ureg.kg / ureg.m**3  # Density
v = 5 * ureg.m / ureg.s  # Velocity
g = 9.81 * ureg.m / ureg.s**2  # Gravity
h = 10 * ureg.m  # Height

term1 = P1 / rho  # Pressure head
term2 = 0.5 * v**2  # Velocity head
term3 = g * h  # Elevation head

print("Checking Bernoulli equation dimensional consistency:")
print(f"P/ρ: {term1.dimensionality}")
print(f"v²/2: {term2.dimensionality}")
print(f"gh: {term3.dimensionality}")
print(f"Consistent: {check_dimensional_consistency([term1, term2, term3])}")


def derive_units(formula_str, variable_units):
    """
    Derive the units of a result from a formula.

    Parameters:
    -----------
    formula_str : str
        Formula as string (e.g., 'P * V / (n * R * T)')
    variable_units : dict
        Dictionary mapping variable names to pint units

    Returns:
    --------
    pint.Unit : Resulting unit
    """
    # Evaluate the formula with unit quantities
    result = eval(formula_str, {"__builtins__": {}}, variable_units)
    return result.units


# Example: Ideal gas law PV = nRT, solve for P
vars_units = {
    'n': 1 * ureg.mol,
    'R': 8.314 * ureg.J / (ureg.mol * ureg.K),
    'T': 298 * ureg.K,
    'V': 0.001 * ureg.m**3
}

result_units = derive_units('n * R * T / V', vars_units)
print(f"\nPressure units from PV=nRT: {result_units}")
```

---

## 2. Engineering Constants

### 2.1 Physical Constants

```python
import numpy as np
from dataclasses import dataclass

@dataclass
class PhysicalConstants:
    """Collection of fundamental physical constants."""

    # Gravitational acceleration
    g_standard = 9.80665  # m/s²
    g_metric = 9.81  # m/s² (commonly used)

    # Universal gas constant
    R_universal = 8.314462618  # J/(mol·K)
    R_specific_air = 287.05  # J/(kg·K) for dry air

    # Standard atmospheric conditions
    P_atm = 101325  # Pa
    T_standard = 273.15  # K (0°C)
    T_stp = 273.15  # K (STP: 0°C, 1 atm)
    T_ntp = 293.15  # K (NTP: 20°C, 1 atm)

    # Water properties at STP
    rho_water_4C = 1000.0  # kg/m³ at 4°C (maximum density)
    rho_water_20C = 998.2  # kg/m³ at 20°C
    mu_water_20C = 1.002e-3  # Pa·s at 20°C

    # Air properties at STP
    rho_air_stp = 1.293  # kg/m³
    mu_air_20C = 1.825e-5  # Pa·s at 20°C

    # Mathematical constants
    pi = np.pi
    e = np.e

    # Conversion factors
    hp_to_W = 745.7
    BTU_to_J = 1055.06
    cal_to_J = 4.184


# Create global instance
PHYS = PhysicalConstants()

# Usage examples
print(f"Standard gravity: {PHYS.g_standard} m/s²")
print(f"Universal gas constant: {PHYS.R_universal} J/(mol·K)")
print(f"Water density at 20°C: {PHYS.rho_water_20C} kg/m³")
```

### 2.2 Standard Conditions

```python
class StandardConditions:
    """Standard condition definitions for various industries."""

    # STP - Standard Temperature and Pressure (IUPAC)
    STP = {
        'T': 273.15,  # K (0°C)
        'P': 100000,  # Pa (1 bar)
        'description': 'IUPAC Standard Temperature and Pressure'
    }

    # NTP - Normal Temperature and Pressure
    NTP = {
        'T': 293.15,  # K (20°C)
        'P': 101325,  # Pa (1 atm)
        'description': 'Normal Temperature and Pressure'
    }

    # Standard conditions for gas flow (USA)
    SCF = {
        'T': 288.71,  # K (60°F)
        'P': 101325,  # Pa (14.696 psi)
        'description': 'Standard Cubic Feet conditions'
    }

    # ISO standard reference conditions
    ISO = {
        'T': 288.15,  # K (15°C)
        'P': 101325,  # Pa (1 atm)
        'description': 'ISO standard reference conditions'
    }

    @staticmethod
    def correct_to_standard(Q_actual, T_actual, P_actual,
                           standard='STP'):
        """
        Correct volumetric flow rate to standard conditions.

        Parameters:
        -----------
        Q_actual : float
            Actual volumetric flow rate (m³/s)
        T_actual : float
            Actual temperature (K)
        P_actual : float
            Actual pressure (Pa)
        standard : str
            Standard condition set ('STP', 'NTP', 'SCF', 'ISO')

        Returns:
        --------
        float : Flow rate at standard conditions (Sm³/s)
        """
        std = getattr(StandardConditions, standard)

        # Using ideal gas law: P₁V₁/T₁ = P₂V₂/T₂
        Q_standard = Q_actual * (P_actual / std['P']) * (std['T'] / T_actual)

        return Q_standard


# Example usage
Q_actual = 10.0  # m³/s at operating conditions
T_actual = 350  # K
P_actual = 500000  # Pa

Q_stp = StandardConditions.correct_to_standard(Q_actual, T_actual, P_actual, 'STP')
print(f"Flow rate at STP: {Q_stp:.2f} Sm³/s")
```

### 2.3 Material Properties Database

```python
class MaterialProperties:
    """Common material properties for engineering calculations."""

    FLUIDS = {
        'water': {
            'density': 998.2,  # kg/m³ at 20°C
            'viscosity': 1.002e-3,  # Pa·s at 20°C
            'specific_heat': 4182,  # J/(kg·K)
            'thermal_conductivity': 0.598,  # W/(m·K)
        },
        'air': {
            'density': 1.204,  # kg/m³ at 20°C, 1 atm
            'viscosity': 1.825e-5,  # Pa·s at 20°C
            'specific_heat': 1005,  # J/(kg·K)
            'thermal_conductivity': 0.0257,  # W/(m·K)
        },
        'oil_sae30': {
            'density': 920,  # kg/m³
            'viscosity': 0.29,  # Pa·s at 20°C
            'specific_heat': 1880,  # J/(kg·K)
            'thermal_conductivity': 0.145,  # W/(m·K)
        },
        'glycol_50': {
            'density': 1070,  # kg/m³ (50% ethylene glycol)
            'viscosity': 5.5e-3,  # Pa·s at 20°C
            'specific_heat': 3350,  # J/(kg·K)
            'thermal_conductivity': 0.38,  # W/(m·K)
        }
    }

    METALS = {
        'steel_carbon': {
            'density': 7850,  # kg/m³
            'youngs_modulus': 200e9,  # Pa
            'yield_strength': 250e6,  # Pa
            'thermal_conductivity': 50,  # W/(m·K)
            'thermal_expansion': 12e-6,  # 1/K
        },
        'stainless_304': {
            'density': 8000,  # kg/m³
            'youngs_modulus': 193e9,  # Pa
            'yield_strength': 215e6,  # Pa
            'thermal_conductivity': 16.2,  # W/(m·K)
            'thermal_expansion': 17.3e-6,  # 1/K
        },
        'aluminum_6061': {
            'density': 2700,  # kg/m³
            'youngs_modulus': 69e9,  # Pa
            'yield_strength': 276e6,  # Pa
            'thermal_conductivity': 167,  # W/(m·K)
            'thermal_expansion': 23.6e-6,  # 1/K
        },
        'copper': {
            'density': 8960,  # kg/m³
            'youngs_modulus': 130e9,  # Pa
            'yield_strength': 70e6,  # Pa
            'thermal_conductivity': 401,  # W/(m·K)
            'thermal_expansion': 17e-6,  # 1/K
        }
    }

    @classmethod
    def get_fluid_properties(cls, fluid_name, temperature=20):
        """
        Get fluid properties, optionally adjusted for temperature.

        Parameters:
        -----------
        fluid_name : str
            Name of fluid
        temperature : float
            Temperature in °C (default: 20)

        Returns:
        --------
        dict : Fluid properties
        """
        if fluid_name not in cls.FLUIDS:
            raise ValueError(f"Unknown fluid: {fluid_name}")

        props = cls.FLUIDS[fluid_name].copy()

        # Note: These are reference values at 20°C
        # For accurate temperature-dependent properties, use CoolProp or similar
        if temperature != 20:
            print(f"Warning: Properties are at reference temperature 20°C")

        return props

    @classmethod
    def get_material_properties(cls, material_name):
        """
        Get metal/structural material properties.

        Parameters:
        -----------
        material_name : str
            Name of material

        Returns:
        --------
        dict : Material properties
        """
        if material_name not in cls.METALS:
            raise ValueError(f"Unknown material: {material_name}")

        return cls.METALS[material_name].copy()


# Example usage
water_props = MaterialProperties.get_fluid_properties('water')
print(f"Water density: {water_props['density']} kg/m³")
print(f"Water viscosity: {water_props['viscosity']} Pa·s")

steel_props = MaterialProperties.get_material_properties('steel_carbon')
print(f"Steel yield strength: {steel_props['yield_strength']/1e6:.0f} MPa")
```

---

## 3. Quick Calculators

### 3.1 Reynolds Number

```python
def reynolds_number(velocity, characteristic_length, kinematic_viscosity=None,
                   density=None, dynamic_viscosity=None):
    """
    Calculate Reynolds number for flow analysis.

    Re = ρVL/μ = VL/ν

    Parameters:
    -----------
    velocity : float
        Flow velocity (m/s)
    characteristic_length : float
        Characteristic length (m) - diameter for pipes
    kinematic_viscosity : float, optional
        Kinematic viscosity (m²/s)
    density : float, optional
        Fluid density (kg/m³)
    dynamic_viscosity : float, optional
        Dynamic viscosity (Pa·s)

    Returns:
    --------
    float : Reynolds number (dimensionless)

    Notes:
    ------
    - Re < 2300: Laminar flow (pipes)
    - 2300 < Re < 4000: Transition
    - Re > 4000: Turbulent flow
    """
    if kinematic_viscosity is not None:
        Re = velocity * characteristic_length / kinematic_viscosity
    elif density is not None and dynamic_viscosity is not None:
        Re = density * velocity * characteristic_length / dynamic_viscosity
    else:
        raise ValueError("Must provide either kinematic_viscosity or (density and dynamic_viscosity)")

    return Re


def flow_regime(Re, geometry='pipe'):
    """
    Determine flow regime based on Reynolds number.

    Parameters:
    -----------
    Re : float
        Reynolds number
    geometry : str
        Flow geometry ('pipe', 'flat_plate', 'sphere')

    Returns:
    --------
    str : Flow regime description
    """
    if geometry == 'pipe':
        if Re < 2300:
            return 'Laminar'
        elif Re < 4000:
            return 'Transitional'
        else:
            return 'Turbulent'
    elif geometry == 'flat_plate':
        if Re < 5e5:
            return 'Laminar'
        elif Re < 1e6:
            return 'Transitional'
        else:
            return 'Turbulent'
    elif geometry == 'sphere':
        if Re < 1:
            return 'Stokes (creeping) flow'
        elif Re < 2000:
            return 'Laminar'
        else:
            return 'Turbulent'
    else:
        raise ValueError(f"Unknown geometry: {geometry}")


# Example usage
V = 2.0  # m/s
D = 0.15  # m (pipe diameter)
rho = 998.2  # kg/m³ (water)
mu = 1.002e-3  # Pa·s (water at 20°C)

Re = reynolds_number(V, D, density=rho, dynamic_viscosity=mu)
regime = flow_regime(Re)
print(f"Reynolds number: {Re:.0f}")
print(f"Flow regime: {regime}")
```

### 3.2 Friction Factor (Darcy-Weisbach)

```python
import numpy as np
from scipy.optimize import fsolve

def friction_factor_laminar(Re):
    """
    Friction factor for laminar flow.

    f = 64/Re
    """
    return 64.0 / Re


def friction_factor_turbulent_smooth(Re):
    """
    Friction factor for turbulent flow in smooth pipes (Blasius).

    f = 0.316/Re^0.25  (valid for Re < 10^5)
    """
    if Re < 4000:
        raise ValueError("Blasius equation not valid for laminar flow")

    return 0.316 / (Re ** 0.25)


def friction_factor_colebrook(Re, relative_roughness):
    """
    Colebrook-White equation for friction factor (implicit).

    1/√f = -2 log₁₀(ε/3.7D + 2.51/(Re√f))

    Parameters:
    -----------
    Re : float
        Reynolds number
    relative_roughness : float
        ε/D (absolute roughness / diameter)

    Returns:
    --------
    float : Darcy friction factor
    """
    def colebrook_eq(f):
        return (1.0 / np.sqrt(f) +
                2.0 * np.log10(relative_roughness / 3.7 + 2.51 / (Re * np.sqrt(f))))

    # Initial guess using Swamee-Jain
    f_guess = 0.25 / (np.log10(relative_roughness / 3.7 + 5.74 / (Re ** 0.9))) ** 2

    # Solve implicit equation
    f = fsolve(colebrook_eq, f_guess)[0]

    return f


def friction_factor_swamee_jain(Re, relative_roughness):
    """
    Swamee-Jain explicit approximation of Colebrook equation.

    Accurate to within 1% for:
    - 5000 < Re < 10^8
    - 10^-6 < ε/D < 10^-2

    Parameters:
    -----------
    Re : float
        Reynolds number
    relative_roughness : float
        ε/D (absolute roughness / diameter)

    Returns:
    --------
    float : Darcy friction factor
    """
    numerator = 0.25
    denominator = (np.log10(relative_roughness / 3.7 + 5.74 / (Re ** 0.9))) ** 2

    return numerator / denominator


def friction_factor(Re, relative_roughness=0, method='auto'):
    """
    Calculate Darcy friction factor.

    Parameters:
    -----------
    Re : float
        Reynolds number
    relative_roughness : float
        ε/D, default 0 (smooth pipe)
    method : str
        'auto', 'laminar', 'blasius', 'colebrook', 'swamee-jain'

    Returns:
    --------
    float : Darcy friction factor
    """
    if method == 'auto':
        if Re < 2300:
            method = 'laminar'
        elif relative_roughness == 0 and Re < 1e5:
            method = 'blasius'
        else:
            method = 'swamee-jain'

    if method == 'laminar':
        return friction_factor_laminar(Re)
    elif method == 'blasius':
        return friction_factor_turbulent_smooth(Re)
    elif method == 'colebrook':
        return friction_factor_colebrook(Re, relative_roughness)
    elif method == 'swamee-jain':
        return friction_factor_swamee_jain(Re, relative_roughness)
    else:
        raise ValueError(f"Unknown method: {method}")


# Example usage
Re = 50000
epsilon_D = 0.0002  # Relative roughness for commercial steel

f = friction_factor(Re, epsilon_D)
print(f"Friction factor (Swamee-Jain): {f:.6f}")

f_colebrook = friction_factor(Re, epsilon_D, method='colebrook')
print(f"Friction factor (Colebrook): {f_colebrook:.6f}")
```

### 3.3 NPSH (Net Positive Suction Head)

```python
def npsh_available(P_atm, P_vapor, z_suction, v_suction, rho, g=9.81):
    """
    Calculate available NPSH for pump suction.

    NPSH_a = (P_atm - P_vapor)/ρg + z_suction + v²/(2g) - h_friction

    Parameters:
    -----------
    P_atm : float
        Atmospheric pressure (Pa)
    P_vapor : float
        Vapor pressure of liquid at pumping temperature (Pa)
    z_suction : float
        Height of liquid surface above pump centerline (m)
        Positive if above, negative if below
    v_suction : float
        Velocity in suction pipe (m/s)
    rho : float
        Liquid density (kg/m³)
    g : float
        Gravitational acceleration (m/s²), default 9.81

    Returns:
    --------
    float : NPSH available (m)

    Notes:
    ------
    NPSH_a must exceed NPSH_r (required) to avoid cavitation.
    Typically maintain NPSH_a > NPSH_r + 1m safety margin.
    """
    pressure_head = (P_atm - P_vapor) / (rho * g)
    elevation_head = z_suction
    velocity_head = (v_suction ** 2) / (2 * g)

    NPSH_a = pressure_head + elevation_head + velocity_head

    return NPSH_a


def npsh_check(NPSH_a, NPSH_r, safety_margin=1.0):
    """
    Check if available NPSH is adequate.

    Parameters:
    -----------
    NPSH_a : float
        Available NPSH (m)
    NPSH_r : float
        Required NPSH from pump curve (m)
    safety_margin : float
        Additional safety margin (m), default 1.0

    Returns:
    --------
    dict : Check results with status and margin
    """
    margin = NPSH_a - NPSH_r
    adequate = margin >= safety_margin

    return {
        'adequate': adequate,
        'margin': margin,
        'status': 'OK' if adequate else 'INSUFFICIENT',
        'message': (f"NPSH margin: {margin:.2f} m "
                   f"(required margin: {safety_margin:.2f} m)")
    }


# Example usage - Water pump at sea level, 20°C
P_atm = 101325  # Pa (1 atm)
P_vapor = 2339  # Pa (vapor pressure of water at 20°C)
z = 2.0  # m (2m above pump)
v = 1.5  # m/s
rho = 998.2  # kg/m³
NPSH_r = 3.5  # m (from pump curve)

NPSH_a = npsh_available(P_atm, P_vapor, z, v, rho)
print(f"NPSH available: {NPSH_a:.2f} m")

check = npsh_check(NPSH_a, NPSH_r)
print(f"Status: {check['status']}")
print(check['message'])
```

### 3.4 Specific Speed

```python
def specific_speed_pump(Q, H, N, units='SI'):
    """
    Calculate pump specific speed.

    SI: N_s = N√Q / H^(3/4)
    US: N_s = N√GPM / H^(3/4)

    Parameters:
    -----------
    Q : float
        Flow rate (m³/s for SI, GPM for US)
    H : float
        Head (m for SI, ft for US)
    N : float
        Rotational speed (rpm)
    units : str
        'SI' or 'US'

    Returns:
    --------
    float : Specific speed

    Notes:
    ------
    Typical ranges:
    - N_s < 1000 (SI): Centrifugal (radial)
    - 1000-2000 (SI): Mixed flow
    - N_s > 2000 (SI): Axial flow
    """
    N_s = N * np.sqrt(Q) / (H ** 0.75)

    return N_s


def pump_type_from_specific_speed(N_s, units='SI'):
    """
    Determine pump type from specific speed.

    Parameters:
    -----------
    N_s : float
        Specific speed
    units : str
        'SI' or 'US'

    Returns:
    --------
    str : Pump type recommendation
    """
    if units == 'SI':
        if N_s < 500:
            return "Radial flow (centrifugal) - Low flow, high head"
        elif N_s < 1000:
            return "Radial flow (centrifugal) - Normal range"
        elif N_s < 2000:
            return "Mixed flow - Medium flow, medium head"
        elif N_s < 4000:
            return "Axial flow - High flow, low head"
        else:
            return "Axial flow (propeller) - Very high flow, very low head"
    else:  # US units
        if N_s < 2000:
            return "Radial flow (centrifugal) - Low flow, high head"
        elif N_s < 4000:
            return "Radial flow (centrifugal) - Normal range"
        elif N_s < 7000:
            return "Mixed flow - Medium flow, medium head"
        elif N_s < 15000:
            return "Axial flow - High flow, low head"
        else:
            return "Axial flow (propeller) - Very high flow, very low head"


# Example usage
Q = 0.050  # m³/s
H = 50  # m
N = 1750  # rpm

N_s = specific_speed_pump(Q, H, N, units='SI')
pump_type = pump_type_from_specific_speed(N_s, units='SI')

print(f"Specific speed: {N_s:.0f} (SI)")
print(f"Recommended pump type: {pump_type}")
```

### 3.5 Affinity Laws

```python
class AffinityLaws:
    """
    Pump and fan affinity laws for scaling performance.

    The affinity laws relate changes in:
    - Speed (N)
    - Impeller diameter (D)
    - Flow rate (Q)
    - Head/Pressure (H)
    - Power (P)
    """

    @staticmethod
    def scale_by_speed(N1, N2, Q1=None, H1=None, P1=None):
        """
        Scale pump/fan performance with speed change.

        Laws:
        Q₂/Q₁ = N₂/N₁
        H₂/H₁ = (N₂/N₁)²
        P₂/P₁ = (N₂/N₁)³

        Parameters:
        -----------
        N1 : float
            Original speed (rpm)
        N2 : float
            New speed (rpm)
        Q1 : float, optional
            Original flow rate
        H1 : float, optional
            Original head
        P1 : float, optional
            Original power

        Returns:
        --------
        dict : Scaled values (Q2, H2, P2)
        """
        ratio = N2 / N1

        result = {
            'N1': N1,
            'N2': N2,
            'speed_ratio': ratio
        }

        if Q1 is not None:
            result['Q1'] = Q1
            result['Q2'] = Q1 * ratio

        if H1 is not None:
            result['H1'] = H1
            result['H2'] = H1 * (ratio ** 2)

        if P1 is not None:
            result['P1'] = P1
            result['P2'] = P1 * (ratio ** 3)

        return result

    @staticmethod
    def scale_by_diameter(D1, D2, Q1=None, H1=None, P1=None):
        """
        Scale pump performance with impeller diameter change.

        Laws:
        Q₂/Q₁ = (D₂/D₁)³
        H₂/H₁ = (D₂/D₁)²
        P₂/P₁ = (D₂/D₁)⁵

        Parameters:
        -----------
        D1 : float
            Original diameter
        D2 : float
            New diameter
        Q1 : float, optional
            Original flow rate
        H1 : float, optional
            Original head
        P1 : float, optional
            Original power

        Returns:
        --------
        dict : Scaled values (Q2, H2, P2)
        """
        ratio = D2 / D1

        result = {
            'D1': D1,
            'D2': D2,
            'diameter_ratio': ratio
        }

        if Q1 is not None:
            result['Q1'] = Q1
            result['Q2'] = Q1 * (ratio ** 3)

        if H1 is not None:
            result['H1'] = H1
            result['H2'] = H1 * (ratio ** 2)

        if P1 is not None:
            result['P1'] = P1
            result['P2'] = P1 * (ratio ** 5)

        return result

    @staticmethod
    def find_speed_for_flow(N1, Q1, Q_target):
        """
        Find required speed to achieve target flow rate.

        Parameters:
        -----------
        N1 : float
            Original speed (rpm)
        Q1 : float
            Original flow rate
        Q_target : float
            Target flow rate

        Returns:
        --------
        float : Required speed (rpm)
        """
        N2 = N1 * (Q_target / Q1)
        return N2

    @staticmethod
    def find_speed_for_head(N1, H1, H_target):
        """
        Find required speed to achieve target head.

        Parameters:
        -----------
        N1 : float
            Original speed (rpm)
        H1 : float
            Original head
        H_target : float
            Target head

        Returns:
        --------
        float : Required speed (rpm)
        """
        N2 = N1 * np.sqrt(H_target / H1)
        return N2


# Example usage: Speed change
print("=== Speed Change Example ===")
result = AffinityLaws.scale_by_speed(
    N1=1750, N2=1450,
    Q1=100, H1=50, P1=15
)
print(f"Original: {result['N1']} rpm, {result['Q1']} m³/h, {result['H1']} m, {result['P1']} kW")
print(f"New: {result['N2']} rpm, {result['Q2']:.1f} m³/h, {result['H2']:.1f} m, {result['P2']:.1f} kW")

# Example usage: Impeller trim
print("\n=== Impeller Diameter Change ===")
result = AffinityLaws.scale_by_diameter(
    D1=250, D2=225,  # mm
    Q1=100, H1=50, P1=15
)
print(f"Original diameter: {result['D1']} mm")
print(f"Trimmed diameter: {result['D2']} mm")
print(f"New flow: {result['Q2']:.1f} m³/h ({result['Q2']/result['Q1']*100:.1f}% of original)")
print(f"New head: {result['H2']:.1f} m ({result['H2']/result['H1']*100:.1f}% of original)")
print(f"New power: {result['P2']:.1f} kW ({result['P2']/result['P1']*100:.1f}% of original)")

# Example: Find speed for target flow
print("\n=== Find Required Speed ===")
N_required = AffinityLaws.find_speed_for_flow(N1=1750, Q1=100, Q_target=120)
print(f"To increase flow from 100 to 120 m³/h:")
print(f"Increase speed from 1750 to {N_required:.0f} rpm")
```

---

## 4. Resource Listers

### 4.1 Package Enumeration

```python
import sys
import pkg_resources
import importlib.util

def list_installed_packages():
    """
    List all installed Python packages.

    Returns:
    --------
    list : List of (package_name, version) tuples
    """
    installed_packages = [(d.project_name, d.version)
                         for d in pkg_resources.working_set]
    installed_packages.sort()

    return installed_packages


def check_package_available(package_name):
    """
    Check if a package is available for import.

    Parameters:
    -----------
    package_name : str
        Name of package to check

    Returns:
    --------
    dict : Availability status and version info
    """
    spec = importlib.util.find_spec(package_name)
    available = spec is not None

    result = {
        'package': package_name,
        'available': available,
        'version': None
    }

    if available:
        try:
            module = importlib.import_module(package_name)
            if hasattr(module, '__version__'):
                result['version'] = module.__version__
        except:
            pass

    return result


def check_engineering_packages():
    """
    Check availability of common engineering packages.

    Returns:
    --------
    dict : Status of engineering packages
    """
    packages = [
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'sympy',
        'pint',
        'CoolProp',
        'fluids',
        'thermo',
        'ht',
        'chemicals'
    ]

    results = {}
    for package in packages:
        results[package] = check_package_available(package)

    return results


def print_package_status(status_dict):
    """Print formatted package status."""
    print("\n{:<15} {:<12} {:<15}".format("Package", "Available", "Version"))
    print("-" * 45)

    for package, info in status_dict.items():
        available = "✓" if info['available'] else "✗"
        version = info['version'] if info['version'] else "N/A"
        print(f"{package:<15} {available:<12} {version:<15}")


# Example usage
print("=== Engineering Package Status ===")
status = check_engineering_packages()
print_package_status(status)

# List all packages (first 10)
print("\n=== Sample of Installed Packages ===")
all_packages = list_installed_packages()
for name, version in all_packages[:10]:
    print(f"{name:<30} {version}")
print(f"... and {len(all_packages) - 10} more packages")
```

### 4.2 Database Connection Checker

```python
import os

def check_database_connections():
    """
    Check for database connection configuration.

    Returns:
    --------
    dict : Database availability status
    """
    results = {}

    # Check SQLite (always available in Python)
    try:
        import sqlite3
        results['sqlite'] = {
            'available': True,
            'version': sqlite3.sqlite_version,
            'module_version': sqlite3.version
        }
    except:
        results['sqlite'] = {'available': False}

    # Check PostgreSQL (psycopg2)
    try:
        import psycopg2
        results['postgresql'] = {
            'available': True,
            'module': 'psycopg2',
            'version': psycopg2.__version__
        }
    except:
        results['postgresql'] = {'available': False}

    # Check MySQL
    try:
        import pymysql
        results['mysql'] = {
            'available': True,
            'module': 'pymysql',
            'version': pymysql.__version__
        }
    except:
        try:
            import mysql.connector
            results['mysql'] = {
                'available': True,
                'module': 'mysql.connector',
                'version': mysql.connector.__version__
            }
        except:
            results['mysql'] = {'available': False}

    # Check MongoDB
    try:
        import pymongo
        results['mongodb'] = {
            'available': True,
            'module': 'pymongo',
            'version': pymongo.__version__
        }
    except:
        results['mongodb'] = {'available': False}

    return results


def test_sqlite_connection(db_path=':memory:'):
    """
    Test SQLite connection.

    Parameters:
    -----------
    db_path : str
        Path to database file, or ':memory:' for in-memory

    Returns:
    --------
    dict : Connection test results
    """
    import sqlite3

    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT sqlite_version()")
        version = cursor.fetchone()[0]
        conn.close()

        return {
            'success': True,
            'database': db_path,
            'version': version
        }
    except Exception as e:
        return {
            'success': False,
            'database': db_path,
            'error': str(e)
        }


# Example usage
print("=== Database Module Status ===")
db_status = check_database_connections()

for db, info in db_status.items():
    if info['available']:
        version = info.get('version', 'N/A')
        module = info.get('module', db)
        print(f"{db:<15} ✓ ({module} v{version})")
    else:
        print(f"{db:<15} ✗ (not available)")

# Test SQLite
print("\n=== SQLite Connection Test ===")
test_result = test_sqlite_connection()
if test_result['success']:
    print(f"SQLite connection successful (v{test_result['version']})")
else:
    print(f"SQLite connection failed: {test_result['error']}")
```

### 4.3 Software Installation Verifier

```python
import subprocess
import shutil

def check_system_tools():
    """
    Check for system tools and software.

    Returns:
    --------
    dict : Tool availability status
    """
    tools = {
        'git': 'git --version',
        'python': 'python --version',
        'pip': 'pip --version',
        'gcc': 'gcc --version',
        'gfortran': 'gfortran --version',
        'make': 'make --version'
    }

    results = {}

    for tool, command in tools.items():
        # Check if executable exists in PATH
        executable = shutil.which(tool)

        if executable is None:
            results[tool] = {
                'available': False,
                'path': None,
                'version': None
            }
        else:
            results[tool] = {
                'available': True,
                'path': executable,
                'version': None
            }

            # Try to get version
            try:
                output = subprocess.check_output(
                    command.split(),
                    stderr=subprocess.STDOUT,
                    timeout=5
                ).decode('utf-8').split('\n')[0]
                results[tool]['version'] = output
            except:
                pass

    return results


def print_tool_status(status_dict):
    """Print formatted tool status."""
    print("\n{:<15} {:<12} {:<50}".format("Tool", "Available", "Version/Path"))
    print("-" * 80)

    for tool, info in status_dict.items():
        available = "✓" if info['available'] else "✗"
        detail = info['version'] if info['version'] else info['path']
        if detail is None:
            detail = "Not found"
        print(f"{tool:<15} {available:<12} {detail}")


# Example usage
print("=== System Tool Status ===")
tool_status = check_system_tools()
print_tool_status(tool_status)
```

---

## 5. Context Initialization

### 5.1 Session Setup

```python
def init_engineering_session(verbose=True):
    """
    Initialize an engineering calculation session with common imports.

    Parameters:
    -----------
    verbose : bool
        Print status messages

    Returns:
    --------
    dict : Dictionary of imported modules and setup status
    """
    modules = {}

    # Core scientific computing
    try:
        import numpy as np
        modules['np'] = np
        if verbose:
            print(f"✓ NumPy {np.__version__}")
    except ImportError:
        if verbose:
            print("✗ NumPy not available")

    try:
        import scipy
        modules['scipy'] = scipy
        from scipy import optimize, integrate, interpolate
        modules['optimize'] = optimize
        modules['integrate'] = integrate
        modules['interpolate'] = interpolate
        if verbose:
            print(f"✓ SciPy {scipy.__version__}")
    except ImportError:
        if verbose:
            print("✗ SciPy not available")

    # Plotting
    try:
        import matplotlib.pyplot as plt
        modules['plt'] = plt
        plt.style.use('seaborn-v0_8-darkgrid')
        if verbose:
            print(f"✓ Matplotlib (style: seaborn-darkgrid)")
    except ImportError:
        if verbose:
            print("✗ Matplotlib not available")

    # Data analysis
    try:
        import pandas as pd
        modules['pd'] = pd
        pd.set_option('display.precision', 4)
        pd.set_option('display.max_rows', 20)
        if verbose:
            print(f"✓ Pandas {pd.__version__}")
    except ImportError:
        if verbose:
            print("✗ Pandas not available")

    # Symbolic math
    try:
        import sympy as sym
        modules['sym'] = sym
        sym.init_printing()
        if verbose:
            print(f"✓ SymPy {sym.__version__}")
    except ImportError:
        if verbose:
            print("✗ SymPy not available")

    # Units
    try:
        import pint
        ureg = pint.UnitRegistry()
        modules['ureg'] = ureg
        if verbose:
            print(f"✓ Pint {pint.__version__}")
    except ImportError:
        if verbose:
            print("✗ Pint not available")

    # Engineering libraries
    try:
        import CoolProp
        modules['CoolProp'] = CoolProp
        if verbose:
            print(f"✓ CoolProp {CoolProp.__version__}")
    except ImportError:
        if verbose:
            print("✗ CoolProp not available")

    try:
        import fluids
        modules['fluids'] = fluids
        if verbose:
            print(f"✓ fluids {fluids.__version__}")
    except ImportError:
        if verbose:
            print("✗ fluids not available")

    if verbose:
        print(f"\n✓ Engineering session initialized with {len(modules)} modules")

    return modules


# Example usage
print("=== Initializing Engineering Session ===")
session = init_engineering_session(verbose=True)

# Access modules from session dictionary
if 'np' in session:
    np = session['np']
    print(f"\nNumPy available: {np.pi:.6f}")
```

### 5.2 Default Units Configuration

```python
class EngineeringContext:
    """
    Manage engineering calculation context with default units.
    """

    def __init__(self):
        """Initialize with SI units as default."""
        self.unit_system = 'SI'
        self.units = self._get_si_units()
        self.precision = 4

    def _get_si_units(self):
        """Get SI unit definitions."""
        return {
            'length': 'm',
            'mass': 'kg',
            'time': 's',
            'temperature': 'K',
            'pressure': 'Pa',
            'energy': 'J',
            'power': 'W',
            'force': 'N',
            'volume': 'm³',
            'flow_rate': 'm³/s',
            'velocity': 'm/s',
            'acceleration': 'm/s²',
            'density': 'kg/m³',
            'viscosity_dynamic': 'Pa·s',
            'viscosity_kinematic': 'm²/s'
        }

    def _get_imperial_units(self):
        """Get Imperial/US unit definitions."""
        return {
            'length': 'ft',
            'mass': 'lb',
            'time': 's',
            'temperature': 'degF',
            'pressure': 'psi',
            'energy': 'BTU',
            'power': 'hp',
            'force': 'lbf',
            'volume': 'ft³',
            'flow_rate': 'gpm',
            'velocity': 'ft/s',
            'acceleration': 'ft/s²',
            'density': 'lb/ft³',
            'viscosity_dynamic': 'lbf·s/ft²',
            'viscosity_kinematic': 'ft²/s'
        }

    def set_unit_system(self, system):
        """
        Set unit system.

        Parameters:
        -----------
        system : str
            'SI' or 'Imperial'
        """
        if system.upper() == 'SI':
            self.unit_system = 'SI'
            self.units = self._get_si_units()
        elif system.upper() in ['IMPERIAL', 'US']:
            self.unit_system = 'Imperial'
            self.units = self._get_imperial_units()
        else:
            raise ValueError(f"Unknown unit system: {system}")

    def get_unit(self, quantity):
        """
        Get default unit for a quantity.

        Parameters:
        -----------
        quantity : str
            Physical quantity name

        Returns:
        --------
        str : Unit string
        """
        if quantity not in self.units:
            raise ValueError(f"Unknown quantity: {quantity}")

        return self.units[quantity]

    def set_precision(self, digits):
        """Set default precision for output."""
        self.precision = digits

    def format_value(self, value, quantity=None):
        """
        Format value with appropriate units.

        Parameters:
        -----------
        value : float
            Numerical value
        quantity : str, optional
            Physical quantity for unit lookup

        Returns:
        --------
        str : Formatted string with units
        """
        formatted = f"{value:.{self.precision}f}"

        if quantity and quantity in self.units:
            formatted += f" {self.units[quantity]}"

        return formatted

    def summary(self):
        """Print context summary."""
        print(f"Unit System: {self.unit_system}")
        print(f"Precision: {self.precision} decimal places")
        print("\nDefault Units:")
        for quantity, unit in self.units.items():
            print(f"  {quantity:<25} {unit}")


# Example usage
print("=== Engineering Context ===")
context = EngineeringContext()
context.summary()

print("\n=== Formatting Examples ===")
print(f"Pressure: {context.format_value(101325, 'pressure')}")
print(f"Flow rate: {context.format_value(0.05, 'flow_rate')}")
print(f"Temperature: {context.format_value(298.15, 'temperature')}")

print("\n=== Switch to Imperial ===")
context.set_unit_system('Imperial')
print(f"Pressure: {context.format_value(14.7, 'pressure')}")
print(f"Flow rate: {context.format_value(100, 'flow_rate')}")
```

### 5.3 Plotting Configuration

```python
import matplotlib.pyplot as plt
import numpy as np

def setup_engineering_plots(style='default'):
    """
    Configure matplotlib for engineering plots.

    Parameters:
    -----------
    style : str
        'default', 'presentation', 'publication', 'dark'
    """
    # Base configuration
    plt.rcParams['figure.figsize'] = (10, 6)
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'

    # Font settings
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 10

    # Grid
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['grid.linestyle'] = '--'

    # Lines
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['lines.markersize'] = 6

    if style == 'presentation':
        # Larger fonts for presentations
        plt.rcParams['font.size'] = 14
        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams['axes.titlesize'] = 18
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['lines.linewidth'] = 3

    elif style == 'publication':
        # Publication-ready settings
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['font.serif'] = ['Times New Roman']
        plt.rcParams['mathtext.fontset'] = 'stix'
        plt.rcParams['figure.figsize'] = (8, 6)
        plt.rcParams['lines.linewidth'] = 1.5

    elif style == 'dark':
        # Dark theme
        plt.style.use('dark_background')
        plt.rcParams['grid.alpha'] = 0.2

    print(f"✓ Plot settings configured for: {style}")


def engineering_plot_template(x_data, y_data,
                              xlabel, ylabel, title,
                              grid=True, legend_labels=None):
    """
    Create a standard engineering plot.

    Parameters:
    -----------
    x_data : array or list of arrays
        X-axis data
    y_data : array or list of arrays
        Y-axis data
    xlabel : str
        X-axis label with units
    ylabel : str
        Y-axis label with units
    title : str
        Plot title
    grid : bool
        Show grid
    legend_labels : list, optional
        Legend labels for multiple datasets

    Returns:
    --------
    fig, ax : matplotlib figure and axes objects
    """
    fig, ax = plt.subplots()

    # Handle multiple datasets
    if isinstance(y_data, list):
        for i, (x, y) in enumerate(zip(x_data, y_data)):
            label = legend_labels[i] if legend_labels else None
            ax.plot(x, y, marker='o', label=label)
        if legend_labels:
            ax.legend()
    else:
        ax.plot(x_data, y_data, marker='o')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if grid:
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, ax


# Example usage
print("=== Configuring Engineering Plots ===")
setup_engineering_plots(style='default')

# Create sample plot
x = np.linspace(0, 10, 50)
y1 = np.exp(-x/5) * np.cos(2*x)
y2 = np.exp(-x/5) * np.sin(2*x)

fig, ax = engineering_plot_template(
    [x, x], [y1, y2],
    xlabel='Time (s)',
    ylabel='Response (mm)',
    title='System Response',
    legend_labels=['Damped Cosine', 'Damped Sine']
)

# plt.savefig('system_response.png')
plt.close()
print("✓ Example plot created")
```

---

## 6. Error Handling

### 6.1 Numerical Issue Detection

```python
import numpy as np
import warnings

def check_numerical_issues(value, name="value"):
    """
    Check for common numerical issues.

    Parameters:
    -----------
    value : float or array
        Value(s) to check
    name : str
        Variable name for error messages

    Returns:
    --------
    dict : Dictionary of detected issues
    """
    issues = {
        'is_nan': False,
        'is_inf': False,
        'is_negative': False,
        'is_zero': False,
        'out_of_range': False,
        'messages': []
    }

    if isinstance(value, (list, tuple)):
        value = np.array(value)

    # Check for NaN
    if np.any(np.isnan(value)):
        issues['is_nan'] = True
        issues['messages'].append(f"{name} contains NaN values")

    # Check for infinity
    if np.any(np.isinf(value)):
        issues['is_inf'] = True
        issues['messages'].append(f"{name} contains infinite values")

    # Check for negative (when not expected)
    if np.any(value < 0):
        issues['is_negative'] = True

    # Check for zero
    if np.any(value == 0):
        issues['is_zero'] = True

    # Check for very small or large values
    if not issues['is_inf'] and not issues['is_nan']:
        if np.any(np.abs(value) < 1e-100):
            issues['out_of_range'] = True
            issues['messages'].append(f"{name} contains very small values (<1e-100)")
        if np.any(np.abs(value) > 1e100):
            issues['out_of_range'] = True
            issues['messages'].append(f"{name} contains very large values (>1e100)")

    return issues


def safe_divide(numerator, denominator, default=np.nan):
    """
    Safe division with handling of zero denominator.

    Parameters:
    -----------
    numerator : float or array
        Numerator
    denominator : float or array
        Denominator
    default : float
        Value to return when denominator is zero

    Returns:
    --------
    float or array : Result of division
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.where(denominator != 0,
                         numerator / denominator,
                         default)
    return result


def safe_log(value, base=np.e, default=np.nan):
    """
    Safe logarithm with handling of non-positive values.

    Parameters:
    -----------
    value : float or array
        Input value
    base : float
        Logarithm base (default: e)
    default : float
        Value to return for non-positive inputs

    Returns:
    --------
    float or array : Logarithm result
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if base == np.e:
            result = np.where(value > 0, np.log(value), default)
        elif base == 10:
            result = np.where(value > 0, np.log10(value), default)
        else:
            result = np.where(value > 0,
                            np.log(value) / np.log(base),
                            default)
    return result


def safe_sqrt(value, default=np.nan):
    """
    Safe square root with handling of negative values.

    Parameters:
    -----------
    value : float or array
        Input value
    default : float
        Value to return for negative inputs

    Returns:
    --------
    float or array : Square root result
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = np.where(value >= 0, np.sqrt(value), default)
    return result


# Example usage
print("=== Numerical Issue Detection ===")

# Test various problematic values
test_values = [1.0, 0.0, -1.0, np.nan, np.inf, 1e-150, 1e150]

for val in test_values:
    issues = check_numerical_issues(val, name=f"value={val}")
    if issues['messages']:
        print(f"{val}: {', '.join(issues['messages'])}")
    else:
        print(f"{val}: OK")

# Test safe operations
print("\n=== Safe Operations ===")
print(f"Safe divide 10/0: {safe_divide(10, 0, default=0)}")
print(f"Safe log of -1: {safe_log(-1, default=0)}")
print(f"Safe sqrt of -4: {safe_sqrt(-4, default=0)}")
```

### 6.2 Physical Validity Checks

```python
class PhysicalValidator:
    """Validate physical parameters for engineering calculations."""

    @staticmethod
    def validate_temperature(T, min_T=0, max_T=1000, units='K'):
        """
        Validate temperature value.

        Parameters:
        -----------
        T : float
            Temperature value
        min_T : float
            Minimum valid temperature
        max_T : float
            Maximum valid temperature
        units : str
            Temperature units ('K', 'C', 'F')

        Returns:
        --------
        dict : Validation result
        """
        result = {
            'valid': True,
            'value': T,
            'units': units,
            'warnings': [],
            'errors': []
        }

        # Convert to Kelvin for absolute checks
        if units == 'C':
            T_K = T + 273.15
        elif units == 'F':
            T_K = (T - 32) * 5/9 + 273.15
        else:
            T_K = T

        # Check absolute zero
        if T_K < 0:
            result['valid'] = False
            result['errors'].append(
                f"Temperature {T} {units} is below absolute zero"
            )
        elif T_K < 1:
            result['warnings'].append(
                f"Temperature {T} {units} is very close to absolute zero"
            )

        # Check range
        if T < min_T or T > max_T:
            result['warnings'].append(
                f"Temperature {T} {units} is outside typical range "
                f"[{min_T}, {max_T}] {units}"
            )

        return result

    @staticmethod
    def validate_pressure(P, min_P=0, max_P=1e8, allow_vacuum=False):
        """
        Validate pressure value.

        Parameters:
        -----------
        P : float
            Pressure value (Pa)
        min_P : float
            Minimum valid pressure (Pa)
        max_P : float
            Maximum valid pressure (Pa)
        allow_vacuum : bool
            Allow pressures below atmospheric

        Returns:
        --------
        dict : Validation result
        """
        result = {
            'valid': True,
            'value': P,
            'units': 'Pa',
            'warnings': [],
            'errors': []
        }

        # Check absolute pressure
        if P < 0:
            result['valid'] = False
            result['errors'].append(
                f"Pressure {P/1000:.1f} kPa is negative (absolute pressure cannot be negative)"
            )

        # Check vacuum conditions
        if not allow_vacuum and P < 101325:
            result['warnings'].append(
                f"Pressure {P/1000:.1f} kPa is below atmospheric (vacuum conditions)"
            )

        # Check range
        if P < min_P or P > max_P:
            result['warnings'].append(
                f"Pressure {P/1e6:.2f} MPa is outside typical range "
                f"[{min_P/1e6:.2f}, {max_P/1e6:.2f}] MPa"
            )

        # Check for very high pressure
        if P > 100e6:  # > 100 MPa
            result['warnings'].append(
                f"Very high pressure {P/1e6:.1f} MPa - ensure equipment rating"
            )

        return result

    @staticmethod
    def validate_flow_rate(Q, min_Q=0, max_Q=1000):
        """
        Validate flow rate.

        Parameters:
        -----------
        Q : float
            Flow rate (m³/s)
        min_Q : float
            Minimum valid flow rate
        max_Q : float
            Maximum valid flow rate

        Returns:
        --------
        dict : Validation result
        """
        result = {
            'valid': True,
            'value': Q,
            'units': 'm³/s',
            'warnings': [],
            'errors': []
        }

        if Q < 0:
            result['valid'] = False
            result['errors'].append(
                f"Flow rate {Q} m³/s is negative"
            )

        if Q == 0:
            result['warnings'].append("Flow rate is zero (no flow condition)")

        if Q < min_Q or Q > max_Q:
            result['warnings'].append(
                f"Flow rate {Q} m³/s is outside typical range "
                f"[{min_Q}, {max_Q}] m³/s"
            )

        return result

    @staticmethod
    def validate_reynolds_number(Re):
        """
        Validate and classify Reynolds number.

        Parameters:
        -----------
        Re : float
            Reynolds number

        Returns:
        --------
        dict : Validation result with flow regime
        """
        result = {
            'valid': True,
            'value': Re,
            'warnings': [],
            'errors': []
        }

        if Re < 0:
            result['valid'] = False
            result['errors'].append(
                f"Reynolds number {Re} is negative"
            )
            return result

        # Classify flow regime
        if Re < 2300:
            result['regime'] = 'Laminar'
        elif Re < 4000:
            result['regime'] = 'Transitional'
            result['warnings'].append(
                f"Re={Re:.0f} is in transition zone - flow may be unstable"
            )
        else:
            result['regime'] = 'Turbulent'

        # Check for very high Re
        if Re > 1e8:
            result['warnings'].append(
                f"Very high Reynolds number (Re={Re:.2e}) - "
                "verify correlations are valid"
            )

        return result


# Example usage
print("=== Physical Validity Checks ===")

validator = PhysicalValidator()

# Temperature check
T_check = validator.validate_temperature(25, units='C')
print(f"\nTemperature 25°C:")
print(f"  Valid: {T_check['valid']}")
if T_check['warnings']:
    for warning in T_check['warnings']:
        print(f"  Warning: {warning}")

# Pressure check
P_check = validator.validate_pressure(50000)  # 50 kPa (vacuum)
print(f"\nPressure 50 kPa:")
print(f"  Valid: {P_check['valid']}")
if P_check['warnings']:
    for warning in P_check['warnings']:
        print(f"  Warning: {warning}")

# Reynolds number check
Re_check = validator.validate_reynolds_number(3000)
print(f"\nReynolds number 3000:")
print(f"  Valid: {Re_check['valid']}")
print(f"  Regime: {Re_check['regime']}")
if Re_check['warnings']:
    for warning in Re_check['warnings']:
        print(f"  Warning: {warning}")
```

### 6.3 Unit Verification

```python
from pint import UnitRegistry

ureg = UnitRegistry()

class UnitValidator:
    """Validate units in engineering calculations."""

    def __init__(self):
        self.ureg = UnitRegistry()

    def check_unit_compatibility(self, value1, unit1, value2, unit2):
        """
        Check if two units are compatible (same dimensionality).

        Parameters:
        -----------
        value1 : float
            First value
        unit1 : str
            First unit
        value2 : float
            Second value
        unit2 : str
            Second unit

        Returns:
        --------
        dict : Compatibility check results
        """
        try:
            q1 = value1 * self.ureg(unit1)
            q2 = value2 * self.ureg(unit2)

            compatible = q1.dimensionality == q2.dimensionality

            result = {
                'compatible': compatible,
                'value1': f"{value1} {unit1}",
                'value2': f"{value2} {unit2}",
                'dimensionality1': str(q1.dimensionality),
                'dimensionality2': str(q2.dimensionality)
            }

            if compatible:
                # Convert to same units for comparison
                q2_converted = q2.to(q1.units)
                result['value2_converted'] = f"{q2_converted.magnitude:.4f} {unit1}"
            else:
                result['error'] = (
                    f"Incompatible dimensions: {q1.dimensionality} vs {q2.dimensionality}"
                )

            return result

        except Exception as e:
            return {
                'compatible': False,
                'error': str(e)
            }

    def validate_equation_units(self, lhs_value, lhs_unit, rhs_value, rhs_unit):
        """
        Validate dimensional consistency of an equation.

        Parameters:
        -----------
        lhs_value : float
            Left-hand side value
        lhs_unit : str
            Left-hand side unit
        rhs_value : float
            Right-hand side value
        rhs_unit : str
            Right-hand side unit

        Returns:
        --------
        dict : Validation results
        """
        compatibility = self.check_unit_compatibility(
            lhs_value, lhs_unit, rhs_value, rhs_unit
        )

        result = {
            'dimensionally_consistent': compatibility['compatible'],
            'lhs': f"{lhs_value} {lhs_unit}",
            'rhs': f"{rhs_value} {rhs_unit}"
        }

        if compatibility['compatible']:
            # Check if values are approximately equal
            lhs = lhs_value * self.ureg(lhs_unit)
            rhs = rhs_value * self.ureg(rhs_unit)
            rhs_converted = rhs.to(lhs.units)

            relative_error = abs(lhs.magnitude - rhs_converted.magnitude) / abs(lhs.magnitude)
            result['relative_error'] = relative_error
            result['numerically_consistent'] = relative_error < 1e-6

            if not result['numerically_consistent']:
                result['warning'] = (
                    f"Dimensionally consistent but numerically different: "
                    f"{relative_error*100:.2f}% error"
                )
        else:
            result['error'] = compatibility.get('error', 'Incompatible dimensions')

        return result

    def suggest_unit_conversion(self, value, from_unit, to_unit):
        """
        Suggest unit conversion and check validity.

        Parameters:
        -----------
        value : float
            Value to convert
        from_unit : str
            Source unit
        to_unit : str
            Target unit

        Returns:
        --------
        dict : Conversion results
        """
        try:
            quantity = value * self.ureg(from_unit)
            converted = quantity.to(self.ureg(to_unit))

            return {
                'success': True,
                'original': f"{value} {from_unit}",
                'converted': f"{converted.magnitude:.6g} {to_unit}",
                'conversion_factor': converted.magnitude / value
            }
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'message': f"Cannot convert {from_unit} to {to_unit}"
            }


# Example usage
print("=== Unit Validation ===")

validator = UnitValidator()

# Check unit compatibility
print("\n--- Unit Compatibility ---")
check1 = validator.check_unit_compatibility(100, 'psi', 6.895, 'bar')
print(f"100 psi vs 6.895 bar:")
print(f"  Compatible: {check1['compatible']}")
if check1['compatible']:
    print(f"  Converted: {check1['value2_converted']}")

# Check incompatible units
check2 = validator.check_unit_compatibility(100, 'm/s', 50, 'kg')
print(f"\n100 m/s vs 50 kg:")
print(f"  Compatible: {check2['compatible']}")
if not check2['compatible']:
    print(f"  Error: {check2['error']}")

# Validate equation
print("\n--- Equation Validation ---")
# Bernoulli: P1/ρ + v1²/2 = P2/ρ + v2²/2 (simplified)
eq_check = validator.validate_equation_units(
    1000, 'Pa',
    1000, 'J/m³'  # Pa = J/m³
)
print("1000 Pa = 1000 J/m³:")
print(f"  Dimensionally consistent: {eq_check['dimensionally_consistent']}")
print(f"  Numerically consistent: {eq_check['numerically_consistent']}")

# Suggest conversion
print("\n--- Unit Conversion ---")
conversion = validator.suggest_unit_conversion(100, 'gallon/minute', 'liter/second')
print(f"Convert 100 gpm to L/s:")
if conversion['success']:
    print(f"  Result: {conversion['converted']}")
    print(f"  Factor: {conversion['conversion_factor']:.6f}")
```

---

## Summary

This document provides a comprehensive collection of engineering helper functions and utilities covering:

1. **Unit Conversions**: Using pint library and custom converters for engineering units
2. **Engineering Constants**: Physical constants, standard conditions, and material properties
3. **Quick Calculators**: Reynolds number, friction factors, NPSH, specific speed, and affinity laws
4. **Resource Listers**: Tools to check available packages, databases, and system software
5. **Context Initialization**: Session setup with common imports and default configurations
6. **Error Handling**: Numerical issue detection, physical validity checks, and unit verification

These utilities form a foundation for engineering calculations and can be extended based on specific application needs. They promote:

- **Consistency**: Standardized approaches to common calculations
- **Reliability**: Built-in error checking and validation
- **Efficiency**: Quick access to frequently used functions
- **Maintainability**: Clear, documented code with examples

For specific applications, these helpers can be combined with specialized engineering libraries such as:
- **CoolProp**: Thermodynamic and transport properties
- **fluids**: Fluid mechanics calculations
- **ht**: Heat transfer calculations
- **chemicals**: Chemical property databases
- **thermo**: Chemical thermodynamics
