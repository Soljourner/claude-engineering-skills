"""
Engineering Context Initialization Script
==========================================

This script sets up a standard engineering environment with:
- Physical constants
- Unit registry
- Common imports
- Plotting configuration
- Precision settings

Usage:
    exec(open('init-script.py').read())

Or copy-paste the entire contents into your Python session.
"""

import sys
import warnings

# ============================================================================
# CORE IMPORTS
# ============================================================================

print("=" * 70)
print("Initializing Engineering Context...")
print("=" * 70)

# NumPy - Essential for numerical computing
try:
    import numpy as np
    print("✓ NumPy imported as 'np'")
except ImportError:
    print("✗ NumPy not available - some features will not work")
    np = None

# Pint - Unit registry for dimensional analysis
try:
    from pint import UnitRegistry
    ureg = UnitRegistry()
    print("✓ Pint imported - Unit registry available as 'ureg'")
except ImportError:
    print("✗ Pint not available - unit conversions will not work")
    ureg = None

# ============================================================================
# SCIENTIFIC COMPUTING LIBRARIES
# ============================================================================

# SciPy - Scientific computing functions
try:
    import scipy
    from scipy import optimize, integrate, interpolate, special
    print("✓ SciPy imported (optimize, integrate, interpolate, special)")
except ImportError:
    print("✗ SciPy not available")
    scipy = None

# Matplotlib - Plotting and visualization
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    print("✓ Matplotlib imported as 'plt'")
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    print("✗ Matplotlib not available - plotting will not work")
    plt = None
    mpl = None
    MATPLOTLIB_AVAILABLE = False

# Pandas - Data manipulation (optional)
try:
    import pandas as pd
    print("✓ Pandas imported as 'pd'")
except ImportError:
    print("  Pandas not available (optional)")
    pd = None

print()

# ============================================================================
# STANDARD PHYSICAL CONSTANTS
# ============================================================================

print("Setting up physical constants...")

# Gravitational acceleration (m/s²)
g = 9.80665

# Universal gas constant (J/(mol·K))
R = 8.314462618

# Specific gas constant for air (J/(kg·K))
R_air = 287.05

# Specific heat ratio for air (dimensionless)
gamma_air = 1.4

# Specific heat at constant pressure for air (J/(kg·K))
cp_air = 1005.0

# Specific heat at constant volume for air (J/(kg·K))
cv_air = 718.0

# Standard atmospheric pressure (Pa)
P_atm = 101325.0

# Standard temperature (K) - 15°C
T_std = 288.15

# Standard air density (kg/m³)
rho_air_std = 1.225

# Dynamic viscosity of air at standard conditions (Pa·s)
mu_air_std = 1.789e-5

# Kinematic viscosity of air at standard conditions (m²/s)
nu_air_std = 1.460e-5

# Speed of light in vacuum (m/s)
c_light = 299792458.0

# Stefan-Boltzmann constant (W/(m²·K⁴))
sigma_sb = 5.670374419e-8

# Boltzmann constant (J/K)
k_B = 1.380649e-23

# Avogadro's number (1/mol)
N_A = 6.02214076e23

# Planck constant (J·s)
h_planck = 6.62607015e-34

# Pi (from numpy if available, otherwise use math)
if np is not None:
    pi = np.pi
else:
    import math
    pi = math.pi

print("✓ Physical constants defined")
print()

# ============================================================================
# NUMPY CONFIGURATION
# ============================================================================

if np is not None:
    print("Configuring NumPy...")

    # Set print options for better readability
    np.set_printoptions(
        precision=6,      # 6 decimal places
        suppress=True,    # Suppress scientific notation for small numbers
        linewidth=100,    # Line width for array printing
        threshold=1000    # Threshold for summarization
    )

    print("✓ NumPy print precision set to 6 decimals")
    print()

# ============================================================================
# MATPLOTLIB CONFIGURATION
# ============================================================================

if MATPLOTLIB_AVAILABLE:
    print("Configuring Matplotlib...")

    # Set default figure size
    mpl.rcParams['figure.figsize'] = (10, 6)

    # Enable grid by default
    mpl.rcParams['axes.grid'] = True
    mpl.rcParams['grid.alpha'] = 0.3

    # Font settings
    mpl.rcParams['font.size'] = 11
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 14
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['legend.fontsize'] = 10

    # Line settings
    mpl.rcParams['lines.linewidth'] = 1.5
    mpl.rcParams['lines.markersize'] = 6

    # Try to use LaTeX rendering if available
    try:
        mpl.rcParams['text.usetex'] = False  # Set to True if you have LaTeX installed
        mpl.rcParams['font.family'] = 'serif'
    except:
        pass

    # Color cycle - professional colors
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=[
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ])

    print("✓ Matplotlib style configured")
    print("  - Figure size: 10x6 inches")
    print("  - Grid enabled with alpha=0.3")
    print("  - Professional color scheme applied")
    print()

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

print("Defining helper functions...")

def show_constants():
    """Display all defined constants"""
    print("\nStandard Physical Constants:")
    print(f"  g           = {g} m/s²")
    print(f"  R           = {R} J/(mol·K)")
    print(f"  R_air       = {R_air} J/(kg·K)")
    print(f"  gamma_air   = {gamma_air}")
    print(f"  cp_air      = {cp_air} J/(kg·K)")
    print(f"  cv_air      = {cv_air} J/(kg·K)")
    print(f"  P_atm       = {P_atm} Pa")
    print(f"  T_std       = {T_std} K ({T_std-273.15}°C)")
    print(f"  rho_air_std = {rho_air_std} kg/m³")
    print(f"  mu_air_std  = {mu_air_std} Pa·s")
    print(f"  nu_air_std  = {nu_air_std} m²/s")
    print(f"  c_light     = {c_light} m/s")
    print(f"  sigma_sb    = {sigma_sb} W/(m²·K⁴)")
    print(f"  k_B         = {k_B} J/K")
    print(f"  N_A         = {N_A} 1/mol")
    print(f"  h_planck    = {h_planck} J·s")
    print(f"  pi          = {pi}")

def atm_properties(altitude_m):
    """
    Calculate atmospheric properties at a given altitude using ISA model.

    Parameters:
        altitude_m: Altitude in meters (valid up to ~11000m for troposphere)

    Returns:
        dict with Temperature (K), Pressure (Pa), and Density (kg/m³)
    """
    # Temperature lapse rate (K/m)
    lapse_rate = 0.0065

    # Calculate temperature
    T = T_std - lapse_rate * altitude_m

    # Calculate pressure
    P = P_atm * (T / T_std) ** (g / (lapse_rate * R_air))

    # Calculate density
    rho = P / (R_air * T)

    return {
        'Temperature_K': T,
        'Temperature_C': T - 273.15,
        'Pressure_Pa': P,
        'Pressure_kPa': P / 1000,
        'Density_kg_m3': rho
    }

def reynolds_number(velocity, length, fluid='air'):
    """
    Calculate Reynolds number.

    Parameters:
        velocity: Flow velocity (m/s)
        length: Characteristic length (m)
        fluid: 'air' or provide kinematic viscosity (m²/s)

    Returns:
        Reynolds number (dimensionless)
    """
    if fluid == 'air':
        nu = nu_air_std
    else:
        nu = fluid

    Re = velocity * length / nu
    return Re

def mach_number(velocity, temperature=T_std, gamma=gamma_air, R_specific=R_air):
    """
    Calculate Mach number.

    Parameters:
        velocity: Flow velocity (m/s)
        temperature: Static temperature (K)
        gamma: Specific heat ratio
        R_specific: Specific gas constant (J/(kg·K))

    Returns:
        Mach number (dimensionless)
    """
    if np is not None:
        a = np.sqrt(gamma * R_specific * temperature)
    else:
        import math
        a = math.sqrt(gamma * R_specific * temperature)

    M = velocity / a
    return M

print("✓ Helper functions defined:")
print("  - show_constants(): Display all physical constants")
print("  - atm_properties(altitude_m): Calculate atmospheric properties")
print("  - reynolds_number(velocity, length, fluid='air'): Calculate Re")
print("  - mach_number(velocity, temperature, gamma, R_specific): Calculate M")
print()

# ============================================================================
# SUMMARY
# ============================================================================

print("=" * 70)
print("Engineering Context Initialized Successfully!")
print("=" * 70)
print("\nAvailable imports:")
print("  np       - NumPy" + (" ✓" if np is not None else " ✗"))
print("  scipy    - SciPy" + (" ✓" if scipy is not None else " ✗"))
print("  plt      - Matplotlib" + (" ✓" if plt is not None else " ✗"))
print("  pd       - Pandas" + (" ✓" if pd is not None else " ✗"))
print("  ureg     - Pint Unit Registry" + (" ✓" if ureg is not None else " ✗"))

print("\nStandard constants available:")
print("  g, R, R_air, gamma_air, cp_air, cv_air")
print("  P_atm, T_std, rho_air_std, mu_air_std, nu_air_std")
print("  c_light, sigma_sb, k_B, N_A, h_planck, pi")

print("\nHelper functions:")
print("  show_constants(), atm_properties(), reynolds_number(), mach_number()")

print("\nTip: Run show_constants() to see all values")
print("=" * 70)
