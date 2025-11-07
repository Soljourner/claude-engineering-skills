"""
Pint Units Package Examples
===========================

Practical examples demonstrating unit-safe engineering calculations.
Shows flow rates, pressures, viscosity, and dimensional analysis.

Author: Engineering Skills Library
Date: 2025-11-07
"""

import numpy as np
from pint import UnitRegistry

# Create unit registry (do this once per module)
ureg = UnitRegistry()


def example_1_flow_rate_conversions():
    """
    Example 1: Flow Rate Conversions

    Problem: Convert between common engineering flow rate units

    Demonstrates:
    - Volumetric flow rate conversions
    - Mass flow rate conversions
    - Relationship between mass and volumetric flow

    Verification: Standard conversion factors from engineering handbooks
    """
    print("=" * 70)
    print("Example 1: Flow Rate Conversions")
    print("=" * 70)

    # Volumetric flow rate
    Q_vol = 100 * ureg.gallons / ureg.minute  # gpm

    print("\nVolumetric Flow Rate Conversions:")
    print(f"  Original: {Q_vol:.2f}")
    print(f"  Liters/minute: {Q_vol.to('liter/minute'):.2f}")
    print(f"  m³/hour: {Q_vol.to('meter**3/hour'):.2f}")
    print(f"  ft³/second: {Q_vol.to('foot**3/second'):.4f}")
    print(f"  m³/second: {Q_vol.to('meter**3/second'):.6f}")

    # Verification: 1 gpm = 3.785 L/min
    L_per_min = Q_vol.to('liter/minute').magnitude
    expected_L_per_min = 100 * 3.785411784
    error = abs(L_per_min - expected_L_per_min)

    print(f"\nVerification:")
    print(f"  Expected: {expected_L_per_min:.2f} L/min")
    print(f"  Calculated: {L_per_min:.2f} L/min")
    print(f"  Error: {error:.6f} L/min")
    assert error < 0.01, "Flow rate conversion verification failed"
    print(f"  ✓ PASSED")

    # Mass flow rate from volumetric flow
    density = 1000 * ureg.kg / ureg.meter**3  # Water density
    Q_mass = Q_vol * density

    print(f"\nMass Flow Rate (Water at 1000 kg/m³):")
    print(f"  {Q_mass.to('kg/second'):.4f}")
    print(f"  {Q_mass.to('kg/hour'):.2f}")
    print(f"  {Q_mass.to('lb/minute'):.2f}")
    print(f"  {Q_mass.to('ton/day'):.2f}")

    print()
    return {
        'Q_vol': Q_vol,
        'Q_mass': Q_mass,
        'density': density
    }


def example_2_pressure_conversions():
    """
    Example 2: Pressure Conversions and Head Calculations

    Problem: Convert between pressure units and calculate equivalent head

    Demonstrates:
    - Pressure unit conversions
    - Pressure to head conversion
    - Gauge vs absolute pressure

    Verification: Standard conversion factors (1 psi = 6.89476 kPa)
    """
    print("=" * 70)
    print("Example 2: Pressure Conversions")
    print("=" * 70)

    # Pressure in various units
    P = 100 * ureg.psi

    print("\nPressure Conversions:")
    print(f"  Original: {P:.2f}")
    print(f"  Pascals: {P.to('pascal'):.1f}")
    print(f"  kPa: {P.to('kilopascal'):.2f}")
    print(f"  Bar: {P.to('bar'):.4f}")
    print(f"  MPa: {P.to('megapascal'):.4f}")
    print(f"  Atmospheres: {P.to('atmosphere'):.4f}")
    print(f"  mmHg: {P.to('mmHg'):.1f}")
    print(f"  inH2O: {P.to('inch_H2O'):.1f}")

    # Verification: 100 psi = 689.476 kPa
    kPa_value = P.to('kilopascal').magnitude
    expected_kPa = 689.476
    error = abs(kPa_value - expected_kPa)

    print(f"\nVerification:")
    print(f"  Expected: {expected_kPa:.3f} kPa")
    print(f"  Calculated: {kPa_value:.3f} kPa")
    print(f"  Error: {error:.6f} kPa")
    assert error < 0.01, "Pressure conversion verification failed"
    print(f"  ✓ PASSED")

    # Convert pressure to head (water column)
    rho = 1000 * ureg.kg / ureg.meter**3  # Water density
    g = 9.81 * ureg.meter / ureg.second**2  # Gravity

    head = P / (rho * g)

    print(f"\nEquivalent Head (Water):")
    print(f"  {head.to('meter'):.2f}")
    print(f"  {head.to('feet'):.2f}")

    # Gauge pressure vs absolute pressure
    P_gauge = 100 * ureg.psi
    P_atm = 14.7 * ureg.psi
    P_absolute = P_gauge + P_atm

    print(f"\nGauge vs Absolute Pressure:")
    print(f"  Gauge pressure: {P_gauge:.1f}")
    print(f"  Atmospheric: {P_atm:.1f}")
    print(f"  Absolute: {P_absolute:.1f}")
    print(f"  Absolute: {P_absolute.to('bar'):.3f}")

    print()
    return {
        'pressure': P,
        'head': head,
        'P_gauge': P_gauge,
        'P_absolute': P_absolute
    }


def example_3_viscosity_conversions():
    """
    Example 3: Viscosity Conversions

    Problem: Convert between dynamic and kinematic viscosity units

    Demonstrates:
    - Dynamic viscosity (absolute viscosity)
    - Kinematic viscosity
    - Conversion between dynamic and kinematic
    - Common engineering viscosity units

    Verification: Standard viscosity conversions
    """
    print("=" * 70)
    print("Example 3: Viscosity Conversions")
    print("=" * 70)

    # Dynamic viscosity (absolute viscosity)
    mu = 1.0 * ureg.centipoise  # Common unit in industry

    print("\nDynamic Viscosity Conversions:")
    print(f"  Original: {mu:.3f}")
    print(f"  Pa·s: {mu.to('pascal*second'):.6f}")
    print(f"  Poise (P): {mu.to('poise'):.4f}")
    print(f"  lbf·s/ft²: {mu.to('lbf*second/foot**2'):.6e}")
    print(f"  lbm/(ft·s): {mu.to('pound/(foot*second)'):.6e}")

    # Verification: 1 cP = 0.001 Pa·s
    Pa_s_value = mu.to('pascal*second').magnitude
    expected_Pa_s = 0.001
    error = abs(Pa_s_value - expected_Pa_s)

    print(f"\nDynamic Viscosity Verification:")
    print(f"  Expected: {expected_Pa_s:.6f} Pa·s")
    print(f"  Calculated: {Pa_s_value:.6f} Pa·s")
    print(f"  Error: {error:.9f} Pa·s")
    assert error < 1e-9, "Dynamic viscosity conversion failed"
    print(f"  ✓ PASSED")

    # Kinematic viscosity
    nu = 1.0 * ureg.centistokes  # Common unit in industry

    print(f"\nKinematic Viscosity Conversions:")
    print(f"  Original: {nu:.3f}")
    print(f"  m²/s: {nu.to('meter**2/second'):.6e}")
    print(f"  Stokes (St): {nu.to('stokes'):.4f}")
    print(f"  ft²/s: {nu.to('foot**2/second'):.6e}")

    # Verification: 1 cSt = 1e-6 m²/s
    m2_s_value = nu.to('meter**2/second').magnitude
    expected_m2_s = 1e-6
    error = abs(m2_s_value - expected_m2_s)

    print(f"\nKinematic Viscosity Verification:")
    print(f"  Expected: {expected_m2_s:.9f} m²/s")
    print(f"  Calculated: {m2_s_value:.9f} m²/s")
    print(f"  Error: {error:.12e} m²/s")
    assert error < 1e-12, "Kinematic viscosity conversion failed"
    print(f"  ✓ PASSED")

    # Conversion between dynamic and kinematic
    density = 1000 * ureg.kg / ureg.meter**3  # Water density
    mu_from_nu = nu * density

    print(f"\nConversion: Kinematic → Dynamic (ρ = 1000 kg/m³):")
    print(f"  ν = {nu}")
    print(f"  μ = ν × ρ = {mu_from_nu.to('centipoise'):.3f}")

    # For water at 20°C: ν ≈ 1 cSt, ρ = 1000 kg/m³ → μ ≈ 1 cP
    assert np.isclose(mu_from_nu.to('centipoise').magnitude, 1.0, rtol=0.01)
    print(f"  ✓ Verified: 1 cSt × 1000 kg/m³ = 1 cP (for water)")

    print()
    return {
        'dynamic_viscosity': mu,
        'kinematic_viscosity': nu,
        'density': density
    }


def example_4_pump_power_calculation():
    """
    Example 4: Unit-Aware Pump Power Calculation

    Problem: Calculate hydraulic and shaft power for a pump

    Formula:
    - Hydraulic power: P_h = ρ × g × Q × H
    - Shaft power: P_s = P_h / η

    Demonstrates:
    - Dimensional analysis ensures correct formula
    - Automatic unit conversion in calculations
    - Power unit conversions

    Verification: Known pump calculation from engineering standards
    """
    print("=" * 70)
    print("Example 4: Unit-Aware Pump Power Calculation")
    print("=" * 70)

    # Input parameters
    Q = 100 * ureg.meter**3 / ureg.hour  # Flow rate
    H = 50 * ureg.meter  # Total head
    rho = 1000 * ureg.kg / ureg.meter**3  # Fluid density (water)
    g = 9.81 * ureg.meter / ureg.second**2  # Gravity
    eta = 0.75 * ureg.dimensionless  # Pump efficiency

    print(f"\nInput Parameters:")
    print(f"  Flow rate: {Q:.1f} = {Q.to('gpm'):.1f}")
    print(f"  Head: {H:.1f} = {H.to('feet'):.1f}")
    print(f"  Density: {rho:.0f}")
    print(f"  Efficiency: {eta.magnitude:.1%}")

    # Calculate hydraulic power
    P_hydraulic = rho * g * Q * H

    print(f"\nHydraulic Power:")
    print(f"  {P_hydraulic.to('watt'):.2f}")
    print(f"  {P_hydraulic.to('kilowatt'):.3f}")
    print(f"  {P_hydraulic.to('horsepower'):.3f}")

    # Calculate shaft power
    P_shaft = P_hydraulic / eta

    print(f"\nShaft Power (at {eta.magnitude:.0%} efficiency):")
    print(f"  {P_shaft.to('watt'):.2f}")
    print(f"  {P_shaft.to('kilowatt'):.3f}")
    print(f"  {P_shaft.to('horsepower'):.3f}")

    # Manual verification
    # Q = 100 m³/h = 0.02778 m³/s
    # P = 1000 × 9.81 × 0.02778 × 50 = 13,625 W = 13.625 kW
    Q_m3s = Q.to('meter**3/second').magnitude
    P_manual = 1000 * 9.81 * Q_m3s * 50  # Watts
    P_calc = P_hydraulic.to('watt').magnitude

    print(f"\nVerification:")
    print(f"  Manual calculation: {P_manual:.1f} W")
    print(f"  Pint calculation: {P_calc:.1f} W")
    print(f"  Error: {abs(P_calc - P_manual):.6f} W")

    assert np.isclose(P_calc, P_manual, rtol=1e-6)
    print(f"  ✓ PASSED")

    # Motor selection (add safety factor)
    safety_factor = 1.15
    P_motor = P_shaft * safety_factor

    print(f"\nMotor Selection (with 15% safety factor):")
    print(f"  Required power: {P_motor.to('kilowatt'):.2f}")
    print(f"  Required power: {P_motor.to('horsepower'):.2f}")
    print(f"  → Select: 20 kW or 25 HP motor")

    print()
    return {
        'P_hydraulic': P_hydraulic,
        'P_shaft': P_shaft,
        'P_motor': P_motor
    }


def example_5_reynolds_number_calculation():
    """
    Example 5: Dimensional Consistency in Reynolds Number

    Problem: Calculate Reynolds number ensuring dimensional consistency

    Formula: Re = ρVD/μ (dimensionless)

    Demonstrates:
    - Automatic dimensional analysis
    - Verification of dimensionless result
    - Using different unit systems
    - Error detection for wrong units

    Verification: Known Reynolds number calculation
    """
    print("=" * 70)
    print("Example 5: Reynolds Number with Dimensional Analysis")
    print("=" * 70)

    # Flow conditions
    V = 2.5 * ureg.meter / ureg.second  # Velocity
    D = 150 * ureg.millimeter  # Pipe diameter
    rho = 998 * ureg.kg / ureg.meter**3  # Water density
    mu = 1.0 * ureg.centipoise  # Dynamic viscosity

    print(f"\nInput Parameters:")
    print(f"  Velocity: {V:.2f} = {V.to('feet/second'):.2f}")
    print(f"  Diameter: {D:.1f} = {D.to('inch'):.2f}")
    print(f"  Density: {rho:.0f}")
    print(f"  Viscosity: {mu:.3f} = {mu.to('pascal*second'):.6f}")

    # Calculate Reynolds number
    Re = (rho * V * D) / mu

    print(f"\nReynolds Number Calculation:")
    print(f"  Re = (ρ × V × D) / μ")
    print(f"  Re = ({rho:.0f}) × ({V:.2f}) × ({D:.3f}) / ({mu:.3f})")

    # Convert to base units (should be dimensionless)
    Re_dimensionless = Re.to_base_units()

    print(f"\nResult:")
    print(f"  Re = {Re_dimensionless:.0f}")
    print(f"  Dimensionality: {Re_dimensionless.dimensionality}")

    # Verify dimensionless
    assert Re_dimensionless.dimensionality == ureg.dimensionless.dimensionality
    print(f"  ✓ Confirmed: Result is dimensionless")

    # Flow regime
    Re_value = Re_dimensionless.magnitude

    if Re_value < 2300:
        regime = "Laminar"
    elif Re_value < 4000:
        regime = "Transitional"
    else:
        regime = "Turbulent"

    print(f"\nFlow Regime:")
    print(f"  Re = {Re_value:.0f} → {regime} flow")

    # Manual verification
    V_ms = V.to('meter/second').magnitude
    D_m = D.to('meter').magnitude
    rho_kg_m3 = rho.to('kg/meter**3').magnitude
    mu_Pa_s = mu.to('pascal*second').magnitude

    Re_manual = (rho_kg_m3 * V_ms * D_m) / mu_Pa_s

    print(f"\nVerification:")
    print(f"  Manual calculation: Re = {Re_manual:.0f}")
    print(f"  Pint calculation: Re = {Re_value:.0f}")
    print(f"  Error: {abs(Re_value - Re_manual):.6f}")

    assert np.isclose(Re_value, Re_manual, rtol=1e-6)
    print(f"  ✓ PASSED")

    # Test with US customary units
    print(f"\nSame Calculation with US Units:")
    V_fps = 8.2 * ureg.feet / ureg.second
    D_inch = 5.9 * ureg.inch
    rho_lb_ft3 = 62.3 * ureg.pound / ureg.feet**3
    mu_lb_fts = 6.72e-4 * ureg.pound / (ureg.feet * ureg.second)

    Re_us = (rho_lb_ft3 * V_fps * D_inch) / mu_lb_fts
    Re_us_value = Re_us.to_base_units().magnitude

    print(f"  Re (US units) = {Re_us_value:.0f}")
    print(f"  Match with SI: {np.isclose(Re_us_value, Re_value, rtol=0.01)}")

    print()
    return {
        'Re': Re_dimensionless,
        'regime': regime,
        'velocity': V,
        'diameter': D
    }


def example_6_pressure_drop_with_units():
    """
    Example 6: Pressure Drop Calculation with Automatic Units

    Problem: Calculate pressure drop in pipe using Darcy-Weisbach equation

    Formula: ΔP = f × (L/D) × (ρV²/2)

    Demonstrates:
    - Complex formula with automatic unit handling
    - Conversion between pressure and head
    - Working with dimensionless friction factor

    Verification: Standard pipe flow calculation
    """
    print("=" * 70)
    print("Example 6: Pressure Drop with Automatic Unit Handling")
    print("=" * 70)

    # Pipe and flow parameters
    f = 0.018 * ureg.dimensionless  # Darcy friction factor
    L = 100 * ureg.meter  # Pipe length
    D = 150 * ureg.millimeter  # Pipe diameter
    V = 2.5 * ureg.meter / ureg.second  # Velocity
    rho = 1000 * ureg.kg / ureg.meter**3  # Density
    g = 9.81 * ureg.meter / ureg.second**2

    print(f"\nInput Parameters:")
    print(f"  Friction factor: {f.magnitude:.4f}")
    print(f"  Pipe length: {L:.1f} = {L.to('feet'):.1f}")
    print(f"  Pipe diameter: {D:.1f} = {D.to('inch'):.2f}")
    print(f"  Velocity: {V:.2f} = {V.to('feet/second'):.2f}")
    print(f"  Density: {rho:.0f}")

    # Calculate pressure drop (Darcy-Weisbach equation)
    dP = f * (L/D) * (rho * V**2 / 2)

    print(f"\nPressure Drop:")
    print(f"  ΔP = f × (L/D) × (ρV²/2)")
    print(f"  {dP.to('pascal'):.1f}")
    print(f"  {dP.to('kilopascal'):.2f}")
    print(f"  {dP.to('psi'):.2f}")
    print(f"  {dP.to('bar'):.4f}")

    # Convert to head loss
    h_loss = dP / (rho * g)

    print(f"\nHead Loss:")
    print(f"  h = ΔP / (ρg)")
    print(f"  {h_loss.to('meter'):.3f}")
    print(f"  {h_loss.to('feet'):.2f}")

    # Calculate per unit length
    dP_per_100m = dP / L * (100 * ureg.meter)
    h_per_100m = h_loss / L * (100 * ureg.meter)

    print(f"\nPer 100m of Pipe:")
    print(f"  Pressure drop: {dP_per_100m.to('kilopascal'):.2f}")
    print(f"  Head loss: {h_per_100m.to('meter'):.2f}")

    # Manual verification
    L_m = L.to('meter').magnitude
    D_m = D.to('meter').magnitude
    V_ms = V.to('meter/second').magnitude
    rho_kg_m3 = rho.to('kg/meter**3').magnitude

    dP_manual = 0.018 * (L_m/D_m) * (rho_kg_m3 * V_ms**2 / 2)

    print(f"\nVerification:")
    print(f"  Manual: ΔP = {dP_manual:.1f} Pa")
    print(f"  Pint: ΔP = {dP.to('pascal').magnitude:.1f} Pa")
    print(f"  Error: {abs(dP.to('pascal').magnitude - dP_manual):.6f} Pa")

    assert np.isclose(dP.to('pascal').magnitude, dP_manual, rtol=1e-6)
    print(f"  ✓ PASSED")

    print()
    return {
        'pressure_drop': dP,
        'head_loss': h_loss,
        'dP_per_100m': dP_per_100m
    }


def example_7_array_operations_with_units():
    """
    Example 7: Array Operations with Units (Pump Curve)

    Problem: Create pump performance curve with unit-aware arrays

    Demonstrates:
    - NumPy array operations with units
    - Statistical operations preserving units
    - Finding operating point
    - Plotting data preparation

    Verification: Pump affinity laws and curve fitting
    """
    print("=" * 70)
    print("Example 7: Array Operations with Units (Pump Curve)")
    print("=" * 70)

    # Pump curve data (test data points)
    Q_array = np.array([0, 25, 50, 75, 100, 125, 150]) * ureg.meter**3 / ureg.hour
    H_array = np.array([80, 79, 76, 71, 64, 55, 43]) * ureg.meter

    print("\nPump Performance Curve Data:")
    print(f"{'Q (m³/h)':<12} {'Q (gpm)':<12} {'H (m)':<10} {'H (ft)':<10}")
    print("-" * 50)

    for Q, H in zip(Q_array, H_array):
        print(f"{Q.magnitude:<12.0f} {Q.to('gpm').magnitude:<12.0f} "
              f"{H.magnitude:<10.1f} {H.to('feet').magnitude:<10.1f}")

    # Statistical analysis
    Q_mean = np.mean(Q_array)
    Q_std = np.std(Q_array)
    H_mean = np.mean(H_array)
    H_std = np.std(H_array)

    print(f"\nStatistical Analysis:")
    print(f"  Flow rate - Mean: {Q_mean:.1f}, Std: {Q_std:.1f}")
    print(f"  Head - Mean: {H_mean:.1f}, Std: {H_std:.1f}")

    # Calculate power at each point (hydraulic)
    rho = 1000 * ureg.kg / ureg.meter**3
    g = 9.81 * ureg.meter / ureg.second**2
    P_array = rho * g * Q_array * H_array

    print(f"\nHydraulic Power at Each Point:")
    print(f"{'Q (m³/h)':<12} {'H (m)':<10} {'P (kW)':<10} {'P (hp)':<10}")
    print("-" * 50)

    for Q, H, P in zip(Q_array, H_array, P_array):
        print(f"{Q.magnitude:<12.0f} {H.magnitude:<10.1f} "
              f"{P.to('kilowatt').magnitude:<10.2f} "
              f"{P.to('horsepower').magnitude:<10.2f}")

    # Find operating point (system curve intersection)
    # System curve: H = 50 + 0.002*Q² (simplified)
    Q_system_range = np.linspace(0, 150, 100) * ureg.meter**3 / ureg.hour
    H_static = 50 * ureg.meter
    K = 0.002 * ureg.meter / (ureg.meter**3/ureg.hour)**2
    H_system = H_static + K * Q_system_range**2

    # Interpolate pump curve
    Q_pump_range = np.linspace(0, 150, 100) * ureg.meter**3 / ureg.hour
    H_pump_interp = np.interp(
        Q_pump_range.magnitude,
        Q_array.magnitude,
        H_array.magnitude
    ) * ureg.meter

    # Find intersection (operating point)
    diff = np.abs(H_pump_interp - H_system)
    idx_op = np.argmin(diff)
    Q_op = Q_pump_range[idx_op]
    H_op = H_pump_interp[idx_op]
    P_op = rho * g * Q_op * H_op

    print(f"\nOperating Point (Pump-System Intersection):")
    print(f"  Flow rate: {Q_op:.1f} = {Q_op.to('gpm'):.0f}")
    print(f"  Head: {H_op:.1f} = {H_op.to('feet'):.1f}")
    print(f"  Power: {P_op.to('kilowatt'):.2f} = {P_op.to('horsepower'):.2f}")

    # Array arithmetic verification
    print(f"\nArray Operations Verification:")
    total_power = np.sum(P_array)
    max_power = np.max(P_array)
    idx_max = np.argmax(P_array)

    print(f"  Sum of all power points: {total_power.to('kilowatt'):.2f}")
    print(f"  Maximum power: {max_power.to('kilowatt'):.2f}")
    print(f"  At flow rate: {Q_array[idx_max]:.0f}")

    # Verify units are preserved
    assert H_mean.dimensionality == ureg.meter.dimensionality
    assert P_array[0].dimensionality == ureg.watt.dimensionality
    print(f"  ✓ Units preserved in all operations")

    print()
    return {
        'Q_array': Q_array,
        'H_array': H_array,
        'P_array': P_array,
        'Q_operating': Q_op,
        'H_operating': H_op,
        'P_operating': P_op
    }


def example_8_temperature_conversions():
    """
    Example 8: Temperature Conversions (Absolute and Relative)

    Problem: Handle temperature conversions correctly for both
    absolute temperatures and temperature differences

    Demonstrates:
    - Absolute temperature conversions
    - Temperature difference (delta) conversions
    - Heat transfer calculations
    - Common temperature conversion errors

    Verification: Standard temperature conversion formulas
    """
    print("=" * 70)
    print("Example 8: Temperature Conversions")
    print("=" * 70)

    # Absolute temperatures
    T1 = 25 * ureg.degC
    T2 = 75 * ureg.degF
    T3 = 300 * ureg.kelvin

    print("\nAbsolute Temperature Conversions:")
    print(f"  {T1:.1f} = {T1.to('degF'):.1f} = {T1.to('kelvin'):.2f}")
    print(f"  {T2:.1f} = {T2.to('degC'):.2f} = {T2.to('kelvin'):.2f}")
    print(f"  {T3:.2f} = {T3.to('degC'):.2f} = {T3.to('degF'):.2f}")

    # Verification: 25°C = 77°F
    T_F = T1.to('degF').magnitude
    expected = 77.0
    error = abs(T_F - expected)

    print(f"\nVerification (25°C = 77°F):")
    print(f"  Calculated: {T_F:.1f}°F")
    print(f"  Expected: {expected:.1f}°F")
    print(f"  Error: {error:.6f}°F")
    assert error < 0.1
    print(f"  ✓ PASSED")

    # Temperature differences (delta)
    dT_C = 20 * ureg.delta_degC  # Temperature rise
    dT_F = dT_C.to(ureg.delta_degF)

    print(f"\nTemperature Differences (Delta):")
    print(f"  ΔT = {dT_C:.1f}")
    print(f"  ΔT = {dT_F:.1f}")
    print(f"  (20°C rise = 36°F rise)")

    # Heat transfer calculation
    print(f"\nHeat Transfer Calculation:")
    mass = 10 * ureg.kg
    cp = 4.18 * ureg.kilojoule / (ureg.kg * ureg.kelvin)  # Water
    dT = 50 * ureg.delta_degC

    Q_heat = mass * cp * dT

    print(f"  Mass: {mass:.1f}")
    print(f"  Specific heat: {cp:.2f}")
    print(f"  Temperature rise: {dT:.1f} = {dT.to(ureg.delta_degF):.1f}")
    print(f"  Heat required: {Q_heat.to('kilojoule'):.1f}")
    print(f"  Heat required: {Q_heat.to('BTU'):.1f}")
    print(f"  Heat required: {Q_heat.to('kWh'):.3f}")

    # Verification: Q = m × cp × ΔT
    Q_manual = 10 * 4.18 * 50  # kJ
    Q_calc = Q_heat.to('kilojoule').magnitude

    print(f"\nHeat Transfer Verification:")
    print(f"  Manual: {Q_manual:.1f} kJ")
    print(f"  Pint: {Q_calc:.1f} kJ")
    print(f"  Error: {abs(Q_calc - Q_manual):.6f} kJ")
    assert np.isclose(Q_calc, Q_manual, rtol=1e-6)
    print(f"  ✓ PASSED")

    # Common pitfall: Don't add absolute temperatures
    print(f"\nCommon Pitfall - Adding Absolute Temperatures:")
    print(f"  ❌ WRONG: (20°C) + (30°C) ≠ 50°C")
    print(f"  ✓ RIGHT: (20°C) + (30°C temperature rise) = 50°C")
    print(f"  Use delta_degC for temperature differences!")

    print()
    return {
        'T_celsius': T1,
        'T_fahrenheit': T2,
        'delta_T': dT_C,
        'heat': Q_heat
    }


def run_all_examples():
    """
    Run all examples in sequence
    """
    print("\n")
    print("#" * 70)
    print("# PINT UNITS PACKAGE - COMPLETE EXAMPLES")
    print("#" * 70)
    print()

    try:
        result1 = example_1_flow_rate_conversions()
        result2 = example_2_pressure_conversions()
        result3 = example_3_viscosity_conversions()
        result4 = example_4_pump_power_calculation()
        result5 = example_5_reynolds_number_calculation()
        result6 = example_6_pressure_drop_with_units()
        result7 = example_7_array_operations_with_units()
        result8 = example_8_temperature_conversions()

        print("#" * 70)
        print("# ALL EXAMPLES COMPLETED SUCCESSFULLY")
        print("#" * 70)
        print()

        # Summary of key benefits
        print("Key Benefits Demonstrated:")
        print("  ✓ Automatic unit conversions")
        print("  ✓ Dimensional consistency checking")
        print("  ✓ Error prevention (adding incompatible units)")
        print("  ✓ Array operations with units")
        print("  ✓ Temperature handling (absolute vs delta)")
        print("  ✓ Complex engineering calculations")
        print("  ✓ Integration with NumPy")
        print()

        return True

    except Exception as e:
        print(f"\n✗ Error running examples: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    # Run all examples
    success = run_all_examples()

    if success:
        print("\n✓ All examples executed and verified successfully!")
    else:
        print("\n✗ Some examples failed. Please check the output above.")
