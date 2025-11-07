"""
CoolProp Database Examples
==========================

Working examples for querying thermodynamic properties using CoolProp.
These examples demonstrate real-world engineering calculations.

Installation:
    pip install CoolProp

Author: Claude Engineering Skills Library
"""

from CoolProp.CoolProp import PropsSI, Props1SI
import numpy as np


# =============================================================================
# Example 1: Water Properties at Various Temperatures
# =============================================================================

def water_properties_vs_temperature():
    """
    Calculate water properties across a temperature range at atmospheric pressure.

    This is useful for:
    - Heat transfer calculations
    - Process design
    - Understanding temperature effects on fluid properties
    """
    print("=" * 70)
    print("Example 1: Water Properties vs Temperature at 1 atm")
    print("=" * 70)

    P_atm = 101325  # Pa (1 atmosphere)
    temperatures_C = [0, 10, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]

    print(f"\n{'T(°C)':>6} {'ρ(kg/m³)':>10} {'μ(mPa·s)':>11} {'Cp(kJ/kg·K)':>12} {'k(W/m·K)':>10} {'Pr':>8}")
    print("-" * 70)

    for T_C in temperatures_C:
        T = T_C + 273.15  # Convert to Kelvin

        try:
            density = PropsSI('D', 'T', T, 'P', P_atm, 'Water')  # kg/m³
            viscosity = PropsSI('V', 'T', T, 'P', P_atm, 'Water') * 1000  # mPa·s
            cp = PropsSI('C', 'T', T, 'P', P_atm, 'Water') / 1000  # kJ/kg·K
            conductivity = PropsSI('L', 'T', T, 'P', P_atm, 'Water')  # W/m·K
            prandtl = PropsSI('PRANDTL', 'T', T, 'P', P_atm, 'Water')

            print(f"{T_C:6.0f} {density:10.2f} {viscosity:11.4f} {cp:12.4f} {conductivity:10.4f} {prandtl:8.3f}")

        except Exception as e:
            print(f"{T_C:6.0f}  Error: {e}")

    print("\nNote: Viscosity decreases significantly with temperature")
    print("      Thermal conductivity increases slightly with temperature")
    print()


# =============================================================================
# Example 2: Refrigerant Properties and Cycle Analysis
# =============================================================================

def refrigerant_saturation_properties():
    """
    Calculate saturation properties for common refrigerants.

    Essential for:
    - Refrigeration cycle design
    - HVAC system sizing
    - Heat pump analysis
    """
    print("=" * 70)
    print("Example 2: Refrigerant Saturation Properties at 25°C")
    print("=" * 70)

    T_sat = 25 + 273.15  # 25°C in Kelvin
    refrigerants = ['R134a', 'R410A', 'R32', 'R404A', 'R407C', 'R717']  # R717 = Ammonia

    print(f"\n{'Fluid':>8} {'P_sat(bar)':>12} {'ρ_liq(kg/m³)':>14} {'ρ_vap(kg/m³)':>14} {'h_fg(kJ/kg)':>13}")
    print("-" * 70)

    for fluid in refrigerants:
        try:
            # Saturation pressure
            P_sat = PropsSI('P', 'T', T_sat, 'Q', 0, fluid) / 1e5  # bar

            # Liquid properties (Q=0)
            rho_liquid = PropsSI('D', 'T', T_sat, 'Q', 0, fluid)
            h_liquid = PropsSI('H', 'T', T_sat, 'Q', 0, fluid)

            # Vapor properties (Q=1)
            rho_vapor = PropsSI('D', 'T', T_sat, 'Q', 1, fluid)
            h_vapor = PropsSI('H', 'T', T_sat, 'Q', 1, fluid)

            # Latent heat of vaporization
            h_fg = (h_vapor - h_liquid) / 1000  # kJ/kg

            print(f"{fluid:>8} {P_sat:12.2f} {rho_liquid:14.2f} {rho_vapor:14.2f} {h_fg:13.2f}")

        except Exception as e:
            print(f"{fluid:>8}  Error: {e}")

    print("\nNote: Higher saturation pressure indicates higher system pressures")
    print("      Larger h_fg means more cooling capacity per unit mass")
    print()


def refrigeration_cycle_analysis():
    """
    Analyze a complete vapor-compression refrigeration cycle.

    Cycle states:
    1. Compressor inlet (saturated vapor)
    2. Compressor outlet (superheated vapor)
    3. Condenser outlet (saturated liquid)
    4. Evaporator inlet (two-phase after expansion)
    """
    print("=" * 70)
    print("Example 3: R134a Refrigeration Cycle Analysis")
    print("=" * 70)

    fluid = 'R134a'

    # Operating conditions
    T_evap = -10 + 273.15  # Evaporator temperature: -10°C
    T_cond = 40 + 273.15   # Condenser temperature: 40°C
    eta_comp = 0.80        # Compressor isentropic efficiency

    # State 1: Saturated vapor leaving evaporator
    P1 = PropsSI('P', 'T', T_evap, 'Q', 1, fluid)
    h1 = PropsSI('H', 'T', T_evap, 'Q', 1, fluid)
    s1 = PropsSI('S', 'T', T_evap, 'Q', 1, fluid)

    print(f"\nState 1 (Evaporator outlet - saturated vapor):")
    print(f"  Temperature: {T_evap-273.15:.1f}°C")
    print(f"  Pressure: {P1/1e5:.2f} bar")
    print(f"  Enthalpy: {h1/1000:.2f} kJ/kg")
    print(f"  Entropy: {s1/1000:.4f} kJ/kg·K")

    # State 2s: Isentropic compression
    P2 = PropsSI('P', 'T', T_cond, 'Q', 0, fluid)
    h2s = PropsSI('H', 'P', P2, 'S', s1, fluid)

    # State 2: Actual compression
    h2 = h1 + (h2s - h1) / eta_comp
    T2 = PropsSI('T', 'P', P2, 'H', h2, fluid)
    s2 = PropsSI('S', 'P', P2, 'H', h2, fluid)

    print(f"\nState 2 (Compressor outlet - superheated vapor):")
    print(f"  Temperature: {T2-273.15:.1f}°C")
    print(f"  Pressure: {P2/1e5:.2f} bar")
    print(f"  Enthalpy: {h2/1000:.2f} kJ/kg")
    print(f"  Entropy: {s2/1000:.4f} kJ/kg·K")

    # State 3: Saturated liquid leaving condenser
    P3 = P2
    h3 = PropsSI('H', 'T', T_cond, 'Q', 0, fluid)
    s3 = PropsSI('S', 'T', T_cond, 'Q', 0, fluid)

    print(f"\nState 3 (Condenser outlet - saturated liquid):")
    print(f"  Temperature: {T_cond-273.15:.1f}°C")
    print(f"  Pressure: {P3/1e5:.2f} bar")
    print(f"  Enthalpy: {h3/1000:.2f} kJ/kg")
    print(f"  Entropy: {s3/1000:.4f} kJ/kg·K")

    # State 4: After expansion valve (isenthalpic)
    P4 = P1
    h4 = h3  # Isenthalpic expansion
    T4 = T_evap
    Q4 = PropsSI('Q', 'P', P4, 'H', h4, fluid)

    print(f"\nState 4 (After expansion - two-phase):")
    print(f"  Temperature: {T4-273.15:.1f}°C")
    print(f"  Pressure: {P4/1e5:.2f} bar")
    print(f"  Enthalpy: {h4/1000:.2f} kJ/kg")
    print(f"  Quality: {Q4:.4f}")

    # Performance calculations
    q_evap = h1 - h4  # Cooling capacity
    w_comp = h2 - h1  # Compressor work
    q_cond = h2 - h3  # Heat rejected
    COP = q_evap / w_comp  # Coefficient of performance

    print(f"\n--- Cycle Performance ---")
    print(f"Cooling capacity: {q_evap/1000:.2f} kJ/kg")
    print(f"Compressor work: {w_comp/1000:.2f} kJ/kg")
    print(f"Heat rejected: {q_cond/1000:.2f} kJ/kg")
    print(f"COP (Coefficient of Performance): {COP:.2f}")
    print()


# =============================================================================
# Example 4: Viscosity Calculations for Fluid Flow
# =============================================================================

def viscosity_temperature_dependence():
    """
    Calculate viscosity for different fluids across temperature range.

    Applications:
    - Pipe flow calculations
    - Pump sizing
    - Heat exchanger design
    """
    print("=" * 70)
    print("Example 4: Dynamic Viscosity vs Temperature")
    print("=" * 70)

    P = 101325  # Atmospheric pressure
    temperatures_C = np.linspace(0, 100, 11)
    fluids = {'Water': 'Water', 'Air': 'Air'}

    for fluid_name, fluid in fluids.items():
        print(f"\n{fluid_name} at 1 atm:")
        print(f"{'T(°C)':>6} {'μ(mPa·s)':>12} {'ν(mm²/s)':>12} {'Re*':>10}")
        print("-" * 45)

        for T_C in temperatures_C:
            T = T_C + 273.15

            try:
                mu = PropsSI('V', 'T', T, 'P', P, fluid) * 1000  # mPa·s
                rho = PropsSI('D', 'T', T, 'P', P, fluid)
                nu = mu / rho * 1e6  # mm²/s (kinematic viscosity)

                # Reynolds number for 1 m/s flow in 0.05m diameter pipe
                V = 1.0  # m/s
                D = 0.05  # m
                Re = rho * V * D / (mu / 1000)

                print(f"{T_C:6.0f} {mu:12.4f} {nu:12.4f} {Re:10.0f}")

            except Exception as e:
                print(f"{T_C:6.0f}  Error: {e}")

    print("\n*Reynolds number for V=1 m/s, D=0.05m")
    print("Re < 2300: Laminar flow")
    print("Re > 4000: Turbulent flow")
    print()


# =============================================================================
# Example 5: Two-Phase Properties
# =============================================================================

def two_phase_properties():
    """
    Explore two-phase region properties with varying quality.

    Quality (Q):
    - Q = 0: Saturated liquid
    - 0 < Q < 1: Two-phase mixture
    - Q = 1: Saturated vapor
    """
    print("=" * 70)
    print("Example 5: Two-Phase Properties - Quality Effects")
    print("=" * 70)

    fluid = 'Water'
    P = 1e5  # 1 bar

    # Find saturation temperature at this pressure
    T_sat = PropsSI('T', 'P', P, 'Q', 0, fluid)

    print(f"\n{fluid} at {P/1e5:.1f} bar (T_sat = {T_sat-273.15:.2f}°C)")
    print(f"\n{'Quality':>8} {'T(°C)':>8} {'ρ(kg/m³)':>12} {'h(kJ/kg)':>11} {'s(kJ/kg·K)':>13} {'v(m³/kg)':>11}")
    print("-" * 75)

    qualities = np.linspace(0, 1, 11)

    for Q in qualities:
        T = PropsSI('T', 'P', P, 'Q', Q, fluid)
        rho = PropsSI('D', 'P', P, 'Q', Q, fluid)
        h = PropsSI('H', 'P', P, 'Q', Q, fluid) / 1000  # kJ/kg
        s = PropsSI('S', 'P', P, 'Q', Q, fluid) / 1000  # kJ/kg·K
        v = 1 / rho  # Specific volume

        print(f"{Q:8.2f} {T-273.15:8.2f} {rho:12.3f} {h:11.2f} {s:13.4f} {v:11.6f}")

    print("\nNote: In two-phase region, temperature remains constant at T_sat")
    print("      Density decreases (specific volume increases) with quality")
    print("      Enthalpy and entropy increase linearly with quality")
    print()


def steam_table_example():
    """
    Generate classic steam table data using CoolProp.

    Demonstrates:
    - Saturation properties vs pressure
    - Subcooled and superheated regions
    """
    print("=" * 70)
    print("Example 6: Steam Table - Saturation Properties vs Pressure")
    print("=" * 70)

    pressures_bar = [0.01, 0.1, 1, 5, 10, 20, 50, 100, 200]

    print(f"\n{'P(bar)':>7} {'T_sat(°C)':>11} {'v_f(m³/kg)':>13} {'v_g(m³/kg)':>13} {'h_f(kJ/kg)':>12} {'h_g(kJ/kg)':>12} {'h_fg(kJ/kg)':>13}")
    print("-" * 95)

    for P_bar in pressures_bar:
        P = P_bar * 1e5  # Convert to Pa

        try:
            # Saturation temperature
            T_sat = PropsSI('T', 'P', P, 'Q', 0, 'Water')

            # Saturated liquid (f)
            rho_f = PropsSI('D', 'P', P, 'Q', 0, 'Water')
            v_f = 1 / rho_f
            h_f = PropsSI('H', 'P', P, 'Q', 0, 'Water') / 1000

            # Saturated vapor (g)
            rho_g = PropsSI('D', 'P', P, 'Q', 1, 'Water')
            v_g = 1 / rho_g
            h_g = PropsSI('H', 'P', P, 'Q', 1, 'Water') / 1000

            # Latent heat
            h_fg = h_g - h_f

            print(f"{P_bar:7.2f} {T_sat-273.15:11.2f} {v_f:13.6f} {v_g:13.4f} {h_f:12.2f} {h_g:12.2f} {h_fg:13.2f}")

        except Exception as e:
            print(f"{P_bar:7.2f}  Error: {e}")

    print("\nLegend: f=saturated liquid, g=saturated vapor, fg=evaporation")
    print()


# =============================================================================
# Example 7: Critical and Triple Point Properties
# =============================================================================

def critical_triple_point_data():
    """
    Retrieve critical and triple point data for various fluids.

    Critical point: Boundary above which distinct liquid/gas phases don't exist
    Triple point: Conditions where solid, liquid, and gas coexist
    """
    print("=" * 70)
    print("Example 7: Critical and Triple Point Properties")
    print("=" * 70)

    fluids = ['Water', 'CO2', 'Nitrogen', 'Oxygen', 'Methane', 'Propane',
              'R134a', 'R410A', 'Ammonia', 'Air']

    print(f"\n{'Fluid':>10} {'T_crit(°C)':>12} {'P_crit(bar)':>13} {'T_triple(°C)':>14} {'M(g/mol)':>11}")
    print("-" * 70)

    for fluid in fluids:
        try:
            T_crit = Props1SI('Tcrit', fluid) - 273.15  # °C
            P_crit = Props1SI('Pcrit', fluid) / 1e5  # bar
            T_triple = Props1SI('Ttriple', fluid) - 273.15  # °C
            M = Props1SI('M', fluid) * 1000  # g/mol

            print(f"{fluid:>10} {T_crit:12.2f} {P_crit:13.2f} {T_triple:14.2f} {M:11.3f}")

        except Exception as e:
            print(f"{fluid:>10}  Error: {e}")

    print()


# =============================================================================
# Example 8: Enthalpy-Entropy (Mollier) Diagram Data
# =============================================================================

def mollier_diagram_data():
    """
    Generate data for plotting on h-s (Mollier) diagram.

    Useful for:
    - Turbine analysis
    - Compressor performance
    - Thermodynamic cycle visualization
    """
    print("=" * 70)
    print("Example 8: Enthalpy-Entropy Data for R134a (Mollier Diagram)")
    print("=" * 70)

    fluid = 'R134a'

    # Generate data along saturation curve
    print("\nSaturation curve data:")
    print(f"{'State':>12} {'T(°C)':>8} {'P(bar)':>9} {'h(kJ/kg)':>11} {'s(kJ/kg·K)':>13}")
    print("-" * 60)

    T_range = np.linspace(-40, 100, 15) + 273.15  # -40°C to 100°C

    for T in T_range:
        try:
            # Saturated liquid
            P = PropsSI('P', 'T', T, 'Q', 0, fluid)
            h_f = PropsSI('H', 'T', T, 'Q', 0, fluid) / 1000
            s_f = PropsSI('S', 'T', T, 'Q', 0, fluid) / 1000

            # Saturated vapor
            h_g = PropsSI('H', 'T', T, 'Q', 1, fluid) / 1000
            s_g = PropsSI('S', 'T', T, 'Q', 1, fluid) / 1000

            print(f"{'Sat. Liquid':>12} {T-273.15:8.1f} {P/1e5:9.2f} {h_f:11.2f} {s_f:13.4f}")
            print(f"{'Sat. Vapor':>12} {T-273.15:8.1f} {P/1e5:9.2f} {h_g:11.2f} {s_g:13.4f}")

        except Exception as e:
            pass

    print()


# =============================================================================
# Example 9: Thermodynamic Process Calculations
# =============================================================================

def isentropic_compression():
    """
    Calculate properties for an isentropic (constant entropy) compression process.

    Common in:
    - Ideal compressor analysis
    - Turbine efficiency calculations
    - Cycle performance
    """
    print("=" * 70)
    print("Example 9: Isentropic Compression of Air")
    print("=" * 70)

    fluid = 'Air'

    # Initial state
    T1 = 25 + 273.15  # 25°C
    P1 = 1e5  # 1 bar

    h1 = PropsSI('H', 'T', T1, 'P', P1, fluid)
    s1 = PropsSI('S', 'T', T1, 'P', P1, fluid)
    rho1 = PropsSI('D', 'T', T1, 'P', P1, fluid)

    print(f"\nInitial state:")
    print(f"  T1 = {T1-273.15:.1f}°C")
    print(f"  P1 = {P1/1e5:.1f} bar")
    print(f"  h1 = {h1/1000:.2f} kJ/kg")
    print(f"  s1 = {s1/1000:.4f} kJ/kg·K")

    # Compression ratios
    pressure_ratios = [2, 5, 10, 20]

    print(f"\nIsentropic compression results:")
    print(f"{'PR':>4} {'P2(bar)':>10} {'T2(°C)':>10} {'h2(kJ/kg)':>12} {'w_comp(kJ/kg)':>16}")
    print("-" * 60)

    for PR in pressure_ratios:
        P2 = P1 * PR

        # Final state (isentropic: s2 = s1)
        T2 = PropsSI('T', 'P', P2, 'S', s1, fluid)
        h2 = PropsSI('H', 'P', P2, 'S', s1, fluid)

        # Compression work
        w_comp = (h2 - h1) / 1000  # kJ/kg

        print(f"{PR:4.0f} {P2/1e5:10.1f} {T2-273.15:10.1f} {h2/1000:12.2f} {w_comp:16.2f}")

    print("\nNote: Isentropic process represents ideal (reversible) compression")
    print()


# =============================================================================
# Example 10: Error Handling and Validation
# =============================================================================

def error_handling_examples():
    """
    Demonstrate proper error handling when querying CoolProp.
    """
    print("=" * 70)
    print("Example 10: Error Handling and Input Validation")
    print("=" * 70)

    # Example 1: Invalid fluid name
    print("\n1. Invalid fluid name:")
    try:
        PropsSI('D', 'T', 300, 'P', 101325, 'InvalidFluid')
    except Exception as e:
        print(f"   Error caught: {e}")

    # Example 2: Temperature out of range
    print("\n2. Temperature out of range:")
    try:
        PropsSI('D', 'T', 10000, 'P', 101325, 'Water')  # 10000 K is too high
    except Exception as e:
        print(f"   Error caught: {type(e).__name__}")

    # Example 3: Invalid input pair in two-phase region
    print("\n3. Invalid input pair (T,P) in two-phase region:")
    try:
        # Water at 100°C and 1 atm is on saturation curve
        PropsSI('D', 'T', 373.15, 'P', 101325, 'Water')
    except Exception as e:
        print(f"   Error caught: {type(e).__name__}")

    # Example 4: Correct way - use quality for saturation
    print("\n4. Correct approach - using quality:")
    try:
        rho = PropsSI('D', 'T', 373.15, 'Q', 1, 'Water')  # Q=1 for saturated vapor
        print(f"   Success: ρ = {rho:.4f} kg/m³")
    except Exception as e:
        print(f"   Error: {e}")

    # Example 5: Checking fluid limits
    print("\n5. Checking fluid temperature limits:")
    fluids_to_check = ['Water', 'R134a', 'CO2']
    for fluid in fluids_to_check:
        T_min = Props1SI('Tmin', fluid) - 273.15
        T_max = Props1SI('Tmax', fluid) - 273.15
        print(f"   {fluid}: {T_min:.1f}°C to {T_max:.1f}°C")

    print()


# =============================================================================
# Main Execution
# =============================================================================

if __name__ == "__main__":
    """
    Run all examples to demonstrate CoolProp capabilities.

    Comment out examples you don't need, or run specific functions individually.
    """

    # Basic properties
    water_properties_vs_temperature()

    # Refrigeration
    refrigerant_saturation_properties()
    refrigeration_cycle_analysis()

    # Transport properties
    viscosity_temperature_dependence()

    # Two-phase
    two_phase_properties()
    steam_table_example()

    # Critical points
    critical_triple_point_data()

    # Advanced
    mollier_diagram_data()
    isentropic_compression()

    # Error handling
    error_handling_examples()

    print("=" * 70)
    print("All examples completed successfully!")
    print("=" * 70)
    print("\nTo run individual examples:")
    print("  python examples.py")
    print("\nOr import and use specific functions:")
    print("  from examples import water_properties_vs_temperature")
    print("  water_properties_vs_temperature()")
    print()
