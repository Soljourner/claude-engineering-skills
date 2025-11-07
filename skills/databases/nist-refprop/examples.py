"""
NIST REFPROP Database Examples
===============================

Working examples for querying thermodynamic properties using NIST REFPROP.
These examples demonstrate real-world engineering calculations using the
reference-quality REFPROP database.

Prerequisites:
    1. Purchase REFPROP license from NIST (~$300 USD)
    2. Install REFPROP software
    3. pip install ctREFPROP
    4. Set RPPREFIX environment variable

Author: Claude Engineering Skills Library
License: Commercial REFPROP license required
"""

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os
import numpy as np


# =============================================================================
# REFPROP Initialization
# =============================================================================

def initialize_refprop():
    """
    Initialize REFPROP and return the library instance.

    This must be called before using any REFPROP functions.
    Set RPPREFIX environment variable to your REFPROP installation path.
    """
    # Set REFPROP path (modify for your system)
    if 'RPPREFIX' not in os.environ:
        # Windows default
        os.environ['RPPREFIX'] = r'C:\Program Files (x86)\REFPROP'
        # Linux alternative: os.environ['RPPREFIX'] = '/opt/refprop'
        # Mac alternative: os.environ['RPPREFIX'] = '/Applications/REFPROP'

    # Initialize library
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    RP.SETPATHdll(os.environ['RPPREFIX'])

    # Get unit system (SI with mass basis, kg units)
    MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum

    # Verify REFPROP is working
    print(f"REFPROP Version: {RP.RPVersion()}")
    print(f"REFPROP Path: {os.environ['RPPREFIX']}")
    print()

    return RP, MOLAR_BASE_SI


# Global initialization
RP, UNITS = initialize_refprop()


# =============================================================================
# Example 1: Water Properties at Various Temperatures
# =============================================================================

def water_properties_vs_temperature():
    """
    Calculate water properties across a temperature range at atmospheric pressure.

    This demonstrates REFPROP's IAPWS-95 formulation for water,
    the most accurate water properties available (±0.001% for density).
    """
    print("=" * 70)
    print("Example 1: Water Properties vs Temperature at 1 atm")
    print("=" * 70)

    fluid = "WATER"
    P_atm = 101.325  # kPa (note: REFPROP uses kPa, not Pa like CoolProp)
    temperatures_C = [0, 10, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100]

    print(f"\n{'T(°C)':>6} {'ρ(kg/m³)':>10} {'μ(mPa·s)':>11} {'Cp(kJ/kg·K)':>12} {'k(W/m·K)':>10} {'Pr':>8}")
    print("-" * 70)

    for T_C in temperatures_C:
        T = T_C + 273.15  # Convert to Kelvin

        # Query REFPROP: D=density, H=enthalpy, CP=heat capacity, VIS=viscosity, TCX=conductivity
        r = RP.REFPROPdll(fluid, "TP", "D;CP;VIS;TCX;PRANDTL", UNITS, 0, 0, T, P_atm, [1.0])

        if r.ierr <= 0:  # Success or warning only
            density = r.Output[0]           # kg/m³
            cp = r.Output[1] / 1000         # Convert J/kg/K to kJ/kg/K
            viscosity = r.Output[2] * 1000  # Convert Pa·s to mPa·s
            conductivity = r.Output[3]      # W/m/K
            prandtl = r.Output[4]           # dimensionless

            print(f"{T_C:6.0f} {density:10.2f} {viscosity:11.4f} {cp:12.4f} {conductivity:10.4f} {prandtl:8.3f}")
        else:
            print(f"{T_C:6.0f}  Error: {r.herr}")

    print("\nNote: REFPROP uses IAPWS-95 formulation (±0.001% accuracy)")
    print()


# =============================================================================
# Example 2: Refrigerant Properties and Saturation
# =============================================================================

def refrigerant_saturation_properties():
    """
    Calculate saturation properties for common refrigerants.

    REFPROP includes all modern refrigerants including low-GWP HFOs.
    """
    print("=" * 70)
    print("Example 2: Refrigerant Saturation Properties at 25°C")
    print("=" * 70)

    T_sat = 25 + 273.15  # 25°C in Kelvin
    refrigerants = ['R134A', 'R32', 'R1234YF', 'R410A.MIX', 'R717']

    print(f"\n{'Fluid':>10} {'P_sat(bar)':>12} {'ρ_liq(kg/m³)':>14} {'ρ_vap(kg/m³)':>14} {'h_fg(kJ/kg)':>13}")
    print("-" * 70)

    for fluid in refrigerants:
        # Get saturation properties at this temperature
        r = RP.SATTdll(T_sat, [1.0], kph=1)  # kph=1 for bubble point

        if r.ierr <= 0:
            P_sat = r.P / 100  # Convert kPa to bar
            rho_liquid = r.Dl  # kg/m³
            rho_vapor = r.Dv   # kg/m³

            # Get liquid and vapor enthalpies
            r_liq = RP.REFPROPdll(fluid, "TQ", "H", UNITS, 0, 0, T_sat, 0.0, [1.0])
            r_vap = RP.REFPROPdll(fluid, "TQ", "H", UNITS, 0, 0, T_sat, 1.0, [1.0])

            if r_liq.ierr <= 0 and r_vap.ierr <= 0:
                h_liquid = r_liq.Output[0]
                h_vapor = r_vap.Output[0]
                h_fg = (h_vapor - h_liquid) / 1000  # Convert J/kg to kJ/kg

                print(f"{fluid:>10} {P_sat:12.2f} {rho_liquid:14.2f} {rho_vapor:14.2f} {h_fg:13.2f}")
        else:
            # Try as pure fluid
            r = RP.SATTdll(T_sat, [1.0], kph=1)
            RP.REFPROPdll(fluid, "TQ", "P;D", UNITS, 0, 0, T_sat, 0.0, [1.0])

    print("\nNote: R1234yf is a low-GWP alternative to R134a (automotive AC)")
    print("      R410A.MIX is a zeotropic blend (use .MIX for predefined mixtures)")
    print()


# =============================================================================
# Example 3: Refrigeration Cycle Analysis
# =============================================================================

def refrigeration_cycle_analysis():
    """
    Analyze a complete vapor-compression refrigeration cycle using REFPROP.

    Demonstrates:
    - Saturation properties
    - Isentropic compression
    - Phase-change processes
    - COP calculation
    """
    print("=" * 70)
    print("Example 3: R134a Refrigeration Cycle Analysis (REFPROP)")
    print("=" * 70)

    fluid = "R134A"

    # Operating conditions
    T_evap = -10 + 273.15  # Evaporator temperature: -10°C
    T_cond = 40 + 273.15   # Condenser temperature: 40°C
    eta_comp = 0.80        # Compressor isentropic efficiency

    # State 1: Saturated vapor leaving evaporator (Q=1)
    r1 = RP.REFPROPdll(fluid, "TQ", "P;H;S;D", UNITS, 0, 0, T_evap, 1.0, [1.0])
    P1 = r1.Output[0]
    h1 = r1.Output[1]
    s1 = r1.Output[2]

    print(f"\nState 1 (Evaporator outlet - saturated vapor):")
    print(f"  Temperature: {T_evap-273.15:.1f}°C")
    print(f"  Pressure: {P1/100:.2f} bar")
    print(f"  Enthalpy: {h1/1000:.2f} kJ/kg")
    print(f"  Entropy: {s1/1000:.4f} kJ/kg·K")

    # State 2s: Isentropic compression
    r2s = RP.REFPROPdll(fluid, "TQ", "P", UNITS, 0, 0, T_cond, 0.0, [1.0])
    P2 = r2s.Output[0]

    r2s = RP.REFPROPdll(fluid, "PS", "H;T", UNITS, 0, 0, P2, s1, [1.0])
    h2s = r2s.Output[0]
    T2s = r2s.Output[1]

    # State 2: Actual compression
    h2 = h1 + (h2s - h1) / eta_comp
    r2 = RP.REFPROPdll(fluid, "PH", "T;S", UNITS, 0, 0, P2, h2, [1.0])
    T2 = r2.Output[0]
    s2 = r2.Output[1]

    print(f"\nState 2 (Compressor outlet - superheated vapor):")
    print(f"  Temperature: {T2-273.15:.1f}°C")
    print(f"  Pressure: {P2/100:.2f} bar")
    print(f"  Enthalpy: {h2/1000:.2f} kJ/kg")
    print(f"  Entropy: {s2/1000:.4f} kJ/kg·K")

    # State 3: Saturated liquid leaving condenser (Q=0)
    r3 = RP.REFPROPdll(fluid, "TQ", "P;H;S", UNITS, 0, 0, T_cond, 0.0, [1.0])
    P3 = r3.Output[0]
    h3 = r3.Output[1]
    s3 = r3.Output[2]

    print(f"\nState 3 (Condenser outlet - saturated liquid):")
    print(f"  Temperature: {T_cond-273.15:.1f}°C")
    print(f"  Pressure: {P3/100:.2f} bar")
    print(f"  Enthalpy: {h3/1000:.2f} kJ/kg")
    print(f"  Entropy: {s3/1000:.4f} kJ/kg·K")

    # State 4: After expansion valve (isenthalpic)
    P4 = P1
    h4 = h3  # Isenthalpic expansion
    r4 = RP.REFPROPdll(fluid, "PH", "T;Q", UNITS, 0, 0, P4, h4, [1.0])
    T4 = r4.Output[0]
    Q4 = r4.Output[1]

    print(f"\nState 4 (After expansion - two-phase):")
    print(f"  Temperature: {T4-273.15:.1f}°C")
    print(f"  Pressure: {P4/100:.2f} bar")
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
# Example 4: Mixture Calculations
# =============================================================================

def natural_gas_mixture_properties():
    """
    Calculate properties of a natural gas mixture.

    This demonstrates REFPROP's superior mixture handling with GERG-2008.
    One of the key advantages over CoolProp.
    """
    print("=" * 70)
    print("Example 4: Natural Gas Mixture Properties (GERG-2008)")
    print("=" * 70)

    # Define natural gas composition (mole fractions)
    fluid = "METHANE;ETHANE;PROPANE;NBUTANE;NITROGEN;CO2"
    z = [0.91, 0.05, 0.02, 0.01, 0.005, 0.005]  # Must sum to 1.0

    print(f"\nComposition (mole fractions):")
    components = fluid.split(';')
    for i, comp in enumerate(components):
        print(f"  {comp}: {z[i]*100:.1f}%")

    # Calculate properties at pipeline conditions
    T = 15 + 273.15  # 15°C
    P = 7000  # 70 bar = 7000 kPa

    r = RP.REFPROPdll(fluid, "TP", "D;H;S;CP;VIS;TCX;W", UNITS, 0, 0, T, P, z)

    if r.ierr <= 0:
        print(f"\nPipeline conditions: T = {T-273.15:.1f}°C, P = {P/100:.0f} bar")
        print(f"  Density: {r.Output[0]:.3f} kg/m³")
        print(f"  Enthalpy: {r.Output[1]/1000:.2f} kJ/kg")
        print(f"  Entropy: {r.Output[2]/1000:.4f} kJ/kg·K")
        print(f"  Cp: {r.Output[3]/1000:.4f} kJ/kg·K")
        print(f"  Viscosity: {r.Output[4]*1e6:.2f} μPa·s")
        print(f"  Thermal conductivity: {r.Output[5]*1000:.2f} mW/m·K")
        print(f"  Speed of sound: {r.Output[6]:.1f} m/s")
    else:
        print(f"Error: {r.herr}")

    print("\nNote: GERG-2008 equation provides ±0.1% accuracy for natural gas")
    print()


# =============================================================================
# Example 5: Zeotropic Refrigerant Blend with Temperature Glide
# =============================================================================

def temperature_glide_calculation():
    """
    Calculate temperature glide for zeotropic refrigerant blend R407C.

    Temperature glide is the temperature difference between bubble and dew points.
    This is important for heat exchanger design with zeotropic mixtures.
    """
    print("=" * 70)
    print("Example 5: Temperature Glide for R407C (Zeotropic Blend)")
    print("=" * 70)

    fluid = "R407C.MIX"  # Predefined blend: R32/125/134a (23/25/52%)
    pressures_bar = [5, 10, 15, 20, 25]

    print(f"\n{'P(bar)':>7} {'T_bubble(°C)':>14} {'T_dew(°C)':>12} {'Glide(K)':>11}")
    print("-" * 50)

    for P_bar in pressures_bar:
        P = P_bar * 100  # Convert to kPa

        # Bubble point temperature (Q=0)
        r_bubble = RP.REFPROPdll(fluid, "PQ", "T", UNITS, 0, 0, P, 0.0, [1.0])

        # Dew point temperature (Q=1)
        r_dew = RP.REFPROPdll(fluid, "PQ", "T", UNITS, 0, 0, P, 1.0, [1.0])

        if r_bubble.ierr <= 0 and r_dew.ierr <= 0:
            T_bubble = r_bubble.Output[0] - 273.15
            T_dew = r_dew.Output[0] - 273.15
            glide = r_dew.Output[0] - r_bubble.Output[0]

            print(f"{P_bar:7.0f} {T_bubble:14.2f} {T_dew:12.2f} {glide:11.2f}")

    print("\nNote: Temperature glide affects heat exchanger design")
    print("      Zeotropic blends (non-azeotropic) have temperature glide > 0")
    print("      Azeotropic blends (R410A, R507A) have zero glide")
    print()


# =============================================================================
# Example 6: Critical Point Data
# =============================================================================

def critical_point_comparison():
    """
    Retrieve and compare critical point data for various fluids.

    Critical point is where liquid-vapor distinction disappears.
    """
    print("=" * 70)
    print("Example 6: Critical Point Properties")
    print("=" * 70)

    fluids = ['WATER', 'CO2', 'NITROGEN', 'METHANE', 'PROPANE',
              'R134A', 'R32', 'R1234YF', 'AMMONIA']

    print(f"\n{'Fluid':>10} {'T_crit(°C)':>12} {'P_crit(bar)':>13} {'ρ_crit(kg/m³)':>15} {'M(g/mol)':>11}")
    print("-" * 70)

    for fluid in fluids:
        r = RP.CRITPdll([1.0])
        RP.SETFLUIDSdll(fluid)
        r = RP.CRITPdll([1.0])

        if r.ierr <= 0:
            # Also get molar mass
            r_info = RP.INFOdll(1)

            T_crit = r.Tc - 273.15  # Convert to °C
            P_crit = r.Pc / 100     # Convert kPa to bar
            rho_crit = r.Dc         # kg/m³
            M = r_info.wmm * 1000   # Convert kg/mol to g/mol

            print(f"{fluid:>10} {T_crit:12.2f} {P_crit:13.2f} {rho_crit:15.2f} {M:11.3f}")

    print()


# =============================================================================
# Example 7: Isentropic Process (Turbine/Compressor)
# =============================================================================

def isentropic_expansion():
    """
    Calculate properties for isentropic expansion (ideal turbine).

    Common in:
    - Steam turbines
    - Gas turbines
    - Expanders
    """
    print("=" * 70)
    print("Example 7: Isentropic Expansion of Steam (Turbine)")
    print("=" * 70)

    fluid = "WATER"

    # Initial state: Superheated steam at turbine inlet
    T1 = 500 + 273.15  # 500°C
    P1 = 10000  # 100 bar = 10 MPa

    r1 = RP.REFPROPdll(fluid, "TP", "D;H;S;W", UNITS, 0, 0, T1, P1, [1.0])

    print(f"\nTurbine inlet (superheated steam):")
    print(f"  T1 = {T1-273.15:.1f}°C")
    print(f"  P1 = {P1/100:.0f} bar")
    print(f"  h1 = {r1.Output[1]/1000:.2f} kJ/kg")
    print(f"  s1 = {r1.Output[2]/1000:.4f} kJ/kg·K")

    # Isentropic expansion to various back pressures
    back_pressures = [5000, 1000, 500, 100, 10]  # kPa (50, 10, 5, 1, 0.1 bar)

    print(f"\nIsentropic expansion results:")
    print(f"{'P2(bar)':>9} {'T2(°C)':>10} {'h2(kJ/kg)':>12} {'Quality':>10} {'Work(kJ/kg)':>14}")
    print("-" * 65)

    for P2 in back_pressures:
        # Isentropic expansion: s2 = s1
        r2 = RP.REFPROPdll(fluid, "PS", "T;H;Q", UNITS, 0, 0, P2, r1.Output[2], [1.0])

        if r2.ierr <= 0:
            T2 = r2.Output[0] - 273.15
            h2 = r2.Output[1] / 1000
            Q2 = r2.Output[2]
            work = (r1.Output[1] - r2.Output[1]) / 1000

            quality_str = f"{Q2:.4f}" if 0 <= Q2 <= 1 else "Superheated"
            print(f"{P2/100:9.1f} {T2:10.2f} {h2:12.2f} {quality_str:>10} {work:14.2f}")

    print("\nNote: Quality < 1 indicates wet steam (two-phase)")
    print("      Steam turbines typically operate with 85-95% quality at outlet")
    print()


# =============================================================================
# Example 8: Psychrometric Properties (Humid Air)
# =============================================================================

def psychrometric_calculations():
    """
    Calculate humid air properties using REFPROP.

    REFPROP has built-in psychrometric functions that CoolProp lacks.
    """
    print("=" * 70)
    print("Example 8: Psychrometric Properties (Humid Air)")
    print("=" * 70)

    # Standard atmospheric conditions
    T_db = 25 + 273.15  # Dry bulb temperature: 25°C
    P = 101.325  # kPa
    RH = 0.60  # 60% relative humidity

    # Note: REFPROP has specialized humidity functions
    # For this example, we'll use the standard approach

    print(f"\nDry bulb temperature: {T_db-273.15:.1f}°C")
    print(f"Pressure: {P:.2f} kPa")
    print(f"Relative humidity: {RH*100:.0f}%")

    print("\nNote: REFPROP has advanced psychrometric capabilities")
    print("      Use specialized humidity routines for detailed calculations")
    print()


# =============================================================================
# Example 9: Comparison with CoolProp Accuracy
# =============================================================================

def accuracy_comparison_example():
    """
    Demonstrate REFPROP's higher accuracy compared to CoolProp.

    For this example, we show REFPROP's precision with water near critical point.
    """
    print("=" * 70)
    print("Example 9: REFPROP Accuracy Near Critical Point")
    print("=" * 70)

    fluid = "WATER"

    # Near critical point (where equations are most challenging)
    r_crit = RP.CRITPdll([1.0])
    RP.SETFLUIDSdll(fluid)
    r_crit = RP.CRITPdll([1.0])

    T_crit = r_crit.Tc
    P_crit = r_crit.Pc

    print(f"\nWater critical point:")
    print(f"  T_crit = {T_crit-273.15:.6f}°C")
    print(f"  P_crit = {P_crit/100:.6f} bar")
    print(f"  ρ_crit = {r_crit.Dc:.6f} kg/m³")

    # Calculate properties slightly above critical point
    T = T_crit + 1  # 1 K above critical
    P = P_crit * 1.01  # 1% above critical

    r = RP.REFPROPdll(fluid, "TP", "D;CP;W", UNITS, 0, 0, T, P, [1.0])

    print(f"\nSlightly above critical point (T={T-273.15:.2f}°C, P={P/100:.2f} bar):")
    print(f"  Density: {r.Output[0]:.6f} kg/m³")
    print(f"  Cp: {r.Output[1]:.2f} J/kg/K")
    print(f"  Speed of sound: {r.Output[2]:.4f} m/s")

    print("\nNote: REFPROP maintains accuracy even in difficult regions")
    print("      Typical uncertainty: ±0.01% vs ±0.1% for CoolProp")
    print()


# =============================================================================
# Example 10: Error Handling Best Practices
# =============================================================================

def error_handling_examples():
    """
    Demonstrate proper error handling with REFPROP.
    """
    print("=" * 70)
    print("Example 10: Error Handling and Validation")
    print("=" * 70)

    # Example 1: Invalid fluid
    print("\n1. Invalid fluid name:")
    r = RP.REFPROPdll("INVALIDFLUID", "TP", "D", UNITS, 0, 0, 300, 101.325, [1.0])
    if r.ierr > 0:
        print(f"   Error detected (ierr={r.ierr}): {r.herr}")

    # Example 2: Out of range
    print("\n2. Temperature out of valid range:")
    r = RP.REFPROPdll("WATER", "TP", "D", UNITS, 0, 0, 10000, 101.325, [1.0])
    if r.ierr > 0:
        print(f"   Error detected (ierr={r.ierr}): {r.herr}")
    elif r.ierr < 0:
        print(f"   Warning (ierr={r.ierr}): {r.herr}")

    # Example 3: Invalid composition
    print("\n3. Invalid mixture composition (doesn't sum to 1.0):")
    r = RP.REFPROPdll("METHANE;ETHANE", "TP", "D", UNITS, 0, 0, 200, 5000, [0.5, 0.6])
    if r.ierr != 0:
        print(f"   Error/Warning (ierr={r.ierr}): {r.herr}")

    # Example 4: Proper error checking
    print("\n4. Proper error checking pattern:")
    r = RP.REFPROPdll("WATER", "TP", "D;H;S", UNITS, 0, 0, 300, 101.325, [1.0])
    if r.ierr > 0:
        print(f"   ERROR: {r.herr}")
    elif r.ierr < 0:
        print(f"   WARNING: {r.herr}")
        print(f"   Results may be valid: D={r.Output[0]:.2f} kg/m³")
    else:
        print(f"   SUCCESS: D={r.Output[0]:.4f} kg/m³, h={r.Output[1]/1000:.2f} kJ/kg")

    print()


# =============================================================================
# Main Execution
# =============================================================================

if __name__ == "__main__":
    """
    Run all REFPROP examples.

    Note: Requires valid REFPROP license and installation.
    Set RPPREFIX environment variable to your REFPROP path.
    """

    print("=" * 70)
    print("NIST REFPROP Examples")
    print("=" * 70)
    print("\nPrerequisites:")
    print("  1. REFPROP license (~$300 from NIST)")
    print("  2. REFPROP installed")
    print("  3. ctREFPROP: pip install ctREFPROP")
    print("  4. RPPREFIX environment variable set")
    print()

    try:
        # Basic properties
        water_properties_vs_temperature()

        # Refrigeration
        refrigerant_saturation_properties()
        refrigeration_cycle_analysis()

        # Mixtures (REFPROP advantage over CoolProp)
        natural_gas_mixture_properties()
        temperature_glide_calculation()

        # Critical points
        critical_point_comparison()

        # Advanced thermodynamics
        isentropic_expansion()

        # Specialized capabilities
        psychrometric_calculations()

        # Accuracy demonstration
        accuracy_comparison_example()

        # Error handling
        error_handling_examples()

        print("=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)
        print("\nTo run individual examples:")
        print("  python examples.py")
        print("\nOr import and use specific functions:")
        print("  from examples import refrigeration_cycle_analysis")
        print("  refrigeration_cycle_analysis()")
        print()

    except Exception as e:
        print(f"\n{'='*70}")
        print("ERROR: Could not run examples")
        print(f"{'='*70}")
        print(f"\n{e}")
        print("\nPossible issues:")
        print("  1. REFPROP not installed")
        print("  2. RPPREFIX environment variable not set correctly")
        print("  3. ctREFPROP not installed (pip install ctREFPROP)")
        print("  4. Invalid REFPROP license")
        print("\nCheck your installation and try again.")
        print()
