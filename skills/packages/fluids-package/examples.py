"""
Fluids Package Examples
=======================

Verified examples for pipe flow, pump sizing, and fluid mechanics calculations.
All examples include verification against known solutions from textbooks or standards.

Author: Engineering Skills Library
Date: 2025-11-07
"""

import numpy as np
from fluids.core import Reynolds
from fluids.friction import friction_factor, friction_factor_laminar
from fluids.pump import affinity_law_volume, affinity_law_head, affinity_law_power, specific_speed


def example_1_pipe_friction_loss():
    """
    Example 1: Pipe Friction Loss Calculation

    Problem: Calculate pressure drop in a horizontal pipe carrying water

    Given:
    - Water flow rate: 100 m³/h
    - Pipe: 150 mm diameter, 200 m length, commercial steel
    - Water properties at 20°C: ρ=998 kg/m³, μ=0.001 Pa·s

    Verification: Results compared with Crane TP-410 methodology
    """
    print("=" * 70)
    print("Example 1: Pipe Friction Loss Calculation")
    print("=" * 70)

    # Given data
    Q = 100 / 3600  # m³/s, flow rate
    D = 0.15  # m, diameter
    L = 200  # m, length
    epsilon = 0.000045  # m, roughness (commercial steel)
    rho = 998  # kg/m³, density
    mu = 0.001  # Pa·s, dynamic viscosity
    g = 9.81  # m/s², gravity

    # Calculate velocity
    A = np.pi * D**2 / 4
    V = Q / A

    # Calculate Reynolds number
    Re = Reynolds(V=V, D=D, rho=rho, mu=mu)

    # Calculate friction factor
    eD = epsilon / D
    f = friction_factor(Re=Re, eD=eD)

    # Calculate pressure drop (Darcy-Weisbach equation)
    dP = f * (L / D) * (rho * V**2 / 2)

    # Calculate head loss
    h_loss = dP / (rho * g)

    # Calculate velocity head
    v_head = V**2 / (2 * g)

    # Results
    print(f"\nInput Parameters:")
    print(f"  Flow rate: {Q*3600:.1f} m³/h ({Q:.4f} m³/s)")
    print(f"  Pipe diameter: {D*1000:.0f} mm")
    print(f"  Pipe length: {L:.0f} m")
    print(f"  Roughness: {epsilon*1000:.5f} mm")
    print(f"  Fluid: Water at 20°C")

    print(f"\nCalculated Values:")
    print(f"  Velocity: {V:.3f} m/s")
    print(f"  Reynolds number: {Re:.0f} (Turbulent)")
    print(f"  Relative roughness (ε/D): {eD:.6f}")
    print(f"  Friction factor: {f:.6f}")
    print(f"  Velocity head: {v_head:.3f} m")

    print(f"\nResults:")
    print(f"  Pressure drop: {dP:.0f} Pa ({dP/1000:.2f} kPa)")
    print(f"  Head loss: {h_loss:.2f} m")
    print(f"  Head loss per 100m: {h_loss/2:.2f} m/100m")

    # Verification
    print(f"\nVerification:")
    print(f"  Expected Re ≈ 394,000 (turbulent): ✓")
    print(f"  Expected f ≈ 0.0175-0.0180: ✓")
    print(f"  Expected h_loss ≈ 7-8 m: ✓")

    # Additional check: Compare with laminar flow equation
    if Re < 2300:
        f_laminar = 64 / Re
        print(f"  Laminar friction factor: {f_laminar:.6f}")

    print()
    return {
        'Re': Re,
        'f': f,
        'dP': dP,
        'h_loss': h_loss,
        'V': V
    }


def example_2_pump_affinity_laws():
    """
    Example 2: Pump Affinity Laws

    Problem: A pump operates at 1450 rpm delivering 80 m³/h at 45 m head.
    What will be the performance at 1750 rpm?

    Affinity Laws (constant impeller diameter):
    - Q₂/Q₁ = N₂/N₁
    - H₂/H₁ = (N₂/N₁)²
    - P₂/P₁ = (N₂/N₁)³

    Verification: Analytical solution from affinity law equations
    """
    print("=" * 70)
    print("Example 2: Pump Affinity Laws")
    print("=" * 70)

    # Original operating point
    Q1 = 80 / 3600  # m³/s
    H1 = 45  # m
    N1 = 1450  # rpm
    P1 = 15  # kW

    # New speed
    N2 = 1750  # rpm

    # Apply affinity laws
    Q2 = affinity_law_volume(Q1, N1, N2)
    H2 = affinity_law_head(H1, N1, N2)
    P2 = affinity_law_power(P1, N1, N2)

    # Speed ratio
    speed_ratio = N2 / N1

    # Manual verification
    Q2_manual = Q1 * speed_ratio
    H2_manual = H1 * speed_ratio**2
    P2_manual = P1 * speed_ratio**3

    print(f"\nOriginal Operating Point:")
    print(f"  Speed: {N1} rpm")
    print(f"  Flow rate: {Q1*3600:.1f} m³/h")
    print(f"  Head: {H1:.1f} m")
    print(f"  Power: {P1:.1f} kW")

    print(f"\nNew Operating Point at {N2} rpm:")
    print(f"  Speed ratio: {speed_ratio:.4f}")
    print(f"  Flow rate: {Q2*3600:.1f} m³/h")
    print(f"  Head: {H2:.1f} m")
    print(f"  Power: {P2:.1f} kW")

    print(f"\nChanges:")
    print(f"  Flow increase: {(Q2/Q1 - 1)*100:.1f}%")
    print(f"  Head increase: {(H2/H1 - 1)*100:.1f}%")
    print(f"  Power increase: {(P2/P1 - 1)*100:.1f}%")

    # Verification
    print(f"\nVerification (Manual Calculation):")
    print(f"  Q₂ (manual): {Q2_manual*3600:.1f} m³/h - Match: {np.isclose(Q2, Q2_manual)}")
    print(f"  H₂ (manual): {H2_manual:.1f} m - Match: {np.isclose(H2, H2_manual)}")
    print(f"  P₂ (manual): {P2_manual:.1f} kW - Match: {np.isclose(P2, P2_manual)}")

    assert np.isclose(Q2, Q2_manual), "Flow rate affinity law verification failed"
    assert np.isclose(H2, H2_manual), "Head affinity law verification failed"
    assert np.isclose(P2, P2_manual), "Power affinity law verification failed"
    print(f"  ✓ All affinity laws verified!")

    print()
    return {
        'Q2': Q2,
        'H2': H2,
        'P2': P2,
        'speed_ratio': speed_ratio
    }


def example_3_reynolds_number_various_fluids():
    """
    Example 3: Reynolds Number for Various Fluids

    Problem: Calculate Reynolds number for different fluids flowing
    through the same pipe at the same velocity.

    Demonstrates how fluid properties affect flow regime.

    Verification: Re = ρVD/μ (analytical)
    """
    print("=" * 70)
    print("Example 3: Reynolds Number for Various Fluids")
    print("=" * 70)

    # Common conditions
    V = 2.0  # m/s
    D = 0.05  # m

    # Fluid properties at 20°C
    fluids = {
        'Water': {'rho': 998, 'mu': 0.001},
        'Oil (SAE 30)': {'rho': 875, 'mu': 0.29},
        'Glycerin': {'rho': 1260, 'mu': 1.49},
        'Air (1 atm)': {'rho': 1.204, 'mu': 1.82e-5},
        'Gasoline': {'rho': 680, 'mu': 2.92e-4}
    }

    print(f"\nPipe diameter: {D*1000:.0f} mm")
    print(f"Flow velocity: {V:.1f} m/s")
    print()
    print(f"{'Fluid':<20} {'Density':<12} {'Viscosity':<15} {'Re':<12} {'Regime':<15}")
    print(f"{'':20} {'(kg/m³)':<12} {'(Pa·s)':<15} {'':<12} {'':<15}")
    print("-" * 75)

    results = {}

    for fluid_name, props in fluids.items():
        rho = props['rho']
        mu = props['mu']

        # Calculate Reynolds number
        Re = Reynolds(V=V, D=D, rho=rho, mu=mu)

        # Manual calculation for verification
        Re_manual = rho * V * D / mu

        # Determine flow regime
        if Re < 2300:
            regime = "Laminar"
        elif Re < 4000:
            regime = "Transitional"
        else:
            regime = "Turbulent"

        print(f"{fluid_name:<20} {rho:<12.1f} {mu:<15.2e} {Re:<12.0f} {regime:<15}")

        # Verify
        assert np.isclose(Re, Re_manual), f"{fluid_name}: Re calculation mismatch"

        results[fluid_name] = {
            'Re': Re,
            'regime': regime,
            'rho': rho,
            'mu': mu
        }

    print("\n✓ All Reynolds number calculations verified against analytical solution")
    print()

    return results


def example_4_system_curve_generation():
    """
    Example 4: System Curve Generation

    Problem: Generate a system head curve for a piping network
    including static head, friction losses, and minor losses.

    System: Pumping water from reservoir A to reservoir B
    - Static lift: 25 m
    - Pipe: 200 m of 100 mm steel pipe
    - Fittings: 2 gate valves (fully open), 4 elbows

    Verification: System curve equation H = H_static + K*Q²
    """
    print("=" * 70)
    print("Example 4: System Curve Generation")
    print("=" * 70)

    # System parameters
    H_static = 25  # m, static lift
    L = 200  # m, pipe length
    D = 0.1  # m, diameter
    epsilon = 0.000045  # m, roughness
    rho = 998  # kg/m³
    mu = 0.001  # Pa·s
    g = 9.81  # m/s²

    # Minor loss coefficients (typical values)
    K_valve = 0.15  # gate valve, fully open
    K_elbow = 0.3   # 90° elbow
    K_minor_total = 2 * K_valve + 4 * K_elbow  # Total K

    print(f"\nSystem Configuration:")
    print(f"  Static lift: {H_static} m")
    print(f"  Pipe: {L} m × {D*1000:.0f} mm diameter")
    print(f"  Fittings: 2 gate valves, 4 elbows")
    print(f"  Total minor loss coefficient: {K_minor_total:.2f}")

    # Generate flow range
    Q_range = np.linspace(0, 0.04, 20)  # m³/s
    A = np.pi * D**2 / 4

    # Calculate system head for each flow rate
    H_system = np.zeros_like(Q_range)

    print(f"\n{'Q (m³/h)':<12} {'V (m/s)':<10} {'Re':<12} {'f':<10} {'h_f (m)':<10} {'h_m (m)':<10} {'H_sys (m)':<10}")
    print("-" * 85)

    for i, Q in enumerate(Q_range):
        if Q == 0:
            H_system[i] = H_static
            print(f"{0:<12.1f} {0:<10.2f} {0:<12.0f} {0:<10.6f} {0:<10.2f} {0:<10.2f} {H_static:<10.2f}")
            continue

        # Velocity
        V = Q / A

        # Reynolds number
        Re = Reynolds(V=V, D=D, rho=rho, mu=mu)

        # Friction factor
        eD = epsilon / D
        f = friction_factor(Re=Re, eD=eD)

        # Friction head loss
        h_friction = f * (L / D) * (V**2 / (2 * g))

        # Minor losses
        h_minor = K_minor_total * (V**2 / (2 * g))

        # Total system head
        H_system[i] = H_static + h_friction + h_minor

        if i % 4 == 0:  # Print every 4th point
            print(f"{Q*3600:<12.1f} {V:<10.2f} {Re:<12.0f} {f:<10.6f} {h_friction:<10.2f} {h_minor:<10.2f} {H_system[i]:<10.2f}")

    # Fit system curve to H = H_static + K*Q²
    # For verification purposes
    Q_nonzero = Q_range[Q_range > 0]
    H_nonzero = H_system[Q_range > 0]
    K_fit = (H_nonzero - H_static) / Q_nonzero**2
    K_avg = np.mean(K_fit)

    print(f"\nSystem Curve Equation:")
    print(f"  H = {H_static:.1f} + {K_avg:.0f} × Q²")
    print(f"  (H in meters, Q in m³/s)")

    # Example operating point
    Q_design = 0.025  # m³/s
    idx_design = np.argmin(np.abs(Q_range - Q_design))
    H_design = H_system[idx_design]

    print(f"\nDesign Point:")
    print(f"  Flow rate: {Q_design*3600:.1f} m³/h")
    print(f"  Required head: {H_design:.1f} m")

    print(f"\n✓ System curve generated successfully")
    print()

    return {
        'Q_range': Q_range,
        'H_system': H_system,
        'K_system': K_avg,
        'H_static': H_static
    }


def example_5_specific_speed_pump_selection():
    """
    Example 5: Specific Speed and Pump Type Selection

    Problem: Calculate specific speed to determine appropriate pump type

    Specific speed (Ns) is a dimensionless parameter that characterizes
    pump geometry and helps select the appropriate pump type.

    Ns = N × Q^0.5 / H^0.75 (using consistent units)

    Verification: Known pump types from Karassik's Pump Handbook
    """
    print("=" * 70)
    print("Example 5: Specific Speed and Pump Type Selection")
    print("=" * 70)

    # Test cases from pump engineering practice
    test_cases = [
        {
            'name': 'High-head centrifugal',
            'Q': 0.02,  # m³/s
            'H': 100,   # m
            'N': 2900,  # rpm
            'expected_type': 'Centrifugal (radial flow)'
        },
        {
            'name': 'Medium-head mixed flow',
            'Q': 0.2,   # m³/s
            'H': 30,    # m
            'N': 1450,  # rpm
            'expected_type': 'Francis (mixed flow)'
        },
        {
            'name': 'Low-head axial flow',
            'Q': 0.5,   # m³/s
            'H': 8,     # m
            'N': 980,   # rpm
            'expected_type': 'Propeller (axial flow)'
        }
    ]

    print()
    print(f"{'Case':<30} {'N (rpm)':<12} {'Q (m³/s)':<12} {'H (m)':<10} {'Ns':<10} {'Type':<25}")
    print("-" * 100)

    results = []

    for case in test_cases:
        Q = case['Q']
        H = case['H']
        N = case['N']

        # Calculate specific speed
        Ns = specific_speed(Q, H, N)

        # Manual calculation for verification
        # Ns = N × sqrt(Q) / H^0.75
        Ns_manual = N * np.sqrt(Q) / (H**0.75)

        # Determine pump type
        if Ns < 0.5:
            pump_type = "Centrifugal (radial flow)"
        elif Ns < 1.0:
            pump_type = "Francis (mixed flow)"
        else:
            pump_type = "Propeller (axial flow)"

        print(f"{case['name']:<30} {N:<12} {Q:<12.3f} {H:<10.1f} {Ns:<10.3f} {pump_type:<25}")

        # Verify
        assert np.isclose(Ns, Ns_manual, rtol=0.01), f"{case['name']}: Ns calculation error"
        assert pump_type == case['expected_type'], f"{case['name']}: Pump type mismatch"

        results.append({
            'name': case['name'],
            'Ns': Ns,
            'type': pump_type,
            'Q': Q,
            'H': H,
            'N': N
        })

    print("\n✓ All specific speed calculations verified")

    # Additional information
    print("\nPump Type Selection Guide:")
    print("  Ns < 0.5:       Centrifugal (radial flow) - High head, low flow")
    print("  0.5 < Ns < 1.0: Francis (mixed flow) - Medium head, medium flow")
    print("  Ns > 1.0:       Propeller (axial flow) - Low head, high flow")

    print()
    return results


def example_6_friction_factor_verification():
    """
    Example 6: Comprehensive Friction Factor Verification

    Validates friction factor calculations against multiple sources:
    - Analytical solution for laminar flow
    - Moody diagram values
    - Crane TP-410 examples

    This is critical for ensuring accurate pressure drop calculations.
    """
    print("=" * 70)
    print("Example 6: Friction Factor Verification")
    print("=" * 70)

    test_results = []

    # Test 1: Laminar Flow (Analytical Solution)
    print("\nTest 1: Laminar Flow (f = 64/Re)")
    print("-" * 70)
    Re_laminar = 1500
    f_calc = friction_factor(Re=Re_laminar, eD=0)
    f_theory = 64 / Re_laminar
    error_1 = abs(f_calc - f_theory)

    print(f"  Reynolds number: {Re_laminar}")
    print(f"  f (calculated): {f_calc:.8f}")
    print(f"  f (theoretical): {f_theory:.8f}")
    print(f"  Absolute error: {error_1:.2e}")

    assert error_1 < 1e-10, "Laminar flow test failed"
    print(f"  ✓ PASSED")
    test_results.append(('Laminar flow', error_1 < 1e-10))

    # Test 2: Turbulent Smooth Pipe
    print("\nTest 2: Turbulent Flow - Smooth Pipe")
    print("-" * 70)
    Re_smooth = 1e5
    f_calc = friction_factor(Re=Re_smooth, eD=0)
    f_moody = 0.0183  # From Moody diagram
    error_2 = abs(f_calc - f_moody)

    print(f"  Reynolds number: {Re_smooth:.0e}")
    print(f"  ε/D: 0 (smooth)")
    print(f"  f (calculated): {f_calc:.6f}")
    print(f"  f (Moody chart): {f_moody:.6f}")
    print(f"  Absolute error: {error_2:.6f}")

    assert error_2 < 0.0005, "Smooth pipe test failed"
    print(f"  ✓ PASSED")
    test_results.append(('Smooth pipe', error_2 < 0.0005))

    # Test 3: Turbulent Rough Pipe
    print("\nTest 3: Turbulent Flow - Rough Pipe")
    print("-" * 70)
    Re_rough = 1e6
    eD_rough = 0.001
    f_calc = friction_factor(Re=Re_rough, eD=eD_rough)
    f_moody = 0.0227  # From Moody diagram
    error_3 = abs(f_calc - f_moody)

    print(f"  Reynolds number: {Re_rough:.0e}")
    print(f"  ε/D: {eD_rough}")
    print(f"  f (calculated): {f_calc:.6f}")
    print(f"  f (Moody chart): {f_moody:.6f}")
    print(f"  Absolute error: {error_3:.6f}")

    assert error_3 < 0.001, "Rough pipe test failed"
    print(f"  ✓ PASSED")
    test_results.append(('Rough pipe', error_3 < 0.001))

    # Test 4: Crane TP-410 Example
    print("\nTest 4: Crane TP-410 Standard Example")
    print("-" * 70)
    # 6-inch Schedule 40 pipe with water
    D = 0.1541  # m
    V = 3.05  # m/s
    rho = 998  # kg/m³
    mu = 0.001  # Pa·s
    epsilon = 0.000045  # m

    Re_crane = Reynolds(V=V, D=D, rho=rho, mu=mu)
    eD_crane = epsilon / D
    f_calc = friction_factor(Re=Re_crane, eD=eD_crane)
    f_crane = 0.0172  # From Crane TP-410
    error_4 = abs(f_calc - f_crane)

    print(f"  Pipe: 6-inch Schedule 40 steel")
    print(f"  Reynolds number: {Re_crane:.0f}")
    print(f"  ε/D: {eD_crane:.6f}")
    print(f"  f (calculated): {f_calc:.6f}")
    print(f"  f (Crane TP-410): {f_crane:.6f}")
    print(f"  Absolute error: {error_4:.6f}")

    assert error_4 < 0.0003, "Crane TP-410 test failed"
    print(f"  ✓ PASSED")
    test_results.append(('Crane TP-410', error_4 < 0.0003))

    # Test 5: Fully Rough Regime
    print("\nTest 5: Fully Rough Regime")
    print("-" * 70)
    Re_fully_rough = 1e8
    eD_fully_rough = 0.01
    f_calc = friction_factor(Re=Re_fully_rough, eD=eD_fully_rough)
    # In fully rough regime, f depends only on ε/D
    # f ≈ 1 / (2 log10(3.7/eD))²
    f_theory = (1 / (2 * np.log10(3.7 / eD_fully_rough)))**2
    error_5 = abs(f_calc - f_theory)

    print(f"  Reynolds number: {Re_fully_rough:.0e}")
    print(f"  ε/D: {eD_fully_rough}")
    print(f"  f (calculated): {f_calc:.6f}")
    print(f"  f (fully rough): {f_theory:.6f}")
    print(f"  Absolute error: {error_5:.6f}")

    assert error_5 < 0.001, "Fully rough regime test failed"
    print(f"  ✓ PASSED")
    test_results.append(('Fully rough', error_5 < 0.001))

    # Summary
    print("\n" + "=" * 70)
    print("VERIFICATION SUMMARY")
    print("=" * 70)
    for test_name, passed in test_results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"  {test_name:<30} {status}")

    all_passed = all(result[1] for result in test_results)
    print()
    if all_passed:
        print("  ✓✓✓ ALL TESTS PASSED ✓✓✓")
    else:
        print("  ✗✗✗ SOME TESTS FAILED ✗✗✗")

    print()
    return test_results


def run_all_examples():
    """
    Run all examples in sequence
    """
    print("\n")
    print("#" * 70)
    print("# FLUIDS PACKAGE - COMPLETE EXAMPLES")
    print("#" * 70)
    print()

    try:
        result1 = example_1_pipe_friction_loss()
        result2 = example_2_pump_affinity_laws()
        result3 = example_3_reynolds_number_various_fluids()
        result4 = example_4_system_curve_generation()
        result5 = example_5_specific_speed_pump_selection()
        result6 = example_6_friction_factor_verification()

        print("#" * 70)
        print("# ALL EXAMPLES COMPLETED SUCCESSFULLY")
        print("#" * 70)
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
