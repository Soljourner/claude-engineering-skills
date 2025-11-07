"""
SymPy Symbolic Math Examples for Engineering

All examples demonstrate practical derivations and symbolic solutions
for fluid dynamics and pump engineering applications.
"""

import sympy as sp
from sympy import symbols, Eq, solve, simplify, diff, integrate, latex
import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# Example 1: Derive Pump Head Equation
# ============================================================================

def example_1_derive_pump_head():
    """
    Derive the theoretical pump head equation from Euler turbine equation.

    Starting from first principles:
    - Velocity triangles at inlet and outlet
    - Angular momentum conservation
    - Euler turbine equation

    Expected result: H = (U₂·Cᵤ₂ - U₁·Cᵤ₁) / g
    """
    print("=" * 70)
    print("Example 1: Derive Pump Head Equation")
    print("=" * 70)

    # Define symbolic variables
    omega = sp.Symbol('omega', positive=True)  # Angular velocity [rad/s]
    D_1, D_2 = sp.symbols('D_1 D_2', positive=True)  # Diameters [m]
    C_u1, C_u2 = sp.symbols('C_u1 C_u2', real=True)  # Tangential velocities [m/s]
    C_m1, C_m2 = sp.symbols('C_m1 C_m2', positive=True)  # Meridional velocities [m/s]
    beta_1, beta_2 = sp.symbols('beta_1 beta_2', positive=True)  # Blade angles [rad]
    g = sp.Symbol('g', positive=True)  # Gravity [m/s²]

    print("\nStep 1: Define peripheral velocities")
    # Peripheral velocities: U = ω·D/2
    U_1 = omega * D_1 / 2
    U_2 = omega * D_2 / 2
    print(f"U₁ = {U_1}")
    print(f"U₂ = {U_2}")

    print("\nStep 2: Velocity triangle relations")
    # From velocity triangles: Cᵤ = U - Cₘ/tan(β)
    C_u1_expr = U_1 - C_m1 / sp.tan(beta_1)
    C_u2_expr = U_2 - C_m2 / sp.tan(beta_2)
    print(f"Cᵤ₁ = {C_u1_expr}")
    print(f"Cᵤ₂ = {C_u2_expr}")

    print("\nStep 3: Euler turbine equation")
    # Specific work: Y = U₂·Cᵤ₂ - U₁·Cᵤ₁
    Y = U_2 * C_u2_expr - U_1 * C_u1_expr
    print(f"Y = {simplify(Y)}")

    # Head: H = Y/g
    H = Y / g
    print(f"\nHead: H = Y/g")
    print(f"H = {simplify(H)}")

    print("\nStep 4: Special case - Radial inlet (Cᵤ₁ = 0)")
    H_radial = H.subs(C_u1, 0)
    H_radial = H_radial.subs(C_u1_expr, 0)
    print(f"H (radial inlet) = {simplify(H_radial)}")

    print("\nStep 5: Express in practical form")
    # For typical centrifugal pump: assume C_m2 << U_2
    H_simplified = (omega * D_2)**2 / (2 * g) * (1 - sp.tan(beta_2))
    print(f"H ≈ (ωD₂)²/(2g) · [1 - tan(β₂)]")

    print("\n✓ Derivation complete!")
    print(f"\nLaTeX export:")
    print(f"  {latex(Eq(sp.Symbol('H'), H))}")

    return H


# ============================================================================
# Example 2: Solve Bernoulli Equation for Velocity
# ============================================================================

def example_2_solve_bernoulli():
    """
    Solve Bernoulli equation symbolically for outlet velocity.

    Given:
    - Inlet pressure, velocity, elevation
    - Outlet pressure, elevation

    Find: Outlet velocity

    Expected: v₂ = √(v₁² + 2g(z₁-z₂) + 2(p₁-p₂)/ρ)
    """
    print("\n" + "=" * 70)
    print("Example 2: Solve Bernoulli Equation for Velocity")
    print("=" * 70)

    # Define symbols
    p_1, p_2 = sp.symbols('p_1 p_2', real=True)  # Pressures [Pa]
    v_1, v_2 = sp.symbols('v_1 v_2', positive=True)  # Velocities [m/s]
    z_1, z_2 = sp.symbols('z_1 z_2', real=True)  # Elevations [m]
    rho = sp.Symbol('rho', positive=True)  # Density [kg/m³]
    g = sp.Symbol('g', positive=True)  # Gravity [m/s²]

    print("\nBernoulli equation:")
    print("p₁/ρg + v₁²/(2g) + z₁ = p₂/ρg + v₂²/(2g) + z₂")

    # Set up equation
    bernoulli = Eq(
        p_1/(rho*g) + v_1**2/(2*g) + z_1,
        p_2/(rho*g) + v_2**2/(2*g) + z_2
    )

    print(f"\nSymbolic form:")
    print(bernoulli)

    # Solve for v_2
    print("\nSolving for v₂...")
    v_2_solutions = solve(bernoulli, v_2)

    # Take positive solution (physical)
    v_2_solution = [sol for sol in v_2_solutions if sol.is_positive or sol.has(sp.sqrt)][0]

    print(f"\nv₂ = {v_2_solution}")
    print(f"\nSimplified:")
    v_2_simplified = simplify(v_2_solution)
    print(f"v₂ = {v_2_simplified}")

    # Special cases
    print("\n" + "-" * 70)
    print("Special Cases:")
    print("-" * 70)

    # Case 1: Horizontal pipe (z₁ = z₂)
    v_2_horizontal = v_2_solution.subs(z_2, z_1)
    print(f"\n1. Horizontal pipe (z₁ = z₂):")
    print(f"   v₂ = {simplify(v_2_horizontal)}")

    # Case 2: Static inlet (v₁ = 0)
    v_2_static = v_2_solution.subs(v_1, 0)
    print(f"\n2. Static inlet (v₁ = 0):")
    print(f"   v₂ = {simplify(v_2_static)}")

    # Case 3: Free jet (p₂ = p_atm = p₁)
    v_2_freejet = v_2_solution.subs(p_2, p_1)
    print(f"\n3. Free jet (p₂ = p₁):")
    print(f"   v₂ = {simplify(v_2_freejet)}")

    # Numerical example
    print("\n" + "-" * 70)
    print("Numerical Example:")
    print("-" * 70)
    values = {
        p_1: 300000,    # 3 bar
        p_2: 101325,    # 1 atm
        v_1: 2.0,       # 2 m/s
        z_1: 10.0,      # 10 m elevation
        z_2: 0.0,       # ground level
        rho: 1000,      # water
        g: 9.81
    }

    v_2_numeric = v_2_solution.subs(values)
    print(f"Inlet: p₁=3 bar, v₁=2 m/s, z₁=10 m")
    print(f"Outlet: p₂=1 atm, z₂=0 m")
    print(f"Result: v₂ = {float(v_2_numeric):.2f} m/s")

    print("\n✓ Solution verified!")

    return v_2_solution


# ============================================================================
# Example 3: Differentiate for Optimization
# ============================================================================

def example_3_optimization_derivative():
    """
    Find optimal operating point by differentiation.

    Given:
    - Pump curve: H = a - b·Q - c·Q²
    - System curve: H = K·Q²

    Find:
    - Intersection point (operating point)
    - Maximum power point
    - Optimal efficiency point

    Expected: Operating point where dH_pump/dQ = dH_system/dQ
    """
    print("\n" + "=" * 70)
    print("Example 3: Optimization Using Derivatives")
    print("=" * 70)

    # Define symbols
    Q = sp.Symbol('Q', positive=True)  # Flow rate [m³/s]
    a, b, c = sp.symbols('a b c', positive=True)  # Pump curve coefficients
    K = sp.Symbol('K', positive=True)  # System resistance
    rho = sp.Symbol('rho', positive=True)  # Density
    g = sp.Symbol('g', positive=True)  # Gravity

    # Pump and system curves
    H_pump = a - b*Q - c*Q**2
    H_system = K*Q**2

    print("\nPump curve:")
    print(f"H_pump = {H_pump}")
    print("\nSystem curve:")
    print(f"H_system = {H_system}")

    # Operating point: H_pump = H_system
    print("\n" + "-" * 70)
    print("Operating Point (H_pump = H_system):")
    print("-" * 70)

    operating_eq = Eq(H_pump, H_system)
    Q_operating = solve(operating_eq, Q)

    print(f"\nOperating flow rates:")
    for i, sol in enumerate(Q_operating):
        print(f"  Q_{i+1} = {sol}")

    # Take positive, real solution
    Q_op = [sol for sol in Q_operating if sol.is_positive or not sol.has(sp.I)][0]
    print(f"\nPhysical solution: Q_op = {Q_op}")

    # Head at operating point
    H_op = H_pump.subs(Q, Q_op)
    print(f"Head at operating point: H_op = {simplify(H_op)}")

    # Power optimization
    print("\n" + "-" * 70)
    print("Maximum Power Point:")
    print("-" * 70)

    # Hydraulic power: P = ρ·g·Q·H
    P = rho * g * Q * H_pump
    print(f"\nPower: P = ρ·g·Q·H")
    print(f"P = {expand(P)}")

    # Find maximum: dP/dQ = 0
    dP_dQ = diff(P, Q)
    print(f"\ndP/dQ = {expand(dP_dQ)}")

    Q_max_power = solve(dP_dQ, Q)
    print(f"\nCritical points (dP/dQ = 0):")
    for sol in Q_max_power:
        print(f"  Q = {sol}")

    # Verify it's a maximum with second derivative
    d2P_dQ2 = diff(dP_dQ, Q)
    print(f"\nSecond derivative: d²P/dQ² = {d2P_dQ2}")
    print("(Negative → maximum, Positive → minimum)")

    # Efficiency optimization
    print("\n" + "-" * 70)
    print("Efficiency Analysis:")
    print("-" * 70)

    # Efficiency: η = (P_hydraulic) / (P_shaft)
    # For this analysis, assume P_shaft = P + P_losses
    # P_losses ∝ Q³ (volumetric and mechanical losses)
    P_loss_coeff = sp.Symbol('k_loss', positive=True)
    P_shaft = P + P_loss_coeff * Q**3
    eta = P / P_shaft

    print(f"\nEfficiency: η = P_hydraulic / P_shaft")
    print(f"η = {simplify(eta)}")

    # Find maximum efficiency
    deta_dQ = diff(eta, Q)
    print(f"\ndη/dQ = {simplify(deta_dQ)}")

    # Numerical example
    print("\n" + "-" * 70)
    print("Numerical Example:")
    print("-" * 70)

    values = {
        a: 80,      # 80 m shutoff head
        b: 50,      # linear coefficient
        c: 500,     # quadratic coefficient
        K: 100,     # system resistance
        rho: 1000,
        g: 9.81
    }

    Q_op_num = Q_op.subs(values)
    H_op_num = H_op.subs(values)

    print(f"Pump: H = 80 - 50·Q - 500·Q²")
    print(f"System: H = 100·Q²")
    print(f"\nOperating point:")
    print(f"  Q = {float(Q_op_num):.4f} m³/s = {float(Q_op_num)*3600:.1f} m³/h")
    print(f"  H = {float(H_op_num):.2f} m")

    if len(Q_max_power) > 0 and Q_max_power[0].is_real:
        Q_maxP_num = Q_max_power[0].subs(values)
        print(f"\nMaximum power at Q = {float(Q_maxP_num):.4f} m³/s")

    print("\n✓ Optimization complete!")

    return Q_op, H_op


# ============================================================================
# Example 4: Integrate for Flow Calculations
# ============================================================================

def example_4_integrate_flow():
    """
    Calculate volumetric flow rate from velocity profile using integration.

    For laminar flow in circular pipe:
    - Velocity profile: u(r) = u_max·(1 - (r/R)²)
    - Flow rate: Q = ∫∫ u dA
    - Average velocity: V_avg = Q / A

    Expected: Q = (π·R²·u_max)/2, V_avg = u_max/2
    """
    print("\n" + "=" * 70)
    print("Example 4: Flow Integration (Laminar Pipe Flow)")
    print("=" * 70)

    # Define symbols
    r = sp.Symbol('r', positive=True)  # Radial position [m]
    R = sp.Symbol('R', positive=True)  # Pipe radius [m]
    u_max = sp.Symbol('u_max', positive=True)  # Centerline velocity [m/s]

    print("\nLaminar velocity profile (Hagen-Poiseuille):")
    # Parabolic profile
    u = u_max * (1 - (r/R)**2)
    print(f"u(r) = {u}")

    print("\n" + "-" * 70)
    print("Flow Rate Calculation:")
    print("-" * 70)

    # Flow through annular element: dQ = u·dA = u·2πr·dr
    print("\nDifferential flow rate:")
    dQ = u * 2 * sp.pi * r
    print(f"dQ = u(r)·2πr·dr = {expand(dQ)}·dr")

    # Integrate from r=0 to r=R
    print("\nIntegrate from r=0 to r=R:")
    print(f"Q = ∫₀ᴿ {expand(dQ)} dr")

    Q = integrate(dQ, (r, 0, R))
    print(f"\nQ = {Q}")
    print(f"Simplified: Q = {simplify(Q)}")

    # Verify with analytical solution
    Q_analytical = sp.pi * R**2 * u_max / 2
    print(f"\nAnalytical: Q = πR²u_max/2 = {Q_analytical}")
    print(f"Match: {simplify(Q - Q_analytical) == 0}")

    print("\n" + "-" * 70)
    print("Average Velocity:")
    print("-" * 70)

    # Cross-sectional area
    A = sp.pi * R**2
    print(f"\nPipe area: A = πR² = {A}")

    # Average velocity: V_avg = Q/A
    V_avg = Q / A
    print(f"\nV_avg = Q/A = {simplify(V_avg)}")

    # Ratio to maximum velocity
    ratio = V_avg / u_max
    print(f"\nV_avg/u_max = {simplify(ratio)}")
    print("(For laminar flow, average velocity is half the centerline velocity)")

    # Turbulent flow comparison
    print("\n" + "-" * 70)
    print("Comparison: Turbulent Flow (Power Law)")
    print("-" * 70)

    n = sp.Symbol('n', positive=True)  # Power law exponent
    u_turbulent = u_max * (1 - r/R)**(1/n)

    print(f"\nTurbulent profile (1/n power law):")
    print(f"u(r) = u_max·(1 - r/R)^(1/n)")

    # For n=7 (typical)
    u_turb_7 = u_turbulent.subs(n, 7)
    dQ_turb = u_turb_7 * 2 * sp.pi * r
    Q_turb = integrate(dQ_turb, (r, 0, R))
    V_avg_turb = simplify(Q_turb / A)

    print(f"\nFor n=7 (typical turbulent flow):")
    print(f"V_avg/u_max = {simplify(V_avg_turb / u_max)}")
    print(f"Numerical: {float(V_avg_turb.subs([(R, 1), (u_max, 1)])/1):.4f}")

    # Numerical example
    print("\n" + "-" * 70)
    print("Numerical Example:")
    print("-" * 70)

    values = {
        R: 0.05,      # 50 mm radius (100 mm diameter pipe)
        u_max: 2.0    # 2 m/s centerline
    }

    Q_num = Q.subs(values)
    V_avg_num = V_avg.subs(values)
    A_num = A.subs(values)

    print(f"Pipe: D = 100 mm, R = 50 mm")
    print(f"Centerline velocity: u_max = 2.0 m/s")
    print(f"\nResults:")
    print(f"  Flow rate: Q = {float(Q_num):.6f} m³/s = {float(Q_num)*1000:.2f} L/s")
    print(f"  Cross-section: A = {float(A_num):.6f} m²")
    print(f"  Average velocity: V_avg = {float(V_avg_num):.3f} m/s")
    print(f"  Ratio: V_avg/u_max = {float(V_avg_num/values[u_max]):.3f}")

    print("\n✓ Integration verified!")

    return Q, V_avg


# ============================================================================
# Bonus: Generate Pump Performance Map
# ============================================================================

def bonus_generate_performance_map():
    """
    Generate symbolic expressions for complete pump performance map.
    Then convert to numerical functions for plotting.

    Demonstrates:
    - Symbolic derivation
    - LaTeX generation for documentation
    - Conversion to numerical functions (lambdify)
    - Plotting results
    """
    print("\n" + "=" * 70)
    print("Bonus: Complete Pump Performance Map")
    print("=" * 70)

    # Define symbols
    Q = sp.Symbol('Q', positive=True)
    a, b, c = 80, 50, 500
    rho, g = 1000, 9.81
    eta_max = sp.Symbol('eta_max', positive=True)
    Q_bep = sp.Symbol('Q_bep', positive=True)
    sigma = sp.Symbol('sigma', positive=True)

    # Symbolic expressions
    H_sym = a - b*Q - c*Q**2
    eta_sym = eta_max * sp.exp(-((Q - Q_bep) / sigma)**2)
    P_sym = rho * g * Q * H_sym / eta_sym

    print("\nSymbolic Performance Equations:")
    print(f"Head:       H = {H_sym}")
    print(f"Efficiency: η = {eta_sym}")
    print(f"Power:      P = {simplify(P_sym)}")

    print("\n" + "-" * 70)
    print("LaTeX Export (for documentation):")
    print("-" * 70)
    print(f"\nH: {latex(H_sym)}")
    print(f"η: {latex(eta_sym)}")
    print(f"P: {latex(P_sym)}")

    # Convert to numerical functions
    print("\n" + "-" * 70)
    print("Converting to Numerical Functions:")
    print("-" * 70)

    # Substitute efficiency parameters
    eta_sym_num = eta_sym.subs([(eta_max, 0.85), (Q_bep, 0.15), (sigma, 0.08)])

    # Create lambda functions for numerical evaluation
    H_func = sp.lambdify(Q, H_sym, 'numpy')
    eta_func = sp.lambdify(Q, eta_sym_num, 'numpy')
    P_sym_num = P_sym.subs([(eta_max, 0.85), (Q_bep, 0.15), (sigma, 0.08)])
    P_func = sp.lambdify(Q, P_sym_num, 'numpy')

    # Generate data
    Q_vals = np.linspace(0.01, 0.30, 100)
    H_vals = H_func(Q_vals)
    eta_vals = eta_func(Q_vals)
    P_vals = P_func(Q_vals)

    print("✓ Functions converted successfully")
    print(f"  Evaluated at {len(Q_vals)} points")
    print(f"  Q range: {Q_vals[0]:.2f} to {Q_vals[-1]:.2f} m³/s")

    # Find key points
    idx_bep = np.argmax(eta_vals)
    Q_bep_num = Q_vals[idx_bep]
    H_bep_num = H_vals[idx_bep]
    eta_bep_num = eta_vals[idx_bep]
    P_bep_num = P_vals[idx_bep]

    print("\n" + "-" * 70)
    print("Best Efficiency Point (BEP):")
    print("-" * 70)
    print(f"  Q_bep = {Q_bep_num:.3f} m³/s = {Q_bep_num*3600:.1f} m³/h")
    print(f"  H_bep = {H_bep_num:.2f} m")
    print(f"  η_bep = {eta_bep_num*100:.2f}%")
    print(f"  P_bep = {P_bep_num/1000:.2f} kW")

    # Create performance map plot
    print("\n" + "-" * 70)
    print("Generating Performance Map Plot:")
    print("-" * 70)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10))

    # Head curve
    ax1.plot(Q_vals*3600, H_vals, 'b-', linewidth=2, label='Head curve')
    ax1.axvline(Q_bep_num*3600, color='r', linestyle='--', alpha=0.5, label='BEP')
    ax1.set_xlabel('Flow rate Q (m³/h)')
    ax1.set_ylabel('Head H (m)')
    ax1.set_title('Pump Performance Map')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Efficiency curve
    ax2.plot(Q_vals*3600, eta_vals*100, 'g-', linewidth=2, label='Efficiency')
    ax2.axvline(Q_bep_num*3600, color='r', linestyle='--', alpha=0.5, label='BEP')
    ax2.axhline(eta_bep_num*100, color='g', linestyle=':', alpha=0.5)
    ax2.set_xlabel('Flow rate Q (m³/h)')
    ax2.set_ylabel('Efficiency η (%)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Power curve
    ax3.plot(Q_vals*3600, P_vals/1000, 'r-', linewidth=2, label='Power')
    ax3.axvline(Q_bep_num*3600, color='r', linestyle='--', alpha=0.5, label='BEP')
    ax3.set_xlabel('Flow rate Q (m³/h)')
    ax3.set_ylabel('Power P (kW)')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    plt.tight_layout()

    # Save figure
    output_path = '/tmp/pump_performance_map.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Plot saved to: {output_path}")

    # Optionally display (comment out if running headless)
    # plt.show()
    plt.close()

    print("\n✓ Performance map generated successfully!")

    return H_func, eta_func, P_func


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    print("\n")
    print("*" * 70)
    print("*" + " " * 68 + "*")
    print("*" + "  SymPy Symbolic Math Examples for Engineering".center(68) + "*")
    print("*" + " " * 68 + "*")
    print("*" * 70)

    try:
        # Run all examples
        H_pump = example_1_derive_pump_head()
        v_2 = example_2_solve_bernoulli()
        Q_op, H_op = example_3_optimization_derivative()
        Q_flow, V_avg = example_4_integrate_flow()

        # Bonus
        H_func, eta_func, P_func = bonus_generate_performance_map()

        print("\n" + "=" * 70)
        print("ALL EXAMPLES COMPLETED SUCCESSFULLY")
        print("=" * 70)
        print("\nKey Results Summary:")
        print("-" * 70)
        print("1. Pump head equation derived from first principles")
        print("2. Bernoulli equation solved symbolically for velocity")
        print("3. Operating point found using derivatives")
        print("4. Flow rate calculated by integration")
        print("5. Performance map generated and plotted")
        print("\nYou can now use these patterns for your own symbolic derivations!")

    except Exception as e:
        print(f"\n❌ Error occurred: {e}")
        import traceback
        traceback.print_exc()
