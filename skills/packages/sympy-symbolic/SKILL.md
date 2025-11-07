---
name: sympy-symbolic
description: "Derive and solve fluid dynamics and pump equations symbolically"
category: packages
domain: general
complexity: intermediate
dependencies:
  - sympy
---

# SymPy Symbolic Math Skill

## Overview

SymPy is a Python library for symbolic mathematics. Unlike numerical libraries (NumPy, SciPy), SymPy works with exact symbolic expressions, allowing you to:

- Derive equations analytically
- Simplify complex mathematical expressions
- Solve equations and systems of equations symbolically
- Perform calculus operations (derivatives, integrals, limits)
- Verify dimensional consistency and unit relationships
- Generate LaTeX equations for documentation

For engineering applications, SymPy is essential for:
- Deriving pump and turbine performance equations from first principles
- Solving flow equations symbolically (Bernoulli, continuity, momentum)
- Simplifying complex expressions before numerical implementation
- Verifying analytical solutions against numerical results
- Performing dimensional analysis and unit verification
- Generating documentation with proper mathematical notation

## Installation

```bash
pip install sympy
```

For engineering workflows with visualization:

```bash
pip install sympy matplotlib numpy
```

## Core Concepts

### Symbols vs Numbers

SymPy works with symbolic variables, not numeric values:

```python
import sympy as sp

# Define symbolic variables
x, y, z = sp.symbols('x y z')

# Symbolic expression (not evaluated numerically)
expr = x**2 + 2*x*y + y**2
print(expr)  # x**2 + 2*x*y + y**2

# Simplify
simplified = sp.simplify(expr)
print(simplified)  # (x + y)**2

# Substitute numerical values
result = expr.subs([(x, 3), (y, 2)])
print(result)  # 25
```

**Key Differences from NumPy:**
- SymPy: Exact symbolic computation (√2, π, fractions)
- NumPy: Approximate floating-point computation (1.414..., 3.14159...)

### Creating Symbols

```python
import sympy as sp

# Single symbol
x = sp.Symbol('x')

# Multiple symbols
a, b, c = sp.symbols('a b c')

# With assumptions (important for engineering)
D = sp.Symbol('D', positive=True)  # Diameter > 0
rho = sp.Symbol('rho', real=True, positive=True)  # Density > 0
n = sp.Symbol('n', integer=True)  # Integer

# With Greek letters (common in engineering)
alpha, beta, gamma = sp.symbols('alpha beta gamma')
omega, phi, psi = sp.symbols('omega phi psi')

# Subscripted symbols
v_1, v_2 = sp.symbols('v_1 v_2')  # v₁, v₂
```

### Basic Operations

```python
import sympy as sp

x, y = sp.symbols('x y')

# Arithmetic
expr1 = x + y
expr2 = x * y
expr3 = x**2 + 2*x*y + y**2

# Simplification
sp.simplify(x**2 + 2*x*y + y**2)  # (x + y)**2

# Expansion
sp.expand((x + y)**3)  # x**3 + 3*x**2*y + 3*x*y**2 + y**3

# Factorization
sp.factor(x**2 - y**2)  # (x - y)*(x + y)

# Substitution
expr = x**2 + y
expr.subs(x, 2)  # y + 4
expr.subs([(x, 2), (y, 3)])  # 7
```

## Engineering Applications

### 1. Euler Turbine Equation Derivation

Derive the fundamental equation for turbomachinery:

```python
import sympy as sp

# Define symbols
U_1, U_2 = sp.symbols('U_1 U_2', positive=True)  # Peripheral velocities
C_u1, C_u2 = sp.symbols('C_u1 C_u2', real=True)  # Tangential velocities
g = sp.Symbol('g', positive=True)  # Gravity

# Euler turbine equation (energy transfer per unit mass)
# H = (U₂·Cᵤ₂ - U₁·Cᵤ₁) / g
H = (U_2 * C_u2 - U_1 * C_u1) / g

print("Euler Turbine Equation:")
print(f"H = {H}")

# For radial inlet (no prerotation): C_u1 = 0
H_radial = H.subs(C_u1, 0)
print(f"\nWith radial inlet (C_u1 = 0):")
print(f"H = {H_radial}")
print(f"Simplified: H = {sp.simplify(H_radial)}")

# Express in terms of rotational speed and diameter
omega, D_2 = sp.symbols('omega D_2', positive=True)
U_2_expr = omega * D_2 / 2  # U = ω·D/2

# Velocity triangle relations
beta_2, C_m2 = sp.symbols('beta_2 C_m2', positive=True)
C_u2_expr = U_2_expr - C_m2 / sp.tan(beta_2)

# Substitute
H_expanded = H_radial.subs(U_2, U_2_expr).subs(C_u2, C_u2_expr)
print(f"\nIn terms of omega, D, and velocity triangle:")
print(f"H = {H_expanded}")

# Generate LaTeX for documentation
print(f"\nLaTeX: {sp.latex(H_expanded)}")
```

### 2. Bernoulli Equation Solving

Solve Bernoulli equation for different unknowns:

```python
import sympy as sp

# Define symbols
p_1, p_2 = sp.symbols('p_1 p_2', real=True)  # Pressures
v_1, v_2 = sp.symbols('v_1 v_2', positive=True)  # Velocities
z_1, z_2 = sp.symbols('z_1 z_2', real=True)  # Elevations
rho = sp.Symbol('rho', positive=True)  # Density
g = sp.Symbol('g', positive=True)  # Gravity

# Bernoulli equation: p₁/ρg + v₁²/2g + z₁ = p₂/ρg + v₂²/2g + z₂
bernoulli = sp.Eq(
    p_1/(rho*g) + v_1**2/(2*g) + z_1,
    p_2/(rho*g) + v_2**2/(2*g) + z_2
)

print("Bernoulli Equation:")
print(bernoulli)

# Solve for velocity v_2
v_2_solution = sp.solve(bernoulli, v_2)
print(f"\nSolving for v_2:")
for sol in v_2_solution:
    print(f"v_2 = {sol}")

# Solve for pressure p_2
p_2_solution = sp.solve(bernoulli, p_2)[0]
print(f"\nSolving for p_2:")
print(f"p_2 = {p_2_solution}")
print(f"Simplified: p_2 = {sp.simplify(p_2_solution)}")

# Special case: horizontal pipe (z_1 = z_2), static inlet (v_1 = 0)
p_2_horizontal = p_2_solution.subs([(z_1, z_2), (v_1, 0)])
print(f"\nHorizontal pipe, static inlet:")
print(f"p_2 = {sp.simplify(p_2_horizontal)}")
```

### 3. Simplifying Complex Expressions

Simplify head loss equations:

```python
import sympy as sp

# Define symbols
f, L, D, V, g = sp.symbols('f L D V g', positive=True)
epsilon, Re = sp.symbols('epsilon Re', positive=True)

# Darcy-Weisbach equation
h_f = f * (L / D) * (V**2 / (2*g))

print("Darcy-Weisbach head loss:")
print(f"h_f = {h_f}")

# Combine with Hagen-Poiseuille for laminar flow (f = 64/Re)
f_laminar = 64 / Re
h_f_laminar = h_f.subs(f, f_laminar)
print(f"\nLaminar flow (f = 64/Re):")
print(f"h_f = {sp.simplify(h_f_laminar)}")

# Reynolds number: Re = ρVD/μ
mu = sp.Symbol('mu', positive=True)
Re_expr = (sp.Symbol('rho', positive=True) * V * D) / mu
h_f_expanded = h_f_laminar.subs(Re, Re_expr)
print(f"\nExpanded (Re = ρVD/μ):")
print(f"h_f = {sp.simplify(h_f_expanded)}")

# Verify dimensions
print(f"\nVerify: All terms should have dimension [length]")
```

### 4. Solving Systems of Equations

Solve pipe network equations:

```python
import sympy as sp

# Three-pipe junction
Q_1, Q_2, Q_3 = sp.symbols('Q_1 Q_2 Q_3', real=True)
h_1, h_2, h_3 = sp.symbols('h_1 h_2 h_3', positive=True)
K_1, K_2, K_3 = sp.symbols('K_1 K_2 K_3', positive=True)

# Conservation of mass at junction
mass_balance = sp.Eq(Q_1, Q_2 + Q_3)

# Head loss equations: h = K·Q²
head_loss_1 = sp.Eq(h_1, K_1 * Q_1**2)
head_loss_2 = sp.Eq(h_2, K_2 * Q_2**2)
head_loss_3 = sp.Eq(h_3, K_3 * Q_3**2)

# Loop equation: h_2 = h_3 (parallel pipes)
loop_equation = sp.Eq(h_2, h_3)

print("Pipe Network Equations:")
print(f"Mass balance: {mass_balance}")
print(f"Head loss 1: {head_loss_1}")
print(f"Head loss 2: {head_loss_2}")
print(f"Head loss 3: {head_loss_3}")
print(f"Loop: {loop_equation}")

# Solve for Q_2 and Q_3 in terms of Q_1
solution = sp.solve([mass_balance, head_loss_2, head_loss_3, loop_equation],
                   [Q_2, Q_3, h_2, h_3])

print("\nSolution:")
for var, expr in solution.items():
    print(f"{var} = {expr}")
```

### 5. Calculus: Derivatives for Optimization

Find optimal pump operating point:

```python
import sympy as sp

# Define symbols
Q = sp.Symbol('Q', positive=True)
a, b, c = sp.symbols('a b c', real=True)
rho, g = sp.symbols('rho g', positive=True)

# Pump curve: H = a - b·Q - c·Q²
H = a - b*Q - c*Q**2

# System curve: H_sys = K·Q²
K = sp.Symbol('K', positive=True)
H_sys = K * Q**2

# Operating point: H_pump = H_sys
operating_point = sp.Eq(H, H_sys)

# Solve for Q
Q_operating = sp.solve(operating_point, Q)
print("Operating Point:")
print(f"Q = {Q_operating}")

# Power: P = ρ·g·Q·H
P = rho * g * Q * H

# Find maximum efficiency point (dP/dQ = 0)
dP_dQ = sp.diff(P, Q)
print(f"\nPower derivative:")
print(f"dP/dQ = {dP_dQ}")

# Critical points
critical_points = sp.solve(dP_dQ, Q)
print(f"\nCritical points:")
for cp in critical_points:
    print(f"Q = {cp}")

# Second derivative test
d2P_dQ2 = sp.diff(dP_dQ, Q)
print(f"\nSecond derivative:")
print(f"d²P/dQ² = {d2P_dQ2}")
```

### 6. Calculus: Integrals for Flow Calculations

Calculate average velocity from profile:

```python
import sympy as sp

# Define symbols
r, R = sp.symbols('r R', positive=True)
u_max = sp.Symbol('u_max', positive=True)

# Laminar velocity profile in circular pipe
# u(r) = u_max·(1 - (r/R)²)
u = u_max * (1 - (r/R)**2)

print("Velocity profile:")
print(f"u(r) = {u}")

# Flow rate: Q = ∫∫ u dA = ∫₀ᴿ u(r)·2πr dr
integrand = u * 2 * sp.pi * r
print(f"\nIntegrand: {integrand}")

# Integrate from 0 to R
Q = sp.integrate(integrand, (r, 0, R))
print(f"\nFlow rate Q = {Q}")
print(f"Simplified: Q = {sp.simplify(Q)}")

# Average velocity: V_avg = Q / A
A = sp.pi * R**2
V_avg = Q / A
print(f"\nAverage velocity:")
print(f"V_avg = {sp.simplify(V_avg)}")
print(f"Ratio: V_avg/u_max = {sp.simplify(V_avg/u_max)}")
```

### 7. Unit Verification and Dimensional Analysis

Verify dimensional consistency:

```python
import sympy as sp

# Define dimensions as symbols
M, L, T = sp.symbols('M L T')  # Mass, Length, Time

# Define physical quantities with dimensions
quantities = {
    'rho': M / L**3,      # Density [kg/m³]
    'V': L / T,           # Velocity [m/s]
    'D': L,               # Diameter [m]
    'mu': M / (L * T),    # Dynamic viscosity [Pa·s]
    'g': L / T**2,        # Gravity [m/s²]
    'p': M / (L * T**2),  # Pressure [Pa]
}

# Reynolds number: Re = ρVD/μ
Re_dimensions = (quantities['rho'] * quantities['V'] * quantities['D']) / quantities['mu']
Re_simplified = sp.simplify(Re_dimensions)

print("Reynolds Number Dimensional Analysis:")
print(f"Re = (ρVD)/μ = {Re_simplified}")
print(f"Dimensionless: {Re_simplified == 1}")

# Darcy-Weisbach: Δp = f·(L/D)·(ρV²/2)
f = 1  # Friction factor (dimensionless)
dp_dimensions = f * (L/L) * (quantities['rho'] * quantities['V']**2)
dp_simplified = sp.simplify(dp_dimensions)

print(f"\nDarcy-Weisbach pressure drop:")
print(f"Δp dimensions = {dp_simplified}")
print(f"Expected: {quantities['p']}")
print(f"Match: {sp.simplify(dp_simplified - quantities['p']) == 0}")

# Head (energy per unit weight): H = p/(ρg)
H_dimensions = quantities['p'] / (quantities['rho'] * quantities['g'])
H_simplified = sp.simplify(H_dimensions)

print(f"\nHead (pressure/density/gravity):")
print(f"H dimensions = {H_simplified}")
print(f"Has dimension [length]: {H_simplified == L}")
```

### 8. Generate LaTeX for Documentation

Create publication-ready equations:

```python
import sympy as sp

# Define symbols with proper formatting
rho = sp.Symbol('rho')
omega = sp.Symbol('omega')
D = sp.Symbol('D')
Q = sp.Symbol('Q')
H = sp.Symbol('H')
P = sp.Symbol('P')
eta = sp.Symbol('eta')
g = sp.Symbol('g')

# Pump equations
pump_head = (omega * D)**2 / (8 * g)
pump_power = rho * g * Q * H / eta
specific_speed = omega * sp.sqrt(Q) / H**(sp.Rational(3,4))

print("Pump Equations in LaTeX:")
print("\nHead:")
print(sp.latex(sp.Eq(H, pump_head)))

print("\nPower:")
print(sp.latex(sp.Eq(P, pump_power)))

print("\nSpecific Speed:")
print(sp.latex(sp.Eq(sp.Symbol('N_s'), specific_speed)))

# For markdown documents
print("\n\nMarkdown format:")
print(f"$$H = {sp.latex(pump_head)}$$")
print(f"$$P = {sp.latex(pump_power)}$$")
print(f"$$N_s = {sp.latex(specific_speed)}$$")
```

## Advanced Features

### Matrix Operations

Solve linear systems symbolically:

```python
import sympy as sp

# Define matrix elements
A = sp.Matrix([
    [1, 2, 3],
    [4, 5, 6],
    [7, 8, 10]
])

b = sp.Matrix([3, 6, 9])

# Solve Ax = b
x = A.solve(b)
print(f"Solution: x = {x}")

# Verify
print(f"Verification: A·x = {A * x}")

# Determinant
print(f"Determinant: det(A) = {A.det()}")

# Inverse
A_inv = A.inv()
print(f"Inverse: A⁻¹ = {A_inv}")

# Eigenvalues
eigenvals = A.eigenvals()
print(f"Eigenvalues: {eigenvals}")
```

### Piecewise Functions

Model pump curves with different regimes:

```python
import sympy as sp

Q = sp.Symbol('Q', positive=True)

# Pump curve with different behaviors
H = sp.Piecewise(
    (85 - 10*Q, Q <= 3),         # Low flow (stable)
    (100 - 20*Q, Q <= 5),        # Mid flow
    (50, Q > 5)                   # High flow (surge)
)

print("Piecewise pump curve:")
print(H)

# Evaluate at different flows
for q in [2, 4, 6]:
    print(f"H(Q={q}) = {H.subs(Q, q)}")

# Derivative (may be discontinuous)
dH_dQ = sp.diff(H, Q)
print(f"\ndH/dQ = {dH_dQ}")
```

### Series Expansions

Taylor series for approximations:

```python
import sympy as sp

x = sp.Symbol('x')

# Original function
f = sp.exp(x)

# Taylor series expansion around x=0
taylor = sp.series(f, x, 0, n=5)
print(f"Taylor series of e^x:")
print(taylor)

# For engineering: small angle approximation
theta = sp.Symbol('theta')
sin_approx = sp.series(sp.sin(theta), theta, 0, n=4)
print(f"\nSmall angle: sin(θ) ≈ {sin_approx.removeO()}")

cos_approx = sp.series(sp.cos(theta), theta, 0, n=4)
print(f"Small angle: cos(θ) ≈ {cos_approx.removeO()}")
```

## Best Practices

1. **Use assumptions** - Declare if variables are positive, real, integer, etc.
   ```python
   D = sp.Symbol('D', positive=True)  # Diameter is always positive
   n = sp.Symbol('n', integer=True)   # RPM is integer
   ```

2. **Simplify expressions** - Always simplify before substituting numbers
   ```python
   expr = sp.simplify(complex_expr)  # Symbolic simplification
   result = expr.subs(values)        # Then substitute
   ```

3. **Use rational numbers** - Avoid floating point in symbolic math
   ```python
   # Good
   expr = sp.Rational(1, 3) * x**2

   # Bad
   expr = 0.333333 * x**2  # Loses exactness
   ```

4. **Check dimensions** - Verify equations are dimensionally consistent
   ```python
   # Define dimension variables and verify
   ```

5. **Generate code** - Convert symbolic expressions to NumPy functions
   ```python
   # Symbolic
   expr = x**2 + sp.sin(x)

   # Convert to numerical function
   f = sp.lambdify(x, expr, 'numpy')

   # Use with NumPy arrays
   import numpy as np
   x_vals = np.linspace(0, 10, 100)
   y_vals = f(x_vals)  # Fast numerical evaluation
   ```

6. **Document with LaTeX** - Export equations for reports
   ```python
   print(sp.latex(equation))  # Copy to LaTeX document
   ```

## Integration with Engineering Workflows

SymPy integrates well with:

- **NumPy**: Convert symbolic expressions to numerical functions
- **Matplotlib**: Plot symbolic expressions
- **Pandas**: Generate formula columns
- **Jupyter**: Display formatted equations
- **Pint**: Verify unit consistency symbolically

Example workflow:

```python
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# 1. Derive equation symbolically
Q, H = sp.symbols('Q H')
a, b, c = 80, 50, 500
H_expr = a - b*Q - c*Q**2

# 2. Simplify and verify
print(f"Pump curve: H = {H_expr}")

# 3. Convert to numerical function
H_func = sp.lambdify(Q, H_expr, 'numpy')

# 4. Evaluate numerically
Q_vals = np.linspace(0, 0.3, 100)
H_vals = H_func(Q_vals)

# 5. Plot
plt.plot(Q_vals, H_vals)
plt.xlabel('Flow rate (m³/s)')
plt.ylabel('Head (m)')
plt.title(f'Pump Curve: $H = {sp.latex(H_expr)}$')
plt.grid(True)
plt.show()
```

## Common Pitfalls

1. **Mixing symbolic and numeric** - Keep them separate
   ```python
   # Bad
   x = sp.Symbol('x')
   result = x + 3.14159  # Mixes symbolic and numeric

   # Good
   result = x + sp.pi    # Both symbolic
   ```

2. **Not simplifying** - Always simplify complex expressions
   ```python
   expr = complex_calculation()
   expr = sp.simplify(expr)  # Much cleaner
   ```

3. **Assuming evaluation** - Symbols don't auto-substitute
   ```python
   x = sp.Symbol('x')
   expr = x**2
   # expr is still symbolic, not a number!
   value = expr.subs(x, 5)  # Now it's 25
   ```

4. **Order of operations** - Be explicit with parentheses
   ```python
   # Ambiguous
   expr = a + b / c * d

   # Clear
   expr = a + (b / c) * d
   expr = a + b / (c * d)
   ```

## References

- SymPy Documentation: https://docs.sympy.org/
- SymPy Tutorial: https://docs.sympy.org/latest/tutorial/
- SymPy Physics Module: https://docs.sympy.org/latest/modules/physics/
- Symbolic Math in Engineering: https://www.sympy.org/scipy-2017-codegen-tutorial/

## Further Learning

- **Intro Tutorial**: https://docs.sympy.org/latest/tutorial/intro.html
- **Solving Equations**: https://docs.sympy.org/latest/tutorial/solvers.html
- **Calculus**: https://docs.sympy.org/latest/tutorial/calculus.html
- **Matrices**: https://docs.sympy.org/latest/tutorial/matrices.html
- **Code Generation**: https://docs.sympy.org/latest/modules/utilities/lambdify.html
