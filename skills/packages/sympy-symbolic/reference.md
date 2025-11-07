# SymPy Symbolic Math Reference Guide

## Quick Symbol Creation

```python
import sympy as sp

# Basic symbols
x = sp.Symbol('x')
x, y, z = sp.symbols('x y z')

# With assumptions
D = sp.Symbol('D', positive=True)      # D > 0
n = sp.Symbol('n', integer=True)        # n is integer
theta = sp.Symbol('theta', real=True)   # θ is real

# Multiple symbols with same assumption
a, b, c = sp.symbols('a b c', positive=True, real=True)
```

## Common Engineering Symbols

### Fluid Dynamics

```python
import sympy as sp

# Pressures and forces
p = sp.Symbol('p', real=True)              # Pressure [Pa]
p_1, p_2 = sp.symbols('p_1 p_2', real=True)  # Inlet/outlet pressures
F = sp.Symbol('F', real=True)              # Force [N]
tau = sp.Symbol('tau', real=True)          # Shear stress [Pa]

# Velocities and flow rates
v, V = sp.symbols('v V', positive=True)    # Velocity [m/s]
u = sp.Symbol('u', real=True)              # x-velocity [m/s]
v = sp.Symbol('v', real=True)              # y-velocity [m/s]
w = sp.Symbol('w', real=True)              # z-velocity [m/s]
Q = sp.Symbol('Q', positive=True)          # Flow rate [m³/s]
V_dot = sp.Symbol('V_dot', positive=True)  # Volume flow rate

# Geometric parameters
D = sp.Symbol('D', positive=True)          # Diameter [m]
r, R = sp.symbols('r R', positive=True)    # Radius [m]
L = sp.Symbol('L', positive=True)          # Length [m]
A = sp.Symbol('A', positive=True)          # Area [m²]
z = sp.Symbol('z', real=True)              # Elevation [m]

# Fluid properties
rho = sp.Symbol('rho', positive=True)      # Density [kg/m³]
mu = sp.Symbol('mu', positive=True)        # Dynamic viscosity [Pa·s]
nu = sp.Symbol('nu', positive=True)        # Kinematic viscosity [m²/s]

# Dimensionless numbers
Re = sp.Symbol('Re', positive=True)        # Reynolds number
Fr = sp.Symbol('Fr', positive=True)        # Froude number
We = sp.Symbol('We', positive=True)        # Weber number
Ma = sp.Symbol('Ma', positive=True)        # Mach number

# Constants
g = sp.Symbol('g', positive=True)          # Gravity [m/s²] = 9.81
```

### Turbomachinery

```python
import sympy as sp

# Rotational parameters
omega = sp.Symbol('omega', positive=True)  # Angular velocity [rad/s]
N = sp.Symbol('N', positive=True)          # Rotational speed [RPM]
n = sp.Symbol('n', positive=True)          # Specific speed

# Velocities
U = sp.Symbol('U', positive=True)          # Peripheral velocity [m/s]
C = sp.Symbol('C', positive=True)          # Absolute velocity [m/s]
W = sp.Symbol('W', positive=True)          # Relative velocity [m/s]
C_u = sp.Symbol('C_u', real=True)          # Tangential component [m/s]
C_m = sp.Symbol('C_m', positive=True)      # Meridional component [m/s]

# Angles
alpha = sp.Symbol('alpha', real=True)      # Absolute flow angle [rad/deg]
beta = sp.Symbol('beta', real=True)        # Relative flow angle [rad/deg]

# Dimensions
D_1, D_2 = sp.symbols('D_1 D_2', positive=True)  # Inlet/outlet diameter
b_1, b_2 = sp.symbols('b_1 b_2', positive=True)  # Inlet/outlet width
Z = sp.Symbol('Z', integer=True, positive=True)  # Number of blades

# Performance
H = sp.Symbol('H', positive=True)          # Head [m]
P = sp.Symbol('P', positive=True)          # Power [W]
T = sp.Symbol('T', positive=True)          # Torque [N·m]
eta = sp.Symbol('eta', positive=True)      # Efficiency [0-1]

# Coefficients
phi = sp.Symbol('phi', positive=True)      # Flow coefficient
psi = sp.Symbol('psi', positive=True)      # Head coefficient
lambda_c = sp.Symbol('lambda', positive=True)  # Power coefficient
```

### Thermodynamics

```python
import sympy as sp

# Temperature
T = sp.Symbol('T', positive=True)          # Temperature [K]
T_1, T_2 = sp.symbols('T_1 T_2', positive=True)

# Energy
E = sp.Symbol('E', real=True)              # Energy [J]
h = sp.Symbol('h', real=True)              # Enthalpy [J/kg]
u = sp.Symbol('u', real=True)              # Internal energy [J/kg]
s = sp.Symbol('s', real=True)              # Entropy [J/(kg·K)]

# Heat and work
Q = sp.Symbol('Q', real=True)              # Heat [J]
W = sp.Symbol('W', real=True)              # Work [J]
Q_dot = sp.Symbol('Q_dot', real=True)      # Heat rate [W]
W_dot = sp.Symbol('W_dot', real=True)      # Power [W]

# Gas properties
R = sp.Symbol('R', positive=True)          # Gas constant [J/(kg·K)]
gamma = sp.Symbol('gamma', positive=True)  # Specific heat ratio
c_p = sp.Symbol('c_p', positive=True)      # Specific heat const P
c_v = sp.Symbol('c_v', positive=True)      # Specific heat const V
```

### Greek Letters (Common in Engineering)

```python
import sympy as sp

# Lowercase
alpha = sp.Symbol('alpha')     # α - angle, thermal diffusivity
beta = sp.Symbol('beta')       # β - angle, compressibility
gamma = sp.Symbol('gamma')     # γ - specific heat ratio
delta = sp.Symbol('delta')     # δ - small change
epsilon = sp.Symbol('epsilon') # ε - roughness, strain
zeta = sp.Symbol('zeta')       # ζ - loss coefficient
eta = sp.Symbol('eta')         # η - efficiency
theta = sp.Symbol('theta')     # θ - angle
kappa = sp.Symbol('kappa')     # κ - thermal conductivity
lambda_ = sp.Symbol('lambda')  # λ - wavelength, friction factor
mu = sp.Symbol('mu')           # μ - dynamic viscosity
nu = sp.Symbol('nu')           # ν - kinematic viscosity
xi = sp.Symbol('xi')           # ξ - coordinate
rho = sp.Symbol('rho')         # ρ - density
sigma = sp.Symbol('sigma')     # σ - stress, surface tension
tau = sp.Symbol('tau')         # τ - shear stress, time constant
phi = sp.Symbol('phi')         # φ - flow coefficient, angle
chi = sp.Symbol('chi')         # χ - quality
psi = sp.Symbol('psi')         # ψ - head coefficient, stream function
omega = sp.Symbol('omega')     # ω - angular velocity

# Uppercase
Delta = sp.Symbol('Delta')     # Δ - change, difference
Phi = sp.Symbol('Phi')         # Φ - potential
Psi = sp.Symbol('Psi')         # Ψ - stream function
Omega = sp.Symbol('Omega')     # Ω - solid angle
```

## Essential SymPy Functions

### Simplification

```python
import sympy as sp

expr = x**2 + 2*x*y + y**2

sp.simplify(expr)           # General simplification
sp.expand(expr)             # Expand products
sp.factor(expr)             # Factorize
sp.collect(expr, x)         # Collect terms by x
sp.cancel(expr)             # Cancel common factors
sp.trigsimp(expr)           # Simplify trig expressions
sp.powsimp(expr)            # Simplify powers
sp.radsimp(expr)            # Rationalize denominators
```

### Solving Equations

```python
import sympy as sp

# Single equation
sp.solve(x**2 - 4, x)                    # Returns [2, -2]

# With equation object
eq = sp.Eq(x**2 + y, 10)
sp.solve(eq, y)                          # Solve for y

# System of equations
eq1 = sp.Eq(x + y, 5)
eq2 = sp.Eq(x - y, 1)
sp.solve([eq1, eq2], [x, y])             # {x: 3, y: 2}

# Nonlinear system
sp.solve([x**2 + y**2 - 1, x - y], [x, y])

# With domain restrictions
sp.solve(x**2 - 4, x, domain=sp.S.Reals)  # Real solutions only
```

### Calculus

```python
import sympy as sp

x = sp.Symbol('x')

# Derivatives
sp.diff(x**2, x)                         # 2*x
sp.diff(x**3*sp.sin(x), x)               # 3*x**2*sin(x) + x**3*cos(x)
sp.diff(x**2 + y**2, x, y)               # Mixed partial: 0

# Multiple derivatives
sp.diff(x**4, x, 3)                      # Third derivative: 24*x

# Integrals
sp.integrate(x**2, x)                    # x**3/3
sp.integrate(x**2, (x, 0, 1))            # Definite: 1/3
sp.integrate(sp.sin(x), (x, 0, sp.pi))   # 2

# Multiple integrals
sp.integrate(x*y, (x, 0, 1), (y, 0, 1))  # 1/4

# Limits
sp.limit(sp.sin(x)/x, x, 0)              # 1
sp.limit(1/x, x, sp.oo)                  # 0

# Series expansion
sp.series(sp.exp(x), x, 0, n=5)          # 1 + x + x²/2 + x³/6 + x⁴/24 + O(x⁵)
```

### Substitution

```python
import sympy as sp

expr = x**2 + 2*x*y + y**2

# Single substitution
expr.subs(x, 3)                          # y**2 + 6*y + 9

# Multiple substitutions
expr.subs([(x, 2), (y, 3)])              # 25

# Substitution with expression
expr.subs(x, sp.sqrt(z))                 # z + 2*sqrt(z)*y + y**2

# Replace sub-expression
expr.subs(x + y, z)                      # z**2
```

### Matrix Operations

```python
import sympy as sp

# Create matrix
A = sp.Matrix([[1, 2], [3, 4]])
B = sp.Matrix([[5, 6], [7, 8]])

# Basic operations
A + B                                    # Matrix addition
A * B                                    # Matrix multiplication
A**2                                     # Matrix power
A.T                                      # Transpose

# Properties
A.det()                                  # Determinant
A.inv()                                  # Inverse
A.eigenvals()                            # Eigenvalues
A.eigenvects()                           # Eigenvectors

# Solving linear systems
# Ax = b
x = A.solve(b)                           # Solution vector

# Row operations
A.rref()                                 # Row reduced echelon form
```

### Special Functions

```python
import sympy as sp

# Trigonometric
sp.sin(x), sp.cos(x), sp.tan(x)
sp.asin(x), sp.acos(x), sp.atan(x)
sp.sinh(x), sp.cosh(x), sp.tanh(x)

# Exponential and logarithmic
sp.exp(x)                                # e^x
sp.log(x)                                # Natural log
sp.log(x, 10)                            # Log base 10

# Square root and powers
sp.sqrt(x)                               # √x
sp.cbrt(x)                               # ∛x
sp.Pow(x, sp.Rational(1, 3))             # x^(1/3)

# Absolute value and sign
sp.Abs(x)                                # |x|
sp.sign(x)                               # Sign function

# Rounding
sp.floor(x)                              # Floor
sp.ceiling(x)                            # Ceiling

# Piecewise functions
sp.Piecewise((x, x >= 0), (-x, x < 0))   # |x|

# Special numbers
sp.pi                                    # π
sp.E                                     # e
sp.I                                     # √(-1)
sp.oo                                    # ∞
```

## Converting to Numerical Functions

### Using lambdify

```python
import sympy as sp
import numpy as np

# Define symbolic expression
x = sp.Symbol('x')
expr = x**2 + sp.sin(x)

# Convert to numerical function
f = sp.lambdify(x, expr, 'numpy')

# Use with NumPy
x_vals = np.linspace(0, 10, 100)
y_vals = f(x_vals)  # Fast numerical evaluation

# Multiple variables
x, y = sp.symbols('x y')
expr = x**2 + y**2
f = sp.lambdify([x, y], expr, 'numpy')
result = f(3, 4)  # 25

# With modules
expr = x**2 * sp.exp(-x)
f_numpy = sp.lambdify(x, expr, 'numpy')    # Uses numpy.exp
f_math = sp.lambdify(x, expr, 'math')      # Uses math.exp
```

## LaTeX Export

```python
import sympy as sp

x, y = sp.symbols('x y')
expr = x**2 + 2*x*y + y**2

# Generate LaTeX
latex_str = sp.latex(expr)
print(latex_str)  # x^{2} + 2 x y + y^{2}

# For equations
eq = sp.Eq(sp.Symbol('H'), expr)
print(sp.latex(eq))  # H = x^{2} + 2 x y + y^{2}

# Pretty printing in terminal
sp.pprint(expr)
# 2          2
# x  + 2⋅x⋅y + y

# For Jupyter notebooks (automatic rendering)
from IPython.display import display
display(expr)  # Renders as formatted math
```

## Dimensional Analysis Template

```python
import sympy as sp

# Define dimension symbols
M, L, T = sp.symbols('M L T')  # Mass, Length, Time

# Define quantities with dimensions
dimensions = {
    'pressure': M / (L * T**2),      # [Pa] = kg/(m·s²)
    'velocity': L / T,                # [m/s]
    'density': M / L**3,              # [kg/m³]
    'viscosity': M / (L * T),         # [Pa·s]
    'force': M * L / T**2,            # [N]
    'energy': M * L**2 / T**2,        # [J]
    'power': M * L**2 / T**3,         # [W]
    'area': L**2,                     # [m²]
    'volume': L**3,                   # [m³]
    'flow_rate': L**3 / T,            # [m³/s]
}

# Check if expression is dimensionless
def check_dimensionless(expr_dimensions):
    """Verify if expression simplifies to 1 (dimensionless)"""
    simplified = sp.simplify(expr_dimensions)
    return simplified == 1

# Example: Reynolds number
Re_dim = (dimensions['density'] * dimensions['velocity'] * L) / dimensions['viscosity']
print(f"Re dimensionless: {check_dimensionless(Re_dim)}")
```

## Common Engineering Equations

### Fluid Mechanics

```python
import sympy as sp

# Define symbols
p, v, z, rho, g, H = sp.symbols('p v z rho g H', real=True)
D, L, f, Q, A = sp.symbols('D L f Q A', positive=True)

# Bernoulli equation
bernoulli = p/(rho*g) + v**2/(2*g) + z  # = constant

# Continuity equation
continuity = sp.Eq(Q, A * v)  # Q = A·v

# Darcy-Weisbach equation
head_loss = f * (L/D) * (v**2/(2*g))

# Reynolds number
mu = sp.Symbol('mu', positive=True)
Re = (rho * v * D) / mu

# Hagen-Poiseuille (laminar flow)
Q_laminar = (sp.pi * D**4 * (p_1 - p_2)) / (128 * mu * L)
```

### Pump Equations

```python
import sympy as sp

# Symbols
Q, H, P, eta = sp.symbols('Q H P eta', positive=True)
rho, g = sp.symbols('rho g', positive=True)
omega, D = sp.symbols('omega D', positive=True)

# Hydraulic power
P_hydraulic = rho * g * Q * H

# Shaft power
P_shaft = P_hydraulic / eta

# Specific speed
N_s = omega * sp.sqrt(Q) / H**sp.Rational(3, 4)

# Affinity laws
# Q ∝ n·D³
# H ∝ n²·D²
# P ∝ n³·D⁵
```

## Quick Reference Card

```python
# Basic operations
sp.Symbol('x')                           # Create symbol
sp.symbols('x y z')                      # Multiple symbols
sp.simplify(expr)                        # Simplify
sp.expand(expr)                          # Expand
sp.factor(expr)                          # Factor

# Solving
sp.solve(eq, x)                          # Solve for x
sp.solve([eq1, eq2], [x, y])             # System of equations

# Calculus
sp.diff(expr, x)                         # Derivative
sp.integrate(expr, x)                    # Integral
sp.integrate(expr, (x, a, b))            # Definite integral
sp.limit(expr, x, a)                     # Limit

# Substitution
expr.subs(x, value)                      # Substitute
expr.subs([(x, 1), (y, 2)])              # Multiple

# Conversion
sp.lambdify(x, expr, 'numpy')            # To numerical function
sp.latex(expr)                           # To LaTeX
sp.pprint(expr)                          # Pretty print

# Special values
sp.pi, sp.E, sp.I, sp.oo                 # π, e, i, ∞
sp.Rational(1, 3)                        # Exact fraction 1/3
```

## SymPy Resources

### Official Documentation

- **Main Documentation**: https://docs.sympy.org/
- **Tutorial**: https://docs.sympy.org/latest/tutorial/
- **API Reference**: https://docs.sympy.org/latest/reference/
- **Examples**: https://docs.sympy.org/latest/modules/index.html

### Specific Topics

- **Solvers**: https://docs.sympy.org/latest/tutorial/solvers.html
- **Calculus**: https://docs.sympy.org/latest/tutorial/calculus.html
- **Matrices**: https://docs.sympy.org/latest/tutorial/matrices.html
- **Simplification**: https://docs.sympy.org/latest/tutorial/simplification.html
- **Physics Module**: https://docs.sympy.org/latest/modules/physics/

### Advanced Features

- **Code Generation**: https://docs.sympy.org/latest/modules/utilities/lambdify.html
- **Assumptions**: https://docs.sympy.org/latest/guides/assumptions.html
- **Units**: https://docs.sympy.org/latest/modules/physics/units/
- **LaTeX Printing**: https://docs.sympy.org/latest/modules/printing.html

### Learning Resources

- **SymPy Tutorial (SciPy)**: https://scipy-lectures.org/packages/sympy.html
- **SymPy GitHub**: https://github.com/sympy/sympy
- **Stack Overflow**: Tag `sympy`
- **SymPy Live**: https://live.sympy.org/ (try in browser)

### Books and Papers

- "Instant SymPy Starter" by Ronan Lamy
- "SymPy Documentation" (comprehensive online manual)
- SymPy paper: Meurer et al., "SymPy: symbolic computing in Python", PeerJ Computer Science (2017)

### Engineering Applications

- Fluid mechanics symbolic derivations
- Control systems transfer functions
- Thermodynamic cycle analysis
- Structural mechanics equations
- Heat transfer equations
- Chemical reaction kinetics

### Related Tools

- **Mathematica**: Commercial alternative
- **Maxima**: Open-source CAS
- **SageMath**: Python-based mathematics system
- **SymEngine**: Fast symbolic manipulation engine (C++ backend for SymPy)

---

*This reference guide provides common symbols and quick access to SymPy functions for engineering applications. For detailed examples, see SKILL.md and examples.py.*
