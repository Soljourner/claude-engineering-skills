# Python Engineering Packages Guide

A comprehensive guide to essential Python packages for engineering applications, organized by domain with practical examples and best practices.

---

## 1. Numerical Computing

### NumPy

**Installation:**
```bash
pip install numpy
```

**Key Functions and Classes:**
- `numpy.array()` - Core N-dimensional array object
- `numpy.linspace()`, `numpy.arange()` - Array generation
- `numpy.linalg` - Linear algebra operations
- `numpy.fft` - Fast Fourier Transform
- `numpy.polynomial` - Polynomial operations
- Broadcasting, vectorization, and universal functions (ufuncs)

**Engineering Applications:**
- Matrix operations for structural analysis
- Signal processing and frequency analysis
- Statistical analysis of experimental data
- Numerical integration and differentiation
- Solving systems of linear equations

**Code Examples:**

```python
import numpy as np

# Example 1: Structural analysis - solving system of equations
# For a truss structure: [K]{d} = {F}
K = np.array([[100, -50, 0],
              [-50, 150, -100],
              [0, -100, 100]])  # Stiffness matrix (kN/m)
F = np.array([0, 1000, 0])  # Force vector (N)

displacements = np.linalg.solve(K, F)
print(f"Node displacements: {displacements} m")

# Example 2: FFT analysis for vibration data
time = np.linspace(0, 1, 1000)
signal = 2 * np.sin(2 * np.pi * 50 * time) + 0.5 * np.sin(2 * np.pi * 120 * time)
noise = 0.2 * np.random.randn(len(time))
signal_noisy = signal + noise

fft_result = np.fft.fft(signal_noisy)
frequencies = np.fft.fftfreq(len(time), time[1] - time[0])
dominant_freq = frequencies[np.argmax(np.abs(fft_result[:len(frequencies)//2]))]
print(f"Dominant frequency: {dominant_freq} Hz")

# Example 3: Stress tensor transformation
stress_tensor = np.array([[100, 30, 0],
                          [30, 60, 0],
                          [0, 0, 40]])  # MPa
eigenvalues, eigenvectors = np.linalg.eig(stress_tensor)
principal_stresses = np.sort(eigenvalues)[::-1]
print(f"Principal stresses: {principal_stresses} MPa")
```

**Documentation:**
- https://numpy.org/doc/stable/

**Common Gotchas:**
- Array copies vs. views: Use `.copy()` when needed to avoid unintended modifications
- Integer division: Use floats for precise calculations
- Broadcasting rules can be confusing with mismatched dimensions
- Matrix multiplication: Use `@` operator or `np.matmul()`, not `*`

---

### SciPy

**Installation:**
```bash
pip install scipy
```

**Key Functions and Classes:**
- `scipy.integrate` - Numerical integration (ODEs, quadrature)
- `scipy.interpolate` - Interpolation methods
- `scipy.signal` - Signal processing
- `scipy.stats` - Statistical functions
- `scipy.sparse` - Sparse matrix operations
- `scipy.linalg` - Extended linear algebra

**Engineering Applications:**
- Solving differential equations (heat transfer, fluid dynamics)
- Curve fitting and data interpolation
- Filter design and signal analysis
- Statistical process control
- Optimization problems

**Code Examples:**

```python
from scipy import integrate, interpolate, signal, stats
import numpy as np

# Example 1: Solving ODE - cooling of an object (Newton's law of cooling)
def cooling_rate(t, T, k, T_ambient):
    """dT/dt = -k(T - T_ambient)"""
    return -k * (T - T_ambient)

k = 0.1  # cooling constant (1/min)
T_ambient = 20  # ambient temperature (°C)
T0 = 100  # initial temperature (°C)
t_span = (0, 50)  # time range (min)
t_eval = np.linspace(0, 50, 100)

sol = integrate.solve_ivp(cooling_rate, t_span, [T0],
                          args=(k, T_ambient), t_eval=t_eval)
print(f"Temperature after 30 min: {np.interp(30, sol.t, sol.y[0]):.2f}°C")

# Example 2: Interpolation for experimental data
# Pressure vs. temperature data from experiment
T_data = np.array([20, 40, 60, 80, 100])  # °C
P_data = np.array([101.3, 107.4, 119.9, 147.1, 181.4])  # kPa

cubic_interp = interpolate.CubicSpline(T_data, P_data)
T_interpolated = np.linspace(20, 100, 50)
P_interpolated = cubic_interp(T_interpolated)
print(f"Interpolated pressure at 75°C: {cubic_interp(75):.2f} kPa")

# Example 3: Low-pass filter for sensor noise reduction
sampling_rate = 1000  # Hz
nyquist = sampling_rate / 2
cutoff_freq = 50  # Hz
order = 4

b, a = signal.butter(order, cutoff_freq / nyquist, btype='low')
# Apply filter to noisy sensor data
time = np.linspace(0, 1, sampling_rate)
noisy_data = np.sin(2 * np.pi * 5 * time) + 0.5 * np.random.randn(len(time))
filtered_data = signal.filtfilt(b, a, noisy_data)
print(f"Signal-to-noise improvement: {np.std(noisy_data) / np.std(filtered_data):.2f}x")
```

**Documentation:**
- https://docs.scipy.org/doc/scipy/

**Common Gotchas:**
- `solve_ivp` requires function signature `f(t, y, *args)`, note the order
- `filtfilt` applies filter twice (forward and backward) for zero phase shift
- Many functions require sorted input arrays
- Sparse matrices have different methods than dense arrays

---

### SymPy

**Installation:**
```bash
pip install sympy
```

**Key Functions and Classes:**
- `sympy.symbols()` - Define symbolic variables
- `sympy.diff()`, `sympy.integrate()` - Calculus operations
- `sympy.solve()`, `sympy.dsolve()` - Equation solving
- `sympy.simplify()`, `sympy.expand()` - Expression manipulation
- `sympy.Matrix` - Symbolic matrices
- `sympy.lambdify()` - Convert to numerical functions

**Engineering Applications:**
- Deriving equations of motion
- Symbolic differentiation for sensitivity analysis
- Exact solutions to differential equations
- Symbolic matrix operations for kinematics
- Automatic code generation from mathematical models

**Code Examples:**

```python
import sympy as sp
import numpy as np

# Example 1: Beam deflection equation derivation
x, L, E, I, w = sp.symbols('x L E I w', positive=True, real=True)

# Simply supported beam with uniform load
# Fourth derivative of deflection equals load/EI
V = w * (L * x / 2 - x**2 / 2)  # Shear force
M = sp.integrate(V, x)  # Bending moment
theta = sp.integrate(M / (E * I), x)  # Slope
y = sp.integrate(theta, x)  # Deflection

# Apply boundary conditions: y(0) = 0, y(L) = 0
C1, C2 = sp.symbols('C1 C2')
y = y + C1 * x + C2
eq1 = y.subs(x, 0)
eq2 = y.subs(x, L)
constants = sp.solve([eq1, eq2], [C1, C2])
y_final = y.subs(constants)
print("Beam deflection equation:")
sp.pprint(sp.simplify(y_final))

# Example 2: Kinematics - velocity and acceleration from position
t = sp.Symbol('t', real=True)
x_pos = 10 * sp.cos(2 * t) + 5  # Position equation (m)
velocity = sp.diff(x_pos, t)
acceleration = sp.diff(velocity, t)

print(f"Velocity: {velocity}")
print(f"Acceleration: {acceleration}")

# Convert to numerical function for plotting
v_func = sp.lambdify(t, velocity, 'numpy')
time_array = np.linspace(0, 2 * np.pi, 100)
velocities = v_func(time_array)

# Example 3: Thermodynamics - Ideal gas law manipulations
P, V, n, R, T = sp.symbols('P V n R T', positive=True)
ideal_gas = sp.Eq(P * V, n * R * T)

# Solve for different variables
T_expr = sp.solve(ideal_gas, T)[0]
print(f"Temperature: {T_expr}")

# Partial derivatives for Maxwell relations
dP_dT = sp.diff(P * V / (n * R), T)
print(f"(∂P/∂T)_V = {dP_dT}")
```

**Documentation:**
- https://docs.sympy.org/

**Common Gotchas:**
- Symbolic computation is much slower than numerical
- Must declare assumptions (positive, real, integer) for proper simplification
- Use `lambdify()` to convert to NumPy functions for numerical evaluation
- `solve()` returns a list, even for single solutions
- Complex expressions may not simplify automatically; try different simplification functions

---

## 2. Fluid Mechanics

### fluids

**Installation:**
```bash
pip install fluids
```

**Key Functions and Classes:**
- `fluids.friction_factor()` - Darcy friction factor calculation
- `fluids.pressure_drop()` - Pipe pressure drop
- `fluids.Reynolds()` - Reynolds number
- `fluids.Cv_to_Kv()`, `fluids.Kv_to_Cv()` - Valve coefficient conversions
- `fluids.drag_sphere()` - Drag force on spheres
- `fluids.pump` - Pump sizing and analysis

**Engineering Applications:**
- Pipeline design and analysis
- Pump and valve sizing
- Flow measurement
- Heat exchanger design
- HVAC system calculations

**Code Examples:**

```python
from fluids import *
import numpy as np

# Example 1: Pipeline pressure drop calculation
L = 100  # pipe length (m)
D = 0.1  # pipe diameter (m)
roughness = 0.000045  # pipe roughness (m) - commercial steel
rho = 1000  # fluid density (kg/m³)
mu = 0.001  # dynamic viscosity (Pa·s)
V = 2.0  # fluid velocity (m/s)

Re = Reynolds(V=V, D=D, rho=rho, mu=mu)
print(f"Reynolds number: {Re:.0f}")

fd = friction_factor(Re=Re, eD=roughness/D)
dP = dP_from_V(V=V, L=L, D=D, rho=rho, fd=fd)
print(f"Pressure drop: {dP/1000:.2f} kPa")

# Example 2: Valve sizing (Cv calculation)
Q = 100  # flow rate (m³/h)
dP_valve = 50000  # pressure drop across valve (Pa)
rho = 1000  # density (kg/m³)
SG = rho / 1000  # specific gravity

Kv = Q / np.sqrt(dP_valve / (1e5 * SG))  # Kv in m³/h
Cv = Kv_to_Cv(Kv)
print(f"Required valve Cv: {Cv:.2f}")

# Example 3: Drag force on a settling particle
D_particle = 0.001  # particle diameter (m)
V_terminal = 0.05  # settling velocity (m/s)
rho_fluid = 1.2  # air density (kg/m³)
mu_fluid = 1.8e-5  # air viscosity (Pa·s)

Re_particle = Reynolds(V=V_terminal, D=D_particle, rho=rho_fluid, mu=mu_fluid)
Cd = drag_sphere(Re_particle)
A_particle = np.pi * (D_particle/2)**2
F_drag = 0.5 * Cd * rho_fluid * V_terminal**2 * A_particle
print(f"Drag force: {F_drag*1e6:.3f} µN")
print(f"Drag coefficient: {Cd:.3f}")
```

**Documentation:**
- https://fluids.readthedocs.io/

**Common Gotchas:**
- Inconsistent units - must use SI units (m, Pa, kg, etc.)
- Some functions have multiple calculation methods; specify method if needed
- Reynolds number transitions affect friction factor significantly
- Valve coefficients (Cv, Kv) have different definitions and units

---

### thermo

**Installation:**
```bash
pip install thermo
```

**Key Functions and Classes:**
- `Chemical` - Access thermodynamic properties
- `Mixture` - Multi-component mixture properties
- `Stream` - Process stream calculations
- Property methods: `Cp`, `H`, `S`, `rho`, etc.
- Flash calculations: `T`, `P`, `VF` (vapor fraction)

**Engineering Applications:**
- Process design and simulation
- Energy balance calculations
- Phase equilibrium
- Property estimation
- Refrigeration cycles

**Code Examples:**

```python
from thermo import Chemical, Mixture, Stream
import numpy as np

# Example 1: Steam properties for power cycle analysis
water = Chemical('water')
water.calculate(T=500+273.15, P=3e6)  # 500°C, 30 bar

print(f"Enthalpy: {water.H/1000:.2f} kJ/kg")
print(f"Entropy: {water.S/1000:.2f} kJ/(kg·K)")
print(f"Density: {water.rho:.2f} kg/m³")

# Isentropic expansion to 0.1 bar
P2 = 1e4  # Pa
water.calculate(S=water.S, P=P2)
h_isentropic = water.H
print(f"Isentropic discharge enthalpy: {h_isentropic/1000:.2f} kJ/kg")

# Example 2: Refrigerant properties for HVAC
refrigerant = Chemical('R134a')
refrigerant.calculate(T=25+273.15, P=1e5)

print(f"R134a at 25°C, 1 bar:")
print(f"  Phase: {refrigerant.phase}")
print(f"  Cp: {refrigerant.Cp:.2f} J/(kg·K)")
print(f"  Thermal conductivity: {refrigerant.k:.4f} W/(m·K)")

# Example 3: Mixture properties for combustion calculations
# Air composition (mole fractions)
air = Mixture(['nitrogen', 'oxygen', 'argon'],
              zs=[0.78, 0.21, 0.01], T=298.15, P=101325)

print(f"Air properties at 25°C:")
print(f"  Molecular weight: {air.MW:.2f} g/mol")
print(f"  Density: {air.rho:.3f} kg/m³")
print(f"  Cp: {air.Cp:.2f} J/(kg·K)")
print(f"  Viscosity: {air.mu*1e6:.2f} µPa·s")
```

**Documentation:**
- https://thermo.readthedocs.io/

**Common Gotchas:**
- Temperature must be in Kelvin, pressure in Pascals
- Must call `.calculate()` with state variables before accessing properties
- Some chemicals have limited data; check availability
- Flash calculations may not converge for all conditions
- Property methods return mass-based values (per kg), not molar

---

### CoolProp

**Installation:**
```bash
pip install CoolProp
```

**Key Functions and Classes:**
- `CoolProp.PropsSI()` - Property evaluation
- `CoolProp.PhaseSI()` - Phase determination
- `CoolProp.HAPropsSI()` - Humid air properties
- Fluid names: 'Water', 'Air', 'R134a', etc.
- Input/output parameters: 'T', 'P', 'H', 'S', 'D', 'Q'

**Engineering Applications:**
- HVAC and refrigeration design
- Steam power cycles
- Cryogenic systems
- Heat pump analysis
- Humid air psychrometrics

**Code Examples:**

```python
from CoolProp.CoolProp import PropsSI, PhaseSI, HAPropsSI
import numpy as np

# Example 1: Rankine cycle analysis
P_high = 10e6  # High pressure (Pa) - 100 bar
P_low = 10e3   # Low pressure (Pa) - 0.1 bar
eta_pump = 0.85
eta_turbine = 0.88

# State 1: Saturated liquid at low pressure
h1 = PropsSI('H', 'P', P_low, 'Q', 0, 'Water')
s1 = PropsSI('S', 'P', P_low, 'Q', 0, 'Water')
T1 = PropsSI('T', 'P', P_low, 'Q', 0, 'Water')

# State 2: Compressed liquid after pump
v1 = 1 / PropsSI('D', 'P', P_low, 'Q', 0, 'Water')
w_pump_ideal = v1 * (P_high - P_low)
h2 = h1 + w_pump_ideal / eta_pump

# State 3: Superheated vapor at boiler exit
T3 = 500 + 273.15  # 500°C
h3 = PropsSI('H', 'P', P_high, 'T', T3, 'Water')
s3 = PropsSI('S', 'P', P_high, 'T', T3, 'Water')

# State 4: After turbine expansion
s4s = s3  # Isentropic
h4s = PropsSI('H', 'P', P_low, 'S', s4s, 'Water')
h4 = h3 - eta_turbine * (h3 - h4s)

# Cycle efficiency
q_in = h3 - h2
w_net = (h3 - h4) - (h2 - h1)
eta_cycle = w_net / q_in

print(f"Rankine Cycle Performance:")
print(f"  Thermal efficiency: {eta_cycle*100:.2f}%")
print(f"  Turbine work: {(h3-h4)/1000:.2f} kJ/kg")
print(f"  Heat input: {q_in/1000:.2f} kJ/kg")

# Example 2: Refrigeration cycle (vapor-compression)
T_evap = -10 + 273.15  # Evaporator temperature
T_cond = 40 + 273.15   # Condenser temperature
fluid = 'R134a'

P_evap = PropsSI('P', 'T', T_evap, 'Q', 0, fluid)
P_cond = PropsSI('P', 'T', T_cond, 'Q', 0, fluid)

# State 1: Saturated vapor leaving evaporator
h1 = PropsSI('H', 'P', P_evap, 'Q', 1, fluid)
s1 = PropsSI('S', 'P', P_evap, 'Q', 1, fluid)

# State 2: After isentropic compression
h2 = PropsSI('H', 'P', P_cond, 'S', s1, fluid)

# State 3: Saturated liquid leaving condenser
h3 = PropsSI('H', 'P', P_cond, 'Q', 0, fluid)

# State 4: After expansion valve (isenthalpic)
h4 = h3

# Performance metrics
q_evap = h1 - h4  # Cooling capacity
w_comp = h2 - h1  # Compressor work
COP = q_evap / w_comp

print(f"\nRefrigeration Cycle (R134a):")
print(f"  COP: {COP:.2f}")
print(f"  Cooling capacity: {q_evap/1000:.2f} kJ/kg")
print(f"  Compressor work: {w_comp/1000:.2f} kJ/kg")

# Example 3: Psychrometric calculations
T_db = 25 + 273.15  # Dry bulb temperature (K)
RH = 0.60           # Relative humidity (60%)
P_atm = 101325      # Atmospheric pressure (Pa)

# Calculate humid air properties
h = HAPropsSI('H', 'T', T_db, 'P', P_atm, 'R', RH)
w = HAPropsSI('W', 'T', T_db, 'P', P_atm, 'R', RH)
T_wb = HAPropsSI('T_wb', 'T', T_db, 'P', P_atm, 'R', RH)
T_dp = HAPropsSI('T_dp', 'T', T_db, 'P', P_atm, 'R', RH)

print(f"\nHumid Air Properties:")
print(f"  Enthalpy: {h/1000:.2f} kJ/kg_da")
print(f"  Humidity ratio: {w*1000:.2f} g_w/kg_da")
print(f"  Wet bulb temp: {T_wb-273.15:.2f}°C")
print(f"  Dew point: {T_dp-273.15:.2f}°C")
```

**Documentation:**
- http://www.coolprop.org/

**Common Gotchas:**
- All properties in SI units (Pa, K, J/kg, etc.)
- Requires exactly 2 state properties as input
- Quality 'Q' ranges from 0 (sat. liquid) to 1 (sat. vapor)
- Not all fluid/property combinations are valid
- Case-sensitive fluid and parameter names

---

## 3. Optimization

### scipy.optimize

**Installation:**
```bash
pip install scipy  # scipy.optimize is part of scipy
```

**Key Functions and Classes:**
- `scipy.optimize.minimize()` - Unconstrained and constrained optimization
- `scipy.optimize.root()` - Root finding
- `scipy.optimize.least_squares()` - Nonlinear least squares
- `scipy.optimize.curve_fit()` - Curve fitting
- `scipy.optimize.differential_evolution()` - Global optimization
- `scipy.optimize.linprog()` - Linear programming

**Engineering Applications:**
- Design optimization
- Parameter estimation from experimental data
- Process optimization
- Control system tuning
- Resource allocation
- Trajectory optimization

**Code Examples:**

```python
from scipy.optimize import minimize, curve_fit, differential_evolution, linprog
import numpy as np

# Example 1: Optimize beam cross-section for minimum weight
def beam_weight(x):
    """Minimize weight of rectangular beam
    x[0] = width, x[1] = height
    """
    width, height = x
    area = width * height
    length = 5.0  # m
    density = 7850  # kg/m³ (steel)
    return area * length * density

def beam_constraints(x):
    """Constraints: moment of inertia and stress limits"""
    width, height = x
    I = width * height**3 / 12  # Second moment of area
    M = 50000  # Applied moment (N·m)
    sigma_max = 250e6  # Allowable stress (Pa)

    # Stress constraint: M*y/I <= sigma_max
    stress_constraint = sigma_max - (M * height/2) / I

    # Minimum I constraint
    I_min = 1e-5  # m⁴
    I_constraint = I - I_min

    return [stress_constraint, I_constraint]

# Define constraints for optimizer
from scipy.optimize import NonlinearConstraint
nlc = NonlinearConstraint(beam_constraints, [0, 0], [np.inf, np.inf])

# Bounds on width and height
bounds = [(0.05, 0.5), (0.1, 1.0)]  # meters

# Initial guess
x0 = [0.2, 0.4]

result = minimize(beam_weight, x0, method='SLSQP', bounds=bounds, constraints=nlc)
print(f"Optimal beam dimensions:")
print(f"  Width: {result.x[0]*1000:.1f} mm")
print(f"  Height: {result.x[1]*1000:.1f} mm")
print(f"  Weight: {result.fun:.2f} kg")

# Example 2: Curve fitting for experimental data
# Fitting heat transfer correlation: Nu = C * Re^m * Pr^n
def nusselt_correlation(X, C, m, n):
    """Nusselt number correlation"""
    Re, Pr = X
    return C * Re**m * Pr**n

# Experimental data
Re_exp = np.array([5000, 10000, 20000, 50000, 100000])
Pr_exp = np.array([5.0, 5.2, 4.8, 5.1, 4.9])
Nu_exp = np.array([45, 82, 155, 350, 650])

# Fit the correlation
popt, pcov = curve_fit(nusselt_correlation, (Re_exp, Pr_exp), Nu_exp,
                       p0=[0.1, 0.8, 0.33])
C_fit, m_fit, n_fit = popt

print(f"\nHeat transfer correlation: Nu = {C_fit:.4f} * Re^{m_fit:.3f} * Pr^{n_fit:.3f}")
print(f"Standard errors: {np.sqrt(np.diag(pcov))}")

# Calculate R² value
Nu_pred = nusselt_correlation((Re_exp, Pr_exp), *popt)
ss_res = np.sum((Nu_exp - Nu_pred)**2)
ss_tot = np.sum((Nu_exp - np.mean(Nu_exp))**2)
r_squared = 1 - (ss_res / ss_tot)
print(f"R² = {r_squared:.4f}")

# Example 3: Production planning optimization (linear programming)
# Maximize profit: produce products A and B
# Profit: $50/unit for A, $60/unit for B
# Constraints: labor hours, material, machine time

c = [-50, -60]  # Negative for maximization

# Inequality constraints (A*x <= b)
# Labor: 2*A + 3*B <= 100 hours
# Material: 4*A + 2*B <= 120 units
# Machine: 1*A + 2*B <= 60 hours
A_ub = [[2, 3],
        [4, 2],
        [1, 2]]
b_ub = [100, 120, 60]

# Bounds: non-negative production
bounds = [(0, None), (0, None)]

result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')

print(f"\nOptimal production plan:")
print(f"  Product A: {result.x[0]:.1f} units")
print(f"  Product B: {result.x[1]:.1f} units")
print(f"  Maximum profit: ${-result.fun:.2f}")
```

**Documentation:**
- https://docs.scipy.org/doc/scipy/reference/optimize.html

**Common Gotchas:**
- `minimize()` performs minimization; negate objective for maximization
- Initial guess significantly affects convergence for nonlinear problems
- Constraint functions should return positive values for satisfied constraints
- Check `result.success` before using results
- Different methods work better for different problems (try 'SLSQP', 'trust-constr')
- Scale variables to similar magnitudes for better convergence

---

### Genetic Algorithms (DEAP)

**Installation:**
```bash
pip install deap
```

**Key Functions and Classes:**
- `deap.creator.create()` - Define fitness and individuals
- `deap.base.Toolbox()` - Register genetic operators
- `deap.algorithms.eaSimple()` - Simple evolutionary algorithm
- `deap.tools` - Selection, crossover, mutation operators
- Multi-objective optimization support (NSGA-II)

**Engineering Applications:**
- Complex optimization with multiple local optima
- Multi-objective design optimization
- Topology optimization
- Parameter tuning for complex systems
- Schedule optimization

**Code Examples:**

```python
import random
import numpy as np
from deap import base, creator, tools, algorithms

# Example 1: Truss structure optimization
# Minimize weight while satisfying stress constraints

# Define the problem
creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessMin)

def eval_truss(individual):
    """Evaluate truss design
    individual = list of member cross-sectional areas (m²)
    """
    areas = np.array(individual)

    # Truss geometry and loading (simplified 10-bar truss)
    length = 9.14  # member length (m)
    E = 68.9e9  # Young's modulus (Pa) - aluminum
    rho = 2770  # density (kg/m³)
    stress_max = 172e6  # allowable stress (Pa)

    # Calculate weight
    weight = np.sum(areas * length * rho * 9.81)

    # Simple stress calculation (simplified)
    # In real application, perform FEA
    force = 444800  # Applied force (N)
    max_stress = force / np.min(areas)

    # Penalty for constraint violation
    penalty = 0
    if max_stress > stress_max:
        penalty = 1e6 * (max_stress - stress_max)

    return (weight + penalty,)

# Setup genetic algorithm
toolbox = base.Toolbox()
toolbox.register("attr_area", random.uniform, 0.0001, 0.01)  # m²
toolbox.register("individual", tools.initRepeat, creator.Individual,
                 toolbox.attr_area, n=10)  # 10 truss members
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", eval_truss)
toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.001, indpb=0.2)
toolbox.register("select", tools.selTournament, tournsize=3)

# Run optimization
population = toolbox.population(n=50)
NGEN = 40
CXPB = 0.7  # Crossover probability
MUTPB = 0.2  # Mutation probability

# Evolution
for gen in range(NGEN):
    offspring = algorithms.varAnd(population, toolbox, CXPB, MUTPB)
    fits = toolbox.map(toolbox.evaluate, offspring)
    for fit, ind in zip(fits, offspring):
        ind.fitness.values = fit
    population = toolbox.select(offspring, k=len(population))

best_ind = tools.selBest(population, k=1)[0]
print(f"Best truss design weight: {best_ind.fitness.values[0]:.2f} N")
print(f"Optimal areas (mm²): {[a*1e6 for a in best_ind]}")

# Example 2: Multi-objective optimization (Pareto front)
# Optimize heat exchanger: minimize cost, maximize heat transfer

creator.create("FitnessMulti", base.Fitness, weights=(-1.0, 1.0))
creator.create("IndividualMulti", list, fitness=creator.FitnessMulti)

def eval_heat_exchanger(individual):
    """Evaluate heat exchanger design
    individual = [tube_diameter, tube_length, num_tubes]
    """
    D, L, N = individual

    # Cost model
    material_cost = 50  # $/kg
    density = 8000  # kg/m³
    thickness = 0.002  # m
    volume = np.pi * D * thickness * L * N
    cost = volume * density * material_cost

    # Heat transfer model (simplified)
    h = 5000  # heat transfer coefficient (W/m²K)
    A = np.pi * D * L * N  # heat transfer area
    Q = h * A * 50  # heat transfer rate (W), 50K temp difference

    return (cost, Q)

toolbox_mo = base.Toolbox()
toolbox_mo.register("attr_d", random.uniform, 0.01, 0.05)  # diameter (m)
toolbox_mo.register("attr_l", random.uniform, 0.5, 3.0)    # length (m)
toolbox_mo.register("attr_n", random.randint, 10, 200)     # number of tubes

toolbox_mo.register("individual", tools.initCycle, creator.IndividualMulti,
                    (toolbox_mo.attr_d, toolbox_mo.attr_l, toolbox_mo.attr_n), n=1)
toolbox_mo.register("population", tools.initRepeat, list, toolbox_mo.individual)

toolbox_mo.register("evaluate", eval_heat_exchanger)
toolbox_mo.register("mate", tools.cxBlend, alpha=0.5)
toolbox_mo.register("mutate", tools.mutGaussian, mu=0, sigma=0.2, indpb=0.3)
toolbox_mo.register("select", tools.selNSGA2)  # Multi-objective selection

pop_mo = toolbox_mo.population(n=100)
algorithms.eaSimple(pop_mo, toolbox_mo, cxpb=0.7, mutpb=0.3, ngen=30, verbose=False)

# Extract Pareto front
pareto = tools.sortNondominated(pop_mo, len(pop_mo), first_front_only=True)[0]
print(f"\nPareto front: {len(pareto)} optimal solutions")
print("Sample trade-offs (Cost, Heat Transfer):")
for ind in pareto[:3]:
    print(f"  ${ind.fitness.values[0]:.0f}, {ind.fitness.values[1]/1000:.1f} kW")

# Example 3: Process parameter optimization with constraints
def eval_process(individual):
    """Optimize manufacturing process parameters
    Minimize: production time + defect rate
    Variables: [temperature, pressure, flow_rate]
    """
    temp, pressure, flow = individual

    # Process model (empirical)
    time = 100 / (temp - 150) + 50 / flow
    defect_rate = 0.01 * abs(pressure - 5) + 0.02 * abs(temp - 200)

    # Constraints
    penalty = 0
    if temp < 160 or temp > 240:
        penalty += 1000
    if pressure < 3 or pressure > 7:
        penalty += 1000
    if flow < 1 or flow > 10:
        penalty += 1000

    return (time + defect_rate * 100 + penalty,)

# This uses the same pattern as Example 1
print("\nGenetic algorithms are powerful for complex optimization problems")
print("where gradient information is unavailable or unreliable.")
```

**Documentation:**
- https://deap.readthedocs.io/

**Common Gotchas:**
- Must create fitness and individual classes before use
- Fitness weights: positive for maximization, negative for minimization
- Population size and generations significantly affect solution quality
- Constraints require penalty functions or repair mechanisms
- Random seed affects reproducibility
- Multi-objective optimization returns Pareto front, not single solution

---

## 4. Visualization

### Matplotlib

**Installation:**
```bash
pip install matplotlib
```

**Key Functions and Classes:**
- `matplotlib.pyplot` - MATLAB-style interface
- `plt.plot()`, `plt.scatter()`, `plt.bar()` - Basic plots
- `plt.subplot()`, `plt.subplots()` - Multiple plots
- `plt.contour()`, `plt.contourf()` - Contour plots
- `plt.imshow()` - Image display
- `plt.savefig()` - Export figures

**Engineering Applications:**
- Data visualization and analysis
- Results presentation
- Control charts and quality control
- Frequency response plots (Bode, Nyquist)
- Stress-strain curves
- Time series analysis

**Code Examples:**

```python
import matplotlib.pyplot as plt
import numpy as np

# Example 1: Stress-strain curve with engineering annotations
strain = np.array([0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2])
stress = np.array([0, 200, 400, 800, 1200, 1400, 1500, 1450, 1400, 1350])  # MPa

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(strain * 100, stress, 'b-', linewidth=2, label='Test data')

# Mark yield point
yield_idx = 4
ax.plot(strain[yield_idx] * 100, stress[yield_idx], 'ro', markersize=10,
        label=f'Yield point ({stress[yield_idx]} MPa)')

# Calculate and plot elastic modulus
elastic_strain = strain[:3]
elastic_stress = stress[:3]
E = np.polyfit(elastic_strain, elastic_stress, 1)[0]
ax.plot(elastic_strain * 100, elastic_stress, 'r--', linewidth=2,
        label=f'Elastic region (E = {E/1000:.0f} GPa)')

ax.set_xlabel('Strain (%)', fontsize=12)
ax.set_ylabel('Stress (MPa)', fontsize=12)
ax.set_title('Tensile Test Results - Steel Specimen', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig('/tmp/stress_strain.png', dpi=300)
print("Stress-strain curve saved to /tmp/stress_strain.png")

# Example 2: Temperature distribution contour plot
x = np.linspace(0, 1, 50)
y = np.linspace(0, 1, 50)
X, Y = np.meshgrid(x, y)

# Simulate 2D heat distribution
T = 100 * np.exp(-((X-0.5)**2 + (Y-0.5)**2) / 0.1) + 20

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Filled contour plot
levels = np.linspace(20, 120, 11)
cf = ax1.contourf(X, Y, T, levels=levels, cmap='hot')
ax1.contour(X, Y, T, levels=levels, colors='black', linewidths=0.5, alpha=0.3)
cbar1 = plt.colorbar(cf, ax=ax1)
cbar1.set_label('Temperature (°C)', fontsize=11)
ax1.set_xlabel('x (m)', fontsize=11)
ax1.set_ylabel('y (m)', fontsize=11)
ax1.set_title('Temperature Distribution - Contour Plot', fontsize=12, fontweight='bold')
ax1.set_aspect('equal')

# 3D surface plot
ax2.remove()
ax2 = fig.add_subplot(122, projection='3d')
surf = ax2.plot_surface(X, Y, T, cmap='hot', edgecolor='none', alpha=0.8)
cbar2 = plt.colorbar(surf, ax=ax2, shrink=0.5)
cbar2.set_label('Temperature (°C)', fontsize=10)
ax2.set_xlabel('x (m)', fontsize=10)
ax2.set_ylabel('y (m)', fontsize=10)
ax2.set_zlabel('Temperature (°C)', fontsize=10)
ax2.set_title('Temperature Distribution - 3D View', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('/tmp/temperature_field.png', dpi=300)
print("Temperature distribution saved to /tmp/temperature_field.png")

# Example 3: Bode plot for control system analysis
frequency = np.logspace(-1, 3, 500)  # Hz
omega = 2 * np.pi * frequency  # rad/s

# Transfer function: H(s) = K / (s^2 + 2*zeta*wn*s + wn^2)
K = 100
wn = 2 * np.pi * 10  # Natural frequency (rad/s)
zeta = 0.3  # Damping ratio

s = 1j * omega
H = K / (s**2 + 2*zeta*wn*s + wn**2)

magnitude_dB = 20 * np.log10(np.abs(H))
phase_deg = np.angle(H, deg=True)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# Magnitude plot
ax1.semilogx(frequency, magnitude_dB, 'b-', linewidth=2)
ax1.axhline(y=-3, color='r', linestyle='--', alpha=0.5, label='-3 dB line')
ax1.grid(True, which='both', alpha=0.3)
ax1.set_ylabel('Magnitude (dB)', fontsize=11)
ax1.set_title('Bode Plot - Second Order System', fontsize=12, fontweight='bold')
ax1.legend()

# Phase plot
ax2.semilogx(frequency, phase_deg, 'r-', linewidth=2)
ax2.axhline(y=-180, color='b', linestyle='--', alpha=0.5, label='-180° line')
ax2.grid(True, which='both', alpha=0.3)
ax2.set_xlabel('Frequency (Hz)', fontsize=11)
ax2.set_ylabel('Phase (degrees)', fontsize=11)
ax2.legend()

plt.tight_layout()
plt.savefig('/tmp/bode_plot.png', dpi=300)
print("Bode plot saved to /tmp/bode_plot.png")

plt.close('all')  # Close all figures
```

**Documentation:**
- https://matplotlib.org/stable/contents.html

**Common Gotchas:**
- Default figure sizes may be too small; specify `figsize` parameter
- Use `plt.tight_layout()` to prevent label clipping
- Save figures before `plt.show()` in scripts
- Use vectorized operations instead of loops for better performance
- Color cycle repeats after 10 colors
- `plt` (pyplot) maintains state; use object-oriented interface for complex plots

---

### PyVista

**Installation:**
```bash
pip install pyvista
```

**Key Functions and Classes:**
- `pyvista.Plotter()` - 3D visualization environment
- `pv.read()` - Read mesh files (VTK, STL, etc.)
- `pv.PolyData()`, `pv.UnstructuredGrid()` - Mesh creation
- `.plot()`, `.show()` - Display visualizations
- `.warp_by_scalar()`, `.contour()` - Data transformations
- `.slice()`, `.clip()` - Mesh operations

**Engineering Applications:**
- FEA results visualization
- CAD model visualization
- Flow field visualization
- 3D mesh manipulation
- Structural deformation plots
- Slice and contour analysis

**Code Examples:**

```python
import pyvista as pv
import numpy as np

# Example 1: Visualize stress distribution on a beam
# Create beam mesh
length, width, height = 2.0, 0.1, 0.2
beam = pv.Box(bounds=[0, length, 0, width, 0, height])

# Generate synthetic stress data
points = beam.points
x_coords = points[:, 0]

# Stress varies parabolically along length (bending stress)
max_stress = 200e6  # Pa
stress = max_stress * (x_coords / length) * (1 - x_coords / length) * 4
beam['von_mises_stress'] = stress / 1e6  # Convert to MPa

# Create plotter and add beam
plotter = pv.Plotter()
plotter.add_mesh(beam, scalars='von_mises_stress', cmap='jet',
                 show_edges=True, scalar_bar_args={'title': 'Stress (MPa)'})
plotter.add_text('Cantilever Beam Stress Distribution', font_size=12)
plotter.camera_position = 'iso'
plotter.show(screenshot='/tmp/beam_stress.png')
print("Beam stress visualization saved to /tmp/beam_stress.png")

# Example 2: Slice through a 3D temperature field
# Create a 3D grid
grid = pv.ImageData(dimensions=(50, 50, 50))
grid.spacing = (0.1, 0.1, 0.1)  # 5m x 5m x 5m volume

# Generate temperature distribution (heat source at center)
x, y, z = grid.points.T
x0, y0, z0 = 2.5, 2.5, 2.5  # Heat source location
r = np.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
T_ambient = 20  # °C
T_source = 200  # °C
temperature = T_ambient + (T_source - T_ambient) * np.exp(-r**2 / 2)
grid['Temperature'] = temperature

# Create visualization with multiple slices
plotter = pv.Plotter()
plotter.add_mesh(grid.outline(), color='black')

# Add three orthogonal slices
slices = grid.slice_orthogonal()
plotter.add_mesh(slices, scalars='Temperature', cmap='hot',
                 scalar_bar_args={'title': 'Temperature (°C)'})

plotter.add_text('3D Temperature Field with Slices', font_size=12)
plotter.show(screenshot='/tmp/temperature_slices.png')
print("Temperature field slices saved to /tmp/temperature_slices.png")

# Example 3: Deformed shape visualization
# Create a plate mesh
plate = pv.Plane(i_resolution=20, j_resolution=20,
                 i_size=1.0, j_size=1.0)

# Generate displacement field (simulated bending)
points = plate.points
x, y = points[:, 0], points[:, 1]
displacement_z = 0.1 * x * (1 - x) * y * (1 - y) * 16  # Maximum 0.1m deflection
displacement = np.zeros_like(points)
displacement[:, 2] = displacement_z

plate['Displacement'] = np.linalg.norm(displacement, axis=1)
plate['Displacement_vec'] = displacement

# Create plotter showing both original and deformed
plotter = pv.Plotter(shape=(1, 2))

# Original shape
plotter.subplot(0, 0)
plotter.add_mesh(plate, color='lightblue', show_edges=True)
plotter.add_text('Original Shape', font_size=10)
plotter.camera_position = 'iso'

# Deformed shape with color map
plotter.subplot(0, 1)
warped = plate.warp_by_vector('Displacement_vec', factor=1.0)
plotter.add_mesh(warped, scalars='Displacement', cmap='rainbow',
                 show_edges=True, scalar_bar_args={'title': 'Displacement (m)'})
plotter.add_text('Deformed Shape', font_size=10)
plotter.camera_position = 'iso'

plotter.show(screenshot='/tmp/deformed_shape.png')
print("Deformed shape comparison saved to /tmp/deformed_shape.png")
```

**Documentation:**
- https://docs.pyvista.org/

**Common Gotchas:**
- Requires VTK backend; may have installation issues on some systems
- Screenshots require `screenshot` parameter in `.show()`
- Scalar data must match number of points or cells
- Some operations modify meshes in-place
- Use `.copy()` to avoid unintended modifications
- Interactive plots don't work well in notebooks without proper backend

---

### Plotly

**Installation:**
```bash
pip install plotly
```

**Key Functions and Classes:**
- `plotly.graph_objects` - Low-level plotting interface
- `plotly.express` - High-level plotting interface
- `go.Figure()` - Create figure objects
- `px.scatter()`, `px.line()`, etc. - Quick plots
- `go.Surface()`, `go.Mesh3d()` - 3D plotting
- Interactive features: zoom, pan, hover

**Engineering Applications:**
- Interactive dashboards
- Web-based data exploration
- Dynamic reporting
- Multi-dimensional data visualization
- Real-time monitoring interfaces

**Code Examples:**

```python
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import pandas as pd

# Example 1: Interactive failure envelope (Mohr-Coulomb)
theta = np.linspace(0, 2*np.pi, 100)
sigma_1_vals = np.linspace(0, 100, 50)  # MPa
sigma_3_vals = np.linspace(0, 50, 50)   # MPa

# Material properties
c = 10  # Cohesion (MPa)
phi = 30 * np.pi / 180  # Friction angle (rad)

# Mohr-Coulomb failure criterion
S1, S3 = np.meshgrid(sigma_1_vals, sigma_3_vals)
failure_function = (S1 - S3) - 2*c*np.cos(phi) - (S1 + S3)*np.sin(phi)

fig = go.Figure(data=[
    go.Contour(
        z=failure_function,
        x=sigma_1_vals,
        y=sigma_3_vals,
        colorscale='RdYlGn',
        contours=dict(
            start=-50,
            end=50,
            size=10,
            showlabels=True
        ),
        colorbar=dict(title="Failure Index"),
        hovertemplate='σ1: %{x:.1f} MPa<br>σ3: %{y:.1f} MPa<br>Index: %{z:.1f}<extra></extra>'
    )
])

fig.add_contour(
    z=failure_function,
    x=sigma_1_vals,
    y=sigma_3_vals,
    contours=dict(
        start=0,
        end=0,
        size=1,
        coloring='lines'
    ),
    line=dict(color='red', width=3),
    showscale=False,
    name='Failure Envelope'
)

fig.update_layout(
    title='Mohr-Coulomb Failure Envelope (Interactive)',
    xaxis_title='Major Principal Stress σ1 (MPa)',
    yaxis_title='Minor Principal Stress σ3 (MPa)',
    width=800,
    height=600,
    hovermode='closest'
)

fig.write_html('/tmp/failure_envelope.html')
print("Interactive failure envelope saved to /tmp/failure_envelope.html")

# Example 2: 3D surface plot - pressure distribution on a wing
# Generate wing geometry and pressure data
x = np.linspace(0, 10, 100)  # Chord direction (m)
y = np.linspace(-5, 5, 80)   # Span direction (m)
X, Y = np.meshgrid(x, y)

# Simulate pressure distribution (simplified)
# Negative pressure (suction) on top
chord = 10
pressure = -5000 * np.exp(-((X-chord*0.3)**2) / (chord*0.1)**2) * \
           np.exp(-(Y**2) / 4) + 101325  # Pa

# Convert to coefficient of pressure
q_inf = 0.5 * 1.225 * 50**2  # Dynamic pressure at 50 m/s
Cp = (pressure - 101325) / q_inf

fig = go.Figure(data=[
    go.Surface(
        x=X,
        y=Y,
        z=Cp,
        colorscale='RdBu',
        colorbar=dict(title="Cp"),
        hovertemplate='x: %{x:.2f}m<br>y: %{y:.2f}m<br>Cp: %{z:.3f}<extra></extra>'
    )
])

fig.update_layout(
    title='Pressure Coefficient Distribution on Wing',
    scene=dict(
        xaxis_title='Chord (m)',
        yaxis_title='Span (m)',
        zaxis_title='Cp',
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.2))
    ),
    width=900,
    height=700
)

fig.write_html('/tmp/wing_pressure.html')
print("Wing pressure distribution saved to /tmp/wing_pressure.html")

# Example 3: Time-series analysis with multiple y-axes
# Simulate pump performance data
time = np.linspace(0, 24, 1000)  # 24 hours
base_flow = 100  # m³/h
flow = base_flow + 10*np.sin(2*np.pi*time/12) + 2*np.random.randn(len(time))

base_power = 50  # kW
power = base_power * (flow/base_flow)**3 + 1*np.random.randn(len(time))

efficiency = (flow * 10 * 1000 * 9.81) / (power * 1000 * 3600) * 100  # %

# Create figure with secondary y-axis
fig = go.Figure()

fig.add_trace(go.Scatter(
    x=time, y=flow,
    name='Flow Rate',
    line=dict(color='blue', width=2),
    hovertemplate='Time: %{x:.1f}h<br>Flow: %{y:.1f} m³/h<extra></extra>'
))

fig.add_trace(go.Scatter(
    x=time, y=power,
    name='Power',
    line=dict(color='red', width=2),
    yaxis='y2',
    hovertemplate='Time: %{x:.1f}h<br>Power: %{y:.1f} kW<extra></extra>'
))

fig.add_trace(go.Scatter(
    x=time, y=efficiency,
    name='Efficiency',
    line=dict(color='green', width=2, dash='dash'),
    yaxis='y3',
    hovertemplate='Time: %{x:.1f}h<br>Efficiency: %{y:.1f}%<extra></extra>'
))

# Create axis objects
fig.update_layout(
    title='Pump Performance Monitoring (24 hours)',
    xaxis=dict(title='Time (hours)'),
    yaxis=dict(
        title='Flow Rate (m³/h)',
        titlefont=dict(color='blue'),
        tickfont=dict(color='blue')
    ),
    yaxis2=dict(
        title='Power (kW)',
        titlefont=dict(color='red'),
        tickfont=dict(color='red'),
        anchor='free',
        overlaying='y',
        side='right',
        position=0.95
    ),
    yaxis3=dict(
        title='Efficiency (%)',
        titlefont=dict(color='green'),
        tickfont=dict(color='green'),
        anchor='free',
        overlaying='y',
        side='right',
        position=1.0
    ),
    width=1000,
    height=600,
    hovermode='x unified'
)

fig.write_html('/tmp/pump_monitoring.html')
print("Pump monitoring dashboard saved to /tmp/pump_monitoring.html")
```

**Documentation:**
- https://plotly.com/python/

**Common Gotchas:**
- Large datasets can make interactive plots slow
- Default color scales may not be appropriate for engineering data
- Multiple y-axes require careful layout configuration
- HTML files can be large with complex plots
- Some features require `plotly-orca` for static image export
- Use `hovertemplate` for custom hover information

---

## 5. Units

### Pint

**Installation:**
```bash
pip install pint
```

**Key Functions and Classes:**
- `pint.UnitRegistry()` - Create unit registry
- `Q_()` - Create quantities with units
- `.to()` - Unit conversion
- `.magnitude`, `.units` - Access value and units
- Context managers for different unit systems
- Automatic unit checking in calculations

**Engineering Applications:**
- Unit conversion and validation
- Preventing unit-related errors
- Multi-unit system calculations
- Dimensional analysis
- Documentation through code

**Code Examples:**

```python
import pint
import numpy as np

# Create unit registry
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# Example 1: Pipe flow calculations with automatic unit handling
# Calculate pressure drop in a pipeline

diameter = Q_(4, 'inch').to('meter')
length = Q_(100, 'feet').to('meter')
flow_rate = Q_(500, 'gallon/minute').to('meter**3/second')
roughness = Q_(0.0018, 'inch').to('meter')
viscosity = Q_(1, 'centipoise').to('pascal*second')
density = Q_(62.4, 'pound/foot**3').to('kilogram/meter**3')

# Calculate velocity
area = np.pi * (diameter/2)**2
velocity = flow_rate / area

# Reynolds number
Re = (density * velocity * diameter / viscosity).to_base_units()
print(f"Reynolds number: {Re:.0f}")

# Friction factor (Colebrook equation - simplified)
eD = roughness / diameter
fd = 0.02  # Initial guess, would iterate normally

# Pressure drop
dP = (fd * length * density * velocity**2) / (2 * diameter)
dP_psi = dP.to('psi')
dP_bar = dP.to('bar')

print(f"Pressure drop: {dP_psi:.2f} = {dP_bar:.3f}")
print(f"Velocity: {velocity.to('feet/second'):.2f}")

# Example 2: Heat exchanger calculations
# Calculate heat transfer and validate energy balance

m_hot = Q_(5000, 'kg/hour')  # Hot fluid mass flow rate
T_hot_in = Q_(90, 'degC')
T_hot_out = Q_(60, 'degC')
Cp_hot = Q_(4.18, 'kJ/(kg*K)')

m_cold = Q_(150, 'lb/minute')  # Cold fluid mass flow rate
T_cold_in = Q_(68, 'degF')
Cp_cold = Q_(1.0, 'BTU/(lb*degF)')

# Heat given up by hot fluid
Q_hot = m_hot * Cp_hot * (T_hot_in - T_hot_out)
Q_hot_kW = Q_hot.to('kilowatt')

print(f"\nHeat transfer rate: {Q_hot_kW:.2f}")

# Calculate cold fluid outlet temperature
# Q_hot = Q_cold, so T_cold_out = T_cold_in + Q_hot/(m_cold*Cp_cold)
T_cold_out = T_cold_in + Q_hot / (m_cold * Cp_cold)
T_cold_out_C = T_cold_out.to('degC')

print(f"Cold fluid outlet temperature: {T_cold_out_C:.1f}")

# Verify energy balance
Q_cold = m_cold * Cp_cold * (T_cold_out - T_cold_in)
Q_cold_kW = Q_cold.to('kilowatt')
error = abs((Q_hot_kW - Q_cold_kW) / Q_hot_kW * 100)

print(f"Energy balance error: {error:.2f}%")

# Example 3: Mechanical power and efficiency
# Motor driving a pump

# Pump parameters
flow = Q_(100, 'm**3/hour')
head = Q_(50, 'meter')
rho = Q_(1000, 'kg/m**3')
g = Q_(9.81, 'm/s**2')

# Hydraulic power
P_hydraulic = rho * g * flow * head
P_hydraulic_kW = P_hydraulic.to('kilowatt')

print(f"\nPump hydraulic power: {P_hydraulic_kW:.2f}")

# Motor power (electrical)
voltage = Q_(460, 'volt')
current = Q_(25, 'ampere')
power_factor = 0.85
efficiency = 0.92

P_electrical = np.sqrt(3) * voltage * current * power_factor
P_shaft = P_electrical * efficiency

print(f"Motor electrical power: {P_electrical.to('kilowatt'):.2f}")
print(f"Motor shaft power: {P_shaft.to('kilowatt'):.2f}")
print(f"Pump efficiency: {(P_hydraulic/P_shaft*100):.1f}%")

# Convert to horsepower
P_hp = P_shaft.to('horsepower')
print(f"Motor power: {P_hp:.1f}")

# Example 4: Context for different unit systems
# Work with multiple unit systems seamlessly

with ureg.context('US'):
    # US customary units
    force = Q_(100, 'lbf')
    area = Q_(2, 'inch**2')
    stress = (force / area).to('psi')
    print(f"\nStress (US units): {stress:.0f}")

with ureg.context('SI'):
    # SI units
    stress_SI = stress.to('MPa')
    print(f"Stress (SI units): {stress_SI:.2f}")
```

**Documentation:**
- https://pint.readthedocs.io/

**Common Gotchas:**
- Must use same unit registry throughout code
- Some units have multiple names (e.g., 'liter' vs 'litre')
- Temperature conversions: use 'degC' for Celsius, 'degF' for Fahrenheit
- Offset units (temperatures) require special handling
- Array operations need NumPy arrays: `Q_(np.array([1,2,3]), 'meter')`
- Unit strings are case-sensitive

---

### astropy.units

**Installation:**
```bash
pip install astropy
```

**Key Functions and Classes:**
- `astropy.units` (import as `u`) - Unit definitions
- `u.Quantity()` - Create quantities
- `.to()` - Unit conversion
- `.si`, `.cgs` - Convert to unit systems
- Equivalencies for special conversions
- Physical constants from `astropy.constants`

**Engineering Applications:**
- Scientific and aerospace calculations
- Physical constants with units
- Spectroscopy and radiation
- Advanced dimensional analysis
- Electromagnetic calculations

**Code Examples:**

```python
from astropy import units as u
from astropy import constants as const
import numpy as np

# Example 1: Orbital mechanics calculations
# Calculate orbital period and velocity for a satellite

altitude = 400 * u.km
R_earth = 6371 * u.km
r_orbit = R_earth + altitude

# Using physical constants
M_earth = const.M_earth
G = const.G

# Orbital velocity
v_orbit = np.sqrt(G * M_earth / r_orbit).to(u.km/u.s)
print(f"Orbital velocity: {v_orbit:.2f}")

# Orbital period
T_orbit = (2 * np.pi * np.sqrt(r_orbit**3 / (G * M_earth))).to(u.minute)
print(f"Orbital period: {T_orbit:.1f}")

# Orbital energy per unit mass
E_specific = -G * M_earth / (2 * r_orbit)
print(f"Specific orbital energy: {E_specific.to(u.MJ/u.kg):.2f}")

# Example 2: Electromagnetic calculations
# Antenna power and signal strength

# Transmitter parameters
P_transmit = 100 * u.W
frequency = 2.4 * u.GHz
wavelength = (const.c / frequency).to(u.cm)

print(f"\nWavelength: {wavelength:.2f}")

# Antenna gain
gain_dBi = 20 * u.dB  # Decibels are dimensionless
gain_linear = 10**(gain_dBi.value / 10) * u.dimensionless_unscaled

# Received power at distance
distance = 1 * u.km
# Friis transmission equation
P_received = (P_transmit * gain_linear**2 * wavelength**2 /
              (16 * np.pi**2 * distance**2)).to(u.mW)

print(f"Received power at {distance}: {P_received:.3e}")
print(f"Received power: {10*np.log10(P_received.value):.1f} dBm")

# Example 3: Radiation heat transfer
# Stefan-Boltzmann law for thermal radiation

# Object parameters
emissivity = 0.9 * u.dimensionless_unscaled
area = 2 * u.m**2
T_object = (500 + 273.15) * u.K
T_surroundings = (20 + 273.15) * u.K

# Stefan-Boltzmann constant
sigma = const.sigma_sb

# Radiated power
P_radiated = emissivity * sigma * area * T_object**4
P_absorbed = emissivity * sigma * area * T_surroundings**4
P_net = (P_radiated - P_absorbed).to(u.kW)

print(f"\nNet radiation heat transfer: {P_net:.2f}")

# Time to cool (simplified, assuming constant properties)
mass = 50 * u.kg
Cp = 500 * u.J / (u.kg * u.K)
dT = 100 * u.K  # Temperature drop

energy_loss = mass * Cp * dT
time_to_cool = (energy_loss / P_net).to(u.hour)
print(f"Approximate cooling time: {time_to_cool:.1f}")

# Example 4: Spectroscopy and photon energy
# Calculate photon energy and number of photons

# Laser parameters
wavelength_laser = 532 * u.nm  # Green laser
power = 5 * u.mW

# Photon energy
E_photon = (const.h * const.c / wavelength_laser).to(u.eV)
print(f"\nPhoton energy at {wavelength_laser}: {E_photon:.3f}")

# Number of photons per second
E_photon_J = E_photon.to(u.J)
n_photons = (power / E_photon_J).to(1/u.s)
print(f"Photon flux: {n_photons:.2e}")

# Energy conversion using equivalencies
# Convert between wavelength and frequency
freq = wavelength_laser.to(u.THz, equivalencies=u.spectral())
print(f"Frequency: {freq:.1f}")

# Example 5: Pressure conversions with physical contexts
# Various pressure units in engineering

pressure_psi = 100 * u.psi
pressure_bar = pressure_psi.to(u.bar)
pressure_Pa = pressure_psi.to(u.Pa)
pressure_atm = pressure_psi.to(u.atm)
pressure_torr = pressure_psi.to(u.Torr)

print(f"\nPressure conversions:")
print(f"  {pressure_psi:.1f}")
print(f"  {pressure_bar:.2f}")
print(f"  {pressure_Pa:.0f}")
print(f"  {pressure_atm:.2f}")
print(f"  {pressure_torr:.1f}")
```

**Documentation:**
- https://docs.astropy.org/en/stable/units/

**Common Gotchas:**
- Some units require equivalencies for conversion (wavelength ↔ frequency)
- Physical constants include uncertainties
- Use `.value` to get dimensionless number
- Can't mix astropy and pint units directly
- Temperature: use Kelvin by default, conversions can be tricky
- Logarithmic units (dB, magnitudes) have special handling

---

## 6. Network Analysis

### NetworkX

**Installation:**
```bash
pip install networkx
```

**Key Functions and Classes:**
- `nx.Graph()`, `nx.DiGraph()` - Create graphs
- `nx.add_node()`, `nx.add_edge()` - Build networks
- `nx.shortest_path()` - Path finding
- `nx.connected_components()` - Component analysis
- `nx.pagerank()`, `nx.betweenness_centrality()` - Centrality measures
- `nx.draw()` - Basic visualization

**Engineering Applications:**
- Pipeline and utility network analysis
- Traffic flow optimization
- Supply chain networks
- Electrical grid analysis
- Structural connectivity
- Project scheduling (PERT/CPM)

**Code Examples:**

```python
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

# Example 1: Pipeline network analysis
# Find optimal routing and critical paths

# Create directed graph for water distribution
G = nx.DiGraph()

# Add nodes (junctions) with elevation data
nodes = {
    'Reservoir': {'elevation': 100, 'demand': 0},
    'J1': {'elevation': 95, 'demand': 50},
    'J2': {'elevation': 90, 'demand': 75},
    'J3': {'elevation': 85, 'demand': 60},
    'J4': {'elevation': 80, 'demand': 40},
    'Tank': {'elevation': 75, 'demand': 0}
}

for node, attrs in nodes.items():
    G.add_node(node, **attrs)

# Add edges (pipes) with length and diameter
pipes = [
    ('Reservoir', 'J1', {'length': 500, 'diameter': 300, 'roughness': 0.15}),
    ('J1', 'J2', {'length': 400, 'diameter': 250, 'roughness': 0.15}),
    ('J1', 'J3', {'length': 600, 'diameter': 200, 'roughness': 0.15}),
    ('J2', 'J4', {'length': 350, 'diameter': 200, 'roughness': 0.15}),
    ('J3', 'J4', {'length': 450, 'diameter': 200, 'roughness': 0.15}),
    ('J4', 'Tank', {'length': 300, 'diameter': 250, 'roughness': 0.15}),
]

for source, target, attrs in pipes:
    G.add_edge(source, target, **attrs)

# Find all paths from reservoir to tank
all_paths = list(nx.all_simple_paths(G, 'Reservoir', 'Tank'))
print(f"Number of paths from Reservoir to Tank: {len(all_paths)}")
for i, path in enumerate(all_paths, 1):
    path_length = sum(G[path[j]][path[j+1]]['length'] for j in range(len(path)-1))
    print(f"  Path {i}: {' -> '.join(path)} (Length: {path_length}m)")

# Calculate hydraulic gradient
shortest_path = nx.shortest_path(G, 'Reservoir', 'Tank', weight='length')
total_length = sum(G[shortest_path[i]][shortest_path[i+1]]['length']
                   for i in range(len(shortest_path)-1))
elevation_drop = nodes['Reservoir']['elevation'] - nodes['Tank']['elevation']
hydraulic_gradient = elevation_drop / total_length

print(f"\nShortest path: {' -> '.join(shortest_path)}")
print(f"Hydraulic gradient: {hydraulic_gradient:.4f} m/m")

# Example 2: Electrical power grid reliability analysis
# Identify critical nodes and edges

# Create power grid network
power_grid = nx.Graph()

# Substations
substations = ['Main', 'North', 'South', 'East', 'West', 'Central']
power_grid.add_nodes_from(substations)

# Transmission lines with capacity (MW)
lines = [
    ('Main', 'North', 500),
    ('Main', 'Central', 800),
    ('Central', 'North', 300),
    ('Central', 'South', 400),
    ('Central', 'East', 350),
    ('Central', 'West', 400),
    ('North', 'East', 250),
    ('South', 'West', 300),
]

for source, target, capacity in lines:
    power_grid.add_edge(source, target, capacity=capacity)

# Calculate centrality measures
betweenness = nx.betweenness_centrality(power_grid, weight='capacity')
closeness = nx.closeness_centrality(power_grid, distance='capacity')

print(f"\nCritical substations (by betweenness centrality):")
for node, centrality in sorted(betweenness.items(), key=lambda x: x[1], reverse=True):
    print(f"  {node}: {centrality:.3f}")

# Find edge connectivity (minimum number of edges to disconnect graph)
edge_connectivity = nx.edge_connectivity(power_grid)
print(f"\nEdge connectivity: {edge_connectivity}")

# Identify critical edges (bridges that would disconnect the network)
edge_betweenness = nx.edge_betweenness_centrality(power_grid, weight='capacity')
print(f"\nCritical transmission lines:")
for edge, centrality in sorted(edge_betweenness.items(),
                               key=lambda x: x[1], reverse=True)[:3]:
    capacity = power_grid[edge[0]][edge[1]]['capacity']
    print(f"  {edge[0]} - {edge[1]}: {centrality:.3f} (Capacity: {capacity} MW)")

# Example 3: Project scheduling with CPM (Critical Path Method)
# Create project network and find critical path

project = nx.DiGraph()

# Activities: (name, duration, predecessors)
activities = [
    ('Start', 0, []),
    ('A_Design', 10, ['Start']),
    ('B_Procurement', 5, ['A_Design']),
    ('C_Site_Prep', 7, ['Start']),
    ('D_Foundation', 12, ['C_Site_Prep']),
    ('E_Structure', 15, ['D_Foundation', 'B_Procurement']),
    ('F_MEP', 10, ['E_Structure']),
    ('G_Finishes', 8, ['F_MEP']),
    ('H_Commissioning', 5, ['G_Finishes']),
    ('End', 0, ['H_Commissioning']),
]

# Build network
for activity, duration, predecessors in activities:
    project.add_node(activity, duration=duration)
    if not predecessors:
        predecessors = []
    for pred in predecessors:
        project.add_edge(pred, activity)

# Calculate earliest start (ES) and earliest finish (EF)
ES = {'Start': 0}
EF = {}

for node in nx.topological_sort(project):
    if node == 'Start':
        EF[node] = ES[node] + project.nodes[node]['duration']
    else:
        predecessors = list(project.predecessors(node))
        ES[node] = max(EF[pred] for pred in predecessors)
        EF[node] = ES[node] + project.nodes[node]['duration']

# Calculate latest start (LS) and latest finish (LF)
LF = {'End': EF['End']}
LS = {}

for node in reversed(list(nx.topological_sort(project))):
    if node == 'End':
        LS[node] = LF[node] - project.nodes[node]['duration']
    else:
        successors = list(project.successors(node))
        LF[node] = min(LS[succ] for succ in successors)
        LS[node] = LF[node] - project.nodes[node]['duration']

# Identify critical path (activities with zero slack)
critical_activities = []
for node in project.nodes():
    slack = LS[node] - ES[node]
    if slack == 0 and node not in ['Start', 'End']:
        critical_activities.append(node)

print(f"\nProject Schedule Analysis:")
print(f"Total project duration: {EF['End']} days")
print(f"\nCritical path activities: {' -> '.join(critical_activities)}")

print(f"\nActivity schedule:")
for node in nx.topological_sort(project):
    if node not in ['Start', 'End']:
        duration = project.nodes[node]['duration']
        slack = LS[node] - ES[node]
        critical = '*' if slack == 0 else ' '
        print(f"  {critical} {node:20s} Duration: {duration:2d} days, "
              f"ES: {ES[node]:2d}, EF: {EF[node]:2d}, Slack: {slack:2d}")
```

**Documentation:**
- https://networkx.org/documentation/stable/

**Common Gotchas:**
- `Graph()` is undirected, `DiGraph()` is directed; choose appropriately
- Node and edge attributes stored as dictionaries
- Some algorithms only work on directed or undirected graphs
- Large graphs can be memory-intensive
- Default drawing is basic; use specialized libraries for better visualization
- Shortest path algorithms require careful weight selection

---

## 7. FEA/CFD

### FEniCS (dolfin)

**Installation:**
```bash
# FEniCS requires specific installation
# Using conda (recommended):
conda create -n fenics -c conda-forge fenics
conda activate fenics

# Or using Docker:
# docker run -ti -v $(pwd):/home/fenics/shared quay.io/fenicsproject/stable
```

**Key Functions and Classes:**
- `dolfin.Mesh()`, `dolfin.UnitSquareMesh()` - Mesh generation
- `dolfin.FunctionSpace()` - Define function spaces
- `dolfin.TrialFunction()`, `dolfin.TestFunction()` - Variational formulation
- `dolfin.solve()` - Solve PDEs
- `dolfin.File()` - Save results
- Support for complex PDEs (heat, elasticity, Navier-Stokes)

**Engineering Applications:**
- Structural mechanics (FEA)
- Heat transfer analysis
- Fluid dynamics (CFD)
- Electromagnetic simulation
- Multi-physics coupling

**Code Examples:**

```python
# Note: FEniCS may not be installed in all environments
# These examples demonstrate the syntax and approach

try:
    from dolfin import *
    import numpy as np

    # Example 1: 2D heat conduction with steady-state
    # Solve: -∇²T = f on a square domain

    # Create mesh and function space
    mesh = UnitSquareMesh(32, 32)
    V = FunctionSpace(mesh, 'P', 1)  # Linear Lagrange elements

    # Define boundary conditions
    def boundary(x, on_boundary):
        return on_boundary

    T_boundary = Constant(0.0)  # Temperature at boundary
    bc = DirichletBC(V, T_boundary, boundary)

    # Define variational problem
    T = TrialFunction(V)
    v = TestFunction(V)
    f = Expression('10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)',
                   degree=2)  # Heat source
    k = Constant(1.0)  # Thermal conductivity

    a = k * dot(grad(T), grad(v)) * dx
    L = f * v * dx

    # Solve
    T_solution = Function(V)
    solve(a == L, T_solution, bc)

    # Save solution
    vtkfile = File('/tmp/temperature.pvd')
    vtkfile << T_solution

    print("Heat conduction solution saved to /tmp/temperature.pvd")
    print(f"Maximum temperature: {T_solution.vector().get_local().max():.2f}")

    # Example 2: Linear elasticity (2D stress analysis)
    # Solve for displacement under applied load

    # Create mesh and function space (vector)
    mesh_elastic = RectangleMesh(Point(0, 0), Point(1, 0.2), 40, 8)
    V_elastic = VectorFunctionSpace(mesh_elastic, 'P', 1)

    # Mark boundaries
    def left_boundary(x, on_boundary):
        return on_boundary and near(x[0], 0, 1e-10)

    def right_boundary(x, on_boundary):
        return on_boundary and near(x[0], 1, 1e-10)

    # Boundary conditions: fixed on left, traction on right
    bc_left = DirichletBC(V_elastic, Constant((0, 0)), left_boundary)

    # Material parameters (steel)
    E = 200e9  # Young's modulus (Pa)
    nu = 0.3   # Poisson's ratio
    mu = E / (2 * (1 + nu))
    lambda_ = E * nu / ((1 + nu) * (1 - 2 * nu))

    # Stress tensor
    def epsilon(u):
        return 0.5 * (nabla_grad(u) + nabla_grad(u).T)

    def sigma(u):
        return lambda_ * nabla_div(u) * Identity(2) + 2 * mu * epsilon(u)

    # Variational problem
    u = TrialFunction(V_elastic)
    v = TestFunction(V_elastic)
    f = Constant((0, 0))  # Body force
    T = Constant((0, -1e6))  # Traction on right edge (Pa)

    a = inner(sigma(u), epsilon(v)) * dx
    L = dot(f, v) * dx + dot(T, v) * ds

    # Solve
    u_solution = Function(V_elastic)
    solve(a == L, u_solution, bc_left)

    # Save solution
    vtkfile_elastic = File('/tmp/displacement.pvd')
    vtkfile_elastic << u_solution

    # Calculate von Mises stress
    s = sigma(u_solution)
    von_Mises = sqrt(3./2 * inner(dev(s), dev(s)))
    V_stress = FunctionSpace(mesh_elastic, 'P', 1)
    von_Mises_proj = project(von_Mises, V_stress)

    vtkfile_stress = File('/tmp/von_mises_stress.pvd')
    vtkfile_stress << von_Mises_proj

    print("\nElasticity solution saved to /tmp/displacement.pvd")
    print(f"Maximum displacement: {u_solution.vector().get_local().max()*1000:.3f} mm")
    print(f"Maximum von Mises stress: {von_Mises_proj.vector().get_local().max()/1e6:.2f} MPa")

    # Example 3: Transient heat equation
    # Solve: ρcp ∂T/∂t - k∇²T = f

    mesh_trans = UnitSquareMesh(32, 32)
    V_trans = FunctionSpace(mesh_trans, 'P', 1)

    # Time parameters
    dt = 0.01
    t_end = 2.0
    num_steps = int(t_end / dt)

    # Material properties
    rho = 8000.0  # kg/m³
    cp = 500.0    # J/(kg·K)
    k = 50.0      # W/(m·K)

    # Initial condition
    T_n = interpolate(Constant(20.0), V_trans)  # Initial temperature

    # Boundary condition
    bc_trans = DirichletBC(V_trans, Constant(20.0),
                          lambda x, on_boundary: on_boundary)

    # Variational formulation
    T = TrialFunction(V_trans)
    v = TestFunction(V_trans)
    f = Expression('100*exp(-t)', t=0, degree=1)  # Time-dependent source

    F = (rho*cp*T*v/dt + k*dot(grad(T), grad(v))) * dx - \
        (rho*cp*T_n*v/dt + f*v) * dx
    a, L = lhs(F), rhs(F)

    # Time-stepping
    T_solution = Function(V_trans)
    t = 0

    # Save initial condition
    vtkfile_trans = File('/tmp/transient_heat.pvd')
    vtkfile_trans << (T_n, t)

    for n in range(num_steps):
        t += dt
        f.t = t

        solve(a == L, T_solution, bc_trans)
        T_n.assign(T_solution)

        # Save every 10 steps
        if n % 10 == 0:
            vtkfile_trans << (T_solution, t)

    print("\nTransient heat solution saved to /tmp/transient_heat.pvd")
    print(f"Final time: {t:.2f} s")
    print(f"Final maximum temperature: {T_solution.vector().get_local().max():.2f}°C")

except ImportError:
    print("FEniCS (dolfin) is not installed in this environment.")
    print("FEniCS requires special installation via conda or Docker.")
    print("\nTo install:")
    print("  conda create -n fenics -c conda-forge fenics")
    print("  conda activate fenics")
    print("\nFor more information: https://fenicsproject.org/download/")
```

**Documentation:**
- https://fenicsproject.org/documentation/

**Common Gotchas:**
- Installation can be complex; use conda or Docker
- Requires understanding of variational formulation
- Mesh quality significantly affects results
- Boundary conditions must be properly specified
- Results saved in .pvd format (use ParaView to visualize)
- Time-dependent problems require careful time-stepping
- Memory usage can be significant for 3D problems

---

## 8. Aerospace

### poliastro

**Installation:**
```bash
pip install poliastro
```

**Key Functions and Classes:**
- `poliastro.bodies` - Planetary bodies (Earth, Mars, etc.)
- `poliastro.twobody.Orbit` - Define orbits
- `poliastro.maneuver` - Orbital maneuvers (Hohmann, bi-elliptic)
- `poliastro.plotting` - Orbit visualization
- `poliastro.iod` - Initial orbit determination
- Propagation and perturbation analysis

**Engineering Applications:**
- Satellite orbit design
- Mission planning
- Trajectory optimization
- Orbit determination
- Launch window analysis
- Interplanetary transfers

**Code Examples:**

```python
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro import constants
from astropy import units as u
from astropy.time import Time
import numpy as np

# Example 1: Low Earth Orbit (LEO) satellite analysis

# Define orbit parameters
altitude = 400 * u.km
a = Earth.R + altitude  # Semi-major axis
ecc = 0.001 * u.one     # Eccentricity (nearly circular)
inc = 51.6 * u.deg      # Inclination (ISS-like)
raan = 0 * u.deg        # Right ascension of ascending node
argp = 0 * u.deg        # Argument of periapsis
nu = 0 * u.deg          # True anomaly

# Create orbit
epoch = Time("2024-01-01 12:00:00", scale="utc")
leo = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu, epoch=epoch)

print("LEO Satellite Orbit:")
print(f"  Altitude: {altitude}")
print(f"  Period: {leo.period.to(u.minute):.2f}")
print(f"  Velocity: {leo.v.value[0]:.2f} km/s")
print(f"  Energy: {leo.energy.to(u.MJ/u.kg):.2f}")

# Calculate ground track parameters
angular_velocity_earth = 360 * u.deg / (23.9344696 * u.hour)
orbital_velocity_angular = 360 * u.deg / leo.period
drift_per_orbit = (orbital_velocity_angular - angular_velocity_earth) * leo.period
print(f"  Ground track drift: {drift_per_orbit.to(u.deg):.2f} per orbit")

# Example 2: Hohmann transfer orbit
# Transfer from LEO to GEO

# Initial orbit (LEO)
r1 = Earth.R + 400 * u.km
orbit_leo = Orbit.circular(Earth, alt=400 * u.km, epoch=epoch)

# Target orbit (GEO)
r2 = 42164 * u.km  # GEO altitude from center of Earth
orbit_geo = Orbit.circular(Earth, alt=(r2 - Earth.R), epoch=epoch)

print("\nHohmann Transfer (LEO to GEO):")

# Calculate Hohmann transfer velocities
v1 = np.sqrt(Earth.k / r1)  # Velocity in LEO
a_transfer = (r1 + r2) / 2  # Transfer orbit semi-major axis
v_transfer_1 = np.sqrt(Earth.k * (2/r1 - 1/a_transfer))  # Velocity at perigee
v_transfer_2 = np.sqrt(Earth.k * (2/r2 - 1/a_transfer))  # Velocity at apogee
v2 = np.sqrt(Earth.k / r2)  # Velocity in GEO

# Delta-V requirements
dv1 = v_transfer_1 - v1
dv2 = v2 - v_transfer_2
total_dv = dv1 + dv2

print(f"  First burn (LEO): {dv1.to(u.km/u.s):.3f}")
print(f"  Second burn (GEO): {dv2.to(u.km/u.s):.3f}")
print(f"  Total delta-V: {total_dv.to(u.km/u.s):.3f}")

# Transfer time
t_transfer = np.pi * np.sqrt(a_transfer**3 / Earth.k)
print(f"  Transfer time: {t_transfer.to(u.hour):.2f}")

# Propellant mass calculation (rocket equation)
Isp = 300 * u.s  # Specific impulse
g0 = constants.g0
m_initial = 1000 * u.kg  # Initial spacecraft mass

# Mass ratio
mass_ratio = np.exp(total_dv / (Isp * g0))
m_final = m_initial / mass_ratio
m_propellant = m_initial - m_final

print(f"  Propellant required: {m_propellant.to(u.kg):.1f} ({m_propellant/m_initial*100:.1f}% of initial mass)")

# Example 3: Interplanetary trajectory (Earth to Mars)

# Define departure and arrival orbits
orbit_earth = Orbit.circular(Sun, alt=Earth.R, epoch=epoch)
orbit_mars = Orbit.circular(Sun, alt=Mars.R, epoch=epoch)

print("\nEarth to Mars Transfer:")

# Approximate Hohmann transfer (simplified)
r_earth = 1.0 * u.AU
r_mars = 1.524 * u.AU

a_transfer_mars = (r_earth + r_mars) / 2
v_earth_orbit = np.sqrt(Sun.k / r_earth)
v_transfer_earth = np.sqrt(Sun.k * (2/r_earth - 1/a_transfer_mars))
v_mars_orbit = np.sqrt(Sun.k / r_mars)
v_transfer_mars = np.sqrt(Sun.k * (2/r_mars - 1/a_transfer_mars))

# Delta-V at Earth
dv_earth = abs(v_transfer_earth - v_earth_orbit)
# Delta-V at Mars
dv_mars = abs(v_mars_orbit - v_transfer_mars)

print(f"  Delta-V at Earth: {dv_earth.to(u.km/u.s):.2f}")
print(f"  Delta-V at Mars: {dv_mars.to(u.km/u.s):.2f}")

# Transfer time
t_transfer_mars = np.pi * np.sqrt(a_transfer_mars**3 / Sun.k)
print(f"  Transfer time: {t_transfer_mars.to(u.day):.0f} ({t_transfer_mars.to(u.month):.1f})")

# Synodic period (launch window frequency)
T_earth = 2 * np.pi * np.sqrt(r_earth**3 / Sun.k)
T_mars = 2 * np.pi * np.sqrt(r_mars**3 / Sun.k)
T_synodic = abs(1 / (1/T_earth - 1/T_mars))
print(f"  Launch window frequency: {T_synodic.to(u.day):.0f} ({T_synodic.to(u.year):.2f})")

# Example 4: Satellite coverage analysis
# Calculate field of view and access time

satellite_altitude = 600 * u.km
ground_station_lat = 40 * u.deg  # Latitude
min_elevation = 10 * u.deg       # Minimum elevation angle for communication

# Calculate maximum range
earth_radius = Earth.R.to(u.km)
sat_radius = earth_radius + satellite_altitude

# Maximum slant range
max_range = np.sqrt(sat_radius**2 - (earth_radius * np.cos(min_elevation))**2) - \
            earth_radius * np.sin(min_elevation)

# Maximum Earth central angle
max_angle = np.arccos(earth_radius * np.cos(min_elevation) / sat_radius)

print(f"\nSatellite Coverage Analysis:")
print(f"  Satellite altitude: {satellite_altitude}")
print(f"  Maximum slant range: {max_range:.0f}")
print(f"  Maximum Earth central angle: {(max_angle * u.rad).to(u.deg):.2f}")

# Coverage area
coverage_area = 2 * np.pi * earth_radius**2 * (1 - np.cos(max_angle))
earth_area = 4 * np.pi * earth_radius**2
coverage_percent = (coverage_area / earth_area * 100)

print(f"  Coverage area: {coverage_area.to(u.km**2):.2e}")
print(f"  Percentage of Earth: {coverage_percent:.2f}%")

# Access time per orbit (approximate)
orbit_period = leo.period.to(u.minute)
access_angle = 2 * max_angle * u.rad
orbit_angle = 2 * np.pi * u.rad
access_time = (access_angle / orbit_angle) * orbit_period

print(f"  Access time per pass: {access_time:.2f}")
```

**Documentation:**
- https://docs.poliastro.space/

**Common Gotchas:**
- Requires astropy units; always use units with parameters
- Coordinate systems and reference frames must be consistent
- Orbital elements have singularities (e.g., circular orbits)
- Perturbations (J2, drag) require additional propagators
- Time scales (UTC, TDB) matter for precision calculations
- Some calculations are approximate for complex scenarios

---

## 9. High Precision

### mpmath

**Installation:**
```bash
pip install mpmath
```

**Key Functions and Classes:**
- `mpmath.mp.dps` - Set decimal precision
- `mpmath.mpf()` - Arbitrary precision float
- `mpmath.matrix()` - Arbitrary precision matrices
- Standard math functions with arbitrary precision
- Special functions (Bessel, gamma, etc.)
- Numerical integration and differentiation

**Engineering Applications:**
- Ill-conditioned numerical problems
- Verification of numerical algorithms
- Accumulation of rounding errors
- Symbolic-numeric computations
- Special function evaluation
- Critical calculations requiring high accuracy

**Code Examples:**

```python
from mpmath import mp, mpf, matrix, exp, sin, cos, pi, sqrt, quad
import numpy as np

# Example 1: Solving ill-conditioned linear system
# Hilbert matrix is notoriously ill-conditioned

# Set precision to 50 decimal places
mp.dps = 50

n = 10  # Matrix size

# Create Hilbert matrix with high precision
H = matrix(n, n)
for i in range(n):
    for j in range(n):
        H[i,j] = mpf(1) / (mpf(i + 1) + mpf(j + 1) - mpf(1))

# Right-hand side vector
b = matrix(n, 1)
for i in range(n):
    b[i] = mpf(1)

# Solve using high precision
x = mp.lu_solve(H, b)

print("High-Precision Linear System Solution:")
print(f"Matrix size: {n}x{n}")
print(f"Precision: {mp.dps} decimal places")
print(f"First solution component: {x[0]}")

# Verify solution by computing residual
residual = H * x - b
residual_norm = sqrt(sum(r**2 for r in residual))
print(f"Residual norm: {residual_norm}")

# Compare with standard numpy (will show numerical issues)
H_numpy = np.array([[1.0/(i+j+1) for j in range(n)] for i in range(n)])
b_numpy = np.ones(n)
try:
    x_numpy = np.linalg.solve(H_numpy, b_numpy)
    residual_numpy = np.linalg.norm(H_numpy @ x_numpy - b_numpy)
    print(f"NumPy residual norm: {residual_numpy:.2e} (shows loss of precision)")
except np.linalg.LinAlgError:
    print("NumPy failed to solve due to ill-conditioning")

# Example 2: Calculating critical engineering constants with high precision
# Sometimes intermediate calculations require high precision

mp.dps = 100  # Increase precision to 100 decimal places

# Calculate gravitational parameter for orbital mechanics
# μ = G * M, where precision matters for long-term orbit prediction

G = mpf('6.67430e-11')  # Gravitational constant (m³/kg/s²)
M_earth = mpf('5.972e24')  # Earth mass (kg)
mu_earth = G * M_earth

print(f"\nHigh-Precision Gravitational Parameter:")
print(f"G = {G}")
print(f"M_earth = {M_earth}")
print(f"μ_earth = {mu_earth} m³/s²")

# Calculate orbital period with high precision
a = mpf('42164000')  # GEO altitude from center (m)
T = mpf(2) * pi * sqrt(a**3 / mu_earth)

print(f"GEO orbital period: {T} seconds")
print(f"GEO orbital period: {T/3600} hours")

# Small differences matter for long-term predictions
# Calculate position error after 1 year due to precision
time_year = mpf('31557600')  # seconds in a year
n_orbits = time_year / T
error_accumulation = n_orbits * mpf('1e-10')  # Hypothetical error per orbit

print(f"Orbits per year: {n_orbits}")
print(f"Accumulated error: {error_accumulation} m")

# Example 3: Evaluating special functions for engineering formulas
# Bessel functions in heat transfer, vibrations, electromagnetics

mp.dps = 30

# Heat transfer in a cylindrical fin
# Temperature distribution involves Bessel functions
# T(r) = T0 * I0(m*r) / I0(m*R)

from mpmath import besseli, besselk

# Parameters
m = mpf('10')  # Parameter related to fin geometry and material
r_values = [mpf(i)/mpf(10) for i in range(11)]  # Radial positions

print(f"\nBessel Function Application (Cylindrical Fin):")
print(f"m = {m}")
print(f"\nr/R\tI0(m*r/R)/I0(m)")

for r in r_values:
    ratio = besseli(0, m * r) / besseli(0, m)
    print(f"{float(r):.1f}\t{float(ratio):.6f}")

# Example 4: Numerical integration with adaptive precision
# Calculate work done by variable force

mp.dps = 30

# Force varies as F(x) = A * exp(-x²/L²) * sin(x)
A = mpf('1000')  # N
L = mpf('0.5')   # m

def force(x):
    return A * exp(-(x**2)/(L**2)) * sin(x)

# Integrate from 0 to 2 meters
x_start = mpf('0')
x_end = mpf('2')

work, error = quad(force, [x_start, x_end], error=True)

print(f"\nNumerical Integration (Work calculation):")
print(f"Work done: {work} J")
print(f"Estimated error: {error} J")
print(f"Relative error: {float(error/work)*100:.2e}%")

# Example 5: Catastrophic cancellation demonstration
# Subtracting nearly equal numbers loses precision

mp.dps = 50

# Calculate (1 + x) - 1 for very small x
# This demonstrates why high precision matters

x_small = mpf('1e-20')

# Method 1: Direct calculation (catastrophic cancellation with limited precision)
result_direct = (mpf(1) + x_small) - mpf(1)

# Method 2: Using high precision
result_precise = x_small

print(f"\nCatastrophic Cancellation Example:")
print(f"x = {x_small}")
print(f"(1 + x) - 1 = {result_direct}")
print(f"x (direct) = {result_precise}")
print(f"Difference preserved: {result_direct == result_precise}")

# Real engineering example: stress concentration factor
# K_t = 1 + 2*sqrt(a/r) where a is crack length, r is radius
# For a/r << 1, need careful calculation

a = mpf('1e-6')  # 1 micrometer crack
r = mpf('1e-3')  # 1 millimeter radius

K_t = mpf(1) + mpf(2) * sqrt(a/r)

print(f"\nStress Concentration Factor:")
print(f"Crack length: {a} m")
print(f"Radius: {r} m")
print(f"K_t = {K_t}")
print(f"K_t - 1 = {K_t - mpf(1)} (small but important)")
```

**Documentation:**
- https://mpmath.org/doc/current/

**Common Gotchas:**
- Set `mp.dps` before calculations; doesn't affect existing values
- Much slower than standard floating point (trade precision for speed)
- Convert to float/int for output if needed
- String input for exact decimal values: `mpf('1.23')`
- Mixed precision operations reduce to lower precision
- Not all NumPy functions have mpmath equivalents
- Memory usage increases with precision

---

## Summary

This guide covers the essential Python packages for engineering applications:

1. **Numerical Computing**: NumPy, SciPy, SymPy - Foundation for all numerical work
2. **Fluid Mechanics**: fluids, thermo, CoolProp - Specialized thermodynamic and flow calculations
3. **Optimization**: scipy.optimize, DEAP - Design optimization and parameter estimation
4. **Visualization**: Matplotlib, PyVista, Plotly - From basic plots to interactive 3D
5. **Units**: Pint, astropy.units - Prevent unit errors and document assumptions
6. **Network Analysis**: NetworkX - Infrastructure and system connectivity
7. **FEA/CFD**: FEniCS - Finite element analysis for PDEs
8. **Aerospace**: poliastro - Orbital mechanics and mission planning
9. **High Precision**: mpmath - Critical calculations requiring arbitrary precision

**General Best Practices:**

- Always validate results with hand calculations or known solutions
- Use version control for code and document package versions
- Check units and dimensional consistency
- Verify convergence for iterative methods
- Consider numerical stability and precision requirements
- Use appropriate visualization for results verification
- Document assumptions and limitations
- Test edge cases and limiting conditions

**Resources:**

- Python Package Index (PyPI): https://pypi.org/
- Conda Forge: https://conda-forge.org/
- Engineering Stack on Stack Overflow
- SciPy Conference proceedings
- Individual package documentation (linked above)

This documentation provides working code examples that can be adapted for specific engineering problems. Always refer to the official documentation for the most up-to-date information and advanced features.
