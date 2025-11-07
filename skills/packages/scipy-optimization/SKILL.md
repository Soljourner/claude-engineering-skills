---
name: scipy-optimization
description: "Optimize pump designs and system parameters using scipy.optimize"
category: packages
domain: general
complexity: intermediate
dependencies:
  - scipy
  - numpy
---

# SciPy Optimization for Engineering Design

## Overview

The `scipy.optimize` module provides a comprehensive suite of optimization algorithms for solving engineering design problems. This skill focuses on applying these methods to pump design optimization, system parameter tuning, and performance analysis.

Key capabilities:
- **Unconstrained optimization**: Find optimal parameters without restrictions
- **Constrained optimization**: Optimize with design constraints and bounds
- **Curve fitting**: Fit models to experimental data
- **Global optimization**: Find global minima in multi-modal problems
- **Multi-objective optimization**: Balance competing design objectives

## Optimization Methods

### 1. minimize() - Unconstrained and Constrained Optimization

The workhorse function for most optimization tasks.

```python
from scipy.optimize import minimize
import numpy as np

# Basic usage
def objective(x):
    """Objective function to minimize"""
    return x[0]**2 + x[1]**2

result = minimize(objective, x0=[1.0, 1.0])
print(f"Optimal point: {result.x}")
print(f"Optimal value: {result.fun}")
```

**Common algorithms:**
- `'Nelder-Mead'`: Derivative-free, robust but slower
- `'BFGS'`: Quasi-Newton method, fast for smooth functions
- `'L-BFGS-B'`: BFGS with bounds
- `'SLSQP'`: Sequential Least Squares with constraints
- `'trust-constr'`: General constrained optimization (preferred)

### 2. least_squares() - Curve Fitting and Residual Minimization

Specialized for problems of the form: minimize sum(residuals^2)

```python
from scipy.optimize import least_squares

def residuals(params, x_data, y_data):
    """Calculate residuals between model and data"""
    a, b, c = params
    y_model = a * x_data**2 + b * x_data + c
    return y_model - y_data

# Fit quadratic to data
x_data = np.array([1, 2, 3, 4, 5])
y_data = np.array([2.1, 3.9, 9.2, 15.8, 25.1])
result = least_squares(residuals, x0=[1, 1, 1], args=(x_data, y_data))
```

**Key features:**
- Handles bounds on variables
- Robust loss functions (handles outliers)
- Efficient for large residual problems
- Returns Jacobian information

### 3. differential_evolution() - Global Optimization

Finds global minimum using evolutionary algorithms. Essential for multi-modal problems.

```python
from scipy.optimize import differential_evolution

def complex_function(x):
    """Function with multiple local minima"""
    return np.sin(x[0]) * np.cos(x[1]) + (x[0] - 1)**2 + (x[1] + 2)**2

# Search over bounds
bounds = [(-5, 5), (-5, 5)]
result = differential_evolution(complex_function, bounds)
```

**When to use:**
- Multiple local minima expected
- No gradient information available
- Need robust global solution
- Can afford longer computation time

### 4. Constrained Optimization

Optimize subject to equality and inequality constraints.

```python
from scipy.optimize import minimize, NonlinearConstraint, LinearConstraint

def objective(x):
    return x[0]**2 + x[1]**2

# Equality constraint: x[0] + x[1] = 1
def eq_constraint(x):
    return x[0] + x[1] - 1

# Inequality constraint: x[0]^2 + x[1]^2 >= 0.5
def ineq_constraint(x):
    return x[0]**2 + x[1]**2 - 0.5

constraints = [
    {'type': 'eq', 'fun': eq_constraint},
    {'type': 'ineq', 'fun': ineq_constraint}
]

result = minimize(objective, x0=[0.5, 0.5],
                 method='SLSQP', constraints=constraints)
```

## Engineering Applications

### 1. Pump Efficiency Maximization

Optimize pump geometry (impeller diameter, blade angle, width) to maximize efficiency at design point.

**Problem formulation:**
- **Objective**: Maximize efficiency η(D, β, b)
- **Constraints**:
  - Minimum flow rate requirement
  - Maximum power consumption
  - Geometric constraints (e.g., D_min ≤ D ≤ D_max)
  - NPSH requirements

**Approach:**
```python
def pump_efficiency(params):
    """Calculate pump efficiency (negative for minimization)"""
    D, beta, b = params  # diameter, blade angle, width

    # Simplified efficiency model
    Q = calculate_flow_rate(D, beta, b)
    H = calculate_head(D, beta, b)
    P = calculate_power(D, beta, b)

    eta = (rho * g * Q * H) / P
    return -eta  # Negative for maximization

bounds = [(0.1, 0.5), (20, 45), (0.02, 0.1)]  # D, beta, b
result = minimize(pump_efficiency, x0=[0.3, 30, 0.05],
                 method='L-BFGS-B', bounds=bounds)
```

### 2. Impeller Geometry Optimization

Optimize blade profile, shroud contour, and hub geometry for specific duty point.

**Key considerations:**
- Velocity triangles at inlet/outlet
- Blade loading distribution
- Cavitation performance (NPSH)
- Structural integrity (stress limits)

### 3. System Cost Minimization

Find optimal pump size and operating speed to minimize total lifecycle cost.

**Cost components:**
- Capital cost (pump + motor + VFD)
- Energy cost (operating hours × efficiency)
- Maintenance cost
- Replacement cost (NPV)

**Formulation:**
```python
def total_cost(params):
    """Total lifecycle cost"""
    pump_size, speed_rpm = params

    # Capital cost
    capital = cost_pump(pump_size) + cost_motor(pump_size, speed_rpm)

    # Operating cost (20 year life)
    annual_energy = operating_hours * power(pump_size, speed_rpm) * electricity_rate
    operating = annual_energy * 20

    # Maintenance
    maintenance = maintenance_rate * capital * 20

    return capital + operating + maintenance
```

### 4. Multi-Objective Optimization

Balance competing objectives (e.g., efficiency vs cost, flow vs head).

**Pareto front approach:**
```python
def weighted_objectives(params, weight):
    """Combine objectives with weighting"""
    eta = pump_efficiency(params)
    cost = pump_cost(params)

    # Normalize and combine
    return -weight * eta + (1 - weight) * cost

# Sweep weights to find Pareto front
pareto_solutions = []
for w in np.linspace(0, 1, 21):
    result = minimize(lambda x: weighted_objectives(x, w),
                     x0=[0.3, 30, 0.05], bounds=bounds)
    pareto_solutions.append(result.x)
```

### 5. Curve Fitting to Experimental Data

Fit pump characteristic curves (H-Q, η-Q, P-Q) to test data.

**Applications:**
- Develop performance models from test data
- Extrapolate performance to off-design conditions
- Identify empirical coefficients in analytical models

**Example:**
```python
from scipy.optimize import curve_fit

def pump_curve(Q, a, b, c):
    """Quadratic head-flow relationship"""
    return a - b * Q - c * Q**2

# Fit to test data
Q_test = np.array([0, 50, 100, 150, 200])  # m³/h
H_test = np.array([120, 115, 105, 90, 70])  # m

params, covariance = curve_fit(pump_curve, Q_test, H_test)
print(f"H = {params[0]:.2f} - {params[1]:.4f}*Q - {params[2]:.6f}*Q²")
```

## Best Practices

### 1. Problem Formulation
- **Scale variables**: Normalize to similar magnitudes (0-1 or -1 to 1)
- **Choose good initial guess**: Start near expected solution
- **Define clear bounds**: Prevent physically impossible solutions
- **Formulate smooth objectives**: Avoid discontinuities

### 2. Algorithm Selection
- **Smooth, gradient-available**: Use `'BFGS'` or `'L-BFGS-B'`
- **With constraints**: Use `'SLSQP'` or `'trust-constr'`
- **Noisy/discontinuous**: Use `'Nelder-Mead'` or derivative-free methods
- **Multiple local minima**: Use `differential_evolution` first, then refine
- **Curve fitting**: Use `curve_fit` or `least_squares`

### 3. Convergence Tips
- Try multiple initial guesses
- Use reasonable tolerances (not too tight)
- Provide analytical gradients if possible
- Check optimization success flag
- Validate physical feasibility of results

### 4. Constraint Formulation
- Express as g(x) ≥ 0 for inequality constraints
- Express as h(x) = 0 for equality constraints
- Keep constraints smooth and differentiable
- Avoid redundant constraints

## Common Pitfalls

1. **Poor scaling**: Variables of vastly different magnitudes
2. **Bad initial guess**: Starting far from feasible region
3. **Wrong algorithm**: Using gradient-based methods on non-smooth functions
4. **Tight tolerances**: Asking for more precision than achievable
5. **Ignoring bounds**: Not constraining to physical limits
6. **No validation**: Not checking if solution makes engineering sense

## Performance Optimization Tips

1. **Vectorize operations**: Use NumPy arrays instead of loops
2. **Provide Jacobians**: Analytical gradients are faster and more accurate
3. **Use appropriate method**: Don't use global optimization if local is sufficient
4. **Cache expensive computations**: Store intermediate results
5. **Start with simple model**: Refine complexity iteratively

## Integration with Pump Design Workflow

```
1. Define design variables (D, β, b, ω, etc.)
2. Establish objective function (maximize η, minimize cost)
3. Define constraints (Q_min, NPSH, stress limits)
4. Set bounds (physical limits)
5. Choose optimization algorithm
6. Run optimization with multiple initial guesses
7. Validate results (CFD, analytical checks)
8. Perform sensitivity analysis
9. Generate documentation and plots
```

## Example Workflow

```python
import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib.pyplot as plt

# Step 1: Define problem
def objective(x):
    D, beta = x
    return -calculate_efficiency(D, beta)

def constraint_flow(x):
    D, beta = x
    return calculate_flow_rate(D, beta) - Q_min

def constraint_power(x):
    D, beta = x
    return P_max - calculate_power(D, beta)

# Step 2: Set up optimization
bounds = [(0.2, 0.6), (20, 50)]
constraints = [
    {'type': 'ineq', 'fun': constraint_flow},
    {'type': 'ineq', 'fun': constraint_power}
]

# Step 3: Solve with global optimizer first
result_global = differential_evolution(
    objective, bounds, constraints=constraints
)

# Step 4: Refine with local optimizer
result_local = minimize(
    objective, x0=result_global.x,
    method='SLSQP', bounds=bounds,
    constraints=constraints
)

# Step 5: Report results
print(f"Optimal diameter: {result_local.x[0]:.3f} m")
print(f"Optimal blade angle: {result_local.x[1]:.1f}°")
print(f"Maximum efficiency: {-result_local.fun:.2%}")
```

## References

- [SciPy Optimization Documentation](https://docs.scipy.org/doc/scipy/reference/optimize.html)
- [Optimization Methods Tutorial](https://scipy-lectures.org/advanced/mathematical_optimization/)
- Engineering optimization textbooks for problem formulation
- Pump design handbooks for realistic constraints and models
