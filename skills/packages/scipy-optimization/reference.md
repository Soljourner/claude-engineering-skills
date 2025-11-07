# SciPy Optimization Reference Guide

## Quick Algorithm Selection Guide

### Decision Tree

```
START
│
├─ Do you have explicit objective function?
│  ├─ NO → Use curve_fit() or least_squares() for data fitting
│  └─ YES ↓
│
├─ Are there multiple local minima?
│  ├─ YES → Use differential_evolution() or basinhopping()
│  └─ NO ↓
│
├─ Do you have constraints?
│  ├─ NO → Unconstrained optimization
│  │  ├─ Gradient available → minimize() with 'BFGS'
│  │  └─ No gradient → minimize() with 'Nelder-Mead'
│  └─ YES ↓
│
└─ Constrained optimization
   ├─ Only bounds → minimize() with 'L-BFGS-B'
   ├─ Linear constraints → minimize() with 'trust-constr'
   └─ Nonlinear constraints → minimize() with 'SLSQP' or 'trust-constr'
```

## Algorithm Comparison Table

| Algorithm | Type | Constraints | Gradient | Speed | Robustness | Use Case |
|-----------|------|-------------|----------|-------|------------|----------|
| `Nelder-Mead` | Local | Bounds only | No | Medium | High | Noisy/discontinuous functions |
| `Powell` | Local | Bounds only | No | Medium | Medium | Separable problems |
| `BFGS` | Local | None | Yes* | Fast | Medium | Smooth unconstrained |
| `L-BFGS-B` | Local | Bounds | Yes* | Fast | Medium | Large-scale with bounds |
| `SLSQP` | Local | All types | Yes* | Fast | Medium | General constrained |
| `trust-constr` | Local | All types | Yes* | Medium | High | Preferred for constraints |
| `differential_evolution` | Global | Bounds | No | Slow | Very High | Multi-modal problems |
| `basinhopping` | Global | Bounds | No | Slow | High | Complex landscapes |
| `least_squares` | Local | Bounds | Yes* | Fast | High | Sum of squares |
| `curve_fit` | Local | Bounds | Yes* | Fast | High | Data fitting |

\* Can use numerical gradients if analytical not provided

## Detailed Algorithm Descriptions

### 1. Nelder-Mead

**When to use:**
- No gradient information available
- Function is noisy or has discontinuities
- Small to medium dimensionality (< 20 variables)
- Quick prototyping

**Strengths:**
- Very robust
- No gradient required
- Good for noisy functions

**Weaknesses:**
- Slow convergence
- May fail in high dimensions
- No constraint handling (except bounds in special versions)

**Example:**
```python
result = minimize(objective, x0, method='Nelder-Mead',
                 options={'maxiter': 1000, 'xatol': 1e-6, 'fatol': 1e-6})
```

### 2. BFGS (Broyden-Fletcher-Goldfarb-Shanno)

**When to use:**
- Smooth, unconstrained problems
- Gradient available or can be approximated
- Medium scale (10-1000 variables)

**Strengths:**
- Fast convergence for smooth functions
- Superlinear convergence rate
- Widely used and reliable

**Weaknesses:**
- Requires gradient (numerical gradient is slow)
- No constraint handling
- Memory intensive for very large problems

**Example:**
```python
result = minimize(objective, x0, method='BFGS', jac=gradient)
```

### 3. L-BFGS-B (Limited-memory BFGS with Bounds)

**When to use:**
- Large-scale optimization (1000+ variables)
- Simple bounds on variables
- Gradient available

**Strengths:**
- Memory efficient
- Handles bounds
- Good for large problems

**Weaknesses:**
- Only handles box constraints
- Requires smooth objective

**Example:**
```python
bounds = [(0, 1), (0, 2), (-np.inf, np.inf)]
result = minimize(objective, x0, method='L-BFGS-B', bounds=bounds)
```

### 4. SLSQP (Sequential Least Squares Programming)

**When to use:**
- Small to medium problems with constraints
- Mix of equality and inequality constraints
- Quick implementation needed

**Strengths:**
- Handles all constraint types
- Relatively fast
- Easy to use

**Weaknesses:**
- Less robust than trust-constr
- Can fail with tight tolerances
- Sensitive to scaling

**Example:**
```python
constraints = [
    {'type': 'eq', 'fun': lambda x: x[0] + x[1] - 1},
    {'type': 'ineq', 'fun': lambda x: x[0] - 0.5}
]
result = minimize(objective, x0, method='SLSQP', constraints=constraints)
```

### 5. trust-constr (Trust Region Constrained)

**When to use:**
- Complex constrained optimization
- Need robust convergence
- Can provide Hessian or Jacobian

**Strengths:**
- Most robust for constrained problems
- Handles all constraint types
- Good convergence properties

**Weaknesses:**
- Slower than SLSQP
- More complex to set up
- Requires recent scipy version

**Example:**
```python
from scipy.optimize import NonlinearConstraint

constraint = NonlinearConstraint(lambda x: x[0]**2 + x[1]**2, 0.5, 2.0)
result = minimize(objective, x0, method='trust-constr', constraints=[constraint])
```

### 6. differential_evolution

**When to use:**
- Global optimization needed
- Multiple local minima suspected
- No gradient available
- Robust solution required

**Strengths:**
- Finds global minimum reliably
- Very robust
- No gradient needed
- Good for discontinuous functions

**Weaknesses:**
- Slow (many function evaluations)
- Requires bounds
- Stochastic (different runs give different results)

**Example:**
```python
bounds = [(0, 10), (0, 10)]
result = differential_evolution(objective, bounds, seed=42, maxiter=100)
```

### 7. least_squares

**When to use:**
- Fitting models to data
- Minimizing sum of squared residuals
- Robust regression needed

**Strengths:**
- Optimized for residual problems
- Robust loss functions available
- Efficient for sparse Jacobians
- Handles bounds

**Weaknesses:**
- Only for least squares problems
- No general constraints

**Example:**
```python
def residuals(params, x, y):
    return model(params, x) - y

result = least_squares(residuals, x0, args=(x_data, y_data),
                      loss='soft_l1', bounds=(lb, ub))
```

## Convergence Tips

### 1. Scaling Variables

**Problem:** Variables with different magnitudes (e.g., D=0.3m, omega=1800 RPM)

**Solution:** Normalize to similar ranges

```python
# Bad: Different scales
x = [D, omega]  # [0.3, 1800]

# Good: Normalized
x_normalized = [D/0.5, omega/3000]  # [0.6, 0.6]

def objective_normalized(x_norm):
    D = x_norm[0] * 0.5
    omega = x_norm[1] * 3000
    return objective([D, omega])
```

### 2. Initial Guess

**Tips:**
- Start near expected solution
- Use engineering judgment
- Try multiple starting points
- Use prior designs as starting point

```python
# Try multiple starting points
initial_guesses = [
    [0.25, 1800, 0.05],  # Conservative design
    [0.30, 1500, 0.06],  # Larger, slower
    [0.20, 2400, 0.04],  # Smaller, faster
]

results = [minimize(objective, x0=ig, method='SLSQP') for ig in initial_guesses]
best_result = min(results, key=lambda r: r.fun)
```

### 3. Tolerance Settings

**Default tolerances may be too tight or too loose**

```python
# Tight tolerance (slow, high precision)
options = {'ftol': 1e-12, 'gtol': 1e-12}

# Loose tolerance (fast, engineering accuracy)
options = {'ftol': 1e-6, 'gtol': 1e-6}

# For most engineering: 1e-6 to 1e-9 is sufficient
result = minimize(objective, x0, method='BFGS', options={'ftol': 1e-8})
```

### 4. Providing Gradients

**Analytical gradients are much faster and more accurate**

```python
def objective(x):
    return x[0]**2 + 2*x[1]**2

def gradient(x):
    return np.array([2*x[0], 4*x[1]])

# Much faster with analytical gradient
result = minimize(objective, x0, jac=gradient, method='BFGS')
```

### 5. Checking Convergence

**Always verify optimization succeeded**

```python
result = minimize(objective, x0, method='SLSQP')

if not result.success:
    print(f"Optimization failed: {result.message}")
    print(f"Status: {result.status}")

# Check optimality conditions
if hasattr(result, 'jac'):
    grad_norm = np.linalg.norm(result.jac)
    print(f"Gradient norm: {grad_norm:.2e}")
    if grad_norm > 1e-4:
        print("Warning: Large gradient at solution")
```

### 6. Dealing with Failed Convergence

**Common issues and fixes:**

| Issue | Symptom | Solution |
|-------|---------|----------|
| Poor scaling | Slow convergence | Normalize variables |
| Bad initial guess | Fails immediately | Try different x0 |
| Tight constraints | Infeasible | Check constraint formulation |
| Discontinuities | Erratic behavior | Use gradient-free method |
| Local minimum | Sub-optimal result | Use global optimizer |
| Numerical noise | No convergence | Increase tolerance |

## Constraint Formulation

### Types of Constraints

**1. Bounds (Box Constraints)**
```python
# Simple variable limits
bounds = [
    (D_min, D_max),     # 0.1 ≤ D ≤ 0.5
    (omega_min, omega_max),  # 1000 ≤ ω ≤ 3000
]
```

**2. Linear Constraints**
```python
from scipy.optimize import LinearConstraint

# a·x ≤ b  (inequality)
# Example: x[0] + 2*x[1] ≤ 10
A = np.array([[1, 2]])
lb = -np.inf
ub = 10
constraint = LinearConstraint(A, lb, ub)
```

**3. Nonlinear Equality Constraints**
```python
# h(x) = 0
def equality_constraint(x):
    return x[0]**2 + x[1]**2 - 1  # Circle: x² + y² = 1

constraint = {'type': 'eq', 'fun': equality_constraint}
```

**4. Nonlinear Inequality Constraints**
```python
# g(x) ≥ 0  (Note: ≥ 0, not ≤ 0)
def inequality_constraint(x):
    return x[0] - x[1]**2  # x[0] ≥ x[1]²

constraint = {'type': 'ineq', 'fun': inequality_constraint}
```

### Best Practices for Constraints

**1. Formulate as g(x) ≥ 0**

```python
# Requirement: Flow rate must be at least 100 L/s
# Bad
def constraint(x):
    Q = calculate_flow(x)
    return 100 - Q  # Wrong sign!

# Good
def constraint(x):
    Q = calculate_flow(x)
    return Q - 100  # Correct: Q ≥ 100
```

**2. Make Constraints Smooth**

```python
# Bad: Discontinuous constraint
def constraint(x):
    if x[0] > 0.5:
        return 1.0
    else:
        return -1.0

# Good: Smooth constraint
def constraint(x):
    return x[0] - 0.5
```

**3. Avoid Redundant Constraints**

```python
# Redundant
constraints = [
    {'type': 'ineq', 'fun': lambda x: x[0] - 1},
    {'type': 'ineq', 'fun': lambda x: x[0] - 0.5},  # Redundant if first satisfied
]

# Better: Only keep the tightest constraint
constraints = [
    {'type': 'ineq', 'fun': lambda x: x[0] - 1},
]
```

**4. Normalize Constraints**

```python
# Bad: Different magnitudes
def constraint_power(x):
    return 50000 - calculate_power(x)  # ~50kW

def constraint_efficiency(x):
    return calculate_efficiency(x) - 0.75  # ~0.75

# Good: Similar magnitudes
def constraint_power(x):
    return (50000 - calculate_power(x)) / 50000  # Normalized

def constraint_efficiency(x):
    return (calculate_efficiency(x) - 0.75) / 0.25  # Normalized
```

## Common Pump Design Constraints

### 1. Performance Constraints

```python
# Minimum flow rate
def constraint_flow_min(x):
    Q = calculate_flow(x)
    return Q - Q_min

# Head requirement
def constraint_head(x):
    H = calculate_head(x)
    return H - H_required

# Maximum power
def constraint_power_max(x):
    P = calculate_power(x)
    return P_max - P

# Efficiency requirement
def constraint_efficiency(x):
    eta = calculate_efficiency(x)
    return eta - eta_min
```

### 2. Geometric Constraints

```python
# Diameter limits
bounds_D = (0.1, 0.5)  # meters

# Width to diameter ratio
def constraint_width_ratio(x):
    D, b = x[0], x[1]
    return b / D - 0.05  # b/D ≥ 0.05

# Blade angle limits
bounds_beta = (15, 50)  # degrees
```

### 3. Cavitation (NPSH)

```python
def constraint_npsh(x):
    NPSH_required = calculate_npsh_required(x)
    NPSH_available = get_npsh_available()
    return NPSH_available - NPSH_required - safety_margin
```

### 4. Stress Limits

```python
def constraint_stress(x):
    D, omega = x[0], x[1]
    rho_material = 7850  # kg/m³ (steel)
    sigma_yield = 250e6  # Pa

    # Centrifugal stress at rim
    sigma = rho_material * (omega * D/2)**2

    # Factor of safety = 2
    return sigma_yield / 2 - sigma
```

### 5. Manufacturing Constraints

```python
# Minimum blade thickness
def constraint_thickness(x):
    t = calculate_blade_thickness(x)
    return t - t_min  # t ≥ 2mm typically

# Maximum diameter (lathe capacity)
def constraint_diameter_max(x):
    D = x[0]
    return D_max_manufacturing - D
```

## Optimization Workflow Template

```python
import numpy as np
from scipy.optimize import minimize, differential_evolution

# ==============================================================================
# 1. DEFINE PROBLEM
# ==============================================================================

# Design variables
# x = [D (m), omega (RPM), b (m), beta (deg)]

# Constants
rho = 1000  # kg/m³
g = 9.81    # m/s²

# Requirements
Q_required = 0.15  # m³/s
H_required = 50    # m

# ==============================================================================
# 2. DEFINE OBJECTIVE FUNCTION
# ==============================================================================

def objective(x):
    """Objective to minimize (e.g., cost, negative efficiency)"""
    D, omega_rpm, b, beta = x

    # Calculate performance
    eta = calculate_efficiency(x)
    cost = calculate_cost(x)

    # Choose objective
    return -eta  # Maximize efficiency
    # return cost  # Minimize cost
    # return -eta + 0.1*cost  # Multi-objective

# ==============================================================================
# 3. DEFINE CONSTRAINTS
# ==============================================================================

def constraint_flow(x):
    Q = calculate_flow(x)
    return Q - Q_required * 0.95  # Allow 5% tolerance

def constraint_head(x):
    H = calculate_head(x)
    return H - H_required * 0.95

constraints = [
    {'type': 'ineq', 'fun': constraint_flow},
    {'type': 'ineq', 'fun': constraint_head},
]

# ==============================================================================
# 4. DEFINE BOUNDS
# ==============================================================================

bounds = [
    (0.15, 0.40),   # D (m)
    (1000, 3000),   # omega (RPM)
    (0.02, 0.10),   # b (m)
    (20, 45),       # beta (deg)
]

# ==============================================================================
# 5. INITIAL GUESS
# ==============================================================================

x0 = [0.25, 1800, 0.05, 28]  # Engineering estimate

# ==============================================================================
# 6. SOLVE
# ==============================================================================

# Option A: Local optimization (fast)
result = minimize(objective, x0=x0, method='SLSQP',
                 bounds=bounds, constraints=constraints,
                 options={'ftol': 1e-8, 'maxiter': 200})

# Option B: Global optimization (robust but slow)
# result = differential_evolution(objective, bounds,
#                                 constraints=constraints,
#                                 seed=42, maxiter=100)

# ==============================================================================
# 7. VERIFY AND REPORT
# ==============================================================================

if result.success:
    print("Optimization successful!")
    print(f"Optimal design: {result.x}")
    print(f"Objective value: {result.fun}")

    # Check constraints
    for i, c in enumerate(constraints):
        val = c['fun'](result.x)
        print(f"Constraint {i}: {val:.6f} (should be ≥ 0)")
else:
    print(f"Optimization failed: {result.message}")

# ==============================================================================
# 8. SENSITIVITY ANALYSIS (OPTIONAL)
# ==============================================================================

# Perturb each variable by ±1% and re-optimize
# This shows which variables most affect the objective
```

## Performance Tips

### 1. Vectorization

```python
# Slow: Python loop
def objective(x):
    result = 0
    for i in range(len(x)):
        result += x[i]**2
    return result

# Fast: NumPy vectorization
def objective(x):
    return np.sum(x**2)
```

### 2. Caching Expensive Computations

```python
from functools import lru_cache

@lru_cache(maxsize=1000)
def expensive_calculation(D, omega):
    # CFD simulation, etc.
    time.sleep(0.1)
    return result

def objective(x):
    D, omega = x
    # Cached - won't recalculate if called with same arguments
    return expensive_calculation(round(D, 6), round(omega, 2))
```

### 3. Parallel Function Evaluations

```python
# For population-based methods like differential_evolution
result = differential_evolution(
    objective, bounds,
    workers=-1,  # Use all CPU cores
    updating='deferred'
)
```

## Resources and Documentation

### Official Documentation
- **SciPy Optimize Tutorial**: https://docs.scipy.org/doc/scipy/tutorial/optimize.html
- **API Reference**: https://docs.scipy.org/doc/scipy/reference/optimize.html
- **SciPy Lecture Notes**: https://scipy-lectures.org/advanced/mathematical_optimization/

### Algorithm Details
- **SLSQP**: Kraft, D. (1988). A software package for sequential quadratic programming
- **L-BFGS-B**: Zhu et al. (1997). Algorithm 778: L-BFGS-B
- **Trust Region**: Conn et al. (2000). Trust-Region Methods
- **Differential Evolution**: Storn & Price (1997). Differential Evolution

### Books
- Nocedal & Wright: "Numerical Optimization" (2006) - Comprehensive optimization theory
- Boyd & Vandenberghe: "Convex Optimization" (2004) - Convex problems
- Rao: "Engineering Optimization" (2009) - Engineering applications

### Pump-Specific Resources
- Gülich: "Centrifugal Pumps" - Performance equations and design constraints
- Stepanoff: "Centrifugal and Axial Flow Pumps" - Classical design methods
- Karassik et al.: "Pump Handbook" - Comprehensive pump engineering

### Online Tools
- **SciPy Cookbook**: Practical examples and recipes
- **Stack Overflow**: scipy + optimization tag
- **GitHub**: Search for scipy optimization pump design examples

## Quick Reference Card

```python
# Unconstrained
minimize(f, x0, method='BFGS')

# With bounds only
minimize(f, x0, method='L-BFGS-B', bounds=bounds)

# With constraints
minimize(f, x0, method='SLSQP', bounds=bounds, constraints=cons)

# Global optimization
differential_evolution(f, bounds)

# Curve fitting
curve_fit(model, x_data, y_data)

# Least squares
least_squares(residuals, x0, args=(x_data, y_data))

# With gradient
minimize(f, x0, jac=grad, method='BFGS')

# Multiple starting points
results = [minimize(f, x0=guess) for guess in guesses]
best = min(results, key=lambda r: r.fun)
```

---

*This reference guide is designed for pump engineers using scipy.optimize for design optimization. For additional help, refer to the SKILL.md overview and examples.py implementations.*
