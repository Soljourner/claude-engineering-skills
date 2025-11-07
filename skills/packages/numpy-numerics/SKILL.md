---
name: numpy-numerics
description: "Numerical array operations for velocity fields, pump curves, and matrix calculations"
category: packages
domain: general
complexity: basic
dependencies:
  - numpy
---

# NumPy Numerics Skill

## Overview

NumPy is the fundamental package for numerical computing in Python. It provides:

- Powerful N-dimensional array objects
- Broadcasting functions for element-wise operations
- Linear algebra, Fourier transform, and random number capabilities
- Tools for integrating C/C++ and Fortran code
- High-performance operations on large arrays

For engineering applications, NumPy is essential for:
- Processing velocity field data from CFD simulations
- Handling pump curve data and interpolation
- Solving linear systems (pipe networks, structural analysis)
- Statistical analysis of experimental data
- Numerical integration and differentiation
- Matrix operations for finite element methods

## Installation

```bash
pip install numpy
```

For scientific computing with additional tools:

```bash
pip install numpy scipy matplotlib
```

## Core Concepts

### Arrays vs Lists

NumPy arrays are more efficient than Python lists for numerical operations:

```python
import numpy as np

# Python list - slow for numerical operations
python_list = [1, 2, 3, 4, 5]
result = [x * 2 for x in python_list]

# NumPy array - fast vectorized operations
numpy_array = np.array([1, 2, 3, 4, 5])
result = numpy_array * 2  # 10-100x faster
```

**Key Differences**:
- Arrays have fixed size and homogeneous type
- Arrays support vectorized operations
- Arrays use contiguous memory (cache-friendly)
- Arrays integrate with C/Fortran libraries

### Array Creation

```python
import numpy as np

# From Python list
a = np.array([1, 2, 3, 4, 5])

# Create ranges
x = np.arange(0, 10, 0.5)  # Start, stop, step
x = np.linspace(0, 10, 100)  # Start, stop, num_points

# Initialize with zeros/ones
data = np.zeros((10, 10))  # 10x10 array of zeros
grid = np.ones((5, 5))     # 5x5 array of ones

# Identity matrix
I = np.eye(3)  # 3x3 identity matrix

# Create meshgrid for 2D functions
x = np.linspace(0, 1, 50)
y = np.linspace(0, 1, 50)
X, Y = np.meshgrid(x, y)
```

### Array Indexing and Slicing

```python
import numpy as np

# 1D array indexing
a = np.array([10, 20, 30, 40, 50])
print(a[0])      # 10 (first element)
print(a[-1])     # 50 (last element)
print(a[1:4])    # [20, 30, 40] (slice)

# 2D array indexing
data = np.array([[1, 2, 3],
                 [4, 5, 6],
                 [7, 8, 9]])
print(data[0, 0])    # 1 (row 0, col 0)
print(data[:, 0])    # [1, 4, 7] (first column)
print(data[0, :])    # [1, 2, 3] (first row)
print(data[1:, 1:])  # [[5, 6], [8, 9]] (subarray)

# Boolean indexing
a = np.array([1, 2, 3, 4, 5])
mask = a > 2
print(a[mask])  # [3, 4, 5] (elements > 2)
```

## Engineering Applications

### 1. Velocity Field Calculations

Computing velocity fields from flow data:

```python
import numpy as np

# Create a 2D velocity field mesh (e.g., from CFD output)
x = np.linspace(0, 1, 50)  # x-coordinates (m)
y = np.linspace(0, 1, 50)  # y-coordinates (m)
X, Y = np.meshgrid(x, y)

# Velocity components (example: potential flow around cylinder)
r = np.sqrt((X - 0.5)**2 + (Y - 0.5)**2)
theta = np.arctan2(Y - 0.5, X - 0.5)
R = 0.1  # Cylinder radius

# Velocity field
u = 1.0 * (1 - R**2 / r**2) * np.cos(theta)
v = -1.0 * (1 - R**2 / r**2) * np.sin(theta)

# Compute velocity magnitude
velocity_magnitude = np.sqrt(u**2 + v**2)

# Find maximum velocity
max_velocity = np.max(velocity_magnitude)
max_location = np.unravel_index(
    np.argmax(velocity_magnitude),
    velocity_magnitude.shape
)

print(f"Maximum velocity: {max_velocity:.2f} m/s")
print(f"At location: ({X[max_location]:.2f}, {Y[max_location]:.2f})")

# Calculate flow statistics
mean_velocity = np.mean(velocity_magnitude)
std_velocity = np.std(velocity_magnitude)
print(f"Mean velocity: {mean_velocity:.2f} ± {std_velocity:.2f} m/s")
```

### 2. Pump Curve Data Handling

Working with pump performance data:

```python
import numpy as np

# Experimental pump data (flow rate, head, efficiency)
Q_data = np.array([0, 20, 40, 60, 80, 100, 120])  # m³/h
H_data = np.array([85, 84, 82, 78, 72, 64, 52])    # m
eta_data = np.array([0, 45, 68, 78, 80, 75, 65])   # %

# Fit pump curve to polynomial: H = a + b*Q + c*Q²
coeffs = np.polyfit(Q_data, H_data, 2)  # 2nd degree polynomial
print(f"Pump curve: H = {coeffs[0]:.6f}*Q² + {coeffs[1]:.4f}*Q + {coeffs[2]:.2f}")

# Create smooth curve for plotting
Q_smooth = np.linspace(0, 120, 100)
H_smooth = np.polyval(coeffs, Q_smooth)

# Find best efficiency point (BEP)
bep_index = np.argmax(eta_data)
Q_bep = Q_data[bep_index]
H_bep = H_data[bep_index]
eta_bep = eta_data[bep_index]

print(f"\nBest Efficiency Point:")
print(f"  Flow rate: {Q_bep:.1f} m³/h")
print(f"  Head: {H_bep:.1f} m")
print(f"  Efficiency: {eta_bep:.1f}%")

# Interpolate head at specific flow rate
Q_target = 65  # m³/h
H_target = np.interp(Q_target, Q_data, H_data)
print(f"\nAt Q = {Q_target} m³/h: H = {H_target:.1f} m")
```

### 3. Matrix Operations for Linear Systems

Solving pipe network equations:

```python
import numpy as np

# Pipe network with 4 nodes, 5 pipes
# Conservation of mass at each node: A * Q = b
# A: connectivity matrix
# Q: flow rates (unknown)
# b: external flows (boundary conditions)

# Example: Simple pipe network
# Node 1: Q1 + Q2 = 10 (inlet 10 L/s)
# Node 2: Q3 - Q1 = 0
# Node 3: Q4 - Q2 - Q3 = 0
# Node 4: -Q4 = -10 (outlet 10 L/s)

A = np.array([
    [1,  1,  0,  0],  # Node 1
    [-1, 0,  1,  0],  # Node 2
    [0, -1, -1,  1],  # Node 3
    [0,  0,  0, -1]   # Node 4
])

b = np.array([10, 0, 0, -10])  # L/s

# Solve linear system
Q = np.linalg.solve(A, b)

print("Flow rates in pipes:")
for i, q in enumerate(Q, 1):
    print(f"  Pipe {i}: {q:.2f} L/s")

# Verify solution
residual = np.dot(A, Q) - b
print(f"\nResidual (should be ~0): {np.max(np.abs(residual)):.2e}")

# Condition number (stability indicator)
cond = np.linalg.cond(A)
print(f"Condition number: {cond:.2f}")
if cond > 1000:
    print("Warning: Matrix is ill-conditioned")
```

### 4. Interpolation of Experimental Data

Interpolating pressure drop data:

```python
import numpy as np

# Experimental data: flow rate vs pressure drop
Q_exp = np.array([0, 10, 20, 30, 40, 50, 60])      # m³/h
dP_exp = np.array([0, 5, 18, 38, 65, 98, 138])     # kPa

# Linear interpolation (fast, simple)
Q_query = 25  # m³/h
dP_linear = np.interp(Q_query, Q_exp, dP_exp)
print(f"Linear interpolation at Q={Q_query}: dP={dP_linear:.1f} kPa")

# Polynomial interpolation for smooth curve
# Fit to quadratic: dP = k*Q²
coeffs = np.polyfit(Q_exp, dP_exp, 2)
dP_poly = np.polyval(coeffs, Q_query)
print(f"Polynomial fit at Q={Q_query}: dP={dP_poly:.1f} kPa")

# For more complex interpolation, use scipy.interpolate
# This is shown as reference (requires scipy)
from scipy import interpolate
f_spline = interpolate.interp1d(Q_exp, dP_exp, kind='cubic')
Q_dense = np.linspace(0, 60, 100)
dP_spline = f_spline(Q_dense)
```

### 5. Numerical Integration

Calculating flow rate from velocity profile:

```python
import numpy as np

# Velocity profile in circular pipe (laminar flow)
# u(r) = u_max * (1 - (r/R)²)
R = 0.05  # m, pipe radius
u_max = 2.0  # m/s, centerline velocity

# Create radial points
r = np.linspace(0, R, 100)

# Velocity distribution
u = u_max * (1 - (r/R)**2)

# Calculate flow rate using trapezoidal rule
# Q = ∫∫ u dA = ∫₀ᴿ u(r) * 2πr dr
integrand = u * 2 * np.pi * r
Q_trapz = np.trapz(integrand, r)

# Analytical solution: Q = (π*R²*u_max)/2
Q_analytical = np.pi * R**2 * u_max / 2

print(f"Flow rate (numerical): {Q_trapz:.6f} m³/s")
print(f"Flow rate (analytical): {Q_analytical:.6f} m³/s")
print(f"Error: {abs(Q_trapz - Q_analytical)/Q_analytical * 100:.2f}%")

# Simpson's rule (more accurate for smooth functions)
# Requires scipy
from scipy import integrate
Q_simps = integrate.simpson(integrand, r)
print(f"Flow rate (Simpson): {Q_simps:.6f} m³/s")
```

### 6. Numerical Differentiation

Computing pressure gradient from pressure field:

```python
import numpy as np

# Pressure distribution along pipe
x = np.linspace(0, 10, 100)  # m, position
P = 500000 - 2000 * x  # Pa, linear pressure drop

# Numerical gradient (central difference)
dP_dx = np.gradient(P, x)

print(f"Pressure gradient: {dP_dx[0]:.1f} Pa/m")
print(f"Expected: -2000 Pa/m")

# For 2D pressure field
x = np.linspace(0, 1, 50)
y = np.linspace(0, 1, 50)
X, Y = np.meshgrid(x, y)

# Example pressure field
P_field = 100000 - 5000 * X - 3000 * Y

# Compute gradients
dP_dx, dP_dy = np.gradient(P_field, x, y)

print(f"\nPressure gradient field:")
print(f"  ∂P/∂x at (0.5, 0.5): {dP_dx[25, 25]:.1f} Pa/m")
print(f"  ∂P/∂y at (0.5, 0.5): {dP_dy[25, 25]:.1f} Pa/m")
```

### 7. Statistical Analysis of Test Data

Analyzing pump test measurements:

```python
import numpy as np

# Multiple test runs for pump efficiency
test_runs = np.array([
    [78.5, 79.2, 78.8, 79.0, 78.6],  # Run 1-5
    [77.8, 78.5, 78.0, 78.3, 77.9],  # Repeat tests
    [79.1, 78.9, 79.3, 78.7, 79.0]
])

# Statistical measures
mean_efficiency = np.mean(test_runs)
std_efficiency = np.std(test_runs, ddof=1)  # Sample std dev
sem = std_efficiency / np.sqrt(test_runs.size)  # Standard error

print(f"Mean efficiency: {mean_efficiency:.2f}%")
print(f"Standard deviation: {std_efficiency:.2f}%")
print(f"Standard error: {sem:.2f}%")
print(f"95% confidence interval: {mean_efficiency:.2f} ± {1.96*sem:.2f}%")

# Detect outliers (beyond 3 standard deviations)
z_scores = np.abs((test_runs - mean_efficiency) / std_efficiency)
outliers = test_runs[z_scores > 3]
if len(outliers) > 0:
    print(f"\nOutliers detected: {outliers}")
else:
    print("\nNo outliers detected")

# Range and percentiles
print(f"\nData range: {np.min(test_runs):.2f}% to {np.max(test_runs):.2f}%")
print(f"Median: {np.median(test_runs):.2f}%")
print(f"25th percentile: {np.percentile(test_runs, 25):.2f}%")
print(f"75th percentile: {np.percentile(test_runs, 75):.2f}%")
```

## Performance Tips for Large Arrays

### 1. Use Vectorization Instead of Loops

```python
import numpy as np
import time

n = 1000000

# Slow: Python loop
start = time.time()
result = []
for i in range(n):
    result.append(i ** 2)
loop_time = time.time() - start

# Fast: NumPy vectorization
start = time.time()
result = np.arange(n) ** 2
vectorized_time = time.time() - start

print(f"Loop: {loop_time:.3f} s")
print(f"Vectorized: {vectorized_time:.3f} s")
print(f"Speedup: {loop_time/vectorized_time:.0f}x")
```

### 2. Preallocate Arrays

```python
import numpy as np

# Bad: Growing array in loop
result = np.array([])
for i in range(1000):
    result = np.append(result, i)  # Creates new array each time

# Good: Preallocate
result = np.zeros(1000)
for i in range(1000):
    result[i] = i  # In-place assignment

# Best: Use vectorization
result = np.arange(1000)
```

### 3. Use In-Place Operations

```python
import numpy as np

a = np.random.rand(1000000)

# Creates new array
b = a * 2
b = a + 5

# In-place operations (faster, less memory)
a *= 2
a += 5
```

### 4. Choose Appropriate Data Types

```python
import numpy as np

# Default: float64 (8 bytes per element)
a = np.array([1.0, 2.0, 3.0])
print(f"Size: {a.nbytes} bytes")

# Use float32 if precision allows (4 bytes per element)
a = np.array([1.0, 2.0, 3.0], dtype=np.float32)
print(f"Size: {a.nbytes} bytes (50% reduction)")

# For integers, use smallest appropriate type
counts = np.array([1, 2, 3], dtype=np.uint8)  # 0-255
```

### 5. Use Memory Views for Large Arrays

```python
import numpy as np

# Original array
data = np.arange(1000000)

# View (no copy, fast)
view = data[::2]  # Every other element
view[0] = 999  # Modifies original array

# Copy (slower, independent)
copy = data[::2].copy()
copy[0] = 888  # Does not modify original
```

### 6. Leverage Broadcasting

```python
import numpy as np

# Broadcasting: automatic array shape matching
a = np.array([[1, 2, 3],
              [4, 5, 6],
              [7, 8, 9]])

b = np.array([10, 20, 30])

# Adds b to each row of a (no loop needed)
result = a + b
# [[11, 22, 33],
#  [14, 25, 36],
#  [17, 28, 39]]
```

## Common Pitfalls

### 1. Array vs Scalar Division

```python
import numpy as np

# Integer division
a = np.array([1, 2, 3])
result = a / 2  # [0.5, 1.0, 1.5] (float result)

# Floor division
result = a // 2  # [0, 1, 1] (integer result)
```

### 2. Copying vs Viewing

```python
import numpy as np

a = np.array([1, 2, 3])

# View (references same data)
b = a
b[0] = 999  # Modifies a!

# Copy (independent data)
c = a.copy()
c[0] = 888  # Does not modify a
```

### 3. Integer Overflow

```python
import numpy as np

# int8 can only hold -128 to 127
a = np.array([100, 50], dtype=np.int8)
result = a + a  # [-56, 100] (overflow!)

# Use appropriate dtype
a = np.array([100, 50], dtype=np.int32)
result = a + a  # [200, 100] (correct)
```

## Best Practices

1. **Vectorize operations** - Avoid Python loops when possible
2. **Use appropriate dtypes** - Balance precision and memory
3. **Preallocate arrays** - Don't grow arrays in loops
4. **Leverage broadcasting** - Implicit dimension matching
5. **Check array shapes** - Use `.shape` to debug dimension errors
6. **Use built-in functions** - NumPy functions are optimized
7. **Profile code** - Identify bottlenecks before optimizing
8. **Document units** - Add comments for physical quantities

## Integration with Engineering Workflows

NumPy integrates seamlessly with:

- **Pandas**: DataFrames built on NumPy arrays
- **Matplotlib**: Plotting arrays directly
- **SciPy**: Advanced scientific computing (interpolation, optimization)
- **Fluids**: Engineering calculations on array data
- **CoolProp**: Thermodynamic property tables

Example workflow:

```python
import numpy as np
import matplotlib.pyplot as plt
from fluids import Reynolds, friction_factor

# Generate velocity range
V = np.linspace(0.1, 5.0, 50)  # m/s

# Calculate Reynolds numbers
D = 0.1  # m
rho = 1000  # kg/m³
mu = 0.001  # Pa·s

Re = np.array([Reynolds(v, D, rho, mu) for v in V])

# Calculate friction factors
eD = 0.0001
f = np.array([friction_factor(re, eD) for re in Re])

# Plot results
plt.plot(Re, f, 'b-', linewidth=2)
plt.xlabel('Reynolds Number')
plt.ylabel('Friction Factor')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.title('Friction Factor vs Reynolds Number')
plt.show()
```

## References

- NumPy Documentation: https://numpy.org/doc/
- NumPy User Guide: https://numpy.org/doc/stable/user/
- SciPy Lecture Notes: https://scipy-lectures.org/
- "Python for Data Analysis" by Wes McKinney
- "Numerical Python" by Robert Johansson

## Further Learning

- **NumPy Quickstart**: https://numpy.org/doc/stable/user/quickstart.html
- **NumPy for MATLAB Users**: https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
- **Array Programming**: https://numpy.org/doc/stable/user/basics.html
- **Linear Algebra**: https://numpy.org/doc/stable/reference/routines.linalg.html
