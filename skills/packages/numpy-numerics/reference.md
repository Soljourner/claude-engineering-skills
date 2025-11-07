# NumPy Numerics Reference

## Quick Reference Guide

This reference provides essential NumPy functions for engineering applications including array operations, linear algebra, and numerical methods.

---

## Array Creation Functions

### np.array()

```python
np.array(object, dtype=None)
```

**Description**: Create array from Python list or tuple

**Parameters**:
- `object`: Array-like structure (list, tuple, etc.)
- `dtype`: Data type (optional, e.g., `np.float64`, `np.int32`)

**Returns**: NumPy array

**Example**:
```python
a = np.array([1, 2, 3, 4, 5])
b = np.array([[1, 2], [3, 4]], dtype=np.float32)
```

---

### np.arange()

```python
np.arange(start, stop, step)
```

**Description**: Create array with evenly spaced values

**Parameters**:
- `start`: Starting value (inclusive)
- `stop`: End value (exclusive)
- `step`: Spacing between values

**Returns**: 1D array

**Example**:
```python
x = np.arange(0, 10, 0.5)  # [0, 0.5, 1.0, ..., 9.5]
```

**Note**: Similar to Python's `range()` but returns array and accepts floats

---

### np.linspace()

```python
np.linspace(start, stop, num=50)
```

**Description**: Create array with specified number of evenly spaced points

**Parameters**:
- `start`: Starting value (inclusive)
- `stop`: End value (inclusive)
- `num`: Number of points

**Returns**: 1D array

**Example**:
```python
x = np.linspace(0, 1, 100)  # 100 points from 0 to 1
```

**Engineering Use**: Preferred over `arange()` for creating coordinate arrays because you specify number of points (not step size)

---

### np.zeros()

```python
np.zeros(shape, dtype=float)
```

**Description**: Create array filled with zeros

**Parameters**:
- `shape`: Shape tuple (e.g., `(3, 4)` for 3×4 array)
- `dtype`: Data type

**Returns**: Array of zeros

**Example**:
```python
data = np.zeros((10, 10))  # 10×10 array of zeros
```

---

### np.ones()

```python
np.ones(shape, dtype=float)
```

**Description**: Create array filled with ones

**Example**:
```python
grid = np.ones((5, 5))  # 5×5 array of ones
```

---

### np.eye()

```python
np.eye(N, M=None, k=0)
```

**Description**: Create identity matrix

**Parameters**:
- `N`: Number of rows
- `M`: Number of columns (default: N)
- `k`: Diagonal offset (0 = main diagonal)

**Returns**: 2D identity matrix

**Example**:
```python
I = np.eye(3)  # 3×3 identity matrix
```

**Engineering Use**: Initial condition for matrix operations, identity transformations

---

### np.meshgrid()

```python
np.meshgrid(x, y)
```

**Description**: Create coordinate matrices from coordinate vectors

**Parameters**:
- `x`: 1D array of x-coordinates
- `y`: 1D array of y-coordinates

**Returns**: Two 2D arrays (X, Y)

**Example**:
```python
x = np.linspace(0, 1, 50)
y = np.linspace(0, 1, 50)
X, Y = np.meshgrid(x, y)
# Now X[i,j] and Y[i,j] are coordinates of point (i,j)
```

**Engineering Use**: Creating 2D grids for velocity fields, temperature distributions, etc.

---

## Array Manipulation

### Indexing and Slicing

```python
a[i]           # Single element
a[i:j]         # Slice from i to j-1
a[i:j:k]       # Slice with step k
a[:, i]        # Column i
a[i, :]        # Row i
a[mask]        # Boolean indexing
```

**Example**:
```python
a = np.array([10, 20, 30, 40, 50])
print(a[0])      # 10
print(a[1:4])    # [20, 30, 40]
print(a[a > 25]) # [30, 40, 50]
```

---

### np.reshape()

```python
np.reshape(a, newshape)
```

**Description**: Change array shape without changing data

**Example**:
```python
a = np.arange(12)  # [0, 1, 2, ..., 11]
b = a.reshape(3, 4)  # 3×4 array
```

---

### np.concatenate()

```python
np.concatenate((a1, a2, ...), axis=0)
```

**Description**: Join arrays along existing axis

**Example**:
```python
a = np.array([1, 2, 3])
b = np.array([4, 5, 6])
c = np.concatenate((a, b))  # [1, 2, 3, 4, 5, 6]
```

---

## Mathematical Operations

### Element-wise Operations

```python
a + b          # Addition
a - b          # Subtraction
a * b          # Multiplication
a / b          # Division
a ** b         # Power
np.sqrt(a)     # Square root
np.exp(a)      # Exponential
np.log(a)      # Natural logarithm
np.sin(a)      # Sine
np.cos(a)      # Cosine
```

**Broadcasting**: Arrays with different shapes are automatically expanded:
```python
a = np.array([[1, 2, 3],
              [4, 5, 6]])
b = np.array([10, 20, 30])
c = a + b  # b is broadcasted to each row
# [[11, 22, 33],
#  [14, 25, 36]]
```

---

## Statistical Functions

### np.mean()

```python
np.mean(a, axis=None)
```

**Description**: Calculate arithmetic mean

**Parameters**:
- `a`: Input array
- `axis`: Axis along which to compute (None = all elements)

**Returns**: Mean value(s)

**Example**:
```python
data = np.array([1, 2, 3, 4, 5])
print(np.mean(data))  # 3.0
```

---

### np.std()

```python
np.std(a, axis=None, ddof=0)
```

**Description**: Calculate standard deviation

**Parameters**:
- `a`: Input array
- `axis`: Axis along which to compute
- `ddof`: Delta degrees of freedom (use 1 for sample std dev)

**Returns**: Standard deviation

**Example**:
```python
data = np.array([1, 2, 3, 4, 5])
print(np.std(data, ddof=1))  # Sample standard deviation
```

---

### np.median()

```python
np.median(a, axis=None)
```

**Description**: Calculate median value

---

### np.percentile()

```python
np.percentile(a, q, axis=None)
```

**Description**: Calculate q-th percentile

**Parameters**:
- `a`: Input array
- `q`: Percentile (0-100)

**Example**:
```python
data = np.random.rand(1000)
p25 = np.percentile(data, 25)  # 25th percentile
p75 = np.percentile(data, 75)  # 75th percentile
```

---

### np.min(), np.max()

```python
np.min(a, axis=None)
np.max(a, axis=None)
```

**Description**: Find minimum/maximum values

**Example**:
```python
data = np.array([3, 1, 4, 1, 5, 9, 2, 6])
print(np.min(data))  # 1
print(np.max(data))  # 9
```

---

### np.argmin(), np.argmax()

```python
np.argmin(a, axis=None)
np.argmax(a, axis=None)
```

**Description**: Find indices of minimum/maximum values

**Example**:
```python
data = np.array([3, 1, 4, 1, 5, 9, 2, 6])
idx_max = np.argmax(data)  # 5
```

---

## Linear Algebra (np.linalg)

### np.linalg.solve()

```python
np.linalg.solve(A, b)
```

**Description**: Solve linear system Ax = b

**Parameters**:
- `A`: Coefficient matrix (n×n)
- `b`: Right-hand side vector/matrix (n×1 or n×m)

**Returns**: Solution vector/matrix x

**Example**:
```python
# Solve 2x + 3y = 8
#       x -  y = 1
A = np.array([[2, 3],
              [1, -1]])
b = np.array([8, 1])
x = np.linalg.solve(A, b)  # [2.6, 1.6]
```

**Engineering Use**: Pipe networks, structural analysis, circuit analysis

**Requirements**: A must be square and non-singular (det(A) ≠ 0)

---

### np.linalg.det()

```python
np.linalg.det(A)
```

**Description**: Calculate determinant of matrix

**Example**:
```python
A = np.array([[1, 2], [3, 4]])
det = np.linalg.det(A)  # -2.0
```

**Engineering Use**: Check if system is solvable (det ≠ 0)

---

### np.linalg.inv()

```python
np.linalg.inv(A)
```

**Description**: Calculate inverse of matrix

**Example**:
```python
A = np.array([[1, 2], [3, 4]])
A_inv = np.linalg.inv(A)
# Verify: A @ A_inv ≈ I
```

**Warning**: Avoid using inverse for solving Ax=b. Use `np.linalg.solve()` instead (more stable and faster)

---

### np.linalg.cond()

```python
np.linalg.cond(A)
```

**Description**: Calculate condition number of matrix

**Returns**: Condition number (ratio of largest to smallest singular value)

**Interpretation**:
- cond ≈ 1: Well-conditioned (stable)
- cond > 1000: Ill-conditioned (numerical errors likely)
- cond → ∞: Singular or near-singular

**Example**:
```python
A = np.array([[1, 2], [3, 4]])
cond = np.linalg.cond(A)
print(f"Condition number: {cond:.2f}")
```

---

### np.linalg.eig()

```python
np.linalg.eig(A)
```

**Description**: Calculate eigenvalues and eigenvectors

**Returns**: (eigenvalues, eigenvectors)

**Example**:
```python
A = np.array([[1, 2], [2, 1]])
eigenvalues, eigenvectors = np.linalg.eig(A)
```

**Engineering Use**: Stability analysis, vibration modes, principal stresses

---

### np.dot() or @

```python
np.dot(A, B)
A @ B  # Matrix multiplication operator (Python 3.5+)
```

**Description**: Matrix multiplication

**Example**:
```python
A = np.array([[1, 2], [3, 4]])
B = np.array([[5, 6], [7, 8]])
C = A @ B  # Matrix product
```

---

### np.linalg.norm()

```python
np.linalg.norm(x, ord=None)
```

**Description**: Calculate vector or matrix norm

**Parameters**:
- `x`: Input array
- `ord`: Norm order (None = 2-norm, 1 = 1-norm, np.inf = ∞-norm)

**Example**:
```python
v = np.array([3, 4])
magnitude = np.linalg.norm(v)  # 5.0
```

**Engineering Use**: Velocity magnitude, force resultant, error norm

---

## Interpolation

### np.interp()

```python
np.interp(x, xp, fp)
```

**Description**: 1D linear interpolation

**Parameters**:
- `x`: Points to interpolate at
- `xp`: x-coordinates of data points (must be increasing)
- `fp`: y-coordinates of data points

**Returns**: Interpolated values

**Example**:
```python
# Pump curve data
Q = np.array([0, 20, 40, 60, 80, 100])  # m³/h
H = np.array([85, 82, 78, 72, 64, 52])   # m

# Interpolate at Q=50 m³/h
H_interp = np.interp(50, Q, H)  # ≈ 75 m
```

**Engineering Use**: Pump curves, material properties, lookup tables

---

### np.polyfit()

```python
np.polyfit(x, y, deg)
```

**Description**: Fit polynomial of degree `deg` to data

**Parameters**:
- `x`: x-coordinates of data points
- `y`: y-coordinates of data points
- `deg`: Degree of polynomial

**Returns**: Polynomial coefficients (highest degree first)

**Example**:
```python
# Fit quadratic: y = ax² + bx + c
x = np.array([0, 1, 2, 3, 4])
y = np.array([1, 3, 7, 13, 21])
coeffs = np.polyfit(x, y, 2)  # [1.0, 1.0, 1.0]
# y ≈ x² + x + 1
```

**Engineering Use**: Curve fitting for pump curves, test data

---

### np.polyval()

```python
np.polyval(p, x)
```

**Description**: Evaluate polynomial at x

**Parameters**:
- `p`: Polynomial coefficients
- `x`: Evaluation points

**Returns**: Polynomial values

**Example**:
```python
coeffs = [1, 2, 3]  # y = x² + 2x + 3
x = np.array([0, 1, 2])
y = np.polyval(coeffs, x)  # [3, 6, 11]
```

---

## Numerical Integration

### np.trapz()

```python
np.trapz(y, x=None, dx=1.0)
```

**Description**: Integrate using trapezoidal rule

**Parameters**:
- `y`: Values to integrate
- `x`: x-coordinates (optional, default: uniform spacing)
- `dx`: Spacing (if x not provided)

**Returns**: Integral approximation

**Example**:
```python
# ∫₀¹ x² dx = 1/3
x = np.linspace(0, 1, 100)
y = x**2
integral = np.trapz(y, x)  # ≈ 0.333
```

**Engineering Use**: Flow rate from velocity profile, work from force-displacement

**Accuracy**: O(h²) where h is spacing

---

### scipy.integrate.simpson()

```python
from scipy import integrate
integrate.simpson(y, x=None, dx=1.0)
```

**Description**: Integrate using Simpson's rule (requires scipy)

**Accuracy**: O(h⁴) - more accurate than trapezoidal for smooth functions

**Example**:
```python
from scipy import integrate
x = np.linspace(0, 1, 100)
y = x**2
integral = integrate.simpson(y, x)  # ≈ 0.3333
```

---

## Numerical Differentiation

### np.gradient()

```python
np.gradient(f, *varargs)
```

**Description**: Calculate gradient (discrete derivative)

**Parameters**:
- `f`: Array to differentiate
- `varargs`: Spacing (scalar) or coordinates (array)

**Returns**: Gradient array(s)

**Example**:
```python
# df/dx for f = x²
x = np.linspace(0, 10, 100)
f = x**2
df_dx = np.gradient(f, x)  # ≈ 2x
```

**Engineering Use**: Pressure gradient, velocity gradient, slope

**Method**: Central differences (2nd order accurate)

---

### np.diff()

```python
np.diff(a, n=1, axis=-1)
```

**Description**: Calculate discrete differences

**Parameters**:
- `a`: Input array
- `n`: Number of times to difference

**Returns**: Differences (array is one element shorter)

**Example**:
```python
x = np.array([1, 3, 6, 10, 15])
dx = np.diff(x)  # [2, 3, 4, 5]
```

---

## Random Number Generation

### np.random.rand()

```python
np.random.rand(d0, d1, ..., dn)
```

**Description**: Random values in [0, 1) from uniform distribution

---

### np.random.randn()

```python
np.random.randn(d0, d1, ..., dn)
```

**Description**: Random values from standard normal distribution (mean=0, std=1)

---

### np.random.normal()

```python
np.random.normal(loc=0.0, scale=1.0, size=None)
```

**Description**: Random values from normal distribution

**Parameters**:
- `loc`: Mean
- `scale`: Standard deviation
- `size`: Output shape

**Example**:
```python
# Generate test data with measurement error
true_value = 100
noise = np.random.normal(0, 2, 10)  # Mean=0, std=2
measurements = true_value + noise
```

---

### np.random.seed()

```python
np.random.seed(seed)
```

**Description**: Set random seed for reproducibility

**Example**:
```python
np.random.seed(42)  # Always generates same sequence
data = np.random.rand(100)
```

---

## Useful Utilities

### np.where()

```python
np.where(condition, x, y)
```

**Description**: Return elements from x or y depending on condition

**Example**:
```python
a = np.array([1, 2, 3, 4, 5])
result = np.where(a > 3, a, 0)  # [0, 0, 0, 4, 5]
```

---

### np.unravel_index()

```python
np.unravel_index(indices, shape)
```

**Description**: Convert flat index to multi-dimensional index

**Example**:
```python
data = np.random.rand(10, 10)
max_idx_flat = np.argmax(data)
i, j = np.unravel_index(max_idx_flat, data.shape)
# (i, j) are row and column of maximum
```

---

### np.isclose()

```python
np.isclose(a, b, rtol=1e-5, atol=1e-8)
```

**Description**: Element-wise comparison with tolerance

**Example**:
```python
a = 0.1 + 0.2
b = 0.3
print(a == b)  # False (floating point error)
print(np.isclose(a, b))  # True
```

---

### np.allclose()

```python
np.allclose(a, b, rtol=1e-5, atol=1e-8)
```

**Description**: Check if all array elements are close

**Example**:
```python
a = np.array([1.0, 2.0, 3.0])
b = np.array([1.0001, 2.0001, 3.0001])
print(np.allclose(a, b, rtol=1e-3))  # True
```

---

## Performance Tips

### Vectorization

```python
# Slow: Python loop
result = []
for x in data:
    result.append(x**2)

# Fast: NumPy vectorization (10-100x faster)
result = data**2
```

---

### In-place Operations

```python
# Creates new array
a = a + 1

# In-place (faster, less memory)
a += 1
```

---

### Data Types

```python
# Default: float64 (8 bytes)
a = np.array([1.0, 2.0, 3.0])

# Specify smaller type if sufficient
a = np.array([1.0, 2.0, 3.0], dtype=np.float32)  # 4 bytes
```

---

## Common Data Types

| Type | Description | Range | Size |
|------|-------------|-------|------|
| `np.int8` | Signed integer | -128 to 127 | 1 byte |
| `np.int16` | Signed integer | -32,768 to 32,767 | 2 bytes |
| `np.int32` | Signed integer | -2³¹ to 2³¹-1 | 4 bytes |
| `np.int64` | Signed integer | -2⁶³ to 2⁶³-1 | 8 bytes |
| `np.uint8` | Unsigned integer | 0 to 255 | 1 byte |
| `np.float32` | Single precision | ~7 digits | 4 bytes |
| `np.float64` | Double precision | ~16 digits | 8 bytes |
| `np.complex64` | Complex number | (float32, float32) | 8 bytes |
| `np.bool_` | Boolean | True/False | 1 byte |

---

## Array Properties

```python
a.shape      # Dimensions (e.g., (3, 4))
a.ndim       # Number of dimensions
a.size       # Total number of elements
a.dtype      # Data type
a.itemsize   # Size of each element (bytes)
a.nbytes     # Total bytes (size × itemsize)
```

---

## Documentation and Resources

### Official Documentation

- **NumPy Documentation**: https://numpy.org/doc/
- **User Guide**: https://numpy.org/doc/stable/user/
- **API Reference**: https://numpy.org/doc/stable/reference/

### Tutorials

- **NumPy Quickstart**: https://numpy.org/doc/stable/user/quickstart.html
- **NumPy for MATLAB Users**: https://numpy.org/doc/stable/user/numpy-for-matlab-users.html
- **SciPy Lecture Notes**: https://scipy-lectures.org/

### Books

- "Python for Data Analysis" by Wes McKinney
- "Numerical Python" by Robert Johansson
- "Elegant SciPy" by Juan Nunez-Iglesias et al.

### Related Libraries

- **SciPy**: https://scipy.org/ - Advanced scientific computing
- **Pandas**: https://pandas.pydata.org/ - Data analysis
- **Matplotlib**: https://matplotlib.org/ - Plotting
- **scikit-learn**: https://scikit-learn.org/ - Machine learning

---

## Common Gotchas

### 1. View vs Copy

```python
a = np.array([1, 2, 3])
b = a        # View - modifies a
c = a.copy() # Copy - independent
```

### 2. Integer Division

```python
a = np.array([1, 2, 3])
b = a / 2    # Float result: [0.5, 1.0, 1.5]
c = a // 2   # Integer result: [0, 1, 1]
```

### 3. Broadcasting Shapes

```python
# Compatible shapes for broadcasting:
(3, 4) and (4,)    # ✓ OK
(3, 4) and (3, 1)  # ✓ OK
(3, 4) and (3,)    # ✗ Error
```

### 4. Axis Parameter

```python
a = np.array([[1, 2, 3],
              [4, 5, 6]])

np.sum(a, axis=0)  # [5, 7, 9] - sum columns
np.sum(a, axis=1)  # [6, 15] - sum rows
```

---

## Engineering Workflow Integration

NumPy works seamlessly with:

```python
# Read data
import pandas as pd
df = pd.read_csv('data.csv')
data = df['pressure'].values  # Convert to NumPy array

# Calculate
mean_pressure = np.mean(data)
std_pressure = np.std(data, ddof=1)

# Plot
import matplotlib.pyplot as plt
plt.plot(data)
plt.axhline(mean_pressure, color='r', linestyle='--')
plt.show()
```

---

## Version Information

Check NumPy version:
```python
import numpy as np
print(np.__version__)
```

This reference is based on NumPy 1.24+

---

*Last Updated: 2025-11-07*
