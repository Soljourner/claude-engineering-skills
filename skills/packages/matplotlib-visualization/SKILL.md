---
name: matplotlib-visualization
description: "Create pump performance curves, velocity contours, and engineering plots"
category: packages
domain: general
complexity: basic
dependencies:
  - matplotlib
  - numpy
---

# Matplotlib Visualization for Engineering

Create professional engineering plots including pump performance curves, CFD results visualization, velocity contours, and publication-quality figures.

## Installation and Setup

```bash
pip install matplotlib numpy scipy
```

Basic import pattern:

```python
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Set style for professional figures
plt.style.use('seaborn-v0_8-darkgrid')  # or 'classic', 'bmh', 'ggplot'
```

## Common Plot Types for Engineering

### Pump Performance Curves (H-Q, P-Q, η-Q)

Standard pump characteristic curves showing head, power, and efficiency vs flow rate:

```python
import matplotlib.pyplot as plt
import numpy as np

# Sample pump data
Q = np.linspace(0, 120, 50)  # Flow rate [m³/h]
H = 80 - 0.005 * Q**2  # Head [m]
P = 5 + 0.08 * Q + 0.0005 * Q**2  # Power [kW]
eta = (Q * H) / (367 * P) * 100  # Efficiency [%]

# Create figure with three subplots
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

# Head-Flow curve
ax1.plot(Q, H, 'b-', linewidth=2, label='Head curve')
ax1.set_ylabel('Head [m]', fontsize=12)
ax1.set_title('Pump Performance Curves', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Power-Flow curve
ax2.plot(Q, P, 'r-', linewidth=2, label='Power curve')
ax2.set_ylabel('Power [kW]', fontsize=12)
ax2.grid(True, alpha=0.3)
ax2.legend()

# Efficiency-Flow curve
ax3.plot(Q, eta, 'g-', linewidth=2, label='Efficiency curve')
ax3.set_xlabel('Flow Rate [m³/h]', fontsize=12)
ax3.set_ylabel('Efficiency [%]', fontsize=12)
ax3.grid(True, alpha=0.3)
ax3.legend()

plt.tight_layout()
plt.savefig('pump_curves.png', dpi=300, bbox_inches='tight')
plt.show()
```

### System Curves

Plot pump curve with system curve overlay to find operating point:

```python
# Pump curve
Q = np.linspace(0, 120, 50)
H_pump = 80 - 0.005 * Q**2

# System curve: H = H_static + k*Q²
H_static = 20  # Static head [m]
k = 0.003  # System resistance coefficient
H_system = H_static + k * Q**2

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(Q, H_pump, 'b-', linewidth=2, label='Pump curve')
ax.plot(Q, H_system, 'r--', linewidth=2, label='System curve')

# Find and mark operating point
idx = np.argmin(np.abs(H_pump - H_system))
ax.plot(Q[idx], H_pump[idx], 'ko', markersize=10, label=f'Operating point\n({Q[idx]:.1f} m³/h, {H_pump[idx]:.1f} m)')

ax.set_xlabel('Flow Rate [m³/h]', fontsize=12)
ax.set_ylabel('Head [m]', fontsize=12)
ax.set_title('Pump and System Curves', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig('system_curve.png', dpi=300, bbox_inches='tight')
plt.show()
```

### Velocity/Pressure Contours

Visualize CFD results with contour plots:

```python
import matplotlib.pyplot as plt
import numpy as np

# Create mesh grid
x = np.linspace(0, 10, 100)
y = np.linspace(0, 5, 50)
X, Y = np.meshgrid(x, y)

# Example velocity field (flow around obstacle)
U = 2.0 - 1.5 * np.exp(-((X-3)**2 + (Y-2.5)**2) / 2)
V = 0.5 * np.sin(np.pi * X / 5) * np.exp(-((X-3)**2 + (Y-2.5)**2) / 4)
velocity = np.sqrt(U**2 + V**2)

fig, ax = plt.subplots(figsize=(12, 6))

# Filled contours
levels = np.linspace(velocity.min(), velocity.max(), 20)
contourf = ax.contourf(X, Y, velocity, levels=levels, cmap='jet')
cbar = plt.colorbar(contourf, ax=ax, label='Velocity [m/s]')

# Optional: Add contour lines
contour = ax.contour(X, Y, velocity, levels=10, colors='k', linewidths=0.5, alpha=0.4)
ax.clabel(contour, inline=True, fontsize=8, fmt='%.2f')

ax.set_xlabel('x [m]', fontsize=12)
ax.set_ylabel('y [m]', fontsize=12)
ax.set_title('Velocity Contour Plot', fontsize=14, fontweight='bold')
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig('velocity_contour.png', dpi=300, bbox_inches='tight')
plt.show()
```

### Vector Fields

Display flow direction and magnitude:

```python
import matplotlib.pyplot as plt
import numpy as np

# Create coarser grid for vectors
x = np.linspace(0, 10, 20)
y = np.linspace(0, 5, 10)
X, Y = np.meshgrid(x, y)

# Velocity components
U = 2.0 - 1.5 * np.exp(-((X-3)**2 + (Y-2.5)**2) / 2)
V = 0.5 * np.sin(np.pi * X / 5) * np.exp(-((X-3)**2 + (Y-2.5)**2) / 4)
velocity = np.sqrt(U**2 + V**2)

fig, ax = plt.subplots(figsize=(12, 6))

# Quiver plot with color by magnitude
quiver = ax.quiver(X, Y, U, V, velocity, cmap='viridis', scale=30, width=0.003)
cbar = plt.colorbar(quiver, ax=ax, label='Velocity magnitude [m/s]')

ax.set_xlabel('x [m]', fontsize=12)
ax.set_ylabel('y [m]', fontsize=12)
ax.set_title('Velocity Vector Field', fontsize=14, fontweight='bold')
ax.set_aspect('equal')

plt.tight_layout()
plt.savefig('velocity_vectors.png', dpi=300, bbox_inches='tight')
plt.show()
```

### 3D Surface Plots

Visualize pressure or temperature distributions:

```python
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Create mesh
x = np.linspace(-5, 5, 50)
y = np.linspace(-5, 5, 50)
X, Y = np.meshgrid(x, y)

# Example: Pressure distribution
Z = 100 + 10 * np.exp(-0.1 * (X**2 + Y**2)) * np.cos(0.5 * np.sqrt(X**2 + Y**2))

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Surface plot
surf = ax.plot_surface(X, Y, Z, cmap='coolwarm', linewidth=0, antialiased=True, alpha=0.9)

# Customize
ax.set_xlabel('x [m]', fontsize=11)
ax.set_ylabel('y [m]', fontsize=11)
ax.set_zlabel('Pressure [kPa]', fontsize=11)
ax.set_title('Pressure Distribution', fontsize=14, fontweight='bold', pad=20)

# Add colorbar
cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Pressure [kPa]')

# Set viewing angle
ax.view_init(elev=25, azim=45)

plt.tight_layout()
plt.savefig('pressure_3d.png', dpi=300, bbox_inches='tight')
plt.show()
```

### Subplots and Multi-Panel Figures

Complex layouts for comprehensive analysis:

```python
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(14, 10))

# Create 2x2 grid
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# Subplot 1: Line plot
ax1 = fig.add_subplot(gs[0, 0])
x = np.linspace(0, 100, 100)
y1 = 80 - 0.005 * x**2
ax1.plot(x, y1, 'b-', linewidth=2)
ax1.set_xlabel('Flow Rate [m³/h]')
ax1.set_ylabel('Head [m]')
ax1.set_title('(a) Pump Curve')
ax1.grid(True, alpha=0.3)

# Subplot 2: Contour
ax2 = fig.add_subplot(gs[0, 1])
x = np.linspace(0, 10, 50)
y = np.linspace(0, 5, 25)
X, Y = np.meshgrid(x, y)
Z = np.sin(X) * np.cos(Y)
contourf = ax2.contourf(X, Y, Z, levels=15, cmap='RdBu_r')
plt.colorbar(contourf, ax=ax2, label='Value')
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]')
ax2.set_title('(b) Contour Plot')

# Subplot 3: Bar chart
ax3 = fig.add_subplot(gs[1, 0])
categories = ['Case 1', 'Case 2', 'Case 3', 'Case 4']
values = [75, 82, 68, 91]
bars = ax3.bar(categories, values, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
ax3.set_ylabel('Efficiency [%]')
ax3.set_title('(c) Performance Comparison')
ax3.grid(True, axis='y', alpha=0.3)

# Subplot 4: Scatter with trend
ax4 = fig.add_subplot(gs[1, 1])
x_data = np.random.rand(50) * 100
y_data = 30 + 0.5 * x_data + np.random.randn(50) * 5
ax4.scatter(x_data, y_data, alpha=0.6, s=50)
# Add trend line
z = np.polyfit(x_data, y_data, 1)
p = np.poly1d(z)
ax4.plot(x_data, p(x_data), 'r--', linewidth=2, label=f'y={z[0]:.2f}x+{z[1]:.2f}')
ax4.set_xlabel('Parameter X')
ax4.set_ylabel('Parameter Y')
ax4.set_title('(d) Correlation Analysis')
ax4.legend()
ax4.grid(True, alpha=0.3)

plt.suptitle('Multi-Panel Engineering Analysis', fontsize=16, fontweight='bold', y=0.995)

plt.savefig('multipanel.png', dpi=300, bbox_inches='tight')
plt.show()
```

## Styling for Professional Figures

### Using Style Sheets

```python
# Available styles
import matplotlib.pyplot as plt
print(plt.style.available)

# Apply a style
plt.style.use('seaborn-v0_8-whitegrid')

# Or create custom style
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['grid.alpha'] = 0.3
```

### Custom Colors and Styling

```python
# Professional color schemes
colors_engineering = {
    'blue': '#0072BD',
    'red': '#D95319',
    'yellow': '#EDB120',
    'purple': '#7E2F8E',
    'green': '#77AC30',
    'cyan': '#4DBEEE',
    'maroon': '#A2142F'
}

# Apply to plots
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, y1, color=colors_engineering['blue'], linewidth=2.5, label='Series 1')
ax.plot(x, y2, color=colors_engineering['red'], linewidth=2.5, label='Series 2')
```

### LaTeX in Labels

Enable LaTeX rendering for mathematical expressions:

```python
import matplotlib.pyplot as plt
import numpy as np

# Enable LaTeX
plt.rcParams['text.usetex'] = False  # Set True if LaTeX installed
plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern font

fig, ax = plt.subplots(figsize=(10, 6))

x = np.linspace(0, 100, 100)
y = 80 - 0.005 * x**2

ax.plot(x, y, 'b-', linewidth=2)

# Use raw strings with LaTeX syntax
ax.set_xlabel(r'Flow Rate $Q$ [m$^3$/h]', fontsize=13)
ax.set_ylabel(r'Head $H$ [m]', fontsize=13)
ax.set_title(r'Pump Curve: $H = H_0 - k Q^2$', fontsize=14, fontweight='bold')

# Greek letters and equations
ax.text(50, 60, r'$\eta = \frac{Q \cdot H \cdot \rho \cdot g}{P_{shaft}}$',
        fontsize=16, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

### Common Engineering Symbols

```python
# Flow rate: Q, $\dot{m}$, $\dot{V}$
# Pressure: p, P, $\Delta p$
# Temperature: T, $\theta$
# Velocity: v, u, V, $\vec{v}$
# Density: $\rho$
# Viscosity: $\mu$, $\nu$
# Efficiency: $\eta$
# Power: P, $\dot{W}$
# Energy: E, $\dot{E}$
# Reynolds number: $Re = \frac{\rho v D}{\mu}$
```

## Saving High-Resolution Plots

### Different Formats

```python
# PNG - raster, good for web
plt.savefig('figure.png', dpi=300, bbox_inches='tight', facecolor='white')

# PDF - vector, good for publications
plt.savefig('figure.pdf', bbox_inches='tight', facecolor='white')

# SVG - vector, editable in Inkscape/Illustrator
plt.savefig('figure.svg', bbox_inches='tight', facecolor='white')

# EPS - vector, required by some journals
plt.savefig('figure.eps', format='eps', bbox_inches='tight')

# TIFF - raster, required by some journals
plt.savefig('figure.tiff', dpi=600, bbox_inches='tight', facecolor='white')
```

### DPI Recommendations

- **Screen display**: 72-96 DPI
- **Presentations**: 150 DPI
- **Print/Publication**: 300-600 DPI
- **Posters**: 150-300 DPI (depends on viewing distance)

### Figure Size Guidelines

```python
# Single column (journal): 3.5 inches width
fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Double column (journal): 7 inches width
fig, ax = plt.subplots(figsize=(7, 5))

# Presentation (16:9):
fig, ax = plt.subplots(figsize=(12, 6.75))

# Presentation (4:3):
fig, ax = plt.subplots(figsize=(10, 7.5))

# Poster:
fig, ax = plt.subplots(figsize=(14, 10))
```

## Best Practices

1. **Always label axes** with units in square brackets: `[m/s]`, `[kPa]`, `[°C]`
2. **Use consistent colors** across related figures
3. **Add legends** when plotting multiple series
4. **Use grid lines** with low alpha (0.3) for readability
5. **Set appropriate line widths** (2-2.5 for main plots)
6. **Choose appropriate colormaps**:
   - Sequential: `viridis`, `plasma`, `inferno` (general data)
   - Diverging: `RdBu_r`, `coolwarm` (data with critical midpoint)
   - Qualitative: `tab10`, `Set1` (categorical data)
7. **Use `tight_layout()`** or `bbox_inches='tight'` to avoid clipped labels
8. **Save before show()** - showing might clear the figure
9. **Close figures** after saving to free memory: `plt.close()`
10. **For CFD results**, consider logarithmic scales for quantities spanning orders of magnitude

## Quick Reference

```python
# Figure and axes
fig, ax = plt.subplots(figsize=(width, height))

# Line plot
ax.plot(x, y, 'b-', linewidth=2, label='Label', marker='o', markersize=5)

# Scatter plot
ax.scatter(x, y, c=colors, s=sizes, alpha=0.6, cmap='viridis')

# Contour plot
contourf = ax.contourf(X, Y, Z, levels=20, cmap='jet')
plt.colorbar(contourf, ax=ax, label='Label')

# Quiver plot
quiver = ax.quiver(X, Y, U, V, scale=50)

# Labels and title
ax.set_xlabel('Label [unit]', fontsize=12)
ax.set_ylabel('Label [unit]', fontsize=12)
ax.set_title('Title', fontsize=14, fontweight='bold')

# Grid and legend
ax.grid(True, alpha=0.3)
ax.legend(loc='best', fontsize=10)

# Limits
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])

# Log scale
ax.set_xscale('log')
ax.set_yscale('log')

# Save
plt.tight_layout()
plt.savefig('figure.png', dpi=300, bbox_inches='tight')
plt.show()
```

## See Also

- `examples.py` - Complete working examples
- `reference.md` - Detailed customization options
- [Matplotlib Documentation](https://matplotlib.org/stable/contents.html)
- [Matplotlib Gallery](https://matplotlib.org/stable/gallery/index.html)
