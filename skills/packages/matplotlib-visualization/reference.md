# Matplotlib Visualization Reference

Comprehensive reference for customization options, colormaps, and best practices for engineering visualizations.

## Table of Contents

- [Common Customization Options](#common-customization-options)
- [Color Maps for Engineering](#color-maps-for-engineering)
- [Figure Size and DPI Recommendations](#figure-size-and-dpi-recommendations)
- [Line Styles and Markers](#line-styles-and-markers)
- [Font and Text Customization](#font-and-text-customization)
- [Axis Configuration](#axis-configuration)
- [Legend Customization](#legend-customization)
- [Grid and Ticks](#grid-and-ticks)
- [Color Palettes](#color-palettes)

## Common Customization Options

### Figure Creation

```python
# Basic figure with single axes
fig, ax = plt.subplots(figsize=(10, 6))

# Multiple subplots (2x2 grid)
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Subplots with shared axes
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 8))

# Using GridSpec for complex layouts
from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize=(12, 8))
gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
ax1 = fig.add_subplot(gs[0, :])    # Top row, all columns
ax2 = fig.add_subplot(gs[1:, 0])   # Bottom rows, first column
ax3 = fig.add_subplot(gs[1:, 1:])  # Bottom rows, remaining columns
```

### Line Plot Options

```python
ax.plot(x, y,
        color='blue',              # Color: name, hex, RGB tuple
        linestyle='-',             # '-', '--', '-.', ':', ''
        linewidth=2,               # Line width in points
        marker='o',                # Marker style (see below)
        markersize=6,              # Marker size
        markerfacecolor='red',     # Fill color
        markeredgecolor='black',   # Edge color
        markeredgewidth=1,         # Edge width
        alpha=0.7,                 # Transparency (0-1)
        label='Data series',       # Legend label
        zorder=5)                  # Drawing order (higher = on top)
```

### Scatter Plot Options

```python
ax.scatter(x, y,
           s=50,                   # Size (scalar or array)
           c='blue',               # Color (scalar, array, or sequence)
           marker='o',             # Marker style
           cmap='viridis',         # Colormap (if c is array)
           alpha=0.6,              # Transparency
           edgecolors='black',     # Edge color
           linewidths=1,           # Edge width
           vmin=0, vmax=100,       # Color scale limits
           label='Data points')
```

### Contour Plot Options

```python
# Filled contours
contourf = ax.contourf(X, Y, Z,
                       levels=20,          # Number or list of levels
                       cmap='jet',         # Colormap
                       vmin=0, vmax=100,   # Data range
                       extend='both',      # 'neither', 'min', 'max', 'both'
                       alpha=0.9)

# Contour lines
contour = ax.contour(X, Y, Z,
                     levels=10,
                     colors='black',       # Single color or list
                     linewidths=1,
                     linestyles='solid',   # 'solid', 'dashed', 'dashdot', 'dotted'
                     alpha=0.5)

# Add labels to contour lines
ax.clabel(contour, inline=True, fontsize=8, fmt='%.1f')

# Add colorbar
cbar = plt.colorbar(contourf, ax=ax,
                    orientation='vertical',  # 'vertical' or 'horizontal'
                    pad=0.02,                # Space between plot and colorbar
                    shrink=0.8,              # Scale relative to axes
                    aspect=20,               # Ratio of long to short dimensions
                    extend='both',           # Arrow indicators
                    label='Label [units]')
```

### Quiver Plot Options

```python
quiver = ax.quiver(X, Y, U, V,
                   C,                      # Color data (optional)
                   cmap='viridis',
                   scale=50,               # Data units per arrow length unit
                   scale_units='width',    # 'width', 'height', 'dots', 'inches'
                   width=0.003,            # Arrow shaft width
                   headwidth=3,            # Head width as multiple of shaft
                   headlength=5,           # Head length as multiple of shaft
                   headaxislength=4.5,     # Head length at shaft intersection
                   alpha=0.8,
                   pivot='mid')            # 'tail', 'mid', 'middle', 'tip'

# Add a reference arrow
ax.quiverkey(quiver, X=0.9, Y=0.95, U=2,
             label='2 m/s', labelpos='E',
             coordinates='axes')
```

## Color Maps for Engineering

### Sequential Colormaps

Use for data that progresses from low to high values.

**Perceptually uniform (recommended):**
- `viridis` - Default, excellent for most applications
- `plasma` - High contrast, good for presentations
- `inferno` - Dark background compatible
- `magma` - Similar to inferno, slightly different hue
- `cividis` - Optimized for color-vision deficiency

**Traditional (use with caution):**
- `jet` - Classic CFD colormap (not perceptually uniform)
- `hot` - Black-red-yellow-white
- `cool` - Cyan to magenta
- `gray` - Grayscale

```python
# Example usage
contourf = ax.contourf(X, Y, Z, levels=20, cmap='viridis')

# Reverse colormap
contourf = ax.contourf(X, Y, Z, levels=20, cmap='viridis_r')
```

### Diverging Colormaps

Use for data with a meaningful center point (e.g., deviations from nominal, positive/negative values).

- `RdBu_r` - Red (negative) to Blue (positive), reversed
- `coolwarm` - Blue (low) to red (high), perceptually uniform
- `seismic` - Blue-white-red
- `bwr` - Blue-white-red (more saturated than seismic)
- `RdYlBu_r` - Red-yellow-blue, reversed

```python
# Example: Pressure deviation from nominal
deviation = pressure - nominal_pressure
contourf = ax.contourf(X, Y, deviation, levels=20, cmap='RdBu_r',
                       vmin=-10, vmax=10)  # Symmetric limits
```

### Qualitative Colormaps

Use for categorical data or multiple line series.

- `tab10` - 10 distinct colors (default)
- `tab20` - 20 distinct colors
- `Set1` - 9 colors, high saturation
- `Set2` - 8 colors, medium saturation
- `Set3` - 12 colors, low saturation
- `Paired` - 12 colors in pairs

```python
# Example: Multiple pump curves
colors = plt.cm.tab10.colors
for i, data in enumerate(datasets):
    ax.plot(x, data, color=colors[i], label=f'Series {i+1}')
```

### Engineering-Specific Recommendations

| Application | Recommended Colormap | Reason |
|-------------|---------------------|---------|
| Velocity magnitude | `viridis`, `plasma` | Perceptually uniform, shows detail |
| Temperature | `hot`, `inferno` | Intuitive (hot = red) |
| Pressure | `viridis`, `coolwarm` (diverging) | Clear progression or deviation |
| Vorticity | `RdBu_r` | Shows positive/negative rotation |
| Turbulence intensity | `YlOrRd` | Progressive intensity |
| Strain rate | `YlGnBu` | Sequential with good contrast |
| Streamlines | Single color or grayscale | Focus on flow direction |
| Printing (B&W) | `gray`, `viridis` | Converts well to grayscale |

## Figure Size and DPI Recommendations

### Journal Submissions

Most journals specify figure dimensions and resolution:

```python
# Single column (typically 3.5 inches)
fig, ax = plt.subplots(figsize=(3.5, 2.5))
plt.savefig('figure.png', dpi=600)  # 600 DPI for print

# Double column (typically 7 inches)
fig, ax = plt.subplots(figsize=(7, 5))
plt.savefig('figure.png', dpi=600)

# Vector format (preferred by many journals)
plt.savefig('figure.pdf')  # DPI doesn't matter for vector
```

### Common Journal Requirements

| Journal/Publisher | Single Column | Double Column | DPI |
|-------------------|---------------|---------------|-----|
| Elsevier | 90 mm (3.54") | 190 mm (7.48") | 300-500 |
| Springer | 84 mm (3.31") | 174 mm (6.85") | 300-600 |
| IEEE | 3.5" | 7.16" | 300-600 |
| AIP | 3.37" | 6.69" | 300 |
| ASME | 3.5" | 7" | 300-600 |

### Presentations

```python
# 16:9 aspect ratio (standard widescreen)
fig, ax = plt.subplots(figsize=(12, 6.75))
plt.savefig('slide.png', dpi=150)

# 4:3 aspect ratio (older projectors)
fig, ax = plt.subplots(figsize=(10, 7.5))
plt.savefig('slide.png', dpi=150)

# Full HD resolution (1920x1080)
fig, ax = plt.subplots(figsize=(19.2, 10.8))
plt.savefig('slide.png', dpi=100)
```

### Posters

```python
# Large format (24" x 36" poster)
fig, ax = plt.subplots(figsize=(18, 12))
plt.savefig('poster_figure.png', dpi=150)

# Very large format (viewing distance > 2m)
fig, ax = plt.subplots(figsize=(20, 15))
plt.savefig('poster_figure.png', dpi=100)
```

### DPI Guidelines

| Medium | Recommended DPI | Notes |
|--------|----------------|-------|
| Screen/Web | 72-96 | Standard screen resolution |
| Presentation | 150 | Good for projectors |
| Print (draft) | 150-200 | Quick review prints |
| Print (publication) | 300-600 | Standard journal quality |
| Print (high-quality) | 600-1200 | Art books, large format |
| Vector (PDF/SVG) | N/A | Resolution-independent |

## Line Styles and Markers

### Line Styles

```python
# Line style codes
'-'   # solid line (default)
'--'  # dashed line
'-.'  # dash-dot line
':'   # dotted line
''    # no line

# Example
ax.plot(x, y1, '-', label='Solid')
ax.plot(x, y2, '--', label='Dashed')
ax.plot(x, y3, '-.', label='Dash-dot')
ax.plot(x, y4, ':', label='Dotted')
```

### Marker Styles

```python
# Common markers
'o'   # circle
's'   # square
'^'   # triangle up
'v'   # triangle down
'D'   # diamond
'*'   # star
'+'   # plus
'x'   # x
'.'   # point

# Example
ax.plot(x, y, 'o-', markersize=6, label='With markers')
ax.plot(x, y, 'o', markersize=8, label='Markers only')
```

### Combined Format Strings

```python
# Format: [marker][line][color]
ax.plot(x, y, 'ro-')   # red circles with solid line
ax.plot(x, y, 'bs--')  # blue squares with dashed line
ax.plot(x, y, 'g^:')   # green triangles with dotted line
```

## Font and Text Customization

### Font Properties

```python
# Global font settings
plt.rcParams['font.family'] = 'serif'  # 'serif', 'sans-serif', 'monospace'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']
plt.rcParams['font.size'] = 11

# Individual text elements
ax.set_xlabel('Label', fontsize=12, fontweight='bold', fontstyle='italic')
ax.set_title('Title', fontsize=14, fontfamily='sans-serif')
```

### LaTeX Support

```python
# Enable LaTeX (requires LaTeX installation)
plt.rcParams['text.usetex'] = True

# Use math text (no LaTeX required)
plt.rcParams['text.usetex'] = False
plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern font

# Examples
ax.set_xlabel(r'Flow rate $Q$ [m$^3$/h]')
ax.set_ylabel(r'$\eta$ [%]')
ax.set_title(r'$Re = \frac{\rho v D}{\mu}$')
```

### Common Greek Letters

```python
# LaTeX notation (with r-string)
'$\\alpha$'   # α (alpha)
'$\\beta$'    # β (beta)
'$\\gamma$'   # γ (gamma)
'$\\delta$'   # δ (delta)
'$\\epsilon$' # ε (epsilon)
'$\\eta$'     # η (eta)
'$\\theta$'   # θ (theta)
'$\\lambda$'  # λ (lambda)
'$\\mu$'      # μ (mu)
'$\\nu$'      # ν (nu)
'$\\rho$'     # ρ (rho)
'$\\sigma$'   # σ (sigma)
'$\\tau$'     # τ (tau)
'$\\phi$'     # φ (phi)
'$\\omega$'   # ω (omega)

# Capital letters
'$\\Delta$'   # Δ (Delta)
'$\\Omega$'   # Ω (Omega)
```

### Text Annotations

```python
# Add text at specific location
ax.text(x, y, 'Text', fontsize=11,
        horizontalalignment='center',  # 'left', 'center', 'right'
        verticalalignment='bottom',     # 'top', 'center', 'bottom'
        transform=ax.transData,         # Coordinate system
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Annotate with arrow
ax.annotate('Operating point',
            xy=(50, 60),              # Point to annotate
            xytext=(70, 80),          # Text location
            arrowprops=dict(arrowstyle='->', lw=1.5, color='red'),
            fontsize=11,
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
```

## Axis Configuration

### Axis Limits and Scaling

```python
# Set limits
ax.set_xlim([0, 100])
ax.set_ylim([0, 80])

# Automatic with padding
ax.margins(0.1)  # 10% padding

# Logarithmic scale
ax.set_xscale('log')
ax.set_yscale('log')

# Symmetric log scale (handles negative values)
ax.set_yscale('symlog', linthresh=1)

# Invert axis
ax.invert_xaxis()
ax.invert_yaxis()
```

### Aspect Ratio

```python
# Equal aspect (for physical x-y plots)
ax.set_aspect('equal')

# Specific aspect ratio
ax.set_aspect(2.0)  # y-unit is 2x longer than x-unit

# Automatic (fill available space)
ax.set_aspect('auto')
```

### Multiple Y-Axes

```python
fig, ax1 = plt.subplots(figsize=(10, 6))

# First y-axis
ax1.plot(x, y1, 'b-', label='Head')
ax1.set_xlabel('Flow rate [m³/h]')
ax1.set_ylabel('Head [m]', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Second y-axis
ax2 = ax1.twinx()
ax2.plot(x, y2, 'r-', label='Power')
ax2.set_ylabel('Power [kW]', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Combine legends
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
```

### Axis Spine Customization

```python
# Hide top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Move spines to zero
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')

# Change spine properties
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_color('blue')
```

## Legend Customization

```python
ax.legend(loc='best',              # 'best', 'upper right', 'lower left', etc.
          fontsize=10,
          frameon=True,            # Show frame
          framealpha=0.9,          # Frame transparency
          edgecolor='black',       # Frame edge color
          fancybox=True,           # Rounded corners
          shadow=True,             # Drop shadow
          ncol=2,                  # Number of columns
          title='Legend Title',    # Legend title
          title_fontsize=11)

# Legend outside plot area
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Custom legend entries
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='blue', lw=2),
                Line2D([0], [0], color='red', lw=2, linestyle='--')]
ax.legend(custom_lines, ['Label 1', 'Label 2'])
```

## Grid and Ticks

### Grid Configuration

```python
# Basic grid
ax.grid(True)

# Customized grid
ax.grid(True,
        which='major',      # 'major', 'minor', 'both'
        axis='both',        # 'x', 'y', 'both'
        linestyle='-',
        linewidth=0.5,
        alpha=0.3,
        color='gray')

# Minor grid
ax.minorticks_on()
ax.grid(which='minor', linestyle=':', alpha=0.2)
```

### Tick Configuration

```python
# Tick positions
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_yticks(np.arange(0, 101, 10))

# Tick labels
ax.set_xticklabels(['0', '25', '50', '75', '100'])

# Tick parameters
ax.tick_params(axis='both',
               which='major',       # 'major', 'minor', 'both'
               direction='in',      # 'in', 'out', 'inout'
               length=6,
               width=1,
               labelsize=11,
               labelrotation=45,
               pad=10)              # Distance from axis

# Scientific notation
from matplotlib.ticker import ScalarFormatter
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
```

## Color Palettes

### Professional Color Schemes

```python
# Tableau colors (professional, colorblind-friendly)
tableau_colors = {
    'blue': '#1F77B4',
    'orange': '#FF7F0E',
    'green': '#2CA02C',
    'red': '#D62728',
    'purple': '#9467BD',
    'brown': '#8C564B',
    'pink': '#E377C2',
    'gray': '#7F7F7F',
    'olive': '#BCBD22',
    'cyan': '#17BECF'
}

# MATLAB-style colors
matlab_colors = {
    'blue': '#0072BD',
    'red': '#D95319',
    'yellow': '#EDB120',
    'purple': '#7E2F8E',
    'green': '#77AC30',
    'cyan': '#4DBEEE',
    'maroon': '#A2142F'
}

# Colorblind-friendly palette (Wong 2011)
colorblind_palette = {
    'blue': '#0173B2',
    'orange': '#DE8F05',
    'bluish_green': '#029E73',
    'yellow': '#FBAFE4',
    'sky_blue': '#56B4E9',
    'vermillion': '#D55E00',
    'reddish_purple': '#CC79A7',
    'black': '#000000'
}

# Usage
ax.plot(x, y, color=tableau_colors['blue'], linewidth=2)
```

### Creating Custom Colormaps

```python
from matplotlib.colors import LinearSegmentedColormap

# Define colors
colors = ['blue', 'cyan', 'yellow', 'red']
n_bins = 100
cmap = LinearSegmentedColormap.from_list('custom', colors, N=n_bins)

# Use custom colormap
contourf = ax.contourf(X, Y, Z, levels=20, cmap=cmap)
```

## Complete Example with All Options

```python
import matplotlib.pyplot as plt
import numpy as np

# Set global parameters
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = 11
plt.rcParams['axes.linewidth'] = 1.2

# Create figure
fig, ax = plt.subplots(figsize=(10, 6))

# Generate data
x = np.linspace(0, 100, 100)
y1 = 80 - 0.005 * x**2
y2 = 70 - 0.004 * x**2

# Plot with full customization
ax.plot(x, y1, color='#0072BD', linestyle='-', linewidth=2.5,
        marker='o', markevery=10, markersize=6,
        label='Pump A', zorder=3)
ax.plot(x, y2, color='#D95319', linestyle='--', linewidth=2.5,
        marker='s', markevery=10, markersize=6,
        label='Pump B', zorder=3)

# Axis labels with LaTeX
ax.set_xlabel('Flow rate $Q$ [m$^3$/h]', fontsize=13, fontweight='bold')
ax.set_ylabel('Head $H$ [m]', fontsize=13, fontweight='bold')

# Title
ax.set_title('Pump Performance Comparison',
             fontsize=15, fontweight='bold', pad=15)

# Grid
ax.grid(True, which='major', linestyle='-', linewidth=0.5,
        alpha=0.3, zorder=0)

# Legend
ax.legend(loc='upper right', fontsize=11, frameon=True,
          framealpha=0.9, edgecolor='black')

# Limits and ticks
ax.set_xlim([0, 100])
ax.set_ylim([0, 90])
ax.set_xticks(np.arange(0, 101, 20))
ax.set_yticks(np.arange(0, 91, 10))

# Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Tight layout
plt.tight_layout()

# Save
plt.savefig('customized_plot.png', dpi=300, bbox_inches='tight',
            facecolor='white')
plt.savefig('customized_plot.pdf', bbox_inches='tight')

plt.show()
```

## Quick Tips

1. **Use consistent styling** across all figures in a publication
2. **Save figures before calling `show()`** to avoid empty files
3. **Use vector formats (PDF, SVG)** when possible for publications
4. **Close figures** to free memory: `plt.close()` or `plt.close('all')`
5. **Test colorblind accessibility** using online simulators
6. **Check journal requirements** before creating final figures
7. **Use `tight_layout()`** or `bbox_inches='tight'` to avoid clipping
8. **Backup original data** used to create figures for reproducibility
9. **Document all custom settings** in your code
10. **Preview at target size** before finalizing (especially for journals)

## See Also

- [Matplotlib Documentation](https://matplotlib.org/stable/contents.html)
- [Matplotlib Gallery](https://matplotlib.org/stable/gallery/index.html)
- [Choosing Colormaps](https://matplotlib.org/stable/tutorials/colors/colormaps.html)
- [Scientific Visualization Best Practices](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833)
