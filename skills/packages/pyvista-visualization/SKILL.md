---
name: pyvista-visualization
description: "Create 3D visualizations of velocity fields and pump CFD results"
category: packages
domain: visualization
complexity: intermediate
dependencies:
  - pyvista
  - numpy
---

# PyVista 3D Visualization for CFD

Create interactive 3D visualizations of computational fluid dynamics results including velocity fields, pressure contours, streamlines, and pump geometry using PyVista - a high-level 3D visualization library built on VTK.

## Installation and Setup

```bash
pip install pyvista numpy
```

Optional dependencies for enhanced functionality:

```bash
# For mesh processing and additional features
pip install meshio vtk

# For Jupyter notebook support
pip install trame jupyter-server-proxy ipywidgets
```

Basic import pattern:

```python
import pyvista as pv
import numpy as np

# Set plotting theme
pv.set_plot_theme('document')  # or 'default', 'dark', 'paraview'

# Enable off-screen rendering (for scripts without display)
# pv.OFF_SCREEN = True
```

## Core Capabilities for CFD Visualization

### 3D Mesh Visualization

PyVista supports various mesh types essential for CFD:

- **Structured Grids**: Regular grids with curvilinear coordinates
- **Unstructured Grids**: Arbitrary cell types (tetrahedral, hexahedral, etc.)
- **PolyData**: Surface meshes and point clouds
- **Uniform Grids**: Regular Cartesian grids (voxels)

```python
import pyvista as pv
import numpy as np

# Create structured grid (e.g., pump volute domain)
x = np.linspace(0, 1, 20)
y = np.linspace(0, 1, 15)
z = np.linspace(0, 0.5, 10)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Create mesh
grid = pv.StructuredGrid(X, Y, Z)

# Visualize
grid.plot(show_edges=True)
```

### Scalar Field Visualization

Display pressure, temperature, turbulence intensity, and other scalar quantities:

```python
# Add scalar field to mesh (e.g., pressure distribution)
pressure = 101325 + 50000 * np.exp(-((X-0.5)**2 + (Y-0.5)**2 + (Z-0.25)**2) / 0.1)
grid['pressure'] = pressure.flatten(order='F')

# Create plotter with contours
plotter = pv.Plotter()
plotter.add_mesh(grid.contour(isosurfaces=10, scalars='pressure'),
                 cmap='coolwarm',
                 show_scalar_bar=True,
                 scalar_bar_args={'title': 'Pressure [Pa]'})
plotter.show()
```

### Vector Field Visualization (Glyphs)

Visualize velocity vectors using arrows or other glyphs:

```python
# Define velocity field
U = -0.5 * (Y - 0.5)
V = 0.5 * (X - 0.5)
W = 0.2 * np.sin(2 * np.pi * Z)

vectors = np.column_stack([U.flatten(order='F'),
                           V.flatten(order='F'),
                           W.flatten(order='F')])

grid['velocity'] = vectors

# Create arrow glyphs
arrows = grid.glyph(orient='velocity', scale='velocity', factor=0.1)

plotter = pv.Plotter()
plotter.add_mesh(arrows, cmap='viridis', show_scalar_bar=True,
                scalar_bar_args={'title': 'Velocity [m/s]'})
plotter.show()
```

### Streamlines

Trace flow pathlines through velocity fields:

```python
# Create seed points for streamlines
seed_points = pv.Disc(center=(0.5, 0.5, 0), inner=0.05, outer=0.3, normal=(0, 0, 1), r_res=6, c_res=12)

# Generate streamlines
streamlines = grid.streamlines(
    vectors='velocity',
    source_center=(0.5, 0.5, 0),
    source_radius=0.3,
    n_points=50,
    max_time=10.0
)

plotter = pv.Plotter()
plotter.add_mesh(streamlines.tube(radius=0.005), cmap='jet',
                scalar_bar_args={'title': 'Velocity Magnitude [m/s]'})
plotter.add_mesh(seed_points, color='red', point_size=10)
plotter.show()
```

### Volume Rendering

Display 3D scalar fields with transparency:

```python
# Volume rendering for density or temperature fields
plotter = pv.Plotter()
plotter.add_volume(grid, scalars='pressure',
                  cmap='coolwarm',
                  opacity='sigmoid',  # or 'linear', custom array
                  scalar_bar_args={'title': 'Pressure [Pa]'})
plotter.show()
```

### Slicing and Clipping

Extract 2D slices or clip regions from 3D domains:

```python
# Create multiple slices through domain
slices = grid.slice_orthogonal(x=0.5, y=0.5, z=0.25)

plotter = pv.Plotter()
plotter.add_mesh(slices, scalars='pressure', cmap='RdBu_r',
                show_scalar_bar=True,
                scalar_bar_args={'title': 'Pressure [Pa]'})
plotter.show()

# Clip half of domain to see internal flow
clipped = grid.clip(normal='x', value=0.5)

plotter = pv.Plotter()
plotter.add_mesh(clipped, scalars='velocity', cmap='jet',
                show_edges=True, show_scalar_bar=True)
plotter.show()
```

## CFD Application Examples

### Velocity Field Visualization

```python
import pyvista as pv
import numpy as np

# Load or create CFD mesh
grid = pv.ImageData(dimensions=(50, 40, 30))
grid.spacing = (0.02, 0.02, 0.02)
grid.origin = (0, 0, 0)

# Simulate velocity field (e.g., pipe flow with swirl)
x, y, z = grid.points.T
r = np.sqrt((x-0.5)**2 + (y-0.4)**2)
theta = np.arctan2(y-0.4, x-0.5)

u = 2.0 * (1 - (r/0.3)**2) * (r < 0.3)  # Axial velocity
v = 0.5 * r * np.cos(theta) * (r < 0.3)  # Tangential component
w = -0.5 * r * np.sin(theta) * (r < 0.3)

grid['velocity'] = np.column_stack([u, v, w])
grid['speed'] = np.linalg.norm(grid['velocity'], axis=1)

# Visualize with streamlines and contours
plotter = pv.Plotter()

# Add velocity magnitude contours on slices
slice_y = grid.slice(normal='y', origin=(0.5, 0.4, 0.3))
plotter.add_mesh(slice_y, scalars='speed', cmap='jet',
                opacity=0.8, show_scalar_bar=True,
                scalar_bar_args={'title': 'Velocity [m/s]', 'height': 0.7})

# Add streamlines
streamlines = grid.streamlines(
    vectors='velocity',
    source_center=(0.5, 0.4, 0),
    source_radius=0.25,
    n_points=30
)
plotter.add_mesh(streamlines.tube(radius=0.003), color='white')

plotter.show()
```

### Pressure Contours

```python
# Add pressure field (example: stagnation and wake regions)
x_norm = (x - 0.5) / 0.5
y_norm = (y - 0.4) / 0.4
z_norm = (z - 0.3) / 0.3

# Pressure distribution around obstacle
pressure = 101325 + 5000 * (1 - np.sqrt(x_norm**2 + y_norm**2 + z_norm**2))
grid['pressure'] = pressure

# Create pressure isosurfaces
plotter = pv.Plotter()
contours = grid.contour(isosurfaces=10, scalars='pressure')
plotter.add_mesh(contours, cmap='coolwarm', opacity=0.7,
                show_scalar_bar=True,
                scalar_bar_args={'title': 'Pressure [Pa]'})
plotter.show()
```

### Turbulence Visualization

```python
# Turbulent kinetic energy (TKE) visualization
k_turb = 0.01 * grid['speed']**2 * np.random.rand(len(grid['speed']))
grid['TKE'] = k_turb

plotter = pv.Plotter()

# Volume rendering for turbulence
plotter.add_volume(grid, scalars='TKE', cmap='hot',
                  opacity='linear',
                  scalar_bar_args={'title': 'Turbulent Kinetic Energy [m²/s²]'})

# Add outline
plotter.add_mesh(grid.outline(), color='black', line_width=2)

plotter.show()
```

### Pump Geometry Display

```python
# Load pump impeller geometry (assuming STL file)
# impeller = pv.read('pump_impeller.stl')

# Or create simple impeller geometry
def create_simple_impeller():
    """Create simplified pump impeller for demonstration."""
    # Hub (cylinder)
    hub = pv.Cylinder(radius=0.02, height=0.05, center=(0, 0, 0),
                     direction=(0, 0, 1), resolution=30)

    # Blades (create 6 blades)
    blades = pv.PolyData()
    for i in range(6):
        angle = i * 60
        blade = pv.Plane(center=(0.04, 0, 0.025), direction=(0, 0, 1),
                        i_size=0.04, j_size=0.05)
        blade.rotate_z(angle, point=(0, 0, 0.025))
        blades += blade

    impeller = hub + blades
    return impeller

impeller = create_simple_impeller()

# Visualize with lighting
plotter = pv.Plotter()
plotter.add_mesh(impeller, color='lightblue', metallic=0.5,
                roughness=0.5, show_edges=True)
plotter.add_light(pv.Light(position=(1, 1, 1), light_type='scene light'))
plotter.show()
```

### Combined Visualization

```python
# Comprehensive CFD visualization combining multiple techniques
plotter = pv.Plotter()

# Background mesh with velocity magnitude
plotter.add_mesh(grid.outline(), color='black', line_width=2)

# Slice through center showing velocity contours
center_slice = grid.slice(normal='z', origin=(0.5, 0.4, 0.3))
plotter.add_mesh(center_slice, scalars='speed', cmap='jet',
                opacity=0.9, show_scalar_bar=True,
                scalar_bar_args={'title': 'Velocity [m/s]',
                               'vertical': True,
                               'height': 0.7})

# Streamlines showing flow patterns
streamlines = grid.streamlines(
    vectors='velocity',
    source_center=(0.5, 0.4, 0.1),
    source_radius=0.2,
    n_points=25,
    max_time=5.0
)
plotter.add_mesh(streamlines.tube(radius=0.002), color='white', opacity=0.8)

# Pressure isosurface highlighting high-pressure region
high_pressure = grid.threshold(value=103000, scalars='pressure')
plotter.add_mesh(high_pressure, color='red', opacity=0.3)

# Set camera and view
plotter.camera_position = 'xy'
plotter.show()
```

## Interactive Features

### Camera Control

```python
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='pressure')

# Set camera position
plotter.camera_position = [
    (2, 2, 2),    # Camera position
    (0.5, 0.4, 0.3),  # Focal point
    (0, 0, 1)     # View up vector
]

# Or use preset views
# plotter.camera_position = 'xy'  # Top view
# plotter.camera_position = 'xz'  # Front view
# plotter.camera_position = 'yz'  # Side view
# plotter.camera_position = 'iso'  # Isometric

plotter.show()
```

### Multiple Viewports

```python
# Create side-by-side comparison
plotter = pv.Plotter(shape=(1, 2))

# Left: velocity magnitude
plotter.subplot(0, 0)
plotter.add_mesh(grid, scalars='speed', cmap='jet')
plotter.add_text('Velocity Magnitude', font_size=12)

# Right: pressure
plotter.subplot(0, 1)
plotter.add_mesh(grid, scalars='pressure', cmap='coolwarm')
plotter.add_text('Pressure', font_size=12)

plotter.link_views()  # Synchronize camera movements
plotter.show()
```

### Animation

```python
# Animate rotating view
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='speed', cmap='jet')

# Save as GIF or MP4
# plotter.open_gif('rotation.gif')
# plotter.open_movie('rotation.mp4')

# Create rotation animation
path = plotter.generate_orbital_path(n_points=36, shift=0)
plotter.orbit_on_path(path, write_frames=False)
```

## Exporting and Saving

### Save Static Images

```python
plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(grid, scalars='pressure', cmap='coolwarm')
plotter.camera_position = 'iso'

# Save high-resolution image
plotter.screenshot('cfd_pressure.png', transparent_background=False,
                  window_size=[1920, 1080])
```

### Export Meshes

```python
# Save mesh with data for later use
grid.save('cfd_results.vtk')

# Export to other formats
grid.save('cfd_results.vtu')  # VTK Unstructured Grid
grid.save('cfd_results.vtp')  # VTK PolyData
```

### Export to ParaView

PyVista meshes are VTK-compatible and can be opened directly in ParaView for advanced post-processing.

## Best Practices for CFD Visualization

1. **Choose appropriate colormaps**:
   - Velocity: `'jet'`, `'viridis'`, `'plasma'`
   - Pressure/Temperature: `'coolwarm'`, `'RdBu_r'`
   - Turbulence: `'hot'`, `'inferno'`

2. **Use consistent scales**: Set `clim` (color limits) for comparing multiple cases

3. **Add context**: Include geometry outlines, coordinate axes, and scale bars

4. **Optimize for large datasets**:
   - Use decimation: `mesh.decimate(0.5)` to reduce points
   - Enable LOD (Level of Detail): `plotter.enable_eye_dome_lighting()`

5. **Label effectively**: Add titles, scalar bar labels with units

6. **Consider lighting**: Add custom lights for better 3D perception

7. **Use transparency wisely**: Combine opaque surfaces with transparent volumes

## Quick Reference

```python
import pyvista as pv
import numpy as np

# Create mesh
grid = pv.ImageData(dimensions=(50, 40, 30))
grid.spacing = (0.02, 0.02, 0.02)

# Add scalar field
grid['pressure'] = np.random.rand(grid.n_points)

# Add vector field
grid['velocity'] = np.random.rand(grid.n_points, 3)

# Create plotter
plotter = pv.Plotter()

# Add mesh with options
plotter.add_mesh(grid,
                scalars='pressure',    # Scalar field name
                cmap='coolwarm',       # Colormap
                opacity=0.8,           # Transparency
                show_edges=True,       # Show cell edges
                show_scalar_bar=True,  # Show colorbar
                clim=[0, 1])          # Color limits

# Add streamlines
streamlines = grid.streamlines(vectors='velocity', n_points=20)
plotter.add_mesh(streamlines.tube(radius=0.001), color='white')

# Slice
slice_z = grid.slice(normal='z')
plotter.add_mesh(slice_z, scalars='pressure')

# Contour isosurfaces
contours = grid.contour(isosurfaces=10, scalars='pressure')
plotter.add_mesh(contours, opacity=0.5)

# Camera and view
plotter.camera_position = 'iso'
plotter.add_axes()
plotter.show_grid()

# Display
plotter.show()

# Save
plotter.screenshot('output.png')
```

## See Also

- `examples.py` - Complete working CFD visualization examples
- `reference.md` - Detailed PyVista API reference and options
- [PyVista Documentation](https://docs.pyvista.org/)
- [PyVista Examples Gallery](https://docs.pyvista.org/examples/index.html)
- [VTK File Formats](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
