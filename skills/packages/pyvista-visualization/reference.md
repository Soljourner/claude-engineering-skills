# PyVista Visualization Reference

Comprehensive reference for PyVista objects, plotting options, and customization for CFD visualization.

## Table of Contents

- [Key PyVista Objects](#key-pyvista-objects)
- [Mesh Types](#mesh-types)
- [Plotting and Rendering](#plotting-and-rendering)
- [Data Visualization Methods](#data-visualization-methods)
- [Camera and View Control](#camera-and-view-control)
- [Color Maps and Styling](#color-maps-and-styling)
- [Export and Save Options](#export-and-save-options)
- [Performance Optimization](#performance-optimization)
- [Common Patterns for CFD](#common-patterns-for-cfd)

## Key PyVista Objects

### Mesh Objects

PyVista provides several mesh types for different data structures:

```python
import pyvista as pv
import numpy as np

# StructuredGrid - curvilinear coordinates (i,j,k indexing)
# Used for: Pump volutes, blade passages, structured CFD grids
x, y, z = np.meshgrid(...)  # 3D arrays
grid = pv.StructuredGrid(x, y, z)

# ImageData (UniformGrid) - regular Cartesian grid
# Used for: Simple rectangular domains, voxel data
grid = pv.ImageData(dimensions=(nx, ny, nz))
grid.spacing = (dx, dy, dz)
grid.origin = (x0, y0, z0)

# RectilinearGrid - non-uniform but axis-aligned
# Used for: Domains with varying resolution in different directions
grid = pv.RectilinearGrid(x_coords, y_coords, z_coords)

# UnstructuredGrid - arbitrary cell types
# Used for: Tetrahedral/hexahedral meshes from FEA/CFD solvers
grid = pv.UnstructuredGrid(cells, cell_types, points)

# PolyData - surface meshes and point clouds
# Used for: Geometry surfaces, STL files, streamlines
mesh = pv.PolyData(points, faces)
```

### Plotter Object

The main visualization interface:

```python
# Basic plotter
plotter = pv.Plotter()

# Multiple viewports (2x2 grid)
plotter = pv.Plotter(shape=(2, 2))

# Off-screen rendering (for scripts)
plotter = pv.Plotter(off_screen=True)

# Custom window size
plotter = pv.Plotter(window_size=[1920, 1080])

# Notebook mode (for Jupyter)
plotter = pv.Plotter(notebook=True)
```

## Mesh Types

### Creating Meshes

```python
# From arrays
x = np.linspace(0, 10, 50)
y = np.linspace(0, 5, 25)
z = np.linspace(0, 2, 20)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
grid = pv.StructuredGrid(X, Y, Z)

# ImageData (fastest for regular grids)
grid = pv.ImageData(dimensions=(50, 25, 20))
grid.spacing = (0.2, 0.2, 0.1)
grid.origin = (0, 0, 0)

# From file (VTK, VTU, VTP, STL, PLY, OBJ, etc.)
mesh = pv.read('data.vtu')
mesh = pv.read('geometry.stl')

# Primitive shapes
sphere = pv.Sphere(radius=1.0, center=(0, 0, 0))
cylinder = pv.Cylinder(radius=0.5, height=2.0)
box = pv.Box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))
cone = pv.Cone(center=(0, 0, 0), direction=(0, 0, 1))
```

### Adding Data to Meshes

```python
# Scalar data (one value per point or cell)
mesh['pressure'] = pressure_array  # Point data (size = n_points)
mesh['cell_volume'] = volume_array  # Cell data (size = n_cells)

# Vector data (3 components per point or cell)
mesh['velocity'] = velocity_array  # Shape: (n_points, 3)

# Access data
pressure = mesh['pressure']
velocity = mesh['velocity']

# Point vs cell data
mesh.point_data['temperature'] = temp_points
mesh.cell_data['density'] = density_cells

# List available data
print(mesh.point_data.keys())
print(mesh.cell_data.keys())
```

## Plotting and Rendering

### Basic Plotting

```python
plotter = pv.Plotter()

# Add mesh to scene
plotter.add_mesh(mesh,
                scalars='pressure',      # Data to color by
                cmap='coolwarm',         # Colormap
                opacity=0.8,             # Transparency (0-1)
                show_edges=True,         # Show cell edges
                edge_color='black',      # Edge color
                line_width=1,            # Edge thickness
                show_scalar_bar=True,    # Show colorbar
                clim=[pmin, pmax])       # Color limits

# Display
plotter.show()
```

### Mesh Styling Options

```python
plotter.add_mesh(mesh,
    # Color and appearance
    color='blue',                    # Solid color (if no scalars)
    scalars='field_name',           # Field to visualize
    cmap='viridis',                 # Colormap

    # Transparency
    opacity=0.5,                    # Overall opacity (0-1)
    opacity='field_name',           # Variable opacity

    # Surface properties
    metallic=0.5,                   # Metallic appearance (0-1)
    roughness=0.3,                  # Surface roughness (0-1)
    specular=0.5,                   # Specular highlight (0-1)

    # Edges
    show_edges=True,                # Display cell edges
    edge_color='black',             # Edge color
    line_width=2,                   # Edge thickness

    # Lighting
    lighting=True,                  # Enable lighting
    ambient=0.3,                    # Ambient light (0-1)
    diffuse=0.7,                    # Diffuse reflection (0-1)
    specular=0.2,                   # Specular reflection (0-1)

    # Scalar bar
    show_scalar_bar=True,           # Show colorbar
    scalar_bar_args={...},          # Colorbar options

    # Rendering style
    style='surface',                # 'surface', 'wireframe', 'points'
    render_points_as_spheres=True,  # For point clouds
    point_size=5,                   # Point size

    # Color limits
    clim=[vmin, vmax],             # Manual color range

    # Other
    name='mesh1',                   # Actor name for later reference
    pickable=True,                  # Enable picking
    render=True)                    # Render immediately
```

### Scalar Bar (Colorbar) Customization

```python
scalar_bar_args = {
    'title': 'Pressure [Pa]',           # Label
    'vertical': True,                   # Orientation
    'height': 0.7,                      # Height (0-1)
    'width': 0.05,                      # Width (0-1)
    'position_x': 0.85,                 # X position (0-1)
    'position_y': 0.15,                 # Y position (0-1)
    'title_font_size': 12,              # Title size
    'label_font_size': 10,              # Label size
    'n_labels': 5,                      # Number of labels
    'fmt': '%.1f',                      # Number format
    'font_family': 'arial',             # Font
    'color': 'black',                   # Text color
    'shadow': False,                    # Text shadow
    'bold': False,                      # Bold text
    'italic': False,                    # Italic text
}

plotter.add_mesh(mesh, scalars='pressure',
                scalar_bar_args=scalar_bar_args)
```

## Data Visualization Methods

### Slicing

```python
# Single slice
slice_x = mesh.slice(normal='x', origin=(5, 0, 0))
slice_y = mesh.slice(normal='y', origin=(0, 2.5, 0))
slice_z = mesh.slice(normal='z', origin=(0, 0, 1))

# Arbitrary plane
slice_plane = mesh.slice(normal=(1, 1, 0), origin=(0, 0, 0))

# Multiple parallel slices
slices = mesh.slice_along_axis(n=10, axis='x')

# Orthogonal slices (x, y, z at once)
slices = mesh.slice_orthogonal(x=5, y=2.5, z=1)
```

### Contouring (Isosurfaces)

```python
# Multiple isosurfaces
contours = mesh.contour(isosurfaces=10,      # Number of levels
                       scalars='pressure',   # Field to contour
                       method='marching_cubes')  # Algorithm

# Specific values
contours = mesh.contour(isosurfaces=[100, 200, 300],
                       scalars='temperature')

# Single isosurface
iso = mesh.contour(isosurfaces=1, scalars='pressure',
                  rng=[p_value, p_value])
```

### Thresholding

```python
# Keep values above threshold
high_p = mesh.threshold(value=105000, scalars='pressure')

# Keep values in range
mid_p = mesh.threshold(value=[100000, 110000], scalars='pressure')

# Keep values below threshold
low_p = mesh.threshold(value=100000, scalars='pressure', invert=True)
```

### Clipping

```python
# Clip by plane
clipped = mesh.clip(normal='x', origin=(5, 0, 0))

# Clip by arbitrary plane
clipped = mesh.clip(normal=(1, 1, 1), origin=(0, 0, 0))

# Clip by box
clipped = mesh.clip_box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))

# Keep both sides
clipped, other = mesh.clip(normal='x', return_clipped=True)
```

### Streamlines

```python
# From source surface
seed = pv.Disc(center=(0, 0, 0), inner=0, outer=0.5)
streamlines = mesh.streamlines_from_source(
    source=seed,
    vectors='velocity',            # Vector field name
    max_time=10.0,                 # Maximum integration time
    initial_step_length=0.01,      # Step size
    max_step_length=0.1,           # Maximum step
    integration_direction='both',  # 'forward', 'backward', 'both'
    interpolator_type='point',     # 'point' or 'cell'
    max_steps=2000)                # Maximum number of steps

# From center and radius
streamlines = mesh.streamlines(
    vectors='velocity',
    source_center=(0, 0, 0),
    source_radius=0.5,
    n_points=50,                   # Number of seed points
    max_time=10.0)

# Visualize as tubes
tubes = streamlines.tube(radius=0.01, n_sides=6)
plotter.add_mesh(tubes, scalars='speed', cmap='jet')
```

### Glyphs (Vector Arrows)

```python
# Create arrow glyphs
arrows = mesh.glyph(
    orient='velocity',      # Vector field for orientation
    scale='velocity',       # Field for scaling (or True for orient field)
    factor=0.1,            # Scaling factor
    geom=pv.Arrow(),       # Glyph geometry (Arrow, Cone, Sphere, etc.)
    tolerance=0.0,         # Point merging tolerance
    absolute=False)        # Absolute vs relative scaling

# Other glyph geometries
arrows = mesh.glyph(geom=pv.Arrow())
cones = mesh.glyph(geom=pv.Cone())
spheres = mesh.glyph(geom=pv.Sphere(radius=0.1))

# Custom glyph geometry
custom = pv.Cylinder(radius=0.05, height=0.2)
glyphs = mesh.glyph(orient='velocity', scale=False, geom=custom)
```

### Volume Rendering

```python
# Simple volume rendering
plotter.add_volume(mesh,
                  scalars='density',
                  cmap='viridis',
                  opacity='linear',    # 'linear', 'sigmoid', or array
                  shade=True)          # Enable shading

# Custom opacity transfer function
opacity = np.array([[0, 0],      # [value, opacity]
                   [50, 0.2],
                   [100, 0.8],
                   [150, 1.0]])
plotter.add_volume(mesh, scalars='temperature',
                  opacity=opacity, cmap='hot')

# Preset opacity functions
# 'linear': Linear ramp from 0 to 1
# 'sigmoid': S-shaped curve
# 'sigmoid_3', 'sigmoid_5', etc.: Different sigmoid steepness
```

## Camera and View Control

### Camera Position

```python
# Set camera position, focal point, and view up vector
plotter.camera_position = [
    (x_cam, y_cam, z_cam),     # Camera position
    (x_focal, y_focal, z_focal),  # Focal point
    (vup_x, vup_y, vup_z)      # View up direction
]

# Preset views
plotter.camera_position = 'xy'   # Top view (looking down z)
plotter.camera_position = 'xz'   # Front view (looking along y)
plotter.camera_position = 'yz'   # Side view (looking along x)
plotter.camera_position = 'iso'  # Isometric view

# Get current camera position
position = plotter.camera_position
```

### View Control

```python
# Reset camera to show all actors
plotter.reset_camera()

# Zoom
plotter.camera.zoom(1.5)  # Zoom in
plotter.camera.zoom(0.5)  # Zoom out

# Rotate camera
plotter.camera.azimuth(45)    # Rotate around up vector
plotter.camera.elevation(30)   # Rotate around right vector
plotter.camera.roll(10)        # Rotate around view direction

# Link multiple viewports
plotter.link_views()  # Synchronize camera in all subplots
```

### Viewing Angles

```python
# View from specific direction
plotter.view_xy()      # Top
plotter.view_xz()      # Front
plotter.view_yz()      # Side
plotter.view_isometric()  # Isometric

# Set parallel projection (no perspective)
plotter.enable_parallel_projection()
plotter.disable_parallel_projection()
```

## Color Maps and Styling

### Colormaps for CFD

```python
# Sequential colormaps (good for scalar fields)
'viridis'    # Perceptually uniform, good default
'plasma'     # High contrast
'inferno'    # Dark to bright
'magma'      # Similar to inferno
'jet'        # Classic (not perceptually uniform, but familiar)
'hot'        # Black-red-yellow-white
'cool'       # Cyan to magenta
'gray'       # Grayscale

# Diverging colormaps (good for pressure, temperature differences)
'coolwarm'   # Blue-white-red
'RdBu_r'     # Red-blue reversed
'seismic'    # Blue-white-red (stronger colors)
'bwr'        # Blue-white-red

# Categorical colormaps
'tab10'      # 10 distinct colors
'Set1'       # 9 distinct colors

# Reversed colormaps: add '_r' suffix
'viridis_r', 'jet_r', 'coolwarm_r', etc.
```

### Plot Themes

```python
# Set global theme
pv.set_plot_theme('default')    # Light background
pv.set_plot_theme('dark')       # Dark background
pv.set_plot_theme('document')   # White background, good for publications
pv.set_plot_theme('paraview')   # ParaView-like appearance

# Custom theme
my_theme = pv.themes.DefaultTheme()
my_theme.background = 'white'
my_theme.color = 'black'
my_theme.font.family = 'arial'
my_theme.font.size = 12
my_theme.show_edges = False
pv.global_theme.load_theme(my_theme)
```

### Lighting

```python
# Add custom lights
light = pv.Light(position=(10, 10, 10),
                focal_point=(0, 0, 0),
                color='white',
                intensity=0.8,
                light_type='scene light')  # or 'camera light', 'headlight'
plotter.add_light(light)

# Enable/disable default lights
plotter.enable_lightkit()  # Studio-quality lighting
plotter.disable_lightkit()

# Eye dome lighting (enhances depth perception)
plotter.enable_eye_dome_lighting()
```

## Export and Save Options

### Save Images

```python
# Save screenshot
plotter.screenshot('output.png',
                  transparent_background=False,
                  window_size=[1920, 1080],
                  return_img=False)  # Return as array

# Supported formats: PNG, JPEG, BMP, TIFF, SVG

# High-resolution output
plotter = pv.Plotter(off_screen=True, window_size=[3840, 2160])
plotter.add_mesh(mesh, scalars='pressure')
plotter.screenshot('high_res.png')
```

### Save Animations

```python
# Save as GIF
plotter.open_gif('rotation.gif')
for angle in range(0, 360, 2):
    plotter.camera.azimuth(2)
    plotter.write_frame()
plotter.close()

# Save as MP4 (requires imageio-ffmpeg)
plotter.open_movie('rotation.mp4', framerate=30)
for angle in range(0, 360, 2):
    plotter.camera.azimuth(2)
    plotter.write_frame()
plotter.close()
```

### Save Meshes

```python
# Save in VTK formats
mesh.save('output.vtk')      # Legacy VTK
mesh.save('output.vtu')      # Unstructured grid
mesh.save('output.vtp')      # PolyData
mesh.save('output.vts')      # Structured grid
mesh.save('output.vti')      # Image data

# Save in other formats (requires meshio)
mesh.save('output.stl')      # STL
mesh.save('output.ply')      # PLY
mesh.save('output.obj')      # Wavefront OBJ
```

## Performance Optimization

### Decimation (Reduce Points)

```python
# Reduce mesh density
decimated = mesh.decimate(target_reduction=0.5)  # Keep 50% of points

# Decimate preserving boundaries
decimated = mesh.decimate(0.7, boundary_vertex_deletion=False)
```

### Level of Detail (LOD)

```python
# Enable automatic LOD for large meshes
plotter.enable_lod()  # Reduce quality during interaction

# Disable LOD
plotter.disable_lod()
```

### GPU Rendering

```python
# GPU-accelerated volume rendering
plotter.add_volume(mesh, mapper='gpu')  # vs 'cpu'

# Enable depth peeling for transparency
plotter.enable_depth_peeling(number_of_peels=4)
```

### Memory Management

```python
# Remove actor from scene
plotter.remove_actor(actor)

# Clear all actors
plotter.clear()

# Deep clean
plotter.deep_clean()

# Close and release resources
plotter.close()
```

## Common Patterns for CFD

### Load and Visualize CFD Results

```python
# Load VTK file from CFD solver
mesh = pv.read('cfd_results.vtu')

# Quick look at available fields
print("Point data:", mesh.point_data.keys())
print("Cell data:", mesh.cell_data.keys())

# Basic visualization
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars='Velocity', cmap='jet')
plotter.show()
```

### Velocity Magnitude

```python
# Calculate velocity magnitude if not available
if 'velocity' in mesh.point_data:
    velocity = mesh['velocity']
    mesh['speed'] = np.linalg.norm(velocity, axis=1)
```

### Compare Multiple Cases

```python
# Side-by-side comparison
mesh1 = pv.read('case1.vtu')
mesh2 = pv.read('case2.vtu')

plotter = pv.Plotter(shape=(1, 2))

plotter.subplot(0, 0)
plotter.add_mesh(mesh1, scalars='pressure', clim=[95000, 105000])
plotter.add_text('Case 1', font_size=12)

plotter.subplot(0, 1)
plotter.add_mesh(mesh2, scalars='pressure', clim=[95000, 105000])
plotter.add_text('Case 2', font_size=12)

plotter.link_views()
plotter.show()
```

### Extract Wall Surfaces

```python
# Extract surface mesh
surface = mesh.extract_surface()

# Get wall boundaries (if marked)
wall = mesh.threshold([wall_id, wall_id], scalars='boundary_id')
```

## Documentation Links

### Official Documentation

- **PyVista Documentation**: https://docs.pyvista.org/
- **API Reference**: https://docs.pyvista.org/api/index.html
- **Examples Gallery**: https://docs.pyvista.org/examples/index.html
- **User Guide**: https://docs.pyvista.org/user-guide/index.html

### VTK Resources

- **VTK Documentation**: https://vtk.org/documentation/
- **VTK File Formats**: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
- **VTK Examples**: https://kitware.github.io/vtk-examples/site/

### CFD-Specific Resources

- **ParaView (GUI alternative)**: https://www.paraview.org/
- **ParaView Guide**: https://docs.paraview.org/en/latest/
- **VTK in CFD**: https://www.openfoam.com/documentation/guides/latest/doc/guide-post-processing-vtk.html

### Tutorials and Examples

- **PyVista Tutorial**: https://tutorial.pyvista.org/
- **3D Visualization Course**: https://docs.pyvista.org/examples/index.html
- **CFD Post-Processing**: https://docs.pyvista.org/examples/00-load/README.html

### Related Packages

- **VTK**: https://vtk.org/
- **meshio**: https://github.com/nschloe/meshio (file I/O)
- **trimesh**: https://trimsh.org/ (mesh processing)
- **matplotlib**: https://matplotlib.org/ (2D plotting)

## Quick Command Reference

```python
# Import
import pyvista as pv
import numpy as np

# Create mesh
mesh = pv.ImageData(dimensions=(50, 50, 50))
mesh.spacing = (0.1, 0.1, 0.1)

# Add data
mesh['pressure'] = np.random.rand(mesh.n_points)
mesh['velocity'] = np.random.rand(mesh.n_points, 3)

# Plotter
p = pv.Plotter()

# Add mesh
p.add_mesh(mesh, scalars='pressure', cmap='coolwarm')

# Slice
s = mesh.slice(normal='z')
p.add_mesh(s, scalars='pressure')

# Contour
c = mesh.contour(isosurfaces=10, scalars='pressure')
p.add_mesh(c, opacity=0.5)

# Streamlines
st = mesh.streamlines(vectors='velocity', n_points=20)
p.add_mesh(st.tube(radius=0.005), color='white')

# Camera
p.camera_position = 'iso'
p.show_axes()

# Show
p.show()

# Save
p.screenshot('output.png')
mesh.save('data.vtu')
```

## Troubleshooting

### Common Issues

**Issue**: Mesh appears empty
```python
# Check mesh has points and cells
print(f"Points: {mesh.n_points}, Cells: {mesh.n_cells}")
# Check data exists
print(mesh.point_data.keys())
```

**Issue**: Colorbar not showing
```python
# Ensure scalars are specified and show_scalar_bar=True
plotter.add_mesh(mesh, scalars='pressure', show_scalar_bar=True)
```

**Issue**: Streamlines fail
```python
# Verify vector field exists and has correct shape
print(mesh['velocity'].shape)  # Should be (n_points, 3)
```

**Issue**: Slow rendering
```python
# Reduce mesh size
mesh_reduced = mesh.decimate(0.5)
# Enable LOD
plotter.enable_lod()
# Use smaller window
plotter = pv.Plotter(window_size=[800, 600])
```

## See Also

- `SKILL.md` - PyVista capabilities and CFD applications
- `examples.py` - Complete working examples for CFD visualization
