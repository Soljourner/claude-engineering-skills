"""
PyVista 3D Visualization Examples for CFD Applications

This file contains complete, working examples for creating 3D CFD visualizations
including velocity fields, pressure contours, streamlines, and pump geometry.

Run individual functions to see examples.
"""

import pyvista as pv
import numpy as np


def example_create_3d_mesh():
    """
    Create different types of 3D meshes for CFD domains.

    Demonstrates structured grids, unstructured grids, and uniform grids
    commonly used in CFD simulations.
    """
    print("Creating 3D meshes...")

    # Method 1: Structured Grid (curvilinear coordinates)
    # Perfect for pump volutes, blade passages, etc.
    x = np.linspace(0, 2, 30)
    y = np.linspace(0, 1, 20)
    z = np.linspace(0, 0.5, 15)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Apply curvilinear transformation (simulate volute expansion)
    theta = np.pi * X / 2
    R = 0.5 + 0.3 * X
    X_curved = R * np.cos(theta)
    Y_curved = R * np.sin(theta)
    Z_curved = Z

    structured_grid = pv.StructuredGrid(X_curved, Y_curved, Z_curved)

    # Method 2: Uniform Grid (regular Cartesian)
    # Good for simple rectangular domains
    uniform_grid = pv.ImageData(dimensions=(40, 30, 25))
    uniform_grid.spacing = (0.05, 0.05, 0.05)
    uniform_grid.origin = (0, 0, 0)

    # Visualize both meshes
    plotter = pv.Plotter(shape=(1, 2))

    # Left: Structured grid
    plotter.subplot(0, 0)
    plotter.add_mesh(structured_grid, show_edges=True, color='lightblue',
                    opacity=0.5, line_width=1)
    plotter.add_text('Structured Grid (Curvilinear)', font_size=12)
    plotter.camera_position = 'iso'

    # Right: Uniform grid
    plotter.subplot(0, 1)
    plotter.add_mesh(uniform_grid.outline(), color='black', line_width=2)
    plotter.add_mesh(uniform_grid, show_edges=True, color='lightgreen',
                    opacity=0.3, line_width=0.5)
    plotter.add_text('Uniform Grid (Cartesian)', font_size=12)
    plotter.camera_position = 'iso'

    plotter.link_views()
    plotter.show()

    print("Mesh creation example completed.")
    return structured_grid, uniform_grid


def example_velocity_field_visualization():
    """
    Visualize 3D velocity field with multiple techniques.

    Creates a pipe flow with swirl and visualizes using contours,
    glyphs, and color mapping.
    """
    print("Creating velocity field visualization...")

    # Create domain
    grid = pv.ImageData(dimensions=(50, 50, 80))
    grid.spacing = (0.02, 0.02, 0.015)
    grid.origin = (-0.5, -0.5, 0)

    # Get point coordinates
    x, y, z = grid.points.T

    # Create swirling pipe flow
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    # Velocity components
    # Axial velocity (parabolic profile)
    w = 3.0 * (1 - (r/0.4)**2) * (r <= 0.4)

    # Tangential velocity (forced vortex in core, free vortex outside)
    swirl_strength = 1.5
    v_theta = np.where(r <= 0.2,
                       swirl_strength * r,  # Solid body rotation
                       swirl_strength * 0.2**2 / (r + 0.01))  # Free vortex

    # Convert to Cartesian
    u = -v_theta * np.sin(theta)
    v = v_theta * np.cos(theta)

    # Add perturbations for turbulence
    u += 0.1 * np.random.randn(len(u)) * (r <= 0.4)
    v += 0.1 * np.random.randn(len(v)) * (r <= 0.4)
    w += 0.1 * np.random.randn(len(w)) * (r <= 0.4)

    # Add to grid
    grid['velocity'] = np.column_stack([u, v, w])
    grid['speed'] = np.linalg.norm(grid['velocity'], axis=1)

    # Mask outside pipe
    grid['speed'][r > 0.4] = 0

    # Create visualization
    plotter = pv.Plotter()

    # Add velocity magnitude on axial slices
    slice_positions = [0.2, 0.6, 1.0]
    for z_pos in slice_positions:
        slice_z = grid.slice(normal='z', origin=(0, 0, z_pos))
        plotter.add_mesh(slice_z, scalars='speed', cmap='jet',
                        opacity=0.9, show_edges=False,
                        clim=[0, 3.5])

    # Add pipe wall
    pipe = pv.Cylinder(radius=0.4, height=1.2, center=(0, 0, 0.6),
                      direction=(0, 0, 1), resolution=60)
    plotter.add_mesh(pipe, color='gray', opacity=0.2, show_edges=True)

    # Add velocity vectors on center plane
    sample_plane = grid.slice(normal='y', origin=(0, 0, 0.6))
    sample_plane = sample_plane.glyph(orient='velocity', scale='speed',
                                     factor=0.03, geom=pv.Arrow())
    plotter.add_mesh(sample_plane, cmap='plasma', show_scalar_bar=True,
                    scalar_bar_args={'title': 'Velocity Magnitude [m/s]',
                                   'vertical': True,
                                   'height': 0.7,
                                   'position_x': 0.85,
                                   'position_y': 0.15})

    plotter.add_text('3D Velocity Field - Swirling Pipe Flow', font_size=14)
    plotter.camera_position = [(2.5, 2.5, 1.5), (0, 0, 0.6), (0, 0, 1)]
    plotter.show_axes()
    plotter.show()

    print("Velocity field visualization completed.")
    return grid


def example_streamlines_in_pump():
    """
    Generate and visualize streamlines in a simplified pump geometry.

    Shows flow pathlines from inlet through impeller to outlet,
    demonstrating particle tracking in 3D flow.
    """
    print("Creating pump streamline visualization...")

    # Create pump domain (simplified cylindrical volute)
    grid = pv.ImageData(dimensions=(60, 60, 40))
    grid.spacing = (0.03, 0.03, 0.02)
    grid.origin = (-0.9, -0.9, -0.4)

    # Get coordinates
    x, y, z = grid.points.T
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    # Create velocity field (centrifugal pump pattern)
    # Radial component (outward flow)
    v_r = 2.0 * np.exp(-((r - 0.4) / 0.2)**2) * (r > 0.1)

    # Tangential component (rotation from impeller)
    omega = 100  # rad/s
    v_theta = omega * 0.3 * np.exp(-((r - 0.3) / 0.15)**2)

    # Axial component (flow along impeller)
    v_z = 1.5 * np.exp(-r**2 / 0.1) * (1 - np.abs(z) / 0.4)

    # Convert to Cartesian
    u = v_r * np.cos(theta) - v_theta * np.sin(theta)
    v = v_r * np.sin(theta) + v_theta * np.cos(theta)
    w = v_z

    grid['velocity'] = np.column_stack([u, v, w])
    grid['speed'] = np.linalg.norm(grid['velocity'], axis=1)

    # Create seed points at inlet (center)
    n_seeds = 40
    seed_radius = 0.15
    seed_theta = np.linspace(0, 2*np.pi, n_seeds, endpoint=False)
    seed_r = np.random.uniform(0, seed_radius, n_seeds)
    seed_x = seed_r * np.cos(seed_theta)
    seed_y = seed_r * np.sin(seed_theta)
    seed_z = np.zeros(n_seeds)
    seed_points = np.column_stack([seed_x, seed_y, seed_z])

    # Generate streamlines
    streamlines = grid.streamlines_from_source(
        pv.PolyData(seed_points),
        vectors='velocity',
        max_time=2.0,
        initial_step_length=0.01,
        integration_direction='forward'
    )

    # Create visualization
    plotter = pv.Plotter()

    # Add streamlines as tubes colored by velocity
    if streamlines.n_points > 0:
        tubes = streamlines.tube(radius=0.008)
        plotter.add_mesh(tubes, scalars='speed', cmap='jet',
                        show_scalar_bar=True, clim=[0, 8],
                        scalar_bar_args={'title': 'Velocity [m/s]',
                                       'vertical': True,
                                       'height': 0.7})

    # Add impeller (simplified as cylinder with blades)
    hub = pv.Cylinder(radius=0.12, height=0.3, center=(0, 0, 0),
                     direction=(0, 0, 1), resolution=40)
    plotter.add_mesh(hub, color='silver', metallic=0.5, roughness=0.3)

    # Add blade hints (simple rectangles)
    for i in range(6):
        angle = i * 60
        blade = pv.Plane(center=(0.2, 0, 0), direction=(0, 0, 1),
                        i_size=0.15, j_size=0.25)
        blade.rotate_z(angle, point=(0, 0, 0))
        plotter.add_mesh(blade, color='lightblue', opacity=0.6)

    # Add volute outline
    volute = pv.Cylinder(radius=0.8, height=0.8, center=(0, 0, 0),
                        direction=(0, 0, 1), resolution=60)
    plotter.add_mesh(volute.extract_feature_edges(), color='black',
                    line_width=3, opacity=0.5)

    # Add seed points
    plotter.add_mesh(pv.PolyData(seed_points), color='red',
                    point_size=10, render_points_as_spheres=True)

    plotter.add_text('Streamlines in Centrifugal Pump', font_size=14)
    plotter.camera_position = [(2.5, 2.5, 2.0), (0, 0, 0), (0, 0, 1)]
    plotter.show_axes()
    plotter.show()

    print("Streamline visualization completed.")
    return grid, streamlines


def example_scalar_contours():
    """
    Create pressure contour visualization with isosurfaces.

    Shows multiple pressure isosurfaces representing constant
    pressure levels in a 3D flow domain.
    """
    print("Creating scalar contour visualization...")

    # Create domain
    grid = pv.ImageData(dimensions=(60, 50, 40))
    grid.spacing = (0.03, 0.03, 0.025)
    grid.origin = (0, 0, 0)

    # Get coordinates
    x, y, z = grid.points.T

    # Create pressure field (flow around obstacle with stagnation)
    # Obstacle at center
    obstacle_x, obstacle_y, obstacle_z = 0.9, 0.75, 0.5
    r_obs = np.sqrt((x - obstacle_x)**2 + (y - obstacle_y)**2 + (z - obstacle_z)**2)

    # Stagnation pressure near obstacle, lower pressure in wake
    P_inlet = 150000  # Pa
    P_dynamic = 5000  # Pa

    # Stagnation point
    pressure = P_inlet + P_dynamic * (1 - np.tanh(5 * (x - obstacle_x)))

    # Low pressure in wake
    wake_factor = np.exp(-((x - obstacle_x - 0.3)**2 + (y - obstacle_y)**2) / 0.1)
    pressure -= P_dynamic * 1.5 * wake_factor * (x > obstacle_x)

    # High pressure at stagnation
    stag_factor = np.exp(-r_obs**2 / 0.02)
    pressure += P_dynamic * 2 * stag_factor

    grid['pressure'] = pressure

    # Create velocity field for context
    u = 10.0 * (1 - np.exp(-r_obs**2 / 0.1))
    v = -2.0 * (y - obstacle_y) * np.exp(-r_obs**2 / 0.05)
    w = -2.0 * (z - obstacle_z) * np.exp(-r_obs**2 / 0.05)
    grid['velocity'] = np.column_stack([u, v, w])

    # Create visualization
    plotter = pv.Plotter()

    # Create multiple isosurfaces
    n_contours = 12
    contours = grid.contour(isosurfaces=n_contours, scalars='pressure')

    plotter.add_mesh(contours, scalars='pressure', cmap='coolwarm',
                    opacity=0.6, show_scalar_bar=True,
                    scalar_bar_args={'title': 'Pressure [Pa]',
                                   'vertical': True,
                                   'height': 0.7,
                                   'n_labels': 5,
                                   'fmt': '%.0f'})

    # Add obstacle (sphere)
    obstacle = pv.Sphere(radius=0.15, center=(obstacle_x, obstacle_y, obstacle_z))
    plotter.add_mesh(obstacle, color='gray', opacity=0.9)

    # Add domain outline
    plotter.add_mesh(grid.outline(), color='black', line_width=2)

    # Add a few streamlines for flow context
    seed = pv.Line((0.1, 0.75, 0.5), (0.3, 0.75, 0.5), resolution=5)
    streamlines = grid.streamlines_from_source(seed, vectors='velocity',
                                              max_time=0.5)
    if streamlines.n_points > 0:
        plotter.add_mesh(streamlines.tube(radius=0.005), color='white',
                        opacity=0.7)

    plotter.add_text('Pressure Isosurfaces - Flow Around Obstacle', font_size=14)
    plotter.camera_position = [(4.0, 3.5, 2.5), (0.9, 0.75, 0.5), (0, 0, 1)]
    plotter.show_axes()
    plotter.show()

    print("Scalar contour visualization completed.")
    return grid


def example_interactive_plots():
    """
    Create interactive visualization with sliders and controls.

    Demonstrates interactive parameter adjustment for exploring
    CFD results with different visualization settings.
    """
    print("Creating interactive visualization...")

    # Create domain
    grid = pv.ImageData(dimensions=(50, 50, 50))
    grid.spacing = (0.04, 0.04, 0.04)
    grid.origin = (-1, -1, -1)

    # Get coordinates
    x, y, z = grid.points.T
    r = np.sqrt(x**2 + y**2 + z**2)

    # Create turbulent kinetic energy field
    k_turb = 0.5 * np.exp(-r**2 / 0.5) * (1 + 0.3 * np.sin(5*x) * np.cos(5*y))
    grid['TKE'] = np.maximum(k_turb, 0)

    # Create dissipation rate field
    epsilon = 0.1 * k_turb**(1.5) / 0.1
    grid['epsilon'] = np.maximum(epsilon, 0)

    # Create velocity magnitude
    grid['speed'] = 5.0 * np.exp(-r**2 / 0.8)

    # Create plotter
    plotter = pv.Plotter()

    # Initial visualization
    slice_mesh = grid.slice_orthogonal(x=0, y=0, z=0)
    actor = plotter.add_mesh(slice_mesh, scalars='TKE', cmap='hot',
                            show_scalar_bar=True, clim=[0, 0.5],
                            scalar_bar_args={'title': 'TKE [m²/s²]'})

    plotter.add_text('Interactive CFD Visualization\n(Rotate/Zoom with mouse)',
                    font_size=12, position='upper_left')

    # Add axes and outline
    plotter.add_mesh(grid.outline(), color='black', line_width=2)
    plotter.show_axes()

    # Show with interaction enabled
    plotter.show()

    print("Interactive visualization completed.")
    print("Note: Use mouse to rotate, zoom. Close window to continue.")
    return grid


def example_volume_rendering():
    """
    Demonstrate volume rendering for 3D scalar fields.

    Uses opacity transfer functions to show internal structure
    of temperature or density distributions.
    """
    print("Creating volume rendering visualization...")

    # Create domain
    grid = pv.ImageData(dimensions=(80, 80, 80))
    grid.spacing = (0.025, 0.025, 0.025)
    grid.origin = (-1, -1, -1)

    # Get coordinates
    x, y, z = grid.points.T
    r = np.sqrt(x**2 + y**2 + z**2)

    # Create temperature field (hot core, cool exterior)
    T_core = 500  # K
    T_ambient = 300  # K
    temperature = T_ambient + (T_core - T_ambient) * np.exp(-r**2 / 0.3)

    # Add turbulent fluctuations
    temperature += 20 * np.sin(8*x) * np.cos(8*y) * np.sin(8*z) * np.exp(-r**2 / 0.5)

    grid['temperature'] = temperature

    # Create plotter
    plotter = pv.Plotter()

    # Volume rendering with opacity mapping
    plotter.add_volume(grid, scalars='temperature',
                      cmap='hot',
                      opacity='sigmoid',  # Sigmoid opacity transfer function
                      scalar_bar_args={'title': 'Temperature [K]',
                                     'vertical': True,
                                     'height': 0.7})

    # Add outline for reference
    plotter.add_mesh(grid.outline(), color='white', line_width=2)

    plotter.add_text('Volume Rendering - Temperature Distribution', font_size=14)
    plotter.camera_position = [(4, 4, 4), (0, 0, 0), (0, 0, 1)]
    plotter.show_axes()
    plotter.show()

    print("Volume rendering visualization completed.")
    return grid


def example_slicing_and_clipping():
    """
    Demonstrate slicing and clipping techniques.

    Shows how to extract 2D slices from 3D data and clip
    regions to reveal internal flow structures.
    """
    print("Creating slicing and clipping visualization...")

    # Create domain
    grid = pv.ImageData(dimensions=(60, 60, 60))
    grid.spacing = (0.03, 0.03, 0.03)
    grid.origin = (0, 0, 0)

    # Get coordinates
    x, y, z = grid.points.T

    # Create vorticity magnitude field
    vorticity = 50 * np.exp(-((x-0.9)**2 + (y-0.9)**2 + (z-0.9)**2) / 0.2)
    vorticity += 30 * np.exp(-((x-0.6)**2 + (y-0.6)**2 + (z-0.6)**2) / 0.15)
    grid['vorticity'] = vorticity

    # Create velocity for streamlines
    u = (y - 0.9) * vorticity / 50
    v = -(x - 0.9) * vorticity / 50
    w = 2.0 * np.ones_like(x)
    grid['velocity'] = np.column_stack([u, v, w])

    # Create visualization with multiple viewports
    plotter = pv.Plotter(shape=(2, 2))

    # Subplot 1: Orthogonal slices
    plotter.subplot(0, 0)
    slices = grid.slice_orthogonal(x=0.9, y=0.9, z=0.9)
    plotter.add_mesh(slices, scalars='vorticity', cmap='plasma',
                    show_scalar_bar=True, clim=[0, 50],
                    scalar_bar_args={'title': 'Vorticity [1/s]'})
    plotter.add_text('Orthogonal Slices', font_size=10)
    plotter.camera_position = 'iso'

    # Subplot 2: Multiple parallel slices
    plotter.subplot(0, 1)
    for z_val in [0.3, 0.6, 0.9, 1.2, 1.5]:
        slice_z = grid.slice(normal='z', origin=(0.9, 0.9, z_val))
        plotter.add_mesh(slice_z, scalars='vorticity', cmap='plasma',
                        opacity=0.7, show_scalar_bar=False, clim=[0, 50])
    plotter.add_mesh(grid.outline(), color='black', line_width=2)
    plotter.add_text('Parallel Z-Slices', font_size=10)
    plotter.camera_position = 'xz'

    # Subplot 3: Clipped domain
    plotter.subplot(1, 0)
    clipped = grid.clip(normal=(1, 1, 1), origin=(0.9, 0.9, 0.9))
    plotter.add_mesh(clipped, scalars='vorticity', cmap='plasma',
                    show_edges=True, edge_color='white',
                    opacity=0.9, show_scalar_bar=False, clim=[0, 50])
    plotter.add_text('Clipped Domain', font_size=10)
    plotter.camera_position = 'iso'

    # Subplot 4: Threshold (show only high vorticity regions)
    plotter.subplot(1, 1)
    high_vorticity = grid.threshold(value=25, scalars='vorticity')
    plotter.add_mesh(high_vorticity, scalars='vorticity', cmap='hot',
                    opacity=0.8, show_scalar_bar=False, clim=[25, 50])
    plotter.add_mesh(grid.outline(), color='black', line_width=2)
    plotter.add_text('Threshold (Vorticity > 25)', font_size=10)
    plotter.camera_position = 'iso'

    plotter.link_views()
    plotter.show()

    print("Slicing and clipping visualization completed.")
    return grid


def example_comprehensive_cfd_visualization():
    """
    Comprehensive CFD visualization combining multiple techniques.

    Demonstrates best practices for presenting complex 3D flow
    results with streamlines, contours, and slices.
    """
    print("Creating comprehensive CFD visualization...")

    # Create realistic pump domain
    grid = pv.ImageData(dimensions=(70, 70, 50))
    grid.spacing = (0.025, 0.025, 0.02)
    grid.origin = (-0.875, -0.875, -0.5)

    # Get coordinates
    x, y, z = grid.points.T
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)

    # Create realistic pump velocity field
    # Radial velocity
    v_r = 5.0 * (r / 0.6) * np.exp(-((r - 0.4) / 0.2)**2) * (r > 0.15)

    # Tangential velocity (impeller rotation)
    omega = 150  # rad/s
    v_theta = omega * r * np.exp(-((r - 0.3) / 0.2)**2) * (1 - np.abs(z) / 0.5)

    # Axial velocity
    v_z = 3.0 * np.exp(-r**2 / 0.2) * np.cos(np.pi * z / 1.0)

    # Convert to Cartesian
    u = v_r * np.cos(theta) - v_theta * np.sin(theta)
    v = v_r * np.sin(theta) + v_theta * np.cos(theta)
    w = v_z

    grid['velocity'] = np.column_stack([u, v, w])
    grid['speed'] = np.linalg.norm(grid['velocity'], axis=1)

    # Pressure field
    P_ref = 101325
    pressure = P_ref + 8000 * (1 - r**2 / 0.64) - 2000 * np.abs(z) / 0.5
    grid['pressure'] = pressure

    # Turbulent kinetic energy
    tke = 0.5 * grid['speed']**2 * 0.05  # 5% turbulence intensity
    grid['TKE'] = tke

    # Create comprehensive visualization
    plotter = pv.Plotter(window_size=[1600, 1200])

    # 1. Central slice showing velocity magnitude
    center_slice = grid.slice(normal='z', origin=(0, 0, 0))
    plotter.add_mesh(center_slice, scalars='speed', cmap='jet',
                    opacity=0.95, show_scalar_bar=True,
                    scalar_bar_args={'title': 'Velocity Magnitude [m/s]',
                                   'vertical': True,
                                   'height': 0.6,
                                   'position_x': 0.85,
                                   'position_y': 0.2,
                                   'n_labels': 5,
                                   'fmt': '%.1f'})

    # 2. Streamlines from inlet
    seed_circle = pv.Disc(center=(0, 0, -0.4), inner=0, outer=0.25,
                         normal=(0, 0, 1), r_res=2, c_res=20)
    streamlines = grid.streamlines_from_source(
        seed_circle,
        vectors='velocity',
        max_time=0.3,
        initial_step_length=0.01
    )
    if streamlines.n_points > 0:
        tubes = streamlines.tube(radius=0.008)
        plotter.add_mesh(tubes, color='white', opacity=0.8)

    # 3. Pressure isosurface (high pressure region)
    high_p = grid.threshold(value=105000, scalars='pressure')
    plotter.add_mesh(high_p, color='red', opacity=0.2)

    # 4. Low pressure isosurface
    low_p = grid.threshold(value=[95000, 98000], scalars='pressure')
    plotter.add_mesh(low_p, color='blue', opacity=0.15)

    # 5. Impeller geometry (simplified)
    hub = pv.Cylinder(radius=0.15, height=0.4, center=(0, 0, 0),
                     direction=(0, 0, 1), resolution=40)
    plotter.add_mesh(hub, color='silver', metallic=0.6, roughness=0.4,
                    specular=0.5)

    # Add blades
    for i in range(6):
        angle = i * 60
        blade = pv.Plane(center=(0.25, 0, 0), direction=(0, 0, 1),
                        i_size=0.2, j_size=0.35)
        blade.rotate_z(angle, point=(0, 0, 0))
        plotter.add_mesh(blade, color='lightblue', opacity=0.7, metallic=0.3)

    # 6. Domain outline
    plotter.add_mesh(grid.outline(), color='black', line_width=3)

    # 7. Add coordinate axes
    plotter.show_axes()

    # Add title and labels
    plotter.add_text('Comprehensive Pump CFD Visualization\n' +
                    'Velocity field, Streamlines, Pressure isosurfaces',
                    position='upper_edge', font_size=12, color='black')

    # Set camera position
    plotter.camera_position = [(2.5, 2.5, 1.5), (0, 0, 0), (0, 0, 1)]

    # Enable better lighting
    plotter.enable_anti_aliasing()

    # Show
    plotter.show()

    print("Comprehensive CFD visualization completed.")
    return grid


def main():
    """
    Run all PyVista CFD visualization examples.
    """
    print("=" * 80)
    print("PyVista 3D Visualization Examples for CFD Applications")
    print("=" * 80)
    print()

    print("Example 1: Creating 3D Meshes")
    print("-" * 80)
    example_create_3d_mesh()
    print()

    print("Example 2: Velocity Field Visualization")
    print("-" * 80)
    example_velocity_field_visualization()
    print()

    print("Example 3: Streamlines in Pump")
    print("-" * 80)
    example_streamlines_in_pump()
    print()

    print("Example 4: Scalar Contours (Pressure)")
    print("-" * 80)
    example_scalar_contours()
    print()

    print("Example 5: Interactive Visualization")
    print("-" * 80)
    example_interactive_plots()
    print()

    print("Example 6: Volume Rendering")
    print("-" * 80)
    example_volume_rendering()
    print()

    print("Example 7: Slicing and Clipping")
    print("-" * 80)
    example_slicing_and_clipping()
    print()

    print("Example 8: Comprehensive CFD Visualization")
    print("-" * 80)
    example_comprehensive_cfd_visualization()
    print()

    print("=" * 80)
    print("All PyVista examples completed!")
    print("=" * 80)


if __name__ == "__main__":
    # Set plot theme
    pv.set_plot_theme('document')

    # Run all examples
    main()

    # Or run individual examples:
    # example_create_3d_mesh()
    # example_velocity_field_visualization()
    # example_streamlines_in_pump()
    # example_scalar_contours()
    # example_interactive_plots()
    # example_volume_rendering()
    # example_slicing_and_clipping()
    # example_comprehensive_cfd_visualization()
