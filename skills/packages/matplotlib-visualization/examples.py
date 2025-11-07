"""
Matplotlib Visualization Examples for Engineering Applications

This file contains complete, working examples for creating professional
engineering plots including pump curves, CFD results, and publication-quality figures.

Run individual functions to see examples.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D


def example_pump_performance_curves():
    """
    Complete pump performance curve with H-Q, P-Q, and η-Q subplots.

    Creates a three-panel figure showing head, power, and efficiency
    curves versus flow rate - standard presentation for centrifugal pumps.
    """
    # Generate pump performance data
    Q = np.linspace(0, 120, 50)  # Flow rate [m³/h]

    # Head curve: H = H0 - k*Q²
    H0 = 80  # Shutoff head [m]
    k_h = 0.005
    H = H0 - k_h * Q**2

    # Power curve: P = P0 + a*Q + b*Q²
    P0 = 5  # No-load power [kW]
    a = 0.08
    b = 0.0005
    P = P0 + a * Q + b * Q**2

    # Efficiency: η = (ρ*g*Q*H) / (P*1000) * 100
    # Simplified: η = Q*H / (367*P) * 100
    eta = (Q * H) / (367 * P) * 100
    eta = np.nan_to_num(eta)  # Handle division by zero at Q=0

    # Create figure with three vertically stacked subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # Plot 1: Head vs Flow Rate
    ax1.plot(Q, H, 'b-', linewidth=2.5, label='Head curve')
    ax1.axhline(y=H0, color='b', linestyle='--', alpha=0.3, linewidth=1)
    ax1.set_ylabel('Head $H$ [m]', fontsize=13, fontweight='bold')
    ax1.set_title('Centrifugal Pump Performance Curves',
                  fontsize=15, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax1.legend(loc='upper right', fontsize=11)
    ax1.set_ylim([0, H0 * 1.1])

    # Plot 2: Power vs Flow Rate
    ax2.plot(Q, P, 'r-', linewidth=2.5, label='Power curve')
    ax2.set_ylabel('Power $P$ [kW]', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax2.legend(loc='upper left', fontsize=11)
    ax2.set_ylim([0, P.max() * 1.1])

    # Plot 3: Efficiency vs Flow Rate
    ax3.plot(Q, eta, 'g-', linewidth=2.5, label='Efficiency curve')

    # Mark best efficiency point (BEP)
    bep_idx = np.argmax(eta)
    ax3.plot(Q[bep_idx], eta[bep_idx], 'ro', markersize=10,
             label=f'BEP: {Q[bep_idx]:.1f} m³/h, {eta[bep_idx]:.1f}%')
    ax3.axvline(x=Q[bep_idx], color='r', linestyle='--', alpha=0.3, linewidth=1)

    ax3.set_xlabel('Flow Rate $Q$ [m³/h]', fontsize=13, fontweight='bold')
    ax3.set_ylabel('Efficiency $\eta$ [%]', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax3.legend(loc='upper right', fontsize=11)
    ax3.set_ylim([0, 100])

    plt.tight_layout()
    plt.savefig('pump_performance_curves.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: pump_performance_curves.png")
    plt.show()


def example_pump_and_system_curves():
    """
    Pump curve overlaid with system curve to show operating point.

    Demonstrates how to find the intersection of pump and system curves
    to determine the actual operating conditions.
    """
    # Pump curve
    Q = np.linspace(0, 120, 100)
    H_pump = 80 - 0.005 * Q**2

    # System curves for different scenarios
    H_static = 20  # Static head [m]
    k1 = 0.003  # Original system
    k2 = 0.005  # Higher resistance (fouled pipes)
    k3 = 0.002  # Lower resistance (new pipes)

    H_system1 = H_static + k1 * Q**2
    H_system2 = H_static + k2 * Q**2
    H_system3 = H_static + k3 * Q**2

    # Find operating points
    def find_operating_point(H_pump, H_system):
        idx = np.argmin(np.abs(H_pump - H_system))
        return idx

    op1 = find_operating_point(H_pump, H_system1)
    op2 = find_operating_point(H_pump, H_system2)
    op3 = find_operating_point(H_pump, H_system3)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))

    # Plot curves
    ax.plot(Q, H_pump, 'b-', linewidth=3, label='Pump curve', zorder=5)
    ax.plot(Q, H_system1, 'r--', linewidth=2, label='Original system')
    ax.plot(Q, H_system2, 'orange', linewidth=2, linestyle='--', label='Fouled system')
    ax.plot(Q, H_system3, 'g--', linewidth=2, label='New system')

    # Mark operating points
    ax.plot(Q[op1], H_pump[op1], 'ro', markersize=12, zorder=10,
            label=f'Op. point 1: {Q[op1]:.1f} m³/h, {H_pump[op1]:.1f} m')
    ax.plot(Q[op2], H_pump[op2], 'o', color='orange', markersize=12, zorder=10,
            label=f'Op. point 2: {Q[op2]:.1f} m³/h, {H_pump[op2]:.1f} m')
    ax.plot(Q[op3], H_pump[op3], 'go', markersize=12, zorder=10,
            label=f'Op. point 3: {Q[op3]:.1f} m³/h, {H_pump[op3]:.1f} m')

    # Add vertical lines from operating points
    ax.axvline(x=Q[op1], color='r', linestyle=':', alpha=0.3, linewidth=1)
    ax.axvline(x=Q[op2], color='orange', linestyle=':', alpha=0.3, linewidth=1)
    ax.axvline(x=Q[op3], color='g', linestyle=':', alpha=0.3, linewidth=1)

    # Labels and styling
    ax.set_xlabel('Flow Rate $Q$ [m³/h]', fontsize=13, fontweight='bold')
    ax.set_ylabel('Head $H$ [m]', fontsize=13, fontweight='bold')
    ax.set_title('Pump and System Curves - Operating Point Analysis',
                 fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=10, framealpha=0.9)
    ax.set_xlim([0, 120])
    ax.set_ylim([0, 90])

    # Add text annotation
    ax.text(5, 75, r'$H_{system} = H_{static} + k \cdot Q^2$',
            fontsize=12, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

    plt.tight_layout()
    plt.savefig('pump_system_curves.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: pump_system_curves.png")
    plt.show()


def example_velocity_contour():
    """
    Velocity contour plot simulating CFD results.

    Shows flow around an obstacle with filled contours, contour lines,
    and a colorbar - typical visualization for CFD velocity fields.
    """
    # Create fine mesh grid
    x = np.linspace(0, 10, 150)
    y = np.linspace(0, 5, 75)
    X, Y = np.meshgrid(x, y)

    # Simulate velocity field (flow around circular obstacle at x=3, y=2.5)
    obstacle_x, obstacle_y = 3.0, 2.5
    obstacle_radius = 0.8

    # Distance from obstacle center
    r = np.sqrt((X - obstacle_x)**2 + (Y - obstacle_y)**2)

    # U-velocity component (horizontal)
    U = 2.0 * np.ones_like(X)
    U[r < obstacle_radius * 1.2] *= 0.1  # Reduced velocity near obstacle
    U += -1.5 * np.exp(-((X-obstacle_x)**2 + (Y-obstacle_y)**2) / 1.5)

    # V-velocity component (vertical)
    V = 0.5 * np.sin(2 * np.pi * X / 10) * np.exp(-((X-obstacle_x)**2 + (Y-obstacle_y)**2) / 2)

    # Total velocity magnitude
    velocity = np.sqrt(U**2 + V**2)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 7))

    # Filled contours
    levels = np.linspace(0, velocity.max(), 25)
    contourf = ax.contourf(X, Y, velocity, levels=levels, cmap='jet', extend='both')

    # Add colorbar
    cbar = plt.colorbar(contourf, ax=ax, label='Velocity magnitude [m/s]',
                        orientation='vertical', pad=0.02)
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label('Velocity magnitude [m/s]', fontsize=12, fontweight='bold')

    # Add contour lines
    contour = ax.contour(X, Y, velocity, levels=12, colors='k',
                         linewidths=0.5, alpha=0.4)
    ax.clabel(contour, inline=True, fontsize=8, fmt='%.2f')

    # Draw obstacle
    circle = plt.Circle((obstacle_x, obstacle_y), obstacle_radius,
                        color='white', ec='black', linewidth=2, zorder=10)
    ax.add_patch(circle)

    # Add streamlines
    skip = (slice(None, None, 3), slice(None, None, 3))
    ax.streamplot(X[skip], Y[skip], U[skip], V[skip],
                  color='white', linewidth=0.5, density=1,
                  arrowsize=1, alpha=0.5)

    # Labels and styling
    ax.set_xlabel('x [m]', fontsize=13, fontweight='bold')
    ax.set_ylabel('y [m]', fontsize=13, fontweight='bold')
    ax.set_title('Velocity Contour Plot - Flow Around Obstacle',
                 fontsize=15, fontweight='bold')
    ax.set_aspect('equal')
    ax.set_xlim([0, 10])
    ax.set_ylim([0, 5])

    plt.tight_layout()
    plt.savefig('velocity_contour.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: velocity_contour.png")
    plt.show()


def example_vector_field():
    """
    Vector field visualization with quiver plot.

    Shows velocity vectors colored by magnitude - useful for understanding
    flow direction and local velocity in CFD results.
    """
    # Create coarser grid for vectors (too many arrows gets cluttered)
    x = np.linspace(0, 10, 25)
    y = np.linspace(0, 5, 13)
    X, Y = np.meshgrid(x, y)

    # Velocity field - vortex pattern
    cx, cy = 5, 2.5  # Center of vortex
    dx = X - cx
    dy = Y - cy
    r = np.sqrt(dx**2 + dy**2)

    # Tangential velocity (rotational flow)
    U = -dy / (r + 0.5) + 1.5
    V = dx / (r + 0.5)

    # Velocity magnitude
    velocity = np.sqrt(U**2 + V**2)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 7))

    # Quiver plot (arrows)
    quiver = ax.quiver(X, Y, U, V, velocity,
                       cmap='viridis',
                       scale=25,
                       width=0.004,
                       headwidth=4,
                       headlength=5,
                       alpha=0.9)

    # Add colorbar
    cbar = plt.colorbar(quiver, ax=ax, label='Velocity magnitude [m/s]',
                        orientation='vertical', pad=0.02)
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label('Velocity magnitude [m/s]', fontsize=12, fontweight='bold')

    # Add contour background for context
    x_fine = np.linspace(0, 10, 100)
    y_fine = np.linspace(0, 5, 50)
    X_fine, Y_fine = np.meshgrid(x_fine, y_fine)
    dx_fine = X_fine - cx
    dy_fine = Y_fine - cy
    r_fine = np.sqrt(dx_fine**2 + dy_fine**2)
    U_fine = -dy_fine / (r_fine + 0.5) + 1.5
    V_fine = dx_fine / (r_fine + 0.5)
    velocity_fine = np.sqrt(U_fine**2 + V_fine**2)

    ax.contour(X_fine, Y_fine, velocity_fine, levels=15,
               colors='gray', linewidths=0.5, alpha=0.3)

    # Labels and styling
    ax.set_xlabel('x [m]', fontsize=13, fontweight='bold')
    ax.set_ylabel('y [m]', fontsize=13, fontweight='bold')
    ax.set_title('Velocity Vector Field - Vortex Flow Pattern',
                 fontsize=15, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)

    # Add reference vector
    ax.quiverkey(quiver, 0.85, 0.95, 2, '2 m/s',
                 labelpos='E', coordinates='axes',
                 color='black', fontproperties={'size': 11, 'weight': 'bold'})

    plt.tight_layout()
    plt.savefig('velocity_vector_field.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: velocity_vector_field.png")
    plt.show()


def example_3d_surface():
    """
    3D surface plot for pressure or temperature distributions.

    Creates an interactive 3D visualization suitable for displaying
    spatial distributions of scalar quantities.
    """
    # Create mesh
    x = np.linspace(-5, 5, 60)
    y = np.linspace(-5, 5, 60)
    X, Y = np.meshgrid(x, y)

    # Pressure distribution - example: wave interference pattern
    r1 = np.sqrt((X + 2)**2 + (Y + 1)**2)
    r2 = np.sqrt((X - 2)**2 + (Y - 1)**2)
    Z = 100 + 10 * np.cos(r1) + 10 * np.cos(r2) + 5 * np.exp(-0.05 * (X**2 + Y**2))

    # Create figure with 3D axes
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot
    surf = ax.plot_surface(X, Y, Z, cmap='coolwarm',
                           linewidth=0, antialiased=True,
                           alpha=0.9, edgecolor='none',
                           vmin=Z.min(), vmax=Z.max())

    # Optional: Add contour projection on bottom
    ax.contour(X, Y, Z, levels=15, offset=Z.min()-5,
               cmap='coolwarm', linewidths=1, alpha=0.5)

    # Customize axes
    ax.set_xlabel('x [m]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_ylabel('y [m]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_zlabel('Pressure [kPa]', fontsize=12, fontweight='bold', labelpad=10)
    ax.set_title('3D Pressure Distribution', fontsize=15, fontweight='bold', pad=20)

    # Set z-axis limits
    ax.set_zlim([Z.min()-5, Z.max()+5])

    # Add colorbar
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.1,
                        label='Pressure [kPa]')
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('Pressure [kPa]', fontsize=11, fontweight='bold')

    # Set viewing angle
    ax.view_init(elev=25, azim=135)

    # Improve grid
    ax.grid(True, alpha=0.3)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    plt.tight_layout()
    plt.savefig('pressure_3d_surface.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("Saved: pressure_3d_surface.png")
    plt.show()


def example_publication_quality():
    """
    Publication-quality multi-panel figure with all styling applied.

    Demonstrates best practices: proper sizing, fonts, colors, labels,
    legends, and layout for submission to technical journals.
    """
    # Set publication-quality defaults
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['grid.alpha'] = 0.3
    plt.rcParams['lines.linewidth'] = 1.8

    # Professional color palette
    colors = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E', '#77AC30']

    # Create figure for journal (double column width = 7 inches)
    fig = plt.figure(figsize=(7, 8))

    # Use GridSpec for flexible layout
    gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.35,
                          left=0.1, right=0.95, top=0.95, bottom=0.08)

    # Panel (a): Pump curves comparison
    ax1 = fig.add_subplot(gs[0, :])
    Q = np.linspace(0, 100, 50)
    for i, rpm in enumerate([1450, 1750, 2000]):
        H = 50 * (rpm/1750)**2 - 0.003 * Q**2 * (rpm/1750)**2
        ax1.plot(Q, H, color=colors[i], linewidth=2,
                label=f'{rpm} RPM', marker='o', markevery=8, markersize=5)

    ax1.set_xlabel('Flow rate, $Q$ [m$^3$/h]', fontsize=11)
    ax1.set_ylabel('Head, $H$ [m]', fontsize=11)
    ax1.set_title('(a) Effect of rotational speed on pump performance',
                  fontsize=11, fontweight='bold', loc='left')
    ax1.legend(frameon=True, fontsize=9, loc='upper right')
    ax1.grid(True, linestyle='-', linewidth=0.5, alpha=0.3)
    ax1.set_xlim([0, 100])

    # Panel (b): Efficiency comparison
    ax2 = fig.add_subplot(gs[1, 0])
    Q = np.linspace(0, 100, 50)
    for i, label in enumerate(['Design A', 'Design B', 'Design C']):
        eta = 85 * np.exp(-((Q - 50 - i*10)**2) / 800)
        ax2.plot(Q, eta, color=colors[i], linewidth=2, label=label)

    ax2.set_xlabel('Flow rate, $Q$ [m$^3$/h]', fontsize=11)
    ax2.set_ylabel('Efficiency, $\eta$ [%]', fontsize=11)
    ax2.set_title('(b) Efficiency comparison', fontsize=11, fontweight='bold', loc='left')
    ax2.legend(frameon=True, fontsize=9)
    ax2.grid(True, linestyle='-', linewidth=0.5, alpha=0.3)
    ax2.set_xlim([0, 100])
    ax2.set_ylim([0, 100])

    # Panel (c): Power consumption
    ax3 = fig.add_subplot(gs[1, 1])
    designs = ['A', 'B', 'C', 'D']
    power = [45, 52, 41, 48]
    bars = ax3.bar(designs, power, color=colors[:4], edgecolor='black', linewidth=1.2)

    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.0f}', ha='center', va='bottom', fontsize=9)

    ax3.set_xlabel('Design variant', fontsize=11)
    ax3.set_ylabel('Power, $P$ [kW]', fontsize=11)
    ax3.set_title('(c) Power consumption', fontsize=11, fontweight='bold', loc='left')
    ax3.grid(True, axis='y', linestyle='-', linewidth=0.5, alpha=0.3)
    ax3.set_ylim([0, 60])

    # Panel (d): Scatter plot with correlation
    ax4 = fig.add_subplot(gs[2, :])
    np.random.seed(42)
    flow_exp = np.random.rand(30) * 100
    flow_sim = flow_exp + np.random.randn(30) * 5

    ax4.scatter(flow_exp, flow_sim, s=50, alpha=0.6,
               color=colors[0], edgecolors='black', linewidth=0.5)

    # Add perfect correlation line
    ax4.plot([0, 100], [0, 100], 'k--', linewidth=1.5,
            alpha=0.5, label='Perfect correlation')

    # Add linear fit
    z = np.polyfit(flow_exp, flow_sim, 1)
    p = np.poly1d(z)
    ax4.plot(flow_exp, p(flow_exp), color=colors[1], linewidth=2,
            label=f'Linear fit: $y={z[0]:.2f}x{z[1]:+.2f}$')

    # Calculate R²
    residuals = flow_sim - p(flow_exp)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((flow_sim - np.mean(flow_sim))**2)
    r_squared = 1 - (ss_res / ss_tot)

    ax4.text(0.05, 0.95, f'$R^2 = {r_squared:.3f}$',
            transform=ax4.transAxes, fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax4.set_xlabel('Experimental flow rate [m$^3$/h]', fontsize=11)
    ax4.set_ylabel('Simulated flow rate [m$^3$/h]', fontsize=11)
    ax4.set_title('(d) Validation: experimental vs. numerical results',
                  fontsize=11, fontweight='bold', loc='left')
    ax4.legend(frameon=True, fontsize=9, loc='lower right')
    ax4.grid(True, linestyle='-', linewidth=0.5, alpha=0.3)
    ax4.set_xlim([0, 100])
    ax4.set_ylim([0, 100])
    ax4.set_aspect('equal')

    # Save in multiple formats
    plt.savefig('publication_figure.png', dpi=600, bbox_inches='tight', facecolor='white')
    plt.savefig('publication_figure.pdf', bbox_inches='tight', facecolor='white')
    print("Saved: publication_figure.png and publication_figure.pdf")
    plt.show()

    # Reset rcParams to defaults
    plt.rcParams.update(plt.rcParamsDefault)


def main():
    """
    Run all examples.
    """
    print("=" * 70)
    print("Matplotlib Visualization Examples for Engineering")
    print("=" * 70)
    print()

    print("Example 1: Pump Performance Curves")
    print("-" * 70)
    example_pump_performance_curves()
    print()

    print("Example 2: Pump and System Curves")
    print("-" * 70)
    example_pump_and_system_curves()
    print()

    print("Example 3: Velocity Contour Plot")
    print("-" * 70)
    example_velocity_contour()
    print()

    print("Example 4: Vector Field Visualization")
    print("-" * 70)
    example_vector_field()
    print()

    print("Example 5: 3D Surface Plot")
    print("-" * 70)
    example_3d_surface()
    print()

    print("Example 6: Publication-Quality Figure")
    print("-" * 70)
    example_publication_quality()
    print()

    print("=" * 70)
    print("All examples completed!")
    print("=" * 70)


if __name__ == "__main__":
    # Run all examples
    main()

    # Or run individual examples:
    # example_pump_performance_curves()
    # example_pump_and_system_curves()
    # example_velocity_contour()
    # example_vector_field()
    # example_3d_surface()
    # example_publication_quality()
