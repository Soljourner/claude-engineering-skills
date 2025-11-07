#!/usr/bin/env python3
"""
ANSYS Fluent - Centrifugal Pump CFD Setup

This script automates the setup of a centrifugal pump CFD simulation
including rotating reference frame, sliding mesh, or frozen rotor approaches.

Requirements:
    - ANSYS Fluent license
    - PyFluent: pip install ansys-fluent-core
    - Mesh file with named zones for impeller, volute, inlet, outlet

Author: ANSYS Automation
License: Commercial use requires valid ANSYS license
"""

import os
from pathlib import Path
from ansys.fluent.core import launch_fluent


class PumpCFDSetup:
    """Automated setup for centrifugal pump CFD simulation."""

    def __init__(self, mesh_file, precision='double', num_procs=8):
        """
        Initialize pump CFD setup.

        Args:
            mesh_file: Path to mesh file (.cas or .msh)
            precision: 'single' or 'double'
            num_procs: Number of processors for parallel execution
        """
        self.mesh_file = mesh_file
        self.precision = precision
        self.num_procs = num_procs
        self.solver = None

        # Pump operating conditions (default values)
        self.rpm = 1500  # Rotational speed (RPM)
        self.flow_rate = 0.05  # Volume flow rate (m^3/s)
        self.fluid = 'water-liquid'  # Working fluid

        # Zone names (should match mesh)
        self.impeller_zone = 'impeller-fluid'
        self.volute_zone = 'volute-fluid'
        self.inlet_zone = 'inlet'
        self.outlet_zone = 'outlet'
        self.impeller_walls = 'impeller-blades'
        self.volute_walls = 'volute-walls'

    def launch_fluent(self):
        """Launch Fluent session."""
        print(f"Launching Fluent with {self.num_procs} processors...")
        self.solver = launch_fluent(
            precision=self.precision,
            processor_count=self.num_procs,
            dimension=3,
            show_gui=False
        )
        print("Fluent session started")

    def read_mesh(self):
        """Read mesh file."""
        print(f"Reading mesh file: {self.mesh_file}")
        self.solver.tui.file.read_case(self.mesh_file)

        # Check mesh
        print("Checking mesh quality...")
        self.solver.tui.mesh.check()

        # Scale mesh if needed (assuming mesh is in mm, convert to m)
        print("Scaling mesh from mm to m...")
        self.solver.tui.mesh.scale(0.001, 0.001, 0.001)

        # Reorder domain for better performance
        print("Reordering domain...")
        self.solver.tui.mesh.reorder.reorder_domain()

    def setup_models(self):
        """Configure physical models for pump simulation."""
        print("Setting up physical models...")

        # Enable energy equation
        print("  - Enabling energy equation")
        self.solver.tui.define.models.energy('yes')

        # Set turbulence model (k-omega SST recommended for pumps)
        print("  - Setting turbulence model: k-omega SST")
        self.solver.tui.define.models.viscous.kw_sst('yes')

        print("Models configured successfully")

    def setup_materials(self):
        """Configure material properties."""
        print(f"Setting up material: {self.fluid}")

        # Copy water properties from database
        self.solver.tui.define.materials.copy('fluid', self.fluid)

        print("Materials configured")

    def setup_cell_zones(self, motion_type='rotating_reference_frame'):
        """
        Configure cell zones with rotation.

        Args:
            motion_type: 'rotating_reference_frame', 'sliding_mesh', or 'frozen_rotor'
        """
        print(f"Setting up cell zones with {motion_type}...")

        # Convert RPM to rad/s
        omega = self.rpm * 2 * 3.14159265359 / 60

        if motion_type == 'rotating_reference_frame':
            # Impeller zone - rotating reference frame
            print(f"  - Impeller: rotating reference frame ({self.rpm} RPM)")
            self.solver.tui.define.boundary_conditions.fluid(
                self.impeller_zone,
                'yes',
                self.fluid,
                'yes',  # Frame motion
                'no',   # Relative to cell zone
                'yes',  # Rotation
                omega,  # Angular velocity (rad/s)
                0, 0, 1,  # Rotation axis (Z-axis)
                0, 0, 0,  # Rotation origin
                'no',   # Porous zone
                'no'    # Source terms
            )

            # Volute zone - stationary
            print("  - Volute: stationary reference frame")
            self.solver.tui.define.boundary_conditions.fluid(
                self.volute_zone,
                'yes',
                self.fluid,
                'no',  # No frame motion
                'no',
                'no',
                'no',
                'no'
            )

        elif motion_type == 'sliding_mesh':
            print(f"  - Impeller: sliding mesh ({self.rpm} RPM)")
            # Enable sliding mesh motion
            self.solver.tui.define.dynamic_mesh.dynamic_mesh('yes')

            # Set up rotation for impeller zone
            self.solver.tui.define.dynamic_mesh.zones.create(
                self.impeller_zone,
                'rigid-body',
                'yes',  # Rotation
                omega,
                0, 0, 1,  # Rotation axis
                0, 0, 0   # Rotation origin
            )

        print("Cell zones configured")

    def setup_boundary_conditions(self):
        """Configure boundary conditions."""
        print("Setting up boundary conditions...")

        # Calculate inlet velocity from flow rate
        # Assuming circular inlet with diameter = 0.1 m
        inlet_area = 3.14159265359 * (0.05 ** 2)  # m^2
        inlet_velocity = self.flow_rate / inlet_area

        # Inlet - mass flow inlet (more stable for pumps)
        print(f"  - Inlet: mass flow inlet ({self.flow_rate} m^3/s)")
        density = 998.2  # Water density kg/m^3
        mass_flow = self.flow_rate * density

        self.solver.tui.define.boundary_conditions.mass_flow_inlet(
            self.inlet_zone,
            'yes',
            'no',   # Not multiphase
            'yes',  # Mass flow rate
            'yes',  # Normal to boundary
            'no',   # Not components
            mass_flow,
            'no',   # Turbulent intensity
            0.05,
            'no',   # Hydraulic diameter
            0.1,
            'no',   # Temperature
            300,
            'no'
        )

        # Outlet - pressure outlet
        print("  - Outlet: pressure outlet (0 Pa gauge)")
        self.solver.tui.define.boundary_conditions.pressure_outlet(
            self.outlet_zone,
            'yes',
            'no',
            0,      # Gauge pressure = 0
            'no',
            300,
            'no',
            'yes'
        )

        # Impeller walls - rotating wall
        print("  - Impeller walls: rotating wall")
        omega = self.rpm * 2 * 3.14159265359 / 60
        self.solver.tui.define.boundary_conditions.wall(
            self.impeller_walls,
            'yes',
            'no',   # Shear condition
            'yes',  # Moving wall
            'yes',  # Rotation
            omega,  # Angular velocity
            0, 0, 1,  # Rotation axis
            0, 0, 0,  # Origin
            'no',   # Heat flux
            'yes',  # Temperature
            'no',   # Fixed temp
            300
        )

        # Volute walls - stationary wall
        print("  - Volute walls: stationary wall")
        self.solver.tui.define.boundary_conditions.wall(
            self.volute_walls,
            'yes',
            'no',
            'no',   # Stationary wall
            'no',
            0,
            'no',
            'yes',
            'no',
            300
        )

        print("Boundary conditions configured")

    def setup_solution_methods(self):
        """Configure solution methods and controls."""
        print("Setting up solution methods...")

        # Pressure-velocity coupling - SIMPLE for steady-state
        print("  - P-V coupling: SIMPLE")
        self.solver.tui.solve.set.p_v_coupling(20)

        # Discretization schemes
        print("  - Discretization: Second-order")
        self.solver.tui.solve.set.discretization_scheme.pressure(12)  # Second-order
        self.solver.tui.solve.set.discretization_scheme.mom(1)        # Second-order upwind
        self.solver.tui.solve.set.discretization_scheme.k(1)
        self.solver.tui.solve.set.discretization_scheme.omega(1)
        self.solver.tui.solve.set.discretization_scheme.temperature(1)

        # Under-relaxation factors (conservative for pumps)
        print("  - Under-relaxation factors")
        self.solver.tui.solve.set.under_relaxation.pressure(0.3)
        self.solver.tui.solve.set.under_relaxation.mom(0.5)
        self.solver.tui.solve.set.under_relaxation.k(0.6)
        self.solver.tui.solve.set.under_relaxation.omega(0.6)
        self.solver.tui.solve.set.under_relaxation.turb_viscosity(0.8)

        print("Solution methods configured")

    def setup_monitors(self):
        """Set up convergence and performance monitors."""
        print("Setting up monitors...")

        # Residual convergence criteria
        self.solver.tui.solve.monitors.residual.convergence_criteria(
            '1e-5',  # Continuity
            '1e-5',  # x-velocity
            '1e-5',  # y-velocity
            '1e-5',  # z-velocity
            '1e-6',  # Energy
            '1e-5',  # k
            '1e-5'   # omega
        )

        # Mass flow rate monitor at outlet
        print("  - Mass flow monitor at outlet")
        self.solver.tui.solve.monitors.surface.set_monitor(
            'mass-flow-rate',
            self.outlet_zone,
            'Mass Flow Rate',
            'yes',  # Print
            'yes',  # Plot
            'yes',  # Write
            'outlet_massflow.out'
        )

        # Pressure monitor at inlet and outlet for head calculation
        print("  - Pressure monitors for head calculation")
        self.solver.tui.solve.monitors.surface.set_monitor(
            'area-weighted-avg',
            self.inlet_zone,
            'Pressure',
            'yes',
            'yes',
            'yes',
            'inlet_pressure.out'
        )

        self.solver.tui.solve.monitors.surface.set_monitor(
            'area-weighted-avg',
            self.outlet_zone,
            'Pressure',
            'yes',
            'yes',
            'yes',
            'outlet_pressure.out'
        )

        print("Monitors configured")

    def initialize_solution(self):
        """Initialize flow field."""
        print("Initializing solution...")

        # Hybrid initialization (recommended for complex geometries)
        print("  - Using hybrid initialization")
        self.solver.tui.solve.initialize.hyb_initialization()

        print("Solution initialized")

    def run_calculation(self, iterations=2000, save_frequency=500):
        """
        Run steady-state calculation.

        Args:
            iterations: Number of iterations
            save_frequency: Save case/data every N iterations
        """
        print(f"Running calculation for {iterations} iterations...")

        # Auto-save settings
        self.solver.tui.file.auto_save.data_frequency(save_frequency)

        # Run iterations
        self.solver.tui.solve.iterate(iterations)

        print("Calculation completed")

    def postprocess_results(self, output_dir='./pump_results'):
        """
        Extract pump performance data.

        Args:
            output_dir: Directory to save results
        """
        print("Post-processing results...")

        Path(output_dir).mkdir(exist_ok=True)

        # Calculate head rise
        print("  - Calculating pump head...")

        # Report mass flow rate
        print("  - Reporting mass flow rate...")
        self.solver.tui.report.surface_integrals.mass_flow_rate(
            self.outlet_zone, '()'
        )

        # Calculate torque on impeller
        print("  - Calculating torque on impeller blades...")
        self.solver.tui.report.forces.wall_forces(
            'yes',
            self.impeller_walls,
            '()',
            'yes',  # Direction vector
            'no',   # Not absolute
            'yes',  # Compute moments
            0, 0, 0,  # Moment center
            0, 0, 1   # Moment axis (Z)
        )

        # Export pressure distribution on impeller
        print("  - Exporting pressure distribution...")
        pressure_file = f"{output_dir}/impeller_pressure.csv"
        self.solver.tui.file.export.ascii(
            pressure_file,
            self.impeller_walls,
            '()',
            'yes',
            'pressure',
            'velocity-magnitude',
            '()'
        )

        print(f"Results exported to {output_dir}")

    def save_case_data(self, output_file):
        """Save case and data file."""
        print(f"Saving case and data to {output_file}")
        self.solver.tui.file.write_case_data(output_file)

    def cleanup(self):
        """Exit Fluent cleanly."""
        if self.solver and self.solver.is_alive():
            print("Exiting Fluent...")
            self.solver.exit()

    def run_complete_setup(self, output_case='pump_setup.cas.h5'):
        """
        Execute complete pump CFD setup workflow.

        Args:
            output_case: Output case/data file
        """
        try:
            self.launch_fluent()
            self.read_mesh()
            self.setup_models()
            self.setup_materials()
            self.setup_cell_zones(motion_type='rotating_reference_frame')
            self.setup_boundary_conditions()
            self.setup_solution_methods()
            self.setup_monitors()
            self.initialize_solution()
            self.save_case_data(output_case)

            print("\n" + "="*60)
            print("Pump CFD setup completed successfully!")
            print(f"Case saved to: {output_case}")
            print("Ready to run calculation with run_calculation()")
            print("="*60)

        except Exception as e:
            print(f"Error during setup: {e}")
            raise
        finally:
            pass  # Don't cleanup yet if running calculation next


def main():
    """Example usage of PumpCFDSetup."""

    # Configuration
    mesh_file = 'pump_mesh.cas'
    output_case = 'pump_setup.cas.h5'

    # Create setup object
    pump = PumpCFDSetup(
        mesh_file=mesh_file,
        precision='double',
        num_procs=8
    )

    # Set pump parameters
    pump.rpm = 1500  # RPM
    pump.flow_rate = 0.05  # m^3/s

    # Run complete setup
    pump.run_complete_setup(output_case)

    # Run calculation
    pump.run_calculation(iterations=2000, save_frequency=500)

    # Post-process
    pump.postprocess_results(output_dir='./pump_results')

    # Save final results
    pump.save_case_data('pump_final.cas.h5')

    # Cleanup
    pump.cleanup()


if __name__ == '__main__':
    main()
