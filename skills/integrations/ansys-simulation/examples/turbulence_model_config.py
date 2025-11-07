#!/usr/bin/env python3
"""
ANSYS Fluent - Turbulence Model Configuration and Comparison

This script provides utilities for configuring various turbulence models
in ANSYS Fluent and running comparison studies to evaluate model performance.

Turbulence Models Supported:
    - RANS: k-epsilon (Standard, RNG, Realizable)
    - RANS: k-omega (Standard, SST)
    - RANS: Spalart-Allmaras
    - Scale-Resolving: LES, DES, DDES, SAS

Requirements:
    - ANSYS Fluent license
    - PyFluent: pip install ansys-fluent-core

Author: ANSYS Automation
License: Commercial use requires valid ANSYS license
"""

import os
from pathlib import Path
from ansys.fluent.core import launch_fluent
from typing import Dict, List
import json


class TurbulenceModelConfig:
    """Configure and manage turbulence models in ANSYS Fluent."""

    # Turbulence model definitions
    RANS_MODELS = {
        'laminar': 'Laminar (no turbulence)',
        'spalart-allmaras': 'Spalart-Allmaras (1-equation)',
        'k-epsilon-standard': 'k-epsilon Standard',
        'k-epsilon-rng': 'k-epsilon RNG',
        'k-epsilon-realizable': 'k-epsilon Realizable',
        'k-omega-standard': 'k-omega Standard',
        'k-omega-sst': 'k-omega SST',
        'reynolds-stress': 'Reynolds Stress Model (RSM)',
    }

    SCALE_RESOLVING = {
        'les': 'Large Eddy Simulation',
        'des': 'Detached Eddy Simulation',
        'ddes': 'Delayed DES',
        'sas': 'Scale-Adaptive Simulation',
    }

    def __init__(self, solver=None):
        """
        Initialize turbulence configuration.

        Args:
            solver: Existing Fluent solver session (or None to create new)
        """
        self.solver = solver
        self.current_model = None

    def set_laminar(self):
        """Set laminar flow (no turbulence model)."""
        print("Setting flow to laminar...")
        self.solver.tui.define.models.viscous.laminar()
        self.current_model = 'laminar'
        print("Laminar flow configured")

    def set_spalart_allmaras(self, production='vorticity', version='standard'):
        """
        Configure Spalart-Allmaras turbulence model.

        Args:
            production: 'vorticity' or 'strain-vorticity'
            version: 'standard' or 'des'
        """
        print("Setting Spalart-Allmaras turbulence model...")

        if version == 'des':
            # DES variant
            self.solver.tui.define.models.viscous.spalart_allmaras(
                'yes',  # Enable
                'yes'   # DES
            )
            self.current_model = 'spalart-allmaras-des'
        else:
            # Standard
            self.solver.tui.define.models.viscous.spalart_allmaras('yes')
            self.current_model = 'spalart-allmaras'

        print(f"Spalart-Allmaras ({version}) configured")

    def set_k_epsilon(self, variant='realizable', near_wall='standard'):
        """
        Configure k-epsilon turbulence model.

        Args:
            variant: 'standard', 'rng', or 'realizable'
            near_wall: 'standard', 'non-equilibrium', 'enhanced', or 'scalable'
        """
        print(f"Setting k-epsilon {variant} turbulence model...")

        if variant == 'standard':
            self.solver.tui.define.models.viscous.ke_standard('yes')
            self.current_model = 'k-epsilon-standard'
        elif variant == 'rng':
            self.solver.tui.define.models.viscous.ke_rng('yes')
            self.current_model = 'k-epsilon-rng'
        elif variant == 'realizable':
            self.solver.tui.define.models.viscous.ke_realizable('yes')
            self.current_model = 'k-epsilon-realizable'
        else:
            raise ValueError(f"Unknown k-epsilon variant: {variant}")

        # Configure near-wall treatment
        if near_wall == 'enhanced':
            print("  - Enabling enhanced wall treatment")
            self.solver.tui.define.models.viscous.near_wall_treatment.enhanced_wall_treatment('yes')

        print(f"k-epsilon {variant} configured")

    def set_k_omega(self, variant='sst', options=None):
        """
        Configure k-omega turbulence model.

        Args:
            variant: 'standard' or 'sst'
            options: Dictionary of additional options
                - low_re_corrections: bool (Low Reynolds number corrections)
                - compressibility: bool (Compressibility effects)
                - shear_flow: bool (Shear flow corrections for SST)
        """
        print(f"Setting k-omega {variant} turbulence model...")

        if options is None:
            options = {}

        if variant == 'standard':
            self.solver.tui.define.models.viscous.kw_standard('yes')
            self.current_model = 'k-omega-standard'
        elif variant == 'sst':
            self.solver.tui.define.models.viscous.kw_sst('yes')
            self.current_model = 'k-omega-sst'

            # SST-specific options
            if options.get('low_re_corrections', False):
                print("  - Enabling low-Re corrections")
                # Access via TUI
                self.solver.tui.define.models.viscous.turbulence_expert.low_re_corrections('yes')

            if options.get('compressibility', False):
                print("  - Enabling compressibility effects")
                self.solver.tui.define.models.viscous.turbulence_expert.compressibility_effects('yes')

        else:
            raise ValueError(f"Unknown k-omega variant: {variant}")

        print(f"k-omega {variant} configured")

    def set_les(self, subgrid_model='smagorinsky', filter_width='cube-root'):
        """
        Configure Large Eddy Simulation (LES).

        Args:
            subgrid_model: 'smagorinsky', 'wale', 'wmles', or 'kinetic-energy-transport'
            filter_width: 'cube-root', 'equiangle-skew', or 'user-defined'
        """
        print(f"Setting LES with {subgrid_model} subgrid model...")

        # Enable LES
        self.solver.tui.define.models.viscous.les('yes')

        # Set subgrid-scale model
        if subgrid_model == 'smagorinsky':
            self.solver.tui.define.models.viscous.les_model.smagorinsky_lilly()
        elif subgrid_model == 'wale':
            self.solver.tui.define.models.viscous.les_model.wale()
        elif subgrid_model == 'wmles':
            self.solver.tui.define.models.viscous.les_model.wmles()
        elif subgrid_model == 'kinetic-energy-transport':
            self.solver.tui.define.models.viscous.les_model.kinetic_energy_transport()
        else:
            raise ValueError(f"Unknown LES subgrid model: {subgrid_model}")

        self.current_model = f'les-{subgrid_model}'
        print(f"LES configured with {subgrid_model} model")

    def set_des(self, variant='ddes', rans_model='sst'):
        """
        Configure Detached Eddy Simulation (DES).

        Args:
            variant: 'des', 'ddes' (Delayed DES), or 'iddes' (Improved DDES)
            rans_model: 'spalart-allmaras' or 'sst'
        """
        print(f"Setting {variant.upper()} with {rans_model} RANS model...")

        if rans_model == 'spalart-allmaras':
            # Spalart-Allmaras based DES
            self.solver.tui.define.models.viscous.spalart_allmaras(
                'yes',  # Enable
                'yes'   # DES
            )

            if variant == 'ddes':
                self.solver.tui.define.models.viscous.des_options.ddes('yes')
                self.current_model = 'ddes-sa'
            elif variant == 'iddes':
                self.solver.tui.define.models.viscous.des_options.iddes('yes')
                self.current_model = 'iddes-sa'
            else:
                self.current_model = 'des-sa'

        elif rans_model == 'sst':
            # k-omega SST based DES
            self.solver.tui.define.models.viscous.kw_sst('yes')
            self.solver.tui.define.models.viscous.des_options.sst_des('yes')

            if variant == 'ddes':
                self.solver.tui.define.models.viscous.des_options.ddes('yes')
                self.current_model = 'ddes-sst'
            else:
                self.current_model = 'des-sst'

        else:
            raise ValueError(f"Unknown RANS model for DES: {rans_model}")

        print(f"{variant.upper()} configured")

    def set_reynolds_stress(self, model_type='linear-pressure-strain'):
        """
        Configure Reynolds Stress Model (RSM).

        Args:
            model_type: 'linear-pressure-strain', 'quadratic-pressure-strain', or 'stress-omega'
        """
        print(f"Setting Reynolds Stress Model ({model_type})...")

        if model_type == 'linear-pressure-strain':
            self.solver.tui.define.models.viscous.reynolds_stress('yes', 'linear-pressure-strain')
        elif model_type == 'quadratic-pressure-strain':
            self.solver.tui.define.models.viscous.reynolds_stress('yes', 'quadratic-pressure-strain')
        elif model_type == 'stress-omega':
            self.solver.tui.define.models.viscous.reynolds_stress('yes', 'stress-omega')
        else:
            raise ValueError(f"Unknown RSM type: {model_type}")

        self.current_model = f'rsm-{model_type}'
        print("Reynolds Stress Model configured")

    def get_model_recommendations(self, flow_type):
        """
        Get turbulence model recommendations for flow type.

        Args:
            flow_type: Type of flow ('external-aero', 'internal-flow', 'separation',
                      'rotating', 'heat-transfer', 'multiphase', 'free-shear')

        Returns:
            dict: Recommended models and notes
        """
        recommendations = {
            'external-aero': {
                'primary': 'k-omega-sst',
                'alternatives': ['spalart-allmaras', 'k-epsilon-realizable'],
                'notes': 'SST excellent for adverse pressure gradients and boundary layers'
            },
            'internal-flow': {
                'primary': 'k-epsilon-realizable',
                'alternatives': ['k-omega-sst', 'k-epsilon-rng'],
                'notes': 'Realizable k-epsilon good for pipes, channels, pumps'
            },
            'separation': {
                'primary': 'k-omega-sst',
                'alternatives': ['ddes-sst', 'les'],
                'notes': 'SST or scale-resolving methods for flow separation'
            },
            'rotating': {
                'primary': 'k-omega-sst',
                'alternatives': ['rsm', 'ddes-sst'],
                'notes': 'SST with curvature correction for rotating flows'
            },
            'heat-transfer': {
                'primary': 'k-omega-sst',
                'alternatives': ['k-epsilon-realizable', 'rsm'],
                'notes': 'Good near-wall resolution required, consider enhanced wall treatment'
            },
            'multiphase': {
                'primary': 'k-epsilon-realizable',
                'alternatives': ['k-omega-sst', 'rsm'],
                'notes': 'Realizable k-epsilon commonly used for multiphase'
            },
            'free-shear': {
                'primary': 'k-epsilon-rng',
                'alternatives': ['k-omega-sst', 'spalart-allmaras'],
                'notes': 'RNG k-epsilon good for jets, wakes, mixing layers'
            }
        }

        return recommendations.get(flow_type, {
            'primary': 'k-omega-sst',
            'alternatives': ['k-epsilon-realizable'],
            'notes': 'General purpose: SST is robust choice'
        })

    def configure_turbulence_numerics(self, accuracy='high'):
        """
        Set recommended numerical schemes for turbulence.

        Args:
            accuracy: 'low' (first-order), 'medium' (bounded second-order),
                     'high' (second-order)
        """
        print(f"Configuring turbulence numerics for {accuracy} accuracy...")

        if accuracy == 'low':
            # First-order upwind (very stable, less accurate)
            self.solver.tui.solve.set.discretization_scheme.k(0)
            self.solver.tui.solve.set.discretization_scheme.omega(0)
            self.solver.tui.solve.set.discretization_scheme.epsilon(0)

        elif accuracy == 'medium':
            # Bounded second-order (good balance)
            # Note: In newer versions, use second-order with bounded option
            self.solver.tui.solve.set.discretization_scheme.k(1)
            self.solver.tui.solve.set.discretization_scheme.omega(1)
            self.solver.tui.solve.set.discretization_scheme.epsilon(1)

        elif accuracy == 'high':
            # Pure second-order upwind (most accurate, may be less stable)
            self.solver.tui.solve.set.discretization_scheme.k(1)
            self.solver.tui.solve.set.discretization_scheme.omega(1)
            self.solver.tui.solve.set.discretization_scheme.epsilon(1)

        print("Turbulence numerics configured")


class TurbulenceModelComparison:
    """Run comparison studies across different turbulence models."""

    def __init__(self, base_case, models_to_test, output_dir='./turbulence_comparison'):
        """
        Initialize comparison study.

        Args:
            base_case: Base case file (.cas)
            models_to_test: List of turbulence models to compare
            output_dir: Output directory for results
        """
        self.base_case = base_case
        self.models_to_test = models_to_test
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.results = {}

    def run_comparison(self, iterations=1000, num_procs=4):
        """
        Run simulations with each turbulence model.

        Args:
            iterations: Number of iterations per model
            num_procs: Number of processors
        """
        print("="*60)
        print("TURBULENCE MODEL COMPARISON STUDY")
        print("="*60)

        for model_name in self.models_to_test:
            print(f"\n--- Testing Model: {model_name} ---")

            # Launch Fluent
            solver = launch_fluent(precision='double', processor_count=num_procs)

            # Read base case
            solver.tui.file.read_case(self.base_case)

            # Configure turbulence model
            turb_config = TurbulenceModelConfig(solver)

            if model_name == 'k-omega-sst':
                turb_config.set_k_omega('sst')
            elif model_name == 'k-epsilon-realizable':
                turb_config.set_k_epsilon('realizable')
            elif model_name == 'spalart-allmaras':
                turb_config.set_spalart_allmaras()
            elif model_name == 'k-epsilon-standard':
                turb_config.set_k_epsilon('standard')
            # Add more models as needed

            # Initialize and solve
            solver.tui.solve.initialize.initialize_flow()
            solver.tui.solve.iterate(iterations)

            # Save results
            output_file = self.output_dir / f"result_{model_name}.cas.h5"
            solver.tui.file.write_case_data(str(output_file))

            # Store result path
            self.results[model_name] = str(output_file)

            # Cleanup
            solver.exit()

            print(f"Model {model_name} completed")

        print("\n" + "="*60)
        print("COMPARISON STUDY COMPLETED")
        print("="*60)

        # Save summary
        self.save_summary()

    def save_summary(self):
        """Save comparison summary to JSON."""
        summary_file = self.output_dir / 'comparison_summary.json'

        summary = {
            'base_case': self.base_case,
            'models_tested': self.models_to_test,
            'results': self.results
        }

        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        print(f"\nSummary saved to: {summary_file}")


def example_flow_over_airfoil():
    """Example: Configure turbulence for flow over airfoil."""
    print("Example: Flow over airfoil - Turbulence configuration")

    # Launch Fluent
    solver = launch_fluent(precision='double', processor_count=4)

    # Read mesh
    solver.tui.file.read_case('airfoil_mesh.cas')

    # Create turbulence configurator
    turb = TurbulenceModelConfig(solver)

    # Get recommendations
    recommendations = turb.get_model_recommendations('external-aero')
    print(f"\nRecommendation: {recommendations['primary']}")
    print(f"Notes: {recommendations['notes']}")

    # Configure k-omega SST (best for external aero)
    turb.set_k_omega('sst', options={'compressibility': True})

    # Set appropriate numerics
    turb.configure_turbulence_numerics(accuracy='high')

    # Complete setup and save
    # ... (boundary conditions, solution methods, etc.)

    solver.tui.file.write_case('airfoil_setup.cas.h5')
    solver.exit()

    print("Airfoil turbulence setup complete")


def example_turbulence_comparison():
    """Example: Compare multiple turbulence models."""
    print("Example: Turbulence model comparison study")

    # Models to test
    models = [
        'k-omega-sst',
        'k-epsilon-realizable',
        'spalart-allmaras',
        'k-epsilon-standard'
    ]

    # Run comparison
    comparison = TurbulenceModelComparison(
        base_case='pipe_flow.cas',
        models_to_test=models,
        output_dir='./pipe_turbulence_study'
    )

    comparison.run_comparison(iterations=1500, num_procs=4)

    print("\nComparison complete. Review results in ./pipe_turbulence_study")


def main():
    """Main function demonstrating turbulence model configuration."""

    print("ANSYS Fluent - Turbulence Model Configuration")
    print("="*60)

    # Example 1: Flow over airfoil
    # example_flow_over_airfoil()

    # Example 2: Model comparison study
    # example_turbulence_comparison()

    # Example 3: Print recommendations for different flow types
    turb = TurbulenceModelConfig()

    flow_types = ['external-aero', 'internal-flow', 'separation',
                  'rotating', 'heat-transfer', 'multiphase']

    print("\nTurbulence Model Recommendations by Flow Type:")
    print("="*60)

    for flow_type in flow_types:
        rec = turb.get_model_recommendations(flow_type)
        print(f"\n{flow_type.upper()}:")
        print(f"  Primary: {rec['primary']}")
        print(f"  Alternatives: {', '.join(rec['alternatives'])}")
        print(f"  Notes: {rec['notes']}")


if __name__ == '__main__':
    main()
