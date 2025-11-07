#!/usr/bin/env python3
"""
ANSYS Fluent - Batch Simulation Framework

This script provides a framework for running multiple ANSYS Fluent simulations
in batch mode, including parametric studies, design of experiments (DOE),
and automated post-processing.

Features:
    - Parametric study automation
    - Design of experiments (DOE)
    - Batch processing with journal files
    - Resource management and job scheduling
    - Automated post-processing and data extraction
    - Result aggregation and reporting

Requirements:
    - ANSYS Fluent license with appropriate HPC licenses
    - PyFluent: pip install ansys-fluent-core
    - NumPy: pip install numpy
    - Pandas: pip install pandas (for result analysis)

Author: ANSYS Automation
License: Commercial use requires valid ANSYS license
"""

import os
import sys
import json
import subprocess
import time
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import numpy as np
except ImportError:
    print("Warning: NumPy not available. Some features may be limited.")
    np = None

try:
    import pandas as pd
except ImportError:
    print("Warning: Pandas not available. Result analysis features limited.")
    pd = None

from ansys.fluent.core import launch_fluent


class FluentBatchSimulation:
    """Framework for running ANSYS Fluent simulations in batch mode."""

    def __init__(self, base_case, output_dir='./batch_results', num_procs=4):
        """
        Initialize batch simulation framework.

        Args:
            base_case: Base case file path
            output_dir: Directory for simulation results
            num_procs: Number of processors per simulation
        """
        self.base_case = Path(base_case)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.num_procs = num_procs
        self.simulations = []
        self.results = {}

        # Create subdirectories
        (self.output_dir / 'cases').mkdir(exist_ok=True)
        (self.output_dir / 'logs').mkdir(exist_ok=True)
        (self.output_dir / 'data').mkdir(exist_ok=True)

    def add_simulation(self, name: str, parameters: Dict[str, Any]):
        """
        Add a simulation to the batch queue.

        Args:
            name: Simulation identifier
            parameters: Dictionary of parameters to modify
        """
        sim_config = {
            'name': name,
            'parameters': parameters,
            'status': 'queued',
            'created': datetime.now().isoformat()
        }
        self.simulations.append(sim_config)
        print(f"Added simulation: {name}")

    def generate_parametric_study(self, param_ranges: Dict[str, tuple],
                                   num_samples: int = 10, method='linspace'):
        """
        Generate parametric study simulations.

        Args:
            param_ranges: Dict mapping parameter names to (min, max) tuples
            num_samples: Number of samples per parameter
            method: 'linspace' or 'logspace' or 'random'
        """
        if np is None:
            raise ImportError("NumPy required for parametric studies")

        print(f"Generating parametric study with {len(param_ranges)} parameters...")

        # Generate parameter combinations
        param_values = {}
        for param_name, (min_val, max_val) in param_ranges.items():
            if method == 'linspace':
                values = np.linspace(min_val, max_val, num_samples)
            elif method == 'logspace':
                values = np.logspace(np.log10(min_val), np.log10(max_val), num_samples)
            elif method == 'random':
                values = np.random.uniform(min_val, max_val, num_samples)
            else:
                raise ValueError(f"Unknown sampling method: {method}")

            param_values[param_name] = values

        # Create combinations (full factorial or one-at-a-time)
        base_params = {name: values[0] for name, values in param_values.items()}

        sim_count = 0
        for param_name, values in param_values.items():
            for value in values:
                params = base_params.copy()
                params[param_name] = value

                sim_name = f"param_{param_name}_{value:.4f}"
                self.add_simulation(sim_name, params)
                sim_count += 1

        print(f"Generated {sim_count} parametric simulations")

    def generate_doe_study(self, factors: Dict[str, tuple], doe_type='full_factorial'):
        """
        Generate Design of Experiments (DOE) study.

        Args:
            factors: Dict mapping factor names to (min, max) tuples
            doe_type: 'full_factorial', 'fractional_factorial', or 'latin_hypercube'
        """
        if np is None:
            raise ImportError("NumPy required for DOE studies")

        print(f"Generating {doe_type} DOE with {len(factors)} factors...")

        if doe_type == 'full_factorial':
            # Full factorial with 2 levels
            levels = [-1, 1]  # Low and high levels
            combinations = np.array(np.meshgrid(*[levels for _ in factors])).T.reshape(-1, len(factors))

            # Map to actual values
            factor_list = list(factors.keys())
            for i, combo in enumerate(combinations):
                params = {}
                for j, factor_name in enumerate(factor_list):
                    min_val, max_val = factors[factor_name]
                    # Map -1 to min, 1 to max
                    params[factor_name] = min_val if combo[j] == -1 else max_val

                sim_name = f"doe_run_{i+1:03d}"
                self.add_simulation(sim_name, params)

        elif doe_type == 'latin_hypercube':
            # Latin Hypercube Sampling
            n_samples = 20  # Adjust as needed
            n_factors = len(factors)

            # Simple LHS implementation
            samples = np.zeros((n_samples, n_factors))
            for i in range(n_factors):
                samples[:, i] = (np.random.permutation(n_samples) + 0.5) / n_samples

            # Scale to factor ranges
            factor_list = list(factors.keys())
            for i, sample in enumerate(samples):
                params = {}
                for j, factor_name in enumerate(factor_list):
                    min_val, max_val = factors[factor_name]
                    params[factor_name] = min_val + sample[j] * (max_val - min_val)

                sim_name = f"lhs_run_{i+1:03d}"
                self.add_simulation(sim_name, params)

        print(f"Generated {len(self.simulations)} DOE simulations")

    def run_single_simulation(self, sim_config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run a single Fluent simulation.

        Args:
            sim_config: Simulation configuration dictionary

        Returns:
            Result dictionary with simulation outcomes
        """
        name = sim_config['name']
        parameters = sim_config['parameters']

        print(f"\n{'='*60}")
        print(f"Running simulation: {name}")
        print(f"Parameters: {parameters}")
        print(f"{'='*60}\n")

        start_time = time.time()

        try:
            # Launch Fluent
            solver = launch_fluent(
                precision='double',
                processor_count=self.num_procs,
                show_gui=False
            )

            # Read base case
            solver.tui.file.read_case(str(self.base_case))

            # Apply parameters (this is problem-specific)
            self._apply_parameters(solver, parameters)

            # Initialize and solve
            solver.tui.solve.initialize.initialize_flow()
            solver.tui.solve.iterate(1000)  # Adjust iterations as needed

            # Extract results
            results = self._extract_results(solver, name)

            # Save case and data
            case_file = self.output_dir / 'cases' / f"{name}.cas.h5"
            solver.tui.file.write_case_data(str(case_file))

            # Cleanup
            solver.exit()

            elapsed_time = time.time() - start_time

            return {
                'name': name,
                'status': 'completed',
                'parameters': parameters,
                'results': results,
                'elapsed_time': elapsed_time,
                'case_file': str(case_file)
            }

        except Exception as e:
            elapsed_time = time.time() - start_time
            print(f"Error in simulation {name}: {e}")

            return {
                'name': name,
                'status': 'failed',
                'parameters': parameters,
                'error': str(e),
                'elapsed_time': elapsed_time
            }

    def _apply_parameters(self, solver, parameters: Dict[str, Any]):
        """
        Apply parameters to Fluent case (problem-specific).

        Args:
            solver: Fluent solver session
            parameters: Parameters to apply

        Note: This is a template - customize for your specific case
        """
        # Example: Modify inlet velocity
        if 'inlet_velocity' in parameters:
            velocity = parameters['inlet_velocity']
            solver.tui.define.boundary_conditions.velocity_inlet(
                'inlet', 'yes', 'no', 'yes', 'yes', 'no',
                velocity, 'no', 0.05, 'no', 0.1, 'no', 300, 'no'
            )

        # Example: Modify material properties
        if 'density' in parameters:
            density = parameters['density']
            # Modify material density
            # (Implementation depends on your setup)

        # Example: Modify temperature
        if 'temperature' in parameters:
            temp = parameters['temperature']
            solver.tui.define.boundary_conditions.velocity_inlet(
                'inlet', 'yes', 'no', 'no', 'no',
                temp, 'no'
            )

        # Add more parameter modifications as needed

    def _extract_results(self, solver, sim_name: str) -> Dict[str, Any]:
        """
        Extract results from simulation (problem-specific).

        Args:
            solver: Fluent solver session
            sim_name: Simulation name

        Returns:
            Dictionary of extracted results

        Note: This is a template - customize for your specific case
        """
        results = {}

        # Example: Extract mass flow rate
        # Note: Actual implementation depends on your case setup
        try:
            # Export surface data
            data_file = self.output_dir / 'data' / f"{sim_name}_surface.csv"
            solver.tui.file.export.ascii(
                str(data_file),
                'outlet',
                '()',
                'yes',
                'pressure',
                'velocity-magnitude',
                '()'
            )
            results['data_file'] = str(data_file)
        except:
            pass

        # Add more result extraction as needed
        results['converged'] = True  # Check actual convergence
        results['timestamp'] = datetime.now().isoformat()

        return results

    def run_sequential(self):
        """Run all simulations sequentially."""
        print(f"\nRunning {len(self.simulations)} simulations sequentially...")

        for i, sim_config in enumerate(self.simulations, 1):
            print(f"\nSimulation {i}/{len(self.simulations)}")
            result = self.run_single_simulation(sim_config)
            self.results[result['name']] = result

            # Save intermediate results
            self._save_results()

        print("\nAll sequential simulations completed")

    def run_parallel(self, max_workers: Optional[int] = None):
        """
        Run simulations in parallel.

        Args:
            max_workers: Maximum number of parallel workers (default: CPU count / num_procs)
        """
        if max_workers is None:
            max_workers = max(1, multiprocessing.cpu_count() // self.num_procs)

        print(f"\nRunning {len(self.simulations)} simulations in parallel")
        print(f"Max workers: {max_workers}")

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs
            future_to_sim = {
                executor.submit(self.run_single_simulation, sim): sim
                for sim in self.simulations
            }

            # Process completed simulations
            for i, future in enumerate(as_completed(future_to_sim), 1):
                sim = future_to_sim[future]
                try:
                    result = future.result()
                    self.results[result['name']] = result
                    print(f"\nCompleted {i}/{len(self.simulations)}: {result['name']}")

                    # Save intermediate results
                    self._save_results()

                except Exception as e:
                    print(f"\nSimulation {sim['name']} generated an exception: {e}")

        print("\nAll parallel simulations completed")

    def run_with_journal(self, journal_template: str):
        """
        Run simulations using journal file template.

        Args:
            journal_template: Path to journal file template
        """
        print(f"Running simulations with journal file: {journal_template}")

        for sim_config in self.simulations:
            name = sim_config['name']
            parameters = sim_config['parameters']

            # Generate journal file from template
            journal_file = self._generate_journal(journal_template, name, parameters)

            # Run Fluent in batch mode
            log_file = self.output_dir / 'logs' / f"{name}.log"
            success = self._run_fluent_batch(journal_file, log_file)

            # Store result
            self.results[name] = {
                'name': name,
                'status': 'completed' if success else 'failed',
                'parameters': parameters,
                'log_file': str(log_file)
            }

            self._save_results()

        print("All journal-based simulations completed")

    def _generate_journal(self, template: str, name: str,
                         parameters: Dict[str, Any]) -> Path:
        """
        Generate journal file from template with parameter substitution.

        Args:
            template: Journal template file path
            name: Simulation name
            parameters: Parameters to substitute

        Returns:
            Path to generated journal file
        """
        with open(template, 'r') as f:
            content = f.read()

        # Simple parameter substitution
        for param_name, param_value in parameters.items():
            placeholder = f"{{{param_name}}}"
            content = content.replace(placeholder, str(param_value))

        # Save generated journal
        journal_file = self.output_dir / 'cases' / f"{name}.jou"
        with open(journal_file, 'w') as f:
            f.write(content)

        return journal_file

    def _run_fluent_batch(self, journal_file: Path, log_file: Path) -> bool:
        """
        Run Fluent in batch mode with journal file.

        Args:
            journal_file: Journal file path
            log_file: Log file path

        Returns:
            True if successful, False otherwise
        """
        cmd = [
            'fluent',
            '3ddp',  # 3D double precision
            '-g',    # No GUI
            f'-t{self.num_procs}',
            '-i', str(journal_file)
        ]

        try:
            with open(log_file, 'w') as log:
                result = subprocess.run(
                    cmd,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=3600  # 1 hour timeout
                )
            return result.returncode == 0
        except Exception as e:
            print(f"Error running Fluent batch: {e}")
            return False

    def _save_results(self):
        """Save results to JSON file."""
        results_file = self.output_dir / 'results_summary.json'
        with open(results_file, 'w') as f:
            json.dump(self.results, f, indent=2)

    def generate_report(self) -> str:
        """
        Generate summary report of batch simulations.

        Returns:
            Report as formatted string
        """
        total = len(self.results)
        completed = sum(1 for r in self.results.values() if r['status'] == 'completed')
        failed = sum(1 for r in self.results.values() if r['status'] == 'failed')

        total_time = sum(r.get('elapsed_time', 0) for r in self.results.values())

        report = f"""
{'='*70}
ANSYS FLUENT BATCH SIMULATION REPORT
{'='*70}

Base Case: {self.base_case}
Output Directory: {self.output_dir}
Processors per Simulation: {self.num_procs}

SUMMARY
-------
Total Simulations: {total}
Completed: {completed}
Failed: {failed}
Total Compute Time: {total_time/3600:.2f} hours

RESULTS
-------
"""
        for name, result in self.results.items():
            status = result['status']
            params = result['parameters']
            report += f"\n{name}:\n"
            report += f"  Status: {status}\n"
            report += f"  Parameters: {params}\n"
            if 'elapsed_time' in result:
                report += f"  Time: {result['elapsed_time']/60:.2f} minutes\n"

        report += f"\n{'='*70}\n"
        report += f"Report generated: {datetime.now().isoformat()}\n"
        report += f"{'='*70}\n"

        # Save report
        report_file = self.output_dir / 'batch_report.txt'
        with open(report_file, 'w') as f:
            f.write(report)

        return report


def example_parametric_study():
    """Example: Parametric study varying inlet velocity."""
    print("Example: Parametric study - varying inlet velocity")

    # Initialize batch framework
    batch = FluentBatchSimulation(
        base_case='pipe_flow.cas',
        output_dir='./velocity_study',
        num_procs=4
    )

    # Define parameter ranges
    param_ranges = {
        'inlet_velocity': (1.0, 10.0),  # m/s
    }

    # Generate parametric study
    batch.generate_parametric_study(param_ranges, num_samples=10, method='linspace')

    # Run sequentially
    batch.run_sequential()

    # Generate report
    report = batch.generate_report()
    print(report)


def example_doe_study():
    """Example: Design of Experiments study."""
    print("Example: DOE study - full factorial design")

    batch = FluentBatchSimulation(
        base_case='heat_exchanger.cas',
        output_dir='./doe_study',
        num_procs=4
    )

    # Define factors
    factors = {
        'inlet_velocity': (1.0, 5.0),
        'inlet_temperature': (300, 400),
        'wall_temperature': (350, 450)
    }

    # Generate DOE
    batch.generate_doe_study(factors, doe_type='full_factorial')

    # Run in parallel
    batch.run_parallel(max_workers=2)

    # Generate report
    report = batch.generate_report()
    print(report)


def main():
    """Main function demonstrating batch simulation capabilities."""
    print("ANSYS Fluent - Batch Simulation Framework")
    print("="*70)

    # Uncomment to run examples:
    # example_parametric_study()
    # example_doe_study()

    print("\nBatch simulation framework ready.")
    print("Use FluentBatchSimulation class to set up your study.")


if __name__ == '__main__':
    main()
