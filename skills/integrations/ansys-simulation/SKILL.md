---
name: ansys-simulation
description: "Automate ANSYS Fluent CFD simulations via Python scripting and journal files"
category: integrations
domain: fluids
complexity: advanced
dependencies: []
---

# ANSYS Simulation Integration

Comprehensive guide for automating ANSYS Fluent CFD (Computational Fluid Dynamics) simulations using Python scripting, journal files, and the PyAnsys ecosystem.

## Overview

ANSYS Fluent is a leading commercial CFD solver used for simulating fluid flow, heat transfer, and related phenomena in complex geometries. This skill covers automation approaches for:

- **ANSYS Workbench**: Unified platform for CAD integration, meshing, and workflow management
- **ANSYS Fluent**: Advanced CFD solver with extensive physics models
- **PyAnsys**: Python API ecosystem for programmatic control
- **Journal Files**: Text-based scripting for batch operations

## ANSYS Workbench and Fluent Architecture

### ANSYS Workbench
- Project-based workflow environment
- Integrates multiple ANSYS tools (Fluent, CFX, Mechanical, etc.)
- Supports parametric studies and design optimization
- Provides geometry and meshing tools (DesignModeler, SpaceClaim, Meshing)

### ANSYS Fluent
- Finite volume-based CFD solver
- Extensive turbulence models (k-epsilon, k-omega SST, LES, DES)
- Multiphase flow capabilities
- Heat transfer and combustion modeling
- User-defined functions (UDFs) in C
- Text User Interface (TUI) for scripting

## Licensing Requirements

### Commercial Licensing
ANSYS products require commercial licenses managed through the ANSYS License Manager:

- **Fluent License Types**:
  - HPC licenses for parallel computing
  - Solver licenses (increment-based)
  - PrepPost licenses for pre/post-processing

- **License Server Setup**:
  ```bash
  # Set license server environment variable
  export ANSYSLMD_LICENSE_FILE=1055@license-server.company.com

  # Or in Windows
  set ANSYSLMD_LICENSE_FILE=1055@license-server.company.com
  ```

- **License Checkout**: Licenses are checked out when starting Fluent and returned on exit
- **HPC Licensing**: Additional licenses required for parallel runs (typically cores/4)

### Academic Licensing
- Student/teaching licenses available with limitations
- Typically limited mesh size and solver capabilities
- Separate license file required

### Important Notes
- Always verify license availability before batch runs
- Consider license usage in automated workflows
- Use `-t` flag to specify number of processors for HPC licensing

## PyAnsys - Python API Introduction

PyAnsys provides a Pythonic interface to ANSYS products, enabling:
- Remote or local session management
- Programmatic control of simulation setup
- Data extraction and post-processing
- Integration with Python data science stack

### Installation
```bash
# Install PyFluent (for Fluent automation)
pip install ansys-fluent-core

# Install other PyAnsys packages as needed
pip install ansys-mapdl-core      # For Mechanical
pip install ansys-dpf-core        # For post-processing
pip install ansys-geometry        # For geometry operations
```

### Basic PyFluent Usage
```python
from ansys.fluent.core import launch_fluent

# Launch Fluent session
session = launch_fluent(precision='double', processor_count=4, mode='solver')

# Access TUI commands
session.tui.file.read_case('case_file.cas')

# Execute commands
session.tui.solve.iterate(100)

# Clean exit
session.exit()
```

### Key PyAnsys Components
- **ansys-fluent-core**: Main Fluent automation package
- **ansys-fluent-visualization**: Post-processing and visualization
- **ansys-geometry**: CAD operations and geometry manipulation
- **ansys-meshing**: Automated meshing workflows

## Journal File Scripting

Journal files are text-based scripts containing TUI commands executed sequentially by Fluent.

### Journal File Format
```scheme
; Comments start with semicolon
; Commands follow TUI menu structure

/file/read-case case_file.cas

/solve/initialize/initialize-flow
/solve/iterate 1000

/file/write-data result.dat

/exit yes
```

### Running Journal Files
```bash
# Run Fluent in batch mode with journal file
fluent 3ddp -g -i journal.jou > fluent.log

# Flags:
# 3ddp: 3D double precision
# -g: Run without GUI
# -i: Input journal file
# -t4: Use 4 processors
```

### Journal File Best Practices
- Use absolute paths for file I/O
- Include error checking where possible
- Add comments for maintainability
- Test interactively before batch execution
- Redirect output to log files for debugging

## Common Workflow

### 1. Geometry Import

**Via PyFluent:**
```python
from ansys.fluent.core import launch_fluent

solver = launch_fluent(precision='double', processor_count=4)

# Import geometry
solver.tui.file.import_.cad_geometry(
    'geometry.scdoc',
    'yes'  # Confirm import
)
```

**Via Journal File:**
```scheme
/file/import/cad-geometry geometry.step yes
```

### 2. Meshing Setup

**Mesh Generation Strategy:**
- Surface mesh extraction from CAD
- Volume mesh generation (tetrahedral, hex, poly)
- Boundary layer refinement for wall-bounded flows
- Mesh quality checks

**PyFluent Meshing:**
```python
from ansys.fluent.core import launch_fluent

meshing = launch_fluent(precision='double', processor_count=4, mode='meshing')

# Import geometry
meshing.workflow.InitializeWorkflow(WorkflowType='Watertight Geometry')

# Set meshing parameters
meshing.meshing.GlobalSettings.set_state({
    'MaxSize': 0.1,
    'MinSize': 0.001,
    'GrowthRate': 1.2
})

# Generate mesh
meshing.meshing.GenerateMesh()

# Switch to solver
solver = meshing.switch_to_solver()
```

**Journal File Meshing:**
```scheme
/mesh/scale 0.001 0.001 0.001
/mesh/check
/mesh/repair-improve/repair
```

### 3. Physics Setup

**Turbulence Models:**
- k-epsilon: General purpose, wall functions
- k-omega SST: Better for adverse pressure gradients, boundary layers
- LES/DES: Large eddy simulation for transient, high-fidelity

**Material Properties:**
- Fluid properties (density, viscosity, thermal conductivity)
- Solid properties for conjugate heat transfer

**PyFluent Physics Setup:**
```python
# Enable energy equation
solver.tui.define.models.energy('yes')

# Set turbulence model (k-omega SST)
solver.tui.define.models.viscous.kw_sst('yes')

# Define material properties
solver.tui.define.materials.change_create(
    'water-liquid',
    'water',
    'yes',
    'constant',
    998.2,  # Density
    'yes',
    'constant',
    0.001003  # Viscosity
)

# Assign material to cell zone
solver.tui.define.boundary_conditions.fluid(
    'fluid-zone',
    'yes',
    'water',
    'no',
    'no',
    'no',
    'no',
    0,
    'no',
    0,
    'no',
    0,
    'no',
    'no'
)
```

**Journal File Physics:**
```scheme
/define/models/energy yes no no no no
/define/models/viscous kw-sst yes
/define/materials/change-create air air yes constant 1.225 yes constant 1.7894e-05
```

### 4. Boundary Conditions

**Common BC Types:**
- Velocity inlet
- Pressure outlet
- Wall (stationary, moving, rotating)
- Symmetry
- Periodic

**PyFluent BC Setup:**
```python
# Velocity inlet
solver.tui.define.boundary_conditions.velocity_inlet(
    'inlet',
    'yes',
    'no',
    'yes',
    'yes',
    'no',
    10,  # Velocity magnitude (m/s)
    'no',
    0,
    'no',
    0,
    'no',
    300,  # Temperature
    'no'
)

# Pressure outlet
solver.tui.define.boundary_conditions.pressure_outlet(
    'outlet',
    'yes',
    'no',
    0,  # Gauge pressure
    'no',
    300,
    'no',
    'yes'
)

# Wall with no-slip
solver.tui.define.boundary_conditions.wall(
    'wall',
    'yes',
    'no',
    'no',
    0,
    'no',
    0,
    'no',
    'yes',
    'no',
    0
)
```

### 5. Solver Settings

**Solution Methods:**
- Pressure-velocity coupling (SIMPLE, SIMPLEC, PISO, Coupled)
- Discretization schemes (first-order, second-order upwind, QUICK)
- Gradient evaluation (least squares, Green-Gauss)

**PyFluent Solver Configuration:**
```python
# Set pressure-velocity coupling
solver.tui.solve.set.p_v_coupling(24)  # 24 = Coupled

# Set discretization schemes
solver.tui.solve.set.discretization_scheme.pressure(12)  # Second order
solver.tui.solve.set.discretization_scheme.mom(1)  # Second order upwind
solver.tui.solve.set.discretization_scheme.k(1)
solver.tui.solve.set.discretization_scheme.omega(1)

# Set under-relaxation factors
solver.tui.solve.set.under_relaxation.pressure(0.3)
solver.tui.solve.set.under_relaxation.mom(0.7)

# Initialize flow field
solver.tui.solve.initialize.initialize_flow()

# Run iterations
solver.tui.solve.iterate(1000)
```

**Convergence Criteria:**
```python
# Set residual convergence criteria
solver.tui.solve.monitors.residual.convergence_criteria(
    '1e-6',  # continuity
    '1e-6',  # x-velocity
    '1e-6',  # y-velocity
    '1e-6',  # z-velocity
    '1e-6',  # energy
    '1e-6',  # k
    '1e-6'   # omega
)
```

### 6. Post-Processing

**Data Extraction:**
```python
# Create surface for reporting
solver.tui.surface.iso_surface(
    'pressure',
    'iso-surface-1',
    'fluid-zone',
    '()',
    101325
)

# Export data
solver.tui.file.export.ascii(
    'surface-data.csv',
    'iso-surface-1',
    '()',
    'yes',
    'pressure',
    'velocity-magnitude',
    '()'
)

# Calculate force on wall
solver.tui.report.forces.wall_forces(
    'yes',
    'wall-zone',
    '()',
    'yes',
    'no',
    'yes'
)
```

**Visualization:**
```python
# Generate contour plot (requires GUI mode or visualization package)
solver.tui.display.set.contours.filled_contours('yes')
solver.tui.display.contour(
    'pressure',
    'pressure',
    0,
    0,
    'surface-1',
    '()'
)

# Save image
solver.tui.display.save_picture('pressure_contour.png')
```

## Batch Mode Execution

### Linux/Unix Batch Execution
```bash
#!/bin/bash
# batch_fluent.sh

# Set license server
export ANSYSLMD_LICENSE_FILE=1055@license-server.com

# Set number of processors
NPROCS=8

# Run Fluent in batch mode
fluent 3ddp -g -t${NPROCS} -i simulation.jou > fluent_${NPROCS}cores.log 2>&1

# Check exit status
if [ $? -eq 0 ]; then
    echo "Simulation completed successfully"
else
    echo "Simulation failed - check log file"
    exit 1
fi
```

### Windows Batch Execution
```batch
@echo off
REM batch_fluent.bat

set ANSYSLMD_LICENSE_FILE=1055@license-server.com
set NPROCS=8

"C:\Program Files\ANSYS Inc\v241\fluent\ntbin\win64\fluent.exe" 3ddp -g -t%NPROCS% -i simulation.jou > fluent.log 2>&1

if %ERRORLEVEL% EQU 0 (
    echo Simulation completed successfully
) else (
    echo Simulation failed - check log file
    exit /b 1
)
```

### Python-Based Batch Execution
```python
import subprocess
import os
from pathlib import Path

def run_fluent_batch(journal_file, num_procs=4, log_file='fluent.log'):
    """
    Run Fluent in batch mode with journal file.

    Args:
        journal_file: Path to journal file
        num_procs: Number of processors
        log_file: Output log file name

    Returns:
        bool: True if successful, False otherwise
    """
    # Set license server
    os.environ['ANSYSLMD_LICENSE_FILE'] = '1055@license-server.com'

    # Build command
    cmd = [
        'fluent',
        '3ddp',           # 3D double precision
        '-g',             # No GUI
        f'-t{num_procs}', # Number of processors
        '-i', journal_file
    ]

    # Run simulation
    with open(log_file, 'w') as log:
        result = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )

    return result.returncode == 0

# Example usage
if __name__ == '__main__':
    success = run_fluent_batch(
        journal_file='simulation.jou',
        num_procs=8,
        log_file='fluent_run.log'
    )

    if success:
        print("Simulation completed successfully")
    else:
        print("Simulation failed - check log file")
```

### Parametric Studies
```python
import numpy as np
from ansys.fluent.core import launch_fluent

def parametric_study(velocities, case_file):
    """
    Run parametric study varying inlet velocity.

    Args:
        velocities: List of inlet velocities to test
        case_file: Base case file
    """
    results = {}

    for velocity in velocities:
        print(f"Running simulation with velocity = {velocity} m/s")

        # Launch Fluent
        solver = launch_fluent(precision='double', processor_count=4)

        # Read base case
        solver.tui.file.read_case(case_file)

        # Modify inlet velocity
        solver.tui.define.boundary_conditions.velocity_inlet(
            'inlet',
            'yes', 'no', 'yes', 'yes', 'no',
            velocity,
            'no', 0, 'no', 0, 'no', 300, 'no'
        )

        # Initialize and solve
        solver.tui.solve.initialize.initialize_flow()
        solver.tui.solve.iterate(500)

        # Extract results
        # (Add result extraction here)

        # Save case and data
        output_file = f'result_v{velocity}.cas.h5'
        solver.tui.file.write_case_data(output_file)

        solver.exit()

        results[velocity] = output_file

    return results

# Run parametric study
velocities = np.linspace(5, 25, 5)  # 5 to 25 m/s
results = parametric_study(velocities, 'base_case.cas.h5')
```

## Best Practices

1. **License Management**: Always check license availability before long batch runs
2. **Mesh Quality**: Verify mesh quality before solving (skewness, aspect ratio)
3. **Convergence Monitoring**: Set appropriate convergence criteria and monitor residuals
4. **File Management**: Use clear naming conventions for parametric studies
5. **Error Handling**: Implement robust error checking in automation scripts
6. **Documentation**: Comment journal files and Python scripts thoroughly
7. **Version Control**: Track journal files and scripts in version control
8. **Testing**: Test automation workflows on small cases before production runs
9. **Resource Management**: Monitor memory and disk usage for large simulations
10. **Reproducibility**: Document all settings for reproducible results

## Troubleshooting

### Common Issues
- **License checkout failures**: Check license server connectivity and availability
- **Mesh issues**: Review mesh quality metrics and refine problematic regions
- **Convergence problems**: Adjust under-relaxation factors, check BC consistency
- **Memory errors**: Reduce mesh size or use more nodes in parallel
- **File I/O errors**: Verify file paths and permissions

### Debugging Tips
- Run interactively first to validate workflow
- Check transcript files for TUI command syntax
- Use `-print` flag in journal files for debugging
- Monitor residual plots for convergence behavior
- Review warning and error messages in log files

## Additional Resources

- **ANSYS Documentation**: Help system within ANSYS products
- **PyAnsys Documentation**: https://docs.pyansys.com
- **ANSYS Learning Hub**: Online training and tutorials
- **ANSYS Community**: User forums and knowledge base
- **Fluent User's Guide**: Comprehensive theory and usage manual
- **Fluent TUI Reference**: Complete TUI command listing
