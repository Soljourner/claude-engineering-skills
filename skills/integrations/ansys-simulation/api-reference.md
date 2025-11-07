# ANSYS Fluent API Reference

Comprehensive reference for PyAnsys Fluent Python API, Journal file scripting, and TUI (Text User Interface) commands.

## PyAnsys Fluent API

### Session Management

#### Launch Fluent
```python
from ansys.fluent.core import launch_fluent

# Launch with default settings
solver = launch_fluent()

# Launch with specific options
solver = launch_fluent(
    precision='double',          # 'single' or 'double'
    processor_count=8,           # Number of processors
    mode='solver',              # 'solver' or 'meshing'
    dimension=3,                # 2 or 3
    show_gui=False,             # Show GUI (default: False)
    product_version='24.1'      # Specific version
)

# Launch meshing mode
meshing = launch_fluent(mode='meshing', processor_count=4)
```

#### Session Control
```python
# Exit cleanly
solver.exit()

# Check if session is alive
if solver.is_alive():
    print("Session is running")

# Get Fluent version
version = solver.get_fluent_version()

# Access TUI
tui = solver.tui

# Access settings (object-oriented API)
settings = solver.settings
```

### File Operations

```python
# Read case file
solver.tui.file.read_case('case.cas')
solver.tui.file.read_case('case.cas.h5')  # HDF5 format

# Read data file
solver.tui.file.read_data('data.dat')

# Read case and data together
solver.tui.file.read_case_data('case_data.cas.h5')

# Write case file
solver.tui.file.write_case('output.cas.h5')

# Write data file
solver.tui.file.write_data('output.dat.h5')

# Write case and data
solver.tui.file.write_case_data('output.cas.h5')

# Import geometry
solver.tui.file.import_.cad_geometry('geometry.step', 'yes')

# Export data
solver.tui.file.export.ascii(
    'results.csv',
    'surface-name',
    '()',
    'yes',
    'pressure',
    'velocity-magnitude',
    '()'
)
```

### Models and Physics

```python
# Energy equation
solver.tui.define.models.energy('yes')  # Enable
solver.tui.define.models.energy('no')   # Disable

# Turbulence models
solver.tui.define.models.viscous.laminar()
solver.tui.define.models.viscous.ke_standard('yes')
solver.tui.define.models.viscous.ke_realizable('yes')
solver.tui.define.models.viscous.kw_standard('yes')
solver.tui.define.models.viscous.kw_sst('yes')
solver.tui.define.models.viscous.spalart_allmaras('yes')
solver.tui.define.models.viscous.les('yes')
solver.tui.define.models.viscous.des('yes')

# Multiphase models
solver.tui.define.models.multiphase.model('vof', 'yes')
solver.tui.define.models.multiphase.model('mixture', 'yes')
solver.tui.define.models.multiphase.model('eulerian', 'yes')

# Radiation
solver.tui.define.models.radiation('p1')
solver.tui.define.models.radiation('do')
solver.tui.define.models.radiation('s2s')

# Species transport
solver.tui.define.models.species.species_transport('yes', 'yes')
```

### Materials

```python
# Create/modify material
solver.tui.define.materials.change_create(
    'water-liquid',     # Material name
    'water',            # Display name
    'yes',              # Density
    'constant',         # Density type
    998.2,              # Density value
    'yes',              # Viscosity
    'constant',         # Viscosity type
    0.001003,           # Viscosity value
    'yes',              # Specific heat
    'constant',
    4182,               # Cp value
    'yes',              # Thermal conductivity
    'constant',
    0.6                 # k value
)

# Copy material from database
solver.tui.define.materials.copy('fluid', 'air')
solver.tui.define.materials.copy('fluid', 'water-liquid')

# Delete material
solver.tui.define.materials.delete('material-name')
```

### Boundary Conditions

```python
# Velocity inlet
solver.tui.define.boundary_conditions.velocity_inlet(
    'inlet',           # Zone name
    'yes',             # Modify?
    'no',              # Set mixture values?
    'yes',             # Velocity magnitude?
    'yes',             # Normal to boundary?
    'no',              # Components?
    10.0,              # Velocity magnitude (m/s)
    'no',              # Turbulent intensity?
    0.05,
    'no',              # Hydraulic diameter?
    0.1,
    'no',              # Temperature?
    300,
    'no'               # More options?
)

# Mass flow inlet
solver.tui.define.boundary_conditions.mass_flow_inlet(
    'inlet',
    'yes',
    'no',
    'yes',
    'yes',
    'no',
    0.5,               # Mass flow rate (kg/s)
    'no',
    'no',
    'no',
    300,
    'no'
)

# Pressure inlet
solver.tui.define.boundary_conditions.pressure_inlet(
    'inlet',
    'yes',
    'no',
    101325,            # Total pressure (Pa)
    'yes',
    'no',
    'no',
    'no',
    300,
    'no'
)

# Pressure outlet
solver.tui.define.boundary_conditions.pressure_outlet(
    'outlet',
    'yes',
    'no',
    0,                 # Gauge pressure (Pa)
    'no',
    300,
    'no',
    'yes'
)

# Wall boundary
solver.tui.define.boundary_conditions.wall(
    'wall',
    'yes',
    'no',              # Shear boundary?
    'no',              # Velocity?
    0,
    'no',              # Rotation?
    0,
    'no',              # Heat flux?
    'yes',             # Temperature?
    'no',
    300
)

# Symmetry
solver.tui.define.boundary_conditions.symmetry('symmetry-plane', 'yes')

# Periodic boundaries
solver.tui.define.boundary_conditions.periodic(
    'periodic-1',
    'periodic-2',
    'yes',
    'yes',
    'translational',
    'yes',
    0, 1, 0            # Translation vector
)
```

### Cell Zones

```python
# Fluid zone
solver.tui.define.boundary_conditions.fluid(
    'fluid-zone',
    'yes',
    'water-liquid',    # Material name
    'no',              # Frame motion?
    'no',              # Porous zone?
    'no',              # Source terms?
    'no',
    0,
    'no',
    0,
    'no',
    0,
    'no',
    'no'
)

# Solid zone
solver.tui.define.boundary_conditions.solid(
    'solid-zone',
    'yes',
    'aluminum',        # Material
    'no',              # Source terms?
    'no'
)
```

### Solution Methods

```python
# Pressure-velocity coupling
solver.tui.solve.set.p_v_coupling(20)  # SIMPLE
solver.tui.solve.set.p_v_coupling(21)  # SIMPLEC
solver.tui.solve.set.p_v_coupling(22)  # PISO
solver.tui.solve.set.p_v_coupling(24)  # Coupled

# Discretization schemes
# Pressure: 10=Standard, 11=Linear, 12=Second-order, 14=PRESTO!
solver.tui.solve.set.discretization_scheme.pressure(12)

# Momentum: 0=First-order upwind, 1=Second-order upwind, 3=QUICK
solver.tui.solve.set.discretization_scheme.mom(1)

# Turbulence
solver.tui.solve.set.discretization_scheme.k(1)
solver.tui.solve.set.discretization_scheme.omega(1)
solver.tui.solve.set.discretization_scheme.epsilon(1)

# Energy
solver.tui.solve.set.discretization_scheme.temperature(1)

# Gradient evaluation
solver.tui.solve.set.gradient_scheme(0)  # Green-Gauss cell-based
solver.tui.solve.set.gradient_scheme(1)  # Green-Gauss node-based
solver.tui.solve.set.gradient_scheme(2)  # Least squares cell-based

# Under-relaxation factors
solver.tui.solve.set.under_relaxation.pressure(0.3)
solver.tui.solve.set.under_relaxation.density(1.0)
solver.tui.solve.set.under_relaxation.body_force(1.0)
solver.tui.solve.set.under_relaxation.mom(0.7)
solver.tui.solve.set.under_relaxation.k(0.8)
solver.tui.solve.set.under_relaxation.omega(0.8)
solver.tui.solve.set.under_relaxation.turb_viscosity(1.0)
solver.tui.solve.set.under_relaxation.temperature(1.0)
```

### Solution Initialization and Execution

```python
# Initialize from inlet
solver.tui.solve.initialize.compute_defaults.velocity_inlet('inlet')
solver.tui.solve.initialize.initialize_flow()

# Initialize with custom values
solver.tui.solve.initialize.set_defaults.pressure(101325)
solver.tui.solve.initialize.set_defaults.x_velocity(0)
solver.tui.solve.initialize.set_defaults.y_velocity(0)
solver.tui.solve.initialize.set_defaults.z_velocity(0)
solver.tui.solve.initialize.set_defaults.temperature(300)
solver.tui.solve.initialize.initialize_flow()

# Hybrid initialization
solver.tui.solve.initialize.hyb_initialization()

# Run iterations
solver.tui.solve.iterate(1000)

# Run until converged (with max iterations)
solver.tui.solve.iterate(5000)
```

### Monitors and Convergence

```python
# Residual monitors
solver.tui.solve.monitors.residual.plot('yes')
solver.tui.solve.monitors.residual.print('yes')

# Set convergence criteria
solver.tui.solve.monitors.residual.convergence_criteria(
    '1e-6',   # Continuity
    '1e-6',   # x-velocity
    '1e-6',   # y-velocity
    '1e-6',   # z-velocity
    '1e-6',   # Energy
    '1e-6',   # k
    '1e-6'    # omega/epsilon
)

# Surface monitors (force, mass flow, etc.)
solver.tui.solve.monitors.surface.set_monitor(
    'mass-flow-rate',
    'outlet',
    'Mass Flow Rate',
    'yes',            # Print?
    'yes',            # Plot?
    'yes',            # Write?
    'outlet_flow.out'
)
```

### Post-Processing

```python
# Create surfaces
solver.tui.surface.plane_surface(
    'x', 0,            # Normal direction and position
    'plane-x0'         # Surface name
)

solver.tui.surface.iso_surface(
    'pressure',        # Field variable
    'iso-p',           # Surface name
    'fluid-zone',      # Cell zone
    '()',
    101325             # Iso-value
)

# Report forces
solver.tui.report.forces.wall_forces(
    'yes',             # Print to console
    'wall',            # Wall zone(s)
    '()',
    'yes',             # Direction vector?
    'no',              # Absolute?
    'yes'              # Total forces?
)

# Report mass flow
solver.tui.report.surface_integrals.mass_flow_rate(
    'inlet', 'outlet', '()'
)

# Report volume integrals
solver.tui.report.volume_integrals.volume('fluid-zone', '()')
solver.tui.report.volume_integrals.volume_avg(
    'fluid-zone', '()', 'pressure'
)

# Export data
solver.tui.file.export.ascii(
    'surface_data.csv',
    'surface-name',
    '()',
    'yes',             # Export?
    'pressure',
    'velocity-magnitude',
    'temperature',
    '()'               # End list
)
```

## Journal File Syntax

Journal files contain TUI commands in text format. Each command follows the menu hierarchy.

### Basic Syntax

```scheme
; Comments start with semicolon
; Commands use forward slashes to navigate menus
; Prompts are answered in sequence

/file/read-case case.cas

; Multi-line commands
/define/boundary-conditions/velocity-inlet inlet yes no yes yes no 10 no 0.05 no 0.1 no 300 no

; List termination with empty parentheses
/solve/monitors/residual/convergence-criteria 1e-6 1e-6 1e-6 1e-6 ()
```

### Common Journal File Commands

#### File I/O
```scheme
; Read files
/file/read-case "case.cas.h5"
/file/read-data "data.dat.h5"
/file/read-case-data "case_data.cas.h5"

; Write files
/file/write-case "output.cas.h5"
/file/write-data "output.dat.h5"
/file/write-case-data "output.cas.h5"

; Import geometry
/file/import/cad-geometry "geometry.step" yes

; Export results
/file/export/ascii "results.csv" surface-1 () yes pressure velocity-magnitude ()
```

#### Models
```scheme
; Energy
/define/models/energy yes no no no no

; Turbulence - k-omega SST
/define/models/viscous kw-sst yes

; Turbulence - k-epsilon realizable
/define/models/viscous ke-realizable yes

; Multiphase - VOF
/define/models/multiphase/model vof yes yes
```

#### Materials
```scheme
; Copy from database
/define/materials/copy fluid water-liquid

; Create custom material
/define/materials/change-create water-liquid water yes constant 998.2 yes constant 0.001003 yes constant 4182 yes constant 0.6 no no no no no
```

#### Boundary Conditions
```scheme
; Velocity inlet
/define/boundary-conditions/velocity-inlet inlet yes no yes yes no 10 no 0.05 no 0.1 no 300 no

; Pressure outlet
/define/boundary-conditions/pressure-outlet outlet yes no 0 no 300 no yes

; Wall
/define/boundary-conditions/wall wall yes no no 0 no 0 no yes no 300
```

#### Solution Setup
```scheme
; Solution methods
/solve/set/p-v-coupling 24
/solve/set/discretization-scheme/pressure 12
/solve/set/discretization-scheme/mom 1
/solve/set/discretization-scheme/k 1
/solve/set/discretization-scheme/omega 1

; Under-relaxation
/solve/set/under-relaxation/pressure 0.3
/solve/set/under-relaxation/mom 0.7

; Initialize
/solve/initialize/initialize-flow

; Iterate
/solve/iterate 1000
```

#### Monitoring
```scheme
; Convergence criteria
/solve/monitors/residual/convergence-criteria 1e-6 1e-6 1e-6 1e-6 1e-6 ()

; Surface monitor
/solve/monitors/surface/set-monitor "Mass Flow Rate" outlet yes yes yes outlet_flow.out
```

#### Reporting
```scheme
; Forces
/report/forces/wall-forces yes wall () yes no yes

; Surface integrals
/report/surface-integrals/mass-flow-rate inlet outlet ()
```

### Control Flow in Journal Files

```scheme
; Exit Fluent
/exit yes

; Execute system command
/system "ls -l"

; Time delay (seconds)
/sleep 10

; Conditional execution (limited)
; Use scheme scripting for complex logic
(if (> 1 2) (ti-menu-load-string "/solve/iterate 100"))
```

## TUI Commands Reference

### Navigation
- `/` - Return to top-level menu
- `q` - Quit current menu level
- `?` - Display available commands
- Tab - Auto-complete command

### Key TUI Menus

```
/define           - Define models, materials, BCs, etc.
  /models         - Physical models
  /materials      - Material properties
  /boundary-conditions - Boundary and cell zone conditions
  /operating-conditions - Reference values, gravity
  /user-defined   - UDFs and custom settings

/solve            - Solution controls and execution
  /set            - Solution methods and controls
  /initialize     - Flow field initialization
  /iterate        - Run iterations
  /monitors       - Convergence and surface monitors

/display          - Visualization (GUI mode)
  /contour        - Contour plots
  /vector         - Vector plots
  /pathline       - Particle tracking

/surface          - Create surfaces for post-processing
  /plane-surface  - Planar surfaces
  /iso-surface    - Iso-value surfaces

/report           - Calculate and report results
  /forces         - Force and moment calculations
  /surface-integrals - Area-weighted averaging
  /volume-integrals - Volume-weighted averaging
  /fluxes         - Mass and energy fluxes

/file             - File operations
  /read-case      - Read case file
  /write-case     - Write case file
  /import         - Import geometry/mesh
  /export         - Export data

/mesh             - Mesh operations
  /check          - Check mesh quality
  /scale          - Scale mesh
  /smooth         - Smooth mesh

/parallel         - Parallel computing controls
```

## Examples for Common Tasks

### Example 1: Complete Flow Setup

```python
from ansys.fluent.core import launch_fluent

# Launch Fluent
solver = launch_fluent(precision='double', processor_count=4)

# Read mesh
solver.tui.file.read_case('mesh.cas.h5')

# Enable models
solver.tui.define.models.energy('yes')
solver.tui.define.models.viscous.kw_sst('yes')

# Set materials
solver.tui.define.materials.copy('fluid', 'water-liquid')

# Boundary conditions
solver.tui.define.boundary_conditions.velocity_inlet(
    'inlet', 'yes', 'no', 'yes', 'yes', 'no', 5.0,
    'no', 0.05, 'no', 0.1, 'no', 300, 'no'
)
solver.tui.define.boundary_conditions.pressure_outlet(
    'outlet', 'yes', 'no', 0, 'no', 300, 'no', 'yes'
)
solver.tui.define.boundary_conditions.wall(
    'walls', 'yes', 'no', 'no', 0, 'no', 0, 'no', 'yes', 'no', 300
)

# Solution methods
solver.tui.solve.set.p_v_coupling(24)
solver.tui.solve.set.discretization_scheme.pressure(12)
solver.tui.solve.set.discretization_scheme.mom(1)

# Initialize and solve
solver.tui.solve.initialize.initialize_flow()
solver.tui.solve.iterate(1000)

# Save results
solver.tui.file.write_case_data('result.cas.h5')

solver.exit()
```

### Example 2: Automated Post-Processing

```python
def extract_results(case_data_file, output_csv):
    """Extract pressure and velocity data from surfaces."""
    solver = launch_fluent(precision='double')
    solver.tui.file.read_case_data(case_data_file)

    # Create plane surface at x=0
    solver.tui.surface.plane_surface('x', 0, 'plane-x0')

    # Export data
    solver.tui.file.export.ascii(
        output_csv,
        'plane-x0',
        '()',
        'yes',
        'pressure',
        'velocity-magnitude',
        'x-velocity',
        'y-velocity',
        'z-velocity',
        'temperature',
        '()'
    )

    # Calculate forces on walls
    solver.tui.report.forces.wall_forces(
        'yes', 'walls', '()', 'yes', 'no', 'yes'
    )

    solver.exit()
    print(f"Results exported to {output_csv}")

# Usage
extract_results('simulation.cas.h5', 'plane_data.csv')
```

### Example 3: Parametric Temperature Study

```python
import numpy as np
from ansys.fluent.core import launch_fluent

def temperature_study(base_case, temperatures, output_dir):
    """Run parametric study with varying inlet temperature."""
    results = {}

    for temp in temperatures:
        print(f"Running case with T = {temp} K")

        solver = launch_fluent(precision='double', processor_count=4)
        solver.tui.file.read_case(base_case)

        # Modify inlet temperature
        solver.tui.define.boundary_conditions.velocity_inlet(
            'inlet', 'yes', 'no', 'yes', 'yes', 'no', 5.0,
            'no', 0.05, 'no', 0.1, 'no', temp, 'no'
        )

        # Solve
        solver.tui.solve.initialize.initialize_flow()
        solver.tui.solve.iterate(500)

        # Save results
        output_file = f"{output_dir}/result_T{temp}.cas.h5"
        solver.tui.file.write_case_data(output_file)

        # Extract outlet temperature
        solver.tui.report.surface_integrals.area_weighted_avg(
            'outlet', '()', 'temperature'
        )

        results[temp] = output_file
        solver.exit()

    return results

# Run study
temps = np.linspace(300, 400, 6)
results = temperature_study('base.cas.h5', temps, './results')
```

### Example 4: Journal File Template

```scheme
; ============================================
; Automated CFD Simulation Journal File
; ============================================

; Read case file
/file/read-case "mesh.cas.h5"

; Enable physical models
/define/models/energy yes no no no no
/define/models/viscous kw-sst yes

; Set material
/define/materials/copy fluid water-liquid

; Boundary conditions
/define/boundary-conditions/velocity-inlet inlet yes no yes yes no 10 no 0.05 no 0.1 no 300 no
/define/boundary-conditions/pressure-outlet outlet yes no 0 no 300 no yes
/define/boundary-conditions/wall walls yes no no 0 no 0 no yes no 300

; Solution methods
/solve/set/p-v-coupling 24
/solve/set/discretization-scheme/pressure 12
/solve/set/discretization-scheme/mom 1
/solve/set/discretization-scheme/k 1
/solve/set/discretization-scheme/omega 1

; Under-relaxation
/solve/set/under-relaxation/pressure 0.3
/solve/set/under-relaxation/mom 0.7

; Convergence criteria
/solve/monitors/residual/convergence-criteria 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 ()

; Initialize
/solve/initialize/initialize-flow

; Solve
/solve/iterate 1000

; Save results
/file/write-case-data "result.cas.h5"

; Report forces
/report/forces/wall-forces yes walls () yes no yes

; Exit
/exit yes
```

### Example 5: Batch Processing Multiple Cases

```python
from pathlib import Path
from ansys.fluent.core import launch_fluent

def batch_process_cases(case_dir, journal_template):
    """
    Process multiple case files with same journal file.

    Args:
        case_dir: Directory containing .cas files
        journal_template: Journal file to apply
    """
    case_files = Path(case_dir).glob('*.cas')

    for case_file in case_files:
        print(f"Processing {case_file.name}")

        solver = launch_fluent(precision='double', processor_count=4)

        # Read case
        solver.tui.file.read_case(str(case_file))

        # Execute journal commands
        with open(journal_template, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith(';'):
                    solver.tui.execute_command(line.strip())

        # Save with _solved suffix
        output_name = case_file.stem + '_solved.cas.h5'
        output_path = case_file.parent / output_name
        solver.tui.file.write_case_data(str(output_path))

        solver.exit()
        print(f"Completed: {output_name}")

# Usage
batch_process_cases('./cases', 'solve_template.jou')
```

## Additional Notes

### HDF5 File Format
- Recommended format: `.cas.h5` and `.dat.h5`
- Benefits: Faster I/O, smaller file size, parallel I/O support
- Use for all modern workflows

### Error Handling
```python
from ansys.fluent.core import launch_fluent

try:
    solver = launch_fluent(precision='double', processor_count=4)
    solver.tui.file.read_case('case.cas.h5')
    solver.tui.solve.iterate(1000)
except Exception as e:
    print(f"Error occurred: {e}")
finally:
    if solver.is_alive():
        solver.exit()
```

### Performance Tips
- Use double precision for better accuracy
- Enable parallel processing with `processor_count`
- Use coupled solver for pressure-velocity coupling when stable
- Increase under-relaxation factors gradually for faster convergence
- Monitor residuals to detect convergence issues early
