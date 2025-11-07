# Engineering Software Integrations

This document provides comprehensive guides for integrating with popular engineering software tools, including setup, scripting interfaces, authentication, and automation examples.

## Table of Contents

1. [OpenFOAM](#openfoam)
2. [ANSYS (Fluent, Workbench)](#ansys)
3. [SolidWorks](#solidworks)
4. [COMSOL Multiphysics](#comsol-multiphysics)
5. [MATLAB](#matlab)
6. [Star-CCM+](#star-ccm)

---

## OpenFOAM

### Overview and Capabilities

OpenFOAM (Open Field Operation and Manipulation) is an open-source computational fluid dynamics (CFD) toolbox written in C++. It provides:

- Incompressible and compressible flow solvers
- Multiphase flow simulation
- Turbulence modeling (RANS, LES, DES)
- Heat transfer and buoyancy-driven flows
- Particle tracking and spray modeling
- Mesh generation (blockMesh, snappyHexMesh)
- Post-processing with ParaView integration

### Setup Instructions

**Linux/Ubuntu Installation:**

```bash
# Add OpenFOAM repository
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.openfoam.org/ubuntu

# Install OpenFOAM
sudo apt-get update
sudo apt-get install openfoam11

# Source the OpenFOAM environment
source /opt/openfoam11/etc/bashrc

# Add to .bashrc for persistent setup
echo "source /opt/openfoam11/etc/bashrc" >> ~/.bashrc
```

**Docker Installation:**

```bash
docker pull openfoam/openfoam11-paraview510
docker run -it -v $PWD:/workspace openfoam/openfoam11-paraview510
```

### API/Scripting Interface

OpenFOAM uses a dictionary-based configuration system and can be scripted using:

**1. Python with PyFoam:**

```bash
pip install PyFoam
```

**2. Python with foamlib:**

```bash
pip install foamlib
```

**3. Shell scripting with OpenFOAM utilities**

### Authentication and Licensing

OpenFOAM is open-source under the GPL license. No authentication or license server required.

- Free for commercial and academic use
- No restrictions on distribution or modification
- Community support via forums and CFD Online

### Example Automation Scripts

**Python Script for Case Setup and Execution:**

```python
#!/usr/bin/env python3
"""
OpenFOAM cavity flow automation using PyFoam
"""

import os
import subprocess
from pathlib import Path

class OpenFOAMCase:
    def __init__(self, case_dir):
        self.case_dir = Path(case_dir)

    def create_case_structure(self):
        """Create standard OpenFOAM case directory structure"""
        dirs = ['0', 'constant', 'system']
        for d in dirs:
            (self.case_dir / d).mkdir(parents=True, exist_ok=True)
        print(f"Created case structure in {self.case_dir}")

    def write_blockMeshDict(self, nx, ny, nz, length, width, height):
        """Generate blockMeshDict for rectangular domain"""
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
    (0 0 0)
    ({length} 0 0)
    ({length} {width} 0)
    (0 {width} 0)
    (0 0 {height})
    ({length} 0 {height})
    ({length} {width} {height})
    (0 {width} {height})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    movingWall
    {{
        type wall;
        faces
        (
            (3 7 6 2)
        );
    }}
    fixedWalls
    {{
        type wall;
        faces
        (
            (0 4 7 3)
            (2 6 5 1)
            (1 5 4 0)
        );
    }}
    frontAndBack
    {{
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }}
);

mergePatchPairs
(
);

// ************************************************************************* //
"""
        file_path = self.case_dir / 'system' / 'blockMeshDict'
        file_path.write_text(content)
        print(f"Written blockMeshDict")

    def write_controlDict(self, solver, end_time, delta_t, write_interval):
        """Generate controlDict for simulation control"""
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     {solver};

startFrom       startTime;

startTime       0;

stopAt          endTime;

endTime         {end_time};

deltaT          {delta_t};

writeControl    timeStep;

writeInterval   {write_interval};

purgeWrite      0;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""
        file_path = self.case_dir / 'system' / 'controlDict'
        file_path.write_text(content)
        print(f"Written controlDict")

    def run_command(self, command):
        """Run OpenFOAM command in case directory"""
        result = subprocess.run(
            command,
            cwd=self.case_dir,
            shell=True,
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"Error: {result.stderr}")
            raise RuntimeError(f"Command failed: {command}")

        print(result.stdout)
        return result

    def generate_mesh(self):
        """Run blockMesh to generate mesh"""
        print("Generating mesh...")
        self.run_command("blockMesh")

    def run_simulation(self, solver="icoFoam"):
        """Run the solver"""
        print(f"Running {solver}...")
        self.run_command(solver)

    def post_process(self, function_object=None):
        """Run post-processing"""
        if function_object:
            print(f"Post-processing with {function_object}...")
            self.run_command(f"postProcess -func {function_object}")

    def extract_forces(self):
        """Extract forces from simulation"""
        forces_file = self.case_dir / 'postProcessing' / 'forces' / '0' / 'force.dat'
        if forces_file.exists():
            data = []
            with open(forces_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        data.append([float(x) for x in line.split()])
            return data
        return None


# Example usage
if __name__ == "__main__":
    # Create and run a cavity case
    case = OpenFOAMCase("cavity_example")

    # Setup case
    case.create_case_structure()
    case.write_blockMeshDict(nx=20, ny=20, nz=1,
                            length=0.1, width=0.1, height=0.01)
    case.write_controlDict(solver="icoFoam", end_time=0.5,
                          delta_t=0.005, write_interval=20)

    # Run simulation
    case.generate_mesh()
    case.run_simulation("icoFoam")

    print("Simulation completed successfully!")
```

**Parametric Study Script:**

```python
#!/usr/bin/env python3
"""
OpenFOAM parametric study automation
"""

import numpy as np
import pandas as pd
from pathlib import Path
import shutil

class ParametricStudy:
    def __init__(self, base_case, study_dir):
        self.base_case = Path(base_case)
        self.study_dir = Path(study_dir)
        self.results = []

    def clone_case(self, case_name):
        """Clone base case for parametric study"""
        new_case = self.study_dir / case_name
        if new_case.exists():
            shutil.rmtree(new_case)
        shutil.copytree(self.base_case, new_case)
        return new_case

    def modify_parameter(self, case_dir, dict_file, parameter, value):
        """Modify a parameter in OpenFOAM dictionary"""
        file_path = case_dir / dict_file
        content = file_path.read_text()

        # Simple parameter replacement (for basic cases)
        # For complex modifications, use PyFoam's parser
        modified = content.replace(f"{parameter}", str(value))
        file_path.write_text(modified)

    def run_study(self, parameter_name, parameter_values):
        """Run parametric study"""
        for i, value in enumerate(parameter_values):
            case_name = f"case_{parameter_name}_{value}"
            print(f"\n{'='*60}")
            print(f"Running case {i+1}/{len(parameter_values)}: {case_name}")
            print(f"{parameter_name} = {value}")
            print(f"{'='*60}")

            # Clone and setup case
            case_dir = self.clone_case(case_name)
            case = OpenFOAMCase(case_dir)

            # Modify parameter
            self.modify_parameter(case_dir, 'system/controlDict',
                                 'PARAMETER_PLACEHOLDER', value)

            # Run simulation
            try:
                case.run_simulation()

                # Extract results
                result = {
                    'case': case_name,
                    parameter_name: value,
                    'converged': True
                }
                self.results.append(result)

            except Exception as e:
                print(f"Error in case {case_name}: {e}")
                self.results.append({
                    'case': case_name,
                    parameter_name: value,
                    'converged': False
                })

        return pd.DataFrame(self.results)

# Example usage
if __name__ == "__main__":
    study = ParametricStudy("base_case", "parametric_study")
    reynolds_numbers = [100, 500, 1000, 5000, 10000]
    results_df = study.run_study("Reynolds", reynolds_numbers)
    results_df.to_csv("parametric_results.csv", index=False)
```

### File I/O and Data Exchange

**Reading OpenFOAM Data with Python:**

```python
import numpy as np
from pathlib import Path

def read_openfoam_scalar_field(file_path):
    """Read OpenFOAM scalar field"""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Find internal field data
    in_data = False
    data = []

    for line in lines:
        if 'internalField' in line and 'uniform' in line:
            # Uniform field
            value = float(line.split()[-1].strip(';'))
            return value
        elif 'internalField' in line:
            in_data = True
            continue

        if in_data:
            if line.strip().startswith('('):
                continue
            elif line.strip().startswith(')'):
                break
            else:
                try:
                    data.append(float(line.strip()))
                except:
                    continue

    return np.array(data)

def read_openfoam_vector_field(file_path):
    """Read OpenFOAM vector field"""
    with open(file_path, 'r') as f:
        lines = f.readlines()

    in_data = False
    data = []

    for line in lines:
        if 'internalField' in line:
            in_data = True
            continue

        if in_data:
            if line.strip().startswith('(') and len(line.strip()) > 1:
                # Vector data line
                vector_str = line.strip().strip('(').strip(')').strip(',')
                try:
                    vector = [float(x) for x in vector_str.split()]
                    if len(vector) == 3:
                        data.append(vector)
                except:
                    continue
            elif line.strip() == ')':
                break

    return np.array(data)

def write_openfoam_scalar_field(file_path, field_name, data,
                               object_class="volScalarField"):
    """Write OpenFOAM scalar field"""
    header = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {object_class};
    object      {field_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   nonuniform List<scalar>
{len(data)}
(
"""

    footer = """)
;

boundaryField
{{
}}

// ************************************************************************* //
"""

    with open(file_path, 'w') as f:
        f.write(header)
        for value in data:
            f.write(f"{value}\n")
        f.write(footer)

# Example: Read pressure field
pressure = read_openfoam_scalar_field("0/p")
velocity = read_openfoam_vector_field("0/U")
```

**VTK Export for ParaView:**

```python
def export_to_vtk(case_dir, times=None):
    """Convert OpenFOAM data to VTK format"""
    import subprocess

    cmd = ["foamToVTK"]
    if times:
        cmd.extend(["-time", ",".join(map(str, times))])

    subprocess.run(cmd, cwd=case_dir)
    print(f"VTK files exported to {case_dir}/VTK")
```

### Best Practices

1. **Case Organization:**
   - Use clear naming conventions for cases
   - Document case setup in a README
   - Use version control (git) for case files
   - Keep base cases separate from parametric studies

2. **Mesh Quality:**
   - Always check mesh quality with `checkMesh`
   - Target non-orthogonality < 70°
   - Keep aspect ratios reasonable (< 100)
   - Use appropriate refinement in critical regions

3. **Solver Settings:**
   - Start with relaxation for steady cases
   - Use appropriate turbulence models for flow regime
   - Monitor residuals and solution convergence
   - Implement adaptive time stepping for transient cases

4. **Performance:**
   - Use parallel decomposition for large cases
   - Profile simulations to identify bottlenecks
   - Use binary format for large fields
   - Clean old time directories to save space

5. **Automation:**
   - Use Allrun scripts for reproducibility
   - Implement proper error handling
   - Log all outputs for debugging
   - Validate results automatically

---

## ANSYS

### Overview and Capabilities

ANSYS provides comprehensive engineering simulation tools:

**ANSYS Fluent (CFD):**
- Advanced turbulence modeling
- Multiphase flows
- Combustion and reacting flows
- Heat transfer with radiation
- Moving/deforming meshes

**ANSYS Workbench (FEA):**
- Structural analysis (static, dynamic, thermal)
- Modal and harmonic analysis
- Fatigue and fracture mechanics
- Coupled multiphysics simulations

### Setup Instructions

**Installation:**

1. Download ANSYS installation package from ANSYS Customer Portal
2. Run the installer:
   ```bash
   # Linux
   sudo ./INSTALL

   # Windows
   setup.exe
   ```
3. Configure license server
4. Set environment variables:
   ```bash
   # Linux
   export ANSYS_ROOT=/ansys_inc/v242
   export PATH=$ANSYS_ROOT/fluent/bin:$PATH
   export PATH=$ANSYS_ROOT/Framework/bin/Linux64:$PATH
   ```

**Python Environment Setup:**

```bash
# ANSYS provides its own Python interpreter
# Typically located at:
# Linux: /ansys_inc/v242/commonfiles/CPython/3_10/linx64/Release/python
# Windows: C:\Program Files\ANSYS Inc\v242\commonfiles\CPython\3_10\winx64\Release\python\python.exe

# Or use PyAnsys packages with your Python
pip install ansys-fluent-core
pip install ansys-mapdl-core
pip install ansys-dpf-core
```

### API/Scripting Interface

ANSYS supports multiple scripting interfaces:

**1. PyFluent (Python API for Fluent):**

```python
import ansys.fluent.core as pyfluent

# Launch Fluent
solver = pyfluent.launch_fluent(
    precision='double',
    processor_count=4,
    mode='solver'
)

# Or connect to existing session
# solver = pyfluent.connect_to_fluent(ip='localhost', port=12345)
```

**2. ANSYS APDL (MAPDL) Scripting:**

```python
from ansys.mapdl.core import launch_mapdl

# Launch MAPDL
mapdl = launch_mapdl()
```

**3. Journal Files (TUI Commands):**

```scheme
; Fluent journal file example
/file/read-case "case.cas"
/solve/initialize/initialize-flow
/solve/iterate 1000
/file/write-data "result.dat"
/exit
```

**4. ACT (ANSYS Customization Toolkit) - Python/IronPython**

**5. Workbench Scripting (Python/JavaScript)**

### Authentication and Licensing

ANSYS uses FlexLM license management:

**License Server Configuration:**

```bash
# Set license server environment variable
# Linux
export ANSYSLMD_LICENSE_FILE=1055@license-server.company.com
export ANSYSLI_SERVERS=2325@license-server.company.com

# Windows (System Environment Variables)
ANSYSLMD_LICENSE_FILE=1055@license-server.company.com
```

**Check License Status:**

```bash
# Linux
$ANSYS_ROOT/licensingclient/linx64/ansysli_util -status

# Windows
"C:\Program Files\ANSYS Inc\v242\licensingclient\winx64\ansysli_util.exe" -status
```

**License Checkout in Python:**

```python
import os

# Set license server before importing ANSYS modules
os.environ['ANSYSLMD_LICENSE_FILE'] = '1055@license-server.company.com'
os.environ['ANSYSLI_SERVERS'] = '2325@license-server.company.com'

# Now import and use ANSYS modules
import ansys.fluent.core as pyfluent
```

### Example Automation Scripts

**Fluent: Complete CFD Simulation Automation:**

```python
#!/usr/bin/env python3
"""
ANSYS Fluent automation for pipe flow simulation
"""

import ansys.fluent.core as pyfluent
from ansys.fluent.core import examples
import numpy as np
import pandas as pd

class FluentAutomation:
    def __init__(self, precision='double', processor_count=4):
        """Initialize Fluent session"""
        self.solver = None
        self.precision = precision
        self.processor_count = processor_count

    def launch(self):
        """Launch Fluent solver"""
        print("Launching ANSYS Fluent...")
        self.solver = pyfluent.launch_fluent(
            precision=self.precision,
            processor_count=self.processor_count,
            mode='solver',
            show_gui=False
        )
        print(f"Fluent launched with {self.processor_count} processors")
        return self.solver

    def import_mesh(self, mesh_file):
        """Import mesh file"""
        print(f"Importing mesh: {mesh_file}")
        self.solver.file.read_case(file_name=mesh_file)

        # Check mesh quality
        self.solver.mesh.check()

        # Get mesh statistics
        stats = self.solver.mesh.get_stats()
        print(f"Mesh statistics: {stats}")

    def setup_physics(self, solver_type='pressure-based',
                     viscous_model='k-epsilon'):
        """Setup physics models"""
        print("Setting up physics models...")

        # Set solver type
        self.solver.setup.general.solver.type = solver_type

        # Enable energy equation if needed
        self.solver.setup.models.energy.enabled = True

        # Set viscous model
        if viscous_model == 'k-epsilon':
            self.solver.setup.models.viscous.model = 'k-epsilon'
            self.solver.setup.models.viscous.k_epsilon_model = 'standard'
        elif viscous_model == 'k-omega':
            self.solver.setup.models.viscous.model = 'k-omega'
        elif viscous_model == 'laminar':
            self.solver.setup.models.viscous.model = 'laminar'

        print(f"Viscous model set to: {viscous_model}")

    def setup_materials(self, material_name='water-liquid'):
        """Setup material properties"""
        print(f"Setting up material: {material_name}")

        # Get material from database
        material = self.solver.setup.materials.database.copy_by_name(
            type='fluid',
            name=material_name
        )

        # Assign material to cell zones
        for zone in self.solver.setup.cell_zone_conditions.get_zones():
            zone.material = material_name

    def setup_boundary_conditions(self, bc_config):
        """Setup boundary conditions

        bc_config: dict with boundary zone names and conditions
        Example: {
            'inlet': {'type': 'velocity-inlet', 'velocity': 1.0},
            'outlet': {'type': 'pressure-outlet', 'pressure': 0},
            'wall': {'type': 'wall', 'wall_motion': 'stationary'}
        }
        """
        print("Setting up boundary conditions...")

        for zone_name, config in bc_config.items():
            zone = self.solver.setup.boundary_conditions[zone_name]

            if config['type'] == 'velocity-inlet':
                zone.type = 'velocity-inlet'
                zone.momentum.velocity = config.get('velocity', 1.0)

            elif config['type'] == 'pressure-outlet':
                zone.type = 'pressure-outlet'
                zone.pressure = config.get('pressure', 0)

            elif config['type'] == 'wall':
                zone.type = 'wall'
                if 'temperature' in config:
                    zone.thermal.temperature = config['temperature']

            print(f"  {zone_name}: {config['type']}")

    def setup_solution_methods(self):
        """Configure solution methods"""
        print("Configuring solution methods...")

        # Set pressure-velocity coupling
        self.solver.solution.methods.p_v_coupling.flow_scheme = 'SIMPLE'

        # Set discretization schemes
        self.solver.solution.methods.discretization.pressure = 'second-order'
        self.solver.solution.methods.discretization.momentum = 'second-order-upwind'
        self.solver.solution.methods.discretization.turbulent_kinetic_energy = 'second-order-upwind'
        self.solver.solution.methods.discretization.turbulent_dissipation_rate = 'second-order-upwind'

    def initialize_solution(self, method='hybrid'):
        """Initialize flow field"""
        print(f"Initializing solution using {method} method...")

        if method == 'hybrid':
            self.solver.solution.initialization.hybrid_initialize()
        else:
            self.solver.solution.initialization.standard_initialize()

    def setup_monitors(self, monitor_config):
        """Setup residual and report monitors

        monitor_config: dict with monitor specifications
        """
        print("Setting up solution monitors...")

        # Residual convergence criteria
        for variable, criterion in monitor_config.get('residuals', {}).items():
            self.solver.solution.monitor.residual.equations[variable] = criterion

        # Force monitors
        if 'forces' in monitor_config:
            for name, config in monitor_config['forces'].items():
                self.solver.solution.report_definitions.surface[name] = {
                    'report_type': 'surface-integrals',
                    'field': 'pressure-coefficient',
                    'surface_names': config['surfaces']
                }

    def run_calculation(self, iterations=1000):
        """Run the calculation"""
        print(f"Running calculation for {iterations} iterations...")

        self.solver.solution.run_calculation.iterate(
            number_of_iterations=iterations
        )

        # Check convergence
        print("Calculation complete")

    def extract_results(self):
        """Extract and return results"""
        print("Extracting results...")

        results = {}

        # Get residuals
        residuals = self.solver.solution.monitor.residual.get_values()
        results['residuals'] = residuals

        # Get report values
        reports = self.solver.solution.report_definitions.get_values()
        results['reports'] = reports

        return results

    def save_results(self, case_file, data_file):
        """Save case and data files"""
        print(f"Saving results to {data_file}")
        self.solver.file.write_case_data(case_file_name=case_file,
                                          data_file_name=data_file)

    def export_data(self, surface_name, filename, variables):
        """Export surface data"""
        print(f"Exporting data from {surface_name} to {filename}")

        self.solver.file.export.ascii(
            file_name=filename,
            surface_names=[surface_name],
            variables=variables
        )

    def close(self):
        """Close Fluent session"""
        if self.solver:
            print("Closing Fluent...")
            self.solver.exit()


# Example usage
if __name__ == "__main__":
    # Initialize automation
    fluent = FluentAutomation(processor_count=4)
    fluent.launch()

    try:
        # Import mesh
        fluent.import_mesh("pipe_flow.msh")

        # Setup physics
        fluent.setup_physics(viscous_model='k-epsilon')

        # Setup materials
        fluent.setup_materials('water-liquid')

        # Setup boundary conditions
        bc_config = {
            'inlet': {'type': 'velocity-inlet', 'velocity': 2.0},
            'outlet': {'type': 'pressure-outlet', 'pressure': 0},
            'wall': {'type': 'wall'}
        }
        fluent.setup_boundary_conditions(bc_config)

        # Configure solution
        fluent.setup_solution_methods()
        fluent.initialize_solution()

        # Setup monitors
        monitor_config = {
            'residuals': {
                'continuity': 1e-4,
                'x-velocity': 1e-4,
                'y-velocity': 1e-4,
                'k': 1e-4,
                'epsilon': 1e-4
            }
        }
        fluent.setup_monitors(monitor_config)

        # Run calculation
        fluent.run_calculation(iterations=1000)

        # Extract and save results
        results = fluent.extract_results()
        fluent.save_results("pipe_flow.cas", "pipe_flow.dat")

        # Export data
        fluent.export_data("outlet", "outlet_data.csv",
                          ['velocity-magnitude', 'pressure'])

        print("Simulation completed successfully!")

    finally:
        fluent.close()
```

**MAPDL: Structural Analysis Automation:**

```python
#!/usr/bin/env python3
"""
ANSYS MAPDL structural analysis automation
"""

from ansys.mapdl.core import launch_mapdl
import numpy as np

class MAPDLAutomation:
    def __init__(self, run_location=None):
        """Initialize MAPDL session"""
        self.mapdl = None
        self.run_location = run_location

    def launch(self, nproc=2):
        """Launch MAPDL"""
        print("Launching ANSYS MAPDL...")
        self.mapdl = launch_mapdl(
            nproc=nproc,
            run_location=self.run_location,
            override=True
        )
        print(f"MAPDL launched with {nproc} processors")
        print(f"Working directory: {self.mapdl.directory}")
        return self.mapdl

    def create_beam_model(self, length=1.0, width=0.1, height=0.05):
        """Create a simple cantilever beam model"""
        print("Creating beam model...")

        # Clear database
        self.mapdl.clear()

        # Set title
        self.mapdl.title("Cantilever Beam Analysis")

        # Define element type
        self.mapdl.prep7()  # Enter preprocessor
        self.mapdl.et(1, "BEAM188")  # 3D beam element

        # Define material properties
        # Steel properties
        E = 2e11  # Young's modulus (Pa)
        nu = 0.3  # Poisson's ratio
        rho = 7850  # Density (kg/m3)

        self.mapdl.mp("EX", 1, E)
        self.mapdl.mp("PRXY", 1, nu)
        self.mapdl.mp("DENS", 1, rho)

        # Define section properties
        self.mapdl.sectype(1, "BEAM", "RECT")
        self.mapdl.secdata(width, height)

        # Create nodes
        n_nodes = 21
        for i in range(n_nodes):
            x = i * length / (n_nodes - 1)
            self.mapdl.n(i + 1, x, 0, 0)

        # Create elements
        self.mapdl.secnum(1)
        self.mapdl.mat(1)
        for i in range(n_nodes - 1):
            self.mapdl.e(i + 1, i + 2)

        # Apply boundary conditions
        # Fix all DOFs at left end
        self.mapdl.d(1, "ALL", 0)

        # Apply load at right end
        force = 1000  # N
        self.mapdl.f(n_nodes, "FY", -force)

        print(f"Model created: {n_nodes} nodes, {n_nodes-1} elements")

        return n_nodes

    def solve_static(self):
        """Solve static analysis"""
        print("Solving static analysis...")

        # Enter solution processor
        self.mapdl.slashsolu()

        # Set analysis type
        self.mapdl.antype("STATIC")

        # Solve
        output = self.mapdl.solve()

        # Check solution status
        self.mapdl.finish()

        print("Static solution complete")
        return output

    def solve_modal(self, n_modes=5):
        """Solve modal analysis"""
        print(f"Solving modal analysis for {n_modes} modes...")

        # Enter solution processor
        self.mapdl.slashsolu()

        # Set analysis type
        self.mapdl.antype("MODAL")

        # Set mode extraction method
        self.mapdl.modopt("LANB", n_modes)

        # Solve
        output = self.mapdl.solve()

        # Check solution status
        self.mapdl.finish()

        print("Modal solution complete")
        return output

    def postprocess_static(self):
        """Post-process static results"""
        print("Post-processing static results...")

        # Enter post-processor
        self.mapdl.post1()

        # Read results
        self.mapdl.set(1, 1)  # Read first load step

        # Get nodal displacements
        displacements = self.mapdl.post_processing.nodal_displacement('ALL')

        # Get nodal stresses
        stresses = self.mapdl.post_processing.nodal_stress('ALL')

        # Get maximum values
        max_disp = np.max(np.abs(displacements))
        max_stress = np.max(np.abs(stresses))

        results = {
            'max_displacement': max_disp,
            'max_stress': max_stress,
            'displacements': displacements,
            'stresses': stresses
        }

        print(f"Maximum displacement: {max_disp:.6e} m")
        print(f"Maximum stress: {max_stress:.6e} Pa")

        return results

    def postprocess_modal(self, n_modes=5):
        """Post-process modal results"""
        print("Post-processing modal results...")

        # Enter post-processor
        self.mapdl.post1()

        # Extract frequencies
        frequencies = []
        for mode in range(1, n_modes + 1):
            self.mapdl.set(1, mode)
            freq = self.mapdl.post_processing.frequency
            frequencies.append(freq)
            print(f"Mode {mode}: {freq:.4f} Hz")

        return frequencies

    def plot_deformation(self, filename="deformation.png"):
        """Plot deformed shape"""
        print(f"Plotting deformation to {filename}")
        self.mapdl.post1()
        self.mapdl.set(1, 1)
        self.mapdl.plnsol("U", "SUM")
        # Save plot
        self.mapdl.show("PNG")
        self.mapdl.savfig(filename)

    def export_results(self, filename="results.txt"):
        """Export results to text file"""
        print(f"Exporting results to {filename}")
        self.mapdl.post1()
        self.mapdl.prnsol("U", "COMP")
        self.mapdl.prnsol("S", "COMP")

    def close(self):
        """Exit MAPDL"""
        if self.mapdl:
            print("Closing MAPDL...")
            self.mapdl.exit()


# Example usage
if __name__ == "__main__":
    # Initialize automation
    mapdl_auto = MAPDLAutomation()
    mapdl_auto.launch(nproc=2)

    try:
        # Create model
        mapdl_auto.create_beam_model(length=1.0, width=0.1, height=0.05)

        # Static analysis
        mapdl_auto.solve_static()
        static_results = mapdl_auto.postprocess_static()

        # Modal analysis
        mapdl_auto.solve_modal(n_modes=5)
        frequencies = mapdl_auto.postprocess_modal(n_modes=5)

        # Plot results
        mapdl_auto.plot_deformation()

        print("\nAnalysis completed successfully!")
        print(f"Natural frequencies: {frequencies}")

    finally:
        mapdl_auto.close()
```

### File I/O and Data Exchange

**Reading ANSYS Results Files:**

```python
from ansys.dpf import core as dpf
from ansys.dpf.core import examples

# Read result file
model = dpf.Model("file.rst")

# Get mesh
mesh = model.metadata.meshed_region

# Get stress results
stress_op = dpf.operators.result.stress()
stress_op.inputs.data_sources(model)
stress = stress_op.outputs.fields_container()

# Get displacement results
disp_op = dpf.operators.result.displacement()
disp_op.inputs.data_sources(model)
displacement = disp_op.outputs.fields_container()

# Export to VTK
export_op = dpf.operators.serialization.vtk_export()
export_op.inputs.mesh(mesh)
export_op.inputs.fields1(stress[0])
export_op.inputs.file_path("output.vtk")
export_op.run()
```

**Converting Between Formats:**

```python
import pandas as pd

def export_mapdl_to_csv(mapdl, filename):
    """Export MAPDL results to CSV"""
    mapdl.post1()
    mapdl.set(1, 1)

    # Get all node numbers
    nodes = mapdl.mesh.nnum

    # Get displacements
    disp = mapdl.post_processing.nodal_displacement('ALL')

    # Create DataFrame
    df = pd.DataFrame({
        'Node': nodes,
        'UX': disp[:, 0],
        'UY': disp[:, 1],
        'UZ': disp[:, 2]
    })

    df.to_csv(filename, index=False)
    print(f"Results exported to {filename}")
```

### Best Practices

1. **Session Management:**
   - Always use try-finally blocks to ensure proper cleanup
   - Close sessions explicitly to release licenses
   - Use context managers when available

2. **Automation:**
   - Validate inputs before running analyses
   - Implement checkpoints for long simulations
   - Log all operations for debugging
   - Handle errors gracefully

3. **Performance:**
   - Use appropriate number of processors for problem size
   - Enable result file writing only when needed
   - Use binary format for large data transfers
   - Leverage GPU acceleration when available

4. **Licensing:**
   - Check license availability before launching
   - Use appropriate license features
   - Release licenses promptly after completion
   - Monitor license usage in batch jobs

5. **Results Management:**
   - Use consistent naming conventions
   - Archive results systematically
   - Document simulation parameters
   - Validate results against analytical solutions

---

## SolidWorks

### Overview and Capabilities

SolidWorks is a 3D CAD software for mechanical design:

- Parametric part and assembly modeling
- Sheet metal design
- Weldments and frame structures
- Surface modeling
- Motion simulation
- Flow simulation (CFD add-in)
- Finite element analysis (Simulation add-in)
- Drawing generation

### Setup Instructions

**Installation:**

1. Download SolidWorks installation manager
2. Run installation (requires Windows)
3. Activate license (standalone or network)

**API Access Setup:**

```bash
# Install Python win32com for COM interface
pip install pywin32

# Or use SolidWorks API via .NET (pythonnet)
pip install pythonnet
```

**Enable API Access in SolidWorks:**

1. Tools → Options → System Options → General
2. Check "Enable API automation"

### API/Scripting Interface

SolidWorks provides several API interfaces:

**1. COM API (VBA, Python via win32com):**

```python
import win32com.client
import pythoncom

# Connect to SolidWorks
swApp = win32com.client.Dispatch("SldWorks.Application")
```

**2. C# .NET API**

**3. VBA (built-in macro editor)**

**4. JavaScript API (for web-based applications)**

### Authentication and Licensing

SolidWorks uses FlexLM or SolidNetWork License (SNL) manager:

**Standalone License:**
- Activated via serial number
- Tied to specific machine

**Network License:**
- License server configuration
- Port: 25734 (SolidNetWork License Manager)

**License Configuration:**

```python
# Check if SolidWorks is running with license
def check_solidworks_license():
    try:
        swApp = win32com.client.Dispatch("SldWorks.Application")
        if swApp is not None:
            print(f"SolidWorks version: {swApp.RevisionNumber()}")
            return True
    except:
        print("Cannot connect to SolidWorks")
        return False
```

### Example Automation Scripts

**Complete Part Creation and Export:**

```python
#!/usr/bin/env python
"""
SolidWorks automation using Python COM interface
Requires: pywin32, Windows, SolidWorks installed
"""

import win32com.client
import pythoncom
import os
import time

class SolidWorksAutomation:
    def __init__(self, visible=True):
        """Initialize SolidWorks connection"""
        self.swApp = None
        self.visible = visible

    def connect(self):
        """Connect to SolidWorks application"""
        print("Connecting to SolidWorks...")
        try:
            # Try to get running instance
            self.swApp = win32com.client.GetActiveObject("SldWorks.Application")
            print("Connected to running SolidWorks instance")
        except:
            # Start new instance
            self.swApp = win32com.client.Dispatch("SldWorks.Application")
            print("Started new SolidWorks instance")

        self.swApp.Visible = self.visible

        # Get constants
        self.swConst = win32com.client.constants

        return self.swApp

    def create_part(self):
        """Create new part document"""
        print("Creating new part...")

        # Get template path
        template_path = self.swApp.GetUserPreferenceStringValue(3)  # swDefaultTemplatePart

        # Create new document
        part = self.swApp.NewDocument(template_path, 0, 0, 0)

        if part is None:
            raise RuntimeError("Failed to create part")

        return part

    def create_rectangle_sketch(self, model, plane="Front",
                               x=0, y=0, width=100, height=50):
        """Create rectangular sketch on specified plane"""
        print(f"Creating rectangle sketch on {plane} plane...")

        # Select plane
        plane_name = f"{plane} Plane"
        model.Extension.SelectByID2(plane_name, "PLANE", 0, 0, 0,
                                   False, 0, None, 0)

        # Create sketch
        model.SketchManager.InsertSketch(True)

        # Create rectangle (in mm)
        x1 = x - width/2
        y1 = y - height/2
        x2 = x + width/2
        y2 = y + height/2

        model.SketchManager.CreateCenterRectangle(
            x/1000, y/1000, 0,  # Center point (convert to meters)
            x2/1000, y2/1000, 0  # Corner point
        )

        return model

    def extrude_feature(self, model, depth=10, reverse=False):
        """Extrude the current sketch"""
        print(f"Extruding sketch by {depth}mm...")

        # Get feature manager
        feat_mgr = model.FeatureManager

        # Extrude (depth in meters)
        feature = feat_mgr.FeatureExtrusion2(
            True,  # Single direction
            False,  # No draft
            reverse,  # Reverse direction
            0,  # Draft type
            0,  # Draft angle
            depth/1000,  # Depth (convert to meters)
            0,  # Second depth
            False,  # Merge result
            False,  # Thin feature
            False,  # Flip side
            False,  # Direction type
            0,  # Direction type
            False  # Both directions
        )

        if feature is None:
            raise RuntimeError("Extrusion failed")

        return feature

    def create_circular_pattern(self, model, feature_name,
                               count=4, angle=360):
        """Create circular pattern of feature"""
        print(f"Creating circular pattern: {count} instances...")

        # Select feature to pattern
        model.Extension.SelectByID2(feature_name, "BODYFEATURE",
                                   0, 0, 0, False, 0, None, 0)

        # Select axis
        model.Extension.SelectByID2("", "AXIS", 0, 0, 0,
                                   True, 0, None, 0)

        # Create pattern
        feat_mgr = model.FeatureManager
        pattern = feat_mgr.FeatureCircularPattern4(
            count,  # Number of instances
            angle * 3.14159/180,  # Angle in radians
            False,  # Equal spacing
            False,  # Geometry pattern
            0,  # Skip instances
            False,  # Staggered
            False,  # Vary sketch
            False  # Vary instances
        )

        return pattern

    def add_fillet(self, model, edges, radius=5):
        """Add fillet to edges"""
        print(f"Adding fillet with radius {radius}mm...")

        # Select edges
        for i, edge in enumerate(edges):
            model.Extension.SelectByID2(edge, "EDGE",
                                       0, 0, 0,
                                       i > 0, 0, None, 0)

        # Create fillet
        feat_mgr = model.FeatureManager
        fillet = feat_mgr.FeatureFillet2(
            0,  # Fillet type (constant radius)
            radius/1000,  # Radius (convert to meters)
            0,  # Second radius
            0,  # Setback radius
            False,  # Overflow type
            False,  # Enable setback
            False,  # Tangent propagation
            False,  # Full preview
            False,  # Overflow type
            False  # Use feature scope
        )

        return fillet

    def apply_material(self, model, material_name="Plain Carbon Steel"):
        """Apply material to part"""
        print(f"Applying material: {material_name}")

        # Get material database path
        material_db = os.path.join(
            os.path.dirname(self.swApp.GetExecutablePath()),
            "lang", "english", "sldmaterials", "solidworks materials.sldmat"
        )

        # Set material
        model.SetMaterialPropertyName2("", "", material_db, material_name)

    def get_mass_properties(self, model):
        """Get mass properties of part"""
        print("Calculating mass properties...")

        # Get mass properties
        mass_props = model.Extension.CreateMassProperty()

        if mass_props is None:
            print("Failed to get mass properties")
            return None

        props = {
            'mass': mass_props.Mass,
            'volume': mass_props.Volume,
            'surface_area': mass_props.SurfaceArea,
            'center_of_mass': mass_props.CenterOfMass,
            'moments_of_inertia': mass_props.GetMomentOfInertia(0)
        }

        print(f"Mass: {props['mass']:.4f} kg")
        print(f"Volume: {props['volume']*1e9:.2f} mm³")
        print(f"Surface Area: {props['surface_area']*1e6:.2f} mm²")

        return props

    def save_part(self, filepath, save_as_version=None):
        """Save part file"""
        print(f"Saving part to {filepath}")

        # Get active document
        model = self.swApp.ActiveDoc

        # Save options
        save_options = 0  # swSaveAsOptions_Silent
        warnings = 0
        errors = 0

        # Save
        if save_as_version:
            success = model.SaveAs4(filepath, save_as_version,
                                   save_options, errors, warnings)
        else:
            success = model.SaveAs3(filepath, save_options, 0)

        if not success:
            raise RuntimeError(f"Failed to save part (errors: {errors})")

        print(f"Part saved successfully")

    def export_step(self, filepath):
        """Export part as STEP file"""
        print(f"Exporting to STEP: {filepath}")

        model = self.swApp.ActiveDoc

        # Set STEP export options
        # 1 = AP203, 2 = AP214, 3 = AP242
        self.swApp.SetUserPreferenceIntegerValue(
            108, 2)  # swStepAP = 108

        # Export
        success = model.SaveAs3(filepath, 0, 2)

        if not success:
            raise RuntimeError("STEP export failed")

        print("STEP export successful")

    def export_stl(self, filepath, quality='fine'):
        """Export part as STL file"""
        print(f"Exporting to STL: {filepath}")

        model = self.swApp.ActiveDoc

        # Set STL quality
        if quality == 'fine':
            # Fine quality settings
            self.swApp.SetUserPreferenceDoubleValue(104, 0.0001)  # Deviation
            self.swApp.SetUserPreferenceDoubleValue(105, 5.0)  # Angle

        # Export
        success = model.SaveAs3(filepath, 0, 2)

        if not success:
            raise RuntimeError("STL export failed")

        print("STL export successful")

    def create_drawing(self, part_path, template="A3"):
        """Create drawing from part"""
        print("Creating drawing...")

        # Get drawing template
        template_path = self.swApp.GetUserPreferenceStringValue(4)  # swDefaultTemplateDrawing

        # Create drawing
        drawing = self.swApp.NewDocument(template_path, 0, 0, 0)

        if drawing is None:
            raise RuntimeError("Failed to create drawing")

        # Get active sheet
        sheet = drawing.GetCurrentSheet()

        # Insert view
        model_view = drawing.CreateDrawViewFromModelView3(
            part_path,  # Model path
            "*Front",  # View orientation
            0.1, 0.2, 0  # X, Y, Z position (meters)
        )

        return drawing

    def close(self):
        """Close SolidWorks connection"""
        print("Closing SolidWorks connection...")
        # Don't close the application if it was already running
        # self.swApp.ExitApp()


# Example usage
if __name__ == "__main__":
    # Initialize automation
    sw = SolidWorksAutomation(visible=True)
    sw.connect()

    try:
        # Create part
        part = sw.create_part()
        model = sw.swApp.ActiveDoc

        # Create base sketch and extrude
        sw.create_rectangle_sketch(model, plane="Front",
                                   width=100, height=50)
        sw.extrude_feature(model, depth=20)

        # Apply material
        sw.apply_material(model, "Plain Carbon Steel")

        # Get mass properties
        props = sw.get_mass_properties(model)

        # Save part
        output_dir = r"C:\Temp\SolidWorks"
        os.makedirs(output_dir, exist_ok=True)

        part_path = os.path.join(output_dir, "test_part.SLDPRT")
        sw.save_part(part_path)

        # Export formats
        sw.export_step(os.path.join(output_dir, "test_part.STEP"))
        sw.export_stl(os.path.join(output_dir, "test_part.STL"))

        print("\nAutomation completed successfully!")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
```

**Batch Processing Script:**

```python
#!/usr/bin/env python
"""
Batch process SolidWorks files
"""

import win32com.client
import os
import glob

class SolidWorksBatchProcessor:
    def __init__(self):
        self.swApp = None

    def connect(self):
        """Connect to SolidWorks"""
        try:
            self.swApp = win32com.client.GetActiveObject("SldWorks.Application")
        except:
            self.swApp = win32com.client.Dispatch("SldWorks.Application")

        self.swApp.Visible = False  # Run in background

    def process_directory(self, input_dir, output_dir, file_types=['SLDPRT']):
        """Process all files in directory"""
        print(f"Processing files in: {input_dir}")

        os.makedirs(output_dir, exist_ok=True)

        results = []

        for file_type in file_types:
            pattern = os.path.join(input_dir, f"*.{file_type}")
            files = glob.glob(pattern)

            for filepath in files:
                print(f"\nProcessing: {os.path.basename(filepath)}")

                try:
                    # Open file
                    errors = 0
                    warnings = 0
                    doc = self.swApp.OpenDoc6(filepath, 1, 0, "",
                                             errors, warnings)

                    if doc is None:
                        print(f"  Failed to open (errors: {errors})")
                        continue

                    # Get mass properties
                    mass_props = doc.Extension.CreateMassProperty()

                    result = {
                        'filename': os.path.basename(filepath),
                        'mass': mass_props.Mass if mass_props else None,
                        'volume': mass_props.Volume if mass_props else None
                    }

                    # Export STEP
                    step_path = os.path.join(output_dir,
                                            os.path.splitext(os.path.basename(filepath))[0] + ".STEP")
                    doc.SaveAs3(step_path, 0, 2)
                    print(f"  Exported: {os.path.basename(step_path)}")

                    # Close document
                    self.swApp.CloseDoc(filepath)

                    results.append(result)

                except Exception as e:
                    print(f"  Error: {e}")
                    results.append({
                        'filename': os.path.basename(filepath),
                        'error': str(e)
                    })

        return results

    def generate_report(self, results, output_file):
        """Generate CSV report"""
        import csv

        with open(output_file, 'w', newline='') as f:
            if results:
                writer = csv.DictWriter(f, fieldnames=results[0].keys())
                writer.writeheader()
                writer.writerows(results)

        print(f"\nReport generated: {output_file}")

# Example usage
if __name__ == "__main__":
    processor = SolidWorksBatchProcessor()
    processor.connect()

    results = processor.process_directory(
        input_dir=r"C:\Projects\Parts",
        output_dir=r"C:\Projects\Export",
        file_types=['SLDPRT', 'SLDASM']
    )

    processor.generate_report(results, "processing_report.csv")
```

### File I/O and Data Exchange

**Reading/Writing Part Parameters:**

```python
def get_custom_properties(model):
    """Get custom properties from SolidWorks model"""
    custom_prop_mgr = model.Extension.CustomPropertyManager("")

    # Get property names
    names = custom_prop_mgr.GetNames()

    properties = {}
    if names:
        for name in names:
            val = ""
            resolved_val = ""
            custom_prop_mgr.Get5(name, False, val, resolved_val, False)
            properties[name] = resolved_val

    return properties

def set_custom_property(model, name, value, prop_type="Text"):
    """Set custom property"""
    custom_prop_mgr = model.Extension.CustomPropertyManager("")

    # Property types: Text, Date, Number, Yes or No
    custom_prop_mgr.Add3(name, 30, value, 2)  # 30=swCustomInfoText
```

**Importing External Data:**

```python
def import_coordinates_from_csv(model, csv_file):
    """Import 3D points from CSV and create sketch points"""
    import csv

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        points = [(float(row['X']), float(row['Y']), float(row['Z']))
                 for row in reader]

    # Create 3D sketch
    model.SketchManager.Insert3DSketch(True)

    # Add points
    for x, y, z in points:
        model.SketchManager.CreatePoint(x/1000, y/1000, z/1000)  # Convert to meters

    model.SketchManager.Insert3DSketch(True)  # Exit sketch
```

### Best Practices

1. **COM Interface:**
   - Always handle COM exceptions
   - Release COM objects properly
   - Use early binding when possible
   - Check return values for None

2. **Performance:**
   - Disable screen updates for batch operations
   - Use invisible mode for background processing
   - Close documents after processing
   - Clear selection sets

3. **Error Handling:**
   - Check for errors in OpenDoc and SaveAs operations
   - Validate inputs before operations
   - Log all operations for debugging
   - Implement retry logic for file operations

4. **Units:**
   - SolidWorks API uses meters internally
   - Convert from mm to meters for API calls
   - Be consistent with units in calculations
   - Document unit conventions

5. **Version Compatibility:**
   - Test with different SolidWorks versions
   - Use version-specific API features carefully
   - Implement fallbacks for older versions
   - Check API compatibility before deployment

---

## COMSOL Multiphysics

### Overview and Capabilities

COMSOL Multiphysics is a finite element analysis (FEA) software for modeling coupled physics:

- Structural mechanics
- Fluid dynamics (CFD)
- Heat transfer
- Electromagnetics
- Acoustics
- Chemical engineering
- Multiphysics coupling
- Optimization and parametric sweeps

### Setup Instructions

**Installation:**

1. Download COMSOL installer
2. Run installer (requires license file)
3. Configure COMSOL Server for remote access (optional)

**COMSOL with MATLAB:**

```matlab
% Add COMSOL MATLAB API to path
addpath('C:\Program Files\COMSOL\COMSOL61\Multiphysics\mli')

% Initialize
mphstart
```

**COMSOL with Java:**

```bash
# Set classpath
export CLASSPATH=$CLASSPATH:/usr/local/comsol61/multiphysics/plugins/*
```

**COMSOL with Python:**

```bash
# COMSOL Python client (via COMSOL Server)
pip install mph
```

### API/Scripting Interface

COMSOL provides several interfaces:

**1. MATLAB API (Native)**

**2. Java API**

**3. Python (via mph package)**

**4. COMSOL Script (JavaScript-like)**

### Authentication and Licensing

COMSOL uses FlexLM license management:

**License Configuration:**

```bash
# Linux
export LMCOMSOL_LICENSE_FILE=1718@license-server.company.com

# Windows
set LMCOMSOL_LICENSE_FILE=1718@license-server.company.com
```

**Check License:**

```matlab
% In MATLAB
mphstart('comsolserver.company.com', 2036)
```

### Example Automation Scripts

**MATLAB: Complete Heat Transfer Simulation:**

```matlab
%% COMSOL Automation Example - Heat Transfer in a Rod
% Demonstrates automated model creation, solving, and post-processing

function comsol_heat_transfer_automation()
    % Initialize COMSOL
    import com.comsol.model.*
    import com.comsol.model.util.*

    % Start COMSOL
    ModelUtil.showProgress(true);

    % Create model
    model = ModelUtil.create('Model');
    model.modelPath('/path/to/models');
    model.label('Heat Transfer in Rod');

    %% Parameters
    % Define parameters for parametric study
    model.param.set('L', '0.5[m]', 'Rod length');
    model.param.set('r', '0.01[m]', 'Rod radius');
    model.param.set('k_thermal', '50[W/(m*K)]', 'Thermal conductivity');
    model.param.set('T_left', '373[K]', 'Left boundary temperature');
    model.param.set('T_right', '293[K]', 'Right boundary temperature');
    model.param.set('h_conv', '10[W/(m^2*K)]', 'Convection coefficient');
    model.param.set('T_ambient', '293[K]', 'Ambient temperature');

    %% Geometry
    % Create geometry component
    model.component.create('comp1', true);
    model.component('comp1').geom.create('geom1', 3);

    % Create cylinder (rod)
    model.component('comp1').geom('geom1').create('cyl1', 'Cylinder');
    model.component('comp1').geom('geom1').feature('cyl1').set('r', 'r');
    model.component('comp1').geom('geom1').feature('cyl1').set('h', 'L');
    model.component('comp1').geom('geom1').feature('cyl1').set('axis', [1 0 0]);

    % Build geometry
    model.component('comp1').geom('geom1').run;

    %% Material
    % Create material
    model.component('comp1').material.create('mat1', 'Common');
    model.component('comp1').material('mat1').label('Steel');

    % Set material properties
    model.component('comp1').material('mat1').propertyGroup('def').set('density', '7850[kg/m^3]');
    model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '475[J/(kg*K)]');
    model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k_thermal' '0' '0' '0' 'k_thermal' '0' '0' '0' 'k_thermal'});

    %% Physics
    % Add heat transfer physics
    model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');

    % Temperature boundary condition (left end)
    model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 2);
    model.component('comp1').physics('ht').feature('temp1').selection.set([1]);
    model.component('comp1').physics('ht').feature('temp1').set('T0', 'T_left');

    % Temperature boundary condition (right end)
    model.component('comp1').physics('ht').create('temp2', 'TemperatureBoundary', 2);
    model.component('comp1').physics('ht').feature('temp2').selection.set([2]);
    model.component('comp1').physics('ht').feature('temp2').set('T0', 'T_right');

    % Convective heat transfer (cylindrical surface)
    model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 2);
    model.component('comp1').physics('ht').feature('hf1').selection.set([3]);
    model.component('comp1').physics('ht').feature('hf1').set('HeatFluxType', 'ConvectiveHeatFlux');
    model.component('comp1').physics('ht').feature('hf1').set('h', 'h_conv');
    model.component('comp1').physics('ht').feature('hf1').set('Text', 'T_ambient');

    %% Mesh
    % Create mesh
    model.component('comp1').mesh.create('mesh1');

    % Set mesh size
    model.component('comp1').mesh('mesh1').autoMeshSize(4); % Fine mesh

    % Build mesh
    model.component('comp1').mesh('mesh1').run;

    %% Study
    % Create stationary study
    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');

    % Run study
    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('s1', 'Stationary');

    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').feature.remove('fcDef');

    % Solve
    model.sol('sol1').runAll;

    %% Post-processing
    % Create temperature plot
    model.result.create('pg1', 'PlotGroup3D');
    model.result('pg1').label('Temperature Distribution');
    model.result('pg1').create('surf1', 'Surface');
    model.result('pg1').feature('surf1').set('expr', 'T');

    % Create temperature along centerline
    model.result.create('pg2', 'PlotGroup1D');
    model.result('pg2').label('Temperature Profile');

    % Extract data along centerline
    model.result.dataset.create('cln1', 'CutLine3D');
    model.result.dataset('cln1').set('data', 'dset1');
    model.result.dataset('cln1').set('pointdata', {{0,0,0}, {L,0,0}});

    model.result('pg2').setIndex('looplevel', [1], 0);
    model.result('pg2').create('lngr1', 'LineGraph');
    model.result('pg2').feature('lngr1').set('xdata', 'expr');
    model.result('pg2').feature('lngr1').set('xdataexpr', 'x');
    model.result('pg2').feature('lngr1').set('ydata', 'expr');
    model.result('pg2').feature('lngr1').set('expr', 'T');
    model.result('pg2').feature('lngr1').set('data', 'cln1');

    %% Export Results
    % Export temperature data
    model.result.export.create('data1', 'Data');
    model.result.export('data1').set('expr', {'T'});
    model.result.export('data1').set('unit', {'K'});
    model.result.export('data1').set('descr', {'Temperature'});
    model.result.export('data1').set('filename', 'temperature_data.txt');
    model.result.export('data1').run;

    % Export image
    model.result('pg1').run;
    model.result.export.create('img1', 'Image');
    model.result.export('img1').set('sourceobject', 'pg1');
    model.result.export('img1').set('filename', 'temperature_plot.png');
    model.result.export('img1').run;

    %% Save Model
    % Save mph file
    mphsave(model, 'heat_transfer_rod.mph');

    fprintf('Simulation completed successfully!\n');
    fprintf('Model saved as: heat_transfer_rod.mph\n');

    %% Parametric Study
    % Run parametric sweep for different conductivities
    conductivities = [10, 50, 100, 200];
    max_temps = zeros(size(conductivities));

    for i = 1:length(conductivities)
        fprintf('Running case %d: k = %d W/(m*K)\n', i, conductivities(i));

        % Set parameter
        model.param.set('k_thermal', sprintf('%d[W/(m*K)]', conductivities(i)));

        % Solve
        model.sol('sol1').runAll;

        % Extract maximum temperature
        max_temps(i) = mphmax(model, 'T', 'surface');

        fprintf('  Maximum temperature: %.2f K\n', max_temps(i));
    end

    % Plot parametric results
    figure;
    plot(conductivities, max_temps, '-o', 'LineWidth', 2);
    xlabel('Thermal Conductivity [W/(m*K)]');
    ylabel('Maximum Temperature [K]');
    title('Maximum Temperature vs Thermal Conductivity');
    grid on;
    saveas(gcf, 'parametric_study.png');

end
```

**Python: COMSOL Client Automation:**

```python
#!/usr/bin/env python3
"""
COMSOL automation using Python mph client
Requires: mph package and running COMSOL Server
"""

import mph
import numpy as np
import matplotlib.pyplot as plt

class COMSOLAutomation:
    def __init__(self, server='localhost', port=2036):
        """Initialize COMSOL client"""
        self.client = None
        self.server = server
        self.port = port

    def connect(self):
        """Connect to COMSOL Server"""
        print(f"Connecting to COMSOL Server at {self.server}:{self.port}")
        self.client = mph.Client(cores=4)
        print("Connected successfully")
        return self.client

    def load_model(self, model_file):
        """Load existing COMSOL model"""
        print(f"Loading model: {model_file}")
        model = self.client.load(model_file)
        return model

    def run_parametric_sweep(self, model, param_name, param_values):
        """Run parametric sweep"""
        print(f"Running parametric sweep for {param_name}")

        results = []

        for value in param_values:
            print(f"  {param_name} = {value}")

            # Set parameter
            model.parameter(param_name, value)

            # Solve
            model.solve()

            # Extract results
            temp_data = model.evaluate('T')
            max_temp = np.max(temp_data)

            results.append({
                'parameter': value,
                'max_temperature': max_temp
            })

            print(f"    Max temperature: {max_temp:.2f} K")

        return results

    def extract_data(self, model, expression, dataset=None):
        """Extract data from model"""
        if dataset:
            data = model.evaluate(expression, dataset=dataset)
        else:
            data = model.evaluate(expression)

        return data

    def export_data(self, model, filename, expressions):
        """Export data to file"""
        print(f"Exporting data to {filename}")

        # Export using COMSOL's export functionality
        model.export(filename, expressions)

    def close(self):
        """Disconnect from server"""
        if self.client:
            print("Disconnecting from COMSOL Server")
            self.client.disconnect()


# Example usage
if __name__ == "__main__":
    # Initialize
    comsol = COMSOLAutomation()
    comsol.connect()

    try:
        # Load model
        model = comsol.load_model('heat_transfer_rod.mph')

        # Run parametric sweep
        conductivities = [10, 50, 100, 200]
        results = comsol.run_parametric_sweep(
            model,
            'k_thermal',
            conductivities
        )

        # Plot results
        params = [r['parameter'] for r in results]
        temps = [r['max_temperature'] for r in results]

        plt.figure(figsize=(10, 6))
        plt.plot(params, temps, '-o', linewidth=2)
        plt.xlabel('Thermal Conductivity [W/(m*K)]')
        plt.ylabel('Maximum Temperature [K]')
        plt.title('Parametric Study Results')
        plt.grid(True)
        plt.savefig('comsol_results.png')

        print("Analysis completed successfully!")

    finally:
        comsol.close()
```

### File I/O and Data Exchange

**Reading COMSOL Results:**

```matlab
% Load model
model = mphload('model.mph');

% Extract solution data
T = mpheval(model, 'T', 'dataset', 'dset1');

% Get coordinates and values
x = T.p(1,:);  % X coordinates
y = T.p(2,:);  % Y coordinates
z = T.p(3,:);  % Z coordinates
temp = T.d1;   % Temperature values

% Export to MAT file
save('results.mat', 'x', 'y', 'z', 'temp');
```

**Importing External Data:**

```matlab
% Import CSV data
data = readtable('input_data.csv');

% Create interpolation function
model.func.create('int1', 'Interpolation');
model.func('int1').set('source', 'file');
model.func('int1').set('filename', 'input_data.csv');

% Use in physics
model.component('comp1').physics('ht').feature('hs1').set('Q0', 'int1(x,y)');
```

### Best Practices

1. **Model Organization:**
   - Use parameters for all numerical values
   - Group related features
   - Add descriptive labels and comments
   - Use version control for .mph files

2. **Performance:**
   - Start with coarse mesh, refine iteratively
   - Use symmetry to reduce problem size
   - Enable parallel solving when possible
   - Store only necessary solution data

3. **Automation:**
   - Validate inputs before solving
   - Implement convergence checks
   - Save intermediate results
   - Log all parameter variations

4. **Results Validation:**
   - Compare with analytical solutions
   - Check energy balance
   - Verify mesh independence
   - Document assumptions

5. **Server Mode:**
   - Use COMSOL Server for remote access
   - Implement proper error handling
   - Monitor license usage
   - Clean up server sessions

---

## MATLAB

### Overview and Capabilities

MATLAB provides extensive engineering toolboxes:

- **Control System Toolbox** - Control design and analysis
- **Signal Processing Toolbox** - Filter design, spectral analysis
- **Optimization Toolbox** - Linear/nonlinear optimization
- **Partial Differential Equation Toolbox** - FEA and PDE solving
- **Simulink** - Dynamic system simulation
- **CFD Toolbox** - Computational fluid dynamics
- **Statistics and Machine Learning Toolbox**

### Setup Instructions

**Installation:**

1. Download MATLAB installer from MathWorks
2. Run installer with license
3. Select desired toolboxes

**Python Integration:**

```bash
# Install MATLAB Engine for Python
cd /usr/local/MATLAB/R2024a/extern/engines/python
python setup.py install
```

**Environment Setup:**

```bash
# Add MATLAB to PATH (Linux)
export PATH=$PATH:/usr/local/MATLAB/R2024a/bin

# Set library path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/MATLAB/R2024a/bin/glnxa64
```

### API/Scripting Interface

MATLAB provides multiple interfaces:

**1. Native M-files**

**2. Python API:**

```python
import matlab.engine

# Start MATLAB engine
eng = matlab.engine.start_matlab()

# Call MATLAB functions
result = eng.sqrt(4.0)
```

**3. Java API**

**4. C/C++ API**

**5. REST API (MATLAB Production Server)**

### Authentication and Licensing

MATLAB uses FlexLM network license manager:

**License Configuration:**

```bash
# Set license file location
export MLM_LICENSE_FILE=27000@license-server.company.com

# Or individual license file
export MLM_LICENSE_FILE=/path/to/license.lic
```

**Check License:**

```matlab
% Check available toolboxes
ver

% Check license
license('inuse')

% Test specific toolbox
license('test', 'optimization_toolbox')
```

### Example Automation Scripts

**Complete Optimization Example:**

```matlab
%% Engineering Optimization Example
% Optimize beam dimensions for minimum weight with strength constraints

function beam_optimization_example()
    %% Problem Definition
    % Design variables: [width, height, length] (all in meters)
    % Objective: Minimize weight
    % Constraints: Maximum stress, deflection limits

    % Material properties (Steel)
    E = 200e9;  % Young's modulus (Pa)
    sigma_max = 250e6;  % Maximum allowable stress (Pa)
    rho = 7850;  % Density (kg/m^3)

    % Loading
    F = 10000;  % Applied force (N)
    L = 2.0;  % Beam length (m)

    %% Optimization Setup
    % Design variable bounds: [width, height]
    lb = [0.01, 0.01];  % Lower bounds (m)
    ub = [0.5, 0.5];    % Upper bounds (m)

    % Initial guess
    x0 = [0.1, 0.1];

    % Optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point', ...
        'MaxIterations', 100, ...
        'OptimalityTolerance', 1e-6, ...
        'StepTolerance', 1e-10);

    %% Define Objective Function
    function weight = objective(x)
        width = x(1);
        height = x(2);

        % Calculate volume and weight
        volume = width * height * L;
        weight = rho * volume;
    end

    %% Define Constraints
    function [c, ceq] = constraints(x)
        width = x(1);
        height = x(2);

        % Moment of inertia
        I = (width * height^3) / 12;

        % Maximum bending moment (cantilever with end load)
        M_max = F * L;

        % Maximum bending stress
        sigma = M_max * (height/2) / I;

        % Maximum deflection
        delta = (F * L^3) / (3 * E * I);
        delta_max = L / 200;  % Deflection limit: L/200

        % Inequality constraints (must be <= 0)
        c = [
            sigma - sigma_max;  % Stress constraint
            delta - delta_max   % Deflection constraint
        ];

        % Equality constraints (none)
        ceq = [];
    end

    %% Run Optimization
    fprintf('Starting optimization...\n\n');

    [x_opt, weight_opt, exitflag, output] = fmincon(...
        @objective, ...
        x0, ...
        [], [], ...  % No linear inequality constraints
        [], [], ...  % No linear equality constraints
        lb, ub, ...
        @constraints, ...
        options);

    %% Display Results
    fprintf('\n=== Optimization Results ===\n');
    fprintf('Exit flag: %d\n', exitflag);
    fprintf('Iterations: %d\n', output.iterations);
    fprintf('Function evaluations: %d\n', output.funcCount);
    fprintf('\nOptimal dimensions:\n');
    fprintf('  Width: %.4f m\n', x_opt(1));
    fprintf('  Height: %.4f m\n', x_opt(2));
    fprintf('  Length: %.4f m\n', L);
    fprintf('\nMinimum weight: %.2f kg\n', weight_opt);

    % Verify constraints at optimal point
    [c_opt, ~] = constraints(x_opt);
    fprintf('\nConstraint values at optimum:\n');
    fprintf('  Stress margin: %.2e Pa (should be <= 0)\n', c_opt(1));
    fprintf('  Deflection margin: %.2e m (should be <= 0)\n', c_opt(2));

    %% Visualization
    % Create design space visualization
    [W, H] = meshgrid(linspace(lb(1), ub(1), 50), ...
                      linspace(lb(2), ub(2), 50));

    % Calculate objective for each point
    Weight = zeros(size(W));
    StressViolation = zeros(size(W));

    for i = 1:numel(W)
        Weight(i) = objective([W(i), H(i)]);
        [c, ~] = constraints([W(i), H(i)]);
        StressViolation(i) = max(c(1), 0);  % Positive if violated
    end

    % Plot
    figure('Position', [100, 100, 1200, 400]);

    subplot(1, 3, 1);
    contourf(W, H, Weight, 20);
    hold on;
    plot(x_opt(1), x_opt(2), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    colorbar;
    xlabel('Width (m)');
    ylabel('Height (m)');
    title('Weight (kg)');

    subplot(1, 3, 2);
    contourf(W, H, StressViolation, 20);
    hold on;
    plot(x_opt(1), x_opt(2), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    colorbar;
    xlabel('Width (m)');
    ylabel('Height (m)');
    title('Stress Constraint Violation (Pa)');

    subplot(1, 3, 3);
    % Plot convergence history
    if isfield(output, 'iterations')
        plot(1:output.iterations, '-o', 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Objective Value');
        title('Convergence History');
        grid on;
    end

    saveas(gcf, 'optimization_results.png');

    fprintf('\nOptimization complete! Results saved to optimization_results.png\n');
end
```

**Python Integration Example:**

```python
#!/usr/bin/env python3
"""
MATLAB Engine for Python - Engineering automation
"""

import matlab.engine
import numpy as np
import matplotlib.pyplot as plt

class MATLABAutomation:
    def __init__(self):
        """Initialize MATLAB engine"""
        self.eng = None

    def start(self):
        """Start MATLAB engine"""
        print("Starting MATLAB engine...")
        self.eng = matlab.engine.start_matlab()
        print("MATLAB engine started")
        return self.eng

    def run_script(self, script_name):
        """Run MATLAB script"""
        print(f"Running MATLAB script: {script_name}")
        self.eng.eval(f"run('{script_name}')", nargout=0)

    def call_function(self, func_name, *args):
        """Call MATLAB function"""
        func = getattr(self.eng, func_name)
        result = func(*args)
        return result

    def transfer_array(self, np_array):
        """Transfer NumPy array to MATLAB"""
        # Convert to MATLAB array
        matlab_array = matlab.double(np_array.tolist())
        return matlab_array

    def transfer_from_matlab(self, matlab_array):
        """Transfer MATLAB array to NumPy"""
        # Convert to NumPy array
        np_array = np.array(matlab_array)
        return np_array

    def solve_ode(self, func, tspan, y0):
        """Solve ODE using MATLAB's ode45"""
        # Define ODE function in MATLAB workspace
        self.eng.workspace['odefun'] = func

        # Convert inputs to MATLAB format
        tspan_matlab = matlab.double(tspan)
        y0_matlab = matlab.double(y0)

        # Solve ODE
        t, y = self.eng.ode45(func, tspan_matlab, y0_matlab, nargout=2)

        # Convert back to NumPy
        t_np = self.transfer_from_matlab(t)
        y_np = self.transfer_from_matlab(y)

        return t_np, y_np

    def fft_analysis(self, signal, fs):
        """Perform FFT analysis"""
        # Transfer signal to MATLAB
        signal_matlab = self.transfer_array(signal)

        # Perform FFT
        fft_result = self.eng.fft(signal_matlab)

        # Calculate frequency vector
        n = len(signal)
        f = self.eng.linspace(0.0, fs/2, n//2)

        # Convert results
        fft_np = self.transfer_from_matlab(fft_result)
        f_np = self.transfer_from_matlab(f)

        # Calculate magnitude
        magnitude = np.abs(fft_np[:n//2])

        return f_np, magnitude

    def design_filter(self, filter_type, order, cutoff, fs):
        """Design digital filter"""
        # Design filter using MATLAB
        if filter_type == 'lowpass':
            b, a = self.eng.butter(order, cutoff/(fs/2), nargout=2)
        elif filter_type == 'highpass':
            b, a = self.eng.butter(order, cutoff/(fs/2), 'high', nargout=2)
        elif filter_type == 'bandpass':
            b, a = self.eng.butter(order, [cutoff[0]/(fs/2), cutoff[1]/(fs/2)],
                                  'bandpass', nargout=2)

        # Convert to NumPy
        b_np = self.transfer_from_matlab(b)
        a_np = self.transfer_from_matlab(a)

        return b_np, a_np

    def optimize(self, objective, x0, bounds):
        """Run optimization using fmincon"""
        # Convert inputs
        x0_matlab = matlab.double(x0)
        lb_matlab = matlab.double([b[0] for b in bounds])
        ub_matlab = matlab.double([b[1] for b in bounds])

        # Define objective in MATLAB workspace
        self.eng.workspace['objective_func'] = objective

        # Run optimization
        x_opt = self.eng.fmincon(
            objective,
            x0_matlab,
            [], [],  # No linear constraints
            [], [],  # No linear equality constraints
            lb_matlab, ub_matlab
        )

        # Convert result
        x_opt_np = self.transfer_from_matlab(x_opt)

        return x_opt_np

    def stop(self):
        """Stop MATLAB engine"""
        if self.eng:
            print("Stopping MATLAB engine...")
            self.eng.quit()
            print("MATLAB engine stopped")


# Example usage
if __name__ == "__main__":
    # Initialize
    matlab = MATLABAutomation()
    matlab.start()

    try:
        # Example 1: FFT Analysis
        print("\n=== FFT Analysis ===")
        fs = 1000  # Sampling frequency
        t = np.linspace(0, 1, fs)
        signal = np.sin(2 * np.pi * 50 * t) + np.sin(2 * np.pi * 120 * t)

        f, magnitude = matlab.fft_analysis(signal, fs)

        plt.figure(figsize=(12, 5))
        plt.subplot(1, 2, 1)
        plt.plot(t[:100], signal[:100])
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.title('Time Domain Signal')
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.plot(f, magnitude)
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Magnitude')
        plt.title('Frequency Domain')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('fft_analysis.png')
        print("FFT analysis complete")

        # Example 2: Filter Design
        print("\n=== Filter Design ===")
        b, a = matlab.design_filter('lowpass', order=4, cutoff=100, fs=1000)
        print(f"Filter coefficients:")
        print(f"  b: {b}")
        print(f"  a: {a}")

        # Example 3: Solve ODE
        print("\n=== ODE Solution ===")
        # Define ODE: dy/dt = -2*y
        matlab.eng.eval("odefun = @(t,y) -2*y;", nargout=0)
        t, y = matlab.solve_ode('odefun', [0, 5], [1.0])

        plt.figure(figsize=(8, 5))
        plt.plot(t, y)
        plt.xlabel('Time')
        plt.ylabel('y(t)')
        plt.title('ODE Solution: dy/dt = -2y, y(0) = 1')
        plt.grid(True)
        plt.savefig('ode_solution.png')
        print("ODE solution complete")

        # Example 4: Call built-in MATLAB functions
        print("\n=== MATLAB Built-in Functions ===")
        result = matlab.call_function('sqrt', 16.0)
        print(f"sqrt(16) = {result}")

        A = np.array([[1, 2], [3, 4]])
        det_A = matlab.call_function('det', matlab.transfer_array(A))
        print(f"det(A) = {det_A}")

    finally:
        matlab.stop()
```

### File I/O and Data Exchange

**MAT-Files:**

```matlab
% Save data
x = 1:100;
y = sin(x);
save('data.mat', 'x', 'y');

% Load data
data = load('data.mat');
x = data.x;
y = data.y;
```

**Python Access:**

```python
import scipy.io as sio

# Read MAT file
mat_data = sio.loadmat('data.mat')
x = mat_data['x']
y = mat_data['y']

# Write MAT file
sio.savemat('output.mat', {'x': x, 'y': y})
```

**Excel Integration:**

```matlab
% Read Excel
data = readtable('input.xlsx');

% Write Excel
writetable(data, 'output.xlsx');

% Read specific range
data = readmatrix('input.xlsx', 'Sheet', 'Sheet1', 'Range', 'A1:D100');
```

### Best Practices

1. **Performance:**
   - Vectorize operations (avoid loops)
   - Preallocate arrays
   - Use appropriate data types
   - Profile code to identify bottlenecks

2. **Code Organization:**
   - Use functions for reusable code
   - Add help text to functions
   - Use meaningful variable names
   - Comment complex algorithms

3. **Python Integration:**
   - Start engine once, reuse for multiple calls
   - Transfer large arrays efficiently
   - Handle MATLAB errors in Python
   - Clean up engine on exit

4. **Toolbox Usage:**
   - Check license availability before use
   - Use toolbox-specific functions when available
   - Consult documentation for best practices
   - Keep toolboxes updated

5. **Debugging:**
   - Use debugging breakpoints
   - Check array dimensions
   - Validate inputs
   - Use try-catch for error handling

---

## Star-CCM+

### Overview and Capabilities

Star-CCM+ is a commercial CFD software by Siemens:

- Advanced CFD for complex flows
- Multiphase flow (VOF, Eulerian, DEM)
- Heat transfer and conjugate heat transfer
- Turbomachinery analysis
- Overset mesh capabilities
- Design exploration and optimization
- Fluid-structure interaction
- Combustion modeling

### Setup Instructions

**Installation:**

1. Download Star-CCM+ installer from Siemens
2. Run installer:
   ```bash
   # Linux
   ./STAR-CCM+17.06.007_01_linux-x86_64.sh

   # Windows
   STAR-CCM+_installer.exe
   ```
3. Configure license server

**Environment Setup:**

```bash
# Linux
export CDLMD_LICENSE_FILE=1999@license-server.company.com
export PATH=$PATH:/opt/Siemens/17.06.007/STAR-CCM+17.06.007/star/bin

# Set working directory
export STAR_WORKING_DIR=/path/to/simulations
```

### API/Scripting Interface

Star-CCM+ provides multiple automation interfaces:

**1. Java Macro API (Primary)**

**2. Python (via Jython/JPype)**

**3. Command Line Interface**

**4. Batch Mode**

### Authentication and Licensing

Star-CCM+ uses FlexLM license management:

**License Configuration:**

```bash
# Set license server
export CDLMD_LICENSE_FILE=1999@license-server.company.com

# Or in license.dat file
SERVER license-server.company.com ANY 1999
VENDOR siemens
```

**Check License:**

```bash
# Check license status
$STAR_ROOT/license/lmutil lmstat -c $CDLMD_LICENSE_FILE -a

# Check available features
$STAR_ROOT/license/lmutil lmstat -c $CDLMD_LICENSE_FILE -f
```

### Example Automation Scripts

**Java Macro - Complete CFD Setup:**

```java
// Star-CCM+ Java macro example
// Automates pipe flow simulation setup

package macro;

import star.common.*;
import star.base.neo.*;
import star.flow.*;
import star.turbulence.*;
import star.energy.*;
import star.segregatedflow.*;

public class PipeFlowAutomation extends StarMacro {

    public void execute() {
        // Get simulation object
        Simulation sim = getActiveSimulation();

        // Print banner
        sim.println("=== Pipe Flow Automation ===");

        // Setup physics
        setupPhysics(sim);

        // Import geometry
        importGeometry(sim, "pipe.stp");

        // Create regions
        createRegion(sim);

        // Setup boundaries
        setupBoundaries(sim);

        // Create mesh
        createMesh(sim);

        // Setup solvers
        setupSolvers(sim);

        // Run simulation
        runSimulation(sim, 1000);

        // Post-process
        postProcess(sim);

        // Export results
        exportResults(sim);

        sim.println("=== Automation Complete ===");
    }

    private void setupPhysics(Simulation sim) {
        sim.println("Setting up physics models...");

        // Create physics continuum
        PhysicsContinuum physics =
            sim.getContinuumManager().createContinuum(PhysicsContinuum.class);
        physics.setPresentationName("Water");

        // Enable models
        physics.enable(ThreeDimensionalModel.class);
        physics.enable(SteadyModel.class);
        physics.enable(SingleComponentLiquidModel.class);
        physics.enable(SegregatedFlowModel.class);
        physics.enable(ConstantDensityModel.class);
        physics.enable(TurbulentModel.class);
        physics.enable(RansTurbulenceModel.class);
        physics.enable(KEpsilonTurbulence.class);
        physics.enable(RkeTwoLayerTurbModel.class);
        physics.enable(KeTwoLayerAllYplusWallTreatment.class);

        sim.println("Physics models configured");
    }

    private void importGeometry(Simulation sim, String filename) {
        sim.println("Importing geometry: " + filename);

        // Import CAD
        sim.get(PartImportManager.class)
            .importCadPart(
                resolvePath(filename),
                "SharpEdges", 30.0, 1, true, 1.0E-5, true, false
            );
    }

    private void createRegion(Simulation sim) {
        sim.println("Creating region...");

        // Get parts
        MeshPartFactory meshPartFactory =
            sim.get(MeshPartFactory.class);

        // Create region from parts
        sim.getRegionManager().newRegionsFromParts(
            new NeoObjectVector(new Object[] {}),
            "OneRegionPerPart",
            null, "OneBoundaryPerPartSurface",
            null, "OneFeatureCurve",
            null, null, "CreateBoundaryLayer"
        );
    }

    private void setupBoundaries(Simulation sim) {
        sim.println("Setting up boundary conditions...");

        // Get region
        Region region = sim.getRegionManager().getRegion("pipe");

        // Inlet boundary
        Boundary inlet = region.getBoundaryManager().getBoundary("inlet");
        VelocityInletProfile velInlet =
            inlet.getValues().get(VelocityInletProfile.class);
        velInlet.getMethod(ConstantScalarProfileMethod.class)
            .getQuantity().setValue(2.0);  // 2 m/s

        // Outlet boundary
        Boundary outlet = region.getBoundaryManager().getBoundary("outlet");
        outlet.setBoundaryType(PressureOutlet.class);

        // Wall boundaries (no-slip)
        Boundary wall = region.getBoundaryManager().getBoundary("wall");
        wall.setBoundaryType(WallBoundary.class);

        sim.println("Boundary conditions configured");
    }

    private void createMesh(Simulation sim) {
        sim.println("Creating mesh...");

        // Get region
        Region region = sim.getRegionManager().getRegion("pipe");

        // Set base size
        region.getValues().get(BaseSize.class).setValue(0.01);  // 10 mm

        // Create automated mesh operation
        AutoMeshOperation meshOp =
            sim.get(MeshOperationManager.class)
                .createAutoMeshOperation(
                    new StringVector(new String[] {}),
                    new NeoObjectVector(new Object[] {region})
                );

        // Configure mesh settings
        meshOp.getMeshers().setMeshers(
            new NeoObjectVector(new Object[] {
                sim.get(MesherManager.class).getMesher(TrimmerMesher.class),
                sim.get(MesherManager.class).getMesher(PrismMesher.class)
            })
        );

        // Set surface size
        meshOp.getDefaultValues().get(BaseSize.class).setValue(0.01);

        // Set prism layer settings
        PrismAutoMesher prismMesher =
            meshOp.getMeshers().getPrismMesher();
        prismMesher.getPrismStretchingFunction()
            .setSelected(PrismStretchingOption.Type.WALL_THICKNESS);

        // Generate mesh
        meshOp.execute();

        // Get mesh statistics
        MeshMetrics metrics =
            region.getMeshRepresentation().getMeshMetrics();
        int cellCount = metrics.getCellCount();

        sim.println("Mesh created: " + cellCount + " cells");
    }

    private void setupSolvers(Simulation sim) {
        sim.println("Configuring solvers...");

        // Get solver
        SegregatedFlowSolver flowSolver =
            sim.getSolverManager().getSolver(SegregatedFlowSolver.class);

        // Set under-relaxation
        flowSolver.getUrf().setValue(0.7);

        // Pressure solver
        PressureSolver pressureSolver =
            sim.getSolverManager().getSolver(PressureSolver.class);
        pressureSolver.getUrf().setValue(0.3);

        sim.println("Solvers configured");
    }

    private void runSimulation(Simulation sim, int iterations) {
        sim.println("Running simulation for " + iterations + " iterations...");

        // Set stopping criterion
        StepStoppingCriterion stepCriterion =
            sim.getSolverStoppingCriterionManager()
                .createSolverStoppingCriterion(StepStoppingCriterion.class);
        stepCriterion.setMaximumNumberSteps(iterations);

        // Run
        sim.getSimulationIterator().run();

        sim.println("Simulation complete");
    }

    private void postProcess(Simulation sim) {
        sim.println("Post-processing...");

        // Create pressure coefficient report
        PressureCoefficientReport pressureReport =
            sim.getReportManager().createReport(PressureCoefficientReport.class);
        pressureReport.setPresentationName("Pressure Drop");

        // Create velocity magnitude report
        MaxReport maxVelReport =
            sim.getReportManager().createReport(MaxReport.class);
        maxVelReport.setPresentationName("Max Velocity");

        sim.println("Post-processing complete");
    }

    private void exportResults(Simulation sim) {
        sim.println("Exporting results...");

        // Save simulation
        String simFile = sim.getSessionDir() + "/pipe_flow.sim";
        sim.saveState(simFile);

        // Export scene images
        for (Scene scene : sim.getSceneManager().getScenes()) {
            scene.printAndWait(
                resolvePath(scene.getPresentationName() + ".png"),
                2, 1920, 1080, true, false
            );
        }

        // Export data tables
        for (XYPlot plot : sim.getPlotManager().getPlots()) {
            plot.export(resolvePath(plot.getPresentationName() + ".csv"), ",");
        }

        sim.println("Results exported");
    }
}
```

**Python Wrapper for Star-CCM+:**

```python
#!/usr/bin/env python3
"""
Star-CCM+ Python automation wrapper
Executes Java macros and processes results
"""

import subprocess
import os
import time
import pandas as pd

class StarCCMAutomation:
    def __init__(self, starccm_path='/opt/Siemens/17.06.007/STAR-CCM+17.06.007/star/bin/starccm+'):
        """Initialize Star-CCM+ automation"""
        self.starccm_path = starccm_path
        self.license_server = os.environ.get('CDLMD_LICENSE_FILE')

        if not os.path.exists(self.starccm_path):
            raise FileNotFoundError(f"Star-CCM+ not found at {self.starccm_path}")

    def run_macro(self, sim_file, macro_file, batch=True,
                  np=4, output_file=None):
        """Run Star-CCM+ with Java macro"""
        print(f"Running Star-CCM+ macro: {macro_file}")
        print(f"Simulation file: {sim_file}")
        print(f"Number of processors: {np}")

        cmd = [
            self.starccm_path,
            '-batch', macro_file if batch else '',
            '-np', str(np),
            sim_file
        ]

        # Remove empty strings
        cmd = [x for x in cmd if x]

        # Set environment
        env = os.environ.copy()
        if self.license_server:
            env['CDLMD_LICENSE_FILE'] = self.license_server

        # Run command
        start_time = time.time()

        if output_file:
            with open(output_file, 'w') as f:
                process = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    env=env,
                    text=True
                )
        else:
            process = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                env=env
            )
            print(process.stdout)

        elapsed = time.time() - start_time

        if process.returncode != 0:
            print(f"Error running Star-CCM+:")
            print(process.stderr)
            raise RuntimeError("Star-CCM+ execution failed")

        print(f"Star-CCM+ completed in {elapsed:.1f} seconds")

        return process

    def create_new_sim(self, sim_file, macro_file=None):
        """Create new simulation file"""
        print(f"Creating new simulation: {sim_file}")

        cmd = [
            self.starccm_path,
            '-new',
            sim_file
        ]

        if macro_file:
            cmd.extend(['-batch', macro_file])

        subprocess.run(cmd)

    def run_parametric_study(self, base_sim, macro_file,
                            parameters, output_dir):
        """Run parametric study"""
        print("Running parametric study...")

        os.makedirs(output_dir, exist_ok=True)
        results = []

        for i, params in enumerate(parameters):
            print(f"\n{'='*60}")
            print(f"Case {i+1}/{len(parameters)}")
            print(f"Parameters: {params}")
            print(f"{'='*60}")

            # Create case directory
            case_dir = os.path.join(output_dir, f"case_{i:03d}")
            os.makedirs(case_dir, exist_ok=True)

            # Copy simulation file
            case_sim = os.path.join(case_dir, os.path.basename(base_sim))
            if not os.path.exists(case_sim):
                import shutil
                shutil.copy(base_sim, case_sim)

            # Create parameter file
            param_file = os.path.join(case_dir, 'parameters.txt')
            with open(param_file, 'w') as f:
                for key, value in params.items():
                    f.write(f"{key}={value}\n")

            # Run simulation
            output_file = os.path.join(case_dir, 'output.log')
            try:
                self.run_macro(case_sim, macro_file,
                              output_file=output_file)

                result = params.copy()
                result['status'] = 'success'
                result['case_dir'] = case_dir

            except Exception as e:
                print(f"Error in case {i}: {e}")
                result = params.copy()
                result['status'] = 'failed'
                result['error'] = str(e)

            results.append(result)

        # Save results summary
        results_df = pd.DataFrame(results)
        results_df.to_csv(os.path.join(output_dir, 'results_summary.csv'),
                         index=False)

        print("\nParametric study complete!")
        print(f"Results saved to: {output_dir}/results_summary.csv")

        return results_df

    def extract_residuals(self, sim_file):
        """Extract residuals from simulation"""
        residual_file = sim_file.replace('.sim', '_residuals.csv')

        if os.path.exists(residual_file):
            df = pd.read_csv(residual_file)
            return df
        else:
            print(f"Residual file not found: {residual_file}")
            return None

    def check_convergence(self, residuals, threshold=1e-4):
        """Check if simulation converged"""
        if residuals is None:
            return False

        # Check if final residuals are below threshold
        final_residuals = residuals.iloc[-1]

        converged = all(final_residuals < threshold)

        return converged


# Example usage
if __name__ == "__main__":
    # Initialize automation
    starccm = StarCCMAutomation()

    # Example 1: Run single simulation
    print("=== Single Simulation ===")
    starccm.run_macro(
        sim_file='pipe_flow.sim',
        macro_file='PipeFlowAutomation.java',
        np=8,
        output_file='simulation.log'
    )

    # Example 2: Parametric study
    print("\n=== Parametric Study ===")
    parameters = [
        {'velocity': 1.0, 'temperature': 300},
        {'velocity': 2.0, 'temperature': 300},
        {'velocity': 3.0, 'temperature': 300},
        {'velocity': 1.0, 'temperature': 350},
        {'velocity': 2.0, 'temperature': 350},
    ]

    results = starccm.run_parametric_study(
        base_sim='base_case.sim',
        macro_file='parametric_macro.java',
        parameters=parameters,
        output_dir='parametric_study'
    )

    print("\nResults:")
    print(results)
```

### File I/O and Data Exchange

**Export Data from Star-CCM+:**

```java
// Export surface data
XYPlot plot = sim.getPlotManager().getPlot("Pressure Distribution");
plot.export(resolvePath("pressure_data.csv"), ",");

// Export field data
FieldFunction pressureField = sim.getFieldFunctionManager()
    .getFunction("Pressure");

Region region = sim.getRegionManager().getRegion("Fluid");

// Create solution view
SolutionView view = sim.get(SolutionViewManager.class)
    .createSolutionView();
view.setPresentationName("Export View");

// Export
view.export(resolvePath("field_data.csv"), ",");
```

**Import Geometry:**

```java
// Import CAD (STEP, IGES, etc.)
sim.get(PartImportManager.class).importCadPart(
    resolvePath("geometry.stp"),
    "SharpEdges",
    30.0,  // Feature angle
    1,     // Create part per body
    true,  // CAD topology
    1.0E-5,// Cleanup tolerance
    true,  // Assign part to regions
    false  // Extract 2D surfaces
);
```

### Best Practices

1. **Macro Development:**
   - Use Star-CCM+ macro recorder to learn API
   - Test macros on simple cases first
   - Add error handling and validation
   - Log all operations for debugging

2. **Performance:**
   - Use appropriate number of processors
   - Enable parallel meshing when possible
   - Optimize solver settings for problem type
   - Monitor memory usage

3. **Batch Processing:**
   - Run in batch mode for automation
   - Capture all output to log files
   - Implement checkpointing for long runs
   - Clean up temporary files

4. **Results Management:**
   - Export results systematically
   - Use consistent naming conventions
   - Archive simulation files with results
   - Document simulation parameters

5. **Licensing:**
   - Check license availability before batch jobs
   - Use appropriate license tokens
   - Release licenses promptly
   - Monitor license usage across team

---

## Summary

This document provides comprehensive integration guides for six major engineering software packages. Key takeaways:

1. **Automation Benefits:**
   - Increased productivity through batch processing
   - Reproducible workflows
   - Parametric studies and optimization
   - Integration with other tools

2. **Common Patterns:**
   - Most tools provide Python/scripting APIs
   - License management via FlexLM
   - File-based data exchange (CSV, VTK, etc.)
   - Batch mode for automation

3. **Best Practices:**
   - Always validate inputs
   - Implement error handling
   - Log all operations
   - Document workflows
   - Use version control

4. **Integration Strategies:**
   - Use native APIs when available
   - Leverage file I/O for tool coupling
   - Implement workflow orchestration
   - Automate post-processing

For specific implementation details, refer to official documentation for each tool.
