---
name: comsol-multiphysics
description: "Set up coupled fluid-structure interaction for pump vibration analysis"
category: integrations
domain: multiphysics
complexity: advanced
dependencies: []
---

# COMSOL Multiphysics Integration

Comprehensive guide for setting up and automating coupled fluid-structure interaction (FSI) simulations using COMSOL Multiphysics, with a focus on pump vibration analysis and related multiphysics applications.

## Overview

COMSOL Multiphysics is a leading commercial simulation platform for modeling and solving complex multiphysics problems. It provides:

- **Unified Environment**: Single platform for multiple physics domains
- **Flexible Coupling**: Built-in tools for coupling different physics
- **Equation-Based Modeling**: Direct access to underlying PDEs
- **Java API**: Programmatic control and automation
- **MATLAB LiveLink**: Integration with MATLAB for pre/post-processing
- **Application Builder**: Create custom simulation apps

COMSOL excels at coupled physics problems including:
- Fluid-structure interaction (FSI)
- Thermal-structural coupling
- Electromagnetics with thermal effects
- Acoustics-structure interaction
- Electrochemistry with transport phenomena

## COMSOL Architecture

### Model Structure
- **Component**: Independent geometry and physics
- **Study**: Solution sequence with one or more steps
- **Physics Interfaces**: Pre-configured PDE systems for specific phenomena
- **Multiphysics Couplings**: Automatic or manual coupling between physics
- **Solver Configurations**: Automatic or customized solver sequences

### Solver Technology
- **COMSOL Solver**: Direct and iterative solvers
- **Fully Coupled Approach**: Simultaneous solution of all equations
- **Segregated Solver**: Physics-based segregation for large models
- **Time-Stepping**: Implicit methods for transient analysis
- **Mesh Adaptation**: Automatic refinement based on solution

## Licensing Requirements

### Commercial Licensing

COMSOL Multiphysics requires commercial licenses for all production use:

- **Base Package**:
  - COMSOL Multiphysics license (required for all users)
  - Includes basic PDE interfaces and computational tools
  - License managed via FlexNet License Server

- **Module Licenses**:
  - Individual modules purchased separately
  - Floating licenses (checked out during use)
  - Batch mode requires license checkout
  - HPC add-ons for parallel computing

- **License Server Setup**:
  ```bash
  # Set license server environment variable (Linux/Mac)
  export LMCOMSOL_LICENSE_FILE=1718@license-server.company.com

  # Or in Windows
  set LMCOMSOL_LICENSE_FILE=1718@license-server.company.com
  ```

- **License Types**:
  - **Floating Network License**: Shared among users
  - **Node-Locked License**: Tied to specific machine
  - **HPC License**: Additional cores for parallel computing
  - **Batch License**: For batch/automated simulations

### Academic Licensing
- Classroom Kit: Limited to educational use
- Research licenses: Available for academic institutions
- Restrictions on commercial applications
- May have feature or model size limitations

### Important Licensing Notes
- **License Check**: Always verify license availability before batch runs
- **License Release**: Properly close COMSOL to release licenses
- **Batch Operations**: Require special batch/HPC licenses
- **Module Dependencies**: Some modules require others (e.g., FSI needs CFD + Structural)
- **API Access**: Requires base license plus relevant modules

## Key Modules for Pump Applications

### CFD Module

Comprehensive computational fluid dynamics for single-phase and multiphase flows:

**Capabilities:**
- **Turbulence Models**:
  - RANS: k-epsilon, k-omega, SST
  - LES and DES for high-fidelity simulations
  - Wall functions and low-Re formulations

- **Flow Types**:
  - Incompressible and compressible
  - Laminar and turbulent
  - Steady-state and transient
  - Rotating machinery (frozen rotor, sliding mesh)

- **Multiphase Flow**:
  - Euler-Euler multiphase
  - Phase field methods
  - Level set tracking
  - Bubbly flow models

**Pump-Specific Features:**
- Rotating reference frames for impellers
- Mixing plane interfaces
- Cavitation modeling
- Pressure pulsation analysis

### Structural Mechanics Module

Comprehensive structural analysis including linear and nonlinear behavior:

**Capabilities:**
- **Analysis Types**:
  - Static structural analysis
  - Eigenfrequency analysis (modal analysis)
  - Frequency response analysis
  - Transient dynamics
  - Prestressed analysis

- **Material Models**:
  - Linear elastic
  - Hyperelastic (rubber, polymers)
  - Plasticity (metals under high stress)
  - Composite materials
  - Contact and friction

- **Dynamic Analysis**:
  - Modal analysis for natural frequencies
  - Harmonic response
  - Time-dependent loading
  - Damping models (Rayleigh, modal)

**Pump-Specific Features:**
- Rotating machinery stress analysis
- Bolt preload and assembly stress
- Fatigue analysis
- Contact between impeller and casing

### Fluid-Structure Interaction (FSI)

Couples CFD and Structural Mechanics for two-way interaction:

**Coupling Approaches:**

1. **One-Way FSI**:
   - Fluid loads applied to structure
   - Structure does not affect fluid
   - Faster computation
   - Suitable for rigid-like structures

2. **Two-Way FSI**:
   - Fully coupled fluid and structure
   - Deforming mesh methods
   - Accounts for large deformations
   - Required for flexible structures

3. **Weak Coupling**:
   - Sequential solution (staggered approach)
   - Fluid → Structure → Fluid iteration
   - Better for loosely coupled problems

4. **Strong Coupling**:
   - Fully implicit simultaneous solution
   - Better convergence for tightly coupled problems
   - Higher computational cost

**FSI Features:**
- Automatic mesh deformation (ALE method)
- Remeshing for large deformations
- Fluid loads computed automatically
- Pressure and viscous forces transferred
- Support for multiple solid bodies

**Pump FSI Applications:**
- Impeller blade vibration under fluid forces
- Casing vibration and acoustic radiation
- Seal deflection under pressure
- Shaft deflection and critical speeds
- Cavitation-induced vibration

### Additional Relevant Modules

**Acoustics Module:**
- Pressure pulsation analysis
- Noise radiation from pump casing
- Structure-borne noise
- Coupled acoustic-structure-fluid problems

**Heat Transfer Module:**
- Thermal loads in pumps
- Coupled thermal-structural analysis
- Conjugate heat transfer (fluid-solid)
- Thermal expansion effects

**Optimization Module:**
- Design optimization
- Topology optimization
- Shape optimization
- Parameter sweeps

## Java API for Automation

COMSOL provides a comprehensive Java API for programmatic model building and automation:

### API Structure

**Core Components:**
- **Model object**: Top-level container
- **Component**: Geometry and physics
- **Physics interfaces**: Add and configure physics
- **Study**: Define solution procedures
- **Results**: Post-processing and visualization

### Basic Java API Workflow

```java
import com.comsol.model.*;
import com.comsol.model.util.*;

public class PumpFSI {
    public static Model run() {
        // Create model
        Model model = ModelUtil.create("PumpFSI");

        // Create component
        model.component().create("comp1", true);

        // Create geometry
        model.component("comp1").geom().create("geom1", 3);

        // Import CAD geometry
        model.component("comp1").geom("geom1").create("imp1", "Import");
        model.component("comp1").geom("geom1").feature("imp1")
            .set("filename", "/path/to/pump_geometry.step");
        model.component("comp1").geom("geom1").run();

        // Add fluid physics (CFD)
        model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");
        model.component("comp1").physics("spf").selection()
            .named("geom1_fluid_domain");

        // Add structural physics
        model.component("comp1").physics().create("solid", "SolidMechanics", "geom1");
        model.component("comp1").physics("solid").selection()
            .named("geom1_solid_domain");

        // Add FSI coupling
        model.component("comp1").multiphysics().create("fsi1", "FluidStructureInteraction", 3);
        model.component("comp1").multiphysics("fsi1")
            .selection().named("geom1_fsi_boundary");

        // Create mesh
        model.component("comp1").mesh().create("mesh1");
        model.component("comp1").mesh("mesh1").automatic(true);
        model.component("comp1").mesh("mesh1").run();

        // Create study
        model.study().create("std1");
        model.study("std1").create("time", "Transient");
        model.study("std1").feature("time").set("tlist", "range(0,0.01,1)");

        // Solve
        model.sol().create("sol1");
        model.sol("sol1").study("std1");
        model.sol("sol1").feature().create("st1", "StudyStep");
        model.sol("sol1").feature().create("v1", "Variables");
        model.sol("sol1").feature().create("t1", "Time");
        model.sol("sol1").attach("std1");
        model.sol("sol1").runAll();

        // Save model
        model.save("/path/to/pump_fsi_model.mph");

        return model;
    }

    public static void main(String[] args) {
        run();
    }
}
```

### Running Java API Scripts

```bash
# Compile Java file
comsol compile PumpFSI.java

# Run with COMSOL
comsol batch -inputfile PumpFSI.class -outputfile results.mph

# Or run directly
java -cp /path/to/comsol/plugins/*:. PumpFSI
```

### API Advantages
- **Reproducibility**: Scripts ensure consistent model building
- **Parametric Studies**: Easy to vary parameters
- **Batch Processing**: Run multiple cases automatically
- **Integration**: Connect with other tools and databases
- **Version Control**: Track model changes in source control

## MATLAB LiveLink

COMSOL integrates seamlessly with MATLAB for enhanced scripting and data processing:

### LiveLink Features

**Model Control from MATLAB:**
- Build and modify COMSOL models
- Run simulations from MATLAB scripts
- Extract results into MATLAB workspace
- Use MATLAB's data analysis tools

**Installation:**
- Requires separate LiveLink for MATLAB license
- Configure MATLAB path to COMSOL installation
- Start COMSOL server or use direct connection

### Basic MATLAB LiveLink Usage

```matlab
% Initialize COMSOL with MATLAB
import com.comsol.model.*
import com.comsol.model.util.*

% Start COMSOL server (if not already running)
mphstart

% Create or load model
model = mphload('pump_model.mph');

% Modify parameters
model.param.set('inlet_velocity', '5[m/s]');
model.param.set('outlet_pressure', '101325[Pa]');

% Run study
model.study('std1').run();

% Extract results
pressure = mpheval(model, 'p', 'dataset', 'dset1');
velocity = mpheval(model, 'u', 'dataset', 'dset1');

% Process in MATLAB
mean_pressure = mean(pressure.d1);
max_velocity = max(sqrt(velocity.d1.^2 + velocity.d2.^2 + velocity.d3.^2));

% Plot using MATLAB
figure;
plot(pressure.p, pressure.d1);
xlabel('Position');
ylabel('Pressure [Pa]');
title('Pressure Distribution');

% Save results
save('pump_results.mat', 'pressure', 'velocity');

% Close COMSOL
ModelUtil.remove('model');
```

### Parametric Study with MATLAB

```matlab
% Parametric study of inlet velocity effects
velocities = 1:1:10;  % m/s
results = struct();

for i = 1:length(velocities)
    fprintf('Running case %d: velocity = %.1f m/s\n', i, velocities(i));

    % Set parameter
    model.param.set('inlet_velocity', sprintf('%f[m/s]', velocities(i)));

    % Solve
    model.study('std1').run();

    % Extract force on impeller
    force = mphint2(model, 'spf.Fp_x', 'surface', 'selection', 5);
    results(i).velocity = velocities(i);
    results(i).force = force;

    % Extract vibration amplitude
    displacement = mphmax(model, 'sqrt(u^2+v^2+w^2)', 'volume', 'selection', 3);
    results(i).max_displacement = displacement;
end

% Plot results
figure;
subplot(2,1,1);
plot([results.velocity], [results.force], '-o');
xlabel('Inlet Velocity [m/s]');
ylabel('Force on Impeller [N]');
grid on;

subplot(2,1,2);
plot([results.velocity], [results.max_displacement]*1e6, '-o');
xlabel('Inlet Velocity [m/s]');
ylabel('Max Displacement [μm]');
grid on;

% Save results
save('parametric_results.mat', 'results');
```

## Common Workflows

### Workflow 1: Pump Casing Vibration Under Fluid Loads

**Application**: Analyze vibration of pump casing due to pressure pulsations from fluid flow.

**Approach**: One-way FSI (fluid loads mapped to structure)

**Steps:**

1. **Geometry Setup**:
   - Import pump casing geometry
   - Create fluid domain inside casing
   - Define FSI boundary (fluid-structure interface)

2. **CFD Setup**:
   - Define inlet and outlet boundaries
   - Set up turbulence model (k-epsilon or k-omega SST)
   - Configure transient solver
   - Apply rotating reference frame if analyzing flow with rotating impeller

3. **Structural Setup**:
   - Define material properties (steel, cast iron)
   - Apply boundary conditions (fixed support at mounting points)
   - Set up transient structural analysis
   - Include damping if known

4. **FSI Coupling**:
   - Map fluid pressure and shear stress to structure
   - One-way coupling (structure assumed rigid relative to fluid)
   - Time synchronization between fluid and structural solvers

5. **Solution**:
   - Solve transient CFD first
   - Extract time-varying pressure loads
   - Apply loads to structural model
   - Solve structural dynamics

6. **Post-Processing**:
   - Extract displacement at critical points
   - FFT analysis for frequency content
   - Compare with natural frequencies
   - Stress concentration analysis

### Workflow 2: Impeller FSI Analysis

**Application**: Coupled analysis of impeller blade deformation under fluid forces.

**Approach**: Two-way FSI with moving mesh

**Steps:**

1. **Geometry and Mesh**:
   - Import impeller geometry
   - Create fluid domain around impeller
   - Define FSI interface at blade surfaces
   - Mesh solid impeller and fluid domain
   - Fine mesh at FSI boundary

2. **Fluid Physics**:
   - Turbulent flow (k-omega SST recommended)
   - Rotating reference frame for impeller
   - Inlet velocity or mass flow boundary
   - Outlet pressure boundary
   - Moving mesh (ALE) enabled

3. **Structural Physics**:
   - Linear elastic material initially
   - Fixed at shaft connection
   - Prestress from centrifugal loading
   - Include rotation effects

4. **Two-Way FSI Coupling**:
   - Fluid applies loads to structure
   - Structure deformation moves mesh
   - Iterative coupling within each time step
   - Convergence criteria for FSI iteration

5. **Prestressed FSI**:
   - First: Stationary structural analysis with centrifugal load
   - Second: Use prestressed state as initial condition for FSI
   - Improves convergence and accuracy

6. **Solution Strategy**:
   - Start with coarse time steps for initial transient
   - Refine time stepping for periodic solution
   - Monitor residuals and FSI iterations
   - Expect longer solve times for strong coupling

7. **Post-Processing**:
   - Blade tip displacement vs. time
   - Stress distribution on blades
   - Effect of deformation on flow field
   - Comparison of rigid vs. flexible results

### Workflow 3: Coupled Thermal-Flow Analysis

**Application**: Thermal effects in pumps handling hot fluids.

**Approach**: Conjugate heat transfer with thermal expansion

**Steps:**

1. **Multi-Domain Setup**:
   - Fluid domain (hot liquid)
   - Solid domains (casing, impeller)
   - Define fluid-solid boundaries

2. **Coupled Physics**:
   - CFD for fluid flow
   - Heat Transfer in Fluids (energy equation)
   - Heat Transfer in Solids (conduction)
   - Thermal Stress in solids
   - Automatic temperature coupling at boundaries

3. **Boundary Conditions**:
   - Inlet: Temperature and velocity
   - Outlet: Pressure outflow
   - External surfaces: Convection or radiation
   - Solid boundaries: Fixed temperature or insulation

4. **Material Properties**:
   - Temperature-dependent viscosity
   - Thermal expansion coefficients
   - Thermal conductivity

5. **Solution Sequence**:
   - Steady-state thermal-flow first
   - Use as initial condition for thermal stress
   - Or fully coupled transient simulation

6. **Post-Processing**:
   - Temperature distribution
   - Thermal stresses
   - Thermal expansion displacements
   - Effect on clearances and fits

### Workflow 4: Modal Analysis with Fluid Loading

**Application**: Natural frequencies of pump casing in contact with fluid.

**Approach**: Eigenfrequency analysis with fluid-structure coupling

**Steps:**

1. **Geometry**:
   - Pump casing structure
   - Surrounding or contained fluid
   - FSI boundaries

2. **Physics Setup**:
   - Pressure Acoustics for fluid
   - Solid Mechanics for structure
   - Acoustic-Structure Boundary coupling

3. **Boundary Conditions**:
   - Structural: Fixed supports
   - Acoustic: Sound hard walls at rigid boundaries
   - FSI: Acceleration coupling

4. **Study**:
   - Eigenfrequency study
   - Search for eigenfrequencies in range of interest
   - Acoustic-structure coupling included automatically

5. **Results**:
   - Natural frequencies with fluid loading
   - Mode shapes
   - Comparison with in-air and in-vacuum modes
   - Fluid added mass effect

## Batch Execution

### Command Line Execution

COMSOL can run in batch mode without GUI:

```bash
# Run existing model file
comsol batch -inputfile pump_model.mph -outputfile results.mph

# Run Java method file
comsol batch -inputfile PumpFSI.java -outputfile results.mph

# Run with specific study
comsol batch -inputfile pump_model.mph -outputfile results.mph -study std1

# Run with parameter override
comsol batch -inputfile pump_model.mph -outputfile results.mph \
  -pname inlet_velocity -plist 5,10,15,20

# Run MATLAB script
comsol batch -inputfile pump_analysis.m -outputfile results.mph

# Specify number of processors
comsol batch -np 8 -inputfile pump_model.mph -outputfile results.mph
```

### Batch Script Example (Linux)

```bash
#!/bin/bash
# batch_comsol.sh

# Set license server
export LMCOMSOL_LICENSE_FILE=1718@license-server.com

# Set number of processors
NPROCS=8

# Input/output files
INPUT_MODEL="pump_fsi_base.mph"
OUTPUT_DIR="results"

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Parameter sweep: inlet velocities
VELOCITIES="2 4 6 8 10"

for VEL in ${VELOCITIES}; do
    echo "Running simulation with inlet velocity = ${VEL} m/s"

    OUTPUT_FILE="${OUTPUT_DIR}/pump_fsi_v${VEL}.mph"
    LOG_FILE="${OUTPUT_DIR}/pump_fsi_v${VEL}.log"

    # Run COMSOL in batch mode
    comsol batch -np ${NPROCS} \
        -inputfile ${INPUT_MODEL} \
        -outputfile ${OUTPUT_FILE} \
        -pname inlet_velocity \
        -plist ${VEL} \
        > ${LOG_FILE} 2>&1

    # Check exit status
    if [ $? -eq 0 ]; then
        echo "  Simulation completed successfully"
    else
        echo "  ERROR: Simulation failed - check ${LOG_FILE}"
        exit 1
    fi
done

echo "All simulations completed"
```

### Python Automation with JPype

Use JPype to call COMSOL Java API from Python:

```python
import jpype
import jpype.imports
from jpype.types import *
import numpy as np

# Start JVM with COMSOL
comsol_root = '/usr/local/comsol56/multiphysics'
jvm_path = jpype.getDefaultJVMPath()

jpype.startJVM(
    jvm_path,
    f"-Djava.library.path={comsol_root}/lib/glnxa64",
    f"-Dcs.lic.path={comsol_root}/license",
    classpath=f"{comsol_root}/plugins/*"
)

# Import COMSOL packages
from com.comsol.model.util import ModelUtil

# Load model
model = ModelUtil.load('/path/to/pump_model.mph')

# Parametric study
velocities = np.linspace(2, 10, 5)
results = []

for vel in velocities:
    print(f"Running velocity = {vel:.1f} m/s")

    # Set parameter
    model.param().set('inlet_velocity', f'{vel}[m/s]')

    # Run study
    model.study('std1').run()

    # Extract results (example)
    # Add actual result extraction here

    # Save case
    model.save(f'results/pump_v{vel:.1f}.mph')

    results.append({
        'velocity': vel,
        # Add extracted results here
    })

# Cleanup
ModelUtil.remove('model')
jpype.shutdownJVM()

print("Parametric study completed")
```

## Best Practices

### Model Setup
1. **Geometry**: Clean, defeatured geometry suitable for meshing
2. **Selections**: Use named selections for boundary conditions
3. **Parameters**: Use global parameters for easy modification
4. **Units**: Consistent unit system throughout model
5. **Symmetry**: Exploit symmetry when possible

### Meshing
1. **Mesh Quality**: Check mesh quality metrics (skewness, growth rate)
2. **Boundary Layers**: Critical for turbulent flow
3. **FSI Interface**: Fine mesh at fluid-structure boundary
4. **Mesh Independence**: Verify solution mesh convergence
5. **Adaptive Meshing**: Use for complex flows

### Physics Setup
1. **Appropriate Models**: Select physics models suitable for application
2. **Material Properties**: Verify material data accuracy
3. **Boundary Conditions**: Ensure physical consistency
4. **Initial Conditions**: Provide good initial guess
5. **Coupling Settings**: Check FSI coupling parameters

### Solver Strategy
1. **Solver Selection**: Automatic vs. manual solver configuration
2. **Convergence Criteria**: Appropriate tolerances
3. **Time Stepping**: Adaptive time stepping for transient
4. **Linearization**: Check nonlinear solver settings
5. **Parallel Computing**: Use multiple cores effectively

### Automation
1. **Version Control**: Track model files and scripts
2. **Documentation**: Comment code and document workflow
3. **Error Handling**: Implement checks and error recovery
4. **Testing**: Validate automation on simple cases first
5. **License Management**: Check license availability in batch jobs

### Post-Processing
1. **Data Extraction**: Export only necessary data
2. **Visualization**: Standardized plots for comparison
3. **Validation**: Compare with experimental or analytical results
4. **Reporting**: Automated report generation
5. **Data Management**: Organize results systematically

## Troubleshooting

### Common Issues

**Licensing Problems:**
- Verify license server connectivity
- Check module availability
- Ensure batch licenses for automated runs

**Meshing Failures:**
- Simplify geometry (remove small features)
- Adjust mesh size parameters
- Use different mesh sequence

**Convergence Issues:**
- Reduce load step or time step
- Improve initial conditions
- Adjust solver tolerances
- Check for poor mesh quality
- Review boundary conditions

**FSI-Specific Issues:**
- Mesh deformation failures: Remesh or adjust smoothing
- Slow FSI convergence: Adjust relaxation factors
- Large deformations: Consider remeshing strategies
- Time step too large: Reduce for stability

**Performance Issues:**
- Use coarser mesh initially
- Simplify physics models
- Exploit symmetry
- Increase parallel cores
- Use segregated solver for large models

## Additional Resources

### Documentation
- **COMSOL Multiphysics User's Guide**: Core functionality
- **CFD Module User's Guide**: Fluid flow physics
- **Structural Mechanics Module User's Guide**: Solid mechanics
- **COMSOL API Reference**: Java API documentation
- **LiveLink for MATLAB User's Guide**: MATLAB integration

### Online Resources
- **COMSOL Website**: https://www.comsol.com
- **COMSOL Blog**: Application examples and tips
- **COMSOL Support Center**: Knowledge base and downloads
- **Application Gallery**: Example models and tutorials
- **COMSOL Forums**: Community discussions

### Training
- **COMSOL Training Courses**: Official training programs
- **Webinar Archive**: Recorded training sessions
- **Learning Center**: Video tutorials and guides

### Support
- **Technical Support**: support@comsol.com
- **Regional Support**: Contact local COMSOL office
- **Application Engineers**: Assistance with complex problems
