# Pump Fluid-Structure Interaction (FSI) Workflow

Complete workflow example for coupled fluid-structure interaction analysis of centrifugal pump impeller vibration.

## Overview

**Application**: Analyze vibration of a centrifugal pump impeller under coupled fluid-structure loading.

**Objectives**:
1. Compute fluid flow field through pump impeller
2. Calculate structural deformation of impeller blades
3. Analyze coupled FSI effects
4. Determine natural frequencies with fluid loading
5. Assess vibration amplitudes and stress levels

**Physics Involved**:
- Turbulent fluid flow (CFD Module)
- Structural dynamics (Structural Mechanics Module)
- Two-way FSI coupling
- Rotating machinery effects

## Problem Setup

### Geometry

**Components**:
- Impeller: 6-blade centrifugal design
- Fluid domain: Volume swept by impeller plus inlet/outlet regions
- Shaft connection: Hub region

**Dimensions** (example):
- Impeller diameter: 200 mm
- Blade height: 40 mm
- Shaft diameter: 30 mm
- Number of blades: 6

**CAD Import**:
- Import impeller geometry (STEP/IGES format)
- Create fluid domain using CFD design tools
- Define FSI interface at blade surfaces

### Operating Conditions

**Flow Parameters**:
- Flow rate: 100 L/s (0.1 m³/s)
- Rotational speed: 1500 rpm (157 rad/s)
- Working fluid: Water at 20°C
- Inlet pressure: Atmospheric (101325 Pa)
- Outlet pressure: 250 kPa (gauge)

**Material Properties**:
- Impeller material: Stainless steel 316
  - Density: 8000 kg/m³
  - Young's modulus: 193 GPa
  - Poisson's ratio: 0.27
  - Yield strength: 290 MPa
- Fluid: Water
  - Density: 998 kg/m³
  - Dynamic viscosity: 0.001 Pa·s

## Complete Java Implementation

```java
import com.comsol.model.*;
import com.comsol.model.util.*;

public class PumpFSIAnalysis {

    public static Model run() {
        // ============================================================
        // 1. MODEL SETUP
        // ============================================================

        Model model = ModelUtil.create("PumpFSI");
        model.label("Centrifugal Pump Impeller FSI Analysis");

        // Create component
        model.component().create("comp1", true);
        model.component("comp1").geom().create("geom1", 3);

        // ============================================================
        // 2. PARAMETERS
        // ============================================================

        model.param().set("rpm", "1500[rpm]", "Rotational speed");
        model.param().set("omega", "rpm*2*pi/60", "Angular velocity");
        model.param().set("Q_flow", "0.1[m^3/s]", "Volume flow rate");
        model.param().set("p_in", "0[Pa]", "Inlet pressure (gauge)");
        model.param().set("p_out", "250000[Pa]", "Outlet pressure (gauge)");

        // Material properties - Fluid
        model.param().set("rho_water", "998[kg/m^3]", "Water density");
        model.param().set("mu_water", "0.001[Pa*s]", "Water viscosity");

        // Material properties - Solid
        model.param().set("rho_steel", "8000[kg/m^3]", "Steel density");
        model.param().set("E_steel", "193e9[Pa]", "Young's modulus");
        model.param().set("nu_steel", "0.27", "Poisson's ratio");

        // Geometry parameters
        model.param().set("D_impeller", "0.2[m]", "Impeller diameter");
        model.param().set("D_shaft", "0.03[m]", "Shaft diameter");

        // ============================================================
        // 3. GEOMETRY
        // ============================================================

        // Import impeller geometry
        model.component("comp1").geom("geom1").create("imp1", "Import");
        model.component("comp1").geom("geom1").feature("imp1")
            .set("type", "cad");
        model.component("comp1").geom("geom1").feature("imp1")
            .set("filename", "/path/to/impeller_geometry.step");

        // Create fluid domain (simplified - actual would be more complex)
        model.component("comp1").geom("geom1").create("inlet_cylinder", "Cylinder");
        model.component("comp1").geom("geom1").feature("inlet_cylinder")
            .set("r", "D_impeller/2")
            .set("h", "0.1")
            .set("pos", new String[]{"0", "0", "-0.1"});

        model.component("comp1").geom("geom1").create("outlet_cylinder", "Cylinder");
        model.component("comp1").geom("geom1").feature("outlet_cylinder")
            .set("r", "D_impeller/2")
            .set("h", "0.1")
            .set("pos", new String[]{"0", "0", "0.05"});

        // Boolean operations to create fluid domain
        model.component("comp1").geom("geom1").create("uni1", "Union");
        model.component("comp1").geom("geom1").feature("uni1")
            .selection("input").set("inlet_cylinder", "outlet_cylinder");

        model.component("comp1").geom("geom1").create("dif1", "Difference");
        model.component("comp1").geom("geom1").feature("dif1")
            .selection("input").set("uni1");
        model.component("comp1").geom("geom1").feature("dif1")
            .selection("input2").set("imp1");

        // Run geometry
        model.component("comp1").geom("geom1").run();

        // ============================================================
        // 4. SELECTIONS (Named selections for boundaries/domains)
        // ============================================================

        // Create selection for fluid domain
        model.component("comp1").selection().create("fluid_domain", "Explicit");
        model.component("comp1").selection("fluid_domain").label("Fluid Domain");
        model.component("comp1").selection("fluid_domain").geom("geom1", 3);
        model.component("comp1").selection("fluid_domain").set(1);  // Adjust index

        // Create selection for solid domain (impeller)
        model.component("comp1").selection().create("solid_domain", "Explicit");
        model.component("comp1").selection("solid_domain").label("Impeller Solid");
        model.component("comp1").selection("solid_domain").geom("geom1", 3);
        model.component("comp1").selection("solid_domain").set(2);  // Adjust index

        // FSI boundary (blade surfaces)
        model.component("comp1").selection().create("fsi_boundary", "Explicit");
        model.component("comp1").selection("fsi_boundary").label("FSI Interface");
        model.component("comp1").selection("fsi_boundary").geom("geom1", 2);
        model.component("comp1").selection("fsi_boundary").set(10, 11, 12, 13, 14, 15);

        // Inlet boundary
        model.component("comp1").selection().create("inlet", "Explicit");
        model.component("comp1").selection("inlet").label("Inlet");
        model.component("comp1").selection("inlet").geom("geom1", 2);
        model.component("comp1").selection("inlet").set(1);

        // Outlet boundary
        model.component("comp1").selection().create("outlet", "Explicit");
        model.component("comp1").selection("outlet").label("Outlet");
        model.component("comp1").selection("outlet").geom("geom1", 2);
        model.component("comp1").selection("outlet").set(2);

        // Fixed boundary (shaft)
        model.component("comp1").selection().create("fixed_shaft", "Explicit");
        model.component("comp1").selection("fixed_shaft").label("Fixed Shaft");
        model.component("comp1").selection("fixed_shaft").geom("geom1", 2);
        model.component("comp1").selection("fixed_shaft").set(20);

        // ============================================================
        // 5. MATERIALS
        // ============================================================

        // Water material
        model.component("comp1").material().create("mat_water", "Common");
        model.component("comp1").material("mat_water").label("Water");
        model.component("comp1").material("mat_water").selection()
            .named("fluid_domain");
        model.component("comp1").material("mat_water").propertyGroup("def")
            .set("density", "rho_water");
        model.component("comp1").material("mat_water").propertyGroup("def")
            .set("dynamicviscosity", "mu_water");

        // Steel material
        model.component("comp1").material().create("mat_steel", "Common");
        model.component("comp1").material("mat_steel").label("Stainless Steel");
        model.component("comp1").material("mat_steel").selection()
            .named("solid_domain");
        model.component("comp1").material("mat_steel").propertyGroup("def")
            .set("density", "rho_steel");
        model.component("comp1").material("mat_steel").propertyGroup("def")
            .set("youngsmodulus", "E_steel");
        model.component("comp1").material("mat_steel").propertyGroup("def")
            .set("poissonsratio", "nu_steel");

        // ============================================================
        // 6. PHYSICS - TURBULENT FLOW (CFD)
        // ============================================================

        // Create turbulent flow physics (k-omega SST)
        model.component("comp1").physics().create("spf", "TurbulentFlowkOmega", "geom1");
        model.component("comp1").physics("spf").selection().named("fluid_domain");

        // Rotating frame (for impeller region)
        model.component("comp1").physics("spf").create("rot1", "RotatingDomain", 3);
        model.component("comp1").physics("spf").feature("rot1").selection()
            .named("fluid_domain");
        model.component("comp1").physics("spf").feature("rot1")
            .set("omega", "omega");
        model.component("comp1").physics("spf").feature("rot1")
            .set("axis", new String[]{"0", "0", "1"});
        model.component("comp1").physics("spf").feature("rot1")
            .set("axispoint", new String[]{"0", "0", "0"});

        // Inlet boundary condition
        model.component("comp1").physics("spf").create("inlet1", "InletBoundary", 2);
        model.component("comp1").physics("spf").feature("inlet1").selection()
            .named("inlet");
        model.component("comp1").physics("spf").feature("inlet1")
            .set("BoundaryCondition", "MeanVelocity");
        model.component("comp1").physics("spf").feature("inlet1")
            .set("Uav", "Q_flow/(pi*(D_impeller/2)^2)");

        // Outlet boundary condition
        model.component("comp1").physics("spf").create("outlet1", "OutletBoundary", 2);
        model.component("comp1").physics("spf").feature("outlet1").selection()
            .named("outlet");
        model.component("comp1").physics("spf").feature("outlet1")
            .set("p0", "p_out");

        // Wall boundary at FSI interface (handled by FSI coupling)
        model.component("comp1").physics("spf").create("wall_fsi", "WallBoundary", 2);
        model.component("comp1").physics("spf").feature("wall_fsi").selection()
            .named("fsi_boundary");

        // Initial values
        model.component("comp1").physics("spf").feature("init1")
            .set("u", "0");
        model.component("comp1").physics("spf").feature("init1")
            .set("p", "p_in");

        // ============================================================
        // 7. PHYSICS - SOLID MECHANICS
        // ============================================================

        // Create solid mechanics physics
        model.component("comp1").physics().create("solid", "SolidMechanics", "geom1");
        model.component("comp1").physics("solid").selection()
            .named("solid_domain");

        // Fixed constraint at shaft
        model.component("comp1").physics("solid").create("fix1", "Fixed", 2);
        model.component("comp1").physics("solid").feature("fix1").selection()
            .named("fixed_shaft");

        // Centrifugal force (body load)
        model.component("comp1").physics("solid").create("cent1", "CentrifugalForce", 3);
        model.component("comp1").physics("solid").feature("cent1").selection()
            .named("solid_domain");
        model.component("comp1").physics("solid").feature("cent1")
            .set("Omega", "omega");
        model.component("comp1").physics("solid").feature("cent1")
            .set("AxisDirection", new int[]{0, 0, 1});
        model.component("comp1").physics("solid").feature("cent1")
            .set("RotationPoint", new String[]{"0", "0", "0"});

        // Boundary load from fluid (handled by FSI coupling)

        // ============================================================
        // 8. MULTIPHYSICS - FLUID-STRUCTURE INTERACTION
        // ============================================================

        model.component("comp1").multiphysics().create("fsi1", "FluidStructureInteraction", 2);
        model.component("comp1").multiphysics("fsi1").selection()
            .named("fsi_boundary");
        model.component("comp1").multiphysics("fsi1")
            .set("SolidModel", "solid");
        model.component("comp1").multiphysics("fsi1")
            .set("FluidModel", "spf");

        // ============================================================
        // 9. MESH
        // ============================================================

        model.component("comp1").mesh().create("mesh1");

        // Fluid mesh
        model.component("comp1").mesh("mesh1").create("ftet_fluid", "FreeTet");
        model.component("comp1").mesh("mesh1").feature("ftet_fluid").selection()
            .geom("geom1", 3);
        model.component("comp1").mesh("mesh1").feature("ftet_fluid").selection()
            .named("fluid_domain");
        model.component("comp1").mesh("mesh1").feature("ftet_fluid")
            .create("size1", "Size");
        model.component("comp1").mesh("mesh1").feature("ftet_fluid").feature("size1")
            .set("hauto", 5);  // Medium mesh

        // Boundary layer at FSI interface
        model.component("comp1").mesh("mesh1").create("bl1", "BoundaryLayer");
        model.component("comp1").mesh("mesh1").feature("bl1").selection()
            .named("fsi_boundary");
        model.component("comp1").mesh("mesh1").feature("bl1")
            .set("numLayers", 5);
        model.component("comp1").mesh("mesh1").feature("bl1")
            .set("thickness", "0.005");

        // Solid mesh
        model.component("comp1").mesh("mesh1").create("ftet_solid", "FreeTet");
        model.component("comp1").mesh("mesh1").feature("ftet_solid").selection()
            .geom("geom1", 3);
        model.component("comp1").mesh("mesh1").feature("ftet_solid").selection()
            .named("solid_domain");
        model.component("comp1").mesh("mesh1").feature("ftet_solid")
            .create("size2", "Size");
        model.component("comp1").mesh("mesh1").feature("ftet_solid").feature("size2")
            .set("hauto", 6);  // Coarser for solid

        // Build mesh
        model.component("comp1").mesh("mesh1").run();

        // ============================================================
        // 10. STUDY 1: STATIONARY (INITIALIZE)
        // ============================================================

        // Stationary study to initialize flow field
        model.study().create("std1");
        model.study("std1").label("Stationary - Initialize Flow");
        model.study("std1").create("stat", "Stationary");

        // Activate only fluid physics for initialization
        model.study("std1").feature("stat").activate("spf", true);
        model.study("std1").feature("stat").activate("solid", false);

        // Solver for stationary study
        model.sol().create("sol1");
        model.sol("sol1").study("std1");
        model.sol("sol1").attach("std1");
        model.sol("sol1").feature().create("st1", "StudyStep");
        model.sol("sol1").feature().create("v1", "Variables");
        model.sol("sol1").feature().create("s1", "Stationary");

        model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");
        model.sol("sol1").feature("s1").feature("fc1").create("i1", "Iterative");
        model.sol("sol1").feature("s1").feature("fc1").feature("i1")
            .set("linsolver", "gmres");

        // Run stationary study
        model.sol("sol1").runAll();

        // ============================================================
        // 11. STUDY 2: PRESTRESSED (CENTRIFUGAL LOAD)
        // ============================================================

        // Stationary structural analysis with centrifugal load
        model.study().create("std2");
        model.study("std2").label("Prestressed - Centrifugal Load");
        model.study("std2").create("stat2", "Stationary");

        // Activate only solid physics
        model.study("std2").feature("stat2").activate("spf", false);
        model.study("std2").feature("stat2").activate("solid", true);

        // Solver for prestress study
        model.sol().create("sol2");
        model.sol("sol2").study("std2");
        model.sol("sol2").attach("std2");
        model.sol("sol2").feature().create("st2", "StudyStep");
        model.sol("sol2").feature().create("v2", "Variables");
        model.sol("sol2").feature().create("s2", "Stationary");

        model.sol("sol2").feature("s2").create("fc2", "FullyCoupled");
        model.sol("sol2").feature("s2").feature("fc2").create("d1", "Direct");

        // Run prestress study
        model.sol("sol2").runAll();

        // ============================================================
        // 12. STUDY 3: TIME-DEPENDENT FSI
        // ============================================================

        model.study().create("std3");
        model.study("std3").label("Transient FSI Analysis");
        model.study("std3").create("time", "Transient");

        // Time range: 0 to 0.5 seconds
        model.study("std3").feature("time")
            .set("tlist", "range(0,0.001,0.5)");

        // Activate both physics
        model.study("std3").feature("time").activate("spf", true);
        model.study("std3").feature("time").activate("solid", true);

        // Use results from previous studies as initial conditions
        model.study("std3").feature("time")
            .set("initstudyhide", "std1");  // Fluid initial condition
        model.study("std3").feature("time")
            .set("initsolhide", "sol1");
        model.study("std3").feature("time")
            .set("solnum", "auto");

        // Solver for FSI study
        model.sol().create("sol3");
        model.sol("sol3").study("std3");
        model.sol("sol3").attach("std3");
        model.sol("sol3").feature().create("st3", "StudyStep");
        model.sol("sol3").feature().create("v3", "Variables");
        model.sol("sol3").feature().create("t1", "Time");

        // Time-dependent solver settings
        model.sol("sol3").feature("t1")
            .set("tlist", "range(0,0.001,0.5)");
        model.sol("sol3").feature("t1")
            .set("rtol", "0.001");
        model.sol("sol3").feature("t1")
            .set("atol", "1e-6");
        model.sol("sol3").feature("t1")
            .set("timemethod", "genalpha");

        // Fully coupled solver with FSI
        model.sol("sol3").feature("t1").create("fc3", "FullyCoupled");
        model.sol("sol3").feature("t1").feature("fc3")
            .set("maxiter", "25");

        model.sol("sol3").feature("t1").feature("fc3").create("i3", "Iterative");
        model.sol("sol3").feature("t1").feature("fc3").feature("i3")
            .set("linsolver", "gmres");
        model.sol("sol3").feature("t1").feature("fc3").feature("i3")
            .set("rhob", "20");

        // Run FSI study
        model.sol("sol3").runAll();

        // ============================================================
        // 13. POST-PROCESSING
        // ============================================================

        // Create datasets
        model.result().dataset().create("dset1", "Solution");
        model.result().dataset("dset1").set("solution", "sol3");

        // Surface dataset for FSI boundary
        model.result().dataset().create("surf_fsi", "Surface");
        model.result().dataset("surf_fsi").set("data", "dset1");
        model.result().dataset("surf_fsi").selection().named("fsi_boundary");

        // ============================================================
        // 14. RESULTS EXTRACTION
        // ============================================================

        // Surface force on impeller
        model.result().numerical().create("force", "IntSurface");
        model.result().numerical("force").set("data", "dset1");
        model.result().numerical("force").selection().named("fsi_boundary");
        model.result().numerical("force")
            .set("expr", new String[]{"spf.Fp_x", "spf.Fp_y", "spf.Fp_z"});
        model.result().numerical("force")
            .set("descr", new String[]{"Force X [N]", "Force Y [N]", "Force Z [N]"});

        // Maximum displacement
        model.result().numerical().create("max_disp", "MaxVolume");
        model.result().numerical("max_disp").set("data", "dset1");
        model.result().numerical("max_disp").selection().named("solid_domain");
        model.result().numerical("max_disp")
            .set("expr", new String[]{"sqrt(u^2+v^2+w^2)"});
        model.result().numerical("max_disp")
            .set("descr", new String[]{"Max Displacement [m]"});

        // Maximum von Mises stress
        model.result().numerical().create("max_stress", "MaxVolume");
        model.result().numerical("max_stress").set("data", "dset1");
        model.result().numerical("max_stress").selection().named("solid_domain");
        model.result().numerical("max_stress")
            .set("expr", new String[]{"solid.mises"});
        model.result().numerical("max_stress")
            .set("descr", new String[]{"Max von Mises Stress [Pa]"});

        // ============================================================
        // 15. VISUALIZATION
        // ============================================================

        // Plot group: Pressure distribution
        model.result().create("pg_pressure", "PlotGroup3D");
        model.result("pg_pressure").label("Pressure Distribution");
        model.result("pg_pressure").set("data", "dset1");
        model.result("pg_pressure").create("surf1", "Surface");
        model.result("pg_pressure").feature("surf1")
            .set("expr", "p");
        model.result("pg_pressure").feature("surf1")
            .set("colortable", "Rainbow");

        // Plot group: Deformation
        model.result().create("pg_deform", "PlotGroup3D");
        model.result("pg_deform").label("Impeller Deformation");
        model.result("pg_deform").set("data", "dset1");
        model.result("pg_deform").create("surf2", "Surface");
        model.result("pg_deform").feature("surf2").selection().named("fsi_boundary");
        model.result("pg_deform").feature("surf2")
            .set("expr", "solid.mises");
        model.result("pg_deform").feature("surf2")
            .set("colortable", "ThermalWave");
        model.result("pg_deform").create("def1", "Deform");
        model.result("pg_deform").feature("def1")
            .set("expr", new String[]{"u", "v", "w"});
        model.result("pg_deform").feature("def1")
            .set("scale", "1000");  // Amplify displacement

        // ============================================================
        // 16. SAVE MODEL
        // ============================================================

        model.save("/path/to/pump_fsi_complete.mph");

        System.out.println("Pump FSI analysis model created successfully");
        return model;
    }

    public static void main(String[] args) {
        run();
    }
}
```

## MATLAB LiveLink Implementation

```matlab
%% Pump FSI Analysis using MATLAB LiveLink
% Initialize COMSOL with MATLAB
import com.comsol.model.*
import com.comsol.model.util.*

mphstart

%% Load or create model
% If running the Java code first, load the model:
model = mphload('/path/to/pump_fsi_complete.mph');

% Otherwise, run the Java method:
% model = PumpFSIAnalysis.run();

%% Parametric Study: Varying Rotational Speed
rpm_values = [1000, 1250, 1500, 1750, 2000];
results = struct();

for i = 1:length(rpm_values)
    fprintf('Running case %d: RPM = %d\n', i, rpm_values(i));

    % Set parameter
    model.param.set('rpm', sprintf('%d[rpm]', rpm_values(i)));

    % Run transient study
    model.study('std3').run();

    % Extract maximum displacement
    max_disp_data = mphmax(model, 'sqrt(u^2+v^2+w^2)', 'volume', ...
        'selection', 'solid_domain', 'dataset', 'dset1');

    % Extract maximum stress
    max_stress_data = mphmax(model, 'solid.mises', 'volume', ...
        'selection', 'solid_domain', 'dataset', 'dset1');

    % Extract force on impeller
    force_data = mphint2(model, {'spf.Fp_x', 'spf.Fp_y', 'spf.Fp_z'}, ...
        'surface', 'selection', 'fsi_boundary', 'dataset', 'dset1');

    % Store results
    results(i).rpm = rpm_values(i);
    results(i).max_displacement = max_disp_data;
    results(i).max_stress = max_stress_data;
    results(i).force_x = force_data(1);
    results(i).force_y = force_data(2);
    results(i).force_z = force_data(3);
    results(i).force_magnitude = sqrt(sum(force_data.^2));

    % Save case
    model.save(sprintf('/path/to/results/pump_fsi_rpm%d.mph', rpm_values(i)));
end

%% Plot Results
figure('Position', [100, 100, 1200, 800]);

% Displacement vs RPM
subplot(2,2,1);
plot([results.rpm], [results.max_displacement]*1e6, '-o', 'LineWidth', 2);
xlabel('Rotational Speed [rpm]');
ylabel('Max Displacement [μm]');
title('Maximum Blade Displacement');
grid on;

% Stress vs RPM
subplot(2,2,2);
plot([results.rpm], [results.max_stress]/1e6, '-o', 'LineWidth', 2);
xlabel('Rotational Speed [rpm]');
ylabel('Max von Mises Stress [MPa]');
title('Maximum Stress in Impeller');
grid on;

% Force magnitude vs RPM
subplot(2,2,3);
plot([results.rpm], [results.force_magnitude], '-o', 'LineWidth', 2);
xlabel('Rotational Speed [rpm]');
ylabel('Force Magnitude [N]');
title('Fluid Force on Impeller');
grid on;

% Safety factor
yield_strength = 290e6;  % Pa
safety_factors = yield_strength ./ [results.max_stress];
subplot(2,2,4);
plot([results.rpm], safety_factors, '-o', 'LineWidth', 2);
hold on;
yline(1.5, 'r--', 'Min SF = 1.5', 'LineWidth', 2);
xlabel('Rotational Speed [rpm]');
ylabel('Safety Factor');
title('Structural Safety Factor');
grid on;
legend('Calculated SF', 'Minimum Acceptable');

%% Save results
save('pump_fsi_parametric_results.mat', 'results');
saveas(gcf, 'pump_fsi_results.png');

%% FFT Analysis of Displacement Time History
% Extract time history at blade tip
blade_tip_point = [0.1, 0, 0.02];  % Coordinates of blade tip
time_data = mphglobal(model, 't', 'dataset', 'dset1');
disp_data = mphinterp(model, 'sqrt(u^2+v^2+w^2)', 'coord', blade_tip_point, ...
    'dataset', 'dset1');

% FFT
Fs = 1/(time_data(2) - time_data(1));  % Sampling frequency
L = length(disp_data);
Y = fft(disp_data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% Plot FFT
figure;
plot(f, P1*1e6, 'LineWidth', 2);
xlabel('Frequency [Hz]');
ylabel('Displacement Amplitude [μm]');
title('Frequency Spectrum of Blade Tip Vibration');
grid on;
xlim([0, 500]);  % Focus on relevant frequency range

% Identify dominant frequencies
[pks, locs] = findpeaks(P1*1e6, f, 'MinPeakHeight', 0.1, 'NPeaks', 5);
fprintf('\nDominant vibration frequencies:\n');
for i = 1:length(pks)
    fprintf('  %.1f Hz: %.2f μm\n', locs(i), pks(i));
end

%% Clean up
fprintf('\nAnalysis complete!\n');
```

## Expected Results

### Flow Field
- **Velocity Distribution**: Higher velocities at blade tips, recirculation at inlet
- **Pressure Rise**: Approximately 250 kPa from inlet to outlet
- **Turbulence**: High turbulence intensity in blade passages

### Structural Response
- **Static Deformation**: Primarily from centrifugal loading
  - Typical blade tip displacement: 50-200 μm
- **Dynamic Deformation**: Small additional vibration from fluid pulsations
  - Vibration amplitude: 1-10 μm
- **Stress Levels**:
  - Maximum at blade root and hub connection
  - Typical range: 50-150 MPa
  - Should be well below yield strength (290 MPa)

### FSI Effects
- **Fluid Loading**: Pressure and shear stress on blades
- **Blade Deformation**: Affects flow passage geometry
- **Coupled Frequencies**: Natural frequencies shifted by fluid added mass
- **Vibration Modes**: Blade bending, torsion, disk modes

### Performance Metrics
- **Safety Factor**: Should exceed 1.5-2.0
- **Resonance Check**: Operating frequency should avoid natural frequencies
- **Fatigue Life**: High-cycle fatigue assessment based on stress range

## Validation and Verification

### Mesh Convergence
1. Run with progressively finer mesh
2. Check convergence of maximum stress and displacement
3. Ensure less than 5% change with mesh refinement

### Time Step Convergence
1. Reduce time step by factor of 2
2. Verify results remain consistent
3. Balance accuracy and computational cost

### Comparison with Simplified Models
1. Compare rigid vs. flexible impeller
2. One-way vs. two-way FSI
3. Assess importance of coupling

### Experimental Validation (if available)
1. Compare pressure rise with pump curve
2. Validate vibration amplitudes with accelerometer data
3. Check natural frequencies with modal testing

## Troubleshooting

### Common Issues and Solutions

**Mesh Deformation Failure:**
- Reduce time step
- Increase mesh smoothing at FSI boundary
- Consider remeshing strategies for large deformations

**FSI Convergence Issues:**
- Adjust FSI relaxation factors
- Start from good initial condition (stationary solution)
- Use smaller time steps during startup transient

**Solver Memory Issues:**
- Use iterative solver instead of direct
- Reduce mesh size
- Exploit symmetry (if applicable)

**Unrealistic Results:**
- Check boundary conditions
- Verify material properties
- Review selection assignments
- Check units consistency

## Best Practices

1. **Initialization**: Always initialize with stationary flow solution
2. **Prestress**: Include centrifugal loading as prestress condition
3. **Mesh Quality**: Fine mesh at FSI interface, boundary layers for turbulence
4. **Time Stepping**: Start with small steps, can increase once steady periodic state reached
5. **Monitoring**: Track residuals, FSI iterations, and solution progress
6. **Validation**: Compare with analytical estimates and experimental data
7. **Documentation**: Save all parameter values and solver settings

## Additional Analysis Options

### Modal Analysis with Fluid Loading
Replace transient FSI with eigenfrequency analysis to find natural frequencies and mode shapes with fluid added mass effect.

### Frequency Response
Apply harmonic excitation and compute frequency response function for forced vibration analysis.

### Cavitation Analysis
Include cavitation model in CFD to assess vibration induced by vapor bubble collapse.

### Multi-Blade Passage
Model full 360° geometry or use periodic boundaries for multiple blades.

### Thermal Effects
Add heat transfer physics for high-temperature applications.
