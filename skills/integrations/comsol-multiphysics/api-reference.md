# COMSOL Java API Reference

Comprehensive reference for automating COMSOL Multiphysics simulations using the Java API.

## API Basics

### Importing Required Packages

```java
import com.comsol.model.*;
import com.comsol.model.util.*;

// For specific physics
import com.comsol.model.physics.*;
import com.comsol.model.mesh.*;
import com.comsol.model.geom.*;
```

### Model Creation and Management

```java
// Create new model
Model model = ModelUtil.create("ModelName");

// Load existing model
Model model = ModelUtil.load("/path/to/model.mph");

// Save model
model.save("/path/to/output.mph");

// Remove model from memory
ModelUtil.remove("ModelName");

// List all models in memory
String[] models = ModelUtil.tags();

// Clear all models
ModelUtil.clear();

// Show model in GUI (if available)
ModelUtil.showProgress(true);
model.show();
```

### Model Parameters

```java
// Create global parameters
model.param().set("inlet_velocity", "5[m/s]");
model.param().set("outlet_pressure", "101325[Pa]");
model.param().set("rotation_speed", "1500[rpm]");
model.param().set("temperature", "293.15[K]");

// Get parameter value
String velocity = model.param().get("inlet_velocity");

// Remove parameter
model.param().remove("parameter_name");

// Create parameter group
model.param().group().create("par1");
model.param("par1").set("diameter", "0.1[m]");
model.param("par1").set("length", "1[m]");

// Parameter expressions (using other parameters)
model.param().set("area", "pi*diameter^2/4");
model.param().set("reynolds", "rho*inlet_velocity*diameter/mu");
```

## Model Building Commands

### Component Management

```java
// Create component (3D by default)
model.component().create("comp1", true);

// Create 2D component
model.component().create("comp1", 2);

// Create 2D axisymmetric component
model.component().create("comp1", "2Daxi");

// Create 1D component
model.component().create("comp1", 1);

// Get component
ModelEntity comp = model.component("comp1");
```

### Geometry Operations

```java
// Create geometry sequence
model.component("comp1").geom().create("geom1", 3);

// Import CAD file
model.component("comp1").geom("geom1").create("imp1", "Import");
model.component("comp1").geom("geom1").feature("imp1")
    .set("filename", "/path/to/geometry.step");
model.component("comp1").geom("geom1").feature("imp1")
    .set("type", "cad");

// Create block
model.component("comp1").geom("geom1").create("blk1", "Block");
model.component("comp1").geom("geom1").feature("blk1")
    .set("size", new String[]{"1", "0.5", "0.3"});
model.component("comp1").geom("geom1").feature("blk1")
    .set("pos", new String[]{"0", "0", "0"});

// Create cylinder
model.component("comp1").geom("geom1").create("cyl1", "Cylinder");
model.component("comp1").geom("geom1").feature("cyl1")
    .set("r", "0.05");
model.component("comp1").geom("geom1").feature("cyl1")
    .set("h", "0.2");

// Create sphere
model.component("comp1").geom("geom1").create("sph1", "Sphere");
model.component("comp1").geom("geom1").feature("sph1")
    .set("r", "0.1");

// Boolean operations
model.component("comp1").geom("geom1").create("uni1", "Union");
model.component("comp1").geom("geom1").feature("uni1")
    .selection("input").set("blk1", "cyl1");

model.component("comp1").geom("geom1").create("dif1", "Difference");
model.component("comp1").geom("geom1").feature("dif1")
    .selection("input").set("blk1");
model.component("comp1").geom("geom1").feature("dif1")
    .selection("input2").set("cyl1");

// Run geometry
model.component("comp1").geom("geom1").run();

// Run specific feature
model.component("comp1").geom("geom1").run("blk1");
```

### Selections

```java
// Create explicit selection (domains)
model.component("comp1").geom("geom1").create("sel1", "ExplicitSelection");
model.component("comp1").geom("geom1").feature("sel1")
    .selection("selection").init(3);
model.component("comp1").geom("geom1").feature("sel1")
    .selection("selection").set("blk1", 1);

// Create box selection
model.component("comp1").geom("geom1").create("boxsel1", "BoxSelection");
model.component("comp1").geom("geom1").feature("boxsel1")
    .set("entitydim", 2);  // Select boundaries
model.component("comp1").geom("geom1").feature("boxsel1")
    .set("xmin", "-0.1");
model.component("comp1").geom("geom1").feature("boxsel1")
    .set("xmax", "0.1");

// Create cumulative selection
model.component("comp1").selection().create("fluid_domain", "Explicit");
model.component("comp1").selection("fluid_domain").geom("geom1", 3);
model.component("comp1").selection("fluid_domain").set(1, 2, 3);

// Boundary selection
model.component("comp1").selection().create("inlet", "Explicit");
model.component("comp1").selection("inlet").geom("geom1", 2);
model.component("comp1").selection("inlet").set(5);
```

### Physics Interfaces

#### CFD Module - Laminar and Turbulent Flow

```java
// Single-phase flow
model.component("comp1").physics().create("spf", "LaminarFlow", "geom1");

// Set domain selection
model.component("comp1").physics("spf").selection()
    .named("geom1_fluid_domain");

// Turbulent flow with k-epsilon
model.component("comp1").physics().create("spf", "TurbulentFlowkEps", "geom1");

// Turbulent flow with k-omega
model.component("comp1").physics().create("spf", "TurbulentFlowkOmega", "geom1");

// Enable compressibility
model.component("comp1").physics("spf").feature("fp1")
    .set("CompressibilityType", "Compressible");

// Rotating reference frame
model.component("comp1").physics("spf").create("rot1", "RotatingDomain", 3);
model.component("comp1").physics("spf").feature("rot1").selection()
    .named("geom1_rotating_domain");
model.component("comp1").physics("spf").feature("rot1")
    .set("omega", "rotation_speed");
model.component("comp1").physics("spf").feature("rot1")
    .set("axis", new String[]{"0", "0", "1"});

// Inlet boundary condition
model.component("comp1").physics("spf").create("inlet1", "InletBoundary", 2);
model.component("comp1").physics("spf").feature("inlet1").selection()
    .named("geom1_inlet");
model.component("comp1").physics("spf").feature("inlet1")
    .set("BoundaryCondition", "Velocity");
model.component("comp1").physics("spf").feature("inlet1")
    .set("U0", "inlet_velocity");

// Outlet boundary condition
model.component("comp1").physics("spf").create("outlet1", "OutletBoundary", 2);
model.component("comp1").physics("spf").feature("outlet1").selection()
    .named("geom1_outlet");
model.component("comp1").physics("spf").feature("outlet1")
    .set("p0", "outlet_pressure");

// Wall boundary condition
model.component("comp1").physics("spf").create("wall1", "WallBoundary", 2);
model.component("comp1").physics("spf").feature("wall1").selection()
    .named("geom1_wall");
model.component("comp1").physics("spf").feature("wall1")
    .set("BoundaryCondition", "NoSlip");

// Symmetry boundary
model.component("comp1").physics("spf").create("sym1", "SymmetryBoundary", 2);
model.component("comp1").physics("spf").feature("sym1").selection()
    .named("geom1_symmetry");

// Initial values
model.component("comp1").physics("spf").feature("init1")
    .set("u", "inlet_velocity");
model.component("comp1").physics("spf").feature("init1")
    .set("p", "outlet_pressure");
```

#### Structural Mechanics Module

```java
// Solid mechanics
model.component("comp1").physics().create("solid", "SolidMechanics", "geom1");

// Set domain selection
model.component("comp1").physics("solid").selection()
    .named("geom1_solid_domain");

// Material assignment
model.component("comp1").physics("solid").feature("lemm1")
    .set("SolidModel", "Linear elastic");
model.component("comp1").physics("solid").feature("lemm1").selection()
    .named("geom1_solid_domain");

// Fixed constraint
model.component("comp1").physics("solid").create("fix1", "Fixed", 2);
model.component("comp1").physics("solid").feature("fix1").selection()
    .named("geom1_fixed_boundary");

// Prescribed displacement
model.component("comp1").physics("solid").create("disp1", "Displacement", 2);
model.component("comp1").physics("solid").feature("disp1").selection()
    .named("geom1_prescribed_boundary");
model.component("comp1").physics("solid").feature("disp1")
    .set("Direction", new int[]{1, 0, 0});
model.component("comp1").physics("solid").feature("disp1")
    .set("U0", "0.001[m]");

// Boundary load (pressure)
model.component("comp1").physics("solid").create("bndl1", "BoundaryLoad", 2);
model.component("comp1").physics("solid").feature("bndl1").selection()
    .named("geom1_loaded_boundary");
model.component("comp1").physics("solid").feature("bndl1")
    .set("LoadType", "Pressure");
model.component("comp1").physics("solid").feature("bndl1")
    .set("Pressure", "1e6[Pa]");

// Body load (gravity)
model.component("comp1").physics("solid").create("bl1", "BodyLoad", 3);
model.component("comp1").physics("solid").feature("bl1").selection()
    .named("geom1_solid_domain");
model.component("comp1").physics("solid").feature("bl1")
    .set("FperVol", new String[]{"0", "0", "-rho_material*g_const"});

// Centrifugal load
model.component("comp1").physics("solid").create("cent1", "CentrifugalForce", 3);
model.component("comp1").physics("solid").feature("cent1").selection()
    .named("geom1_rotating_domain");
model.component("comp1").physics("solid").feature("cent1")
    .set("Omega", "rotation_speed");
model.component("comp1").physics("solid").feature("cent1")
    .set("AxisDirection", new int[]{0, 0, 1});

// Initial stress (prestress)
model.component("comp1").physics("solid").create("is1", "InitialStress", 3);
model.component("comp1").physics("solid").feature("is1").selection()
    .named("geom1_prestressed_domain");
model.component("comp1").physics("solid").feature("is1")
    .set("Sil", "sol2", "comp1.Sl");
```

### Materials

```java
// Create material
model.component("comp1").material().create("mat1", "Common");
model.component("comp1").material("mat1").label("Water");

// Assign to domain
model.component("comp1").material("mat1").selection()
    .named("geom1_fluid_domain");

// Set properties
model.component("comp1").material("mat1").propertyGroup("def")
    .set("density", "998.2[kg/m^3]");
model.component("comp1").material("mat1").propertyGroup("def")
    .set("dynamicviscosity", "0.001003[Pa*s]");

// Create material for solid
model.component("comp1").material().create("mat2", "Common");
model.component("comp1").material("mat2").label("Steel");
model.component("comp1").material("mat2").selection()
    .named("geom1_solid_domain");

model.component("comp1").material("mat2").propertyGroup("def")
    .set("density", "7850[kg/m^3]");
model.component("comp1").material("mat2").propertyGroup("def")
    .set("youngsmodulus", "200e9[Pa]");
model.component("comp1").material("mat2").propertyGroup("def")
    .set("poissonsratio", "0.33");

// Use material from built-in library
model.component("comp1").material().create("mat3", "Common");
model.component("comp1").material("mat3").label("Air");
model.component("comp1").material("mat3").materialModel("def")
    .set("from", "mat", "Common", "Air");
```

### Multiphysics Couplings

```java
// Fluid-Structure Interaction
model.component("comp1").multiphysics().create("fsi1", "FluidStructureInteraction", 2);
model.component("comp1").multiphysics("fsi1").selection()
    .named("geom1_fsi_boundary");

// Pair physics interfaces
model.component("comp1").multiphysics("fsi1")
    .set("SolidModel", "solid");
model.component("comp1").multiphysics("fsi1")
    .set("FluidModel", "spf");

// Thermal-structure coupling
model.component("comp1").multiphysics().create("te1", "ThermalExpansion", 3);
model.component("comp1").multiphysics("te1").selection()
    .named("geom1_solid_domain");
model.component("comp1").multiphysics("te1")
    .set("SolidModel", "solid");
model.component("comp1").multiphysics("te1")
    .set("TemperatureModel", "ht");

// Acoustic-structure interaction
model.component("comp1").multiphysics().create("actd1", "AcousticStructureBoundary", 2);
model.component("comp1").multiphysics("actd1").selection()
    .named("geom1_acoustic_boundary");
```

### Mesh Generation

```java
// Create mesh
model.component("comp1").mesh().create("mesh1");

// Automatic mesh
model.component("comp1").mesh("mesh1").automatic(true);

// Or custom mesh sequence
model.component("comp1").mesh("mesh1").automatic(false);

// Free tetrahedral
model.component("comp1").mesh("mesh1").create("ftet1", "FreeTet");
model.component("comp1").mesh("mesh1").feature("ftet1")
    .set("hauto", 3);  // Mesh size: 1-9 (fine to coarse)

// Size parameters
model.component("comp1").mesh("mesh1").create("size1", "Size");
model.component("comp1").mesh("mesh1").feature("size1")
    .set("hauto", 1);
model.component("comp1").mesh("mesh1").feature("size1")
    .set("hmax", "0.05");
model.component("comp1").mesh("mesh1").feature("size1")
    .set("hmin", "0.001");
model.component("comp1").mesh("mesh1").feature("size1")
    .set("hgrad", "1.3");  // Growth rate
model.component("comp1").mesh("mesh1").feature("size1")
    .set("hcurve", "0.3");  // Curvature factor

// Boundary layer mesh
model.component("comp1").mesh("mesh1").create("bl1", "BoundaryLayer");
model.component("comp1").mesh("mesh1").feature("bl1").selection()
    .named("geom1_wall");
model.component("comp1").mesh("mesh1").feature("bl1")
    .set("numLayers", 5);
model.component("comp1").mesh("mesh1").feature("bl1")
    .set("thickness", "0.01");
model.component("comp1").mesh("mesh1").feature("bl1")
    .set("thicknessAdj", "soft");  // soft or hard

// Swept mesh
model.component("comp1").mesh("mesh1").create("swe1", "Sweep");
model.component("comp1").mesh("mesh1").feature("swe1")
    .set("numLayers", 10);

// Run mesh
model.component("comp1").mesh("mesh1").run();

// Mesh statistics
model.component("comp1").mesh("mesh1").stat();
```

## Solver Settings

### Study Configuration

```java
// Create study
model.study().create("std1");

// Stationary study
model.study("std1").create("stat", "Stationary");

// Time-dependent study
model.study("std1").create("time", "Transient");
model.study("std1").feature("time")
    .set("tlist", "range(0,0.01,1)");  // 0 to 1 with 0.01 step

// Eigenfrequency study
model.study("std1").create("eig", "Eigenfrequency");
model.study("std1").feature("eig")
    .set("neigs", 10);  // Number of eigenfrequencies
model.study("std1").feature("eig")
    .set("shift", "0");  // Search near this frequency

// Frequency domain study
model.study("std1").create("freq", "Frequency");
model.study("std1").feature("freq")
    .set("plist", "range(0,10,1000)");  // 0 to 1000 Hz

// Parametric sweep
model.study("std1").create("param", "Parametric");
model.study("std1").feature("param")
    .set("pname", "inlet_velocity");
model.study("std1").feature("param")
    .set("plistarr", new String[]{"2", "4", "6", "8", "10"});
model.study("std1").feature("param")
    .set("punit", "m/s");

// Activate physics and variables
model.study("std1").feature("time")
    .activate("spf", true);
model.study("std1").feature("time")
    .activate("solid", true);
```

### Solver Configuration

```java
// Create solver
model.sol().create("sol1");
model.sol("sol1").study("std1");

// Attach to study
model.sol("sol1").attach("std1");

// Study step
model.sol("sol1").feature().create("st1", "StudyStep");
model.sol("sol1").feature("st1").set("study", "std1");
model.sol("sol1").feature("st1").set("studystep", "time");

// Variables
model.sol("sol1").feature().create("v1", "Variables");
model.sol("sol1").feature("v1").set("control", "time");

// Time-dependent solver
model.sol("sol1").feature().create("t1", "Time");
model.sol("sol1").feature("t1").set("tlist", "range(0,0.01,1)");
model.sol("sol1").feature("t1").set("rtol", "0.001");  // Relative tolerance
model.sol("sol1").feature("t1").set("atol", "1e-6");   // Absolute tolerance
model.sol("sol1").feature("t1").set("timemethod", "genalpha");  // Time stepping

// Fully coupled solver
model.sol("sol1").feature("t1").create("fc1", "FullyCoupled");
model.sol("sol1").feature("t1").feature("fc1")
    .set("termonres", "auto");  // Termination criterion

// Iterative solver (for large models)
model.sol("sol1").feature("t1").feature("fc1").create("i1", "Iterative");
model.sol("sol1").feature("t1").feature("fc1").feature("i1")
    .set("linsolver", "gmres");  // GMRES, BiCGStab, etc.
model.sol("sol1").feature("t1").feature("fc1").feature("i1")
    .set("rhob", "20");  // Maximum iterations

// Direct solver (MUMPS, PARDISO)
model.sol("sol1").feature("t1").feature("fc1").create("d1", "Direct");
model.sol("sol1").feature("t1").feature("fc1").feature("d1")
    .set("linsolver", "mumps");

// Segregated solver
model.sol("sol1").feature("t1").create("se1", "Segregated");
model.sol("sol1").feature("t1").feature("se1")
    .set("segvar", new String[]{"u", "p", "T"});  // Variables to segregate

// Run solver
model.sol("sol1").runAll();

// Or run specific study
model.study("std1").run();
```

### Advanced Solver Settings

```java
// Adaptive time stepping
model.sol("sol1").feature("t1")
    .set("initialstepgenalpha", "0.001");  // Initial step
model.sol("sol1").feature("t1")
    .set("maxstepgenalpha", "0.1");  // Maximum step
model.sol("sol1").feature("t1")
    .set("minstepgenalpha", "1e-6");  // Minimum step
model.sol("sol1").feature("t1")
    .set("estrat", "exclude");  // Error estimation strategy

// FSI coupling settings
model.sol("sol1").feature("t1").create("fsi", "FluidStructureInteraction");
model.sol("sol1").feature("t1").feature("fsi")
    .set("maxiter", "25");  // Max FSI iterations
model.sol("sol1").feature("t1").feature("fsi")
    .set("reltol", "1e-3");  // FSI convergence tolerance

// Nonlinear solver settings
model.sol("sol1").feature("t1").feature("fc1")
    .set("dtech", "auto");  // Damping technique
model.sol("sol1").feature("t1").feature("fc1")
    .set("maxiter", "25");  // Max nonlinear iterations

// Solver log
model.sol("sol1").feature("t1")
    .set("solfile", "results/solution.sol");
```

## Results Extraction

### Post-Processing Setup

```java
// Create dataset
model.result().dataset().create("dset1", "Solution");
model.result().dataset("dset1").set("solution", "sol1");

// Create surface dataset
model.result().dataset().create("surf1", "Surface");
model.result().dataset("surf1").set("data", "dset1");
model.result().dataset("surf1").selection().named("geom1_surface");

// Create cut plane
model.result().dataset().create("cpl1", "CutPlane");
model.result().dataset("cpl1").set("data", "dset1");
model.result().dataset("cpl1").set("quickplane", "xy");  // xy, xz, yz
model.result().dataset("cpl1").set("quickz", "0");

// Create particle tracing dataset
model.result().dataset().create("ptcl1", "Particle");
model.result().dataset("ptcl1").set("data", "dset1");
model.result().dataset("ptcl1").set("expr", new String[]{"u", "v", "w"});
```

### Numerical Evaluation

```java
// Global evaluation
model.result().numerical().create("gev1", "EvalGlobal");
model.result().numerical("gev1").set("data", "dset1");
model.result().numerical("gev1").set("expr", new String[]{"inlet_velocity", "outlet_pressure"});
double[][] results = model.result().numerical("gev1").getData();

// Surface integration
model.result().numerical().create("int1", "IntSurface");
model.result().numerical("int1").set("data", "dset1");
model.result().numerical("int1").selection().named("geom1_surface");
model.result().numerical("int1").set("expr", new String[]{"spf.Fp_x", "spf.Fp_y", "spf.Fp_z"});
model.result().numerical("int1").set("descr", new String[]{"Force X", "Force Y", "Force Z"});
double[][] forces = model.result().numerical("int1").getData();

// Volume integration
model.result().numerical().create("intv1", "IntVolume");
model.result().numerical("intv1").set("data", "dset1");
model.result().numerical("intv1").selection().named("geom1_domain");
model.result().numerical("intv1").set("expr", new String[]{"0.5*rho*(u^2+v^2+w^2)"});
model.result().numerical("intv1").set("descr", new String[]{"Kinetic Energy"});
double[][] energy = model.result().numerical("intv1").getData();

// Average value
model.result().numerical().create("av1", "AvSurface");
model.result().numerical("av1").set("data", "dset1");
model.result().numerical("av1").selection().named("geom1_outlet");
model.result().numerical("av1").set("expr", new String[]{"p", "spf.U"});
double[][] avg = model.result().numerical("av1").getData();

// Maximum/minimum
model.result().numerical().create("max1", "MaxVolume");
model.result().numerical("max1").set("data", "dset1");
model.result().numerical("max1").set("expr", new String[]{"sqrt(u^2+v^2+w^2)"});
double[][] maxVel = model.result().numerical("max1").getData();
```

### Data Export

```java
// Export to text file
model.result().export().create("data1", "Data");
model.result().export("data1").set("data", "dset1");
model.result().export("data1").set("expr", new String[]{"u", "v", "w", "p"});
model.result().export("data1").set("filename", "/path/to/results.txt");
model.result().export("data1").run();

// Export image
model.result().export().create("img1", "Image");
model.result().export("img1").set("sourceobject", "pg1");  // Plot group
model.result().export("img1").set("filename", "/path/to/plot.png");
model.result().export("img1").set("size", "manualweb");
model.result().export("img1").set("unit", "px");
model.result().export("img1").set("height", "600");
model.result().export("img1").set("width", "800");
model.result().export("img1").run();

// Export animation
model.result().export().create("anim1", "Animation");
model.result().export("anim1").set("sourceobject", "pg1");
model.result().export("anim1").set("filename", "/path/to/animation.gif");
model.result().export("anim1").set("fps", "10");
model.result().export("anim1").run();
```

### Visualization

```java
// Create plot group
model.result().create("pg1", "PlotGroup3D");
model.result("pg1").set("data", "dset1");

// Surface plot
model.result("pg1").create("surf1", "Surface");
model.result("pg1").feature("surf1").set("expr", "p");
model.result("pg1").feature("surf1").set("colortable", "Rainbow");
model.result("pg1").feature("surf1").set("rangecoloractive", true);
model.result("pg1").feature("surf1").set("rangecolormin", "0");
model.result("pg1").feature("surf1").set("rangecolormax", "1e5");

// Contour plot
model.result("pg1").create("con1", "Contour");
model.result("pg1").feature("con1").set("expr", "sqrt(u^2+v^2+w^2)");
model.result("pg1").feature("con1").set("number", "20");

// Streamline plot
model.result("pg1").create("strl1", "Streamline");
model.result("pg1").feature("strl1").set("expr", new String[]{"u", "v", "w"});
model.result("pg1").feature("strl1").set("number", "50");

// Vector plot
model.result("pg1").create("arws1", "ArrowSurface");
model.result("pg1").feature("arws1").set("expr", new String[]{"u", "v", "w"});
model.result("pg1").feature("arws1").set("scale", "0.01");

// Deformation plot
model.result("pg1").create("def1", "Deform");
model.result("pg1").feature("def1").set("expr", new String[]{"u", "v", "w"});
model.result("pg1").feature("def1").set("scale", "1000");  // Amplification factor

// Run plot
model.result("pg1").run();
```

## Documentation Links

### Official COMSOL Documentation

**Main Documentation:**
- COMSOL Multiphysics Reference Manual
- COMSOL Multiphysics User's Guide
- COMSOL API Reference Manual

**Module-Specific:**
- CFD Module User's Guide
  - https://doc.comsol.com/6.0/docserver/#!/com.comsol.help.cfd
- Structural Mechanics Module User's Guide
  - https://doc.comsol.com/6.0/docserver/#!/com.comsol.help.sme
- Acoustics Module User's Guide
  - https://doc.comsol.com/6.0/docserver/#!/com.comsol.help.aco

**Programming:**
- COMSOL Multiphysics Programming Reference Manual
  - Complete Java API documentation
  - https://doc.comsol.com/6.0/docserver/#!/com.comsol.help.comsol/Programming
- LiveLink for MATLAB User's Guide
  - https://doc.comsol.com/6.0/docserver/#!/com.comsol.help.llmatlab

**Application Library:**
- Searchable database of example models
- https://www.comsol.com/models
- Filter by module and application area

### Online Resources

**COMSOL Website:**
- Main site: https://www.comsol.com
- Support Center: https://www.comsol.com/support
- Blog: https://www.comsol.com/blogs
- Video Gallery: https://www.comsol.com/video

**Community:**
- Discussion Forum: https://www.comsol.com/forum
- Model Exchange: https://www.comsol.com/community/exchange

**Training:**
- Learning Center: https://www.comsol.com/learning-center
- Upcoming Events: https://www.comsol.com/events
- Webinar Archive: https://www.comsol.com/video/webinar

### Support Channels

**Technical Support:**
- Email: support@comsol.com
- Phone: Regional offices (see COMSOL website)
- Online case submission via Support Center

**Application Engineering:**
- Consulting on complex multiphysics problems
- Model development assistance
- Custom training sessions

**Regional Offices:**
- North America: Burlington, MA (headquarters)
- Europe: Stockholm, Sweden (headquarters)
- Asia: Multiple locations

### Additional Resources

**Textbooks and References:**
- "Multiphysics Modeling Using COMSOL" (various authors)
- "Introduction to COMSOL Multiphysics" (Tabatabaian)
- Academic papers using COMSOL (Google Scholar)

**Third-Party Tutorials:**
- YouTube channels with COMSOL tutorials
- University course materials
- Blog posts and technical articles

**Version-Specific Documentation:**
Note: API syntax may vary between COMSOL versions. Always refer to documentation matching your installed version.

Current latest version: COMSOL 6.2 (as of 2024)
- Check version: Help > About COMSOL Multiphysics
- Documentation in COMSOL: Help > Documentation
- Built-in API help: Methods Browser in Model Builder

### Quick Reference

**Common Property Names:**
- Geometry: `pos`, `size`, `r` (radius), `h` (height)
- Physics: `U0` (velocity), `p0` (pressure), `T0` (temperature)
- Solver: `rtol` (relative tolerance), `atol` (absolute tolerance)
- Mesh: `hauto` (size), `hmax`, `hmin`, `hgrad` (growth rate)

**Physics Tags:**
- `spf`: Single-Phase Flow
- `solid`: Solid Mechanics
- `ht`: Heat Transfer
- `acpr`: Pressure Acoustics
- `fsi1`: Fluid-Structure Interaction

**Study Tags:**
- `stat`: Stationary
- `time`: Time Dependent
- `eig`: Eigenfrequency
- `freq`: Frequency Domain

**Common Methods:**
- `.create()`: Create new object
- `.set()`: Set property
- `.get()`: Get property
- `.selection()`: Access selection object
- `.run()`: Execute operation
- `.attach()`: Attach to parent object
