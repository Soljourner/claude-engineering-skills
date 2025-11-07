---
name: openfoam-cfd
description: "Set up turbulent flow simulations in OpenFOAM with automated case generation"
category: integrations
domain: fluids
complexity: advanced
dependencies: []
---

# OpenFOAM CFD Integration

## Overview

OpenFOAM (Open Field Operation and Manipulation) is a free, open-source CFD software package that provides a comprehensive suite of tools for solving complex fluid dynamics problems. It uses the finite volume method to discretize and solve partial differential equations describing fluid flow, heat transfer, turbulence, and chemical reactions.

### Key Features

- **Open Source**: Licensed under GPLv3, completely free to use and modify
- **Extensive Physics**: Handles incompressible/compressible flows, multiphase, combustion, heat transfer
- **Parallel Computing**: Built-in support for MPI-based parallel processing
- **Customizable**: Written in C++, highly extensible for custom physics and solvers
- **Industry Standard**: Used in academia and industry for research and commercial applications

### Architecture

OpenFOAM uses a case-based structure where each simulation is a self-contained directory with:
- Initial and boundary conditions
- Physical properties
- Numerical schemes
- Solver settings
- Mesh data

## Installation

### Ubuntu/Debian Installation

#### Option 1: Official Ubuntu Repository (Latest)

```bash
# Add OpenFOAM repository
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get update

# Install OpenFOAM v11 (or latest version)
sudo apt-get install openfoam11

# Source the OpenFOAM environment
echo "source /opt/openfoam11/etc/bashrc" >> ~/.bashrc
source ~/.bashrc
```

#### Option 2: ESI-OpenCFD Version

```bash
# Add repository
sudo sh -c "wget -O - https://dl.openfoam.com/add-debian-repo.sh | bash"

# Install OpenFOAM v2312 (or latest version)
sudo apt-get install openfoam2312-default

# Source environment
echo "source /usr/lib/openfoam/openfoam2312/etc/bashrc" >> ~/.bashrc
source ~/.bashrc
```

#### Option 3: Build from Source

```bash
# Install dependencies
sudo apt-get install build-essential cmake git ca-certificates \
    flex libfl-dev bison zlib1g-dev libboost-system-dev \
    libboost-thread-dev libopenmpi-dev openmpi-bin \
    gnuplot libreadline-dev libncurses-dev libxt-dev \
    qt5-default libqt5x11extras5-dev libqt5help5 qttools5-dev \
    qtxmlpatterns5-dev-tools libqt5opengl5-dev freeglut3-dev \
    libscotch-dev libcgal-dev

# Clone OpenFOAM
git clone https://github.com/OpenFOAM/OpenFOAM-11.git
cd OpenFOAM-11

# Source and compile
source etc/bashrc
./Allwmake -j
```

#### Verify Installation

```bash
# Check OpenFOAM environment
which simpleFoam
foamInstallationTest

# Run tutorial
mkdir -p $FOAM_RUN
cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily $FOAM_RUN/
cd $FOAM_RUN/pitzDaily
blockMesh
simpleFoam
```

### Install ParaView for Post-Processing

```bash
# Install ParaView
sudo apt-get install paraview

# Or use OpenFOAM's reader
paraFoam
```

## Case Directory Structure

Every OpenFOAM case follows a standardized directory structure:

```
caseDirectory/
├── 0/                      # Initial and boundary conditions
│   ├── U                   # Velocity field
│   ├── p                   # Pressure field
│   ├── T                   # Temperature field (if applicable)
│   ├── k                   # Turbulent kinetic energy
│   ├── epsilon             # Turbulent dissipation rate
│   ├── omega               # Specific dissipation rate
│   └── nut                 # Turbulent viscosity
├── constant/               # Physical properties and mesh
│   ├── polyMesh/          # Mesh data
│   │   ├── points         # Vertex coordinates
│   │   ├── faces          # Face definitions
│   │   ├── owner          # Cell ownership
│   │   ├── neighbour      # Cell neighbors
│   │   └── boundary       # Boundary patches
│   ├── transportProperties # Fluid properties
│   ├── turbulenceProperties # Turbulence model selection
│   ├── g                  # Gravity vector
│   └── regionProperties   # Multi-region properties
├── system/                 # Numerical settings
│   ├── controlDict        # Simulation control
│   ├── fvSchemes          # Discretization schemes
│   ├── fvSolution         # Solver settings
│   ├── blockMeshDict      # Mesh generation
│   ├── snappyHexMeshDict  # Advanced meshing
│   ├── decomposeParDict   # Parallel decomposition
│   └── sampleDict         # Sampling/probing
└── [time directories]/     # Results at each time step
    ├── U
    ├── p
    └── ...
```

## Key Dictionaries

### controlDict (system/controlDict)

Controls simulation execution and output:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      controlDict;
}

application     simpleFoam;      // Solver to use

startFrom       startTime;        // latestTime, firstTime
startTime       0;
stopAt          endTime;
endTime         1000;
deltaT          1;                // Time step

writeControl    timeStep;         // runTime, adjustableRunTime, cpuTime
writeInterval   100;              // Write every N steps
purgeWrite      2;                // Keep only last N time directories

writeFormat     ascii;            // binary
writePrecision  6;
writeCompression off;             // on for binary

timeFormat      general;
timePrecision   6;

runTimeModifiable true;           // Re-read dictionaries during run

// Function objects for monitoring and post-processing
functions
{
    forceCoeffs
    {
        type            forceCoeffs;
        libs            ("libforces.so");
        writeControl    timeStep;
        writeInterval   1;

        patches         (wall);
        rho             rhoInf;
        rhoInf          1.225;
        CofR            (0 0 0);
        liftDir         (0 1 0);
        dragDir         (1 0 0);
        pitchAxis       (0 0 1);
        magUInf         10;
        lRef            1;
        Aref            1;
    }

    probes
    {
        type            probes;
        libs            ("libsampling.so");
        writeControl    timeStep;
        writeInterval   1;

        fields          (p U);
        probeLocations
        (
            (0.5 0.5 0.5)
            (1.0 0.5 0.5)
        );
    }
}
```

### fvSchemes (system/fvSchemes)

Defines discretization schemes:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}

// Time derivatives
ddtSchemes
{
    default         steadyState;    // Euler, backward, CrankNicolson
}

// Gradient schemes
gradSchemes
{
    default         Gauss linear;   // cellLimited Gauss linear 1
    grad(p)         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}

// Divergence schemes (convection)
divSchemes
{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}

// Laplacian schemes (diffusion)
laplacianSchemes
{
    default         Gauss linear corrected;
}

// Interpolation schemes
interpolationSchemes
{
    default         linear;
}

// Surface normal gradient schemes
snGradSchemes
{
    default         corrected;
}

// Wall distance calculation
wallDist
{
    method          meshWave;
}
```

### fvSolution (system/fvSolution)

Solver algorithms and tolerances:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}

solvers
{
    p
    {
        solver          GAMG;           // Geometric-Algebraic Multi-Grid
        tolerance       1e-6;
        relTol          0.01;
        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     2;
        cacheAgglomeration on;
        agglomerator    faceAreaPair;
        nCellsInCoarsestLevel 10;
        mergeLevels     1;
    }

    U
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-7;
        relTol          0.1;
        nSweeps         1;
    }

    "(k|epsilon|omega)"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0.1;
        nSweeps         1;
    }
}

// SIMPLE algorithm for steady-state
SIMPLE
{
    nNonOrthogonalCorrectors 0;
    consistent      yes;            // SIMPLEC algorithm

    residualControl
    {
        p               1e-5;
        U               1e-5;
        "(k|epsilon|omega)" 1e-5;
    }
}

// PIMPLE algorithm for transient
PIMPLE
{
    nOuterCorrectors    50;
    nCorrectors         2;
    nNonOrthogonalCorrectors 1;

    residualControl
    {
        U               1e-5;
        p               1e-5;
    }
}

relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        k               0.7;
        epsilon         0.7;
        omega           0.7;
    }
}
```

### transportProperties (constant/transportProperties)

Physical properties of fluids:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}

// Single phase
transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] 1e-05;  // Kinematic viscosity (m²/s)

// For multiphase (interFoam)
phases (water air);

water
{
    transportModel  Newtonian;
    nu              [0 2 -1 0 0 0 0] 1e-06;
    rho             [1 -3 0 0 0 0 0] 1000;
}

air
{
    transportModel  Newtonian;
    nu              [0 2 -1 0 0 0 0] 1.48e-05;
    rho             [1 -3 0 0 0 0 0] 1;
}

sigma           [1 0 -2 0 0 0 0] 0.07;  // Surface tension
```

### turbulenceProperties (constant/turbulenceProperties)

Turbulence model selection:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      turbulenceProperties;
}

simulationType  RAS;  // laminar, LES, RAS (RANS)

// RANS models
RAS
{
    RASModel        kEpsilon;  // kOmegaSST, realizableKE, SpalartAllmaras

    turbulence      on;
    printCoeffs     on;
}

// LES models
LES
{
    LESModel        Smagorinsky;  // kEqn, WALE, dynamicKEqn

    turbulence      on;
    printCoeffs     on;

    delta           cubeRootVol;

    cubeRootVolCoeffs
    {
        deltaCoeff      1;
    }
}
```

## Boundary Conditions

### Common Boundary Types

#### Velocity (U)

**Fixed Value (Inlet)**
```cpp
inlet
{
    type            fixedValue;
    value           uniform (1 0 0);  // 1 m/s in x-direction
}
```

**Zero Gradient (Outlet)**
```cpp
outlet
{
    type            zeroGradient;
}
```

**No Slip Wall**
```cpp
wall
{
    type            noSlip;  // equivalent to fixedValue uniform (0 0 0)
}
```

**Slip Wall**
```cpp
symmetry
{
    type            slip;
}
```

**Inlet with Profile**
```cpp
inlet
{
    type            codedFixedValue;
    value           uniform (0 0 0);
    name            parabolicVelocity;
    code
    #{
        const fvPatch& boundaryPatch = patch();
        const vectorField& Cf = boundaryPatch.Cf();
        vectorField& field = *this;

        forAll(Cf, faceI)
        {
            scalar y = Cf[faceI].y();
            scalar H = 0.1;  // Channel height
            scalar Umax = 1.5;
            field[faceI] = vector(Umax * 4 * y * (H - y) / (H * H), 0, 0);
        }
    #};
}
```

#### Pressure (p)

**Fixed Value (Outlet)**
```cpp
outlet
{
    type            fixedValue;
    value           uniform 0;
}
```

**Zero Gradient (Inlet)**
```cpp
inlet
{
    type            zeroGradient;
}
```

**Wall**
```cpp
wall
{
    type            zeroGradient;
}
```

**Total Pressure (Inlet)**
```cpp
inlet
{
    type            totalPressure;
    p0              uniform 101325;
    U               U;
    phi             phi;
    rho             none;
    psi             none;
    gamma           1;
    value           uniform 101325;
}
```

#### Turbulence (k, epsilon, omega)

**Inlet - Fixed Value**
```cpp
// Turbulent kinetic energy
k_inlet
{
    type            turbulentIntensityKineticEnergyInlet;
    intensity       0.05;  // 5% turbulence intensity
    value           uniform 0.375;
}

// Dissipation rate
epsilon_inlet
{
    type            turbulentMixingLengthDissipationRateInlet;
    mixingLength    0.005;  // 0.5 cm
    value           uniform 14.855;
}

// Specific dissipation rate
omega_inlet
{
    type            turbulentMixingLengthFrequencyInlet;
    mixingLength    0.005;
    value           uniform 440;
}
```

**Wall Functions**
```cpp
// k at wall
k_wall
{
    type            kqRWallFunction;
    value           uniform 0.375;
}

// epsilon at wall
epsilon_wall
{
    type            epsilonWallFunction;
    value           uniform 14.855;
}

// omega at wall
omega_wall
{
    type            omegaWallFunction;
    value           uniform 440;
}

// nut at wall
nut_wall
{
    type            nutkWallFunction;
    value           uniform 0;
}
```

**Outlet**
```cpp
k_outlet
{
    type            zeroGradient;
}

epsilon_outlet
{
    type            zeroGradient;
}

omega_outlet
{
    type            zeroGradient;
}
```

### Calculating Turbulence Parameters

For inlet conditions:

```
Turbulent intensity: I = u'/U = 0.16 * Re^(-1/8)
Turbulent kinetic energy: k = 3/2 * (U * I)²
Turbulent dissipation: ε = C_μ^(3/4) * k^(3/2) / L
Specific dissipation: ω = k^(1/2) / (C_μ^(1/4) * L)

Where:
- U = mean velocity
- Re = Reynolds number
- L = turbulent length scale ≈ 0.07 * characteristic_length
- C_μ = 0.09
```

## Solver Selection

### Incompressible Solvers

**simpleFoam** - Steady-state, incompressible, turbulent
- Applications: External aerodynamics, HVAC, industrial flows
- Algorithm: SIMPLE
- Speed: Fast convergence for steady problems

**pimpleFoam** - Transient, incompressible, turbulent
- Applications: Vortex shedding, transient flows, LES
- Algorithm: PIMPLE (merged PISO-SIMPLE)
- Speed: Good stability, handles large time steps

**icoFoam** - Transient, incompressible, laminar
- Applications: Low Reynolds number, educational
- Algorithm: PISO
- Speed: Fast for laminar flows

**pisoFoam** - Transient, incompressible, turbulent
- Applications: LES, DNS
- Algorithm: PISO
- Speed: Small time steps required

### Multiphase Solvers

**interFoam** - VOF method, two immiscible fluids
- Applications: Free surface flows, sloshing, wave breaking
- Tracks interface between fluids

**multiphaseInterFoam** - VOF for multiple fluids
- Applications: More than two phases

**twoPhaseEulerFoam** - Eulerian-Eulerian multiphase
- Applications: Fluidized beds, bubbly flows

### Compressible Solvers

**rhoSimpleFoam** - Steady-state, compressible, turbulent
- Applications: High-speed flows, compressible effects

**rhoPimpleFoam** - Transient, compressible, turbulent
- Applications: Transient compressible flows

**sonicFoam** - Transient, compressible, supersonic
- Applications: Supersonic flows, shock waves

### Heat Transfer Solvers

**buoyantSimpleFoam** - Steady, buoyancy-driven
- Applications: Natural convection

**buoyantPimpleFoam** - Transient, buoyancy-driven
- Applications: Transient natural convection

**chtMultiRegionFoam** - Conjugate heat transfer
- Applications: Solid-fluid heat transfer

## Mesh Generation

### blockMesh - Structured Hex Meshes

Basic block mesh for rectangular domain:

```cpp
// system/blockMeshDict
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}

convertToMeters 1;

vertices
(
    (0 0 0)      // 0
    (1 0 0)      // 1
    (1 1 0)      // 2
    (0 1 0)      // 3
    (0 0 0.1)    // 4
    (1 0 0.1)    // 5
    (1 1 0.1)    // 6
    (0 1 0.1)    // 7
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (20 20 1) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }
    outlet
    {
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }
    walls
    {
        type wall;
        faces
        (
            (0 1 5 4)
            (3 7 6 2)
        );
    }
    frontAndBack
    {
        type empty;
        faces
        (
            (0 3 2 1)
            (4 5 6 7)
        );
    }
);

mergePatchPairs
(
);
```

Execute: `blockMesh`

### snappyHexMesh - Unstructured Hex-Dominant

For complex geometries from STL files:

```cpp
// system/snappyHexMeshDict (simplified)
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}

castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    geometry.stl
    {
        type triSurfaceMesh;
        name geometry;
    }
}

castellatedMeshControls
{
    maxLocalCells 100000;
    maxGlobalCells 2000000;
    minRefinementCells 0;
    nCellsBetweenLevels 3;

    features
    (
        {
            file "geometry.eMesh";
            level 2;
        }
    );

    refinementSurfaces
    {
        geometry
        {
            level (2 2);
        }
    }

    resolveFeatureAngle 30;
    refinementRegions
    {
    }

    locationInMesh (0.01 0.01 0.01);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 30;
    nRelaxIter 5;
}

addLayersControls
{
    relativeSizes true;
    layers
    {
        geometry
        {
            nSurfaceLayers 3;
        }
    }

    expansionRatio 1.2;
    finalLayerThickness 0.3;
    minThickness 0.1;
    nGrow 0;
    featureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}

meshQualityControls
{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-30;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.02;
    minVolRatio 0.01;
    minTriangleTwist -1;
    nSmoothScale 4;
    errorReduction 0.75;
}

debug 0;
mergeTolerance 1e-6;
```

Execute: `snappyHexMesh -overwrite`

## Post-Processing with ParaView

### Using paraFoam

```bash
# In case directory
paraFoam

# With reconstructed parallel case
paraFoam -builtin

# Touch file for newer versions
touch case.foam
paraview case.foam
```

### Common Visualization Tasks

1. **Velocity Magnitude**
   - Filters → Calculator
   - Result Array Name: UMag
   - Formula: `mag(U)`

2. **Streamlines**
   - Filters → Stream Tracer
   - Seed Type: Point Source or Line Source
   - Vectors: U

3. **Slice/Plane**
   - Filters → Slice
   - Plane Parameters: Origin and Normal

4. **Iso-surfaces**
   - Filters → Contour
   - Select variable and values

5. **Wall Shear Stress**
   - Filters → Calculator
   - Formula: `nu * mag(wallShearStress)`

### Sampling and Data Extraction

Create `system/sampleDict`:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      sampleDict;
}

interpolationScheme cellPoint;

setFormat   raw;

sets
(
    centerline
    {
        type    uniform;
        axis    distance;
        start   (0 0.5 0.05);
        end     (10 0.5 0.05);
        nPoints 100;
    }
);

surfaces
(
    plane
    {
        type    cuttingPlane;
        planeType pointAndNormal;
        pointAndNormalDict
        {
            point   (0 0 0);
            normal  (0 0 1);
        }
        interpolate true;
    }
);

fields  (p U k epsilon);
```

Execute: `postProcess -func sampleDict`

## Parallel Processing

### Decompose Case

Create `system/decomposeParDict`:

```cpp
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}

numberOfSubdomains 4;

method          scotch;  // simple, hierarchical, manual

// For simple decomposition
simpleCoeffs
{
    n               (2 2 1);
    delta           0.001;
}

// For hierarchical decomposition
hierarchicalCoeffs
{
    n               (2 2 1);
    delta           0.001;
    order           xyz;
}
```

### Run Parallel

```bash
# Decompose the mesh
decomposePar

# Run solver in parallel
mpirun -np 4 simpleFoam -parallel

# Reconstruct results
reconstructPar

# Reconstruct only latest time
reconstructPar -latestTime
```

## Workflow Summary

```bash
# 1. Create case structure
mkdir -p myCase/{0,constant,system}

# 2. Copy base configuration
cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily/* myCase/

# 3. Edit dictionaries
# - system/controlDict
# - system/fvSchemes
# - system/fvSolution
# - constant/transportProperties
# - constant/turbulenceProperties
# - 0/U, 0/p, 0/k, 0/epsilon

# 4. Generate mesh
cd myCase
blockMesh

# 5. Check mesh quality
checkMesh

# 6. Run solver
simpleFoam > log.simpleFoam 2>&1

# 7. Monitor residuals
foamMonitor -l postProcessing/residuals/0/residuals.dat

# 8. Post-process
paraFoam

# For parallel:
decomposePar
mpirun -np 4 simpleFoam -parallel
reconstructPar
```

## Best Practices

1. **Mesh Quality**: Always run `checkMesh` before simulation
2. **Residuals**: Monitor convergence (< 1e-4 for engineering accuracy)
3. **Courant Number**: Keep Co < 1 for stability (transient cases)
4. **Y+ Values**: 30-300 for wall functions, < 1 for resolved boundary layers
5. **Grid Independence**: Perform mesh refinement studies
6. **Initialization**: Use `potentialFoam` for better initial velocity field
7. **Backup**: Version control your case files with git

## Troubleshooting

### Common Issues

**Mesh Check Fails**
```bash
checkMesh
# Look for non-orthogonality, skewness
# Increase nNonOrthogonalCorrectors in fvSolution
```

**Divergence**
- Reduce relaxation factors
- Improve initial conditions
- Refine mesh
- Reduce time step (transient)

**Slow Convergence**
- Check boundary conditions
- Adjust solver tolerances
- Use better preconditioners (GAMG for pressure)

**Parallel Issues**
```bash
# Check decomposition
decomposePar -cellDist

# Verify MPI installation
mpirun -np 2 hostname
```

## Resources

- [Official Documentation](https://www.openfoam.com/documentation/)
- [OpenFOAM Wiki](https://openfoamwiki.net/)
- [CFD Online Forums](https://www.cfd-online.com/Forums/openfoam/)
- [User Guide PDF](https://www.openfoam.com/documentation/user-guide)
- [Programmer's Guide](https://www.openfoam.com/documentation/cpp-guide)
- [Tutorials](https://www.openfoam.com/documentation/tutorial-guide)

## See Also

- `examples/pipe-flow/` - Complete turbulent pipe flow case
- `examples/pump-impeller/` - Rotating machinery simulation
- `examples/scripts/` - Python case generation utilities
- `reference.md` - Detailed dictionary and solver reference
