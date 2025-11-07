# OpenFOAM Reference Guide

Comprehensive reference for OpenFOAM dictionaries, solvers, and models.

## Table of Contents
1. [Dictionary Syntax](#dictionary-syntax)
2. [Solver Descriptions](#solver-descriptions)
3. [Turbulence Models](#turbulence-models)
4. [Numerical Schemes](#numerical-schemes)
5. [Linear Solvers](#linear-solvers)
6. [Dimensioned Units](#dimensioned-units)
7. [Function Objects](#function-objects)
8. [Utilities](#utilities)

---

## Dictionary Syntax

### Basic Structure

```cpp
FoamFile
{
    version     2.0;
    format      ascii;  // or binary
    class       dictionary;
    location    "system";
    object      controlDict;
}

// Single line comment

/* Multi-line
   comment */

keyword     value;

dictionaryName
{
    entry1      value1;
    entry2      value2;
}

listName
(
    item1
    item2
    item3
);
```

### Data Types

**Scalars**
```cpp
scalar1     1.5;
scalar2     1.5e-3;
```

**Vectors**
```cpp
vector1     (1 0 0);
vector2     (1.5 -2.3 0.0);
```

**Tensors**
```cpp
tensor1     (1 0 0 0 1 0 0 0 1);  // Identity tensor
```

**Lists**
```cpp
scalarList  (1.0 2.0 3.0 4.0);
vectorList  ((1 0 0) (0 1 0) (0 0 1));
```

**Dimensioned Types**
```cpp
dimensions  [0 1 -1 0 0 0 0];  // Velocity: m/s
value       [0 2 -2 0 0 0 0] 100;  // Pressure: m²/s²
```

### Include Files

```cpp
#include "filename"
#includeIfPresent "optionalFile"
```

### Variable Substitution

```cpp
inlet
{
    type        fixedValue;
    value       $internalField;
}
```

### Conditional Compilation

```cpp
#ifdef VARIABLE
    entry1  value1;
#else
    entry2  value2;
#endif
```

### Dictionary Expansion

```cpp
commonSettings
{
    solver      GAMG;
    tolerance   1e-6;
}

p
{
    $commonSettings;
    relTol      0.01;
}
```

---

## Solver Descriptions

### Incompressible Solvers

#### simpleFoam
- **Type**: Steady-state
- **Compressibility**: Incompressible
- **Turbulence**: RANS/LES
- **Algorithm**: SIMPLE
- **Applications**: External aerodynamics, industrial flows, HVAC
- **Advantages**: Fast convergence, robust
- **Limitations**: Steady-state only

#### pimpleFoam
- **Type**: Transient
- **Compressibility**: Incompressible
- **Turbulence**: RANS/LES
- **Algorithm**: PIMPLE (PISO-SIMPLE)
- **Applications**: Vortex shedding, transient flows, fluid-structure interaction
- **Advantages**: Large time steps, stable
- **Limitations**: Slower than pisoFoam

#### pisoFoam
- **Type**: Transient
- **Compressibility**: Incompressible
- **Turbulence**: RANS/LES/DNS
- **Algorithm**: PISO
- **Applications**: LES, DNS, high-accuracy transient
- **Advantages**: High accuracy
- **Limitations**: Small time steps required

#### icoFoam
- **Type**: Transient
- **Compressibility**: Incompressible
- **Turbulence**: Laminar only
- **Algorithm**: PISO
- **Applications**: Low Reynolds number, educational
- **Advantages**: Simple, fast for laminar
- **Limitations**: Laminar only

#### nonNewtonianIcoFoam
- **Type**: Transient
- **Compressibility**: Incompressible
- **Turbulence**: Laminar
- **Algorithm**: PISO
- **Applications**: Non-Newtonian fluids (polymers, blood)
- **Viscosity Models**: Power law, Cross, Bird-Carreau

#### SRFSimpleFoam
- **Type**: Steady-state
- **Compressibility**: Incompressible
- **Turbulence**: RANS
- **Algorithm**: SIMPLE
- **Applications**: Rotating reference frame (mixers, turbines)
- **Special**: Single rotating frame

### Compressible Solvers

#### rhoSimpleFoam
- **Type**: Steady-state
- **Compressibility**: Compressible
- **Turbulence**: RANS
- **Algorithm**: SIMPLE
- **Applications**: High-speed subsonic, transonic flows
- **Mach Range**: 0.1 - 3.0

#### rhoPimpleFoam
- **Type**: Transient
- **Compressibility**: Compressible
- **Turbulence**: RANS/LES
- **Algorithm**: PIMPLE
- **Applications**: Transient compressible flows
- **Special**: Includes energy equation

#### sonicFoam
- **Type**: Transient
- **Compressibility**: Compressible (supersonic)
- **Turbulence**: RANS/LES
- **Algorithm**: PISO
- **Applications**: Supersonic flows, shock waves
- **Mach Range**: > 1.0

#### rhoCentralFoam
- **Type**: Transient
- **Compressibility**: Compressible (high Mach)
- **Turbulence**: LES/DNS
- **Algorithm**: Central differencing (Kurganov-Tadmor)
- **Applications**: Supersonic/hypersonic flows, shock capturing
- **Mach Range**: 0.5 - 20+

### Multiphase Solvers

#### interFoam
- **Type**: Transient
- **Method**: VOF (Volume of Fluid)
- **Phases**: 2 immiscible
- **Applications**: Free surface, sloshing, wave breaking, dam break
- **Interface**: Sharp interface capturing

#### multiphaseInterFoam
- **Type**: Transient
- **Method**: VOF
- **Phases**: Multiple immiscible
- **Applications**: More than 2 phases
- **Interface**: Multiple interfaces

#### twoPhaseEulerFoam
- **Type**: Transient
- **Method**: Eulerian-Eulerian
- **Phases**: 2 (gas-solid, liquid-gas)
- **Applications**: Fluidized beds, bubbly flows, particle-laden flows
- **Models**: Drag, lift, virtual mass, turbulent dispersion

#### multiphaseEulerFoam
- **Type**: Transient
- **Method**: Eulerian-Eulerian
- **Phases**: Multiple
- **Applications**: Complex multiphase reactors
- **Models**: Inter-phase momentum/heat/mass transfer

#### driftFluxFoam
- **Type**: Transient
- **Method**: Drift-flux model
- **Phases**: 2 (mixture model)
- **Applications**: Sedimentation, flotation
- **Advantages**: Computationally efficient

### Heat Transfer Solvers

#### buoyantSimpleFoam
- **Type**: Steady-state
- **Compressibility**: Compressible
- **Turbulence**: RANS
- **Applications**: Natural convection, buoyancy-driven flows
- **Special**: Boussinesq approximation

#### buoyantPimpleFoam
- **Type**: Transient
- **Compressibility**: Compressible
- **Turbulence**: RANS/LES
- **Applications**: Transient natural convection, fires
- **Special**: Density variation with temperature

#### chtMultiRegionFoam
- **Type**: Steady/Transient
- **Regions**: Multiple (solid and fluid)
- **Applications**: Conjugate heat transfer, electronics cooling
- **Special**: Coupled solid-fluid heat transfer

#### laplacianFoam
- **Type**: Transient
- **Equation**: Laplace equation
- **Applications**: Heat conduction, potential flow
- **Special**: Simple diffusion solver

### Combustion Solvers

#### reactingFoam
- **Type**: Transient
- **Compressibility**: Compressible
- **Chemistry**: Finite-rate chemistry
- **Applications**: Combustion, reacting flows
- **Models**: EDC, PaSR, laminar finite rate

#### fireFoam
- **Type**: Transient
- **Compressibility**: Compressible
- **Applications**: Fire simulation, pyrolysis
- **Special**: Radiation, soot modeling

### Particle Tracking

#### icoUncoupledKinematicParcelFoam
- **Type**: Transient
- **Coupling**: One-way (fluid → particles)
- **Applications**: Spray, aerosol
- **Particle Types**: Kinematic

#### sprayFoam
- **Type**: Transient
- **Coupling**: Two-way
- **Applications**: Spray combustion, fuel injection
- **Special**: Lagrangian particles, breakup models

---

## Turbulence Models

### RANS Models

#### k-epsilon Models

**Standard k-epsilon**
```cpp
RASModel        kEpsilon;
```
- **Applications**: General purpose, industrial flows
- **Advantages**: Robust, computationally efficient
- **Limitations**: Poor for swirl, separation, adverse pressure gradients
- **Y+ requirement**: 30-300 (wall functions)

**Realizable k-epsilon**
```cpp
RASModel        realizableKE;
```
- **Improvements**: Better for separation, jets, mixing layers
- **Advantages**: More accurate than standard k-ε
- **Applications**: Complex flows with strong streamline curvature

**RNG k-epsilon**
```cpp
RASModel        RNGkEpsilon;
```
- **Theory**: Renormalization group theory
- **Improvements**: Low Reynolds number effects, swirl
- **Applications**: Swirling flows, separation

#### k-omega Models

**Standard k-omega**
```cpp
RASModel        kOmega;
```
- **Applications**: Boundary layers, adverse pressure gradients
- **Advantages**: Good near-wall behavior
- **Limitations**: Sensitive to freestream values

**k-omega SST (Shear Stress Transport)**
```cpp
RASModel        kOmegaSST;
```
- **Applications**: Aerodynamics, turbomachinery, separation
- **Advantages**: Best of k-ε and k-ω, accurate for adverse pressure gradients
- **Industry**: Widely used in aerospace
- **Y+ requirement**: < 1 for full resolution, 30-300 for wall functions

#### Reynolds Stress Models

**Launder-Reece-Rodi (LRR)**
```cpp
RASModel        LRR;
```
- **Equations**: 7 (6 stresses + ε)
- **Applications**: Complex 3D flows, strong anisotropy
- **Advantages**: Most accurate RANS
- **Limitations**: Computationally expensive, less stable

#### Spalart-Allmaras

**Spalart-Allmaras**
```cpp
RASModel        SpalartAllmaras;
```
- **Equations**: 1 (modified turbulent viscosity)
- **Applications**: Aerodynamics, boundary layers
- **Advantages**: Robust, fast, good for attached flows
- **Limitations**: Tuned for aerospace applications

### LES Models

#### Smagorinsky

```cpp
LESModel        Smagorinsky;

SmagorinskyCoeffs
{
    Ce              1.048;
    Ck              0.094;
}
```
- **Applications**: General purpose LES
- **Advantages**: Simple, stable
- **Limitations**: Too dissipative near walls

#### WALE (Wall-Adapting Local Eddy-viscosity)

```cpp
LESModel        WALE;

WALECoeffs
{
    Ce              1.048;
    Cw              0.325;
}
```
- **Applications**: Wall-bounded flows
- **Advantages**: Correct near-wall behavior, no wall damping needed
- **Better than**: Smagorinsky for wall-bounded flows

#### Dynamic k-equation

```cpp
LESModel        dynamicKEqn;
```
- **Applications**: Complex geometries
- **Advantages**: Dynamic coefficient calculation
- **Transport**: Solves equation for subgrid k

### Wall Treatment

**High Reynolds Number (Wall Functions)**
- Y+ range: 30-300
- Models: Standard, non-equilibrium, low-Re corrections
- Advantages: Coarse mesh near walls, computationally efficient
- Use with: k-ε, k-ω SST (wall function mode)

**Low Reynolds Number (Resolved)**
- Y+ range: < 1 (ideally < 0.5)
- Models: Low-Re k-ε, k-ω SST
- Advantages: Accurate near-wall prediction
- Disadvantages: Fine mesh required, computationally expensive

### Turbulence Model Selection Guide

| Application | Recommended Model | Alternative |
|-------------|------------------|-------------|
| External aerodynamics | k-ω SST | Spalart-Allmaras |
| Internal flows | k-ε, k-ω SST | Realizable k-ε |
| Separation | k-ω SST | RSM |
| Rotating machinery | k-ω SST | k-ε |
| Heat transfer | k-ω SST, Low-Re | v2-f |
| Free shear flows | k-ε | k-ω SST |
| Buoyancy | k-ε | k-ω SST |
| High accuracy | LES (WALE) | DNS |

---

## Numerical Schemes

### Time Derivatives (ddtSchemes)

**Steady State**
```cpp
default         steadyState;
```

**First Order Implicit (Euler)**
```cpp
default         Euler;
```
- Accuracy: First order
- Stability: Unconditionally stable
- Time step: Moderate

**Second Order Implicit (Backward)**
```cpp
default         backward;
```
- Accuracy: Second order
- Stability: Stable
- Time step: Moderate to large

**Crank-Nicolson**
```cpp
default         CrankNicolson 0.5;  // 0.5 = pure CN, 1.0 = pure Euler
```
- Accuracy: Second order
- Stability: Conditionally stable
- Time step: Large for low values of coefficient

### Gradient Schemes (gradSchemes)

**Gauss Linear (Second Order)**
```cpp
default         Gauss linear;
```
- Accuracy: Second order
- Best for: General purpose

**Cell Limited Gauss Linear**
```cpp
grad(U)         cellLimited Gauss linear 1;
```
- Limiting coefficient: 0-2 (1 = moderate limiting)
- Use for: Velocity gradients, stability

**Least Squares**
```cpp
default         leastSquares;
```
- Accuracy: Second order
- Use for: Unstructured meshes

### Divergence Schemes (divSchemes)

**Linear Upwind (Second Order)**
```cpp
div(phi,U)      Gauss linearUpwind grad(U);
```
- Accuracy: Second order
- Stability: Good
- Use for: Velocity

**Bounded Linear Upwind**
```cpp
div(phi,U)      bounded Gauss linearUpwind grad(U);
```
- Bounded: Yes (prevents overshoots)
- Use for: All convection terms

**Upwind (First Order)**
```cpp
div(phi,k)      Gauss upwind;
```
- Accuracy: First order
- Stability: Very stable
- Use for: Initial convergence, turbulence

**Linear (Central Differencing)**
```cpp
div(phi,U)      Gauss linear;
```
- Accuracy: Second order
- Stability: Can be unstable
- Use for: Low Peclet number, LES

**LUST (LES)**
```cpp
div(phi,U)      Gauss LUST grad(U);
```
- Accuracy: High order
- Use for: LES, DNS

**Schemes for Turbulence**
```cpp
div(phi,k)              bounded Gauss upwind;
div(phi,epsilon)        bounded Gauss upwind;
div(phi,omega)          bounded Gauss upwind;
```

### Laplacian Schemes (laplacianSchemes)

**Linear Corrected (Second Order)**
```cpp
default         Gauss linear corrected;
```
- Accuracy: Second order
- Non-orthogonal correction: Yes
- Use for: Diffusion terms

**Limited (Bounded)**
```cpp
default         Gauss linear limited 0.5;
```
- Limiting: 0-1 (0.5 = moderate)
- Use for: Poor quality meshes

**Uncorrected (First Order)**
```cpp
default         Gauss linear uncorrected;
```
- Accuracy: First order
- Use for: Orthogonal meshes only

### Interpolation Schemes

**Linear (Second Order)**
```cpp
default         linear;
```
- Standard choice

**Upwind (First Order)**
```cpp
interpolate(U)  upwind phi;
```
- More stable

### Surface Normal Gradient (snGradSchemes)

**Corrected (Second Order)**
```cpp
default         corrected;
```
- Non-orthogonal correction: Yes

**Limited**
```cpp
default         limited 0.5;
```
- Use for: Poor meshes

---

## Linear Solvers

### Pressure Solvers

**GAMG (Geometric-Algebraic Multi-Grid)**
```cpp
p
{
    solver          GAMG;
    tolerance       1e-6;
    relTol          0.01;
    smoother        GaussSeidel;
    nPreSweeps      0;
    nPostSweeps     2;
    cacheAgglomeration true;
    nCellsInCoarsestLevel 10;
    agglomerator    faceAreaPair;
    mergeLevels     1;
}
```
- **Best for**: Pressure equation
- **Performance**: Very fast for large meshes
- **Scalability**: Excellent parallel scaling

**PCG (Preconditioned Conjugate Gradient)**
```cpp
p
{
    solver          PCG;
    preconditioner  DIC;
    tolerance       1e-6;
    relTol          0.01;
}
```
- **Best for**: Symmetric positive-definite matrices
- **Preconditioners**: DIC, DILU, FDIC, GAMG

### Velocity and Scalar Solvers

**smoothSolver**
```cpp
U
{
    solver          smoothSolver;
    smoother        GaussSeidel;
    tolerance       1e-7;
    relTol          0.1;
    nSweeps         1;
}
```
- **Best for**: Velocity, turbulence, scalars
- **Smoothers**: GaussSeidel, symGaussSeidel, DIC, DILU

**PBiCGStab (Preconditioned Bi-Conjugate Gradient Stabilized)**
```cpp
U
{
    solver          PBiCGStab;
    preconditioner  DILU;
    tolerance       1e-7;
    relTol          0.1;
}
```
- **Best for**: Asymmetric matrices
- **Preconditioners**: DILU, DIAGONAL

### Preconditioners

| Preconditioner | Description | Use Case |
|----------------|-------------|----------|
| **DIC** | Diagonal Incomplete Cholesky | Symmetric, pressure |
| **DILU** | Diagonal Incomplete LU | Asymmetric, velocity |
| **FDIC** | Faster DIC | Symmetric, faster |
| **GAMG** | Multi-grid | Large systems |
| **DIAGONAL** | Diagonal | Simple, fast |
| **none** | No preconditioning | Small systems |

### Solver Tolerances

**Absolute Tolerance**
```cpp
tolerance       1e-6;
```
- Residual must drop below this value

**Relative Tolerance**
```cpp
relTol          0.01;
```
- Residual must drop by this factor (1% = 0.01)

**Typical Values**
- Pressure: tolerance 1e-6, relTol 0.01
- Velocity: tolerance 1e-7, relTol 0.1
- Turbulence: tolerance 1e-8, relTol 0.1

---

## Dimensioned Units

### Dimension Vector

Format: `[mass length time temperature quantity current luminosity]`

| Physical Quantity | Dimensions | Example Value |
|-------------------|------------|---------------|
| **Length** | [0 1 0 0 0 0 0] | meter |
| **Time** | [0 0 1 0 0 0 0] | second |
| **Mass** | [1 0 0 0 0 0 0] | kilogram |
| **Temperature** | [0 0 0 1 0 0 0] | Kelvin |
| **Velocity** | [0 1 -1 0 0 0 0] | m/s |
| **Pressure** | [0 2 -2 0 0 0 0] | m²/s² (kinematic) |
| **Pressure** | [1 -1 -2 0 0 0 0] | kg/(m·s²) = Pa |
| **Force** | [1 1 -2 0 0 0 0] | kg·m/s² = N |
| **Kinematic Viscosity** | [0 2 -1 0 0 0 0] | m²/s |
| **Dynamic Viscosity** | [1 -1 -1 0 0 0 0] | kg/(m·s) = Pa·s |
| **Density** | [1 -3 0 0 0 0 0] | kg/m³ |
| **Energy** | [1 2 -2 0 0 0 0] | kg·m²/s² = J |
| **Power** | [1 2 -3 0 0 0 0] | kg·m²/s³ = W |
| **Heat Flux** | [1 0 -3 0 0 0 0] | kg/s³ = W/m² |
| **Thermal Conductivity** | [1 1 -3 -1 0 0 0] | kg·m/(s³·K) = W/(m·K) |
| **Specific Heat** | [0 2 -2 -1 0 0 0] | m²/(s²·K) = J/(kg·K) |
| **Turb. Kinetic Energy (k)** | [0 2 -2 0 0 0 0] | m²/s² |
| **Dissipation Rate (ε)** | [0 2 -3 0 0 0 0] | m²/s³ |
| **Spec. Dissipation (ω)** | [0 0 -1 0 0 0 0] | 1/s |
| **Angular Velocity** | [0 0 -1 0 0 0 0] | rad/s |
| **Dimensionless** | [0 0 0 0 0 0 0] | - |

### Common Fluid Properties

**Air (20°C, 1 atm)**
```cpp
nu      [0 2 -1 0 0 0 0] 1.5e-05;    // m²/s
rho     [1 -3 0 0 0 0 0] 1.225;      // kg/m³
mu      [1 -1 -1 0 0 0 0] 1.8e-05;   // Pa·s
cp      [0 2 -2 -1 0 0 0] 1005;      // J/(kg·K)
k       [1 1 -3 -1 0 0 0] 0.0257;    // W/(m·K)
Pr      [0 0 0 0 0 0 0] 0.7;         // Prandtl number
```

**Water (20°C)**
```cpp
nu      [0 2 -1 0 0 0 0] 1.0e-06;    // m²/s
rho     [1 -3 0 0 0 0 0] 998;        // kg/m³
mu      [1 -1 -1 0 0 0 0] 0.001;     // Pa·s
cp      [0 2 -2 -1 0 0 0] 4182;      // J/(kg·K)
k       [1 1 -3 -1 0 0 0] 0.598;     // W/(m·K)
Pr      [0 0 0 0 0 0 0] 7.0;         // Prandtl number
```

---

## Function Objects

Function objects enable in-situ data processing during simulation runtime.

### Forces and Coefficients

**forces**
```cpp
forces
{
    type            forces;
    libs            ("libforces.so");
    writeControl    timeStep;
    writeInterval   1;

    patches         (wall);
    rho             rhoInf;      // For incompressible
    rhoInf          1.225;       // kg/m³
    CofR            (0 0 0);     // Center of rotation
}
```

**forceCoeffs**
```cpp
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
    lRef            1;           // Reference length
    Aref            1;           // Reference area
}
```

### Sampling

**probes**
```cpp
probes
{
    type            probes;
    libs            ("libsampling.so");
    writeControl    timeStep;
    writeInterval   1;

    fields          (p U);
    probeLocations
    (
        (1 0.5 0.5)
        (2 0.5 0.5)
        (3 0.5 0.5)
    );
}
```

### Field Averaging

**fieldAverage**
```cpp
fieldAverage
{
    type            fieldAverage;
    libs            ("libfieldFunctionObjects.so");
    writeControl    writeTime;

    fields
    (
        U
        {
            mean        on;
            prime2Mean  on;
            base        time;
        }
        p
        {
            mean        on;
            prime2Mean  off;
            base        time;
        }
    );
}
```

### Residuals

**residuals**
```cpp
residuals
{
    type            residuals;
    libs            ("libutilityFunctionObjects.so");
    writeControl    timeStep;
    writeInterval   1;
    fields          (p U k epsilon);
}
```

---

## Utilities

### Pre-Processing

| Utility | Description |
|---------|-------------|
| **blockMesh** | Generate structured hex meshes |
| **snappyHexMesh** | Generate unstructured hex-dominant meshes from STL |
| **extrudeMesh** | Extrude 2D to 3D |
| **refineMesh** | Refine existing mesh |
| **checkMesh** | Check mesh quality |
| **transformPoints** | Transform mesh coordinates |
| **createPatch** | Create/modify patches |
| **topoSet** | Create cell/face/point sets |
| **setFields** | Initialize field values |
| **decomposePar** | Decompose for parallel |

### Post-Processing

| Utility | Description |
|---------|-------------|
| **paraFoam** | Launch ParaView |
| **foamToVTK** | Convert to VTK format |
| **sample** | Extract data along lines/planes |
| **postProcess** | Run function objects post-simulation |
| **foamLog** | Extract residuals from log |
| **foamCalc** | Calculate derived fields |
| **reconstructPar** | Reconstruct parallel results |
| **execFlowFunctionObjects** | Execute function objects |

### Manipulation

| Utility | Description |
|---------|-------------|
| **mapFields** | Map results between meshes |
| **moveDynamicMesh** | Move mesh according to dictionary |
| **renumberMesh** | Optimize mesh numbering |
| **stitchMesh** | Stitch two meshes together |
| **mergeMeshes** | Merge multiple meshes |
| **splitMeshRegions** | Split into disconnected regions |

---

## Links to Documentation

### Official Resources

- **OpenFOAM.org**: https://www.openfoam.org/
- **OpenFOAM.com**: https://www.openfoam.com/
- **User Guide**: https://www.openfoam.com/documentation/user-guide
- **Programmer's Guide**: https://www.openfoam.com/documentation/cpp-guide
- **Tutorials**: https://www.openfoam.com/documentation/tutorial-guide

### Community Resources

- **CFD Online**: https://www.cfd-online.com/Forums/openfoam/
- **OpenFOAM Wiki**: https://openfoamwiki.net/
- **GitHub**: https://github.com/OpenFOAM/OpenFOAM-11

### Learning Resources

- **OpenFOAM Course**: http://www.wolfdynamics.com/training.html
- **YouTube Tutorials**: Search "OpenFOAM tutorial"
- **Books**: "The OpenFOAM Technology Primer" by Moukalled et al.

### Validation and Benchmarks

- **NASA TMR**: Turbulence Modeling Resource
- **NPARC Alliance**: Validation Archive
- **ERCOFTAC**: Database of experimental data
