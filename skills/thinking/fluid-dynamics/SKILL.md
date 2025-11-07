---
name: fluid-dynamics-workflow
description: "Systematic workflow for fluid dynamics analysis from problem setup to validation"
category: thinking
domain: fluids
complexity: advanced
dependencies: []
---

# Fluid Dynamics Analysis Workflow

A comprehensive, systematic approach to computational fluid dynamics (CFD) and analytical fluid mechanics problems.

## Overview

This workflow guides you through complete fluid dynamics analysis, from initial problem formulation through validation. It ensures rigorous treatment of governing physics, appropriate model selection, and proper verification of results.

## Complete Analysis Workflow

### 1. Problem Definition

#### 1.1 Geometry and Domain
- **Physical Domain**:
  - Define computational domain extent
  - Identify key geometric features (boundary layers, separation zones, recirculation)
  - Determine symmetry planes to reduce computational cost
  - Establish coordinate system (Cartesian, cylindrical, spherical)
  - Assess 2D vs 3D requirements

- **Domain Size Considerations**:
  - Inlet: 5-10 characteristic lengths upstream
  - Outlet: 20-30 characteristic lengths downstream
  - Lateral boundaries: 5-10 characteristic lengths from object
  - Verify no artificial confinement effects

#### 1.2 Flow Regime Characterization

**Reynolds Number Analysis**:
```
Re = ρVL/μ = VL/ν

Where:
- ρ: density
- V: characteristic velocity
- L: characteristic length
- μ: dynamic viscosity
- ν: kinematic viscosity
```

**Flow Regime Classification**:
- Re < 2300: Laminar flow (pipes)
- 2300 < Re < 4000: Transitional
- Re > 4000: Turbulent (pipes)
- External flows:
  - Re < 5×10⁵: Laminar boundary layer
  - Re > 5×10⁵: Turbulent boundary layer

**Mach Number Analysis**:
```
Ma = V/c

Where:
- V: flow velocity
- c: speed of sound = √(γRT)
```

**Compressibility Classification**:
- Ma < 0.3: Incompressible (ρ = constant)
- 0.3 < Ma < 0.8: Subsonic compressible
- 0.8 < Ma < 1.2: Transonic
- 1.2 < Ma < 5.0: Supersonic
- Ma > 5.0: Hypersonic

**Froude Number** (free surface flows):
```
Fr = V/√(gL)
```
- Fr < 1: Subcritical flow
- Fr = 1: Critical flow
- Fr > 1: Supercritical flow

#### 1.3 Additional Dimensionless Numbers

**Strouhal Number** (unsteady flows):
```
St = fL/V
```
- Characterizes oscillating flow mechanisms
- Vortex shedding: St ≈ 0.2 for cylinders

**Prandtl Number** (heat transfer):
```
Pr = ν/α = μcp/k
```

**Grashof Number** (natural convection):
```
Gr = gβΔTL³/ν²
```

### 2. Governing Equations

#### 2.1 Continuity Equation

**Incompressible**:
```
∇·V = 0
```

**Compressible**:
```
∂ρ/∂t + ∇·(ρV) = 0
```

#### 2.2 Momentum Equations (Navier-Stokes)

**Incompressible, Newtonian fluid**:
```
ρ(∂V/∂t + V·∇V) = -∇p + μ∇²V + ρg
```

**Component form (Cartesian)**:
```
ρ(∂u/∂t + u∂u/∂x + v∂u/∂y + w∂u/∂z) = -∂p/∂x + μ(∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²) + ρgₓ
ρ(∂v/∂t + u∂v/∂x + v∂v/∂y + w∂v/∂z) = -∂p/∂y + μ(∂²v/∂x² + ∂²v/∂y² + ∂²v/∂z²) + ρgᵧ
ρ(∂w/∂t + u∂w/∂x + v∂w/∂y + w∂w/∂z) = -∂p/∂z + μ(∂²w/∂x² + ∂²w/∂y² + ∂²w/∂z²) + ρgᵣ
```

**Compressible form**:
```
ρ(∂V/∂t + V·∇V) = -∇p + ∇·τ + ρg
```

Where stress tensor:
```
τᵢⱼ = μ(∂uᵢ/∂xⱼ + ∂uⱼ/∂xᵢ) + λ(∇·V)δᵢⱼ
```

#### 2.3 Energy Equation

**For compressible flows**:
```
ρcp(∂T/∂t + V·∇T) = ∇·(k∇T) + Φ + Q
```

Where:
- Φ: viscous dissipation
- Q: heat source term

**Total energy form**:
```
∂(ρE)/∂t + ∇·(ρEV) = -∇·(pV) + ∇·(k∇T) + ∇·(τ·V)
```

#### 2.4 Turbulence Modeling

**Reynolds-Averaged Navier-Stokes (RANS)**:
Decompose variables: u = ū + u'

**Reynolds Stress Tensor**:
```
τᵢⱼᴿᴬᴺˢ = -ρu'ᵢu'ⱼ
```

**Closure Models Decision Tree**:

```
START
  │
  ├─ Simple geometry, equilibrium turbulence?
  │  └─ YES → k-ε Standard
  │
  ├─ Wall-bounded flow, adverse pressure gradient?
  │  └─ YES → k-ω SST (Menter)
  │
  ├─ Separated flow, vortex shedding?
  │  └─ YES → k-ω SST or Transition SST
  │
  ├─ Swirling flow, rotation?
  │  └─ YES → Reynolds Stress Model (RSM)
  │
  ├─ High accuracy needed, large-scale unsteadiness?
  │  └─ YES → Large Eddy Simulation (LES)
  │
  └─ Research, turbulence structure critical?
     └─ YES → Direct Numerical Simulation (DNS)
```

**Model Comparison**:

| Model | Advantages | Limitations | Best Use Cases |
|-------|-----------|-------------|----------------|
| k-ε Standard | Robust, well-validated | Poor near walls, separated flows | Industrial flows, simple geometries |
| k-ε Realizable | Better for separated flows | Still needs wall functions | Complex geometries, separation |
| k-ω Standard | Good near walls | Sensitive to freestream | Wall-bounded flows |
| k-ω SST | Excellent all-around | Higher computational cost | Aerodynamics, adverse gradients |
| Spalart-Allmaras | One equation, efficient | Calibrated for aerospace | Airfoils, external aerodynamics |
| RSM | Captures anisotropy | 5-7 equations, complex | Swirling, rotating flows |
| LES | Resolves large eddies | Very expensive | Unsteady flows, acoustics |
| DNS | No modeling errors | Extremely expensive | Fundamental research, Re < 10⁴ |

### 3. Boundary Conditions

#### 3.1 Inlet Boundary Conditions

**Velocity Inlet**:
- Specify: u, v, w (velocity components)
- For turbulent flows, also specify:
  - Turbulent intensity: I = u'/U (typically 1-10%)
  - Turbulent length scale: l = 0.07L
  - Or directly: k, ε (or ω)

**Calculating turbulence quantities**:
```
k = 3/2(UI)²
ε = C_μ^(3/4) k^(3/2) / l
ω = k^(1/2) / (C_μ^(1/4) l)
```
Where C_μ = 0.09

**Pressure Inlet**:
- Specify: total pressure p₀
- Flow direction
- Temperature (compressible)

**Mass Flow Inlet**:
- Specify: ṁ = ρVA
- Flow direction

#### 3.2 Outlet Boundary Conditions

**Pressure Outlet**:
- Specify static pressure
- Zero gradient for all other quantities: ∂φ/∂n = 0
- Most common and stable

**Outflow**:
- Fully developed flow assumption
- ∂u/∂x = 0 (all variables)
- Use when outlet is far downstream

**Outlet Distance Requirement**:
- Minimum 10-15 hydraulic diameters downstream
- For separated flows: 20-30 diameters

#### 3.3 Wall Boundary Conditions

**No-Slip Wall** (standard):
- V = 0 at wall
- Temperature specified (isothermal) or heat flux

**Slip Wall**:
- Tangential velocity ≠ 0
- Normal velocity = 0
- V·n = 0, ∂V_tangential/∂n = 0

**Wall Function Approach**:
For turbulent flows, use wall functions when:
- 30 < y⁺ < 300

```
y⁺ = ρu_τy/μ
u_τ = √(τ_w/ρ)  [friction velocity]
```

**Near-Wall Resolution Requirements**:

| Approach | y⁺ first cell | Total prism layers | Growth rate |
|----------|--------------|-------------------|-------------|
| Wall Functions | 30-300 | 3-5 | 1.2-1.3 |
| Low-Re Model | y⁺ < 1 | 10-20 | 1.1-1.2 |
| Enhanced Wall Treatment | y⁺ < 5 | 8-15 | 1.15-1.25 |

**Moving Wall**:
- Specify wall velocity
- Applications: rotating machinery, conveyors

#### 3.4 Symmetry and Periodic Boundaries

**Symmetry**:
- Normal velocity = 0
- Zero gradient for other quantities
- ∂φ/∂n = 0

**Periodic**:
- Translational: φ(x₁) = φ(x₂)
- Rotational: for rotating geometries
- Use to reduce computational domain

#### 3.5 Free Surface (Multiphase)

**Volume of Fluid (VOF)**:
- Track interface with volume fraction
- α = 1 (phase 1), α = 0 (phase 2)

**Level Set Method**:
- Use signed distance function
- Sharp interface tracking

### 4. Discretization

#### 4.1 Mesh Requirements

**Mesh Types**:

| Mesh Type | Advantages | Disadvantages | Best For |
|-----------|-----------|---------------|----------|
| Structured (Hexahedral) | Efficient, accurate | Complex geometry difficult | Simple domains, channels |
| Unstructured (Tetrahedral) | Flexible geometry | More cells needed | Complex CAD geometries |
| Hybrid (Prism + Tet) | Best of both | Requires careful setup | Most CFD applications |
| Polyhedral | Fewer cells, stable | Specialized mesher | Industrial CFD |

**Mesh Quality Metrics**:
- Aspect ratio < 100 (< 20 preferred)
- Skewness < 0.85 (< 0.5 preferred)
- Orthogonality > 0.15 (> 0.7 preferred)

**Grid Independence Study**:
```
Refinement ratio: r = 1.5 to 2.0
Test at least 3 mesh levels:
- Coarse: N cells
- Medium: r²N cells
- Fine: r⁴N cells

Check convergence of key outputs:
- Forces (drag, lift)
- Pressures at critical points
- Flow rates
```

**Richardson Extrapolation**:
```
f_exact ≈ f_fine + (f_fine - f_coarse)/(r^p - 1)
```
Where p is the order of accuracy

#### 4.2 Y⁺ Values for Turbulence

**Critical y⁺ ranges**:
```
Viscous sublayer:     y⁺ < 5
Buffer layer:         5 < y⁺ < 30
Log-layer:            30 < y⁺ < 300
Outer layer:          y⁺ > 300
```

**First cell height calculation**:
```
y = y⁺μ/(ρu_τ)

For flat plate:
u_τ ≈ U∞√(C_f/2)
C_f = 0.058Re_L^(-0.2)  [turbulent]
```

**Estimating y for desired y⁺**:
```
Example: Airfoil at Re = 1×10⁶, V = 50 m/s
1. Estimate C_f ≈ 0.0037
2. u_τ = 50√(0.0037/2) = 2.15 m/s
3. For y⁺ = 1: y = 1×1.5e-5/(1.225×2.15) ≈ 5.7e-6 m

First cell: 5.7 μm
```

#### 4.3 Time Step Selection (Unsteady)

**CFL Condition** (Courant-Friedrichs-Lewy):
```
CFL = VΔt/Δx ≤ CFL_max

For explicit schemes: CFL < 1
For implicit schemes: CFL < 10-20 (stability)
                      CFL < 1-5 (accuracy)
```

**Time step calculation**:
```
Δt = CFL × Δx_min / V_max
```

**Temporal Resolution for Unsteady Features**:
- Capture vortex shedding: Δt < T/100 (T = shedding period)
- Transient analysis: Δt < τ/50 (τ = characteristic time)

**Dual Time Stepping**:
- Physical time step: for accuracy
- Inner iterations: for convergence

### 5. Solution Strategy

#### 5.1 Solver Selection

**Pressure-Velocity Coupling**:

| Algorithm | Type | Advantages | Use Cases |
|-----------|------|-----------|-----------|
| SIMPLE | Segregated | Memory efficient | Steady state |
| SIMPLEC | Segregated | Faster convergence | Steady/transient |
| PISO | Segregated | Transient accuracy | Highly unsteady |
| Coupled | Fully-coupled | Robust, fast | Compressible, complex |

**Momentum Discretization**:
- First-order upwind: Stable but diffusive
- Second-order upwind: Good accuracy
- QUICK: Higher accuracy, bounded
- Central differencing: Low diffusion, needs stability

**Pressure Interpolation**:
- Standard: Most cases
- PRESTO: Swirling, high-speed rotating flows
- Body-force-weighted: Buoyancy-driven flows

#### 5.2 Convergence Criteria

**Residual Monitoring**:
```
Continuity:   < 1e-4
Velocity:     < 1e-4
Turbulence:   < 1e-4
Energy:       < 1e-6
```

**Additional Convergence Checks**:
- Monitor force coefficients (C_L, C_D)
- Track mass flow balance (in = out)
- Check integrated quantities
- Visualize flow field for reasonableness

**Convergence plots should show**:
- Smooth decrease in residuals
- Asymptotic approach to steady value
- All residuals dropping together

#### 5.3 Under-Relaxation Factors

**Purpose**: Stabilize solution, prevent divergence

**Typical values (SIMPLE)**:
```
Pressure:         0.3
Momentum:         0.7
Turbulence (k):   0.8
Turbulence (ε,ω): 0.8
Energy:           0.9
```

**Adjustment strategy**:
- Lower for unstable problems
- Increase gradually as solution develops
- Coupled solver: higher factors possible

### 6. Post-Processing and Validation

#### 6.1 Flow Visualization

**Essential Visualizations**:

1. **Velocity Contours/Vectors**:
   - Identify recirculation zones
   - Check boundary layer development
   - Verify flow direction

2. **Pressure Distribution**:
   - Surface pressure coefficient: C_p = (p - p_∞)/(½ρV²)
   - Identify pressure gradients
   - Check for unphysical oscillations

3. **Streamlines/Pathlines**:
   - Visualize flow patterns
   - Identify separation/reattachment
   - Check for reverse flow at outlets

4. **Vorticity**:
   - ω = ∇ × V
   - Identify vortex structures
   - Turbulence visualization

5. **Wall Shear Stress**:
   - τ_w = μ(∂u/∂y)_wall
   - Separation: τ_w = 0
   - Friction coefficient: C_f = τ_w/(½ρV²)

6. **Turbulence Quantities**:
   - Turbulent kinetic energy (k)
   - Turbulent viscosity ratio (μ_t/μ)
   - y⁺ distribution

#### 6.2 Force Calculations

**Pressure and Viscous Forces**:
```
F_pressure = ∫∫ p·n dA
F_viscous = ∫∫ τ·n dA
F_total = F_pressure + F_viscous
```

**Non-dimensional Coefficients**:
```
C_D = D/(½ρV²A)    [Drag coefficient]
C_L = L/(½ρV²A)    [Lift coefficient]
C_M = M/(½ρV²AL)   [Moment coefficient]
```

**Pressure Coefficient**:
```
C_p = (p - p_∞)/(½ρV²)
```

**Friction Coefficient**:
```
C_f = τ_w/(½ρV²)
```

#### 6.3 Validation and Verification

**Verification (Solving equations right)**:

1. **Grid Convergence Index (GCI)**:
```
GCI = F_s |ε|/(r^p - 1)

Where:
- F_s = safety factor (1.25 for 3+ grids)
- ε = (f_fine - f_coarse)/f_fine
- r = refinement ratio
- p = order of accuracy
```

2. **Temporal Convergence**:
- Repeat with Δt/2
- Verify time-step independence

3. **Iterative Convergence**:
- Ensure residuals decrease adequately
- Monitor solution variables

**Validation (Solving right equations)**:

1. **Analytical Solutions**:
   - Poiseuille flow (channel)
   - Couette flow
   - Blasius boundary layer
   - Potential flow over cylinder

2. **Experimental Data**:
   - Wind tunnel measurements
   - PIV (Particle Image Velocimetry)
   - Hot-wire anemometry

3. **Benchmark Cases**:
   - Lid-driven cavity
   - Flow over cylinder
   - Backward-facing step
   - Ahmed body

**Validation Metrics**:
```
Mean Absolute Error:
MAE = (1/N)Σ|φ_CFD - φ_exp|

Root Mean Square Error:
RMSE = √[(1/N)Σ(φ_CFD - φ_exp)²]

Coefficient of Determination:
R² = 1 - Σ(φ_exp - φ_CFD)²/Σ(φ_exp - φ̄_exp)²
```

## Decision Trees

### Flow Regime Decision Tree

```
START: Define Problem
  │
  ├─ Calculate Re = ρVL/μ
  │
  ├─ Re < 1 (Stokes flow)
  │  └─ Use: Creeping flow equations
  │
  ├─ 1 < Re < 100 (Low Re)
  │  └─ Use: Laminar NS, full domain
  │
  ├─ Re < 2300 (Internal)
  │  └─ Use: Laminar NS
  │
  ├─ 2300 < Re < 4000 (Transitional)
  │  └─ Use: Transition models (γ-Re_θ)
  │
  ├─ Re > 4000 (Turbulent)
  │  ├─ Calculate Ma = V/c
  │  │
  │  ├─ Ma < 0.3
  │  │  └─ Incompressible RANS
  │  │
  │  ├─ 0.3 < Ma < 1
  │  │  └─ Compressible RANS
  │  │
  │  └─ Ma > 1
  │     └─ Compressible + shock capturing
  │
  └─ Unsteady features?
     ├─ YES → Calculate St = fL/V
     │  ├─ St > 0.1 → Use URANS or LES
     │  └─ St < 0.1 → RANS may suffice
     └─ NO → RANS sufficient
```

### Turbulence Model Selection Tree

```
START: Turbulent Flow Confirmed
  │
  ├─ Computational Budget?
  │
  ├─ VERY HIGH (Research)
  │  └─ Use DNS (if Re allows)
  │
  ├─ HIGH
  │  ├─ Massively separated?
  │  │  └─ YES → LES (WALE, Smagorinsky)
  │  │
  │  └─ NO → URANS with fine mesh
  │
  └─ MODERATE (Engineering)
     │
     ├─ Flow Type?
     │
     ├─ EXTERNAL (Aero/Hydro)
     │  ├─ Pressure gradient?
     │  │  ├─ ADVERSE → k-ω SST
     │  │  └─ FAVORABLE → Spalart-Allmaras
     │  │
     │  └─ Separation likely? → k-ω SST
     │
     ├─ INTERNAL (Pipes, Channels)
     │  ├─ Simple geometry → k-ε Standard
     │  ├─ Complex, bends → k-ε Realizable
     │  └─ Heat transfer → k-ω SST
     │
     ├─ ROTATING/SWIRLING
     │  └─ RSM (Reynolds Stress Model)
     │
     └─ TRANSITIONAL
        └─ γ-Re_θ Transition SST
```

### Mesh Strategy Decision Tree

```
START: Domain Defined
  │
  ├─ Geometry Complexity?
  │
  ├─ SIMPLE (Rectangular, cylindrical)
  │  └─ Structured hexahedral
  │     └─ O-grid for cylinders
  │
  ├─ MODERATE
  │  └─ Hybrid mesh
  │     ├─ Prism layers at walls
  │     └─ Hex/tet in core
  │
  └─ COMPLEX (CAD geometry)
     └─ Unstructured
        ├─ Prism layers at walls (always)
        └─ Tet or poly in core

  │
  ├─ Boundary Layer Resolution?
  │
  ├─ WALL FUNCTIONS (y⁺ = 30-300)
  │  ├─ First cell: y⁺ ≈ 50
  │  ├─ Prism layers: 3-5
  │  └─ Growth rate: 1.2-1.3
  │
  └─ NEAR-WALL RESOLUTION (y⁺ < 1)
     ├─ First cell: y⁺ < 1
     ├─ Prism layers: 10-20
     ├─ Growth rate: 1.1-1.2
     └─ At least 10 cells in boundary layer

  │
  └─ Refinement Zones?
     ├─ Wake regions: 2-3× density
     ├─ Shear layers: graded refinement
     ├─ Recirculation: fine mesh
     └─ Free stream: coarse OK
```

## Common Pitfalls and Solutions

### Problem: Solution Not Converging

**Diagnosis**:
- Residuals oscillating
- Force coefficients fluctuating
- Solution diverging

**Solutions**:
1. Reduce under-relaxation factors
2. Improve mesh quality (aspect ratio, skewness)
3. Use first-order schemes initially
4. Check boundary condition consistency
5. Reduce time step (unsteady)
6. Switch to coupled solver

### Problem: Unphysical Results

**Diagnosis**:
- Negative pressures (absolute)
- Reverse flow at inlet
- Mass imbalance > 1%
- Unrealistic separation

**Solutions**:
1. Verify boundary conditions
2. Check material properties
3. Ensure proper reference values
4. Verify coordinate system
5. Check for mesh errors
6. Review turbulence model selection

### Problem: High y⁺ Values

**Diagnosis**:
- y⁺ > 5 when using low-Re model
- y⁺ < 30 when using wall functions

**Solutions**:
1. Recalculate first cell height
2. Adjust inflation parameters
3. Use enhanced wall treatment
4. Switch between wall function/resolved approaches

## Best Practices Summary

1. **Always start with problem physics**:
   - Calculate dimensionless numbers
   - Identify dominant phenomena
   - Select appropriate simplifications

2. **Build mesh systematically**:
   - Start coarse, refine iteratively
   - Prioritize boundary layer resolution
   - Perform grid independence study

3. **Start simple, add complexity**:
   - Begin with steady RANS
   - Use first-order upwind for stability
   - Progress to higher-order schemes
   - Add unsteadiness if needed

4. **Monitor everything**:
   - Residuals
   - Force coefficients
   - Mass flow balance
   - y⁺ distribution

5. **Validate rigorously**:
   - Compare with experiments/theory
   - Check physical realizability
   - Perform uncertainty quantification
   - Document assumptions and limitations

6. **Report uncertainty**:
   - Grid convergence index
   - Comparison with validation data
   - Sensitivity to boundary conditions
   - Model form uncertainty

## References and Resources

**Foundational Texts**:
- "Computational Fluid Dynamics" by Anderson
- "An Introduction to Computational Fluid Dynamics" by Versteeg & Malalasekera
- "Turbulent Flows" by Pope

**Validation Databases**:
- NASA Turbulence Modeling Resource
- ERCOFTAC Database
- NPARC Alliance Validation Archive

**Best Practice Guidelines**:
- AIAA Guide for Verification and Validation in CFD
- ASME V&V Standards
- ERCOFTAC Best Practice Guidelines
