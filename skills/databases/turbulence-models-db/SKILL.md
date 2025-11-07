---
name: turbulence-models-db
description: "Select and configure turbulence models (k-epsilon, k-omega SST, LES) for CFD"
category: databases
domain: fluids
complexity: advanced
dependencies: []
---

# Turbulence Models Database

A comprehensive reference for selecting and configuring turbulence models in computational fluid dynamics (CFD) simulations.

## Overview of Turbulence Modeling

Turbulence is a chaotic, three-dimensional, time-dependent flow phenomenon characterized by random fluctuations in velocity, pressure, and other flow quantities. Direct numerical simulation (DNS) of turbulent flows is computationally prohibitive for most engineering applications, necessitating turbulence modeling approaches.

### Modeling Approaches Hierarchy

1. **DNS (Direct Numerical Simulation)**: Resolves all turbulent scales, no modeling
2. **LES (Large Eddy Simulation)**: Resolves large scales, models small scales
3. **RANS (Reynolds-Averaged Navier-Stokes)**: Models all turbulent scales
4. **Laminar**: No turbulence modeling

## RANS Turbulence Models

RANS models solve time-averaged equations and model turbulent fluctuations using the Reynolds stress concept. They are the most widely used in industrial CFD due to computational efficiency.

### k-ε Models (k-epsilon)

The k-ε family models turbulent kinetic energy (k) and its dissipation rate (ε).

#### Standard k-ε

**Characteristics:**
- Two-equation model
- Robust and widely validated
- Good for free shear flows and fully turbulent flows
- Poor for flows with strong adverse pressure gradients
- Not suitable for low-Reynolds number flows without modifications

**Best Applications:**
- Fully turbulent flows
- Free shear layers, mixing layers, jets
- Flow in ducts and channels (far from walls)
- Industrial flows with high Reynolds numbers

**Limitations:**
- Overpredicts separation
- Poor near-wall performance without wall functions
- Inaccurate for swirling flows
- Stagnation point anomaly

**Wall Treatment:**
- Requires y+ > 30 (typically 30-300) with wall functions
- Not suitable for wall-resolved simulations

#### RNG k-ε

**Characteristics:**
- Derived using Renormalization Group theory
- Improved performance for swirling flows and streamline curvature
- Better handles low-Reynolds number effects
- Modified ε equation improves accuracy for rapidly strained flows

**Best Applications:**
- Flows with strong streamline curvature
- Swirling and rotating flows
- Transitional flows (with enhanced wall treatment)
- Separated flows (better than standard k-ε)

**Improvements over Standard k-ε:**
- Additional term in ε equation for rapid strain
- Modified turbulent viscosity formula
- Better prediction of near-wall flows

**Wall Treatment:**
- Can use wall functions (y+ > 30)
- Enhanced wall treatment allows y+ ≈ 1

#### Realizable k-ε

**Characteristics:**
- Ensures mathematical realizability constraints
- Variable Cμ coefficient
- Improved prediction of spreading rate for planar and round jets
- Better performance for rotating flows and boundary layers under strong adverse pressure gradients

**Best Applications:**
- Flows with rotation and recirculation
- Boundary layers with strong pressure gradients
- Separated flows
- Jets and mixing layers

**Advantages:**
- More accurate for complex flows than standard k-ε
- Superior prediction of jet spreading rates
- Better captures effects of streamline curvature

**Wall Treatment:**
- Standard wall functions (y+ > 30)
- Enhanced wall treatment available (y+ ≈ 1)

### k-ω Models (k-omega)

The k-ω family models turbulent kinetic energy (k) and specific dissipation rate (ω).

#### Standard k-ω (Wilcox)

**Characteristics:**
- Two-equation model
- Superior near-wall treatment without wall functions
- Accurate for adverse pressure gradients
- Sensitive to freestream values of ω
- Good for transitional flows

**Best Applications:**
- Low-Reynolds number flows
- Transitional flows
- Flows with adverse pressure gradients
- Aerodynamic flows (airfoils, wings)
- Wall-bounded flows

**Limitations:**
- Highly sensitive to freestream ω values
- Less accurate in free shear flows compared to k-ε
- Can be numerically stiff

**Wall Treatment:**
- Integrates to the wall (y+ ≈ 1 required)
- No wall functions needed for near-wall region

#### k-ω SST (Shear Stress Transport)

**Characteristics:**
- Blends k-ω near walls with k-ε in freestream
- Insensitive to freestream values
- Accounts for transport of turbulent shear stress
- Modified turbulent viscosity formulation
- Industry standard for aerodynamics

**Best Applications:**
- Aerodynamic flows (external aerodynamics)
- Flows with adverse pressure gradients and separation
- Transonic flows
- Heat transfer problems
- Turbomachinery

**Advantages:**
- Combines strengths of k-ω (near-wall) and k-ε (far-field)
- Accurate separation prediction
- Not sensitive to freestream turbulence values
- Robust and reliable

**Limitations:**
- Requires fine near-wall mesh (y+ ≈ 1)
- More computationally expensive than standard models
- Can underpredict separation in some cases

**Wall Treatment:**
- Designed for low-Reynolds number (y+ ≈ 1)
- Automatic wall functions available for coarse meshes
- Best results with wall-resolved mesh

### Spalart-Allmaras

**Characteristics:**
- One-equation model (solves for modified turbulent viscosity)
- Designed for aerodynamic flows
- Low computational cost
- Good for wall-bounded flows
- Limited for free shear flows and decaying turbulence

**Best Applications:**
- Aerospace applications
- External aerodynamics (airfoils, wings, fuselages)
- Mild separation and attached flows
- Transonic flows

**Advantages:**
- Computationally efficient (one equation)
- Robust and stable
- Good near-wall behavior
- Well-suited for structured meshes

**Limitations:**
- Not suitable for complex flows with multiple physics
- Limited accuracy for free shear flows
- Not ideal for internal flows
- Poor for flows with large separation regions

**Wall Treatment:**
- Designed for low-Reynolds number (y+ ≈ 1)
- Can use wall functions for coarser meshes

### Reynolds Stress Models (RSM)

**Characteristics:**
- Seven-equation model (6 Reynolds stresses + ε or ω)
- Solves transport equations for each Reynolds stress component
- Accounts for anisotropy of turbulence
- Most complex RANS approach

**Best Applications:**
- Highly swirling flows
- Flows with strong streamline curvature
- Complex 3D flows
- Rotating flows and cyclone separators
- Flows where turbulence anisotropy is critical

**Advantages:**
- Most accurate RANS model for complex flows
- Captures turbulence anisotropy
- No isotropic eddy viscosity assumption

**Limitations:**
- Most computationally expensive RANS model
- Convergence can be challenging
- Requires very good mesh quality
- More sensitive to numerical settings

## Large Eddy Simulation (LES)

### Overview

LES resolves large turbulent eddies directly while modeling small-scale (subgrid-scale) turbulence. Provides time-accurate flow structures.

**Characteristics:**
- Spatially filtered Navier-Stokes equations
- Resolves energy-containing eddies
- Models universal small-scale turbulence
- Requires 3D time-dependent simulation

**Mesh Requirements:**
- Very fine mesh (Δx, Δy, Δz ≈ local turbulent length scale)
- Isotropic or near-isotropic cells in turbulent regions
- y+ < 1 for wall-resolved LES
- Wall-modeled LES: y+ can be 30-100

**Computational Cost:**
- 10-100x more expensive than RANS
- Scales as Re^(9/4) for channel flows
- Requires long simulation times for statistical convergence

### Subgrid-Scale Models

#### Smagorinsky-Lilly
- Classic algebraic model
- Cs ≈ 0.1-0.2 (model constant)
- Overly dissipative near walls

#### Dynamic Smagorinsky
- Computes Cs dynamically
- More accurate than standard Smagorinsky
- Self-adapting to flow conditions

#### WALE (Wall-Adapting Local Eddy-viscosity)
- Better near-wall behavior
- Returns correct y³ scaling near walls
- No dynamic procedure needed

#### Kinetic Energy Subgrid-Scale
- One-equation model for subgrid kinetic energy
- More accurate but more expensive

### Applications of LES

**Best suited for:**
- Acoustics (noise prediction)
- Combustion and reacting flows
- Complex unsteady flows
- Flows with large-scale instabilities
- Vortex shedding and wake flows
- Mixing problems

**Not recommended for:**
- Steady-state problems
- High-Reynolds number wall-bounded flows (prohibitive cost)
- Industrial simulations with limited resources

## Hybrid RANS-LES Methods

### DES (Detached Eddy Simulation)
- RANS near walls, LES in separated regions
- Good for massively separated flows
- More affordable than pure LES

### DDES (Delayed DES)
- Improved shielding of boundary layer
- Prevents premature switch to LES mode

### SDES (Shielded DES)
- Further improvements to RANS-LES interface
- Better suited for attached flows

### SAS (Scale-Adaptive Simulation)
- RANS-based but resolves large unsteady structures
- Automatic adjustment to resolved scales

## Model Selection Criteria

### Flow Type Classification

#### Internal Flows
**Examples:** Pipes, ducts, channels, valves, pumps
**Recommended models:**
- k-ε Realizable (general purpose)
- k-ω SST (with heat transfer or separation)
- Standard k-ε (simple fully turbulent)

#### External Flows
**Examples:** Airfoils, vehicles, buildings, external aerodynamics
**Recommended models:**
- k-ω SST (industry standard)
- Spalart-Allmaras (aerospace)
- Realizable k-ε (blunt bodies)

#### Free Shear Flows
**Examples:** Jets, wakes, mixing layers
**Recommended models:**
- Realizable k-ε
- Standard k-ε
- RNG k-ε

#### Separated Flows
**Examples:** Flow over backward-facing step, airfoil stall
**Recommended models:**
- k-ω SST (best for mild-moderate separation)
- DES/DDES (massive separation)
- LES (if resources available)

#### Rotating/Swirling Flows
**Examples:** Turbomachinery, cyclones, swirl burners
**Recommended models:**
- RNG k-ε
- RSM
- k-ω SST

### Reynolds Number Considerations

**Low Re (Re < 10⁴):**
- Low-Re k-ω or k-ω SST
- Spalart-Allmaras
- May need transitional models

**Moderate Re (10⁴ < Re < 10⁶):**
- Most RANS models applicable
- k-ω SST for aerodynamics
- Realizable k-ε for internal flows

**High Re (Re > 10⁶):**
- Standard k-ε (with wall functions)
- k-ω SST
- Realizable k-ε

**Very High Re (Re > 10⁷):**
- Wall function approaches necessary
- Standard k-ε
- Realizable k-ε

### Mesh Requirements and y+ Values

#### Wall Functions Approach
**y+ range:** 30 < y+ < 300 (ideally 30-100)
**Models:**
- Standard k-ε
- Realizable k-ε
- RNG k-ε (with standard wall functions)

**Advantages:**
- Coarser mesh acceptable
- Lower computational cost
- Suitable for high-Re flows

**Limitations:**
- Less accurate near-wall gradients
- Not suitable for low-Re or transitional flows
- Poor for flows with separation or reattachment

#### Wall-Resolved Approach
**y+ range:** y+ ≈ 1 (first cell)
**Models:**
- k-ω SST
- Standard k-ω
- Spalart-Allmaras
- Low-Re k-ε variants

**Requirements:**
- Very fine near-wall mesh
- At least 10-15 cells in boundary layer
- y+ < 1 for first cell
- Growth ratio ≤ 1.2 near wall

**Advantages:**
- Accurate near-wall resolution
- Captures boundary layer accurately
- Better for heat transfer
- Handles low-Re and transitional flows

**Limitations:**
- High cell count
- Increased computational cost
- Mesh generation more complex

#### Enhanced Wall Treatment
**y+ range:** y+ < 5 or 30 < y+ < 300 (adaptive)
**Models:**
- Realizable k-ε with EWT
- RNG k-ε with EWT

**Advantages:**
- Flexibility in mesh resolution
- Blends wall functions and low-Re formulation
- Handles variable y+ in domain

#### y+ Guidelines by Application

| Application | Target y+ | Model Recommendation |
|-------------|-----------|---------------------|
| External aerodynamics | y+ ≈ 1 | k-ω SST |
| Heat transfer | y+ < 1 | k-ω SST, Low-Re |
| Turbomachinery | y+ ≈ 1-2 | k-ω SST |
| Internal flows (simple) | y+ = 30-100 | Realizable k-ε |
| Separation prediction | y+ ≈ 1 | k-ω SST |
| High-speed flows | y+ ≈ 1 | k-ω SST, SA |
| LES wall-resolved | y+ < 1 | LES with WALE/Dynamic |
| LES wall-modeled | y+ = 30-100 | WMLES |

### Computational Cost Comparison

**Relative cost (normalized to standard k-ε = 1):**

| Model | Relative Cost | Memory | Convergence |
|-------|---------------|---------|-------------|
| Spalart-Allmaras | 0.8 | Low | Good |
| Standard k-ε | 1.0 | Low | Excellent |
| RNG k-ε | 1.1 | Low | Good |
| Realizable k-ε | 1.1 | Low | Good |
| Standard k-ω | 1.2 | Low | Fair |
| k-ω SST | 1.3 | Low | Good |
| RSM | 2.0-2.5 | Medium | Fair-Poor |
| DES/DDES | 5-20 | High | Fair |
| LES | 50-500 | Very High | N/A (time-accurate) |

## Model Constants and Parameters

### Standard k-ε Constants
- Cμ = 0.09
- C1ε = 1.44
- C2ε = 1.92
- σk = 1.0
- σε = 1.3

### RNG k-ε Constants
- Cμ = 0.0845
- C1ε = 1.42
- C2ε = 1.68
- σk = 0.7179
- σε = 0.7179
- η0 = 4.38
- β = 0.012

### Realizable k-ε Constants
- C1ε = 1.44
- C2 = 1.9
- σk = 1.0
- σε = 1.2
- Cμ = variable (function of strain rate and rotation)

### Standard k-ω Constants
- α = 5/9
- β = 0.075
- β* = 0.09
- σk = 2.0
- σω = 2.0

### k-ω SST Constants
**k-ω inner:**
- α1 = 5/9
- β1 = 0.075
- σk1 = 2.0
- σω1 = 2.0

**k-ε outer (transformed):**
- α2 = 0.44
- β2 = 0.0828
- σk2 = 1.0
- σω2 = 1.168

**Other:**
- β* = 0.09
- a1 = 0.31 (SST limiter constant)

### Spalart-Allmaras Constants
- cb1 = 0.1355
- cb2 = 0.622
- σ = 2/3
- κ = 0.41 (von Karman constant)
- cw1 = cb1/κ² + (1 + cb2)/σ
- cw2 = 0.3
- cw3 = 2.0
- cv1 = 7.1

## Wall Functions vs Wall-Resolved

### Standard Wall Functions

**Theory:**
- Based on law of the wall
- Assumes equilibrium boundary layer
- Logarithmic law: u+ = (1/κ)ln(y+) + B

**Requirements:**
- 30 < y+ < 300
- Equilibrium turbulent boundary layer
- No significant pressure gradients

**When to use:**
- High-Re fully turbulent flows
- Simple geometries
- Limited computational resources
- Steady-state simulations

**Limitations:**
- Inaccurate for adverse pressure gradients
- Poor for separation and reattachment
- Not suitable for heat transfer predictions
- Fails in transitional flows

### Scalable Wall Functions

**Improvements:**
- Avoid deterioration for fine meshes
- y+ insensitive formulation
- Better for y+ < 30

### Non-Equilibrium Wall Functions

**Improvements:**
- Account for pressure gradient effects
- Improved separation prediction
- Better suited for complex flows

**When to use:**
- Flows with pressure gradients
- Separation and reattachment
- Complex geometries

### Enhanced Wall Treatment (EWT)

**Characteristics:**
- Two-layer approach
- Blends wall functions (high y+) with low-Re formulation (low y+)
- Adaptive based on local y+

**When to use:**
- Variable mesh resolution
- Uncertainty in y+ values
- Complex geometries with varying resolution

### Wall-Resolved (Low-Re)

**Requirements:**
- y+ ≈ 1 (ideally y+ < 1)
- 10-15+ cells in boundary layer
- Growth ratio ≤ 1.2 near wall
- Integration to wall (no wall functions)

**When to use:**
- Accurate heat transfer required
- Separation prediction critical
- Low-Re or transitional flows
- Aerodynamic design optimization

**Models requiring wall-resolved:**
- k-ω SST (optimal)
- Standard k-ω
- Spalart-Allmaras (optimal)
- LES

## Turbulence Model Selection Decision Tree

```
START: What is your flow problem?
│
├─ Need time-accurate unsteady structures?
│  │
│  ├─ YES: Go to LES/Hybrid
│  │  │
│  │  ├─ Can afford very fine mesh and long run time?
│  │  │  ├─ YES: LES (wall-resolved or wall-modeled)
│  │  │  └─ NO: DES/DDES (massively separated flows) or SAS
│  │  │
│  │  └─ Is flow primarily attached?
│  │     ├─ YES: URANS (k-ω SST or Realizable k-ε)
│  │     └─ NO: DES/DDES
│  │
│  └─ NO: Continue to RANS selection
│
├─ Flow type?
│  │
│  ├─ External aerodynamics (airfoils, vehicles, aircraft)
│  │  └─ k-ω SST (first choice) or Spalart-Allmaras
│  │     Required: y+ ≈ 1, wall-resolved mesh
│  │
│  ├─ Internal flows (pipes, ducts, channels)
│  │  │
│  │  ├─ Heat transfer important?
│  │  │  ├─ YES: k-ω SST (y+ ≈ 1)
│  │  │  └─ NO: Continue
│  │  │
│  │  ├─ Separation or adverse pressure gradients?
│  │  │  ├─ YES: k-ω SST (y+ ≈ 1) or Realizable k-ε (EWT)
│  │  │  └─ NO: Realizable k-ε (wall functions, y+ = 30-100)
│  │  │
│  │  └─ Simple fully turbulent?
│  │     └─ Standard k-ε (wall functions, y+ = 30-100)
│  │
│  ├─ Free shear flows (jets, wakes, mixing)
│  │  └─ Realizable k-ε or RNG k-ε
│  │
│  ├─ Rotating/swirling flows (turbomachinery, cyclones)
│  │  │
│  │  ├─ Simple rotation?
│  │  │  └─ RNG k-ε or k-ω SST
│  │  │
│  │  └─ Complex 3D rotation with strong curvature?
│  │     └─ RSM (if resources available) or RNG k-ε
│  │
│  └─ Separated flows (backward step, airfoil stall)
│     │
│     ├─ Mild-moderate separation?
│     │  └─ k-ω SST (y+ ≈ 1)
│     │
│     └─ Massive separation?
│        └─ DES/DDES or LES (if affordable)
│
├─ Can you achieve y+ ≈ 1 near walls?
│  │
│  ├─ YES: k-ω SST, Spalart-Allmaras, or Low-Re k-ε
│  │
│  └─ NO: Must use wall functions
│     └─ Realizable k-ε or RNG k-ε (with EWT if possible)
│
└─ Special considerations:
   │
   ├─ Transitional flows (low Re)? → k-ω SST + transition model
   ├─ Compressible/high-speed? → k-ω SST or Spalart-Allmaras
   ├─ Buoyancy-driven? → Realizable k-ε or k-ω SST with buoyancy terms
   ├─ Multiphase flows? → Realizable k-ε or k-ω SST
   └─ Limited resources? → Standard k-ε or Spalart-Allmaras
```

## Quick Selection Guide by Application

### Aerospace
**Model:** k-ω SST or Spalart-Allmaras
**y+:** ≈ 1
**Rationale:** Accurate prediction of boundary layers, separation, and pressure distribution

### Automotive (External)
**Model:** k-ω SST
**y+:** ≈ 1
**Rationale:** Separation prediction, drag/lift accuracy

### HVAC / Building Ventilation
**Model:** Realizable k-ε with wall functions
**y+:** 30-100
**Rationale:** Large domains, computational efficiency, adequate accuracy

### Turbomachinery
**Model:** k-ω SST
**y+:** 1-2
**Rationale:** Adverse pressure gradients, rotation, heat transfer

### Combustion
**Model:** Realizable k-ε or LES
**y+:** Depends (wall functions for RANS, y+ < 1 for LES)
**Rationale:** Mixing, turbulence-chemistry interaction

### Heat Exchangers
**Model:** k-ω SST or Realizable k-ε with EWT
**y+:** < 5 or wall-resolved
**Rationale:** Heat transfer accuracy critical

### Mixing / Chemical Reactors
**Model:** Realizable k-ε or RSM
**y+:** 30-100
**Rationale:** Capturing mixing patterns, turbulence anisotropy

### Hydraulics (Dams, Spillways)
**Model:** Realizable k-ε or k-ω SST
**y+:** Variable
**Rationale:** Free surface flows, separation, aeration

### Environmental (Atmospheric flows)
**Model:** Standard k-ε or Realizable k-ε
**y+:** Wall functions
**Rationale:** Large scale, computational cost, atmospheric boundary layer

## Best Practices

### General Guidelines

1. **Start simple:** Begin with simpler models (k-ε, k-ω SST) before trying complex models
2. **Validate mesh:** Ensure y+ values are appropriate for chosen model
3. **Check mesh quality:** Turbulence models are sensitive to mesh quality
4. **Monitor residuals:** Turbulence equations often converge slower than momentum
5. **Use appropriate boundary conditions:** Turbulent intensity and length scale matter
6. **Verify sensitivity:** Test mesh independence and model sensitivity

### Boundary Conditions

**Inlet:**
- Turbulent intensity: I = 0.16 Re^(-1/8) for fully developed pipe flow
- Low turbulence: I = 0.1% - 1%
- Medium turbulence: I = 1% - 5%
- High turbulence: I = 5% - 20%
- Turbulent length scale: l ≈ 0.07 × characteristic length

**Wall:**
- No-slip condition
- Wall functions or wall-resolved (model dependent)
- Roughness can be specified (equivalent sand-grain roughness)

**Outlet:**
- Zero gradient (outflow)
- Fixed pressure

### Convergence Tips

1. **Initialize properly:** Use potential flow or previous solution
2. **Under-relax initially:** Start with URF = 0.3-0.5 for turbulence equations
3. **Gradually increase:** Increase URF as solution stabilizes
4. **Coupled solvers:** Often help with k-ω models
5. **Monitor flow features:** Check separation, reattachment, vortex shedding
6. **Residuals alone insufficient:** Verify force/flux convergence

### Common Pitfalls

1. **Incorrect y+:** Most common error; check and adjust mesh
2. **Poor mesh quality:** High skewness, aspect ratio issues
3. **Inadequate refinement:** In regions of high gradients
4. **Wrong model for application:** Using k-ε for separated flows
5. **Freestream sensitivity:** k-ω without SST correction
6. **Ignoring validation:** Always compare with experiments or benchmarks

## Summary Table: RANS Model Comparison

| Model | Equations | y+ Requirement | Best For | Avoid For | Computational Cost |
|-------|-----------|----------------|----------|-----------|-------------------|
| Standard k-ε | 2 | 30-300 | Simple internal flows, fully turbulent | Separation, low-Re, adverse pressure gradients | Low |
| RNG k-ε | 2 | 30-300 (or EWT) | Swirl, rotation, moderate separation | Simple flows (overkill) | Low-Medium |
| Realizable k-ε | 2 | 30-300 (or EWT) | Complex flows, jets, separation, general purpose | Critical aerodynamics | Low-Medium |
| Standard k-ω | 2 | ≈ 1 | Low-Re, transitional, near-wall | Freestream flows (sensitive to ω∞) | Medium |
| k-ω SST | 2 | ≈ 1 | Aerodynamics, separation, adverse pressure gradients | Simple internal flows (expensive) | Medium |
| Spalart-Allmaras | 1 | ≈ 1 | Aerospace, external aerodynamics, mild separation | Internal flows, free shear, complex physics | Low |
| RSM | 7 | ≈ 1 or 30-300 | Anisotropic turbulence, strong swirl, complex 3D | Simple flows, limited resources | High |

## References

See reference.md for detailed equations, validation cases, and academic references.
