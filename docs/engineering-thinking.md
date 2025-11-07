# Engineering Thinking: Structured Workflows for Technical Analysis

This document provides structured thinking workflows for common engineering analysis tasks. Each workflow includes decision trees, step-by-step methodologies, and actionable guidance for systematic problem-solving.

---

## 1. Fluid Dynamics Analysis Workflow

### Overview
A systematic approach to analyzing fluid flow problems using computational and analytical methods.

### Workflow Steps

#### Step 1: Problem Definition
```
┌─────────────────────────────────────────┐
│     Define Fluid Dynamics Problem      │
└─────────────────┬───────────────────────┘
                  │
                  ▼
        ┌─────────────────────┐
        │ Identify Parameters │
        └─────────┬───────────┘
                  │
    ┌─────────────┼─────────────┐
    │             │             │
    ▼             ▼             ▼
[Geometry]   [Fluid Props]  [Flow Type]
```

**Actionable Checklist:**
- [ ] Define geometry and domain boundaries
- [ ] Identify fluid properties (density, viscosity, temperature)
- [ ] Determine flow regime (laminar vs. turbulent)
- [ ] Specify Reynolds number range
- [ ] Identify compressibility effects (Mach number)
- [ ] Determine if flow is steady or unsteady
- [ ] Identify heat transfer requirements

**Key Decisions:**
- **Compressible vs. Incompressible**: Ma < 0.3 → Incompressible
- **Laminar vs. Turbulent**: Re < 2300 (pipe) → Laminar
- **Steady vs. Unsteady**: Time-dependent boundary conditions → Unsteady

#### Step 2: Governing Equations Selection

```
Flow Type Decision Tree:

                    START
                      │
                      ▼
              ┌───────────────┐
              │ Compressible? │
              └───┬───────┬───┘
                  │       │
              Yes │       │ No
                  │       │
                  ▼       ▼
         [Full Navier-   [Incompressible
          Stokes +       Navier-Stokes]
          Energy Eqn]           │
                                ▼
                        ┌───────────────┐
                        │   Turbulent?  │
                        └───┬───────┬───┘
                            │       │
                        Yes │       │ No
                            │       │
                            ▼       ▼
                     [RANS/LES/   [Direct
                      DNS]        Solution]
```

**Equation Selection Matrix:**

| Flow Characteristics | Governing Equations | Complexity |
|---------------------|--------------------|-----------:|
| Incompressible, Laminar | Navier-Stokes | Low |
| Incompressible, Turbulent | RANS (k-ε, k-ω) | Medium |
| Compressible, Subsonic | Euler + Viscous terms | Medium |
| Compressible, Supersonic | Full NS + Shocks | High |
| Multiphase | VOF/Eulerian-Eulerian | Very High |

**Actionable Steps:**
1. Calculate dimensionless numbers (Re, Ma, Fr, We)
2. Select appropriate simplifications (inviscid, potential flow, boundary layer)
3. Choose turbulence model if Re > 4000
4. Determine if thermal effects must be included
5. Select coordinate system (Cartesian, cylindrical, spherical)

#### Step 3: Boundary Conditions

**Boundary Condition Decision Matrix:**

```
Inlet BC:
    Velocity Known? ──Yes──> [Velocity Inlet]
         │
         No
         │
         ▼
    Mass Flow Known? ──Yes──> [Mass Flow Inlet]
         │
         No
         │
         ▼
    Pressure Known? ──Yes──> [Pressure Inlet]

Outlet BC:
    Pressure Known? ──Yes──> [Pressure Outlet]
         │
         No
         │
         ▼
    Flow Fully Developed? ──Yes──> [Outflow]
         │
         No
         │
         ▼
    [Zero Gradient]

Wall BC:
    Moving Wall? ──Yes──> [Moving Wall + Velocity]
         │
         No
         │
         ▼
    Heat Transfer? ──Yes──> [Isothermal/Heat Flux]
         │
         No
         │
         ▼
    [No-Slip Wall]
```

**Checklist:**
- [ ] Define all inlet conditions (velocity/pressure/temperature)
- [ ] Define all outlet conditions
- [ ] Specify wall boundary conditions (no-slip, slip, moving)
- [ ] Define thermal boundary conditions
- [ ] Specify symmetry planes if applicable
- [ ] Define periodic boundaries if applicable

#### Step 4: Discretization and Meshing

**Meshing Strategy Decision Tree:**

```
                    START
                      │
                      ▼
            ┌─────────────────┐
            │ Complex Geometry?│
            └────┬────────┬────┘
                 │        │
             Yes │        │ No
                 │        │
                 ▼        ▼
         [Unstructured   [Structured
          Tetrahedral]    Hexahedral]
                 │            │
                 ▼            ▼
          ┌──────────┐  ┌──────────┐
          │High Grad?│  │High Grad?│
          └─┬────────┘  └─┬────────┘
            │ Yes          │ Yes
            ▼              ▼
        [Boundary      [Boundary
         Layer Mesh]    Layer Mesh]
```

**Mesh Quality Criteria:**

| Criterion | Target Value | Critical? |
|-----------|--------------|-----------|
| Skewness | < 0.85 | Yes |
| Aspect Ratio | < 100 | Yes |
| Orthogonality | > 0.15 | Yes |
| y+ (wall) | < 1 (LES/DNS) or 30-300 (wall fn) | Yes |

**Actionable Steps:**
1. **Generate base mesh**
   - Start with coarse mesh for initial solution
   - Use structured mesh where possible for accuracy
   - Apply unstructured mesh for complex geometries

2. **Refine critical regions**
   - Boundary layers (inflation layers)
   - High velocity gradient regions
   - Shock waves or flow separation zones
   - Wake regions

3. **Mesh independence study**
   - Generate 3 meshes: coarse, medium, fine
   - Compare results (velocity, pressure, forces)
   - Ensure < 5% difference between medium and fine
   - Select mesh based on accuracy/cost tradeoff

4. **Quality checks**
   - Check minimum cell volume > 0
   - Verify y+ values at walls
   - Check aspect ratios in critical regions
   - Smooth transitions (growth rate < 1.2)

#### Step 5: Solving

**Solver Selection:**

```
                    Steady State?
                      │
            ┌─────────┴─────────┐
            │                   │
           Yes                 No
            │                   │
            ▼                   ▼
    ┌──────────────┐    ┌──────────────┐
    │Pressure-Based│    │  Time Method │
    │   (SIMPLE)   │    └──────┬───────┘
    └──────────────┘           │
                               ▼
                    ┌──────────────────────┐
                    │ Explicit vs Implicit?│
                    └──┬────────────────┬──┘
                       │                │
                  Explicit          Implicit
                       │                │
                       ▼                ▼
                  [Low CFL,        [Large time
                   Small Δt]        steps, stable]
```

**Convergence Monitoring:**

Residuals Checklist:
- [ ] Continuity residual < 1e-3
- [ ] Momentum residuals < 1e-3
- [ ] Energy residual < 1e-6 (if thermal)
- [ ] Turbulence residuals < 1e-3
- [ ] Monitor convergence history (should decrease monotonically)

Monitored Quantities:
- [ ] Drag/Lift coefficients (stable within 1%)
- [ ] Mass flow balance (inlet = outlet within 0.5%)
- [ ] Pressure drop (stable)
- [ ] Velocity at key locations (stable)

**Actionable Solving Procedure:**

1. **Initialization**
   - Use hybrid initialization or FMG
   - Set appropriate initial conditions
   - Check for negative volumes or temperatures

2. **Start with simplified physics**
   - Run inviscid solution first
   - Add viscosity
   - Enable turbulence models
   - Add heat transfer if needed

3. **Under-relaxation factors**
   - Start conservative (0.3-0.5)
   - Gradually increase as solution stabilizes
   - Monitor residuals for oscillations

4. **Adaptive time stepping (transient)**
   - Start with small time steps
   - Use CFL number for guidance
   - Increase as solution stabilizes

#### Step 6: Post-Processing and Verification

**Verification Checklist:**

```
Post-Processing Workflow:

    Results Ready
         │
         ▼
    ┌─────────────────┐
    │ Mass Balance OK?│
    └────┬──────────┬─┘
         │          │
        Yes         No ──> Refine mesh/solution
         │
         ▼
    ┌─────────────────┐
    │ Energy Balance? │
    └────┬──────────┬─┘
         │          │
        Yes         No ──> Check BC/models
         │
         ▼
    ┌─────────────────┐
    │Grid Independent?│
    └────┬──────────┬─┘
         │          │
        Yes         No ──> Refine mesh
         │
         ▼
    ┌─────────────────┐
    │Physical Results?│
    └────┬──────────┬─┘
         │          │
        Yes         No ──> Review setup
         │
         ▼
    [VALIDATED]
```

**Actionable Verification Steps:**

1. **Physical Checks**
   - [ ] Velocity vectors aligned with expected flow
   - [ ] No reverse flow at inlets
   - [ ] Boundary layer profiles reasonable
   - [ ] Pressure gradients physical
   - [ ] Temperature fields physical

2. **Conservation Checks**
   - [ ] Mass: ∑(ṁ_in) - ∑(ṁ_out) < 0.5%
   - [ ] Momentum: Forces balance
   - [ ] Energy: Heat in = Heat out + Work (± 1%)

3. **Comparison with Known Solutions**
   - [ ] Compare with analytical solutions (if available)
   - [ ] Benchmark against experimental data
   - [ ] Compare with published CFD results
   - [ ] Verify dimensionless parameters (Cd, Cl, Nu, etc.)

4. **Sensitivity Analysis**
   - [ ] Vary mesh density (±20%)
   - [ ] Vary time step (transient, ±50%)
   - [ ] Test different turbulence models
   - [ ] Vary inlet conditions (±10%)

**Visualization Guidelines:**
- Velocity contours and vectors
- Pressure contours
- Streamlines or pathlines
- Vorticity contours
- Turbulence quantities (k, ε, ω)
- Wall shear stress
- Heat transfer coefficient distributions

---

## 2. Pump Design Workflow

### Overview
Systematic methodology for designing centrifugal and positive displacement pumps from requirements to validation.

### Workflow Steps

#### Step 1: Requirements Gathering

**Requirements Template:**

```
┌──────────────────────────────────────────┐
│        PUMP SPECIFICATION SHEET          │
├──────────────────────────────────────────┤
│ Operating Conditions:                    │
│  • Flow Rate (Q):         ___ m³/h       │
│  • Head (H):              ___ m          │
│  • NPSH Available:        ___ m          │
│  • Fluid:                 ________       │
│  • Temperature:           ___ °C         │
│  • Density:               ___ kg/m³      │
│  • Viscosity:             ___ cP         │
│  • Vapor Pressure:        ___ kPa        │
│                                          │
│ Performance Requirements:                │
│  • Efficiency Target:     ___ %          │
│  • Operating Range:       ___-___ m³/h   │
│  • Max Power:             ___ kW         │
│  • Speed Range:           ___ RPM        │
│                                          │
│ Special Requirements:                    │
│  • Solids content:        ___ %          │
│  • Corrosivity:           Yes/No         │
│  • Abrasiveness:          Yes/No         │
│  • Seal type required:    ________       │
│  • Installation:          H / V          │
│  • Standards:             ________       │
└──────────────────────────────────────────┘
```

**Checklist:**
- [ ] Define duty point (Q, H)
- [ ] Specify fluid properties
- [ ] Determine operating range (min/nom/max flow)
- [ ] Calculate NPSHa
- [ ] Identify environmental constraints
- [ ] Specify materials requirements
- [ ] Define reliability/maintenance requirements
- [ ] Determine cost constraints
- [ ] Identify applicable standards (API, ANSI, ISO)

#### Step 2: Pump Type Selection

**Selection Decision Tree:**

```
                        START
                          │
                          ▼
                 ┌────────────────┐
                 │  Flow < 50 m³/h│
                 │  AND           │
                 │  H > 100 m?    │
                 └───┬────────┬───┘
                     │        │
                    Yes       No
                     │        │
                     ▼        ▼
              [Consider PD  ┌──────────────┐
               Pump]        │Viscosity High?│
                            │(> 100 cP)    │
                            └──┬────────┬──┘
                               │        │
                              Yes       No
                               │        │
                               ▼        ▼
                        [PD Pump:    ┌──────────┐
                         Gear/Screw/ │Head > 50m│
                         Lobe]       │ AND      │
                                     │Q<500m³/h?│
                                     └─┬──────┬─┘
                                       │      │
                                      Yes     No
                                       │      │
                                       ▼      ▼
                                  [Centrifugal: [Centrifugal:
                                   Multistage]   Single Stage]
```

**Pump Type Selection Matrix:**

| Pump Type | Flow Range | Head Range | Viscosity | Solids | Efficiency |
|-----------|------------|------------|-----------|--------|------------|
| Centrifugal - Single Stage | 10-10000 m³/h | 5-150 m | < 200 cP | < 5% | 70-85% |
| Centrifugal - Multistage | 10-1000 m³/h | 50-2000 m | < 100 cP | < 1% | 65-80% |
| Axial Flow | 500-50000 m³/h | 2-20 m | < 50 cP | < 2% | 75-90% |
| Mixed Flow | 100-10000 m³/h | 10-80 m | < 100 cP | < 3% | 70-85% |
| Gear Pump | 0.1-500 m³/h | 10-200 m | 1-100000 cP | < 1% | 70-80% |
| Screw Pump | 1-1000 m³/h | 10-100 m | 1-500000 cP | < 30% | 50-75% |
| Diaphragm | 0.01-100 m³/h | 10-200 m | Any | < 70% | 30-60% |

**Specific Speed Calculation:**

```
Ns = N × √Q / H^0.75

Where:
  N = Rotational speed (RPM)
  Q = Flow rate (m³/s)
  H = Head (m)

Selection Guide:
  Ns < 20    → Radial flow (low flow, high head)
  20 < Ns < 40 → Francis type
  40 < Ns < 80 → Mixed flow
  Ns > 80     → Axial flow (high flow, low head)
```

#### Step 3: Preliminary Sizing

**For Centrifugal Pumps:**

**Actionable Calculation Sequence:**

1. **Calculate Specific Speed**
   ```
   Ns = N × √Q / H^0.75
   ```
   - Choose optimal speed based on motor availability
   - Target Ns for best efficiency (typically 30-60 for centrifugal)

2. **Estimate Impeller Diameter**
   ```
   D2 = 60 × √(2gH) / (π × N × ψ)

   Where:
     ψ = Head coefficient (0.4-1.2 for centrifugal)
     g = 9.81 m/s²
   ```

3. **Calculate Impeller Width**
   ```
   b2 = Q / (π × D2 × Vm2)

   Where:
     Vm2 = Meridional velocity at outlet (2-5 m/s typical)
   ```

4. **Determine Blade Angles**
   ```
   β2 = 20°-40° (backward curved blades for efficiency)
   β1 = 90° - arctan(Vm1/U1)

   Where:
     U1 = π × D1 × N / 60 (peripheral velocity at inlet)
   ```

5. **Number of Blades**
   ```
   Z = 6.5 × (D2 + D1)/(D2 - D1) × sin((β1 + β2)/2)

   Round to nearest integer: typically 5-9 blades
   ```

**Preliminary Design Checklist:**
- [ ] Specific speed within optimal range
- [ ] Impeller diameter sized
- [ ] Blade angles determined
- [ ] Number of blades selected
- [ ] Volute or diffuser type chosen
- [ ] Shaft diameter estimated (torsion + deflection)
- [ ] Bearing locations determined
- [ ] Seal arrangement selected
- [ ] NPSH required calculated

#### Step 4: Detailed Design

**Design Process Flow:**

```
Preliminary Design
       │
       ▼
┌─────────────────┐
│ Impeller Design │
└────────┬────────┘
         │
         ├──> [Hydraulic Design]
         │      • Blade profile
         │      • Meridional contour
         │      • Blade thickness distribution
         │
         ├──> [Mechanical Design]
         │      • Stress analysis
         │      • Dynamic balance
         │      • Material selection
         │
         ▼
┌─────────────────┐
│  Volute Design  │
└────────┬────────┘
         │
         ├──> [Geometry]
         │      • Cross-section area distribution
         │      • Cutwater design
         │      • Discharge nozzle
         │
         ▼
┌─────────────────┐
│  Shaft Design   │
└────────┬────────┘
         │
         ├──> [Sizing]
         │      • Torque capacity
         │      • Critical speed
         │      • Deflection limits
         │
         ▼
┌─────────────────┐
│ Bearing Select. │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   Seal Design   │
└────────┬────────┘
         │
         ▼
   [Complete Design]
```

**Impeller Detailed Design:**

1. **Meridional Contour Design**
   - [ ] Define hub and shroud profiles
   - [ ] Ensure smooth acceleration through impeller
   - [ ] Check for separation zones
   - [ ] Optimize for minimum losses

2. **Blade Profile Design**
   - [ ] Define blade angles from inlet to outlet
   - [ ] Use conformal mapping or streamline curvature method
   - [ ] Ensure smooth velocity distribution
   - [ ] Check blade loading (Δp across blade)

3. **Blade Thickness**
   - [ ] Structural requirements (stress)
   - [ ] Manufacturing constraints
   - [ ] Blockage effects (< 10% of passage area)
   - [ ] Leading edge: 2-5 mm
   - [ ] Trailing edge: 1-3 mm

4. **Mechanical Considerations**
   - [ ] Stress analysis (centrifugal, bending)
   - [ ] Critical speed > 1.3 × operating speed
   - [ ] Balance to ISO Grade G6.3 or better
   - [ ] Fatigue life > 10^7 cycles

**Volute Design:**

```
Area Distribution Formula:
  A(θ) = (Q × r)/(8.5 × H^0.5) × (θ/360°)

  Where:
    θ = Angle from cutwater (degrees)
    r = Radius to centroid of section
    Q = Flow rate (m³/s)
    H = Head (m)
```

**Checklist:**
- [ ] Cross-sectional area increases proportionally with angle
- [ ] Cutwater positioned correctly (typically 1.05-1.1 × D2)
- [ ] Cutwater angle optimized (3-10° typically)
- [ ] Discharge nozzle sized (velocity 2-4 m/s)
- [ ] Transition smooth to discharge flange

**Shaft Design:**

1. **Diameter from Torsion**
   ```
   d_shaft = ∛(16 × T / (π × τ_allow))

   Where:
     T = Torque = Power / (2π × N/60)
     τ_allow = Allowable shear stress (typically 40-60 MPa for steel)
   ```

2. **Critical Speed Check**
   ```
   N_critical = (60/2π) × √(k/m)

   Operating speed < 0.75 × N_critical (rigid shaft)
   OR
   Operating speed > 1.3 × N_critical (flexible shaft)
   ```

3. **Deflection Check**
   - [ ] Deflection at impeller < 0.05 mm
   - [ ] Deflection at seal < 0.02 mm
   - [ ] Slope at bearings within limits

**Seal Selection:**

| Seal Type | Pressure Limit | Speed Limit | Application |
|-----------|----------------|-------------|-------------|
| Gland Packing | 10 bar | 15 m/s | General service, low cost |
| Mechanical Seal - Single | 25 bar | 25 m/s | Most common, clean fluids |
| Mechanical Seal - Double | 100 bar | 30 m/s | Hazardous fluids |
| Magnetic Drive | 40 bar | 20 m/s | Zero leakage required |

#### Step 5: Performance Prediction

**Analytical Method (for centrifugal pumps):**

```
Performance Prediction Sequence:

1. Theoretical Head:
   H_th = (U2 × Cu2 - U1 × Cu1)/g

   Where:
     U = Peripheral velocity
     Cu = Tangential component of absolute velocity

2. Hydraulic Efficiency:
   η_h = H_actual / H_th ≈ 0.85-0.95

3. Volumetric Efficiency:
   η_v = Q_actual / Q_theoretical ≈ 0.92-0.98

   Account for:
     • Leakage past wear rings
     • Leakage to balance device
     • Disk friction recirculation

4. Mechanical Efficiency:
   η_m = P_hydraulic / P_shaft ≈ 0.96-0.99

   Account for:
     • Bearing friction
     • Seal friction
     • Disk friction

5. Overall Efficiency:
   η_overall = η_h × η_v × η_m
```

**Loss Estimation:**

| Loss Component | Typical Magnitude | Calculation Method |
|----------------|-------------------|-------------------|
| Incidence loss | 0.5-2% H | (C1 - C1,design)² / 2g |
| Friction loss | 2-5% H | f × L/Dh × V²/2g |
| Diffusion loss | 1-3% H | (1-Cp) × V²/2g |
| Recirculation | 1-4% H | (Q_leak/Q) × H |
| Shock loss | 0-10% H | Off-design operation |
| Leakage loss | 2-5% Q | Δp × clearance³ |

**Performance Curve Generation:**

```
Operating Points to Calculate:

Flow Rate (% of Design):
  0%, 25%, 50%, 75%, 100%, 110%, 120%

For each point, calculate:
  1. Velocity triangles
  2. Theoretical head
  3. Losses (hydraulic, volumetric, mechanical)
  4. Actual head
  5. Power consumption
  6. Efficiency
  7. NPSH required

Plot curves:
  • H vs Q
  • P vs Q
  • η vs Q
  • NPSHr vs Q
```

**Actionable Prediction Steps:**

- [ ] Calculate theoretical performance at design point
- [ ] Estimate losses at design point
- [ ] Predict efficiency
- [ ] Generate full performance curve (0-120% Q)
- [ ] Calculate power consumption across range
- [ ] Predict NPSH required
- [ ] Estimate operating range (min continuous flow to max flow)
- [ ] Identify best efficiency point (BEP)
- [ ] Check cavitation limits

#### Step 6: Testing and Validation

**Test Planning:**

```
Test Procedure Workflow:

    Pump Assembly
         │
         ▼
    Pre-Test Checks
    • Alignment
    • Rotation direction
    • Lubrication
    • Instrumentation
         │
         ▼
    Hydrostatic Test
    • 1.5 × design pressure
    • Hold 10 minutes
    • Check for leaks
         │
         ▼
    Performance Test
    • Multiple flow points
    • Measure H, Q, P, T
    • Vibration monitoring
         │
         ▼
    NPSH Test
    • Reduce suction pressure
    • Identify 3% head drop
         │
         ▼
    Endurance Test (optional)
    • 24-72 hour run
    • Monitor wear
         │
         ▼
    Analysis & Report
```

**Required Instrumentation:**

| Parameter | Instrument | Accuracy | Location |
|-----------|------------|----------|----------|
| Flow Rate | Magnetic flowmeter | ±0.5% | Discharge line |
| Suction Pressure | Pressure transducer | ±0.1% | Suction flange |
| Discharge Pressure | Pressure transducer | ±0.1% | Discharge flange |
| Speed | Tachometer | ±0.1% | Shaft |
| Power | Wattmeter | ±0.5% | Motor |
| Temperature | RTD | ±0.5°C | Suction/discharge |
| Vibration | Accelerometer | ±5% | Bearing housings |

**Performance Test Procedure:**

1. **Preparation**
   - [ ] Verify instrumentation calibration
   - [ ] Fill system, vent air
   - [ ] Check baseline vibration
   - [ ] Verify motor rotation
   - [ ] Set speed to design value

2. **Data Collection**
   - [ ] Test at 7 flow points: 0%, 40%, 60%, 80%, 100%, 110%, 120% Q
   - [ ] Allow stabilization (5 minutes minimum)
   - [ ] Record all parameters simultaneously
   - [ ] Take multiple readings (minimum 3)
   - [ ] Monitor vibration continuously

3. **Calculations**
   ```
   Head: H = (Pd - Ps)/(ρ × g) + (Vd² - Vs²)/(2g) + Zd - Zs

   Power: P_hydraulic = ρ × g × Q × H

   Efficiency: η = P_hydraulic / P_shaft

   Specific Speed: Ns = N × √Q / H^0.75
   ```

4. **Acceptance Criteria**
   - [ ] Head at design flow: +0% to -3% of predicted
   - [ ] Efficiency: > 95% of predicted or within 2%
   - [ ] Power: ±5% of predicted
   - [ ] NPSH required: < predicted value
   - [ ] Vibration: < ISO 10816 limits
   - [ ] Temperature rise: < 15°C

**Validation Checklist:**

- [ ] Performance meets specification
- [ ] Efficiency within tolerance
- [ ] No cavitation at operating point
- [ ] Vibration acceptable
- [ ] No excessive temperature rise
- [ ] No leakage at seals
- [ ] Bearing temperatures normal
- [ ] Noise levels acceptable
- [ ] Operating range verified
- [ ] Transient behavior acceptable (startup/shutdown)

**Common Issues and Resolutions:**

| Issue | Possible Cause | Resolution |
|-------|----------------|------------|
| Low head | Wrong impeller diameter | Measure, replace if needed |
| Low efficiency | Excessive clearances | Check wear rings, adjust |
| High vibration | Imbalance, misalignment | Re-balance, re-align |
| Cavitation | NPSH insufficient | Increase suction pressure |
| High power | Wrong impeller trim | Verify diameter |
| Noisy operation | Cavitation, bearing issues | Check NPSHa, inspect bearings |

---

## 3. Thermodynamic Analysis

### Overview
Systematic approach to analyzing thermodynamic systems, cycles, and processes.

### Workflow Steps

#### Step 1: System Definition

**System Boundary Definition:**

```
         ┌─────────────────────────────────┐
         │ Define System Boundary          │
         └──────────────┬──────────────────┘
                        │
         ┌──────────────┴──────────────┐
         │                             │
         ▼                             ▼
  [Closed System]              [Open System]
  • Fixed mass                 • Mass flow
  • No mass transfer           • Energy transfer
  • Energy transfer            • Work & heat
         │                             │
         ▼                             ▼
  [Examples:]                  [Examples:]
  • Piston-cylinder            • Turbines
  • Rigid tanks                • Compressors
  • Sealed vessels             • Heat exchangers
                               • Nozzles
```

**System Classification Decision Tree:**

```
                    Is mass crossing boundary?
                              │
                    ┌─────────┴─────────┐
                    │                   │
                   Yes                  No
                    │                   │
                    ▼                   ▼
            [Open/Control       [Closed System]
             Volume]                    │
                    │                   ▼
                    │            Is volume changing?
                    │                   │
                    │              ┌────┴────┐
                    │              │         │
                    │             Yes        No
                    │              │         │
                    │              ▼         ▼
                    │         [Moving     [Rigid]
                    │          Boundary]
                    │
                    ▼
            Is it steady state?
                    │
              ┌─────┴─────┐
              │           │
             Yes          No
              │           │
              ▼           ▼
         [Steady     [Unsteady
          Flow]       Flow]
```

**Actionable Checklist:**

- [ ] Sketch system and draw boundary
- [ ] Identify all mass flows (in/out)
- [ ] Identify all energy transfers (Q, W)
- [ ] Determine if system is open or closed
- [ ] Determine if process is steady or unsteady
- [ ] List all relevant properties (P, T, v, h, s, etc.)
- [ ] Identify constraints or assumptions
- [ ] Specify reference state if needed

#### Step 2: State Point Analysis

**Property Determination Workflow:**

```
State Point Analysis:

For each state point:
    │
    ▼
┌────────────────────────┐
│ How many properties    │
│ are known?             │
└───┬────────────────┬───┘
    │                │
   Two            < Two
    │                │
    ▼                ▼
[State Fully    [Need more
 Defined]        information]
    │                │
    ▼                ▼
Use property    Apply constraints:
tables/EOS      • Process type
    │           • Conservation laws
    │           • Equilibrium conditions
    │                │
    │                ▼
    │           [Iterate until
    │            two properties
    └─────┬──── per state known]
          │
          ▼
    [Calculate all
     properties]
```

**Property Evaluation Methods:**

| Substance State | Method | Tools |
|----------------|--------|-------|
| Ideal Gas | P×v = R×T | Gas tables, equations |
| Compressed Liquid | h ≈ hf(T), s ≈ sf(T) | Liquid tables |
| Saturated Mixture | Quality x = (v-vf)/(vg-vf) | Saturation tables |
| Superheated Vapor | Interpolation | Superheat tables |
| Real Gas | P-v-T relations | EOS (van der Waals, RK, PR) |

**State Point Documentation Template:**

```
State 1:
  Known: P1 = ___ kPa, T1 = ___ K
  Phase: ________________
  Properties:
    • v1 = ___ m³/kg
    • u1 = ___ kJ/kg
    • h1 = ___ kJ/kg
    • s1 = ___ kJ/(kg·K)
  Source: [Table/EOS/Calculation]

[Repeat for all state points]
```

**Actionable Steps:**

1. **Identify number of state points**
   - Inlet states (all inlet streams)
   - Outlet states (all outlet streams)
   - Intermediate states (within system)

2. **For each state point:**
   - [ ] List known properties (from problem statement)
   - [ ] Determine phase (compressed liquid, saturated, superheated)
   - [ ] Use appropriate property source (tables, EOS)
   - [ ] Calculate all properties (P, T, v, u, h, s)
   - [ ] Document assumptions

3. **Check consistency:**
   - [ ] All states physically realizable
   - [ ] Properties internally consistent
   - [ ] Phase appropriate for conditions

#### Step 3: Cycle Calculations

**General Cycle Analysis Procedure:**

```
Cycle Analysis Workflow:

    Define Cycle
         │
         ▼
    Identify States
    (Label 1, 2, 3, ...)
         │
         ▼
    For Each Process:
         │
    ┌────┴────┐
    │         │
    ▼         ▼
[Heat     [Work
 Transfer] Transfer]
    │         │
    └────┬────┘
         │
         ▼
    Apply Conservation:
    • Energy
    • Mass
    • Entropy
         │
         ▼
    Calculate Net Work
    W_net = W_out - W_in
         │
         ▼
    Calculate Net Heat
    Q_net = Q_in - Q_out
         │
         ▼
    Verify: W_net = Q_net
         │
         ▼
    Calculate Efficiency
```

**Process-Specific Work/Heat Calculations:**

| Process Type | Work | Heat |
|--------------|------|------|
| Isobaric (P=const) | W = P(V2-V1) | Q = mcp(T2-T1) |
| Isochoric (V=const) | W = 0 | Q = mcv(T2-T1) |
| Isothermal (T=const) | W = mRT×ln(V2/V1) | Q = W |
| Adiabatic (Q=0) | W = m(u1-u2) | Q = 0 |
| Polytropic | W = (P2V2-P1V1)/(1-n) | Q = W + m(u2-u1) |

**Power Cycle Analysis (e.g., Rankine, Brayton):**

```
Rankine Cycle Example:

States: 1→2→3→4→1

State 1: Pump inlet (saturated liquid)
  Known: P1, x1=0
  Find: h1, s1

Process 1→2: Pump (isentropic)
  s2 = s1
  P2 = boiler pressure
  Find: h2
  W_pump = h2 - h1

State 3: Turbine inlet (superheated)
  Known: P3=P2, T3
  Find: h3, s3

Process 3→4: Turbine (isentropic)
  s4 = s3
  P4 = P1
  Find: h4, x4
  W_turbine = h3 - h4

Process 4→1: Condenser
  Q_out = h4 - h1

Process 2→3: Boiler
  Q_in = h3 - h2

Net Work:
  W_net = W_turbine - W_pump

Efficiency:
  η_th = W_net / Q_in
```

**Refrigeration Cycle Analysis (e.g., Vapor-Compression):**

```
States: 1→2→3→4→1

Process 1→2: Compressor
  W_comp = h2 - h1

Process 2→3: Condenser
  Q_out = h2 - h3

Process 3→4: Expansion valve
  h4 = h3 (isenthalpic)

Process 4→1: Evaporator
  Q_in = h1 - h4

COP:
  COP_cooling = Q_in / W_comp
  COP_heating = Q_out / W_comp
```

**Actionable Cycle Analysis Steps:**

1. **Sketch T-s or P-h diagram**
   - [ ] Plot all state points
   - [ ] Draw process lines
   - [ ] Label processes
   - [ ] Shade area = work (for T-s diagram)

2. **Calculate properties at all states**
   - [ ] Use state point analysis procedure
   - [ ] Document all values

3. **Calculate work for each component**
   - [ ] Turbines: W_out = ṁ(h_in - h_out)
   - [ ] Compressors/Pumps: W_in = ṁ(h_out - h_in)

4. **Calculate heat transfer for each component**
   - [ ] Boilers/Heaters: Q_in = ṁ(h_out - h_in)
   - [ ] Condensers/Coolers: Q_out = ṁ(h_in - h_out)

5. **Calculate net work and heat**
   - [ ] W_net = ∑W_out - ∑W_in
   - [ ] Q_net = ∑Q_in - ∑Q_out

6. **Verify energy balance**
   - [ ] W_net = Q_net (within rounding error)

#### Step 4: Efficiency Evaluation

**Efficiency Definitions:**

```
Efficiency Types:

                Thermal Systems
                       │
         ┌─────────────┼─────────────┐
         │             │             │
         ▼             ▼             ▼
    [Power Cycle] [Refrigeration] [Heat Engine]
         │             │             │
         ▼             ▼             ▼
    η = W_net/Q_in  COP = Q_L/W   η = W/Q_H
```

**Efficiency Evaluation Decision Tree:**

```
                What type of system?
                       │
        ┌──────────────┼──────────────┐
        │              │              │
        ▼              ▼              ▼
    [Power       [Refrigerator/  [Heat Engine
     Plant]       Heat Pump]      Component]
        │              │              │
        ▼              ▼              ▼
   η_th = W_net    COP_R = Q_L    η_turbine =
        /Q_in           /W         (h1-h2_actual)
        │                          /(h1-h2_ideal)
        │              COP_HP =
        │               Q_H/W
        │
        ▼
    Compare to Carnot:
    η_Carnot = 1 - T_L/T_H
```

**Component Efficiency Definitions:**

| Component | Isentropic Efficiency | Equation |
|-----------|----------------------|----------|
| Turbine | η_T = W_actual / W_isentropic | (h1-h2a)/(h1-h2s) |
| Compressor | η_C = W_isentropic / W_actual | (h2s-h1)/(h2a-h1) |
| Pump | η_P = W_isentropic / W_actual | (h2s-h1)/(h2a-h1) |
| Nozzle | η_N = KE_actual / KE_isentropic | V2a² / V2s² |

**Exergy (Second Law) Efficiency:**

```
Exergy Analysis:

Exergy of a stream:
  ex = (h - h0) - T0(s - s0) + V²/2 + gz

Exergy destroyed:
  I = T0 × S_gen

Second-law efficiency:
  η_II = Exergy_output / Exergy_input
       = 1 - I / Exergy_input
```

**Actionable Efficiency Evaluation:**

1. **Calculate First-Law Efficiency**
   - [ ] Identify energy input (fuel, electricity, heat)
   - [ ] Identify desired energy output (work, cooling, heating)
   - [ ] Calculate η = Output / Input

2. **Calculate Carnot Efficiency (if applicable)**
   - [ ] Identify high temperature reservoir (TH)
   - [ ] Identify low temperature reservoir (TL)
   - [ ] Calculate η_Carnot = 1 - TL/TH

3. **Calculate Component Efficiencies**
   - [ ] For each component, calculate actual vs. ideal
   - [ ] Use isentropic efficiency definitions

4. **Perform Exergy Analysis**
   - [ ] Define dead state (T0, P0)
   - [ ] Calculate exergy at all states
   - [ ] Calculate exergy destruction in each component
   - [ ] Identify major sources of irreversibility

5. **Optimization Opportunities**
   - [ ] Rank components by exergy destruction
   - [ ] Identify processes far from reversible
   - [ ] Calculate potential efficiency improvements

**Efficiency Improvement Checklist:**

- [ ] Reduce temperature difference in heat exchangers
- [ ] Improve component isentropic efficiencies
- [ ] Reduce pressure drops (friction losses)
- [ ] Implement regeneration/recuperation
- [ ] Optimize operating pressures and temperatures
- [ ] Minimize throttling losses
- [ ] Improve insulation (reduce heat loss)

---

## 4. Structural Analysis (FEA)

### Overview
Systematic workflow for finite element structural analysis from problem setup to results validation.

### Workflow Steps

#### Step 1: Load Definition

**Load Classification:**

```
Load Types:

        Loads
          │
    ┌─────┴─────┐
    │           │
    ▼           ▼
[Mechanical] [Thermal]
    │           │
    ├─> Point Force
    ├─> Distributed Load
    ├─> Pressure
    ├─> Body Force (gravity)
    ├─> Torque
    ├─> Displacement (enforced)
    │
    └─> [Static vs Dynamic]
              │
        ┌─────┴─────┐
        │           │
        ▼           ▼
    [Static]   [Dynamic]
                    │
              ┌─────┴─────┐
              │           │
              ▼           ▼
        [Harmonic]  [Transient]
        [Impact]    [Seismic]
```

**Load Definition Checklist:**

- [ ] Identify all applied forces (magnitude, direction, location)
- [ ] Identify pressure loads (uniform, varying)
- [ ] Include body forces (gravity, centrifugal, magnetic)
- [ ] Define thermal loads (temperature distribution)
- [ ] Specify constraints/boundary conditions
- [ ] Determine if loads are static or dynamic
- [ ] Check load combinations (if multiple load cases)
- [ ] Apply safety factors per design codes

**Load Case Documentation:**

```
Load Case 1: [Description]
─────────────────────────────
Mechanical Loads:
  • Force 1: F = ___ N, Direction: ___, Location: ___
  • Pressure: P = ___ Pa, Surface: ___
  • Gravity: g = 9.81 m/s², Direction: -Z

Thermal Loads:
  • Temperature: T = ___ °C
  • Heat flux: q = ___ W/m²

Constraints:
  • Fixed support: Face/Edge ___
  • Frictionless support: Face ___
  • Symmetry: Plane ___

Load Factor: ___
```

**Actionable Steps:**

1. **Create free body diagram**
   - [ ] Sketch component/assembly
   - [ ] Draw all external loads
   - [ ] Show all reactions
   - [ ] Label coordinate system

2. **Quantify loads**
   - [ ] Calculate magnitudes
   - [ ] Define directions (unit vectors)
   - [ ] Specify points of application
   - [ ] Document assumptions

3. **Define load cases**
   - [ ] Normal operating condition
   - [ ] Maximum load condition
   - [ ] Startup/shutdown transients
   - [ ] Emergency/fault conditions
   - [ ] Fatigue loading (if cyclic)

#### Step 2: Material Selection

**Material Selection Decision Tree:**

```
              Material Selection
                      │
        ┌─────────────┼─────────────┐
        │             │             │
        ▼             ▼             ▼
    [Metals]    [Polymers]    [Composites]
        │             │             │
    ┌───┴───┐        │         [Direction-
    │       │        │          dependent
    ▼       ▼        ▼          properties]
 [Steel] [Aluminum] [Plastic]
    │
    ├─> Low Carbon (structural)
    ├─> Alloy Steel (high strength)
    ├─> Stainless (corrosion)
    └─> Tool Steel (hardness)
```

**Material Property Requirements:**

| Property | Symbol | Unit | Required For |
|----------|--------|------|--------------|
| Young's Modulus | E | GPa | Stiffness, deflection |
| Poisson's Ratio | ν | - | Multi-axial stress |
| Yield Strength | σy | MPa | Plastic deformation |
| Ultimate Strength | σu | MPa | Fracture |
| Density | ρ | kg/m³ | Mass, dynamic analysis |
| Thermal Expansion | α | 1/K | Thermal stress |
| Thermal Conductivity | k | W/(m·K) | Heat transfer |
| Specific Heat | cp | J/(kg·K) | Thermal analysis |

**Material Behavior Models:**

```
Material Model Selection:

    Is deformation < 0.2% strain?
              │
        ┌─────┴─────┐
        │           │
       Yes          No
        │           │
        ▼           ▼
    [Linear     Is material ductile?
     Elastic]         │
                 ┌────┴────┐
                 │         │
                Yes        No
                 │         │
                 ▼         ▼
            [Elasto-   [Brittle
             plastic]   Fracture]
```

**Material Model Checklist:**

- [ ] Linear elastic (most common, small deformations)
- [ ] Elasto-plastic (large deformations, yielding)
- [ ] Hyperelastic (rubber, large elastic strains)
- [ ] Viscoelastic (time-dependent, polymers)
- [ ] Orthotropic (composites, wood)
- [ ] Temperature-dependent properties
- [ ] Creep model (high temperature, long duration)
- [ ] Fatigue data (S-N curve for cyclic loading)

**Actionable Material Selection:**

1. **Determine requirements**
   - [ ] Strength requirements (yield, ultimate)
   - [ ] Stiffness requirements (deflection limits)
   - [ ] Environmental (temperature, corrosion)
   - [ ] Manufacturing constraints
   - [ ] Cost constraints

2. **Select material**
   - [ ] Choose material class
   - [ ] Select specific alloy/grade
   - [ ] Verify availability
   - [ ] Check design codes/standards

3. **Gather properties**
   - [ ] Get material data sheet
   - [ ] Input into FEA software
   - [ ] Apply temperature dependencies if needed
   - [ ] Document source of data

#### Step 3: Meshing Strategy

**Mesh Type Selection:**

```
Meshing Decision Tree:

        Geometry Complexity?
                │
        ┌───────┴───────┐
        │               │
    Simple          Complex
        │               │
        ▼               ▼
[Structured    [Unstructured Mesh]
 Hexahedral]            │
        │               ├─> Tetrahedral (3D)
        │               ├─> Triangular (2D)
        │               └─> Hybrid
        │
        └─> [Sweep/Map mesh]
```

**Element Type Selection:**

| Geometry | Element Type | Nodes | Use Case |
|----------|--------------|-------|----------|
| 3D Solid | Hex8 / Tet4 | 8/4 | General 3D |
| 3D Solid | Hex20 / Tet10 | 20/10 | Higher accuracy |
| Thin structures | Shell (quad/tri) | 4/3 | Sheet metal, pressure vessels |
| Long slender | Beam | 2 | Frames, trusses |
| Axisymmetric | 2D axisymmetric | 4 | Rotational symmetry |

**Mesh Quality Metrics:**

| Metric | Target Range | Critical? |
|--------|--------------|-----------|
| Aspect Ratio | < 10 (3:1 ideal) | Yes |
| Skewness | < 0.8 | Yes |
| Jacobian Ratio | > 0.6 | Yes |
| Warpage (shells) | < 10° | Yes |
| Element Quality | > 0.3 | Moderate |

**Meshing Strategy Workflow:**

```
Meshing Process:

    Start
      │
      ▼
    ┌─────────────────┐
    │ Simplify Geometry│
    │ • Remove fillets│
    │ • Remove holes  │
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │ Apply Symmetry  │
    │ • 1/4 model     │
    │ • 1/2 model     │
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │ Choose Element  │
    │ Type            │
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │ Set Global Size │
    │ (coarse)        │
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │ Refine Critical │
    │ Regions         │
    │ • Stress conc.  │
    │ • Contact areas │
    │ • Crack tips    │
    └────────┬────────┘
             │
             ▼
    ┌─────────────────┐
    │ Check Quality   │
    └────────┬────────┘
             │
        ┌────┴────┐
        │         │
       Pass      Fail
        │         │
        ▼         ▼
    [Solve]  [Remesh]
```

**Critical Region Refinement:**

```
Regions Requiring Fine Mesh:

High Stress Concentration:
  • Fillets (r/t < 0.1)
  • Holes
  • Re-entrant corners
  • Load application points
  Target: 8-12 elements across radius

Contact Areas:
  • At least 3 elements across contact width
  • Matched mesh if possible

Crack/Fracture:
  • Fine mesh at crack tip
  • Mesh size << crack length

Thin Features:
  • At least 3 elements through thickness
```

**Actionable Meshing Steps:**

1. **Geometry preparation**
   - [ ] Remove unnecessary details (< 5% of overall size)
   - [ ] Simplify small features
   - [ ] Apply symmetry if applicable
   - [ ] Partition for structured meshing

2. **Initial mesh generation**
   - [ ] Select element type (solid, shell, beam)
   - [ ] Set global element size (start coarse)
   - [ ] Generate mesh
   - [ ] Check element count (target 10k-1M for statics)

3. **Mesh refinement**
   - [ ] Identify high-stress regions
   - [ ] Apply local refinement (sizing controls)
   - [ ] Ensure smooth transitions (size ratio < 2:1)
   - [ ] Refine contact regions

4. **Quality check**
   - [ ] Check element quality metrics
   - [ ] Fix highly distorted elements
   - [ ] Verify aspect ratios in critical regions
   - [ ] Check for free edges/overlaps

5. **Mesh convergence study**
   - [ ] Generate 3 meshes (coarse, medium, fine)
   - [ ] Compare critical results (stress, deflection)
   - [ ] Ensure < 5% change between medium and fine
   - [ ] Select appropriate mesh density

#### Step 4: Solver Selection

**Analysis Type Selection:**

```
Solver Selection Tree:

    Time-dependent loads?
              │
        ┌─────┴─────┐
        │           │
       No          Yes
        │           │
        ▼           ▼
    [Static]   [Dynamic]
        │           │
        │      ┌────┴────┐
        │      │         │
        │      ▼         ▼
        │  [Modal]  [Transient]
        │      │         │
        │      │    ┌────┴────┐
        │      │    │         │
        │      │    ▼         ▼
        │      │ [Explicit] [Implicit]
        │      │
        │      └──> [Harmonic]
        │
        ├─> [Linear]
        └─> [Nonlinear]
                │
           ┌────┴────┐
           │         │
           ▼         ▼
      [Geometric] [Material]
      [Contact]
```

**Solver Type Comparison:**

| Solver | Applications | Time Step | Accuracy | Cost |
|--------|--------------|-----------|----------|------|
| Static Linear | Small deformations, linear materials | N/A | High | Low |
| Static Nonlinear | Large deformations, plasticity, contact | N/A | High | Medium |
| Modal | Natural frequencies, mode shapes | N/A | High | Low |
| Harmonic | Frequency response | N/A | High | Medium |
| Transient Implicit | General dynamics, long duration | Large Δt | High | High |
| Transient Explicit | Impact, crash, high-speed | Small Δt | Medium | Very High |

**Nonlinearity Decision:**

```
Is analysis nonlinear?

Material Nonlinearity:
  • Plasticity? (stress > yield)
  • Hyperelasticity? (rubber)
  • Creep? (time-dependent)

Geometric Nonlinearity:
  • Large deformations? (strain > 5%)
  • Large rotations?
  • Buckling?

Boundary Nonlinearity:
  • Contact? (changing contact area)
  • Friction?
  • Gaps?

If ANY are Yes → [Nonlinear Solver Required]
```

**Solver Settings Checklist:**

**For Static Analysis:**
- [ ] Choose solver (direct vs. iterative)
- [ ] Set convergence criteria (force, displacement)
- [ ] Enable/disable large deflection
- [ ] Set maximum iterations (nonlinear)
- [ ] Configure load stepping (nonlinear)

**For Dynamic Analysis:**
- [ ] Set time step size (Δt < T_min / 20)
- [ ] Set total simulation time
- [ ] Choose integration method (Newmark, HHT)
- [ ] Set damping (Rayleigh coefficients)
- [ ] Configure output frequency

**For Contact:**
- [ ] Define contact pairs
- [ ] Set contact algorithm (penalty, Lagrange, augmented)
- [ ] Define friction coefficient
- [ ] Set normal stiffness
- [ ] Enable/disable bonded, no-separation, frictionless

**Actionable Solver Selection:**

1. **Classify problem type**
   - [ ] Static or dynamic?
   - [ ] Linear or nonlinear?
   - [ ] Identify nonlinearities (material, geometric, contact)

2. **Select solver type**
   - [ ] Match solver to problem type
   - [ ] Check software capabilities
   - [ ] Consider computational resources

3. **Configure solver settings**
   - [ ] Set appropriate tolerances
   - [ ] Configure load stepping (if nonlinear)
   - [ ] Set output requests
   - [ ] Enable necessary physics (thermal coupling, etc.)

4. **Verify setup**
   - [ ] Check units consistency
   - [ ] Verify loads applied correctly
   - [ ] Verify constraints defined
   - [ ] Check for rigid body modes (should be zero for constrained model)

#### Step 5: Results Interpretation

**Results Verification Workflow:**

```
Results Available
       │
       ▼
┌──────────────────┐
│ Sanity Checks    │
│ • Deformation    │
│ • Reaction forces│
│ • Energy balance │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Convergence OK?  │
└────┬────────┬────┘
     │        │
    Yes       No → [Adjust settings,
     │              refine mesh,
     │              check BCs]
     ▼
┌──────────────────┐
│ Stress Evaluation│
│ • von Mises      │
│ • Principal      │
│ • Safety factor  │
└────────┬─────────┘
         │
         ▼
┌──────────────────┐
│ Design Criteria? │
└────┬────────┬────┘
     │        │
    Pass     Fail → [Redesign]
     │
     ▼
┌──────────────────┐
│ Report Results   │
└──────────────────┘
```

**Sanity Checks:**

```
Quick Verification Steps:

1. Deformed Shape:
   Does it make physical sense?
   ┌─────────────────────────┐
   │ ✓ Deflects in expected  │
   │   direction             │
   │ ✓ Maximum at expected   │
   │   location              │
   │ ✗ Penetration through   │
   │   other parts           │
   │ ✗ Unrestrained rigid    │
   │   body motion           │
   └─────────────────────────┘

2. Reaction Forces:
   Sum of reactions = Applied loads?
   ∑Rx = ∑Fx
   ∑Ry = ∑Fy
   ∑Rz = ∑Fz
   ∑M = ∑Applied moments

3. Energy Balance:
   Strain Energy + Kinetic Energy = Work Done
```

**Stress Evaluation:**

| Stress Type | Use Case | Failure Criterion |
|-------------|----------|-------------------|
| von Mises | Ductile materials | σvm < σy / SF |
| Maximum Principal | Brittle materials (tension) | σ1 < σult / SF |
| Minimum Principal | Brittle materials (compression) | |σ3| < σult_comp / SF |
| Tresca | Conservative ductile | σ1 - σ3 < σy / SF |
| Maximum Shear | Shear-critical | τmax < τy / SF |

**Safety Factor Evaluation:**

```
Safety Factor Calculation:

For Ductile Materials:
  SF = σ_yield / σ_von_Mises

For Brittle Materials:
  SF = σ_ultimate / σ_max_principal

Typical Target Safety Factors:
  • Static, well-defined loads: SF ≥ 1.5
  • Dynamic, uncertain loads:  SF ≥ 2.0-3.0
  • Pressure vessels (ASME):   SF ≥ 4.0
  • Fatigue:                   Infinite life design
```

**Critical Results to Extract:**

- [ ] **Displacement**
  - Maximum displacement location and magnitude
  - Compare to allowable deflection limits

- [ ] **Stress**
  - Maximum von Mises stress
  - Location of maximum stress
  - Stress distribution in critical regions
  - Compare to yield/ultimate strength

- [ ] **Safety Factor**
  - Minimum safety factor
  - Distribution of safety factor
  - Identify regions with SF < target

- [ ] **Reaction Forces**
  - Total reactions at constraints
  - Compare to applied loads (equilibrium check)

- [ ] **Strain**
  - Maximum strain
  - Check if in linear elastic range (< 0.2%)

- [ ] **Contact Results** (if applicable)
  - Contact pressure
  - Contact area
  - Sliding distance

**Actionable Interpretation Steps:**

1. **Visual inspection**
   - [ ] View deformed shape (with undeformed overlay)
   - [ ] Check deformation direction makes sense
   - [ ] Identify maximum displacement location
   - [ ] Look for unexpected behavior

2. **Convergence verification**
   - [ ] Check convergence plot (should be decreasing)
   - [ ] Verify residuals < tolerance
   - [ ] Check final load step completed
   - [ ] Review warnings/errors in log

3. **Equilibrium check**
   - [ ] Sum reaction forces in X, Y, Z
   - [ ] Compare to applied forces
   - [ ] Should balance within 1%

4. **Stress analysis**
   - [ ] Plot von Mises stress
   - [ ] Identify maximum stress location
   - [ ] Check if at expected location (geometry, loading)
   - [ ] Extract stress values at critical points
   - [ ] Compare to allowable stress

5. **Mesh sensitivity check**
   - [ ] Verify mesh refined in high-stress regions
   - [ ] Check that max stress not at single node (singularity)
   - [ ] If questionable, perform mesh refinement study

6. **Design verification**
   - [ ] Calculate safety factors
   - [ ] Compare to design requirements
   - [ ] Check displacement limits
   - [ ] Verify against design codes/standards

7. **Document results**
   - [ ] Create images (deformation, stress contours)
   - [ ] Extract data tables
   - [ ] Document maximum values
   - [ ] Note locations of critical stresses
   - [ ] Compare to hand calculations (if available)

**Common Pitfalls:**

| Issue | Symptom | Solution |
|-------|---------|----------|
| Singularity | Very high localized stress | Ignore if at sharp corner, refine mesh, or use stress averaging |
| Rigid body mode | Zero/very small stiffness | Add proper constraints |
| Mesh too coarse | Stress changes drastically with refinement | Refine mesh, perform convergence study |
| Units error | Results off by orders of magnitude | Check unit consistency |
| Contact not converging | No solution or oscillations | Adjust contact settings, initial penetration |

---

## 5. Design Optimization

### Overview
Structured approach to engineering design optimization using systematic methodologies.

### Workflow Steps

#### Step 1: Objective Function Definition

**Objective Function Classification:**

```
Optimization Goal:

        Objective Function
                │
    ┌───────────┼───────────┐
    │           │           │
    ▼           ▼           ▼
[Minimize]  [Maximize]  [Target]
    │           │           │
    │           │           │
    ▼           ▼           ▼
Examples:   Examples:   Examples:
• Mass      • Efficiency • Frequency
• Cost      • Strength   • Temperature
• Stress    • Flow rate  • Dimension
• Energy    • Power
```

**Single vs. Multi-Objective:**

```
Number of Objectives:

    Single Objective
          │
          ▼
    [Standard
     Optimization]
    min/max f(x)

    Multi-Objective
          │
          ▼
    ┌──────────────┐
    │ Trade-offs   │
    │ Pareto Front │
    └──────────────┘
          │
    ┌─────┴─────┐
    │           │
    ▼           ▼
[Weighted   [Pareto
 Sum]        Optimization]
```

**Objective Function Checklist:**

- [ ] **Identify primary objective**
  - What are you trying to achieve?
  - Minimize: mass, cost, stress, drag, energy
  - Maximize: efficiency, stiffness, strength, throughput
  - Target: specific frequency, temperature, dimension

- [ ] **Quantify objective mathematically**
  - Express as function of design variables: f(x1, x2, ..., xn)
  - Ensure objective is calculable/measurable
  - Define units clearly

- [ ] **Handle multiple objectives**
  - If multiple objectives, rank by importance
  - Consider weighted sum: f = w1×f1 + w2×f2 + ...
  - Or use multi-objective method (Pareto optimization)

- [ ] **Verify objective makes engineering sense**
  - Is minimizing mass actually what you want?
  - Or do you want maximize stiffness-to-weight ratio?
  - Consider competing objectives (cost vs. performance)

**Objective Function Examples:**

| Engineering Problem | Objective Function | Type |
|---------------------|-------------------|------|
| Structural design | Minimize: Mass = ρ×V | Single, minimize |
| Heat exchanger | Maximize: Q/(ΔP×Cost) | Single, maximize |
| Pump impeller | Multi: Max η, Min cost | Multi-objective |
| Beam design | Minimize: Deflection subject to mass < M_max | Constrained |
| Aircraft wing | Minimize: Drag subject to Lift = W | Constrained |

**Actionable Objective Definition:**

1. **State design goal clearly**
   ```
   Example: "Design a bracket that minimizes mass
            while supporting 1000 N load with
            maximum stress < 200 MPa and
            deflection < 2 mm"
   ```

2. **Formulate mathematical objective**
   ```
   Minimize: f(x) = Mass(x) = ρ × Volume(x)

   Where x = [thickness, width, height, ...]
   ```

3. **Identify if multi-objective**
   ```
   If multiple goals:
   • Option 1: Weighted sum
     f = w1×Mass + w2×Cost

   • Option 2: Pareto optimization
     Find solutions minimizing both Mass and Cost
   ```

4. **Document objective function**
   - [ ] Write mathematical expression
   - [ ] Define all variables
   - [ ] Specify units
   - [ ] State minimize/maximize/target

#### Step 2: Constraints

**Constraint Classification:**

```
Constraint Types:

        Constraints
              │
    ┌─────────┼─────────┐
    │         │         │
    ▼         ▼         ▼
[Equality] [Inequality] [Bounds]
    │          │           │
    │          │           │
   h(x)=0    g(x)≤0     L≤x≤U
    │          │           │
    ▼          ▼           ▼
Examples:  Examples:   Examples:
• Balance  • Stress    • 5≤t≤20 mm
• Flow     • Deflection• 0≤θ≤90°
• Energy   • Temp      • D>0
```

**Constraint Types:**

| Constraint Type | Form | Example | Handling |
|----------------|------|---------|----------|
| Equality | h(x) = 0 | Mass_total = 50 kg | Lagrange multipliers, penalty |
| Inequality | g(x) ≤ 0 | σ_max ≤ 200 MPa | Active set, interior point |
| Box constraints | L ≤ x ≤ U | 5 mm ≤ t ≤ 20 mm | Simple bounds |
| Linear | Ax ≤ b | x1 + 2x2 ≤ 100 | Linear programming |
| Nonlinear | g(x) ≤ 0 | Volume(x) ≤ V_max | Nonlinear programming |

**Common Engineering Constraints:**

```
Structural:
  • σ_max ≤ σ_allowable
  • Deflection ≤ δ_max
  • Frequency > f_min (avoid resonance)
  • Buckling_load > SF × Applied_load

Thermal:
  • T_max ≤ T_allowable
  • Heat_flux < q_max
  • Temperature_gradient < ΔT_max

Manufacturing:
  • Thickness > t_min (manufacturing limit)
  • Feature_size > Tool_size
  • Draft_angle > α_min

Geometric:
  • Volume ≤ V_max
  • Length < L_max
  • Aspect_ratio < AR_max

Performance:
  • Efficiency > η_min
  • Flow_rate > Q_min
  • Power < P_max

Cost:
  • Cost < Budget
```

**Constraint Formulation Workflow:**

```
For each constraint:

    Identify Requirement
           │
           ▼
    Physical/Engineering Limit
    • Safety factor
    • Code requirement
    • Manufacturing limit
    • Operational requirement
           │
           ▼
    Express Mathematically
    g(x) ≤ 0 or h(x) = 0
           │
           ▼
    ┌──────────────┐
    │ Can you      │
    │ calculate it?│
    └──┬────────┬──┘
       │        │
      Yes       No
       │        │
       ▼        ▼
    [Include   [Surrogate
     directly]  model or
                simplify]
```

**Actionable Constraint Definition:**

1. **List all requirements**
   - [ ] Design code requirements (ASME, AISC, API, etc.)
   - [ ] Safety factors
   - [ ] Manufacturing constraints
   - [ ] Assembly constraints
   - [ ] Operational limits
   - [ ] Cost/budget limits

2. **Formulate mathematically**

   Example: Stress constraint
   ```
   Requirement: "Maximum stress must not exceed
                 allowable stress with SF=2"

   Constraint: σ_max(x) ≤ σ_yield / 2

   Or in standard form: g(x) = σ_max(x) - σ_yield/2 ≤ 0
   ```

3. **Check if constraint is active**
   - Active constraint: Optimal solution lies on constraint boundary
   - Inactive constraint: Constraint not limiting (slack)

   Example: If optimized design has σ_max = 100 MPa
            and σ_allow = 200 MPa, constraint is inactive

4. **Ensure feasible region exists**
   ```
   Check that constraints don't conflict:

   Example of conflict:
   • Constraint 1: Volume < 100 cm³
   • Constraint 2: Volume > 200 cm³
   → No feasible solution!
   ```

5. **Document all constraints**
   ```
   Constraint 1 (Stress):
     g1(x) = σ_max(x) / σ_allow - 1 ≤ 0

   Constraint 2 (Deflection):
     g2(x) = δ_max(x) / δ_allow - 1 ≤ 0

   Constraint 3 (Mass):
     g3(x) = Mass(x) / Mass_max - 1 ≤ 0

   Bounds:
     5 mm ≤ thickness ≤ 50 mm
     0° ≤ angle ≤ 90°
   ```

#### Step 3: Variable Selection

**Design Variable Types:**

```
Design Variables:

           Variables
                │
    ┌───────────┼───────────┐
    │           │           │
    ▼           ▼           ▼
[Continuous] [Discrete] [Integer]
    │           │           │
    │           │           │
    ▼           ▼           ▼
Examples:    Examples:   Examples:
• Dimensions • Material   • Number of
• Angles       choice       ribs
• Thickness  • Bolt size  • Teeth
• Position   • Catalog      count
              part

            ┌────┴────┐
            │         │
            ▼         ▼
        [Categorical] [Ordinal]
        Color,        Small/Med/
        Material      Large
```

**Variable Selection Criteria:**

```
Is this a good design variable?

    ┌────────────────────────┐
    │ Does it significantly  │
    │ affect objective or    │
    │ constraints?           │
    └───┬────────────────┬───┘
        │                │
       Yes               No → Don't include
        │
        ▼
    ┌────────────────────────┐
    │ Can it be changed      │
    │ independently?         │
    └───┬────────────────┬───┘
        │                │
       Yes               No → Derived variable
        │
        ▼
    ┌────────────────────────┐
    │ Is range well-defined? │
    └───┬────────────────┬───┘
        │                │
       Yes               No → Define bounds first
        │
        ▼
    [Good Design Variable]
```

**Variable Parameterization:**

| Design Feature | Poor Parameterization | Good Parameterization |
|----------------|----------------------|----------------------|
| Beam cross-section | 20 individual points | Width, height (2 vars) |
| Fillet radius | Radii at 10 locations | Single radius parameter |
| Hole pattern | Individual x,y for each | Pitch, angle, number |
| Shape | Many control points | Few key dimensions + formula |

**Actionable Variable Selection:**

1. **Identify candidate variables**
   - [ ] Dimensions (lengths, widths, thicknesses)
   - [ ] Angles (orientation, draft, chamfer)
   - [ ] Material properties (if material is variable)
   - [ ] Positions (hole locations, component placement)
   - [ ] Counts (number of features, blades, fins)
   - [ ] Shapes (parameterized via key dimensions)

2. **Sensitivity screening**
   ```
   For each candidate variable:
   1. Change variable by ±10%
   2. Calculate change in objective
   3. Sensitivity = Δf / Δx
   4. Rank variables by sensitivity
   5. Keep high-sensitivity variables
   ```

3. **Reduce variable count**
   - Combine related variables (use ratios, aspect ratios)
   - Use symmetry (reduce by half or quarter)
   - Fix low-sensitivity variables at reasonable values
   - Use design rules to link variables

   ```
   Example: Instead of t1, t2, t3, t4 (4 variables)
            Use: t_base, and ratios r1, r2, r3 (4 variables)
            Or:  t_base, uniform scaling factor (2 variables)
   ```

4. **Define bounds for each variable**
   ```
   Variable 1: Thickness (t)
     • Lower bound: 5 mm (manufacturing limit)
     • Upper bound: 50 mm (weight/cost limit)
     • Initial value: 20 mm

   Variable 2: Width (w)
     • Lower bound: 10 mm
     • Upper bound: 100 mm
     • Initial value: 50 mm
   ```

5. **Document design variables**
   ```
   Design Variable Table:

   | Var | Description | Type | Min | Max | Units | Initial |
   |-----|-------------|------|-----|-----|-------|---------|
   | x1  | Thickness   | Cont | 5   | 50  | mm    | 20      |
   | x2  | Width       | Cont | 10  | 100 | mm    | 50      |
   | x3  | Angle       | Cont | 0   | 90  | deg   | 45      |
   | x4  | Material    | Disc | -   | -   | -     | Steel   |
   | x5  | Num_ribs    | Int  | 3   | 10  | -     | 5       |
   ```

**Typical Variable Counts:**

| Problem Complexity | Number of Variables | Algorithm Suitability |
|-------------------|--------------------|-----------------------|
| Simple | 1-5 | Any method works |
| Moderate | 5-20 | Gradient-based, GA |
| Complex | 20-100 | Gradient-based preferred |
| High-dimensional | 100+ | Advanced methods, surrogate models |

#### Step 4: Optimization Algorithm Choice

**Algorithm Selection Decision Tree:**

```
                 Algorithm Selection
                         │
                         ▼
              ┌──────────────────────┐
              │ Gradient available?  │
              └──┬───────────────┬───┘
                 │               │
                Yes              No
                 │               │
                 ▼               ▼
         ┌───────────────┐  ┌──────────────┐
         │ Many variables│  │ Discrete or  │
         │ (>20)?        │  │ Non-smooth?  │
         └──┬────────┬───┘  └──┬────────┬──┘
            │        │         │        │
           Yes       No       Yes       No
            │        │         │        │
            ▼        ▼         ▼        ▼
      [Gradient-  [SQP,     [GA,      [Pattern
       based:     Interior   SA,       Search,
       BFGS]      Point]     PSO]      Simplex]
```

**Algorithm Comparison:**

| Algorithm | Type | Variables | Speed | Global? | Gradient? |
|-----------|------|-----------|-------|---------|-----------|
| Sequential Quadratic Programming (SQP) | Gradient | <100 | Fast | No | Yes |
| Interior Point | Gradient | <1000 | Fast | No | Yes |
| BFGS/L-BFGS | Gradient | <1000 | Fast | No | Yes |
| Genetic Algorithm (GA) | Evolutionary | <100 | Slow | Yes | No |
| Particle Swarm (PSO) | Swarm | <50 | Medium | Yes | No |
| Simulated Annealing (SA) | Stochastic | <50 | Slow | Yes | No |
| Nelder-Mead Simplex | Direct search | <20 | Medium | No | No |
| Pattern Search | Direct search | <50 | Medium | No | No |
| Response Surface | Surrogate | <20 | Medium | No | No |

**Algorithm Selection Guide:**

```
Problem Characteristics → Recommended Algorithm

Smooth, few variables (<10):
  → SQP, Interior Point

Smooth, many variables (10-100):
  → BFGS, L-BFGS

Smooth, very many variables (>100):
  → Conjugate Gradient, L-BFGS

Non-smooth or discrete:
  → Genetic Algorithm, Simulated Annealing

Need global optimum:
  → Multi-start + local optimizer
  → Genetic Algorithm
  → Particle Swarm

Expensive function evaluations:
  → Response Surface / Kriging
  → Bayesian Optimization
  → Few iterations: Pattern Search

Mixed integer:
  → Branch and Bound
  → Genetic Algorithm
```

**Actionable Algorithm Selection:**

1. **Characterize your problem**
   - [ ] Number of variables: ___
   - [ ] Continuous or discrete: ___
   - [ ] Gradient available: Yes / No / Finite difference
   - [ ] Function evaluation cost: Cheap (<1s) / Moderate (1s-1min) / Expensive (>1min)
   - [ ] Multiple local minima expected: Yes / No
   - [ ] Constraints: None / Simple bounds / Linear / Nonlinear

2. **Select primary algorithm**
   ```
   Decision flow:

   IF (gradient available AND variables < 100)
     → Use gradient-based (SQP, Interior Point)

   ELSE IF (discrete variables OR non-smooth)
     → Use evolutionary (GA, PSO)

   ELSE IF (expensive function evaluations)
     → Use surrogate-based optimization

   ELSE
     → Use pattern search or simplex
   ```

3. **Configure algorithm parameters**

   **For Gradient-Based (SQP):**
   - [ ] Convergence tolerance: 1e-6 (tight) to 1e-3 (loose)
   - [ ] Maximum iterations: 100-1000
   - [ ] Finite difference step (if no analytical gradient): 1e-6
   - [ ] Line search method: Backtracking, cubic interpolation

   **For Genetic Algorithm:**
   - [ ] Population size: 10-20 × number of variables
   - [ ] Generations: 50-500
   - [ ] Crossover probability: 0.7-0.9
   - [ ] Mutation probability: 0.01-0.1
   - [ ] Selection: Tournament, roulette wheel

   **For Response Surface:**
   - [ ] Sampling method: Latin Hypercube, DOE
   - [ ] Number of samples: 10-20 × number of variables
   - [ ] Surrogate type: Polynomial, Kriging, RBF
   - [ ] Infill criteria: Expected improvement

4. **Plan multi-start strategy (if seeking global)**
   ```
   Global Optimization Strategy:

   1. Run global search (GA, PSO) for N generations
   2. Take best M solutions
   3. Use each as starting point for local optimizer (SQP)
   4. Select best result

   Or:

   1. Generate N random starting points (Latin Hypercube)
   2. Run local optimizer from each
   3. Select best result
   ```

5. **Select optimization software/tools**

   | Tool | Strengths | Algorithms Available |
   |------|-----------|---------------------|
   | MATLAB Optimization Toolbox | General purpose, easy | SQP, Interior Point, GA, PSO, Pattern Search |
   | Python scipy.optimize | Free, flexible | SLSQP, COBYLA, L-BFGS-B, Differential Evolution |
   | ANSYS DesignXplorer | Integrated with FEA | Response Surface, MOGA, Screening |
   | modeFRONTIER | Multi-objective | Many algorithms, surrogate models |
   | OpenMDAO | Multidisciplinary | Multiple optimizers, gradient-based |
   | Dakota | General purpose | Many algorithms, UQ |

#### Step 5: Verification of Results

**Optimization Result Verification Workflow:**

```
Optimization Complete
         │
         ▼
┌────────────────────┐
│ Convergence Check  │
│ • Iterations       │
│ • Tolerance met?   │
└────────┬───────────┘
         │
         ▼
┌────────────────────┐
│ Constraint Check   │
│ • All satisfied?   │
│ • Active vs slack  │
└────────┬───────────┘
         │
         ▼
┌────────────────────┐
│ Optimality Check   │
│ • KKT conditions   │
│ • Sensitivity      │
└────────┬───────────┘
         │
         ▼
┌────────────────────┐
│ Physical Validity  │
│ • Makes sense?     │
│ • Manufacturable?  │
└────────┬───────────┘
         │
         ▼
┌────────────────────┐
│ Robustness Check   │
│ • Perturbation test│
│ • Tolerances       │
└────────┬───────────┘
         │
         ▼
    [Validated]
```

**Convergence Verification:**

```
Convergence Criteria:

1. Objective function:
   |f(x_k) - f(x_k-1)| < ε_f
   Typical: ε_f = 1e-6

2. Design variables:
   ||x_k - x_k-1|| < ε_x
   Typical: ε_x = 1e-6

3. Gradient (if available):
   ||∇f(x_k)|| < ε_g
   Typical: ε_g = 1e-6

4. Maximum iterations reached?
   If yes, may not be converged!
```

**Checklist:**
- [ ] Convergence tolerance met
- [ ] Not stopped due to iteration limit
- [ ] Objective function history shows plateau
- [ ] Variable changes are small
- [ ] No oscillations in late iterations

**Constraint Verification:**

```
For each constraint g_i(x*) ≤ 0:

┌────────────────────────────────┐
│ Constraint i:                  │
│   Value: g_i(x*) = ______      │
│   Status: [ ] Active (≈ 0)     │
│           [ ] Slack (< 0)      │
│           [ ] Violated (> 0)   │
└────────────────────────────────┘

Active constraints define the optimum
Violated constraints → infeasible solution!
```

**Checklist:**
- [ ] All constraints satisfied (g_i ≤ 0)
- [ ] Identify active constraints (g_i ≈ 0)
- [ ] Active constraints make engineering sense
- [ ] Check if constraint values are close to bounds
- [ ] Verify no numerical violations due to tolerance

**Optimality Conditions (KKT):**

```
For constrained optimization, at optimum x*:

1. Stationarity:
   ∇f(x*) + ∑λ_i × ∇g_i(x*) = 0

2. Primal feasibility:
   g_i(x*) ≤ 0 for all i

3. Dual feasibility:
   λ_i ≥ 0 for all i

4. Complementary slackness:
   λ_i × g_i(x*) = 0 for all i

If satisfied → Local optimum found
```

**Actionable Optimality Check:**
- [ ] Check gradient of Lagrangian is zero (if available)
- [ ] Verify Lagrange multipliers are non-negative
- [ ] Check complementary slackness (λ_i > 0 only for active constraints)

**Physical Validity Check:**

```
Engineering Sanity Checks:

Does the optimal design:
  [ ] Make physical sense?
  [ ] Have reasonable dimensions?
  [ ] Meet manufacturing constraints?
  [ ] Avoid singularities or extreme values?
  [ ] Result in expected improvement?

Red flags:
  • Variables at bounds (may need to expand range)
  • Unrealistic geometry
  • Very thin features (< manufacturing limit)
  • Extremely high aspect ratios
```

**Actionable Physical Checks:**

1. **Compare to baseline**
   ```
   Baseline design:
     • Objective: f_baseline = ___
     • Mass: ___
     • Stress: ___

   Optimized design:
     • Objective: f_opt = ___
     • Mass: ___
     • Stress: ___

   Improvement: (f_baseline - f_opt) / f_baseline × 100%

   Expected improvement: ___%
   Actual improvement: ___%
   ```

2. **Visualize optimized geometry**
   - [ ] Create CAD model of optimized design
   - [ ] Check for manufacturing issues
   - [ ] Verify assembly compatibility
   - [ ] Check accessibility for machining/assembly

3. **Re-run analysis with optimal values**
   - [ ] Run FEA with optimized dimensions
   - [ ] Verify stress < allowable
   - [ ] Verify deflection < limit
   - [ ] Check all constraints satisfied in full model

**Robustness and Sensitivity:**

```
Sensitivity Analysis:

For each design variable x_i:

  ∂f/∂x_i ≈ [f(x + Δx_i) - f(x)] / Δx_i

High sensitivity → Design sensitive to this variable
                 → Tight tolerances required

Low sensitivity → Variable not critical
                → Loose tolerances acceptable
```

**Robustness Check:**

```
Tolerance Simulation:

1. Define manufacturing tolerances:
   x_i = x_i,opt ± δ_i

2. Monte Carlo or worst-case:
   • Sample N combinations within tolerances
   • Evaluate objective and constraints
   • Check all constraints still satisfied

3. Robust optimization (optional):
   Minimize: f(x) + β × σ_f
   Where σ_f = standard deviation of f due to tolerances
```

**Actionable Robustness Check:**

- [ ] Identify critical variables (high sensitivity)
- [ ] Apply realistic manufacturing tolerances
- [ ] Evaluate objective with worst-case tolerances
- [ ] Ensure constraints still satisfied with tolerances
- [ ] If not robust, consider:
  - Tightening tolerances (cost increase)
  - Adding margin to constraints
  - Re-optimize with robustness objective

**Local vs. Global Optimum:**

```
Is it a global optimum?

    For simple, convex problems:
      → Local = Global

    For non-convex problems:
      [ ] Try multiple starting points
      [ ] Use global optimization algorithm
      [ ] Compare results
      [ ] Best result found is likely global
```

**Actionable Global Check:**

1. **Multi-start verification**
   - [ ] Run optimization from 5-10 random starting points
   - [ ] Compare final objectives
   - [ ] If all converge to same result → likely global
   - [ ] If different results → take best, consider global method

2. **Physical bounds check**
   - [ ] Is solution at a variable bound?
   - [ ] If yes, extend bound and re-optimize
   - [ ] Check if solution moves or stays at bound

**Documentation:**

```
Optimization Results Report:

1. Problem Statement:
   • Objective: Minimize ___
   • Variables: x1, x2, ... (n variables)
   • Constraints: g1, g2, ... (m constraints)

2. Algorithm:
   • Method: ___
   • Settings: ___
   • Convergence tolerance: ___

3. Results:
   • Initial design: f(x0) = ___
   • Optimal design: f(x*) = ___
   • Improvement: ___%
   • Iterations: ___
   • Function evaluations: ___

4. Optimal Design Variables:
   | Variable | Initial | Optimal | Change |
   |----------|---------|---------|--------|
   | x1       |         |         |        |
   | ...      |         |         |        |

5. Constraint Status:
   | Constraint | Value | Status | Margin |
   |------------|-------|--------|--------|
   | g1         |       | Active/Slack | |
   | ...        |       |        |        |

6. Verification:
   • Convergence: Met / Not met
   • Constraints: All satisfied
   • KKT conditions: Satisfied
   • Physical validity: Confirmed
   • Robustness: Acceptable / Needs attention

7. Recommendations:
   • Implement optimal design
   • Manufacturing considerations: ___
   • Tolerances: ___
   • Further work: ___
```

---

## Summary

These five structured workflows provide systematic approaches to common engineering analysis and design tasks:

1. **Fluid Dynamics Analysis**: Problem definition → Equation selection → Boundary conditions → Meshing → Solving → Verification
2. **Pump Design**: Requirements → Type selection → Sizing → Detailed design → Performance prediction → Testing
3. **Thermodynamic Analysis**: System definition → State points → Cycle calculations → Efficiency evaluation
4. **Structural Analysis (FEA)**: Load definition → Material selection → Meshing → Solver setup → Results interpretation
5. **Design Optimization**: Objective definition → Constraints → Variables → Algorithm selection → Results verification

Each workflow includes:
- Decision trees for systematic choices
- Step-by-step checklists
- Actionable procedures
- Common pitfalls and solutions
- Verification methods

Use these workflows as structured guides to ensure thorough, systematic, and verifiable engineering analysis and design.
