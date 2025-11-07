# Fluid Dynamics Workflow Examples

Comprehensive step-by-step examples demonstrating the complete CFD workflow for various practical problems.

## Table of Contents

1. [Example 1: Flow Over a Cylinder](#example-1-flow-over-a-cylinder)
2. [Example 2: Pipe Flow with Heat Transfer](#example-2-pipe-flow-with-heat-transfer)
3. [Example 3: Airfoil Analysis](#example-3-airfoil-analysis)
4. [Example 4: Turbulent Mixing](#example-4-turbulent-mixing)
5. [Example 5: Natural Convection](#example-5-natural-convection)

---

## Example 1: Flow Over a Cylinder

### Problem Statement

**Objective**: Analyze the flow over a circular cylinder and calculate drag coefficient.

**Given**:
- Cylinder diameter: D = 0.1 m
- Freestream velocity: V∞ = 10 m/s
- Fluid: Air at 20°C
- Properties: ρ = 1.225 kg/m³, μ = 1.8×10⁻⁵ Pa·s

### Step 1: Problem Definition

#### Calculate Reynolds Number
```
Re = ρV∞D/μ
Re = (1.225 kg/m³)(10 m/s)(0.1 m)/(1.8×10⁻⁵ Pa·s)
Re = 68,056

Conclusion: Turbulent flow (Re > 2000 for cylinder)
```

#### Calculate Mach Number
```
c = √(γRT) = √(1.4 × 287 × 293) = 343 m/s
Ma = V∞/c = 10/343 = 0.029

Conclusion: Incompressible (Ma < 0.3)
```

#### Expected Flow Features
- Boundary layer separation
- Vortex shedding (Kármán vortex street)
- Unsteady wake

```
Strouhal number estimate:
St ≈ 0.2 (for Re > 1000)
f = StV∞/D = 0.2 × 10/0.1 = 20 Hz
Expected shedding frequency: 20 Hz
```

### Step 2: Governing Equations

**Selected equations**:
- Continuity: ∇·V = 0 (incompressible)
- Momentum: ρ(∂V/∂t + V·∇V) = -∇p + μ∇²V
- Turbulence: k-ω SST (good for external flows with separation)

**Why k-ω SST?**
- External flow with adverse pressure gradient
- Separation expected
- Superior near-wall treatment
- Well-validated for cylinder flows

### Step 3: Computational Domain

**Domain size**:
```
Inlet:    15D upstream   (1.5 m)
Outlet:   30D downstream (3.0 m)
Top:      15D above      (1.5 m)
Bottom:   15D below      (1.5 m)
Span:     4D (3D case)   (0.4 m)

Total domain: 45D × 30D × 4D
```

**Blockage ratio**:
```
BR = D/H = 0.1/(2×1.5) = 3.3%
Acceptable (< 5% for minimal confinement effects)
```

### Step 4: Boundary Conditions

**Inlet** (Velocity inlet):
```
Velocity: u = 10 m/s, v = 0, w = 0
Turbulence intensity: I = 1%
Turbulent length scale: l = 0.07L = 0.07 × 0.1 = 0.007 m

k = 3/2(V∞I)² = 1.5 × (10 × 0.01)² = 0.015 m²/s²
ω = k^0.5/(C_μ^0.25 × l) = √0.015/(0.09^0.25 × 0.007) = 55.8 s⁻¹
```

**Outlet** (Pressure outlet):
```
Gauge pressure: 0 Pa
Backflow conditions: match inlet turbulence
```

**Cylinder surface** (Wall):
```
No-slip: u = v = w = 0
Smooth wall
```

**Top/Bottom** (Symmetry or far-field):
```
Zero normal velocity
Zero gradient for other quantities
```

**Front/Back** (3D case):
```
Option 1: Periodic (if span sufficient)
Option 2: Symmetry (if 2D assumption valid)
```

### Step 5: Mesh Generation

#### Target y⁺ Calculation

**Using wall functions** (y⁺ ≈ 50):
```
Estimate skin friction:
C_f ≈ 0.058Re_D^(-0.2) = 0.058 × 68056^(-0.2) = 0.0088

u_τ = V∞√(C_f/2) = 10√(0.0088/2) = 0.664 m/s

First cell height:
y = y⁺μ/(ρu_τ) = 50 × 1.8×10⁻⁵/(1.225 × 0.664)
y = 1.10×10⁻³ m = 1.1 mm
```

#### Mesh Strategy

**Near-wall region** (structured):
```
First cell height: 1.1 mm (y⁺ ≈ 50)
Number of prism layers: 15
Growth rate: 1.2
Total inflation thickness: ~0.2D
```

**Wake region** (refined):
```
Cell size: D/40 = 2.5 mm
Extend refinement 10D downstream
Width: 4D laterally
```

**Far field** (coarse):
```
Cell size: D/5 = 20 mm
Gradual transition from refined region
```

**Total cell count**:
- 2D: ~150,000 cells
- 3D: ~4,000,000 cells

#### Mesh Quality Targets
```
Aspect ratio: < 100 (< 500 in wake)
Skewness: < 0.8
Orthogonality: > 0.3
```

### Step 6: Solver Settings

**Solver type**: Pressure-based, coupled

**Discretization**:
```
Pressure: Second-order
Momentum: Second-order upwind
Turbulence: Second-order upwind
```

**Transient settings**:
```
Time step: Δt = 0.001 s
  (Captures ~20 time steps per shedding cycle: T = 1/f = 0.05 s)

CFL check:
CFL = VΔt/Δx = 10 × 0.001/0.0025 = 4
Acceptable for implicit solver

Total time: 5 seconds (100 shedding cycles)
Inner iterations: 20 per time step
```

**Convergence criteria**:
```
Residuals: < 1×10⁻⁴ (all variables)
Monitor: C_D, C_L (should oscillate regularly)
```

### Step 7: Solution Strategy

**Initialization**:
1. Start with steady RANS solution (100-200 iterations)
2. Switch to transient
3. Run 1-2 seconds to flush transients
4. Begin data collection

**Under-relaxation** (steady phase):
```
Pressure: 0.3
Momentum: 0.7
Turbulence: 0.8
```

**Monitoring**:
- Plot C_D and C_L vs. time
- Check for periodic shedding
- Monitor mass flow balance
- Visualize vorticity field

### Step 8: Post-Processing

#### Force Calculations

**Drag coefficient**:
```
F_D = F_D,pressure + F_D,viscous
C_D = F_D/(1/2 ρV²∞A)

Where A = D × span (3D) or D × 1 (2D per unit depth)

Expected: C_D ≈ 1.0-1.2 for Re ≈ 68,000
```

**Lift coefficient**:
```
C_L = F_L/(1/2 ρV²∞A)

Should oscillate: C_L,rms ≈ 0.3-0.5
```

#### Flow Visualization

**Velocity contours**:
- Identify separation points (θ ≈ 80°-90° from stagnation)
- Visualize wake structure

**Vorticity contours**:
```
ω_z = ∂v/∂x - ∂u/∂y

Shows vortex street clearly
Measure vortex spacing
```

**Pressure coefficient**:
```
C_p = (p - p∞)/(1/2 ρV²∞)

Expected:
- Stagnation point: C_p = 1.0
- Minimum: C_p ≈ -1.5 to -2.0
- Base: C_p ≈ -1.2 to -0.8
```

**Streamlines**:
- Show recirculation zones
- Visualize separation/reattachment

#### Frequency Analysis

**Extract lift force time series**:
```
Perform FFT (Fast Fourier Transform)
Identify dominant frequency: f_shedding

Calculate Strouhal number:
St = fD/V∞

Compare with experimental: St ≈ 0.2
```

### Step 9: Validation

**Compare with correlations**:
```
Drag coefficient (Schlichting):
C_D ≈ 1.2 for Re = 68,000

Strouhal number:
St ≈ 0.21 for Re > 1000

Roshko correlation:
St = 0.212(1 - 21.2/Re)
St = 0.212(1 - 21.2/68056) = 0.206
```

**Grid independence**:
```
Test three meshes:
- Coarse: 75,000 cells
- Medium: 150,000 cells
- Fine: 300,000 cells

Calculate GCI for C_D:
GCI = 1.25|ε|/(r^p - 1)
Where ε = (C_D,fine - C_D,medium)/C_D,fine
```

### Step 10: Results Summary

**Expected results**:
```
Drag coefficient: C_D = 1.1 ± 0.05
Lift coefficient (rms): C_L,rms = 0.4 ± 0.1
Strouhal number: St = 0.20 ± 0.02
Separation angle: θ_sep = 85° ± 5°
```

---

## Example 2: Pipe Flow with Heat Transfer

### Problem Statement

**Objective**: Analyze fully developed turbulent flow in a heated pipe and calculate heat transfer coefficient.

**Given**:
- Pipe diameter: D = 0.05 m
- Length: L = 2.5 m (50D)
- Mass flow rate: ṁ = 0.5 kg/s
- Inlet temperature: T_in = 300 K
- Wall heat flux: q_w = 10,000 W/m²
- Fluid: Water at 20°C
- Properties: ρ = 998 kg/m³, μ = 1.002×10⁻³ Pa·s, c_p = 4182 J/(kg·K), k = 0.6 W/(m·K)

### Step 1: Problem Definition

#### Calculate Flow Parameters
```
Cross-sectional area:
A = πD²/4 = π(0.05)²/4 = 1.963×10⁻³ m²

Bulk velocity:
V = ṁ/(ρA) = 0.5/(998 × 1.963×10⁻³) = 0.255 m/s

Reynolds number:
Re_D = ρVD/μ = (998)(0.255)(0.05)/(1.002×10⁻³)
Re_D = 12,725

Conclusion: Turbulent flow (Re > 2300)
```

#### Calculate Heat Transfer Parameters
```
Prandtl number:
Pr = μc_p/k = (1.002×10⁻³)(4182)/(0.6) = 6.98

Peclet number:
Pe = Re × Pr = 12,725 × 6.98 = 88,822

High Pe → Convection dominated

Thermal entrance length:
L_t/D ≈ 10 (turbulent)
L_t ≈ 0.5 m
```

#### Expected Results
```
Nusselt number (Dittus-Boelter):
Nu = 0.023Re^0.8 Pr^0.4
Nu = 0.023 × (12,725)^0.8 × (6.98)^0.4
Nu = 91.6

Heat transfer coefficient:
h = Nu × k/D = 91.6 × 0.6/0.05 = 1,099 W/(m²·K)

Outlet temperature estimate:
Q_total = q_w × πDL = 10,000 × π × 0.05 × 2.5 = 3,927 W
ΔT = Q/(ṁc_p) = 3,927/(0.5 × 4182) = 1.88 K
T_out ≈ 301.88 K
```

### Step 2: Governing Equations

**Selected equations**:
- Continuity: ∇·V = 0
- Momentum: ρ(∂V/∂t + V·∇V) = -∇p + μ∇²V
- Energy: ρc_p(∂T/∂t + V·∇T) = ∇·(k∇T)
- Turbulence: k-ε Standard (suitable for pipe flow)

**Why k-ε Standard?**
- Fully developed pipe flow
- No separation or adverse gradients
- Well-validated for this geometry
- Computationally efficient

### Step 3: Computational Domain

**Domain geometry**:
```
Option 1: Full 3D pipe (computational expensive)
Option 2: Axisymmetric (preferred)
Option 3: Periodic 2D slice (for developed flow)

Selected: Axisymmetric
- Reduces 3D to 2D
- Captures all physics
- Much faster solution
```

**Domain size**:
```
Length: L = 2.5 m (50D)
Radius: R = 0.025 m

Sufficient for:
- Hydrodynamic development
- Thermal development
- Fully developed region
```

### Step 4: Boundary Conditions

**Inlet**:
```
Option 1: Velocity inlet
- Mass flow rate: ṁ = 0.5 kg/s
- Temperature: T = 300 K
- Turbulent intensity: I = 5%
- Hydraulic diameter: D_h = D = 0.05 m

k = 3/2(VI)² = 1.5 × (0.255 × 0.05)² = 2.44×10⁻⁴ m²/s²
ε = C_μ^0.75 k^1.5 / l
  where l = 0.07D_h = 0.0035 m
ε = 0.09^0.75 × (2.44×10⁻⁴)^1.5 / 0.0035 = 1.02×10⁻³ m²/s³

Option 2: Fully developed profile (if available)
```

**Outlet**:
```
Pressure outlet: p = 0 Pa (gauge)
Temperature: ∂T/∂x = 0 (fully developed)
```

**Wall**:
```
No-slip: V = 0
Heat flux: q_w = 10,000 W/m²
Smooth wall
```

**Axis** (centerline):
```
Axisymmetric boundary condition
∂φ/∂r = 0 for all variables
```

### Step 5: Mesh Generation

#### Near-Wall Resolution

**Target y⁺ = 1** (resolve thermal boundary layer):
```
For heated wall, need low-Re treatment

Estimate u_τ:
C_f ≈ 0.079Re_D^(-0.25) = 0.079 × (12,725)^(-0.25) = 0.0075
u_τ = V√(C_f/2) = 0.255√(0.0075/2) = 0.0156 m/s

First cell height:
y = y⁺μ/(ρu_τ) = 1.0 × (1.002×10⁻³)/(998 × 0.0156)
y = 6.44×10⁻⁵ m = 64.4 μm
```

#### Mesh Structure

**Radial direction** (structured):
```
First cell: 64.4 μm (y⁺ = 1)
Number of cells: 80
Growth rate: 1.1
Bias toward wall (70% of cells in boundary layer)
```

**Axial direction** (structured):
```
Upstream region (0-0.5 m): 100 cells
Developed region (0.5-2.5 m): 200 cells
Total axial cells: 300

Cell size: ~8-10 mm axially
```

**Total cells**: 300 × 80 = 24,000 (axisymmetric 2D)

#### Mesh Quality
```
Aspect ratio: < 50 (acceptable for pipe flow)
Skewness: < 0.2 (structured mesh)
Orthogonality: > 0.95
```

### Step 6: Solver Settings

**Solver**: Pressure-based, steady

**Discretization**:
```
Pressure: PRESTO! (for high aspect ratio)
Momentum: Second-order upwind
Energy: Second-order upwind
Turbulence: First-order upwind (for stability)
```

**Solution controls**:
```
Under-relaxation:
- Pressure: 0.3
- Momentum: 0.7
- Energy: 0.95
- Turbulence: 0.8
```

**Convergence criteria**:
```
Residuals:
- Continuity: < 1×10⁻⁴
- Velocity: < 1×10⁻⁴
- Energy: < 1×10⁻⁶
- Turbulence: < 1×10⁻⁴

Monitors:
- Outlet temperature
- Pressure drop
- Heat transfer rate
```

### Step 7: Solution Strategy

**Initialization**:
1. Initialize with inlet conditions
2. Run with first-order schemes (50 iterations)
3. Switch to second-order
4. Continue until convergence

**Monitoring**:
```
Check every 50 iterations:
- Residual trends
- Mass flow balance (inlet = outlet)
- Energy balance (Q_in = Q_wall)
- Outlet temperature
```

**Energy balance check**:
```
Q_wall = q_w × A_wall = 10,000 × (πDL)
Q_fluid = ṁc_p(T_out - T_in)

Should match within 1%
```

### Step 8: Post-Processing

#### Temperature Profile

**Radial temperature at various x locations**:
```
Plot T(r) at:
- x = 0.1 m (entrance)
- x = 0.5 m (developing)
- x = 1.5 m (developed)
- x = 2.5 m (outlet)

Show thermal boundary layer growth
```

**Non-dimensional temperature**:
```
θ = (T - T_centerline)/(T_wall - T_centerline)

Should show universal profile in developed region
```

#### Velocity Profile

**Check for fully developed flow**:
```
At x > L_h = 10D = 0.5 m

Power law:
u(r)/U_max = (1 - r/R)^(1/n)
n ≈ 7 for Re ≈ 10⁴

Log law (near wall):
u⁺ = (1/κ)ln(y⁺) + B
```

#### Heat Transfer Coefficient

**Local h(x)**:
```
h(x) = q_w/(T_w(x) - T_bulk(x))

T_bulk(x) = (1/A)∫∫ ρVc_pT dA / (ṁc_p)

Plot h(x):
- Should be high near entrance
- Decrease to constant in developed region
```

**Average Nusselt number**:
```
h_avg = (1/L)∫₀ᴸ h(x) dx

Nu_avg = h_avg D/k

Compare with Dittus-Boelter:
Nu = 0.023Re^0.8 Pr^0.4 = 91.6
```

#### Pressure Drop

**Calculate friction factor**:
```
Δp = f(L/D)(ρV²/2)

f = 2Δp D/(ρV²L)

Compare with Colebrook (smooth pipe):
1/√f = -2.0 log₁₀(2.51/(Re√f))

For Re = 12,725: f ≈ 0.0268
```

### Step 9: Validation

**Heat transfer validation**:
```
Dittus-Boelter equation:
Nu = 0.023Re^0.8 Pr^n
n = 0.4 (heating), 0.3 (cooling)

Nu_predicted = 91.6
Nu_CFD should be within ±10%
```

**Friction factor validation**:
```
Blasius equation (smooth pipe):
f = 0.316Re^(-0.25) = 0.316 × (12,725)^(-0.25)
f = 0.0299

Colebrook equation (more accurate):
f ≈ 0.0268

CFD should match within 5%
```

**Grid independence**:
```
Test meshes:
- Coarse: 150 × 40
- Medium: 300 × 80
- Fine: 600 × 160

Monitor:
- Average Nusselt number
- Pressure drop
- Outlet temperature
```

### Step 10: Results Summary

**Expected outputs**:
```
Outlet temperature: T_out = 301.88 K
Average heat transfer coefficient: h = 1,100 W/(m²·K)
Average Nusselt number: Nu = 91.6
Pressure drop: Δp ≈ 200 Pa
Friction factor: f = 0.027
```

**Result presentation**:
- Contour plots: T(x,r), u(x,r)
- Line plots: T(r) at various x, u(r) profile
- Graphs: h(x), Nu(x), f(x)
- Tables: Comparison with correlations

---

## Example 3: Airfoil Analysis

### Problem Statement

**Objective**: Analyze flow over NACA 0012 airfoil at angle of attack and calculate lift and drag.

**Given**:
- Airfoil: NACA 0012
- Chord: c = 1.0 m
- Angle of attack: α = 4°
- Velocity: V∞ = 50 m/s
- Fluid: Air at sea level
- Properties: ρ = 1.225 kg/m³, μ = 1.789×10⁻⁵ Pa·s

### Step 1: Problem Definition

#### Flow Parameters
```
Reynolds number:
Re_c = ρV∞c/μ = (1.225)(50)(1.0)/(1.789×10⁻⁵)
Re_c = 3.42 × 10⁶

Conclusion: Turbulent boundary layer

Transition location: x/c ≈ 0.05-0.1 (both surfaces)

Mach number:
Ma = V∞/c = 50/343 = 0.146

Conclusion: Incompressible (Ma < 0.3)
```

#### Expected Aerodynamic Characteristics
```
From thin airfoil theory:
C_L,theory = 2πα = 2π(4π/180) = 0.436

From experimental data (NACA 0012):
C_L,exp ≈ 0.44 at α = 4°
C_D,exp ≈ 0.0080
C_M,c/4 ≈ 0 (symmetric airfoil)
```

### Step 2: Governing Equations

**Selected equations**:
- Incompressible RANS
- k-ω SST turbulence model
- Transition SST (optional, for better accuracy)

**Why k-ω SST?**
- External aerodynamics
- Adverse pressure gradient on upper surface
- Excellent pressure gradient handling
- Industry standard for airfoils

### Step 3: Computational Domain

**C-mesh topology** (most common for airfoils):
```
Far-field distance: 20c from airfoil
  (20 m in all directions)

Domain blocks:
- Inner block: structured O-grid around airfoil
- Outer block: structured C-grid to far field
```

**Far-field boundary check**:
```
Verify:
- Velocity perturbation < 1% at boundary
- Pressure coefficient → 0 at boundary
```

### Step 4: Boundary Conditions

**Far-field** (velocity inlet/pressure outlet combined):
```
Velocity components:
u = V∞cosα = 50cos(4°) = 49.88 m/s
v = V∞sinα = 50sin(4°) = 3.49 m/s

Or use far-field boundary condition:
- Mach number: 0.146
- Static pressure: 101,325 Pa
- Static temperature: 288 K
- Flow direction: α = 4°

Turbulent intensity: I = 0.1% (low for wind tunnel)
Turbulent viscosity ratio: μₜ/μ = 1-10
```

**Airfoil surface**:
```
No-slip wall: V = 0
Smooth surface (or specify roughness)
Adiabatic (∂T/∂n = 0, if solving energy)
```

**Upper/lower** (if using rectangular domain):
```
Symmetry or far-field
```

### Step 5: Mesh Generation

#### Surface Mesh

**Leading edge**:
```
Very fine clustering at stagnation point
Cell size: ~0.0001c = 0.1 mm
Spacing ratio: 1.05
```

**Surface distribution**:
```
Number of points on airfoil: 400-600
- Upper surface: 250
- Lower surface: 250
- Leading edge: dense clustering
- Trailing edge: fine resolution

Non-dimensional spacing:
Δs/c < 0.001 at leading/trailing edges
Δs/c ≈ 0.01 on mid-chord regions
```

#### Boundary Layer Mesh

**First cell height** (y⁺ < 1):
```
Estimate C_f:
C_f ≈ 0.058Re_c^(-0.2) = 0.058 × (3.42×10⁶)^(-0.2)
C_f = 0.00264

u_τ = V∞√(C_f/2) = 50√(0.00264/2) = 1.816 m/s

For y⁺ = 1:
y = y⁺μ/(ρu_τ) = 1 × (1.789×10⁻⁵)/(1.225 × 1.816)
y = 8.04×10⁻⁶ m = 8.04 μm

Use: y₁ = 8 μm
```

**Prism layer specifications**:
```
Number of layers: 30-40
First layer: 8 μm (y⁺ ≈ 1)
Growth rate: 1.15
Total thickness: ~0.02c = 20 mm

Ensures:
- y⁺ < 1 everywhere
- 15-20 cells within boundary layer
- Smooth transition to isotropic region
```

#### Off-body Mesh

**Wake refinement**:
```
Extend: 5c downstream
Width: ±0.5c from airfoil
Cell size: 0.01c = 10 mm
```

**Far-field**:
```
Cell size: 1.0c = 1000 mm (coarse)
Gradual coarsening from wake to far-field
```

**Total cell count**:
- 2D: ~100,000-200,000 cells
- 3D (with spanwise): 5-10 million cells

### Step 6: Solver Settings

**Solver**: Pressure-based, coupled

**Discretization**:
```
Pressure: Second-order
Momentum: Second-order upwind
Turbulence: Second-order upwind
```

**Spatial discretization**:
```
Gradient: Least squares cell-based
  or Green-Gauss node-based (higher accuracy)
```

**Solution controls**:
```
Courant number: 5-10 (pseudo-transient)

Under-relaxation (if segregated):
- Pressure: 0.3
- Momentum: 0.7
- Turbulence: 0.8
```

### Step 7: Solution Strategy

**Initialization**:
```
Hybrid initialization (from far-field conditions)
Or initialize from potential flow solution
```

**Solution sequence**:
1. First-order upwind schemes (100 iterations)
2. Switch to second-order (500-1000 iterations)
3. Enable pseudo-transient (if convergence issues)

**Convergence monitoring**:
```
Residuals: < 1×10⁻⁵ (tight for accurate forces)

Force coefficients:
- C_L should stabilize within ±0.001
- C_D should stabilize within ±0.0001
- C_M should stabilize within ±0.001

Monitors:
- Lift coefficient: C_L
- Drag coefficient: C_D
- Moment coefficient: C_M,c/4
```

### Step 8: Post-Processing

#### Force Coefficients

**Lift coefficient**:
```
L = ∫ (p_lower - p_upper)cosθ ds + ∫ τ_w sinθ ds
C_L = L/(1/2 ρV²∞S)

Where S = c × span (per unit span for 2D)

Pressure component: ~95% of lift
Viscous component: ~5%

Expected: C_L ≈ 0.44
```

**Drag coefficient**:
```
D = ∫ (p_lower - p_upper)sinθ ds + ∫ τ_w cosθ ds
C_D = D/(1/2 ρV²∞S)

Components:
C_D = C_D,pressure + C_D,friction

Expected:
C_D,friction ≈ 0.0060 (dominant at low α)
C_D,pressure ≈ 0.0020 (form drag)
C_D,total ≈ 0.0080
```

**Moment coefficient** (about quarter-chord):
```
M_c/4 = ∫ (p-p∞)(x-c/4) ds + ∫ τ_w moments
C_M,c/4 = M/(1/2 ρV²∞Sc)

Expected: C_M,c/4 ≈ 0 (symmetric airfoil)
```

**Aerodynamic center**: x_ac ≈ 0.25c for thin airfoils

#### Pressure Distribution

**Pressure coefficient**:
```
C_p = (p - p∞)/(1/2 ρV²∞)

Plot C_p vs. x/c:
- Upper surface: negative (suction)
- Lower surface: mostly positive (pressure)
- Leading edge suction peak: C_p ≈ -2 to -3
```

**Pressure integration**:
```
Normal force coefficient:
C_N = ∫₀¹ (C_p,lower - C_p,upper) d(x/c)

Axial force coefficient:
C_A ≈ 0.01-0.02 (small at low α)

C_L = C_Ncosα - C_Asinα ≈ C_N (for small α)
C_D = C_Nsinα + C_Acosα ≈ C_A (for small α)
```

#### Skin Friction

**Friction coefficient distribution**:
```
C_f = τ_w/(1/2 ρV²∞)

Plot C_f vs. x/c:
- Identify transition location (C_f minimum)
- Check for separation (C_f = 0)
- Turbulent: C_f ~ x^(-0.2)
```

**Boundary layer properties**:
```
Displacement thickness: δ*(x)
Momentum thickness: θ(x)
Shape factor: H(x) = δ*/θ

H ≈ 2.6 (laminar)
H ≈ 1.3-1.4 (turbulent)
H > 2.5 (approaching separation)
```

#### Flow Visualization

**Velocity contours**:
- Show boundary layer growth
- Identify wake structure

**Streamlines**:
- Flow acceleration over upper surface
- Stagnation streamline

**Vorticity**:
```
ω_z = ∂v/∂x - ∂u/∂y

Concentrated in:
- Boundary layers
- Wake
- Separation bubbles (if any)
```

**y⁺ distribution**:
```
Check that y⁺ ≤ 1 everywhere on surface
Maximum y⁺ usually at leading edge stagnation
```

### Step 9: Validation

**Experimental comparison**:
```
Data sources:
- Abbott & von Doenhoff airfoil data
- NASA wind tunnel tests
- NACA reports

Compare:
1. C_L vs. α curve
2. C_D vs. α curve
3. C_p distribution at α = 4°
```

**Theoretical comparison**:
```
Thin airfoil theory:
C_L = 2π(α + α_L=0)
For NACA 0012: α_L=0 = 0
C_L,theory = 2π(4π/180) = 0.436

CFD should give: C_L ≈ 0.44 (within 1%)
```

**Grid independence**:
```
Three mesh levels:
- Coarse: 50K cells, y₁ = 20 μm
- Medium: 100K cells, y₁ = 10 μm
- Fine: 200K cells, y₁ = 5 μm

Calculate GCI for C_L and C_D:
Should show < 1% change from medium to fine
```

**Turbulence model sensitivity**:
```
Compare:
- k-ω SST
- Spalart-Allmaras
- Transition SST

Usually < 5% difference in C_L
C_D more sensitive (10-20% variation)
```

### Step 10: Results Summary

**Aerodynamic coefficients**:
```
Lift coefficient: C_L = 0.440 ± 0.005
Drag coefficient: C_D = 0.0080 ± 0.0002
  - Friction drag: C_D,f = 0.0060
  - Pressure drag: C_D,p = 0.0020
Moment coefficient: C_M,c/4 = 0.000 ± 0.001
Lift-to-drag ratio: L/D = 55.0
```

**Polar data** (for multiple α):
```
α (deg)  |  C_L  |  C_D   |  L/D
--------------------------------
   0     | 0.000 | 0.0070 |  0
   2     | 0.220 | 0.0072 | 30.6
   4     | 0.440 | 0.0080 | 55.0
   6     | 0.658 | 0.0095 | 69.3
   8     | 0.872 | 0.0118 | 73.9
  10     | 1.080 | 0.0150 | 72.0
```

---

## Example 4: Turbulent Mixing

### Problem Statement

**Objective**: Analyze turbulent mixing of two co-flowing streams in a Y-junction mixer.

**Given**:
- Inlet 1: Water, V₁ = 2 m/s, T₁ = 300 K, ṁ₁ = 1 kg/s
- Inlet 2: Water, V₂ = 1 m/s, T₂ = 350 K, ṁ₂ = 0.5 kg/s
- Pipe diameter: D = 0.05 m each inlet
- Mixing section diameter: D_mix = 0.07 m
- Length: L = 1.0 m (20D_mix)

### Step 1: Problem Definition

#### Flow Characterization
```
Inlet 1:
Re₁ = ρV₁D/μ = 998 × 2 × 0.05 / 1.002×10⁻³ = 99,601

Inlet 2:
Re₂ = ρV₂D/μ = 998 × 1 × 0.05 / 1.002×10⁻³ = 49,800

Both turbulent

Mixed flow:
V_mix = (ṁ₁ + ṁ₂)/(ρA_mix)
V_mix = 1.5/(998 × π × 0.07²/4) = 0.39 m/s

Re_mix = 998 × 0.39 × 0.07 / 1.002×10⁻³ = 27,000
```

#### Expected Outlet Temperature
```
Energy balance (neglecting losses):
T_out = (ṁ₁c_pT₁ + ṁ₂c_pT₂)/(ṁ₁ + ṁ₂)
T_out = (1.0 × 300 + 0.5 × 350)/(1.0 + 0.5)
T_out = 316.7 K
```

#### Mixing Length
```
Turbulent mixing length: L_mix ≈ 10-20 diameters
L_mix ≈ 10 × 0.07 = 0.7 m

Domain of 1.0 m should show good mixing
```

### Step 2: Governing Equations

**Selected model**:
- Continuity, momentum, energy
- k-ε Realizable (good for mixing)
- Or k-ω SST

### Step 3-10: Solution Steps

*(Abbreviated for brevity - follows similar pattern to previous examples)*

**Key aspects**:
1. Two velocity-inlet boundaries
2. Fine mesh in mixing region
3. Monitor temperature and species mixing
4. Calculate mixing efficiency
5. Visualize temperature contours and streamlines

**Expected results**:
- Temperature uniformity at outlet
- Pressure drop calculation
- Mixing length determination
- Turbulent intensity distribution

---

## Example 5: Natural Convection in Enclosure

### Problem Statement

**Objective**: Analyze natural convection in a square cavity with differentially heated walls.

**Given**:
- Cavity size: H × W = 1 m × 1 m
- Hot wall (left): T_h = 310 K
- Cold wall (right): T_c = 290 K
- Top/bottom: Adiabatic
- Fluid: Air
- Pressure: 101,325 Pa

### Step 1: Problem Definition

#### Calculate Rayleigh Number
```
ΔT = T_h - T_c = 20 K
T_avg = 300 K

β = 1/T_avg = 1/300 = 3.33×10⁻³ K⁻¹

ν = 1.5×10⁻⁵ m²/s
α = 2.2×10⁻⁵ m²/s

Grashof number:
Gr = gβΔTH³/ν²
Gr = 9.81 × 3.33×10⁻³ × 20 × 1³ / (1.5×10⁻⁵)²
Gr = 2.9×10⁹

Rayleigh number:
Ra = Gr × Pr = 2.9×10⁹ × (0.71) = 2.06×10⁹

Conclusion: Turbulent natural convection (Ra > 10⁹)
```

#### Expected Flow Pattern
```
For Ra > 10⁶:
- Thin boundary layers on hot/cold walls
- Rising plume along hot wall
- Falling plume along cold wall
- Stratified core region
- Multiple recirculation cells possible
```

### Step 2-10: Solution Steps

*(Abbreviated)*

**Key aspects**:
1. Boussinesq approximation or full buoyancy
2. Very fine mesh near walls (boundary layers)
3. Transient solution (natural convection often unsteady)
4. Monitor heat transfer (Nusselt number)
5. Visualize temperature and velocity fields

**Expected results**:
```
Average Nusselt number:
Nu = hH/k ≈ 0.1Ra^(1/3) = 0.1 × (2.06×10⁹)^(1/3) ≈ 127

Heat transfer:
Q = Nu × k × A × ΔT/H
Q ≈ 127 × 0.026 × 1 × 20/1 ≈ 66 W
```

---

## Summary of Best Practices

### Across All Examples

1. **Always calculate dimensionless numbers first**
   - Determines flow regime
   - Guides model selection
   - Provides validation benchmarks

2. **Start with analytical/empirical estimates**
   - Sanity check for results
   - Initial conditions
   - Validation targets

3. **Mesh systematically**
   - Calculate required y⁺
   - Refine critical regions
   - Perform grid independence study

4. **Monitor convergence rigorously**
   - Residuals
   - Integrated quantities
   - Mass/energy balance

5. **Validate thoroughly**
   - Compare with experiments
   - Compare with theory
   - Check physical realizability

6. **Document assumptions**
   - Geometry simplifications
   - Boundary condition choices
   - Model selections
   - Limitations

### Common Pitfalls to Avoid

1. Insufficient domain size
2. Poor near-wall resolution
3. Inadequate convergence
4. Wrong turbulence model
5. Unrealistic boundary conditions
6. Lack of validation
7. No uncertainty quantification

### Reporting Checklist

- [ ] Problem statement and objectives
- [ ] Dimensionless numbers and regime
- [ ] Domain and boundary conditions
- [ ] Mesh description and quality metrics
- [ ] Solver settings and convergence
- [ ] Results with uncertainty estimates
- [ ] Validation against experiments/theory
- [ ] Conclusions and recommendations

---

These examples demonstrate the complete workflow from problem definition through validation. Adapt the approach based on your specific problem requirements while maintaining the systematic methodology.
