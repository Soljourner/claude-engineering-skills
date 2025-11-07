# Pump Efficiency Optimization - Technical Reference

## Fundamental Equations

### 1. Efficiency Definitions

#### Overall Pump Efficiency

```
η_overall = P_hydraulic / P_shaft

P_hydraulic = ρ × g × Q × H

where:
  ρ = fluid density (kg/m³)
  g = gravitational acceleration (9.81 m/s²)
  Q = flow rate (m³/s)
  H = head (m)
  P_shaft = shaft power (W)
```

#### Component Efficiencies

**Hydraulic Efficiency:**
```
η_h = (g × H) / (u₂ × c_u2)

where:
  u₂ = impeller tip velocity = ω × r₂ (m/s)
  c_u2 = tangential component of absolute velocity at exit (m/s)
  ω = angular velocity (rad/s)
```

**Volumetric Efficiency:**
```
η_v = Q_delivered / (Q_delivered + Q_leakage)

Q_leakage = (π × D × c³ × ΔP) / (12 × μ × L)

where:
  D = seal diameter (m)
  c = radial clearance (m)
  ΔP = pressure differential (Pa)
  μ = dynamic viscosity (Pa·s)
  L = seal length (m)
```

**Mechanical Efficiency:**
```
η_m = P_hydraulic / (P_hydraulic + P_friction)

P_friction = P_bearing + P_seal + P_disk

where:
  P_bearing = bearing friction power
  P_seal = seal friction power
  P_disk = disk friction power
```

**Overall Efficiency Relationship:**
```
η_overall = η_h × η_v × η_m
```

### 2. Euler Turbomachine Equation

The fundamental equation relating head to impeller geometry:

```
g × H = u₂ × c_u2 - u₁ × c_u1

For radial entry (c_u1 = 0):
H = (u₂ × c_u2) / g
```

**Velocity triangles:**
```
u = ω × r                    (tangential velocity)
c_m = Q / (π × D × b)        (meridional velocity)
c_u = u - c_m / tan(β)       (tangential component, absolute frame)

where:
  β = blade angle from tangential direction
```

### 3. Specific Speed

Dimensionless parameter characterizing pump type:

**US Customary Units:**
```
N_s = N × √Q / H^(3/4)

where:
  N = rotational speed (RPM)
  Q = flow rate (GPM)
  H = head (ft)
```

**SI Units:**
```
Ω_s = ω × √Q / (g × H)^(3/4)

where:
  ω = angular velocity (rad/s)
  Q = flow rate (m³/s)
  H = head (m)
```

**Pump type selection:**
- N_s < 2000: Radial flow (centrifugal)
- 2000 < N_s < 7000: Mixed flow
- N_s > 7000: Axial flow

### 4. Affinity Laws

For geometrically similar pumps:

**Speed variation (constant diameter):**
```
Q₂/Q₁ = N₂/N₁
H₂/H₁ = (N₂/N₁)²
P₂/P₁ = (N₂/N₁)³
```

**Diameter variation (constant speed):**
```
Q₂/Q₁ = D₂/D₁
H₂/H₁ = (D₂/D₁)²
P₂/P₁ = (D₂/D₁)³
```

**Combined variation:**
```
Q₂/Q₁ = (N₂/N₁) × (D₂/D₁)
H₂/H₁ = (N₂/N₁)² × (D₂/D₁)²
P₂/P₁ = (N₂/N₁)³ × (D₂/D₁)³
```

## Loss Mechanisms

### 1. Hydraulic Losses

**Friction losses in passages:**
```
h_f = f × (L/D_h) × (v²/2g)

where:
  f = friction factor (function of Reynolds number and roughness)
  L = passage length (m)
  D_h = hydraulic diameter = 4A/P (m)
  v = flow velocity (m/s)
```

**Darcy-Weisbach friction factor:**
```
1/√f = -2.0 × log₁₀(ε/(3.7D) + 2.51/(Re√f))

where:
  ε = surface roughness (m)
  Re = Reynolds number = ρvD/μ
```

**Shock losses (incidence):**
```
h_shock = k × (Δc)²/(2g)

where:
  Δc = difference between flow angle and blade angle
  k = shock coefficient (typically 0.5-1.0)
```

**Separation losses:**
Occur when adverse pressure gradient causes flow detachment. Minimize by:
- Gradual area changes (diffuser angle < 7-10°)
- Smooth transitions
- Avoid sharp corners

### 2. Disk Friction Losses

Power consumed by impeller surface rotation in fluid:

**Enclosed disk (small clearance):**
```
P_disk = k × ρ × ω³ × r⁵ × (Re_r)^(-0.2)

where:
  k ≈ 0.073 (empirical constant)
  Re_r = ρ × ω × r² / μ (rotational Reynolds number)
```

**Open disk (large clearance):**
```
P_disk = 0.5 × ρ × C_f × ω³ × r⁵

where:
  C_f = skin friction coefficient ≈ 0.005-0.01
```

**Reduction methods:**
- Minimize impeller outside diameter
- Use pump-out vanes on shrouds
- Optimize clearances (too small → high friction, too large → recirculation)

### 3. Leakage Losses

**Through annular clearances:**
```
Q_leak = (π × D × c³ × ΔP) / (12 × μ × L) × (1 + 1.5 × (c/L)²)

Power loss:
P_leak = Q_leak × ΔP
```

**Typical clearances:**
- Wear rings: c = 0.010 - 0.025 inches per inch of diameter
- Example: 6" diameter → clearance = 0.060 - 0.150 inches

**Minimization strategies:**
- Close tolerances (limited by thermal expansion, deflection, wear)
- Labyrinth seals (multiple restrictions)
- Hard facings (tungsten carbide, ceramic) for wear resistance
- Proper material selection (minimize galling)

### 4. Recirculation Losses

**Suction recirculation:**
Occurs at low flow (< 60-70% BEP):
```
NPSH_req increases dramatically
Noise, vibration, cavitation damage

Onset flow:
Q_sr ≈ 0.6 × Q_BEP
```

**Discharge recirculation:**
Occurs at high flow (> 120% BEP):
```
Flow separation at impeller exit
Efficiency drop
Possible head instability

Onset flow:
Q_dr ≈ 1.2 × Q_BEP
```

**Prevention:**
- Operate within recommended flow range (70-120% BEP)
- Use inlet guide vanes for low-flow operation
- Variable speed control to match system demand

## Design Optimization Equations

### 1. Impeller Design

**Blade angle calculation:**

**Inlet angle (for shock-free entry):**
```
β₁ = arctan(c_m1 / u₁)

where:
  c_m1 = meridional velocity at inlet
  u₁ = tangential velocity at inlet
```

**Outlet angle (for target head):**
```
tan(β₂) = c_m2 / (u₂ - c_u2)

Given target head H:
c_u2 = (g × H) / u₂

Therefore:
β₂ = arctan(c_m2 / (u₂ - (g × H)/u₂))
```

**Optimal outlet angle range:**
- β₂ = 15-20°: High efficiency, lower head
- β₂ = 20-25°: Balanced design (most common)
- β₂ = 25-35°: Higher head, lower efficiency

**Number of blades:**
```
Z = 6.5 × (D₂ + D₁)/(D₂ - D₁) × sin((β₁ + β₂)/2)

Typical ranges:
  Centrifugal pumps: 5-7 blades
  Mixed flow: 4-6 blades
  Axial flow: 3-5 blades
```

**Impeller width:**
```
b₂ = Q / (π × D₂ × c_m2 × η_v)

Width ratio:
b₂/D₂ ≈ 0.03-0.15 (varies with specific speed)
```

### 2. Volute Design

**Volute throat area:**
```
A_throat = Q / c_throat

where:
  c_throat ≈ 0.15-0.20 × u₂ (design velocity)
```

**Volute expansion angle:**
```
α < 7-10° to avoid flow separation
```

**Cutwater clearance:**
```
Gap ≈ 0.03-0.05 × D₂
```

### 3. Best Efficiency Point (BEP)

**Empirical correlations for efficiency:**

**Hydraulic efficiency:**
```
η_h = 1 - 0.8 / N_s^0.17    (for well-designed pumps)

Typical values:
  Radial (N_s < 2000): 85-92%
  Mixed (2000-7000): 88-94%
  Axial (> 7000): 86-92%
```

**Overall efficiency estimates:**

```
Anderson correlation:
η = 1 - 0.095 / (Q_BEP [m³/s])^0.16

Typical ranges:
  Small pumps (< 10 HP): 30-60%
  Medium pumps (10-100 HP): 60-80%
  Large pumps (> 100 HP): 80-90%
```

## Optimization Algorithms

### 1. Gradient-Based Optimization

**Objective function:**
```
Maximize: η_overall(x) = f(D₂, b₂, β₂, Z, ...)

Subject to:
  g₁(x): H_achieved ≥ H_required
  g₂(x): NPSH_available > NPSH_required
  g₃(x): D₂,min ≤ D₂ ≤ D₂,max
  ... (geometric and operational constraints)
```

**Sequential Quadratic Programming (SQP):**

Solves sequence of quadratic subproblems:
```
Minimize: ∇f(xₖ)ᵀp + 0.5pᵀBₖp

Subject to:
  ∇gᵢ(xₖ)ᵀp + gᵢ(xₖ) ≤ 0

Update:
  xₖ₊₁ = xₖ + αₖpₖ
```

**Advantages:**
- Fast convergence for smooth problems
- Handles constraints effectively
- Good for local optimization

**Limitations:**
- Requires gradient information
- May converge to local optima
- Sensitive to initial guess

### 2. Genetic Algorithms (GA)

**Process:**
1. Initialize random population of designs
2. Evaluate fitness (efficiency) for each
3. Selection (tournament, roulette wheel)
4. Crossover (combine parent designs)
5. Mutation (random modifications)
6. Repeat until convergence

**Encoding example:**
```
Individual = [D₂, b₂, β₂, Z]
            = [0.350, 0.035, 23.5, 7]
```

**Fitness function:**
```
Fitness = η_overall - Σ(penalty_i × constraint_violation_i)
```

**Parameters:**
- Population size: 50-200
- Crossover rate: 0.7-0.9
- Mutation rate: 0.01-0.1
- Generations: 50-500

**Multi-objective GA (NSGA-II):**

Simultaneously optimize multiple objectives:
```
Objectives:
  f₁(x): Maximize efficiency
  f₂(x): Minimize cost
  f₃(x): Maximize NPSH margin
  f₄(x): Minimize size/weight

Result: Pareto front of optimal trade-offs
```

### 3. Particle Swarm Optimization (PSO)

**Update equations:**
```
v_i(t+1) = w × v_i(t) + c₁r₁(p_best,i - x_i) + c₂r₂(g_best - x_i)

x_i(t+1) = x_i(t) + v_i(t+1)

where:
  v_i = velocity of particle i
  x_i = position (design variables)
  p_best,i = best position of particle i
  g_best = global best position
  w = inertia weight (0.4-0.9)
  c₁, c₂ = cognitive and social parameters (≈ 2.0)
  r₁, r₂ = random numbers [0, 1]
```

**Advantages:**
- Simple implementation
- Few parameters to tune
- Good global search capability

### 4. Surrogate-Based Optimization

**Process:**
1. Generate initial samples (Latin Hypercube, Sobol sequence)
2. Evaluate high-fidelity model (CFD, experiment)
3. Build surrogate model (Kriging, RBF, neural network)
4. Optimize surrogate (cheap evaluations)
5. Validate optimal design with high-fidelity model

**Kriging (Gaussian Process) model:**
```
ŷ(x) = μ + r(x)ᵀR⁻¹(y - μ1)

where:
  μ = mean of observations
  r(x) = correlation between x and sample points
  R = correlation matrix of sample points
  y = observations
```

**Expected Improvement criterion:**
```
EI(x) = (μ(x) - y_best) × Φ((μ(x) - y_best)/σ(x)) + σ(x) × φ((μ(x) - y_best)/σ(x))

where:
  Φ = cumulative normal distribution
  φ = probability density function
  σ(x) = prediction uncertainty
```

## Energy Savings Calculations

### 1. Annual Energy Consumption

**Basic calculation:**
```
E_annual = P_electrical × hours_per_year

P_electrical = (ρ × g × Q × H) / (η_pump × η_motor × η_VFD)

Units:
  E_annual: kWh/year
  P_electrical: kW
```

**Example:**
```
Q = 1000 GPM = 0.063 m³/s
H = 100 ft = 30.5 m
η_pump = 0.75
η_motor = 0.95
η_VFD = 0.97
hours = 6000 hr/year

P_hydraulic = 1000 × 100 / 3960 = 25.3 HP = 18.9 kW
P_electrical = 18.9 / (0.75 × 0.95 × 0.97) = 27.4 kW

E_annual = 27.4 × 6000 = 164,400 kWh/year

At $0.10/kWh:
Cost = $16,440/year
```

### 2. Efficiency Improvement Savings

**Savings from efficiency increase:**
```
Savings = Cost_baseline × (1/η_baseline - 1/η_improved)

Percent savings = (1 - η_baseline/η_improved) × 100%

Example:
  70% → 80% efficiency
  Savings = $16,440 × (1/0.70 - 1/0.80)
         = $16,440 × (1.429 - 1.250)
         = $16,440 × 0.179
         = $2,943/year

  Percent = (1 - 0.70/0.80) × 100% = 12.5%
```

### 3. VFD Energy Savings

**Power reduction with speed:**
```
P₂ = P₁ × (N₂/N₁)³

Savings = P₁ × [1 - (N₂/N₁)³] × hours × $/kWh

Example:
  P₁ = 30 HP = 22.4 kW
  N₂/N₁ = 0.80 (80% speed)
  hours = 6000
  cost = $0.10/kWh

P₂ = 22.4 × 0.80³ = 22.4 × 0.512 = 11.5 kW
Savings = (22.4 - 11.5) × 6000 × 0.10
        = $6,540/year (49% savings)
```

**VFD efficiency considerations:**
```
η_VFD ≈ 0.96-0.98 at full load
η_VFD ≈ 0.93-0.96 at 50% load

Net savings = Throttle_loss - VFD_loss
```

### 4. Pump Staging Savings

**Multiple pumps vs. single pump:**

Single 100 HP pump at 50% flow:
```
By affinity laws:
N₂/N₁ ≈ 0.50 (50% flow)
P₂ = 100 × 0.50³ = 12.5 HP

But: Efficiency drops significantly at off-BEP
Actual: ~20-30 HP due to poor efficiency
```

Two 50 HP pumps, one running at BEP:
```
Power: 50 HP at peak efficiency
Savings: 20-50% vs. single large pump at part load
```

### 5. Life Cycle Cost (LCC)

**Present value calculation:**
```
PV = FV / (1 + r)ⁿ

Annual cost present value:
PV_annual = C × [1 - (1 + r)⁻ⁿ] / r

where:
  FV = future value
  PV = present value
  r = discount rate (e.g., 0.05 for 5%)
  n = number of years
  C = annual cost
```

**Life cycle cost:**
```
LCC = C_capital + C_install + PV_energy + PV_maintenance - PV_salvage

Example (20 years, 5% discount):
  Capital: $20,000
  Install: $5,000
  Energy: $15,000/year → PV = $187,185
  Maintenance: $1,500/year → PV = $18,718
  Salvage: $2,000 → PV = $754

LCC = 20,000 + 5,000 + 187,185 + 18,718 - 754
    = $230,149
```

**Annualized cost:**
```
C_annual = LCC × [r / (1 - (1 + r)⁻ⁿ)]

Example:
C_annual = 230,149 × [0.05 / (1 - 1.05⁻²⁰)]
         = 230,149 × 0.0802
         = $18,458/year
```

### 6. Simple Payback Period

```
Payback = (Investment - Rebates) / Annual_savings

Example:
  New high-efficiency pump: $30,000
  Current pump trade-in: $5,000
  Utility rebate: $2,000
  Energy savings: $5,000/year

Payback = (30,000 - 5,000 - 2,000) / 5,000
        = 4.6 years
```

### 7. Return on Investment (ROI)

**Simple ROI:**
```
ROI = (Total_savings - Investment) / Investment × 100%

Over lifetime:
  Total savings = Annual_savings × Years

Example (10 years):
  Investment: $23,000 (net)
  Annual savings: $5,000
  Total savings: $50,000

ROI = (50,000 - 23,000) / 23,000 × 100%
    = 117%
```

**Internal Rate of Return (IRR):**

Discount rate where NPV = 0:
```
0 = -Investment + Σ(Cash_flow_t / (1 + IRR)ᵗ)

Solve iteratively for IRR
```

## Empirical Correlations

### 1. Efficiency Prediction

**Stepanoff correlation:**
```
η = 1 - 0.8 / (N_s × √(D₂ [m]))

For average manufacturing quality
```

**Anderson correlation:**
```
log(1 - η) = -0.165 × log(Q_BEP [m³/s]) - 1.024

Based on statistical analysis of pump data
```

**Karassik correlation (overall efficiency):**
```
η = 0.94 - 0.0896 / N_s^0.17

For well-designed centrifugal pumps
```

### 2. NPSH Required

**Wislicenus correlation:**
```
NPSH_req = (N_s × n / 1000)^(4/3) × (Q / n)^(2/3)

where:
  N_s = specific speed
  n = rotational speed (RPM)
  Q = flow rate
```

**Thoma cavitation parameter:**
```
σ = NPSH_req / H

Typical values:
  Radial pumps: σ = 0.05-0.15
  Mixed flow: σ = 0.15-0.35
  Axial flow: σ = 0.3-0.5
```

### 3. Power Consumption

**Bare shaft power:**
```
P_shaft = (ρ × g × Q × H) / η_pump

For water at 20°C:
P_shaft [kW] = (Q [m³/s] × H [m]) / (102 × η)
P_shaft [HP] = (Q [GPM] × H [ft] × SG) / (3960 × η)
```

**Motor power (with safety factor):**
```
P_motor = P_shaft / η_motor × SF

where:
  SF = safety factor (typically 1.15-1.25)
  η_motor = motor efficiency
```

### 4. Wear Ring Clearance

**Recommended clearances:**
```
c = K × √D

where:
  c = radial clearance (inches)
  D = diameter (inches)
  K = 0.010-0.025 (depending on service)

Examples:
  Clean water: K = 0.010-0.015
  Dirty water: K = 0.015-0.020
  Abrasive slurry: K = 0.020-0.025
```

**Leakage vs. clearance:**
```
Q_leak ∝ c³

Doubling clearance → 8× leakage!
```

## Optimization Constraints

### 1. Geometric Constraints

```
Diameter ratio: 0.3 ≤ D₁/D₂ ≤ 0.7
Width ratio: 0.03 ≤ b₂/D₂ ≤ 0.15
Blade angle: 15° ≤ β₂ ≤ 35°
Blade count: 5 ≤ Z ≤ 9 (centrifugal)
Clearance: 0.001 ≤ c ≤ 0.010 m

Flow velocity: c_m < 10 m/s (to avoid erosion)
Impeller tip speed: u₂ < 50 m/s (structural limit)
```

### 2. Performance Constraints

```
Efficiency: η ≥ η_target
Head: H_achieved ≥ H_required
Flow: Q_min ≤ Q ≤ Q_max
NPSH: NPSH_available > NPSH_required × 1.2 (safety margin)
Speed: N_min ≤ N ≤ N_max
```

### 3. Operating Constraints

```
Flow range: 0.7 × Q_BEP ≤ Q ≤ 1.2 × Q_BEP (avoid recirculation)
Speed (VFD): 0.4 × N_rated ≤ N ≤ N_rated
Parallel pumps: Same head characteristic curve
```

### 4. Material/Manufacturing Constraints

```
Minimum blade thickness: t ≥ 3 mm (castability)
Surface finish: Ra ≤ 6.3 μm (standard), Ra ≤ 1.6 μm (high efficiency)
Tolerance: IT7-IT9 (precision manufacturing)
Clearances achievable: c ≥ 0.002 D (practical limit)
```

## Case Study: Energy Audit Results

### Typical Energy Savings Opportunities

**1. Efficiency Improvement (70% → 80%)**
```
Savings: 12.5% energy reduction
Cost: New impeller or complete pump
Payback: 2-5 years
```

**2. VFD Installation (80% average speed)**
```
Savings: 49% energy reduction
Cost: VFD + installation
Payback: 1-3 years
```

**3. Impeller Trimming (10% diameter reduction)**
```
Savings: 27% power reduction (if oversized)
Cost: Machining + labor
Payback: < 1 year
```

**4. Parallel Pump Optimization**
```
Savings: 20-40% at part load
Cost: Controls + sensors
Payback: 1-2 years
```

**5. System Optimization (reduced resistance)**
```
Savings: 10-30% energy reduction
Cost: Pipe upsizing, valve replacement
Payback: 2-4 years
```

### Real-World Example

**Municipal Water Plant:**
```
Baseline:
  10 pumps × 200 HP each
  Average load: 60%
  Operating: 24/7 (8760 hr/year)
  Efficiency: 72% average
  Power: 10 × 200 × 0.746 × 0.60 / 0.72 = 1242 kW
  Annual energy: 10,884,000 kWh
  Cost: $1,088,400/year ($0.10/kWh)

Improvements:
  1. Install VFD on 5 pumps: $200,000
  2. Replace 3 old pumps with high-efficiency: $300,000
  3. Optimize staging logic: $50,000
  Total investment: $550,000

Results:
  New average efficiency: 81%
  Better load matching
  Average load: 70% of BEP (better efficiency point)
  Power: 960 kW (23% reduction)
  Annual energy: 8,410,000 kWh
  Cost: $841,000/year

Savings:
  Annual: $247,400 (23%)
  Simple payback: 2.2 years
  20-year NPV: $2.5 million
```

## Best Practice Guidelines

### 1. Design Phase

- Select pump type based on specific speed
- Design for BEP operation at most common load
- Allow 10-20% margin for fouling/aging
- Minimize clearances within manufacturing limits
- Use CFD for detailed optimization
- Consider variable speed operation
- Specify high surface finish for critical areas

### 2. Selection Phase

- Wire-to-water efficiency (include motor, VFD)
- Life cycle cost, not just first cost
- Operating range flexibility
- Maintenance accessibility
- Proven reliability in similar service
- Energy rebate eligibility

### 3. Installation Phase

- Proper alignment (< 0.002" indicator runout)
- Minimize suction piping losses
- Avoid sharp bends near pump suction
- Install flow/pressure/power monitoring
- Verify performance against curve
- Commission with full system test

### 4. Operation Phase

- Monitor efficiency trends
- Operate near BEP when possible
- Use VFD instead of throttling
- Optimize pump staging
- Schedule batch processes for off-peak
- Implement demand-based control

### 5. Maintenance Phase

- Regular vibration monitoring
- Track wear ring clearances
- Replace seals before failure
- Maintain smooth internal surfaces
- Balance after impeller work
- Trending analysis (CMMS)

## References

1. Pump Handbook, 4th Edition - Karassik et al. (2008)
2. Centrifugal Pumps: Design and Application - Val S. Lobanoff (1992)
3. ANSI/HI 9.6.7 - Effects of Liquid Viscosity on Rotodynamic Pump Performance
4. DOE Motor Challenge: Pump System Assessment Tool (PSAT)
5. Europump/Hydraulic Institute - Variable Speed Pumping Guide
6. ASHRAE Fundamentals - Fluid Flow
7. ISO 9906 - Rotodynamic Pumps - Hydraulic Performance Acceptance Tests

## Nomenclature

| Symbol | Description | Units |
|--------|-------------|-------|
| η | Efficiency | - |
| ρ | Density | kg/m³ |
| ω | Angular velocity | rad/s |
| Q | Flow rate | m³/s |
| H | Head | m |
| P | Power | W |
| N | Rotational speed | RPM |
| D | Diameter | m |
| b | Width | m |
| β | Blade angle | degrees |
| Z | Number of blades | - |
| u | Tangential velocity | m/s |
| c | Velocity (absolute frame) | m/s |
| w | Velocity (relative frame) | m/s |
| μ | Dynamic viscosity | Pa·s |
| ΔP | Pressure difference | Pa |
| NPSH | Net Positive Suction Head | m |
| BEP | Best Efficiency Point | - |
| N_s | Specific speed | - |
| Re | Reynolds number | - |
| c_m | Meridional velocity | m/s |
| c_u | Tangential component (absolute) | m/s |

## Quick Reference Formulas

**Pump power:**
```
HP = (GPM × ft × SG) / (3960 × η)
kW = (m³/s × m × ρ) / (102 × η)
```

**Efficiency change:**
```
Savings(%) = (1 - η_old/η_new) × 100
```

**VFD savings:**
```
Power_ratio = (Speed_ratio)³
```

**Specific speed:**
```
N_s = RPM × √GPM / (ft)^0.75
```

**Affinity laws:**
```
Q₂/Q₁ = N₂/N₁
H₂/H₁ = (N₂/N₁)²
P₂/P₁ = (N₂/N₁)³
```
