---
name: structural-analysis-workflow
description: "FEA workflow for pump casings and impellers under fluid loads"
category: thinking
domain: structural
complexity: advanced
dependencies: []
---

# Structural Analysis Workflow for Pump Components

This skill provides a systematic approach to performing Finite Element Analysis (FEA) on pump components, specifically casings and impellers under fluid loads.

## FEA Workflow

### 1. Geometry and Material Definition

**Geometry Preparation:**
- Import or create CAD geometry of pump casing or impeller
- Simplify geometry by removing non-structural features (small fillets, bolt holes, cosmetic features)
- Identify critical regions: high-stress areas, pressure boundaries, mounting interfaces
- Define coordinate system aligned with pump axis (Z-axis typically along shaft)

**Material Definition:**
- Assign material properties to all components
- Define isotropic or anisotropic properties as needed
- Input temperature-dependent properties if thermal effects exist
- Consider material orientation for cast or forged components

**Critical Parameters:**
- Young's Modulus (E)
- Poisson's Ratio (ν)
- Yield Strength (σ_y)
- Ultimate Tensile Strength (σ_u)
- Density (ρ)
- Thermal expansion coefficient (α)

### 2. Load Application

**Pressure Loads:**
- Apply internal pressure to casing wetted surfaces
- Use static pressure at design point or maximum operating pressure
- For impellers, apply pressure distribution from CFD analysis
- Consider pressure pulsations for fatigue analysis

**Centrifugal Forces (Impellers):**
- Apply rotational velocity to impeller body
- Centrifugal force: F_c = mω²r
- Ensure correct direction and magnitude
- Consider overspeed conditions (typically 120% of rated speed)

**Hydraulic Forces:**
- Radial thrust on impeller
- Axial thrust from pressure difference
- Tangential forces from fluid momentum change
- Unbalanced forces at off-design conditions

**Thermal Loads:**
- Temperature distribution from thermal analysis
- Differential expansion between components
- Thermal gradients during startup/shutdown

### 3. Boundary Conditions (Constraints)

**Casing Constraints:**
- Fixed supports at mounting feet
- Symmetry planes if applicable
- Pipe connections (nozzles) - typically allow thermal expansion
- Anchor points: fully constrained in all directions

**Impeller Constraints:**
- Hub constrained to shaft (cylindrical constraint)
- Contact between impeller and shaft (may use bonded or frictional contact)
- Symmetry conditions if analyzing only a sector

**Connection Types:**
- Bonded: No relative motion (welded joints)
- Frictional: Coulomb friction model (shrink fits)
- No separation: Contact that can slide but not separate
- Frictionless: Free sliding contact

### 4. Mesh Generation

**Element Types:**
- Tetrahedral (Tet10): General purpose, automatic meshing
- Hexahedral (Hex20): Higher accuracy, structured regions
- Wedge (Prism): Transition between tet and hex meshes
- Shell elements: Thin-walled sections of casings

**Mesh Quality Criteria:**
- Element aspect ratio < 5:1 (ideally < 3:1)
- Skewness < 0.85
- Jacobian ratio > 0.6
- Minimum of 3 elements through thickness
- Finer mesh at stress concentrations (fillets, holes, welds)

**Mesh Refinement:**
- Use sphere of influence around critical features
- Create local mesh controls at:
  * Fillet radii
  * Nozzle-to-shell junctions
  * Impeller eye (inlet)
  * Impeller blade-to-hub junction
  * Wear ring gaps
- Perform mesh convergence study (verify results converge as mesh density increases)

**Recommended Element Sizes:**
- Impeller blade thickness: minimum 4 elements
- Casing wall thickness: minimum 3 elements
- Fillet radius: minimum 6 elements around arc
- Global element size: 5-10% of characteristic dimension

### 5. Solver Selection

**Linear Static Analysis:**
- Use when: Deformations are small, material is elastic, no time-dependent effects
- Suitable for: Initial design verification, factor of safety calculations
- Fast computation, good for parametric studies

**Nonlinear Analysis:**
- Use when: Large deformations, material plasticity, contact nonlinearity
- Required for: Plastic collapse analysis, contact pressure distribution
- Slower computation, requires incremental loading

**Modal Analysis:**
- Determine natural frequencies and mode shapes
- Avoid resonance with pump operating speed and harmonics
- Check: 1st critical frequency > 1.5 × operating frequency

**Harmonic/Transient Analysis:**
- Pressure pulsations at blade passing frequency
- Startup and shutdown transients
- Water hammer events

**Fatigue Analysis:**
- High-cycle fatigue (> 10^4 cycles): Impellers, shafts
- Low-cycle fatigue: Thermal cycling, startup/shutdown
- Use S-N curves or strain-life methods

**Solver Settings:**
- Convergence criteria: typically 0.1% force/moment convergence
- Large deflection: ON for impellers with thin blades
- Contact settings: program-controlled for initial runs
- Stabilization: may be needed for very thin features

### 6. Results Interpretation

**Primary Results:**

**Von Mises Stress (σ_vm):**
- Equivalent stress for ductile materials
- Compare against material yield strength
- σ_vm = √[(σ₁-σ₂)² + (σ₂-σ₃)² + (σ₃-σ₁)²] / √2
- Most commonly used failure criterion

**Principal Stresses (σ₁, σ₂, σ₃):**
- Maximum and minimum normal stresses
- Critical for brittle materials
- Check maximum principal stress against ultimate strength

**Deformation:**
- Total deformation: overall displacement magnitude
- Directional deformation: check clearances (impeller-to-casing)
- Critical gaps: wear ring clearances, axial thrust bearing clearances

**Safety Factor:**
- Factor of Safety (FOS) = σ_allowable / σ_actual
- Color-coded: FOS < 1 (red/failure), FOS > 1 (safe)

**Secondary Results:**

**Shear Stress:**
- Critical at welds, keyways, shaft connections
- Maximum shear stress theory (Tresca criterion)

**Strain:**
- Elastic strain vs. plastic strain
- Strain energy density
- Equivalent plastic strain for ductile failure

**Contact Pressure:**
- Shrink fit interfaces
- Bearing surfaces
- Sealing faces

**Critical Review Checklist:**
1. Do stress concentrations align with geometric features?
2. Are boundary conditions realistic (no rigid body motion)?
3. Is the deformation pattern physically reasonable?
4. Are peak stresses at expected locations?
5. Does the solution satisfy equilibrium (sum of reactions = applied loads)?

### 7. Safety Factor Evaluation

**Design Codes and Standards:**
- ASME Section VIII Division 1: Pressure vessels (casings)
- ASME Section VIII Division 2: Higher pressure, FEA-based design
- API 610: Centrifugal pumps for petroleum industry
- ISO 13709: Centrifugal pumps for petroleum, petrochemical and natural gas industries

**Required Safety Factors:**

**Static Loads:**
- Against yield: FOS ≥ 1.5 (general machinery)
- Against yield: FOS ≥ 2.0 (pressure vessels per ASME)
- Against ultimate: FOS ≥ 2.5 (static loads)
- Against ultimate: FOS ≥ 4.0 (dynamic/impact loads)

**Fatigue Loads:**
- Infinite life design: Stress < endurance limit
- Finite life: Use cumulative damage theory (Miner's rule)
- FOS on fatigue life: typically 2.0 or higher

**Safety Factor Calculation:**

For ductile materials:
```
FOS_yield = σ_y / σ_vm
FOS_ultimate = σ_u / σ_vm
Design FOS = min(FOS_yield, FOS_ultimate)
```

For brittle materials:
```
FOS = σ_u / σ₁ (based on maximum principal stress)
```

**Material-Specific Adjustments:**
- Cast iron: Reduce allowable stress by 25%
- Welded joints: Apply joint efficiency factor (0.7-1.0)
- High temperature: Apply creep reduction factors
- Corrosive environment: Add corrosion allowance

## Applications to Pumps

### Pump Casing Pressure Containment

**Analysis Objectives:**
- Verify casing can withstand maximum allowable working pressure (MAWP)
- Check stress at nozzle-to-shell junctions
- Evaluate flange sealing capability
- Assess casing distortion effects on wear ring clearances

**Critical Regions:**
- Volute cutwater (point of maximum pressure gradient)
- Nozzle reinforcement pads
- Casing split lines (horizontal or vertical split)
- Mounting foot attachments

**Load Cases:**
- Maximum operating pressure
- Hydrostatic test pressure (typically 1.5 × design pressure)
- Pressure + thermal expansion
- Pressure + piping loads (nozzle forces and moments)

**Acceptance Criteria:**
- Von Mises stress < allowable stress
- Flange face distortion < gasket capability
- No yielding at bolt holes
- Maximum deformation < 0.5% of casing diameter

**Common Failure Modes:**
- Casing rupture at cutwater
- Flange leakage due to uneven loading
- Cracking at weld toes
- Bolt failure under pressure cycling

### Impeller Centrifugal Stress

**Analysis Objectives:**
- Calculate maximum stress in blades, shrouds, and hub
- Evaluate stress at blade-to-shroud junctions
- Check deflection to ensure clearance maintenance
- Assess fatigue life at blade passing frequency

**Critical Regions:**
- Blade trailing edge (maximum bending moment)
- Blade-to-hub fillet (stress concentration)
- Blade-to-shroud junction (closed impellers)
- Hub bore (stress concentration from keyway)

**Load Cases:**
- Centrifugal force at rated speed
- Centrifugal force at overspeed (120% of rated)
- Combined centrifugal + hydraulic pressure
- Thermal loads from fluid temperature

**Stress Components:**

**Centrifugal Stress (Hoop Stress):**
```
σ_hoop = ρω²r² (for thin disk)
where:
  ρ = material density (kg/m³)
  ω = angular velocity (rad/s)
  r = radius (m)
```

**Blade Bending Stress:**
- From pressure difference across blade
- Maximum at blade root (hub junction)
- Increases with blade length and pressure differential

**Acceptance Criteria:**
- Von Mises stress at rated speed < 0.5 × σ_y
- Von Mises stress at overspeed < 0.67 × σ_y
- Blade tip deflection < 50% of clearance gap
- No natural frequencies within ±20% of operating speed range

**Design Optimizations:**
- Increase hub diameter to reduce stress
- Taper blade thickness from hub to shroud
- Add fillets at blade-to-hub junction (r ≥ 3mm typical)
- Use materials with higher strength-to-density ratio

### Shaft Deflection

**Analysis Objectives:**
- Calculate shaft deflection under radial and axial loads
- Verify bearing alignment is maintained
- Ensure coupling alignment tolerances are met
- Check shaft critical speeds vs. operating speed

**Critical Regions:**
- Midspan between bearings (maximum deflection)
- Impeller location (affects wear ring clearance)
- Coupling location (misalignment causes vibration)
- Bearing journals (stress concentration at shoulders)

**Load Cases:**
- Radial hydraulic thrust (maximum at shutoff)
- Axial hydraulic thrust
- Impeller weight + fluid weight
- Thermal expansion

**Hydraulic Radial Thrust:**
```
F_radial = K_r × ρ × g × H × D₂ × b₂
where:
  K_r = radial thrust coefficient (0.2-0.4 for single volute)
  H = head (m)
  D₂ = impeller diameter (m)
  b₂ = impeller outlet width (m)
```

**Acceptance Criteria:**
- Maximum shaft deflection < 0.0005 × shaft span
- Deflection at impeller < 25% of wear ring clearance
- Shaft stress < 0.3 × σ_y (for fatigue resistance)
- 1st critical speed > 1.5 × maximum operating speed

**Deflection Calculation:**
- Use beam theory for preliminary estimates
- FEA for complex geometries and load distributions
- Include gyroscopic effects for high-speed pumps

### Fatigue Analysis

**Analysis Objectives:**
- Predict fatigue life in cycles or years
- Identify locations prone to crack initiation
- Evaluate stress concentration factors
- Design for infinite life or safe finite life

**Fatigue-Critical Locations:**
- Impeller blade trailing edges (pressure pulsations)
- Blade-to-hub fillets (stress concentration)
- Shaft keyway (stress concentration factor ~3)
- Shaft shoulder at bearing locations
- Weld toes on casings

**Load Cycles:**
- Startup/shutdown cycles: low-cycle fatigue
- Blade passing frequency: high-cycle fatigue
- Pressure pulsations: high-cycle fatigue
- Rotor imbalance: high-cycle fatigue

**Fatigue Analysis Methods:**

**S-N Curve Method (High-Cycle Fatigue):**
- Applicable for N > 10⁴ cycles
- Use material S-N curves (stress vs. cycles to failure)
- Apply mean stress correction (Goodman or Soderberg)
- Calculate cumulative damage using Miner's rule

**Strain-Life Method (Low-Cycle Fatigue):**
- Applicable for N < 10⁴ cycles
- Uses plastic strain range
- Coffin-Manson equation
- Required for startup/shutdown analysis

**Stress Concentration Factors:**
- Sharp corners: K_t = 2-3
- Keyways: K_t = 2-3
- Fillets: K_t = 1.5-2.5 (radius dependent)
- Threads: K_t = 2-4

**Acceptance Criteria:**
- Infinite life: alternating stress < endurance limit
- Finite life: calculated life > 2 × required life
- Factor of safety on stress amplitude ≥ 1.5
- Safety factor on life ≥ 2.0

**Mean Stress Effects:**
```
Goodman correction:
σ_a / S_e + σ_m / σ_u = 1 / FOS

where:
  σ_a = alternating stress amplitude
  σ_m = mean stress
  S_e = endurance limit
  σ_u = ultimate tensile strength
```

## Material Selection Criteria

### Pump Casing Materials

**Cast Iron (ASTM A48, A278):**
- Applications: Low-pressure water service, non-corrosive fluids
- Properties: σ_u = 150-400 MPa, low cost, good castability
- Limitations: Brittle, poor fatigue resistance, limited to <250°C
- Use for: Municipal water, HVAC, pressures < 20 bar

**Ductile Iron (ASTM A536):**
- Applications: General industrial service, moderate pressures
- Properties: σ_y = 275-550 MPa, σ_u = 400-800 MPa, good machinability
- Grades: 65-45-12, 80-55-06 (σ_u-σ_y-elongation)
- Use for: Process water, oil transfer, pressures < 40 bar

**Carbon Steel (ASTM A216 WCB):**
- Applications: High-pressure, high-temperature service
- Properties: σ_y = 250 MPa, σ_u = 485 MPa, good weldability
- Temperature range: -29°C to 400°C
- Use for: Boiler feed, power plants, pressures < 100 bar

**Stainless Steel 316 (ASTM A743 CF-8M):**
- Applications: Corrosive fluids, seawater, chemicals
- Properties: σ_y = 205 MPa, σ_u = 485 MPa, excellent corrosion resistance
- Temperature range: -196°C to 400°C
- Use for: Chemical processing, marine, food & beverage

**Bronze (ASTM B584):**
- Applications: Seawater service, small pumps
- Properties: Good corrosion resistance, moderate strength
- Use for: Marine, desalination, pump internals

### Impeller Materials

**Ductile Iron (ASTM A536):**
- Best strength-to-cost ratio
- Density: 7,100 kg/m³
- σ_y = 275 MPa minimum
- Use for: General service, water

**316 Stainless Steel:**
- Corrosion resistance
- Density: 8,000 kg/m³
- σ_y = 275 MPa (cast), 290 MPa (investment cast)
- Use for: Chemical, food, pharmaceutical

**Duplex Stainless (CD4MCu):**
- High strength and corrosion resistance
- σ_y = 450 MPa
- Expensive but long-lasting
- Use for: Seawater, aggressive chemicals

**Nickel-Aluminum Bronze (NAB):**
- Excellent cavitation resistance
- σ_y = 240 MPa
- Best for seawater
- Use for: Marine, desalination

**Titanium (Ti-6Al-4V):**
- Highest strength-to-weight ratio
- Density: 4,430 kg/m³
- σ_y = 880 MPa
- Very expensive
- Use for: High-speed, aerospace

### Shaft Materials

**Carbon Steel (AISI 1045, 4140):**
- General purpose
- σ_y = 400-650 MPa
- Good machinability
- Requires corrosion protection

**Stainless Steel 416, 17-4PH:**
- Corrosive environments
- σ_y = 520-1,170 MPa (17-4PH heat treated)
- Good fatigue resistance
- 17-4PH for high loads

**Alloy Steel (AISI 4340):**
- High-power applications
- σ_y = 860-1,380 MPa (heat treated)
- Excellent fatigue properties
- Requires protection coating

### Material Selection Process

**Step 1: Service Conditions**
- Fluid type (corrosive, abrasive, clean)
- Temperature range
- Pressure level
- Environmental exposure

**Step 2: Performance Requirements**
- Required strength (stress levels)
- Fatigue resistance (cyclic loading)
- Wear resistance (sliding contact)
- Thermal expansion compatibility

**Step 3: Manufacturing Considerations**
- Casting vs. machining from bar stock
- Weldability requirements
- Heat treatment feasibility
- Inspection requirements (NDT)

**Step 4: Cost-Benefit Analysis**
- Material cost
- Processing cost
- Expected lifetime
- Maintenance costs

**Step 5: Code Compliance**
- ASME allowable stresses
- Temperature limits
- Impact testing requirements (low temperature)
- Corrosion allowance

## Stress Concentrations

### Theoretical Stress Concentration Factor (K_t)

The stress concentration factor relates the peak local stress to the nominal stress:

```
σ_max = K_t × σ_nominal
```

### Common Stress Concentrations in Pumps

**Fillets:**
```
K_t = 1 + 2√(d/r)
where:
  d = smaller diameter
  r = fillet radius
```

Design rule: r ≥ 0.1 × d for K_t ≤ 2.0

**Shoulder (step change in diameter):**
- r/d = 0.05: K_t ≈ 2.5
- r/d = 0.10: K_t ≈ 2.0
- r/d = 0.20: K_t ≈ 1.6

**Circular hole in plate:**
- Single hole: K_t = 3.0
- Multiple holes: interaction increases K_t

**Keyway on shaft:**
- Profile keyway: K_t = 2.0-2.5
- Sled runner keyway: K_t = 2.5-3.0
- End-milled keyway: K_t = 3.0-3.5

**Threads:**
- Coarse thread: K_t = 2.2-3.0
- Fine thread: K_t = 2.8-4.0
- Stress relief groove reduces K_t

**Welds:**
- Toe of weld: K_t = 2.0-3.0
- Root of weld: K_t = 2.5-4.0
- Grinding flush reduces K_t by 30-50%

### Stress Concentration Mitigation Strategies

**Increase Fillet Radius:**
- Double the radius → reduce K_t by 20-30%
- Optimize using FEA parametric study
- Blended radii better than constant radius

**Add Relief Grooves:**
- Undercut at shoulder changes
- Reduces stress at critical section
- Common at bearing fits

**Eliminate Sharp Corners:**
- All internal corners should have radius
- Minimum radius = 2 × wall thickness
- Check manufacturability

**Relocate Transitions:**
- Move stress concentrations away from high nominal stress regions
- Place diameter changes away from high bending moment locations

**Use Material Locally:**
- Increase local thickness at stress concentration
- Add boss or reinforcement pad
- Ensure smooth transition

**Surface Treatment:**
- Shot peening introduces compressive residual stress
- Improves fatigue life by 2-4x
- Effective for shafts, blades

### Fatigue Notch Factor (K_f)

For fatigue analysis, the effective stress concentration is reduced:

```
K_f = 1 + q(K_t - 1)
where:
  q = notch sensitivity (0 to 1)
  K_t = theoretical stress concentration factor
```

Notch sensitivity depends on:
- Material ductility (lower for ductile materials)
- Notch radius (lower for larger radii)
- Material grain size

Typical values:
- Steel, r = 0.5 mm: q ≈ 0.9
- Steel, r = 2.0 mm: q ≈ 0.95
- Aluminum: q ≈ 0.7-0.8

### Verification with FEA

**Local Mesh Refinement:**
- Mesh size at notch < r/4 (fillet radius)
- Use at least 6-8 elements around fillet arc
- Verify stress convergence with mesh refinement

**Linearization:**
- Extract stress through section thickness
- Separate membrane and bending components
- Use for pressure vessel code compliance

**Submodeling:**
- Global model with coarse mesh
- Local model with fine mesh around stress concentration
- Apply displacement boundary conditions from global to local

---

## Usage Guidelines

This structural analysis workflow should be applied iteratively:

1. **Preliminary Design:** Use hand calculations and simplified FEA
2. **Design Refinement:** Detailed FEA with accurate geometry and loading
3. **Design Verification:** Final analysis with all load cases and safety factors
4. **Fatigue Assessment:** If cyclic loading is significant
5. **Documentation:** Report stress levels, safety factors, and compliance with codes

Always validate FEA results against:
- Theoretical calculations where possible
- Previous designs with field experience
- Industry standards and guidelines
- Physical testing when available

**Remember:** FEA is a tool to aid engineering judgment, not replace it. Always perform sanity checks on results and ensure boundary conditions represent actual installation and operating conditions.
