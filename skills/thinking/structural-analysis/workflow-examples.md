# Structural Analysis Workflow Examples

This document provides detailed step-by-step examples for performing FEA on pump components.

## Example 1: Pump Casing Pressure Analysis

### Problem Statement

Analyze a horizontal split case pump casing for a water supply application:
- Design pressure: 15 bar (1.5 MPa)
- Test pressure: 22.5 bar (1.5 × design pressure)
- Casing material: Ductile iron ASTM A536 Grade 65-45-12
- Casing outside diameter: 500 mm
- Wall thickness: 12 mm
- Operating temperature: 20°C

**Objective:** Verify casing can safely withstand test pressure.

### Step 1: Geometry and Material Definition

**Geometry Preparation:**
```
1. Import CAD model of casing
2. Remove features not affecting structural analysis:
   - Small bolt holes (< M10)
   - Drain plugs
   - Vent holes
   - Nameplates
3. Retain critical features:
   - Suction and discharge nozzles
   - Mounting feet
   - Casing split line
   - Wear ring housing
4. Check geometry for errors:
   - Overlapping surfaces
   - Gaps between components
   - Sliver faces
```

**Material Definition:**
```
Material: ASTM A536 Grade 65-45-12
Young's Modulus (E): 169 GPa
Poisson's Ratio (ν): 0.29
Density (ρ): 7,100 kg/m³
Yield Strength (σ_y): 310 MPa
Tensile Strength (σ_u): 448 MPa
Elongation: 12%
```

**Code Allowable Stress (ASME B31.1):**
```
S = σ_y / 1.5 = 310 / 1.5 = 207 MPa (design condition)
```

### Step 2: Load Application

**Internal Pressure:**
```
Test Pressure: 22.5 bar = 2.25 MPa
Apply to all internal wetted surfaces:
- Volute chamber
- Suction nozzle internal surface
- Discharge nozzle internal surface
- Wear ring housing
```

**Pressure Application in FEA:**
1. Select all internal faces
2. Apply normal pressure: 2.25 MPa
3. Direction: inward (compression on casing wall)

### Step 3: Boundary Conditions

**Constraints:**
```
Mounting Feet (4 locations):
- Fix bottom face of each foot in Y-direction (vertical)
- Fix one foot fully (X, Y, Z) to prevent rigid body motion
- Allow other feet to slide in X and Z (thermal expansion)

Nozzles:
- Frictionless support on nozzle flanges
- Simulates pipe connection without restraint

Split Line:
- Bonded contact (represents bolted joint with sufficient preload)
- Alternative: Include bolts and model contact
```

### Step 4: Mesh Generation

**Element Type:**
```
Primary: 10-node tetrahedral elements (SOLID187)
Refinement zones: 20-node hexahedral elements
```

**Mesh Controls:**
```
Global element size: 15 mm (3% of casing diameter)

Local refinement:
- Volute cutwater: 3 mm element size
- Nozzle-to-casing junction: 5 mm element size
- Mounting foot fillet: 4 mm element size
- Split line region: 5 mm element size

Through-thickness elements: 4 elements (12mm / 4 = 3mm per element)
```

**Mesh Quality Check:**
```
Total elements: ~185,000
Total nodes: ~312,000
Average aspect ratio: 2.1
Maximum aspect ratio: 4.3
Minimum element quality: 0.67
```

### Step 5: Solver Selection

**Analysis Type:**
```
Linear Static Analysis
- Small deformation theory
- Elastic material model
- No time-dependent effects
```

**Solver Settings:**
```
Solver: Direct Sparse
Large deflection: OFF (deformations << dimensions)
Convergence criteria: Default (0.1% force convergence)
```

### Step 6: Results Interpretation

**Von Mises Stress Distribution:**
```
Maximum von Mises stress: 158 MPa
Location: Volute cutwater region

Stress at critical locations:
- Cutwater: 158 MPa
- Discharge nozzle junction: 112 MPa
- Suction nozzle junction: 87 MPa
- Mounting foot fillet: 62 MPa
- Split line: 45 MPa
```

**Deformation:**
```
Maximum total deformation: 0.38 mm
Location: Casing top center (farthest from supports)

Deformation percentage: 0.38 / 500 = 0.076% (acceptable)
```

**Principal Stresses:**
```
Maximum principal stress (σ₁): 165 MPa (tension, inside surface)
Minimum principal stress (σ₃): -2.25 MPa (compression, outer surface)
```

### Step 7: Safety Factor Evaluation

**Factor of Safety Against Yield:**
```
FOS = σ_y / σ_max = 310 / 158 = 1.96
```

**Assessment:**
```
Required FOS per ASME VIII Div 1: ≥ 1.5 at test pressure
Actual FOS: 1.96 > 1.5 ✓ PASSES

At design pressure (15 bar):
σ_max ≈ 158 × (15/22.5) = 105 MPa
FOS = 310 / 105 = 2.95 > 2.0 ✓ PASSES
```

**Conclusion:**
```
The casing design is adequate for the specified service conditions.
Maximum stress occurs at the volute cutwater, which is expected.
Consider adding local reinforcement if pressure increases beyond 18 bar.
```

---

## Example 2: Impeller Centrifugal Stress Analysis

### Problem Statement

Analyze a closed impeller for centrifugal pump:
- Impeller diameter (D₂): 300 mm
- Impeller eye diameter: 100 mm
- Number of blades: 6
- Shroud thickness: 4 mm
- Blade thickness: 5 mm (root) to 3 mm (tip)
- Hub bore: 50 mm
- Operating speed: 2,900 rpm
- Overspeed: 3,480 rpm (120% of rated)
- Material: 316 Stainless Steel (cast)
- Design pressure: 10 bar at impeller outlet

**Objective:** Verify impeller can withstand overspeed condition.

### Step 1: Geometry and Material Definition

**Geometry:**
```
Full 360° impeller model
- Front shroud
- Back shroud
- 6 blades (equal spacing, 60° apart)
- Hub with central bore

Key dimensions:
- Outer diameter: 300 mm
- Hub diameter: 120 mm
- Blade height: 25 mm
- Fillet radius at blade-to-shroud: 2 mm
```

**Material:**
```
316 Stainless Steel (Cast, ASTM A743 CF-8M)
Young's Modulus (E): 193 GPa
Poisson's Ratio (ν): 0.27
Density (ρ): 8,000 kg/m³
Yield Strength (σ_y): 275 MPa
Tensile Strength (σ_u): 485 MPa
```

**Coordinate System:**
```
Origin: Hub center
Z-axis: Along shaft (rotation axis)
Rotation: About Z-axis
```

### Step 2: Load Application

**Rotational Velocity:**
```
Overspeed: 3,480 rpm
ω = 3,480 × 2π / 60 = 364.4 rad/s

Apply to entire impeller body
Direction: About Z-axis (coordinate system dependent)
```

**Pressure Distribution (from CFD or theoretical):**
```
Simplified approach - Linear pressure distribution:
P_eye (inlet): 0 bar (gage)
P_outlet: 10 bar (gage)

Applied to blade pressure side (concave surface)
```

**Pressure on Shrouds:**
```
Front shroud (pressure side): 10 bar at periphery, 0 at eye
Back shroud: Average pressure ≈ 5 bar
```

### Step 3: Boundary Conditions

**Hub Bore Constraint:**
```
Fixed support on hub bore cylindrical surface
- Constrains all translations (X, Y, Z)
- Constrains rotation about X and Y
- Allows rotation about Z (shaft axis)

Alternative: Use cylindrical support
- Radial direction: Fixed
- Axial direction: Fixed at one axial location
- Tangential direction: Free
```

**Symmetry (if using 60° sector):**
```
Not recommended for this analysis
Use full model to capture:
- Manufacturing variations
- Blade interaction
- Unbalanced effects
```

### Step 4: Mesh Generation

**Element Type:**
```
10-node tetrahedral elements (high-order)
Alternative: 20-node hexahedral for thicker sections
```

**Mesh Strategy:**
```
Global element size: 5 mm

Critical refinement zones:
- Blade thickness: 1.5 mm (3 elements through 5mm thickness)
- Blade-to-shroud fillet: 0.5 mm element size
- Blade-to-hub fillet: 0.8 mm element size
- Hub bore: 2 mm element size (stress concentration)
- Blade trailing edge: 2 mm

Total elements: ~425,000
Total nodes: ~650,000
```

**Mesh Quality:**
```
Aspect ratio: < 3:1 for 90% of elements
Minimum orthogonal quality: > 0.3
Through-thickness elements:
- Shroud (4mm): 3 elements
- Blade root (5mm): 4 elements
- Blade tip (3mm): 2-3 elements
```

### Step 5: Solver Selection

**Analysis Type:**
```
Static Structural (with large deflection ON)
Reason: Thin blades may have moderate deformation

Solver: Direct Sparse (preferred for contact)
OR Iterative (for very large models)
```

**Advanced Settings:**
```
Large deflection: ON
Stabilization: None (default)
Weak springs: OFF
Contact settings: N/A (no contact in this model)
```

### Step 6: Results Interpretation

**Von Mises Stress:**
```
Maximum stress locations:
1. Blade-to-hub fillet (trailing edge): 182 MPa
2. Back shroud at periphery: 145 MPa
3. Front shroud at periphery: 138 MPa
4. Hub bore (stress concentration): 95 MPa

Stress distribution:
- Increases radially from hub to periphery
- Maximum at blade trailing edge root
- Higher on back shroud (no balancing holes)
```

**Principal Stresses:**
```
At blade root fillet:
σ₁ (max principal): 195 MPa (tangential tension)
σ₂ (mid principal): 78 MPa
σ₃ (min principal): -12 MPa (radial compression)

At shroud periphery:
σ₁: 152 MPa (hoop stress)
σ₂: 68 MPa (radial stress)
```

**Deformation:**
```
Maximum total deformation: 0.24 mm
Location: Blade tip at periphery

Radial growth at periphery: 0.18 mm
Axial deflection of shroud: 0.12 mm

Critical check - Wear ring clearance:
Original clearance: 0.40 mm
Radial growth: 0.18 mm
Remaining clearance: 0.22 mm (adequate, >50% remains)
```

**Centrifugal Force Verification:**
```
Hand calculation for shroud hoop stress:
σ_hoop = ρ × ω² × r²
σ_hoop = 8,000 × (364.4)² × (0.15)²
σ_hoop = 23.9 MPa

This is membrane stress; FEA includes bending and constraints.
FEA result (145 MPa) includes:
- Hoop stress (~24 MPa)
- Bending from pressure (~35 MPa)
- Stress concentration (~2.5×)
- Secondary stresses

Reasonable agreement when factoring in all effects.
```

### Step 7: Safety Factor Evaluation

**At Overspeed (3,480 rpm):**
```
Maximum von Mises stress: 182 MPa
Material yield strength: 275 MPa

FOS = 275 / 182 = 1.51
```

**Assessment at Overspeed:**
```
Required FOS: ≥ 1.2 (overspeed is transient condition)
Actual FOS: 1.51 > 1.2 ✓ PASSES (marginal)
```

**At Rated Speed (2,900 rpm):**
```
Stress scales with ω²:
σ_rated = 182 × (2,900 / 3,480)² = 126 MPa

FOS = 275 / 126 = 2.18
Required: ≥ 1.5 for continuous operation
Actual: 2.18 > 1.5 ✓ PASSES
```

**Recommendations:**
```
1. Increase fillet radius at blade-to-hub junction to 3mm (from 2mm)
   Expected stress reduction: ~15-20%

2. Taper blade thickness more gradually to reduce bending stress

3. Consider material upgrade to 17-4PH (σ_y = 1,000 MPa heat treated)
   Would provide FOS > 5 at rated speed

4. Optimize shroud thickness (parametric study)

5. Perform fatigue analysis at blade passing frequency
```

---

## Example 3: Shaft Deflection and Critical Speed

### Problem Statement

Analyze pump shaft for a multi-stage pump:
- Shaft length between bearings: 800 mm
- Shaft diameter: 50 mm (main section)
- Impeller diameter: 200 mm
- Impeller weight: 12 kg
- Number of stages: 4
- Operating speed: 3,600 rpm
- Material: AISI 4140 steel (hardened)
- Radial hydraulic load per stage: 1,200 N

**Objectives:**
1. Calculate shaft deflection at impeller locations
2. Verify shaft critical speed > 1.5 × operating speed
3. Check stress due to bending

### Step 1: Geometry and Material Definition

**Geometry:**
```
Shaft model:
- Main section: Ø50 mm × 800 mm
- Impeller hubs: Ø70 mm × 30 mm (4 locations)
- Keyway: 12 mm wide × 8 mm deep × 50 mm long (at each impeller)
- Bearing journals: Ø45 mm × 40 mm (2 locations, at ends)

Impeller locations (from left bearing):
- Stage 1: 200 mm
- Stage 2: 350 mm
- Stage 3: 500 mm
- Stage 4: 650 mm
```

**Material:**
```
AISI 4140 Steel (Hardened and Tempered)
Young's Modulus (E): 205 GPa
Poisson's Ratio (ν): 0.29
Density (ρ): 7,850 kg/m³
Yield Strength (σ_y): 650 MPa (hardened)
Tensile Strength (σ_u): 800 MPa
```

### Step 2: Load Application

**Radial Hydraulic Thrust:**
```
At each impeller location:
Force magnitude: 1,200 N
Direction: Radial (perpendicular to shaft axis)
All forces in same direction (worst case)

Apply as point load at each impeller hub center
```

**Impeller Weight:**
```
At each impeller location:
Force magnitude: 12 kg × 9.81 m/s² = 118 N
Direction: Downward (gravity, -Y direction)

For horizontal shaft - significant
For vertical shaft - less critical for lateral deflection
```

**Self-Weight:**
```
Apply gravity to entire shaft
Acceleration: 9.81 m/s² in -Y direction
```

### Step 3: Boundary Conditions

**Bearing Supports:**
```
Left bearing (Drive end):
- Radial constraint: Fixed in X and Y
- Axial constraint: Fixed in Z (thrust bearing)
- Rotation: Free about Z-axis

Right bearing (Non-drive end):
- Radial constraint: Fixed in X and Y
- Axial constraint: Free (allows thermal expansion)
- Rotation: Free about Z-axis

Model as:
- Cylindrical support on bearing journal surfaces
- Axial constraint at one end only
```

**Impeller Hubs:**
```
No constraints (loads applied here)
Impeller mass represented by remote mass elements
- Mass: 12 kg each
- Attached to impeller hub center
```

### Step 4: Mesh Generation

**Element Type:**
```
20-node hexahedral elements (mapped mesh)
- Better for circular cross-sections
- More accurate stress calculations

Alternative: 10-node tetrahedral if complex geometry
```

**Mesh Controls:**
```
Circumferential divisions: 24 (15° increments)
Radial divisions: 6 (through radius)
Axial divisions: Variable
- Fine mesh at impeller hubs: 2 mm element length
- Coarse mesh at main shaft: 10 mm element length
- Very fine at keyways: 1 mm element length

Total elements: ~68,000
Total nodes: ~285,000
```

**Keyway Mesh:**
```
Critical for stress concentration
Use swept mesh with fine elements
Minimum 3 elements across keyway depth
Fillet radius at bottom: 0.5 mm minimum
```

### Step 5: Solver Selection

**Static Analysis:**
```
Analysis type: Static Structural
Large deflection: OFF (deflections are small)
Solver: Direct Sparse
```

**Modal Analysis (for critical speed):**
```
Analysis type: Modal
Number of modes: 10
Frequency range: 0-200 Hz (0-12,000 rpm)
```

### Step 6: Results Interpretation

**Deflection Results (Static Analysis):**
```
Maximum deflection: 0.089 mm
Location: Between Stage 2 and Stage 3 (midspan)

Deflection at impeller locations:
- Stage 1: 0.052 mm
- Stage 2: 0.071 mm
- Stage 3: 0.076 mm
- Stage 4: 0.058 mm

Allowable deflection:
- Shaft span criterion: 0.0005 × 800 = 0.40 mm
- Actual: 0.089 mm < 0.40 mm ✓ PASSES

Wear ring clearance check:
- Typical clearance: 0.3 mm
- Deflection at impeller: 0.076 mm max
- Remaining clearance: 0.224 mm (75% retained) ✓ GOOD
```

**Stress Results (Static Analysis):**
```
Maximum von Mises stress: 165 MPa
Location: Keyway at Stage 3 (combined bending + stress concentration)

Bending stress at midspan (no keyway):
σ_bending = 48 MPa

Stress concentration factor from FEA:
K_t = 165 / 48 = 3.4 (typical for keyway)
```

**Critical Speed (Modal Analysis):**
```
Mode 1 (1st bending mode): 82.5 Hz = 4,950 rpm
Mode 2 (2nd bending mode): 142.3 Hz = 8,538 rpm
Mode 3 (3rd bending mode): 187.6 Hz = 11,256 rpm

Operating speed: 3,600 rpm
Critical speed margin: 4,950 / 3,600 = 1.38
```

**Mode Shape Analysis:**
```
Mode 1: Single hump at midspan (classical 1st mode)
- Maximum deflection between Stage 2 and 3
- Both ends nodal (at bearings)

Mode 2: Two humps with node near center
- Higher frequency
- Less concern for operating speed
```

### Step 7: Safety Factor Evaluation

**Stress Assessment:**
```
Maximum stress: 165 MPa (at keyway)
Yield strength: 650 MPa

FOS = 650 / 165 = 3.94 > 1.5 ✓ PASSES

For fatigue (rotating bending):
Endurance limit (S_e) ≈ 0.5 × σ_u = 400 MPa (polished)
Keyway reduces this: S_e,actual ≈ 400 / 3.4 = 118 MPa

Alternating stress: σ_a = σ_bending = 48 MPa
Mean stress: σ_m ≈ 0 (pure bending)

FOS_fatigue = 118 / 48 = 2.46 > 2.0 ✓ PASSES
```

**Critical Speed Assessment:**
```
Required: 1st critical > 1.5 × operating speed
Required critical speed: 1.5 × 3,600 = 5,400 rpm
Actual 1st critical: 4,950 rpm

4,950 < 5,400 ✗ FAILS (marginal)

Recommendation: Increase shaft diameter to 55 mm
Expected 1st critical with Ø55: ~5,800 rpm ✓
```

**Coupling Alignment:**
```
Deflection at drive end: 0.045 mm
Typical coupling tolerance: 0.10 mm
Margin: 0.055 mm ✓ ADEQUATE
```

**Recommendations:**
```
1. Increase shaft diameter from 50 mm to 55 mm
   - Improves critical speed to ~5,800 rpm
   - Reduces bending stress by ~30%
   - Reduces deflection by ~35%

2. Optimize bearing positions (move inward by 50 mm)
   - Would increase 1st critical to ~5,500 rpm with Ø50 shaft
   - Reduces span from 800 mm to 700 mm

3. Add central steady bearing if layout permits
   - Dramatically increases critical speed
   - Reduces deflection to < 0.03 mm

4. Improve keyway design:
   - Use sled-runner profile instead of profile keyway
   - Increase fillet radius at bottom
   - Consider spline connection (no stress concentration)

5. Verify dynamic response at blade passing frequency:
   - BPF = 3,600 rpm / 60 × 6 blades = 360 Hz
   - No modal resonances in range ✓
```

---

## Example 4: Fatigue Analysis - Impeller Blade

### Problem Statement

Perform fatigue life assessment for pump impeller blade:
- Operating speed: 1,800 rpm
- Number of blades: 5
- Service: Continuous operation
- Expected life: 20 years
- Pressure pulsation amplitude: ±5% of mean pressure
- Mean pressure difference across blade: 8 bar
- Material: Duplex stainless steel CD4MCu
- Startup/shutdown cycles: 100 per year

**Objective:** Calculate fatigue life and verify > 20 years.

### Fatigue Load Calculation

**High-Cycle Fatigue (Pressure Pulsations):**
```
Operating hours per year: 8,760 hours (continuous)
Operating cycles per year: 1,800 rpm × 60 min/hr × 8,760 hr = 9.47 × 10⁸ cycles

Total cycles in 20 years: 1.89 × 10¹⁰ cycles (high-cycle fatigue)

Mean pressure: 8 bar = 0.8 MPa
Pressure amplitude: ±5% = ±0.04 MPa
```

**Low-Cycle Fatigue (Startup/Shutdown):**
```
Startup/shutdown cycles per year: 100
Total in 20 years: 2,000 cycles (low-cycle fatigue)

During startup:
- Speed ramps from 0 to 1,800 rpm
- Centrifugal stress cycles from 0 to maximum
- Thermal stress from temperature change
```

### Static Analysis Setup

**Load Case 1 - Maximum Operating Condition:**
```
Rotational speed: 1,800 rpm (ω = 188.5 rad/s)
Pressure: 8.4 bar (mean + amplitude)
Temperature: 80°C
```

**Load Case 2 - Minimum Operating Condition:**
```
Rotational speed: 1,800 rpm (ω = 188.5 rad/s)
Pressure: 7.6 bar (mean - amplitude)
Temperature: 80°C
```

**Load Case 3 - Startup:**
```
Rotational speed: 0 rpm
Pressure: 0 bar
Temperature: 20°C
```

### FEA Results

**Load Case 1 (Maximum):**
```
Maximum von Mises stress: 145 MPa
Location: Blade trailing edge at hub fillet

Principal stresses at critical location:
σ₁ = 152 MPa (tension)
σ₂ = 68 MPa
σ₃ = -8 MPa (compression)
```

**Load Case 2 (Minimum):**
```
Maximum von Mises stress: 138 MPa
Location: Same as Load Case 1
```

**Load Case 3 (Startup):**
```
Maximum von Mises stress: 12 MPa
Location: Thermal gradient region
(Residual stresses from shutdown)
```

### Fatigue Analysis - High-Cycle

**Alternating and Mean Stress:**
```
σ_max = 145 MPa (Load Case 1)
σ_min = 138 MPa (Load Case 2)

σ_mean = (σ_max + σ_min) / 2 = 141.5 MPa
σ_alternating = (σ_max - σ_min) / 2 = 3.5 MPa
```

**Material Properties:**
```
CD4MCu Duplex Stainless Steel
Ultimate tensile strength: σ_u = 655 MPa
Endurance limit (rotating beam): S_e = 0.5 × σ_u = 328 MPa (polished)

Surface finish factor (cast): k_a = 0.6
S_e,actual = 0.6 × 328 = 197 MPa
```

**Stress Concentration:**
```
From FEA at fillet:
K_t = σ_local / σ_nominal = 145 / 68 = 2.13

Notch sensitivity for duplex SS:
q = 0.90 (for r = 2 mm fillet)

Fatigue notch factor:
K_f = 1 + q(K_t - 1) = 1 + 0.90(2.13 - 1) = 2.02
```

**Modified Endurance Limit:**
```
S_e,modified = S_e,actual / K_f = 197 / 2.02 = 97.5 MPa
```

**Mean Stress Correction (Goodman):**
```
Equivalent alternating stress:
σ_eq = σ_a / (1 - σ_m / σ_u)
σ_eq = 3.5 / (1 - 141.5 / 655)
σ_eq = 3.5 / 0.784 = 4.46 MPa
```

**Assessment:**
```
σ_eq = 4.46 MPa << S_e,modified = 97.5 MPa

This indicates infinite life (> 10⁶ cycles)
Safety factor on stress: 97.5 / 4.46 = 21.9 ✓ EXCELLENT

Conclusion: High-cycle fatigue is not a concern
```

### Fatigue Analysis - Low-Cycle

**Stress Range (Startup/Shutdown):**
```
σ_max = 145 MPa (operating)
σ_min = 12 MPa (shutdown)

Δσ = 145 - 12 = 133 MPa
σ_a = Δσ / 2 = 66.5 MPa
σ_m = (145 + 12) / 2 = 78.5 MPa
```

**S-N Curve Method:**
```
For CD4MCu at N = 2,000 cycles:
From material data or: S_f = σ_u × (N / 10³)^(-0.12)
S_f = 655 × (2,000 / 1,000)^(-0.12) = 645 MPa

Apply mean stress correction (Goodman):
σ_eq = σ_a / (1 - σ_m / σ_u)
σ_eq = 66.5 / (1 - 78.5 / 655)
σ_eq = 66.5 / 0.880 = 75.6 MPa

Safety factor: 645 / 75.6 = 8.5 > 2.0 ✓ PASSES
```

**Miner's Rule (Cumulative Damage):**
```
High-cycle contribution:
n₁ = 1.89 × 10¹⁰ cycles
N₁ = ∞ (stress below endurance limit)
D₁ = n₁ / N₁ ≈ 0

Low-cycle contribution:
n₂ = 2,000 cycles
N₂ = allowable cycles at σ_eq = 75.6 MPa
From S-N curve: N₂ ≈ 1.2 × 10⁶ cycles
D₂ = 2,000 / 1,200,000 = 0.0017

Total damage: D = D₁ + D₂ ≈ 0.0017
Failure occurs when D = 1.0
Safety factor on life: 1.0 / 0.0017 = 588 ✓ EXCELLENT
```

### Conclusions and Recommendations

**Fatigue Life Assessment:**
```
✓ High-cycle fatigue: Infinite life (stress << endurance limit)
✓ Low-cycle fatigue: 588× required life
✓ Total fatigue damage: Negligible (0.17% in 20 years)

Expected fatigue life: > 11,000 years (588 × 20 years)
```

**Additional Considerations:**
```
1. Corrosion fatigue:
   - Seawater service may reduce endurance limit by 30-50%
   - Even with 50% reduction, still infinite life

2. Surface finish:
   - Assume as-cast surface (k_a = 0.6)
   - Polishing could increase endurance limit by 40%

3. Residual stresses:
   - Shot peening would add compressive stress
   - Could increase fatigue life by 2-4×

4. Stress concentrations:
   - Increase fillet radius if possible
   - Every 1 mm increase reduces K_f by ~10%
```

**Design Validation:**
```
The impeller design has adequate fatigue life for the specified service.
No design changes required from fatigue standpoint.
```

---

## Example 5: Contact Analysis - Shrink Fit

### Problem Statement

Analyze shrink fit between pump shaft and impeller hub:
- Shaft diameter: 80 mm
- Hub bore diameter (at assembly temp): 79.92 mm
- Interference: 0.08 mm (diametral)
- Hub outer diameter: 140 mm
- Hub length: 60 mm
- Shaft material: AISI 4140 steel
- Hub material: 316 SS
- Torque to transmit: 500 Nm
- Required friction coefficient: 0.12

**Objectives:**
1. Calculate contact pressure
2. Verify torque transmission capacity
3. Check stress in hub

### Step 1: Geometry and Material Definition

**Shaft:**
```
Diameter: 80 mm
Length: 200 mm (extended beyond hub for constraints)
Material: AISI 4140
E = 205 GPa, ν = 0.29
```

**Hub:**
```
Inner diameter: 79.92 mm (before interference)
Outer diameter: 140 mm
Length: 60 mm
Material: 316 SS
E = 193 GPa, ν = 0.27
```

**Interference Setup:**
```
Method 1: Offset shaft by +0.04 mm radially
Method 2: Apply thermal load to expand hub, then cool

Using Method 1 (interference offset)
```

### Step 2: Load Application

**Interference:**
```
Radial interference: 0.04 mm (half of diametral)
Applied by offsetting shaft outer surface
```

**Torque:**
```
500 Nm applied to shaft end
Transmitted through friction at interface
```

### Step 3: Boundary Conditions

**Contact:**
```
Contact region: Shaft outer surface to hub bore
Contact type: Frictional
Friction coefficient: μ = 0.12
Contact formulation: Augmented Lagrangian
```

**Constraints:**
```
Shaft: Fixed at one end (opposite of hub)
Hub: Free (constrained only through contact)
Symmetry: None (full 360° model for torque transmission)
```

### Step 4: Mesh Generation

**Contact Mesh:**
```
Element size at contact surfaces: 1 mm
Minimum 3 elements through hub wall
Contact elements: Face-to-face
```

### Step 5: Solver Selection

```
Analysis type: Nonlinear Static
Substeps: 10 (gradual interference application)
Contact solver: Augmented Lagrangian
Large deflection: ON
```

### Step 6: Results Interpretation

**Contact Pressure:**
```
Average contact pressure: 45.2 MPa
Maximum contact pressure: 52.1 MPa (at hub edges)
Minimum contact pressure: 41.8 MPa (at hub center)

Pressure distribution: Relatively uniform with edge effects
```

**Theoretical Validation:**
```
Lame's equation for thick-walled cylinder:
p = E_eff × δ / d

where:
δ = radial interference = 0.04 mm
d = diameter = 80 mm
E_eff = combined modulus

1/E_eff = (1-ν_s²)/E_s + (1-ν_h²)/E_h × (1 + (d_i/d_o)²) / (1 - (d_i/d_o)²)

For shaft (solid): (1-0.29²)/205 = 0.00447 GPa⁻¹
For hub: (1-0.27²)/193 × (1+(80/140)²)/(1-(80/140)²) = 0.00889 GPa⁻¹

1/E_eff = 0.00447 + 0.00889 = 0.01336 GPa⁻¹
E_eff = 74.9 GPa

p = 74,900 × 0.04 / 80 = 37.5 MPa

FEA: 45.2 MPa vs Theory: 37.5 MPa
Difference: 20% (FEA higher due to edge effects and refinement)
```

**Hub Stress:**
```
Maximum von Mises stress in hub: 78 MPa
Location: Inner surface at hub mid-length

Hub material yield: 275 MPa (316 SS)
FOS = 275 / 78 = 3.53 ✓ ADEQUATE
```

**Torque Capacity:**
```
Friction force per unit area: f = μ × p
f = 0.12 × 45.2 = 5.42 MPa

Contact area: A = π × d × L
A = π × 80 × 60 = 15,080 mm²

Total friction force: F = f × A
F = 5.42 × 15,080 = 81,734 N

Torque capacity: T = F × (d/2)
T = 81,734 × (80/2) = 3,269,360 Nmm = 3,269 Nm

Required torque: 500 Nm
Safety factor: 3,269 / 500 = 6.5 ✓ EXCELLENT
```

### Step 7: Assessment

**Results Summary:**
```
✓ Contact pressure: 45 MPa (adequate for torque transmission)
✓ Hub stress: 78 MPa (FOS = 3.5 > 2.0)
✓ Torque capacity: 3,269 Nm (FOS = 6.5 > 1.5)
✓ No yielding or excessive stress

Design is adequate with significant safety margins.
```

**Sensitivity Study:**

**Effect of Interference:**
```
δ = 0.06 mm: p = 34 MPa, T = 2,452 Nm (FOS = 4.9)
δ = 0.08 mm: p = 45 MPa, T = 3,269 Nm (FOS = 6.5) ← Current
δ = 0.10 mm: p = 56 MPa, T = 4,085 Nm (FOS = 8.2)

Recommendation: Current interference is adequate.
Increasing to 0.10 mm provides more margin but higher stress.
```

**Effect of Friction Coefficient:**
```
μ = 0.10: T = 2,724 Nm (FOS = 5.4) ← Minimum realistic
μ = 0.12: T = 3,269 Nm (FOS = 6.5) ← Current assumption
μ = 0.15: T = 4,086 Nm (FOS = 8.2) ← With oil retention

Even with lower friction (μ = 0.10), torque capacity is adequate.
```

---

## Summary

These examples demonstrate the complete FEA workflow for typical pump structural analyses:

1. **Casing Analysis:** Pressure containment and code compliance
2. **Impeller Analysis:** Centrifugal stress and overspeed conditions
3. **Shaft Analysis:** Deflection and critical speed calculations
4. **Fatigue Analysis:** High-cycle and low-cycle fatigue life prediction
5. **Contact Analysis:** Shrink fit interference and torque transmission

Each analysis follows the seven-step workflow:
1. Geometry and material definition
2. Load application
3. Boundary conditions
4. Mesh generation
5. Solver selection
6. Results interpretation
7. Safety factor evaluation

**Key Takeaways:**

- Always validate FEA results with hand calculations where possible
- Perform mesh convergence studies for critical analyses
- Check that deformed shape and stress distribution are physically reasonable
- Apply appropriate safety factors per design codes
- Document assumptions and load cases clearly
- Use parametric studies to optimize designs

**Common Pitfalls to Avoid:**

- Insufficient mesh refinement at stress concentrations
- Unrealistic boundary conditions (over-constrained or under-constrained)
- Ignoring stress concentration factors in fatigue analysis
- Not performing mesh convergence studies
- Neglecting to check units consistency
- Applying safety factors incorrectly (to stress vs. to load)
