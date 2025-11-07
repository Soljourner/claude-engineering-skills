# Structural Analysis Reference

This document provides theoretical background, material properties, and design code references for pump structural analysis.

## Table of Contents

1. [FEA Theory](#fea-theory)
2. [Failure Theories](#failure-theories)
3. [Material Properties](#material-properties)
4. [Design Codes and Standards](#design-codes-and-standards)
5. [Stress Analysis Formulas](#stress-analysis-formulas)
6. [Fatigue Analysis](#fatigue-analysis)
7. [Contact Mechanics](#contact-mechanics)

---

## FEA Theory

### Fundamentals of Finite Element Analysis

**Basic Principle:**

FEA discretizes a continuous domain into finite elements and approximates the solution using shape functions. The governing equation for structural analysis:

```
[K]{u} = {F}
```

where:
- [K] = global stiffness matrix
- {u} = nodal displacement vector
- {F} = applied force vector

**Strain-Displacement Relationship:**

For small deformations:
```
ε_x = ∂u/∂x
ε_y = ∂v/∂y
ε_z = ∂w/∂z
γ_xy = ∂u/∂y + ∂v/∂x
γ_yz = ∂v/∂z + ∂w/∂y
γ_zx = ∂w/∂x + ∂u/∂z
```

**Stress-Strain Relationship (Hooke's Law):**

For isotropic materials:
```
{σ} = [D]{ε}
```

where [D] is the elasticity matrix:

```
[D] = E/((1+ν)(1-2ν)) ×
[
  1-ν    ν      ν      0         0         0
  ν      1-ν    ν      0         0         0
  ν      ν      1-ν    0         0         0
  0      0      0      (1-2ν)/2  0         0
  0      0      0      0         (1-2ν)/2  0
  0      0      0      0         0         (1-2ν)/2
]
```

### Element Types

**1D Elements:**

**Beam Elements:**
- 2 nodes (linear) or 3 nodes (quadratic)
- 6 DOF per node (3 translations, 3 rotations)
- Suitable for long, slender structures (shafts, frames)
- Based on Euler-Bernoulli or Timoshenko beam theory

**2D Elements:**

**Shell Elements:**
- Thin shell: Kirchhoff-Love theory (neglect transverse shear)
- Thick shell: Mindlin-Reissner theory (include transverse shear)
- Suitable for thin-walled pressure vessels, impeller shrouds
- 5-6 DOF per node

**Plane Stress Elements:**
- σ_z = τ_xz = τ_yz = 0
- For thin plates with in-plane loading

**Plane Strain Elements:**
- ε_z = γ_xz = γ_yz = 0
- For long structures with no variation in one direction

**3D Elements:**

**Tetrahedral Elements:**
- 4 nodes (linear, SOLID72) or 10 nodes (quadratic, SOLID187)
- Automatic meshing, complex geometries
- Quadratic elements significantly more accurate

**Hexahedral Elements:**
- 8 nodes (linear, SOLID45) or 20 nodes (quadratic, SOLID186)
- More accurate than tetrahedrals with fewer elements
- Requires structured meshing
- Preferred for regular geometries

**Wedge/Prism Elements:**
- 6 nodes (linear) or 15 nodes (quadratic)
- Transition between hex and tet meshes
- Useful for radial patterns (shaft-hub interfaces)

### Shape Functions

**Linear Tetrahedral (4-node):**
```
N_i = (a_i + b_i×x + c_i×y + d_i×z) / (6V)
```
where V is element volume

**Quadratic Tetrahedral (10-node):**
- Corner nodes: Standard shape functions
- Mid-side nodes: Quadratic interpolation
- Provides better stress accuracy

### Convergence and Accuracy

**Mesh Convergence:**

Stress should converge as element size decreases:
```
σ_converged ≈ σ_n + (σ_n - σ_n-1) / (r^p - 1)
```
where:
- r = refinement ratio (typically 2)
- p = convergence rate (1-2 for linear elements, 2-3 for quadratic)

**Recommended Convergence Criterion:**
```
|σ_n - σ_n-1| / σ_n < 5%
```

**Accuracy Indicators:**
- Aspect ratio < 5:1 (3:1 preferred)
- Skewness < 0.85
- Jacobian ratio > 0.6
- Element quality > 0.3
- Minimum 3 elements through thickness

### Solution Methods

**Direct Solvers:**
- Gaussian elimination on full stiffness matrix
- Accurate, memory intensive
- Best for small to medium problems (< 500,000 DOF)
- Examples: Sparse, Frontal

**Iterative Solvers:**
- Conjugate gradient methods
- Memory efficient
- Best for large problems (> 1,000,000 DOF)
- Examples: PCG (Preconditioned Conjugate Gradient)

**Modal Analysis:**

Eigenvalue problem:
```
([K] - ω²[M]){φ} = 0
```
where:
- [M] = mass matrix
- ω = natural frequency (rad/s)
- {φ} = mode shape

Natural frequency: f = ω / (2π) Hz

---

## Failure Theories

### Von Mises Criterion (Maximum Distortion Energy)

**For ductile materials:**

```
σ_vm = √[(σ₁-σ₂)² + (σ₂-σ₃)² + (σ₃-σ₁)²] / √2
```

Alternative form using stress components:
```
σ_vm = √(σ_x² - σ_x×σ_y + σ_y² + 3τ_xy²)  (2D)

σ_vm = √[(σ_x-σ_y)² + (σ_y-σ_z)² + (σ_z-σ_x)² + 6(τ_xy² + τ_yz² + τ_zx²)] / √2  (3D)
```

**Failure criterion:**
```
σ_vm ≥ σ_y  →  Yielding begins
```

**Safety factor:**
```
FOS = σ_y / σ_vm
```

**Applicability:** Best for ductile materials (steel, aluminum, ductile iron)

### Maximum Principal Stress (Rankine)

**For brittle materials:**

```
Failure when: σ₁ ≥ σ_u  (tension)
           or: σ₃ ≤ -σ_uc (compression)
```

where:
- σ₁ = maximum principal stress (most tensile)
- σ₃ = minimum principal stress (most compressive)
- σ_u = ultimate tensile strength
- σ_uc = ultimate compressive strength

**Safety factor:**
```
FOS = σ_u / σ₁  (if σ₁ > 0)
FOS = σ_uc / |σ₃|  (if σ₃ < 0)
```

**Applicability:** Cast iron, ceramics, concrete

### Maximum Shear Stress (Tresca)

**For ductile materials (conservative):**

```
τ_max = (σ₁ - σ₃) / 2
```

**Failure criterion:**
```
τ_max ≥ σ_y / 2  →  Yielding begins
```

**Safety factor:**
```
FOS = σ_y / (2τ_max)
```

**Relationship to Von Mises:**
- Tresca is more conservative (predicts failure at 15% lower stress)
- Tresca: inscribed hexagon in π-plane
- Von Mises: circle in π-plane

### Mohr-Coulomb Criterion

**For materials with different tensile and compressive strengths:**

```
σ₁/σ_t - σ₃/σ_c = 1
```

where:
- σ_t = tensile strength
- σ_c = compressive strength

**Applicability:** Soils, rocks, some composites

### Comparison of Failure Theories

| Material | Recommended Theory | Conservatism |
|----------|-------------------|--------------|
| Steel | Von Mises | Moderate |
| Stainless Steel | Von Mises | Moderate |
| Aluminum | Von Mises | Moderate |
| Ductile Iron | Von Mises or Tresca | Tresca more conservative |
| Cast Iron | Maximum Principal Stress | - |
| Bronze | Von Mises | Moderate |
| Composites | Tsai-Wu or Tsai-Hill | Varies |

---

## Material Properties

### Common Pump Materials

#### Ferrous Alloys

**Cast Iron (ASTM A48, Grade 30):**
```
Young's Modulus (E): 90-115 GPa (typical: 100 GPa)
Poisson's Ratio (ν): 0.26
Density (ρ): 7,200 kg/m³
Tensile Strength (σ_u): 210 MPa
Compressive Strength: 620 MPa (3× tensile)
Brinell Hardness: 150-200 HB
Thermal Expansion (α): 10.8 × 10⁻⁶ /°C
Thermal Conductivity (k): 50 W/(m·K)
Max Temperature: 250°C

Notes:
- Brittle, no distinct yield point
- Good damping properties
- Low cost
- Use maximum principal stress theory
```

**Ductile Iron (ASTM A536, Grade 65-45-12):**
```
Young's Modulus (E): 169 GPa
Poisson's Ratio (ν): 0.29
Density (ρ): 7,100 kg/m³
Yield Strength (σ_y): 310 MPa (0.2% offset)
Tensile Strength (σ_u): 448 MPa
Elongation: 12%
Brinell Hardness: 140-180 HB
Thermal Expansion (α): 11.0 × 10⁻⁶ /°C
Thermal Conductivity (k): 36 W/(m·K)
Max Temperature: 370°C

Endurance Limit: S_e ≈ 0.4 × σ_u = 179 MPa (polished)

Common Grades:
- 60-40-18: σ_y=276 MPa, σ_u=414 MPa, El=18%
- 65-45-12: σ_y=310 MPa, σ_u=448 MPa, El=12%
- 80-55-06: σ_y=379 MPa, σ_u=552 MPa, El=6%
- 100-70-03: σ_y=483 MPa, σ_u=689 MPa, El=3%
```

**Carbon Steel, Cast (ASTM A216, Grade WCB):**
```
Young's Modulus (E): 200 GPa
Poisson's Ratio (ν): 0.29
Density (ρ): 7,850 kg/m³
Yield Strength (σ_y): 250 MPa
Tensile Strength (σ_u): 485 MPa
Elongation: 22%
Brinell Hardness: 135-175 HB
Thermal Expansion (α): 11.7 × 10⁻⁶ /°C
Thermal Conductivity (k): 52 W/(m·K)
Max Temperature: 400°C

Endurance Limit: S_e ≈ 0.5 × σ_u = 243 MPa (polished, unnotched)

ASME Allowable Stress (Section VIII):
- At 20°C: 138 MPa
- At 100°C: 138 MPa
- At 200°C: 138 MPa
- At 300°C: 135 MPa
```

**Carbon Steel, Wrought (AISI 1045):**
```
Young's Modulus (E): 205 GPa
Poisson's Ratio (ν): 0.29
Density (ρ): 7,850 kg/m³
Yield Strength (σ_y): 310 MPa (as-rolled)
                      530 MPa (Q&T at 205°C)
Tensile Strength (σ_u): 565 MPa (as-rolled)
                        625 MPa (Q&T at 205°C)
Elongation: 16% (as-rolled)
Brinell Hardness: 163 HB (as-rolled), 229 HB (Q&T)
Thermal Expansion (α): 11.5 × 10⁻⁶ /°C
Max Temperature: 425°C

Endurance Limit: S_e ≈ 0.5 × σ_u = 280 MPa (polished)
                                   195 MPa (machined)
```

**Alloy Steel (AISI 4140):**
```
Young's Modulus (E): 205 GPa
Poisson's Ratio (ν): 0.29
Density (ρ): 7,850 kg/m³

Heat Treatment Dependent:
Annealed:
- Yield Strength: 417 MPa
- Tensile Strength: 655 MPa
- Elongation: 25%

Normalized:
- Yield Strength: 655 MPa
- Tensile Strength: 1,020 MPa
- Elongation: 18%

Q&T (Oil quench, 205°C temper):
- Yield Strength: 1,570 MPa
- Tensile Strength: 1,770 MPa
- Elongation: 9%

Typical for shafts (Q&T at 400°C):
- Yield Strength: 910 MPa
- Tensile Strength: 1,080 MPa
- Elongation: 14%
- Brinell Hardness: 321 HB

Endurance Limit: S_e ≈ 0.5 × σ_u = 540 MPa (polished, Q&T)
Thermal Expansion (α): 12.3 × 10⁻⁶ /°C
```

#### Stainless Steels

**316 Stainless Steel, Cast (ASTM A743, Grade CF-8M):**
```
Young's Modulus (E): 193 GPa
Poisson's Ratio (ν): 0.27
Density (ρ): 8,000 kg/m³
Yield Strength (σ_y): 205 MPa (annealed)
Tensile Strength (σ_u): 485 MPa
Elongation: 30%
Brinell Hardness: 140-217 HB
Thermal Expansion (α): 16.5 × 10⁻⁶ /°C
Thermal Conductivity (k): 16 W/(m·K)
Max Temperature: 425°C (continuous)
                 870°C (intermittent)

Endurance Limit: S_e ≈ 0.4 × σ_u = 194 MPa (cast surface)

ASME Allowable Stress:
- At 20°C: 115 MPa
- At 100°C: 100 MPa
- At 200°C: 92 MPa

Corrosion Resistance: Excellent in most environments
```

**316 Stainless Steel, Wrought (ASTM A276):**
```
Young's Modulus (E): 193 GPa
Poisson's Ratio (ν): 0.27
Density (ρ): 8,000 kg/m³
Yield Strength (σ_y): 290 MPa (annealed)
Tensile Strength (σ_u): 580 MPa
Elongation: 40%
Brinell Hardness: 217 HB
Thermal Expansion (α): 16.5 × 10⁻⁶ /°C

Cold Worked (1/4 hard):
- Yield Strength: 515 MPa
- Tensile Strength: 860 MPa

Endurance Limit: S_e ≈ 0.5 × σ_u = 290 MPa (polished)
```

**17-4 PH Stainless Steel (ASTM A564, Condition H900):**
```
Young's Modulus (E): 196 GPa
Poisson's Ratio (ν): 0.27
Density (ρ): 7,800 kg/m³
Yield Strength (σ_y): 1,170 MPa
Tensile Strength (σ_u): 1,310 MPa
Elongation: 10%
Rockwell Hardness: 42-46 HRC
Thermal Expansion (α): 10.8 × 10⁻⁶ /°C
Max Temperature: 315°C

Other Conditions:
- H1025: σ_y=1,000 MPa, σ_u=1,070 MPa, El=12%
- H1150: σ_y=795 MPa, σ_u=930 MPa, El=16%

Endurance Limit: S_e ≈ 0.45 × σ_u = 590 MPa (H900, polished)

Applications: High-strength shafts, impellers for corrosive service
```

**Duplex Stainless Steel, Cast (ASTM A890, Grade 4A / CD4MCu):**
```
Young's Modulus (E): 185 GPa
Poisson's Ratio (ν): 0.28
Density (ρ): 7,800 kg/m³
Yield Strength (σ_y): 450 MPa (0.2% offset)
Tensile Strength (σ_u): 655 MPa
Elongation: 15%
Brinell Hardness: 250 HB
Thermal Expansion (α): 13.0 × 10⁻⁶ /°C
Thermal Conductivity (k): 15 W/(m·K)
Max Temperature: 315°C

Composition: Fe-25Cr-5Ni-3Mo-2Cu
Microstructure: ~50% ferrite, ~50% austenite

Endurance Limit: S_e ≈ 0.5 × σ_u = 328 MPa

Advantages:
- 2× strength of 316 SS
- Excellent corrosion resistance (especially pitting, SCC)
- Good erosion resistance
- Lower thermal expansion than austenitic SS

Applications: Seawater service, high-pressure impellers
```

#### Non-Ferrous Alloys

**Aluminum Bronze (ASTM B148, Alloy C95400):**
```
Young's Modulus (E): 120 GPa
Poisson's Ratio (ν): 0.30
Density (ρ): 7,600 kg/m³
Yield Strength (σ_y): 205 MPa
Tensile Strength (σ_u): 585 MPa
Elongation: 18%
Brinell Hardness: 170 HB
Thermal Expansion (α): 16.2 × 10⁻⁶ /°C
Thermal Conductivity (k): 63 W/(m·K)
Max Temperature: 275°C

Composition: Cu-11Al-4Fe-4Ni
Endurance Limit: S_e ≈ 0.3 × σ_u = 176 MPa

Advantages:
- Excellent corrosion resistance (seawater)
- Good cavitation resistance
- Good wear resistance
```

**Nickel-Aluminum Bronze (NAB) (ASTM B148, Alloy C95800):**
```
Young's Modulus (E): 130 GPa
Poisson's Ratio (ν): 0.32
Density (ρ): 7,650 kg/m³
Yield Strength (σ_y): 240 MPa
Tensile Strength (σ_u): 620 MPa
Elongation: 15%
Brinell Hardness: 180 HB
Thermal Expansion (α): 16.2 × 10⁻⁶ /°C
Max Temperature: 300°C

Composition: Cu-9Al-5Ni-5Fe-1Mn
Endurance Limit: S_e ≈ 0.35 × σ_u = 217 MPa

Advantages:
- Best cavitation resistance of common pump alloys
- Excellent seawater corrosion resistance
- Good strength and toughness
- Self-healing oxide film

Applications: Marine pumps, desalination, seawater injection
```

**Titanium Alloy (Ti-6Al-4V, ASTM B367, Grade 5):**
```
Young's Modulus (E): 114 GPa
Poisson's Ratio (ν): 0.34
Density (ρ): 4,430 kg/m³
Yield Strength (σ_y): 880 MPa (annealed)
Tensile Strength (σ_u): 950 MPa
Elongation: 14%
Rockwell Hardness: 36 HRC
Thermal Expansion (α): 8.6 × 10⁻⁶ /°C
Thermal Conductivity (k): 7.4 W/(m·K)
Max Temperature: 400°C (continuous)

Endurance Limit: S_e ≈ 0.55 × σ_u = 520 MPa

Advantages:
- Highest strength-to-weight ratio (σ/ρ = 198 kN·m/kg)
- Excellent corrosion resistance (forming passive TiO₂ film)
- Low thermal expansion
- Non-magnetic

Disadvantages:
- Very expensive (~$30-50/kg vs $2-5/kg for steel)
- Difficult to machine
- Poor galling resistance (needs surface treatment)

Applications: Aerospace, military, high-speed turbomachinery
```

### Temperature-Dependent Properties

Most structural properties degrade with temperature. For carbon steel:

| Temperature (°C) | E (GPa) | σ_y Reduction | σ_u Reduction |
|------------------|---------|---------------|---------------|
| 20 | 200 | 1.00 | 1.00 |
| 100 | 197 | 0.95 | 0.97 |
| 200 | 193 | 0.85 | 0.90 |
| 300 | 186 | 0.75 | 0.82 |
| 400 | 177 | 0.65 | 0.73 |
| 500 | 165 | 0.50 | 0.60 |

**Creep Effects:**

For temperatures above 0.4 × T_melt (absolute), consider creep:
- Carbon steel: Above ~370°C
- Stainless steel: Above ~425°C
- Aluminum: Above ~150°C

---

## Design Codes and Standards

### ASME Section VIII - Pressure Vessels

**Division 1: Design by Formula**

**Allowable Stress:**
```
S = minimum of:
- σ_u / 3.5
- σ_y / 1.5
```

For cast materials, multiply by quality factor:
- Radiographed castings: 1.0
- Non-radiographed: 0.8

**Cylindrical Shell (Internal Pressure):**

Circumferential (hoop) stress:
```
t = PR / (SE - 0.6P)
```
where:
- t = required thickness (mm)
- P = internal pressure (MPa)
- R = inside radius (mm)
- S = allowable stress (MPa)
- E = joint efficiency (0.7-1.0)

Longitudinal stress (half of hoop stress):
```
t = PR / (2SE + 0.4P)
```

**Spherical Shell:**
```
t = PR / (2SE - 0.2P)
```

**Ellipsoidal Head:**
```
t = PD / (2SE - 0.2P)
```
where D = major diameter

**Torispherical Head (ASME F&D):**
```
t = 0.885PL / (SE - 0.1P)
```
where:
- L = crown radius (≈ D for F&D head)
- M = stress intensification factor

**Nozzle Reinforcement:**

Required reinforcement area:
```
A = d × t × F
```
where:
- d = nozzle opening diameter
- t = shell thickness
- F = correction factor (1.0 for no internal pressure)

Available reinforcement:
- Shell in excess of required thickness
- Nozzle wall in excess of required thickness
- Welds
- Reinforcement pad

**Hydrostatic Test Pressure:**
```
P_test = 1.5 × P_design × (S_design / S_test)
```

Typically: P_test = 1.5 × P_design (if temperature correction negligible)

**Division 2: Design by Analysis**

More rigorous, allows FEA-based design.

**Stress Categories:**

1. **Primary Stresses (P):**
   - Produced by mechanical loads
   - Not self-limiting (can cause failure)

   Types:
   - P_m: Primary membrane (average through thickness)
   - P_b: Primary bending
   - P_L: Primary local membrane

2. **Secondary Stresses (Q):**
   - Self-limiting (caused by constraint)
   - Examples: thermal expansion, stress concentration
   - Can cause fatigue but not immediate failure

3. **Peak Stresses (F):**
   - Local stresses (notches, welds)
   - Only significant for fatigue

**Allowable Stresses (Div 2):**

```
P_m ≤ S_m (design membrane stress intensity)
P_L ≤ 1.5 × S_m
P_L + P_b ≤ 1.5 × S_m
P_L + P_b + Q ≤ 3 × S_m
P_L + P_b + Q + F: Check fatigue
```

where:
```
S_m = minimum of (σ_u / 3.5, σ_y / 1.5)  [Same as Div 1]
```

**Stress Linearization:**

For FEA results, linearize stress through thickness:
```
Membrane stress: σ_m = (1/t) ∫ σ dt
Bending stress: σ_b = (6/t²) ∫ σ × (z - t/2) dt
Peak stress: σ_p = σ_total - σ_m - σ_b
```

**Fatigue Evaluation (Div 2, Part 5):**

1. Calculate stress amplitude: S_alt
2. Determine mean stress: S_mean
3. Apply stress concentration factors
4. Obtain allowable cycles from design fatigue curves
5. Calculate cumulative damage (Miner's rule)

**Design Fatigue Curves:**

ASME provides S-N curves for various materials with safety factors:
- Factor of 2 on stress
- Factor of 20 on cycles (or 10 on cycles for N < 10⁴)

### API 610 - Centrifugal Pumps

**Shaft Deflection Limits:**

```
Deflection at impeller ≤ 0.005 inches (0.127 mm)
OR
Deflection ≤ 50% of wear ring clearance
```

**Shaft Critical Speed:**

```
1st critical speed ≥ 1.4 × maximum continuous speed
OR
1st critical speed ≥ 1.2 × trip speed
```

For variable speed:
```
Operating speed range shall not be within ±20% of any critical speed
```

**Shaft Design:**

Lateral analysis required considering:
- Hydraulic loads (maximum expected)
- Rotor weight
- Coupling forces
- Thermal effects

**Minimum Shaft Factor of Safety:**

```
FOS ≥ 2.0 on yield strength (static loads)
FOS ≥ 1.5 on endurance limit (fatigue)
```

**Casing Design Pressure:**

```
Maximum Allowable Working Pressure (MAWP):
MAWP ≥ 1.5 × shutoff pressure at maximum impeller diameter

OR

MAWP ≥ discharge pressure + 1.5 × rated differential pressure
```

**Hydrostatic Test:**
```
Test pressure = 1.5 × MAWP (minimum)
Hold for 30 minutes, no leakage permitted
```

### ISO 13709 (API 610 equivalent)

Similar requirements to API 610, adopted internationally.

### ASME B31.1 - Power Piping

Used for piping connected to pumps (nozzle loads).

**Allowable Stress:**
```
S = minimum of:
- σ_u / 3.5
- σ_y / 1.5
- creep stress for temperature
```

**Pressure Design:**
```
t = PD / (2SE + Py)
```
Plus corrosion allowance (typically 3 mm)

**Nozzle Loads:**

Maximum allowable nozzle loads transmitted to pump (typical):
```
F_axial: 500-5,000 N (depends on nozzle size)
F_shear: 500-5,000 N
M_bending: 200-2,000 Nm
M_torsion: 100-1,000 Nm
```

Verify pump manufacturer's limits.

### DIN/EN Standards (European)

**EN 13445 - Unfired Pressure Vessels:**
Similar to ASME Section VIII, European equivalent

**DIN 24960 - Centrifugal Pumps:**
German standard for pump design

---

## Stress Analysis Formulas

### Thin-Walled Pressure Vessels

**Criteria:** t/r < 0.1 (t = thickness, r = radius)

**Cylindrical Shell:**
```
Hoop stress: σ_θ = Pr / t
Axial stress: σ_a = Pr / (2t)
Radial stress: σ_r ≈ 0 (thin-wall assumption)

Maximum stress: σ_max = Pr / t (hoop)
```

**Spherical Shell:**
```
σ = Pr / (2t)  (uniform in all directions)
```

### Thick-Walled Cylinders (Lamé Equations)

**Criteria:** t/r ≥ 0.1

Internal pressure P_i, external pressure P_o:

**Radial stress:**
```
σ_r = [P_i×r_i² - P_o×r_o² - (P_i - P_o)×r_i²×r_o²/r²] / (r_o² - r_i²)
```

**Tangential (hoop) stress:**
```
σ_θ = [P_i×r_i² - P_o×r_o² + (P_i - P_o)×r_i²×r_o²/r²] / (r_o² - r_i²)
```

At inner surface (r = r_i):
```
σ_r = -P_i (compressive)
σ_θ = [P_i×(r_o² + r_i²) - 2×P_o×r_o²] / (r_o² - r_i²)
```

For internal pressure only (P_o = 0):
```
σ_θ,inner = P_i × (r_o² + r_i²) / (r_o² - r_i²)

σ_vm,inner ≈ σ_θ,inner  (since σ_r = -P_i << σ_θ)
```

### Rotating Disks and Impellers

**Solid Disk:**
```
σ_r = ρω²[(3+ν)/8] × (r_o² - r²)
σ_θ = ρω²[(3+ν)/8]×r_o² - [(1+3ν)/8]×r²

Maximum stress (at center, r=0):
σ_r = σ_θ = ρω²(3+ν)r_o² / 8
```

**Disk with Central Hole:**
```
At inner radius (r = r_i):
σ_θ = ρω² × [(3+ν)/4] × [r_i² + (1-ν)/(3+ν)×r_o²]
```

**Approximate Hoop Stress (impeller shroud):**
```
σ_hoop ≈ ρ × ω² × r²
```

where:
- ρ = density (kg/m³)
- ω = angular velocity (rad/s)
- r = radius (m)

### Beam Bending

**Euler-Bernoulli Beam Theory:**

Bending stress:
```
σ = -M×y / I
```
where:
- M = bending moment
- y = distance from neutral axis
- I = second moment of area

Deflection:
```
d²y/dx² = M / (EI)
```

**Maximum deflection (simply supported, central load F):**
```
y_max = F×L³ / (48EI)
```

**Maximum deflection (simply supported, uniformly distributed load w):**
```
y_max = 5wL⁴ / (384EI)
```

**Circular shaft bending:**
```
I = πd⁴ / 64
Z = πd³ / 32

σ_max = M / Z = 32M / (πd³)
```

### Shaft Torsion

**Circular shaft:**
```
τ_max = Tr / J
```
where:
- T = torque
- r = outer radius
- J = polar moment of inertia = πd⁴ / 32

Shear stress:
```
τ_max = 16T / (πd³)
```

Angle of twist:
```
θ = TL / (GJ)
```
where G = shear modulus = E / [2(1+ν)]

**Power transmission:**
```
T = P / ω = 9,550 × P / n
```
where:
- P = power (kW)
- ω = angular velocity (rad/s)
- n = speed (rpm)
- T = torque (Nm)

### Contact Stress (Hertzian)

**Cylinder on cylinder (parallel axes):**
```
p_max = √[F(E*) / (πL)] × √[1/r₁ + 1/r₂]
```
where:
- F = normal force
- L = contact length
- E* = effective modulus = 2E₁E₂ / [E₁(1-ν₂²) + E₂(1-ν₁²)]

**Sphere on flat surface:**
```
p_max = [3F / (2πa²)]
where: a³ = 3FR / (4E*)
```

---

## Fatigue Analysis

### S-N Curve (Stress-Life) Method

**For high-cycle fatigue (N > 10⁴ cycles):**

**S-N Relationship:**
```
S = σ'_f × (2N)^b
```
where:
- S = stress amplitude
- N = cycles to failure
- σ'_f = fatigue strength coefficient
- b = fatigue strength exponent (-0.05 to -0.12)

**Simplified (Basquin's equation):**
```
S × N^k = C
```
where k = 0.08-0.12 (typically 0.1 for steels)

**Endurance Limit (S_e):**

For ferrous metals, exists at N ≈ 10⁶ cycles:
```
S_e = 0.5 × σ_u  (wrought steel, polished)
S_e = 0.4 × σ_u  (cast steel)
S_e = 0.3 × σ_u  (cast iron)
```

**Endurance Limit Modifying Factors:**

```
S_e = k_a × k_b × k_c × k_d × k_e × S'_e
```

where:
- k_a = surface finish factor
- k_b = size factor
- k_c = load type factor
- k_d = temperature factor
- k_e = reliability factor
- S'_e = 0.5 × σ_u (base endurance limit)

**Surface Finish Factor (k_a):**
```
k_a = a × (σ_u)^b  [σ_u in MPa]
```

| Surface Finish | a | b |
|----------------|------|-------|
| Ground | 1.58 | -0.085 |
| Machined | 4.51 | -0.265 |
| Cold-rolled | 2.70 | -0.265 |
| Hot-rolled | 57.7 | -0.718 |
| As-forged | 272 | -0.995 |
| As-cast | 0.50 | 0 |

**Size Factor (k_b):**
```
For d ≤ 8 mm: k_b = 1.0
For 8 < d < 250 mm: k_b = 0.85
For d ≥ 250 mm: k_b = 0.75

More accurate: k_b = (d / 7.62)^(-0.107)  [d in mm]
```

**Load Type Factor (k_c):**
```
Bending: k_c = 1.0
Axial: k_c = 0.85
Torsion: k_c = 0.59
```

**Temperature Factor (k_d):**
```
T < 450°C: k_d = 1.0
450-550°C: k_d = 1 - 0.0058(T - 450)  [T in °C]
```

**Reliability Factor (k_e):**

| Reliability | k_e |
|-------------|------|
| 50% | 1.000 |
| 90% | 0.897 |
| 95% | 0.868 |
| 99% | 0.814 |
| 99.9% | 0.753 |

### Mean Stress Effects

**Goodman Relation (most common):**
```
σ_a / S_e + σ_m / σ_u = 1
```

Solving for allowable alternating stress:
```
σ_a,allow = S_e × (1 - σ_m / σ_u)
```

**Gerber Relation (less conservative):**
```
σ_a / S_e + (σ_m / σ_u)² = 1
```

**Soderberg Relation (most conservative):**
```
σ_a / S_e + σ_m / σ_y = 1
```

**Modified Goodman (with safety factor):**
```
σ_a / S_e + σ_m / σ_u = 1 / n
```
where n = desired factor of safety

### Notch Effects on Fatigue

**Fatigue Notch Factor:**
```
K_f = 1 + q(K_t - 1)
```
where:
- K_t = theoretical stress concentration factor
- q = notch sensitivity (0 to 1)

**Notch Sensitivity (q):**

Depends on material and notch radius:
```
q = 1 / [1 + a/√r]
```

where:
- r = notch radius (mm)
- a = material constant

For steels:
```
a ≈ 0.025√(2070/σ_u)  [σ_u in MPa]
```

Typical values:
- Steel (σ_u = 500 MPa), r = 1 mm: q ≈ 0.90
- Steel (σ_u = 1000 MPa), r = 1 mm: q ≈ 0.95
- Aluminum, r = 1 mm: q ≈ 0.80

**Reduced Endurance Limit:**
```
S_e,notch = S_e / K_f
```

### Cumulative Damage (Miner's Rule)

For varying stress amplitudes:
```
D = Σ(n_i / N_i)
```
where:
- n_i = actual number of cycles at stress level i
- N_i = cycles to failure at stress level i
- D = cumulative damage

Failure predicted when D ≥ 1.0

Design typically requires: D ≤ 0.5 (FOS = 2 on life)

### Strain-Life (ε-N) Method

**For low-cycle fatigue (N < 10⁴ cycles):**

**Coffin-Manson Equation:**
```
Δε_total / 2 = Δε_e / 2 + Δε_p / 2

Δε_e / 2 = (σ'_f / E) × (2N)^b  (elastic component)

Δε_p / 2 = ε'_f × (2N)^c  (plastic component)
```

where:
- σ'_f = fatigue strength coefficient ≈ σ_u
- b = fatigue strength exponent ≈ -0.12
- ε'_f = fatigue ductility coefficient ≈ ε_f (true fracture strain)
- c = fatigue ductility exponent ≈ -0.6

**Transition Life (N_t):**

Point where elastic and plastic components are equal:
```
N_t = (ε'_f × E / σ'_f)^(1/(b-c))
```

For N < N_t: plastic strain dominates (low-cycle)
For N > N_t: elastic strain dominates (high-cycle)

### Fatigue Safety Factor

**On stress:**
```
n = S_e / σ_a,eq
```
where σ_a,eq includes mean stress correction

**On life:**
```
n = N_allow / N_required
```

Typically require:
- n ≥ 1.5 on stress
- n ≥ 2.0 on life

---

## Contact Mechanics

### Hertzian Contact Theory

**Two Cylinders in Contact (parallel axes):**

Contact half-width:
```
b = √[4FR / (πLE*)] × √[1/(1/r₁ + 1/r₂)]
```

Maximum contact pressure:
```
p_max = 2F / (πbL)
```

where:
```
1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂  (effective modulus)
1/r = 1/r₁ + 1/r₂  (effective radius)
```

**Two Spheres in Contact:**

Contact radius:
```
a = [3FR / (4E*)]^(1/3)
```
where:
```
1/R = 1/R₁ + 1/R₂
```

Maximum contact pressure:
```
p_max = 3F / (2πa²)
```

**Maximum Shear Stress:**

Occurs below surface at depth z ≈ 0.5a:
```
τ_max ≈ 0.31 × p_max
```

### Shrink Fit Analysis

**Interface Pressure:**

For hub on shaft:
```
p = δ / [d × (1/E_hub × C_hub + 1/E_shaft × C_shaft)]
```

where:
```
C_hub = (d_o² + d_i²) / (d_o² - d_i²) + ν_hub
C_shaft = (1 + ν_shaft)  [for solid shaft]

δ = radial interference
d = nominal diameter
d_i = hub inner diameter
d_o = hub outer diameter
```

**Torque Capacity:**
```
T = μ × p × π × d × L × (d/2)
```
where:
- μ = coefficient of friction
- L = hub length

**Required Interference:**

To transmit torque T:
```
δ = T × [(d × E_hub) / (2πμL)] × [(d_o² - d_i²) / (d_o² + d_i²)]
```

Plus safety factor (typically 1.5-2.0)

**Stress in Hub:**

At inner surface:
```
σ_θ = p × (d_o² + d_i²) / (d_o² - d_i²)
σ_r = -p
```

Von Mises stress:
```
σ_vm ≈ σ_θ - σ_r = p × (d_o² + d_i²) / (d_o² - d_i²) + p
σ_vm = p × [2d_o² / (d_o² - d_i²)]
```

Must satisfy: σ_vm < σ_y

---

## Quick Reference Tables

### Typical Safety Factors

| Application | Yield | Ultimate | Fatigue |
|-------------|--------|----------|---------|
| Pump casing (ASME VIII) | 1.5 | 3.5 | - |
| Impeller (rated speed) | 1.5 | 2.5 | 2.0 |
| Impeller (overspeed 120%) | 1.2 | 2.0 | - |
| Shaft (continuous) | 2.0 | - | 2.0 |
| Bolts (pressure) | 2.5 | 4.0 | - |
| General machinery | 1.5-2.0 | 2.5-3.0 | 2.0 |

### Stress Concentration Factors (K_t)

| Feature | Geometry | K_t |
|---------|----------|-----|
| Hole in tension plate | d/W = 0.25 | 2.5 |
| Hole in tension plate | d/W = 0.50 | 2.3 |
| Shoulder fillet | r/d = 0.05, D/d = 1.5 | 2.3 |
| Shoulder fillet | r/d = 0.10, D/d = 1.5 | 1.9 |
| Shoulder fillet | r/d = 0.20, D/d = 1.5 | 1.5 |
| Keyway (profiled) | - | 2.0-2.5 |
| Keyway (end-milled) | - | 3.0-3.5 |
| Thread (rolled) | - | 2.2-2.8 |
| Thread (cut) | - | 3.0-4.0 |
| Weld toe | - | 2.0-3.0 |

### Material Selection Guide

| Service Condition | Recommended Material |
|-------------------|---------------------|
| Clean water, < 80°C, < 20 bar | Ductile iron, WCB steel |
| Clean water, < 200°C, < 100 bar | Carbon steel (WCB) |
| Seawater, ambient temp | 316 SS, Duplex SS, NAB |
| Chemicals, corrosive | 316 SS, Duplex SS, Hastelloy |
| High temperature (> 400°C) | Chrome-moly steel, stainless |
| Abrasive slurries | Hard iron, duplex SS, chrome steel |
| High-speed impellers | Duplex SS, 17-4PH, titanium |
| Cryogenic (< -50°C) | 316 SS, 9% Ni steel, aluminum |

### Conversion Factors

```
Stress: 1 MPa = 1 N/mm² = 145 psi
        1 ksi = 6.895 MPa
        1 bar = 0.1 MPa = 14.5 psi

Density: 1 g/cm³ = 1,000 kg/m³

Angular velocity: ω (rad/s) = 2πn / 60
                 where n = rpm

Torque-Power: T (Nm) = 9,550 × P (kW) / n (rpm)
              T (lb-ft) = 5,252 × P (hp) / n (rpm)
```

---

## References

1. **ASME Boiler and Pressure Vessel Code, Section VIII, Divisions 1 & 2** - Pressure vessel design rules
2. **API Standard 610** - Centrifugal Pumps for Petroleum, Petrochemical and Natural Gas Industries
3. **ISO 13709** - Centrifugal pumps for petroleum, petrochemical and natural gas industries
4. **Shigley's Mechanical Engineering Design** - Comprehensive machine design reference
5. **Roark's Formulas for Stress and Strain** - Stress analysis formulas
6. **Peterson's Stress Concentration Factors** - Detailed K_t values
7. **Machinery's Handbook** - Material properties and machining data
8. **ASM Metals Handbook** - Comprehensive material properties
9. **Theory of Plates and Shells** (Timoshenko) - Classical theory
10. **ANSYS/ABAQUS Documentation** - FEA software theory manuals

---

*This reference document provides theoretical foundation for structural analysis of pump components. Always verify critical calculations with multiple methods and consult applicable design codes.*
