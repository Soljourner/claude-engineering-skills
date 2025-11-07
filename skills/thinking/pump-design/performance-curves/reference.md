# Pump Performance Curves - Reference Documentation

## Table of Contents

1. [Fundamental Equations](#fundamental-equations)
2. [Affinity Laws](#affinity-laws)
3. [Curve Shape Interpretations](#curve-shape-interpretations)
4. [Industry Standards](#industry-standards)
5. [Dimensionless Parameters](#dimensionless-parameters)
6. [Loss Analysis](#loss-analysis)
7. [References and Literature](#references-and-literature)

---

## Fundamental Equations

### Head-Flow Relationship

**General Polynomial Form:**
```
H(Q) = a₀ + a₁·Q + a₂·Q²
```

For stable operation: `a₂ < 0` (descending curve)

**Alternative Forms:**

**Shutoff-based form:**
```
H(Q) = H₀ - K₁·Q - K₂·Q²
```

**Normalized form:**
```
H/H_BEP = 1 + a·(Q/Q_BEP) + b·(Q/Q_BEP)²
```

**Theoretical basis (Euler equation):**
```
H = (u₂·c_u2 - u₁·c_u1) / g
```

For radial entry (no pre-swirl):
```
H = u₂·c_u2 / g = u₂²/g - (u₂·c_m2)/(g·tan(β₂))
```

### Efficiency Definition

**Hydraulic Efficiency:**
```
η = P_hydraulic / P_shaft = (ρ·g·Q·H) / P_shaft
```

**Component Breakdown:**
```
η = η_hydraulic · η_volumetric · η_mechanical
```

Where:
- `η_hydraulic = H / H_theoretical` (accounts for friction and shock losses)
- `η_volumetric = Q / (Q + Q_leakage)` (accounts for internal leakage)
- `η_mechanical = P_hydraulic / P_input` (accounts for bearing/seal friction)

**Typical Values:**
- η_hydraulic: 0.85-0.95
- η_volumetric: 0.96-0.99
- η_mechanical: 0.95-0.98

### Power Equation

**Shaft Power:**
```
P_shaft = (ρ·g·Q·H) / η
```

**In terms of pressure:**
```
P_shaft = (ΔP·Q) / η
```

Where `ΔP = ρ·g·H`

**Power curve forms:**

**Backward-curved blades:**
```
P ∝ Q  (nearly linear, non-overloading)
```

**Radial blades:**
```
P ∝ Q^n  where n > 1 (overloading characteristic)
```

### NPSH Required

**Definition:**
```
NPSH_req = (P_inlet - P_vapor)/(ρ·g) + v²/(2g)
```

**Empirical correlation:**
```
NPSH_req(Q) = NPSH_BEP · (Q/Q_BEP)^n
```

Where n ≈ 1.5-2.5 (typically 1.8-2.0)

**Thoma Cavitation Parameter:**
```
σ = NPSH_req / H
```

Typical values by pump type:
- Radial (Ns < 50): σ = 0.04-0.12
- Mixed flow (50 < Ns < 150): σ = 0.10-0.20
- Axial (Ns > 150): σ = 0.15-0.35

**Safety margin:**
```
NPSH_available ≥ 1.3 · NPSH_required  (minimum per ANSI/HI 9.6.1)
NPSH_available ≥ 2.0 · NPSH_required  (recommended for critical services)
```

---

## Affinity Laws

### Speed Changes (Constant Diameter, Constant Fluid)

**Flow Rate:**
```
Q₂ / Q₁ = N₂ / N₁
```

**Head:**
```
H₂ / H₁ = (N₂ / N₁)²
```

**Power:**
```
P₂ / P₁ = (N₂ / N₁)³
```

**Efficiency:**
```
η₂ ≈ η₁  (approximately constant)
```

More precisely (accounting for Reynolds number effects):
```
η₂ = 1 - (1 - η₁)·(N₁/N₂)^0.1
```

**NPSH Required:**
```
NPSH₂ / NPSH₁ = (N₂ / N₁)²
```

**Validity:**
- Accurate within ±20% speed change
- Assumes Reynolds number > 10⁶ (fully turbulent)
- Small efficiency improvement at higher speeds due to reduced relative surface roughness

### Diameter Changes (Constant Speed, Constant Fluid)

**Flow Rate:**
```
Q₂ / Q₁ = (D₂ / D₁)
```

**Head:**
```
H₂ / H₁ = (D₂ / D₁)²
```

**Power:**
```
P₂ / P₁ = (D₂ / D₁)³
```

**Efficiency with trimming penalty:**
```
η₂ = η₁ - Δη
```

Where:
```
Δη ≈ 0.02 · (1 - D₂/D₁) / 0.10  (2% per 10% trim)
```

**Practical limits:**
```
0.75 ≤ D₂/D₁ ≤ 1.00  (absolute range)
0.85 ≤ D₂/D₁ ≤ 1.00  (recommended range)
```

**Why efficiency decreases:**
1. Increased clearance gaps (relative to impeller)
2. Increased hydraulic losses
3. Deviation from optimal blade angles
4. Increased recirculation zones

### Combined Speed and Diameter Changes

**Flow Rate:**
```
Q₂ / Q₁ = (N₂/N₁) · (D₂/D₁)
```

**Head:**
```
H₂ / H₁ = (N₂/N₁)² · (D₂/D₁)²
```

**Power:**
```
P₂ / P₁ = (N₂/N₁)³ · (D₂/D₁)³
```

### Density Changes (Different Fluids or Temperatures)

**Important: Head remains constant!**
```
H₂ = H₁  (head is independent of density)
Q₂ = Q₁  (flow rate is independent of density)
```

**Power scales with density:**
```
P₂ / P₁ = ρ₂ / ρ₁ = SG₂ / SG₁
```

**Pressure change:**
```
ΔP₂ / ΔP₁ = ρ₂ / ρ₁
```

**Viscosity effects:**

For significant viscosity changes (Re < 10⁶), use HI correction factors:

```
Q_viscous = Q_water · C_Q
H_viscous = H_water · C_H
η_viscous = η_water · C_η
```

Where C_Q, C_H, C_η are correction factors from HI charts (function of Q, Ns, viscosity)

---

## Curve Shape Interpretations

### H-Q Curve Shape vs Specific Speed

**Low Specific Speed (Ns = 10-30, Radial Flow):**
- Steep, nearly linear H-Q curve
- High shutoff head ratio: H₀/H_BEP = 1.25-1.35
- Stable throughout operating range
- dH/dQ strongly negative

**Medium Specific Speed (Ns = 30-80, Francis Type):**
- Moderate slope H-Q curve
- Moderate shutoff ratio: H₀/H_BEP = 1.15-1.25
- Generally stable
- May show slight dip at low flows

**High Specific Speed (Ns = 80-150, Mixed Flow):**
- Flat H-Q curve
- Low shutoff ratio: H₀/H_BEP = 1.10-1.15
- Risk of instability at low flows
- Potential dip in curve due to recirculation

**Very High Specific Speed (Ns > 150, Axial Flow):**
- Very flat H-Q curve
- Minimal shutoff head increase
- Often shows S-shaped curve with unstable region
- Requires careful system matching

### Stability Analysis

**Stable Curve:**
```
dH/dQ < 0  for all Q in operating range
```

Characteristics:
- Single operating point with system curve
- Predictable performance
- Safe for parallel operation
- No hunting or oscillation

**Unstable Curve (Dip):**
```
dH/dQ > 0  for some Q range
```

Problems:
- Multiple possible operating points
- Hunting between stable points
- Parallel pump flow distribution issues
- Vibration and noise

**Causes of instability:**
1. Recirculation at impeller inlet (low flow)
2. Flow separation in diffuser
3. Rotating stall
4. Poor impeller design (forward-curved blades)

**Solutions:**
1. Minimum flow bypass
2. Impeller redesign
3. Operating range restrictions
4. Variable speed control

### Efficiency Curve Characteristics

**Peak location vs specific speed:**

- **Low Ns (radial):** Broad, flat efficiency peak
  - η_max at Q/Q_shutoff ≈ 0.6-0.7
  - Wide operating range (50-140% Q_BEP)

- **Medium Ns (mixed):** Sharp efficiency peak
  - η_max at Q/Q_shutoff ≈ 0.7-0.8
  - Moderate range (60-130% Q_BEP)

- **High Ns (axial):** Very sharp peak
  - η_max at Q/Q_shutoff ≈ 0.8-0.9
  - Narrow range (70-120% Q_BEP)

**Efficiency vs pump size (Scale effects):**

Gülich correlation:
```
η_opt ≈ 0.94 - 0.0525 / √Ns  (European Ns units)
```

Anderson correlation:
```
η_max ≈ 1 - 0.8 / (Ns/50)^0.25
```

Typical values:
- Small pumps (Q < 50 m³/h): η = 50-75%
- Medium pumps (50-500 m³/h): η = 70-85%
- Large pumps (Q > 500 m³/h): η = 80-92%

### Power Curve Implications

**Non-overloading characteristic (backward-curved blades):**
```
P_max occurs at Q_max (runout)
```

**Advantage:** Motor sized for maximum flow won't be overloaded at any operating point

**Overloading characteristic (radial/forward-curved blades):**
```
P increases sharply with Q
```

**Risk:** Motor can be overloaded at high flow rates

**Motor selection:**
```
P_motor_rated ≥ P_max / η_motor · SF
```

Where:
- P_max = maximum power in operating range
- η_motor = motor efficiency (typically 0.90-0.96)
- SF = service factor (1.10-1.25)

---

## Industry Standards

### ISO 9906 (Performance Testing)

**Scope:** Hydraulic performance acceptance tests for centrifugal, mixed flow, and axial pumps

**Tolerance Grades:**

**Grade 1 (Precision testing):**
- Flow: ±2.5%
- Head: ±3.0%
- Power: ±4.0%
- Efficiency: ±2.0 percentage points

**Grade 2 (Standard testing):**
- Flow: ±4.0%
- Head: ±5.0%
- Power: ±6.5%
- Efficiency: ±3.5 percentage points

**Grade 3 (Routine testing):**
- Flow: ±7.0%
- Head: ±7.5%
- Power: ±9.0%
- Efficiency: ±5.0 percentage points

**Test requirements:**
- Minimum 3 test points including rated point
- Test points between 70-120% of rated flow
- Cavitation tests at 100%, 110%, 120% rated flow

### ANSI/HI 9.6.1 (Pump Acceptance Testing)

**Standard conditions:**
- Atmospheric pressure: 101.3 kPa
- Water temperature: 20°C
- Density: 998 kg/m³

**Acceptance criteria:**
- Must meet or exceed manufacturer's published curves
- Tolerances per contract specification
- NPSH testing required for critical applications

### HI 9.6.3 (NPSH Testing)

**3% head drop criterion:**
```
NPSH_3% = NPSH at which H = 0.97 · H_non-cavitating
```

**0% head drop criterion (NPSH₀):**
```
NPSH_0 = NPSH at incipient cavitation
```

Relationship:
```
NPSH_3% ≈ 1.3 · NPSH_0  (typical)
```

**Testing procedure:**
1. Establish baseline head at high NPSH
2. Progressively reduce suction pressure
3. Record head vs NPSH_available
4. Determine NPSH_3% from curve

### API 610 (Centrifugal Pumps for Petroleum Industry)

**Performance requirements:**
- Preferred operating range: 70-120% BEP
- Allowable operating range: 60-130% BEP (with restrictions)
- NPSH margin: NPSH_avail ≥ NPSH_req + 0.6 m (minimum)

**Curve stability:**
- H-Q curve must be continuously rising (dH/dQ < 0)
- No dips or humps allowed
- Critical for parallel operation

**Testing:**
- Witnessed performance test per API 610 Section 6.2
- NPSH test per API 610 Annex F
- Mechanical run test

---

## Dimensionless Parameters

### Specific Speed

**European (dimensionless) convention:**
```
Ns = N · √Q / H^(3/4)
```

Units:
- N: rpm
- Q: m³/s
- H: m

**US customary units:**
```
Ns_US = N(rpm) · √Q(gpm) / H(ft)^(3/4)
```

Conversion:
```
Ns_US = 51.64 · Ns_European
```

**Physical meaning:**
- Characterizes impeller geometry
- Indicates best pump type for application
- Correlates with efficiency

### Suction Specific Speed

```
S = N · √Q / NPSH_req^(3/4)
```

**Typical limits:**
- Radial pumps: S < 230 (US units), S < 4.5 (SI)
- Mixed flow: S < 280 (US units), S < 5.4 (SI)
- Axial: S < 350 (US units), S < 6.8 (SI)

**High S values:**
- Indicate good suction performance
- Risk of cavitation if NPSH not met

### Flow Coefficient

```
φ = Q / (N·D³) = c_m / u
```

Or:
```
φ = Q / (ω·D³) = Q / (π·N·D³/60)
```

### Head Coefficient

```
ψ = g·H / (N·D)² = g·H / u²
```

Or:
```
ψ = g·H / (ω·D)²
```

### Power Coefficient

```
λ = P / (ρ·N³·D⁵)
```

### Reynolds Number

```
Re = ρ·u·D / μ = ρ·ω·D² / μ
```

**Significance:**
- Re > 10⁶: Fully turbulent, scale-independent performance
- Re < 10⁶: Viscous effects, reduced performance
- Re < 10⁵: Significant efficiency penalty

---

## Loss Analysis

### Major Loss Components

**1. Friction Losses (h_f):**
```
h_f = λ·(L/D)·(v²/2g)
```

Contribution: 30-50% of total losses

**2. Shock Losses (h_shock):**

At design point: minimal
Off-design:
```
h_shock = k·|i|²/(2g)
```

Where i = incidence angle

**3. Leakage Losses:**
```
Q_leak = C_d·A_gap·√(2g·ΔH)
```

Reduces volumetric efficiency:
```
η_vol = Q/(Q + Q_leak)
```

**4. Disk Friction:**
```
P_disk = k_f·ρ·ω³·r⁵
```

Significant in large, high-speed pumps

**5. Recirculation Losses:**

Occur at Q < 0.5·Q_BEP
- Inlet recirculation: low NPSH margin, vibration
- Outlet recirculation: efficiency drop, heating

### Total Loss Coefficient

```
H_theoretical - H_actual = Σh_losses
```

Efficiency:
```
η_h = H_actual / H_theoretical = 1 / (1 + Σk_loss)
```

---

## System Curve and Operating Point

### System Head Equation

**General form:**
```
H_system(Q) = H_static + H_friction + H_minor
```

**Expanded:**
```
H_sys = (z₂ - z₁) + (P₂ - P₁)/(ρ·g) + (8·f·L·Q²)/(π²·g·D⁵) + Σ(K_i·v²/(2g))
```

**Simplified form:**
```
H_sys = H_s + K·Q²
```

Where:
```
K = (8·f·L)/(π²·g·D⁵) + (8·Σ K_minor)/(π²·g·D⁴)
```

### Operating Point Determination

**Graphical method:**
1. Plot pump H-Q curve
2. Plot system H-Q curve
3. Intersection = operating point

**Analytical method:**

Solve:
```
H_pump(Q) = H_system(Q)
```

For quadratic forms:
```
a₀ + a₁·Q + a₂·Q² = H_s + K·Q²
```

Solution:
```
Q_op = [-a₁ ± √(a₁² - 4·(a₂-K)·(a₀-H_s))] / [2·(a₂-K)]
```

Take positive root.

### Multiple Pumps

**Series operation:**
```
H_combined = H₁ + H₂
Q_combined = Q  (same through both)
```

**Parallel operation:**
```
Q_combined = Q₁ + Q₂
H_combined = H  (same for both)
```

**Parallel operating point:**

Each pump must satisfy:
```
H₁(Q₁) = H₂(Q₂) = H_system(Q₁ + Q₂)
```

**Warning:** Unequal flow distribution if pump curves differ!

---

## Curve Fitting Techniques

### Polynomial Regression

**H-Q curve (2nd order):**
```python
import numpy as np
coeffs = np.polyfit(Q_data, H_data, 2)
```

**Quality metrics:**

**Coefficient of determination:**
```
R² = 1 - SS_residual / SS_total
```

Target: R² > 0.99

**Root Mean Square Error:**
```
RMSE = √(Σ(y_i - ŷ_i)² / n)
```

### Weighted Least Squares

Emphasize important regions (e.g., near BEP):

```python
weights = np.exp(-k·(Q - Q_BEP)²)
coeffs = np.polyfit(Q_data, H_data, 2, w=weights)
```

### Physical Constraints

For physically realistic curves:

1. **Head curve:**
   - H(0) > 0 (shutoff head must be positive)
   - dH/dQ < 0 (stable curve)

2. **Efficiency curve:**
   - η(0) = 0 (zero efficiency at zero flow)
   - 0 ≤ η ≤ 1 (physical bounds)
   - Single maximum (BEP)

3. **NPSH curve:**
   - NPSH(Q) > 0 (always positive)
   - dNPSH/dQ > 0 (increases with flow)

---

## References and Literature

### Textbooks

1. **Gülich, J.F.** (2020). *Centrifugal Pumps*, 4th Edition. Springer.
   - Comprehensive reference for pump design
   - Detailed performance prediction methods
   - Chapter 6: Performance characteristics

2. **Karassik, I.J., et al.** (2008). *Pump Handbook*, 4th Edition. McGraw-Hill.
   - Industry standard reference
   - Section 2.1: Centrifugal pump theory
   - Section 2.2: Performance characteristics

3. **Stepanoff, A.J.** (1957). *Centrifugal and Axial Flow Pumps*, 2nd Edition. Wiley.
   - Classic text on pump theory
   - Fundamental performance relationships

4. **Lobanoff, V.S. & Ross, R.R.** (1992). *Centrifugal Pumps: Design and Application*, 2nd Edition. Gulf Professional Publishing.
   - Practical design methods
   - Chapter 4: Performance characteristics

### Standards

1. **ISO 9906:2012** - Rotodynamic pumps - Hydraulic performance acceptance tests - Grades 1, 2 and 3

2. **ANSI/HI 9.6.1-2012** - Rotodynamic Pumps Guideline for NPSH Margin

3. **ANSI/HI 9.6.3-2012** - Rotodynamic Pumps Guideline for Net Positive Suction Head Margin

4. **API 610 (12th Edition)** - Centrifugal Pumps for Petroleum, Petrochemical and Natural Gas Industries

5. **ANSI/HI 9.6.7-2015** - Rotodynamic Pumps Guideline for Effects of Liquid Viscosity on Performance

### Technical Papers

1. **Anderson, H.H.** (1980). "A Design Method for Centrifugal Compressor Impellers." ASME Paper 80-GT-70.
   - Efficiency correlations

2. **Lazarkiewicz, S. & Troskolanski, A.T.** (1965). "Impeller Pumps." Pergamon Press.
   - Velocity triangle analysis

3. **Brennen, C.E.** (1994). "Hydrodynamics of Pumps." Oxford University Press.
   - Cavitation and NPSH theory

4. **Japikse, D., et al.** (1997). "Introduction to Turbomachinery." Concepts ETI.
   - Performance prediction methods

### Online Resources

1. **Hydraulic Institute (HI):** www.pumps.org
   - Standards and guidelines
   - Educational resources

2. **ANSI/HI Pump Standards:**
   - Complete set of pump industry standards

3. **Engineering Toolbox:**
   - Fluid properties
   - Pump calculations

---

## Nomenclature

| Symbol | Description | Units |
|--------|-------------|-------|
| Q | Volumetric flow rate | m³/s, m³/h |
| H | Total head | m |
| P | Power | W, kW |
| η | Efficiency | - (fraction) |
| N | Rotational speed | rpm |
| D | Impeller diameter | m |
| ρ | Density | kg/m³ |
| g | Gravitational acceleration | m/s² |
| Ns | Specific speed | dimensionless |
| NPSH | Net positive suction head | m |
| σ | Thoma cavitation parameter | - |
| u | Blade velocity | m/s |
| c | Absolute fluid velocity | m/s |
| w | Relative fluid velocity | m/s |
| β | Blade angle | degrees, radians |
| f | Friction factor | - |
| Re | Reynolds number | - |
| BEP | Best efficiency point | - |

### Subscripts

- 0: shutoff condition (Q = 0)
- 1: impeller inlet
- 2: impeller outlet
- s: static
- th: theoretical
- req: required
- avail: available
- op: operating point

---

*Document Version: 1.0*
*Last Updated: 2025*
*Cross-referenced with ISO 9906, API 610, and ANSI/HI standards*
