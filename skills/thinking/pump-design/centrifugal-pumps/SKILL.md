---
name: centrifugal-pump-design
description: "Design centrifugal pumps using Euler equations, velocity triangles, and specific speed"
category: thinking
domain: mechanical
complexity: advanced
dependencies:
  - numpy
  - matplotlib
---

# Centrifugal Pump Design

## Overview

This skill provides a comprehensive workflow for designing centrifugal pumps from first principles using Euler turbine equations, velocity triangles, and dimensionless parameters. The methodology follows classical pump design theory as presented in Stepanoff, Gülich, and Karassik.

## Design Workflow

### 1. Requirements Analysis

**Input Parameters:**
- Flow rate: Q (m³/s or m³/h)
- Head: H (m)
- Rotational speed: N (rpm)
- Fluid properties:
  - Density: ρ (kg/m³)
  - Kinematic viscosity: ν (m²/s)
  - NPSH available (m)

**Derived Parameters:**
- Power requirement: P = ρ·g·Q·H / η
- Specific speed: Ns = N·√Q / H^0.75 (dimensionless form)
- Suction specific speed: S = N·√Q / NPSH^0.75

### 2. Specific Speed Calculation and Pump Type Selection

**Specific Speed Definition:**

The specific speed is a dimensionless parameter that characterizes pump geometry:

```
Ns = (N·Q^0.5) / H^0.75
```

Where:
- N = rotational speed (rpm)
- Q = flow rate (m³/s)
- H = head per stage (m)

**Pump Type Selection Based on Ns:**

| Specific Speed (Ns) | Pump Type | Impeller Shape | Typical Efficiency |
|---------------------|-----------|----------------|-------------------|
| 10-30 | Radial flow | Narrow, radial | 70-85% |
| 30-50 | Francis-vane | Medium width | 80-88% |
| 50-80 | Mixed flow | Wide, angled | 85-90% |
| 80-150 | Mixed flow | Very wide | 85-92% |
| 150-300 | Axial flow | Propeller | 82-88% |

**Note:** Specific speed units vary by convention. The European convention uses:
```
nq = N·√Q / H^0.75  [where Q in m³/s, H in m, N in rpm]
```

US convention multiplies by different constants. Always verify which convention is being used.

### 3. Impeller Diameter Estimation

**Stepanoff Correlation for Outlet Diameter:**

```
D2 = 84.6 · (H / N)^0.5  [D2 in mm, H in m, N in rpm]
```

Or in SI units:
```
D2 = K · (H·g / u2²)^0.5 · u2 / (0.5·π·N/60)
```

**Simplified form:**
```
D2 = 60 / π · √(2·g·H / Ku²) / N  [D2 in m]
```

Where Ku is the outlet velocity coefficient (typically 0.95-1.05)

**Outlet Width Estimation:**

```
b2 = K_b · D2 / Ns^0.65
```

Where K_b = 2.5-3.5 for radial pumps

**Eye Diameter Estimation:**

Based on suction specific speed:
```
D_eye = (4·Q / (π·c_m1))^0.5
```

Where c_m1 is the meridional velocity at inlet (typically 2-4 m/s)

### 4. Velocity Triangle Construction

**Velocity Components:**

At any point in the impeller:
- **u** = blade velocity (tangential)
- **c** = absolute fluid velocity
- **w** = relative velocity (fluid velocity relative to blade)

**Velocity Triangle Relations:**

```
c = u + w  (vector addition)
```

**Component breakdown:**
- c_m = meridional (through-flow) component
- c_u = tangential (whirl) component
- u = ω·r = π·N·D/60

**Inlet Triangle (Station 1 - Eye):**

For axial inlet with no pre-rotation:
```
c_u1 = 0  (design condition)
c_m1 = Q / A1 = Q / (π·D_eye²/4)
c1 = c_m1
u1 = π·N·D_eye/60
w1 = √(u1² + c_m1²)
β1 = atan(c_m1 / u1)  [blade angle at inlet]
```

**Outlet Triangle (Station 2):**

```
u2 = π·N·D2/60
c_m2 = Q / A2 = Q / (π·D2·b2)
```

From Euler equation and design choices:
```
c_u2 = g·H / u2  (for backward-curved blades)
c2 = √(c_m2² + c_u2²)
w2 = √((u2 - c_u2)² + c_m2²)
β2 = atan(c_m2 / (u2 - c_u2))  [blade angle at outlet]
```

**Important:** β2 is typically 15-35° for backward-curved blades

### 5. Euler Turbine Equation Application

**Fundamental Equation:**

The theoretical head developed by a pump is:

```
H_th = (u2·c_u2 - u1·c_u1) / g
```

**For radial entry (c_u1 = 0):**

```
H_th = u2·c_u2 / g
```

**In terms of blade angles:**

```
H_th = u2² / g - u2·c_m2 / (g·tan(β2))
```

**Slip Factor Correction:**

Real pumps have finite number of blades, causing slip:

```
σ = 1 - (π·sin(β2)) / Z
```

Where Z is the number of blades (typically 5-9 for radial pumps)

**Actual head:**
```
H = σ·H_th
```

### 6. Blade Angle Calculations

**Inlet Blade Angle (β1):**

```
β1 = atan(c_m1 / u1)
```

Typically: β1 = 15-30° for radial pumps

**Outlet Blade Angle (β2):**

```
β2 = atan(c_m2 / (u2 - c_u2))
```

**Blade Type Selection:**
- **Backward-curved:** β2 < 90° (most common, stable operation)
- **Radial:** β2 = 90° (higher head, less stable)
- **Forward-curved:** β2 > 90° (rarely used, unstable)

**Blade Number Selection:**

Pfleiderer formula:
```
Z = 6.5 · (D2 + D1)/(D2 - D1) · sin((β1 + β2)/2)
```

Round to nearest integer (typically 5-9 blades)

### 7. Efficiency Estimation

**Component Efficiencies:**

Total efficiency:
```
η = η_h · η_vol · η_mech
```

**Hydraulic Efficiency (η_h):**

Accounts for friction and shock losses:
```
η_h = H / H_th ≈ 0.85-0.95
```

More detailed (Gülich):
```
η_h = 1 / (1 + k_friction + k_shock + k_recirculation)
```

**Volumetric Efficiency (η_vol):**

Accounts for leakage:
```
η_vol = Q / (Q + Q_leak) ≈ 0.96-0.99
```

**Mechanical Efficiency (η_mech):**

Accounts for bearing and seal friction:
```
η_mech = (P_hydraulic) / (P_shaft) ≈ 0.95-0.98
```

**Overall Efficiency Estimation:**

For preliminary design (Gülich correlation):
```
η_opt ≈ 0.94 - 0.0525·(Ns)^(-0.5)  [for Ns in European units]
```

Or Anderson correlation:
```
η = 1 - 0.8 / (Ns/50)^0.25  [for well-designed pumps]
```

### 8. Performance Curve Generation

**Head-Flow Characteristic:**

The H-Q curve can be approximated as:
```
H = a·Q² + b·Q + c
```

Where coefficients are determined from:
- Design point (Q, H)
- Shut-off head (Q=0, H_max)
- Slope at design point

**Power Curve:**
```
P = ρ·g·Q·H / η(Q)
```

**NPSH Required:**

```
NPSH_req = (c1²)/(2g) + k·(c1²)/(2g)
```

Where k accounts for acceleration and losses (k ≈ 1.5-2.5)

### 9. Loss Analysis

**Major Loss Sources:**

1. **Friction losses:** Proportional to c²/2g
2. **Shock losses:** At off-design conditions
3. **Leakage losses:** Through clearances
4. **Disk friction:** Proportional to ω³·r⁵
5. **Recirculation:** At low flow rates

**Total head loss coefficient:**
```
k_total = k_friction + k_shock + k_incidence + k_diffusion
```

## Design Rules of Thumb

### General Guidelines

1. **Specific Speed Range:**
   - Single-stage: Ns = 10-150
   - Multi-stage for Ns < 20 (high head, low flow)
   - Axial for Ns > 150 (low head, high flow)

2. **Blade Angles:**
   - Inlet: β1 = 18-28° (typical)
   - Outlet: β2 = 20-30° (backward-curved)
   - Never exceed β2 > 40° (poor efficiency)

3. **Velocities:**
   - Meridional: c_m = 2-5 m/s
   - Outlet absolute: c2 = 10-25 m/s
   - Inlet relative: w1 < 20 m/s (cavitation limit)

4. **Dimensional Ratios:**
   - D2/D1 = 1.5-2.5 (radial pumps)
   - b2/D2 = 0.02-0.08 (decreases with Ns)
   - Blade thickness: t = 3-6 mm (depends on material and size)

5. **Number of Blades:**
   - Radial pumps: Z = 5-7
   - Mixed flow: Z = 5-8
   - Axial flow: Z = 3-5
   - More blades = higher head, lower efficiency

6. **Clearances:**
   - Impeller shroud clearance: 0.2-0.5 mm (small pumps)
   - Wear ring clearance: 0.3-0.8 mm
   - Minimize to reduce leakage losses

### Performance Guidelines

1. **Efficiency:**
   - Small pumps (P < 20 kW): η = 60-75%
   - Medium pumps (20-200 kW): η = 75-85%
   - Large pumps (P > 200 kW): η = 85-92%

2. **Operating Range:**
   - Recommended: 70-120% of design flow
   - Maximum: 50-150% of design flow
   - Avoid continuous operation below 40% (recirculation)

3. **NPSH Safety Margin:**
   - NPSH_available > 1.3 × NPSH_required (minimum)
   - NPSH_available > 2.0 × NPSH_required (recommended)

### Material Selection

1. **Impeller Materials:**
   - Cast iron: Standard water service
   - Bronze: Seawater, corrosive fluids
   - Stainless steel (316): Chemical service
   - Duplex stainless: High corrosion + high strength

2. **Surface Finish:**
   - Ra < 1.6 μm for hydraulic surfaces
   - Ra < 3.2 μm for back shroud
   - Polish for slurries (Ra < 0.8 μm)

## Common Design Mistakes to Avoid

### Critical Errors

1. **Incorrect Specific Speed Application:**
   - ❌ Using total head instead of head per stage
   - ❌ Mixing unit systems (US vs. European Ns)
   - ✅ Always calculate Ns with consistent units
   - ✅ Select pump type based on correct Ns range

2. **Velocity Triangle Errors:**
   - ❌ Assuming c_u1 ≠ 0 without pre-rotation device
   - ❌ Using u1 at hub instead of mean diameter
   - ✅ Calculate velocities at correct radii
   - ✅ Verify vector addition: c = u + w

3. **Excessive Blade Angles:**
   - ❌ β2 > 35° causes instability and low efficiency
   - ❌ β1 < 10° causes manufacturing difficulties
   - ✅ Keep β2 = 20-30° for stable operation
   - ✅ Ensure blade angles are manufacturable

4. **Ignoring Slip:**
   - ❌ Using H_th without slip correction
   - ❌ Neglecting finite blade number effects
   - ✅ Apply slip factor: σ = 1 - π·sin(β2)/Z
   - ✅ Actual head H = σ·H_th

5. **Poor Suction Design:**
   - ❌ Insufficient NPSH margin
   - ❌ Excessive inlet velocity (cavitation)
   - ✅ Keep inlet velocity < 3-4 m/s
   - ✅ NPSH_avail > 1.5 × NPSH_req minimum

6. **Unrealistic Efficiency:**
   - ❌ Assuming η > 0.90 for small pumps
   - ❌ Not accounting for scale effects
   - ✅ Use correlations based on size and Ns
   - ✅ Validate against similar pump data

### Design Optimization Issues

1. **Over-Designing:**
   - Using too many blades (increases friction)
   - Excessive impeller diameter (higher cost, power)
   - Too narrow outlet width (high velocity, losses)

2. **Under-Designing:**
   - Too few blades (slip losses, uneven flow)
   - Insufficient outlet width (cavitation)
   - Inadequate structural analysis (fatigue failure)

3. **Manufacturing Constraints:**
   - Blade angles too complex for casting
   - Insufficient blade thickness for welding
   - Tight tolerances without justification

4. **System Mismatch:**
   - Selecting pump for single duty point
   - Not considering system curve variation
   - Ignoring start-up and shutdown conditions

## Verification and Validation

### Design Verification Checklist

- [ ] Specific speed within acceptable range for pump type
- [ ] Blade angles manufacturable (β1, β2 > 10°)
- [ ] Velocities within recommended limits
- [ ] NPSH available exceeds required by 30%+
- [ ] Efficiency estimate reasonable for pump size
- [ ] Slip factor applied to theoretical head
- [ ] Structural analysis for rotational stresses
- [ ] Leakage paths minimized with proper clearances

### Performance Validation

1. **Dimensional Analysis:**
   - Check Reynolds number: Re = ρ·u2·D2/μ > 10⁶
   - Verify Mach number: Ma = u2/a < 0.3

2. **Affinity Laws:**
   - Q₂/Q₁ = (N₂/N₁) · (D₂/D₁)³
   - H₂/H₁ = (N₂/N₁)² · (D₂/D₁)²
   - P₂/P₁ = (N₂/N₁)³ · (D₂/D₁)⁵

3. **Comparison with Similar Pumps:**
   - Benchmark against published data
   - Compare dimensional ratios
   - Validate efficiency predictions

## Advanced Topics

### Multi-Stage Pumps

For high heads:
```
Number of stages = H_total / H_per_stage
```

Select H_per_stage for optimal Ns (40-80 for best efficiency)

### Variable Speed Operation

Using affinity laws:
```
Q ∝ N
H ∝ N²
P ∝ N³
```

### Cavitation Analysis

Thoma cavitation parameter:
```
σ = NPSH_req / H
```

Typical values:
- σ = 0.06-0.12 (radial pumps)
- σ = 0.15-0.30 (axial pumps)

### Computational Fluid Dynamics (CFD)

Modern design uses CFD for:
- Detailed flow field analysis
- Loss prediction
- Cavitation prediction
- Performance optimization

However, preliminary design still follows classical theory presented here.

## Summary

Centrifugal pump design is a systematic process combining:
1. Dimensional analysis (specific speed)
2. Euler turbine equation (energy transfer)
3. Velocity triangles (kinematics)
4. Empirical correlations (efficiency, losses)

Success requires balancing theoretical calculations with practical experience and validation against existing designs.

## References

See `reference.md` for detailed derivations, correlations, and literature references.
