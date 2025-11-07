---
name: pump-performance-curves
description: "Generate and interpret H-Q curves, apply affinity laws, and predict off-design performance"
category: thinking
domain: mechanical
complexity: intermediate
dependencies:
  - numpy
  - matplotlib
  - scipy
---

# Pump Performance Curves

## Overview

This skill provides comprehensive methods for generating, interpreting, and manipulating pump performance curves. Performance curves are fundamental to pump selection, system design, and operational analysis. The skill covers the four primary curve types (H-Q, η-Q, P-Q, and NPSH-Q), affinity laws for scaling performance, curve fitting techniques, and operating point determination.

## Performance Curve Components

### 1. Head-Flow (H-Q) Curve

The H-Q curve shows the relationship between total dynamic head and volumetric flow rate.

**Characteristics:**
- **Shutoff Head (H₀):** Head at zero flow (Q = 0)
- **Design Point:** Nominal operating condition (Q_design, H_design)
- **Runout Flow:** Maximum flow at zero head
- **Curve Shape:** Typically descending, influenced by impeller design

**Mathematical Representation:**

```
H(Q) = H₀ - a·Q - b·Q²
```

Or polynomial form:
```
H(Q) = c₀ + c₁·Q + c₂·Q²
```

**Physical Interpretation:**
- Steep curve: Radial impeller, low specific speed
- Flat curve: Axial impeller, high specific speed
- Rising curve (unstable): Poor design, avoid

**Stability Criteria:**
```
dH/dQ < 0  (stable operation)
dH/dQ > 0  (unstable - avoid)
```

### 2. Efficiency-Flow (η-Q) Curve

The efficiency curve shows how effectively the pump converts shaft power to hydraulic power.

**Key Points:**
- **Best Efficiency Point (BEP):** Peak of efficiency curve
- **Preferred Operating Range:** 70-120% of BEP flow
- **Acceptable Range:** 50-150% of BEP flow
- **Avoid:** Operation below 40% or above 150% of BEP

**Efficiency Definition:**
```
η = (ρ·g·Q·H) / P_shaft
```

**Typical Efficiency Curve Shape:**
```
η(Q) = η_max · (1 - k₁·(Q/Q_BEP - 1)² - k₂·(Q/Q_BEP - 1)⁴)
```

**Efficiency by Pump Size:**
- Small pumps (P < 20 kW): η = 50-75%
- Medium pumps (20-200 kW): η = 70-85%
- Large pumps (P > 200 kW): η = 80-92%

### 3. Power-Flow (P-Q) Curve

The power curve shows shaft power consumption vs. flow rate.

**Power Calculation:**
```
P_shaft = (ρ·g·Q·H) / η
```

**Curve Shapes by Impeller Type:**

**Backward-Curved Blades (most common):**
- Power increases with flow
- Non-overloading characteristic
- Motor sized for maximum flow

**Radial Blades:**
- Steep power increase
- Can overload motor
- Requires careful motor selection

**Forward-Curved Blades (rare):**
- Very steep power curve
- Highly overloading
- Generally avoided

**Motor Selection:**
```
P_motor = P_max / η_motor · SF
```
Where SF = safety factor (typically 1.10-1.25)

### 4. NPSH-Flow (NPSH_req-Q) Curve

Net Positive Suction Head Required curve shows minimum inlet pressure needed to avoid cavitation.

**NPSH Definition:**
```
NPSH_req = (P_inlet - P_vapor)/(ρ·g) + v²/(2g)
```

**Typical Behavior:**
```
NPSH_req(Q) = NPSH₀ + k·(Q/Q_BEP)^n
```
Where n ≈ 1.5-2.5

**Safety Margin:**
```
NPSH_available > NPSH_required · 1.3  (minimum)
NPSH_available > NPSH_required · 2.0  (recommended)
```

**Thoma Cavitation Parameter:**
```
σ = NPSH_req / H
```

Typical values:
- Radial pumps: σ = 0.06-0.12
- Mixed flow: σ = 0.10-0.20
- Axial pumps: σ = 0.15-0.35

## Curve Characteristics and Interpretation

### Shutoff Head

**Definition:** Maximum head at zero flow.

**Typical Ratios:**
```
H₀/H_BEP = 1.10-1.35  (centrifugal pumps)
```

**Significance:**
- Protection against deadhead operation
- Check valve selection
- System pressure relief sizing

### Best Efficiency Point (BEP)

**Definition:** Operating point with maximum hydraulic efficiency.

**Significance:**
- Minimum wear and vibration
- Longest bearing life
- Lowest energy consumption
- Reference for operating range

**Design Guidelines:**
- Size pump so normal duty is 70-110% of BEP
- Avoid continuous operation far from BEP
- Maximum deviation: ±30% of Q_BEP

**BEP Identification:**
```
dη/dQ = 0  (at BEP)
d²η/dQ² < 0  (maximum, not minimum)
```

### Operating Range

**Preferred Operating Range (POR):**
```
0.7·Q_BEP ≤ Q ≤ 1.2·Q_BEP
```

**Acceptable Operating Range (AOR):**
```
0.5·Q_BEP ≤ Q ≤ 1.5·Q_BEP
```

**Continuous Operation Limits:**
- Minimum: 40% Q_BEP (avoid recirculation)
- Maximum: 150% Q_BEP (avoid cavitation and vibration)

### Stable vs. Unstable Regions

**Stable Curve:**
```
dH/dQ < 0  for all Q
```
- Single intersection with system curve
- Predictable operation
- Desirable characteristic

**Unstable Curve:**
```
dH/dQ > 0  for some Q range
```
- Multiple intersections possible
- Hunting and oscillation
- Parallel pump problems
- Avoid in design

**Dip in H-Q Curve:**
- Common at low specific speeds
- Caused by recirculation
- Creates unstable region
- May require minimum flow bypass

## Affinity Laws

The affinity laws allow scaling pump performance for speed or diameter changes.

### Speed Changes (Constant Diameter)

**Flow Rate:**
```
Q₂/Q₁ = N₂/N₁
```

**Head:**
```
H₂/H₁ = (N₂/N₁)²
```

**Power:**
```
P₂/P₁ = (N₂/N₁)³
```

**Efficiency:**
```
η₂ ≈ η₁  (approximately constant)
```

**NPSH Required:**
```
NPSH₂/NPSH₁ = (N₂/N₁)²
```

**Validity:**
- Accurate for ±20% speed change
- Assumes Reynolds number independence
- Efficiency slightly improves at higher speeds

**Example Application:**
```python
# Original: 1750 RPM, Q=100 m³/h, H=50m, P=20 kW
# New speed: 1450 RPM

N_ratio = 1450/1750  # = 0.829
Q_new = 100 * 0.829 = 82.9 m³/h
H_new = 50 * 0.829² = 34.4 m
P_new = 20 * 0.829³ = 11.4 kW
```

### Impeller Diameter Changes (Constant Speed)

**Flow Rate:**
```
Q₂/Q₁ = (D₂/D₁)
```

**Head:**
```
H₂/H₁ = (D₂/D₁)²
```

**Power:**
```
P₂/P₁ = (D₂/D₁)³
```

**Efficiency:**
```
η₂ ≈ η₁ - Δη
```
Where Δη accounts for increased clearance effects (typically 2-5% loss)

**Practical Limits:**
- Maximum trim: 75% of original diameter
- Below 75%: Efficiency penalty increases significantly
- Minimum recommended: 85% of original diameter

**Trimming Guidelines:**
```
D_min = 0.75·D_original  (absolute minimum)
D_practical = 0.85·D_original  (recommended minimum)
```

**Example Application:**
```python
# Original: D=300mm, Q=200 m³/h, H=80m, P=60 kW
# Trimmed to: D=270mm (90% of original)

D_ratio = 270/300  # = 0.90
Q_new = 200 * 0.90 = 180 m³/h
H_new = 80 * 0.90² = 64.8 m
P_new = 60 * 0.90³ = 43.7 kW
η_penalty = 0.02  # 2% efficiency loss
```

### Combined Changes

For simultaneous speed AND diameter changes:

```
Q₂/Q₁ = (N₂/N₁) · (D₂/D₁)
H₂/H₁ = (N₂/N₁)² · (D₂/D₁)²
P₂/P₁ = (ρ₂/ρ₁) · (N₂/N₁)³ · (D₂/D₁)³
```

### Density Effects

For different fluids or temperatures:

```
H₂ = H₁  (head is independent of density)
P₂/P₁ = ρ₂/ρ₁  (power scales with density)
NPSH₂ = NPSH₁  (head-based, not pressure)
```

**Important:** Specific gravity affects power, not head!

## Curve Fitting from Test Data

### Polynomial Curve Fitting

**H-Q Curve (2nd Order):**
```
H = a₀ + a₁·Q + a₂·Q²
```

**Fitting Methods:**
1. Least squares regression
2. Three-point method (Q=0, Q=Q_BEP, Q=Q_max)
3. Weighted least squares (emphasize BEP region)

**Implementation:**
```python
import numpy as np
from scipy.optimize import curve_fit

# Polynomial model
def head_model(Q, a0, a1, a2):
    return a0 + a1*Q + a2*Q**2

# Fit to data
params, _ = curve_fit(head_model, Q_data, H_data)
```

**Quality Metrics:**
```
R² = 1 - SS_res/SS_tot  (coefficient of determination)
RMSE = √(Σ(H_predicted - H_measured)²/n)
```

Target: R² > 0.99 for good fit

### Efficiency Curve Fitting

**Gaussian-like Model:**
```
η(Q) = η_max · exp(-k·(Q - Q_BEP)²)
```

**Polynomial Model:**
```
η(Q) = b₀ + b₁·Q + b₂·Q² + b₃·Q³
```

**Constraints:**
- η(0) = 0
- η(Q_BEP) = η_max
- η → 0 as Q → ∞

### Power Curve Fitting

**Direct Calculation:**
```
P(Q) = ρ·g·Q·H(Q) / η(Q)
```

**Or Polynomial Fit:**
```
P(Q) = c₀ + c₁·Q + c₂·Q² + c₃·Q³
```

### NPSH Curve Fitting

**Power Law Model:**
```
NPSH_req(Q) = NPSH₀ + k·Q^n
```
Where n ≈ 1.5-2.0

**Exponential Model:**
```
NPSH_req(Q) = a·exp(b·Q) + c
```

### Uncertainty Analysis

**Measurement Uncertainties:**
- Flow: ±2-5%
- Head: ±1-3%
- Power: ±1-2%
- Efficiency: ±2-5% (compound error)

**Propagation:**
```
σ_η = η·√((σ_Q/Q)² + (σ_H/H)² + (σ_P/P)²)
```

## System Curve Intersection (Operating Point)

### System Curve Definition

**System Head Equation:**
```
H_system = H_static + K·Q²
```

Where:
- H_static = elevation difference + pressure head
- K = friction coefficient = (8·f·L)/(π²·g·D⁵)

**Components:**
```
H_static = (z₂ - z₁) + (P₂ - P₁)/(ρ·g)
H_friction = (f·L·v²)/(D·2g) = K·Q²
H_minor = Σ(K_i·v²/(2g))
```

### Operating Point Determination

**Graphical Method:**
1. Plot H-Q curve (pump)
2. Plot H_system curve
3. Intersection = operating point

**Analytical Method:**
```
H_pump(Q) = H_system(Q)
a₀ + a₁·Q + a₂·Q² = H_s + K·Q²
Q_op = (-a₁ + √(a₁² - 4·(a₂-K)·(a₀-H_s)))/(2·(a₂-K))
```

**Numerical Method (Newton-Raphson):**
```python
def find_operating_point(pump_curve, system_curve, Q_initial=1.0):
    f = lambda Q: pump_curve(Q) - system_curve(Q)
    Q_op = newton(f, Q_initial)
    H_op = pump_curve(Q_op)
    return Q_op, H_op
```

### Multiple Pump Configurations

**Pumps in Series:**
```
H_total = H₁ + H₂
Q_total = Q  (same flow through both)
```
Use: High head applications

**Pumps in Parallel:**
```
Q_total = Q₁ + Q₂
H_total = H  (same head for both)
```
Use: High flow applications

**Parallel Operation Stability:**
- Requires H-Q curves with dH/dQ < 0
- Unequal flow distribution if curves differ
- One pump may operate at shutoff

### Variable Speed Operation

**Operating Point Trajectory:**
```
H/H₀ = (Q/Q₀)²  (system curve)
H/H₀ = (N/N₀)²  (affinity law)
```

Combined:
```
Q/Q₀ = N/N₀
```

**Energy Savings:**
```
P₂/P₁ = (Q₂/Q₁)³  (for system curve flow)
```

Reducing flow by 20% saves ~50% power!

## Design Guidelines and Best Practices

### Curve Selection

1. **Select pump with BEP near duty point**
   - Target: Q_duty = 0.8-1.1·Q_BEP
   - Margin for fouling and wear

2. **Verify operating range**
   - All expected flows within 50-150% Q_BEP
   - Consider seasonal variations
   - Account for future expansion

3. **Check curve stability**
   - Ensure dH/dQ < 0 throughout
   - Avoid dips or humps
   - Important for parallel operation

4. **Validate NPSH margin**
   - NPSH_avail > 1.3·NPSH_req (minimum)
   - Higher margin at elevated temperatures
   - Consider transient conditions

### Common Mistakes to Avoid

**1. Ignoring Efficiency:**
- ❌ Selecting pump based only on H-Q intersection
- ✅ Verify operating point near BEP
- ✅ Calculate lifecycle energy costs

**2. Neglecting Operating Range:**
- ❌ Designing for single duty point
- ✅ Consider full operating envelope
- ✅ Verify performance at extremes

**3. Misapplying Affinity Laws:**
- ❌ Extrapolating beyond ±20% speed change
- ❌ Ignoring efficiency changes with trimming
- ✅ Use manufacturer curves when available
- ✅ Account for Reynolds number effects

**4. System Curve Errors:**
- ❌ Using only friction losses (forgetting static head)
- ❌ Assuming constant K with control valves
- ✅ Include all system resistances
- ✅ Account for changing conditions

**5. Parallel Pump Design:**
- ❌ Assuming 2 pumps = 2× flow
- ❌ Ignoring individual pump operating points
- ✅ Plot combined H-Q curve
- ✅ Verify stable operation for all combinations

## Verification Checklist

Before finalizing pump selection:

- [ ] Operating point within 70-120% of BEP flow
- [ ] Efficiency at duty point > 90% of maximum efficiency
- [ ] H-Q curve stable (dH/dQ < 0) throughout operating range
- [ ] NPSH_avail > 1.5·NPSH_req at all operating points
- [ ] Power curve non-overloading (or motor sized appropriately)
- [ ] System curve accounts for all resistances
- [ ] Future operating conditions considered
- [ ] Affinity laws applied correctly (if used)
- [ ] Curve fitting R² > 0.99 (if fitting data)
- [ ] Parallel/series operation verified (if applicable)

## Advanced Topics

### Variable Frequency Drive (VFD) Operation

**Speed Control Advantage:**
```
Power savings = 1 - (N₂/N₁)³
```

**Considerations:**
- Bearing minimum speed limits
- Cooling at reduced flow
- VFD efficiency (typically 95-97%)

### Specific Speed Effects

**Curve Shape vs. Specific Speed:**
```
Ns = N·√Q / H^0.75
```

- Low Ns (10-30): Steep H-Q curve, flat efficiency
- Medium Ns (30-80): Moderate slope, peak efficiency
- High Ns (80-150): Flat H-Q curve, narrow efficiency peak

### Temperature Effects

**Viscosity Impact:**
- Reynolds number: Re = ρ·v·D/μ
- Head reduction: H_viscous < H_water
- Efficiency penalty increases at low Re
- Use HI correction charts for Re < 10⁶

**Vapor Pressure:**
- NPSH_req increases with temperature
- Critical for hot liquids
- May require booster pumps

## Summary

Pump performance curves are essential tools for:
1. Pump selection and sizing
2. Operating point prediction
3. Energy consumption analysis
4. System optimization
5. Troubleshooting and diagnostics

Mastery requires understanding:
- Physical meaning of each curve type
- Proper application of affinity laws
- System-pump interaction
- Operating range considerations
- Impact of off-design operation

## References

See `reference.md` for detailed standards, equations, and literature references.
See `plotter.py` for verified computational examples.
