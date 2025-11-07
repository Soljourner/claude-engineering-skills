---
name: fluids-package
description: "Pipe flow, pump sizing, friction factor, and compressible flow calculations"
category: packages
domain: fluids
complexity: intermediate
dependencies:
  - fluids
  - scipy
---

# Fluids Package Skill

## Overview

The `fluids` library is a comprehensive Python package for mechanical and chemical engineers working with fluid flow problems. It provides validated correlations and functions for:

- Pipe flow and friction factor calculations
- Pump sizing and performance analysis
- Compressible and incompressible flow
- Heat exchanger design
- Two-phase flow
- Pressure drop calculations
- Fluid properties

The library implements over 100 correlations from the literature with extensive validation against published test cases.

## Installation

```bash
pip install fluids
```

For full functionality including optimization routines:

```bash
pip install fluids[complete]
```

## Key Modules

### fluids.core
Core utilities and dimensional analysis functions.

### fluids.friction
Friction factor calculations for pipe flow including:
- Darcy-Weisbach equation
- Colebrook-White correlation
- Moody diagram implementations
- Turbulent and laminar flow regimes

### fluids.pump
Pump performance calculations:
- Affinity laws
- Specific speed
- NPSH calculations
- Pump curves and efficiency

### fluids.compressible
Compressible flow calculations:
- Mach number relationships
- Choked flow conditions
- Isentropic flow
- Normal shock waves

### fluids.fittings
Pressure drop through valves, fittings, and pipe components.

## Common Functions with Engineering Context

### Reynolds Number Calculation

The Reynolds number (Re) determines flow regime and is fundamental to all pipe flow calculations.

```python
from fluids.core import Reynolds

# For pipe flow: Re = ρVD/μ
Re = Reynolds(V=2.5, D=0.05, rho=1000, mu=0.001)
# Result: 125000 (turbulent flow)

# Interpretation:
# Re < 2300: Laminar flow
# 2300 < Re < 4000: Transition
# Re > 4000: Turbulent flow
```

### Friction Factor (Darcy-Weisbach)

The friction factor (f) is used in the Darcy-Weisbach equation: ΔP = f(L/D)(ρV²/2)

```python
from fluids.friction import friction_factor

# Colebrook-White correlation (implicit, most accurate)
f = friction_factor(Re=125000, eD=0.0001)  # eD = roughness/diameter
# Result: ~0.0178

# Moody correlation (explicit approximation)
from fluids.friction import friction_factor_Moody
f_moody = friction_factor_Moody(Re=125000, eD=0.0001)

# For laminar flow (Re < 2300):
f_laminar = friction_factor(Re=1500, eD=0.0001)
# Result: 0.0427 (f = 64/Re)
```

### Head Loss in Pipes

Calculate pressure drop and head loss in piping systems.

```python
from fluids.friction import friction_factor, head_from_P
from fluids.core import Reynolds

# Given: Water flow through steel pipe
D = 0.1  # m, pipe diameter
L = 100  # m, pipe length
V = 2.0  # m/s, velocity
rho = 1000  # kg/m³, density
mu = 0.001  # Pa·s, viscosity
epsilon = 0.000045  # m, roughness (steel)

# Step 1: Calculate Reynolds number
Re = Reynolds(V=V, D=D, rho=rho, mu=mu)

# Step 2: Calculate friction factor
eD = epsilon / D
f = friction_factor(Re=Re, eD=eD)

# Step 3: Calculate pressure drop
# Darcy-Weisbach: ΔP = f(L/D)(ρV²/2)
dP = f * (L/D) * (rho * V**2 / 2)

# Step 4: Convert to head loss
h_loss = head_from_P(dP, rho)  # meters of fluid

print(f"Reynolds: {Re:.0f}")
print(f"Friction factor: {f:.5f}")
print(f"Pressure drop: {dP:.0f} Pa")
print(f"Head loss: {h_loss:.2f} m")
```

### Pump Affinity Laws

Relate pump performance at different speeds and impeller diameters.

```python
from fluids.pump import affinity_law_volume, affinity_law_head, affinity_law_power

# Original pump operating point
Q1 = 100  # m³/h, flow rate
H1 = 50   # m, head
P1 = 20   # kW, power
N1 = 1450 # rpm, speed
D1 = 0.3  # m, impeller diameter

# New speed
N2 = 1750  # rpm

# Affinity laws (constant impeller diameter):
# Q2/Q1 = N2/N1
# H2/H1 = (N2/N1)²
# P2/P1 = (N2/N1)³

Q2 = affinity_law_volume(Q1, N1, N2)
H2 = affinity_law_head(H1, N1, N2)
P2 = affinity_law_power(P1, N1, N2)

print(f"New flow: {Q2:.1f} m³/h")
print(f"New head: {H2:.1f} m")
print(f"New power: {P2:.1f} kW")
```

### Specific Speed Calculations

Specific speed (Ns) characterizes pump type and efficiency.

```python
from fluids.pump import specific_speed

# Pump operating conditions
Q = 0.05  # m³/s, flow rate
H = 40    # m, head
N = 1450  # rpm, rotational speed

# Calculate specific speed (dimensionless)
Ns = specific_speed(Q, H, N)

# Interpretation:
# Ns < 0.5: Centrifugal (radial flow)
# 0.5 < Ns < 1.0: Francis (mixed flow)
# 1.0 < Ns < 4.0: Propeller (axial flow)

print(f"Specific speed: {Ns:.2f}")
if Ns < 0.5:
    pump_type = "Centrifugal (radial flow)"
elif Ns < 1.0:
    pump_type = "Francis (mixed flow)"
else:
    pump_type = "Propeller (axial flow)"
print(f"Recommended pump type: {pump_type}")
```

### Compressible Flow - Mach Number

For gas flow in pipes and nozzles.

```python
from fluids.compressible import Mach

# Calculate Mach number from velocity
V = 200  # m/s, velocity
c = 340  # m/s, speed of sound in air at 15°C

Ma = Mach(V, c)
# Result: 0.588

# Flow classification:
# Ma < 0.3: Incompressible
# 0.3 < Ma < 0.8: Subsonic
# 0.8 < Ma < 1.2: Transonic
# Ma > 1.2: Supersonic
```

### Compressible Flow - Choked Flow

Determine if flow is choked in a nozzle or orifice.

```python
from fluids.compressible import P_critical_flow

# Gas properties
P_upstream = 500000  # Pa, upstream pressure
k = 1.4  # heat capacity ratio (air)

# Critical pressure for choked flow
P_crit = P_critical_flow(P=P_upstream, k=k)
# Result: ~264,000 Pa

# If downstream pressure < P_crit, flow is choked
P_downstream = 200000  # Pa
if P_downstream < P_crit:
    print("Flow is choked - mass flow is at maximum")
    print(f"Critical pressure: {P_crit:.0f} Pa")
else:
    print("Flow is not choked")
```

## Complete Engineering Examples

### Example 1: Pump Selection and System Curve

```python
import numpy as np
from fluids.friction import friction_factor
from fluids.core import Reynolds
import matplotlib.pyplot as plt

def system_curve(Q_range, static_head, pipe_specs):
    """
    Calculate system head curve for a piping system.

    Parameters:
    -----------
    Q_range : array, flow rates (m³/s)
    static_head : float, static lift (m)
    pipe_specs : dict with keys:
        - L: pipe length (m)
        - D: pipe diameter (m)
        - epsilon: roughness (m)
        - rho: fluid density (kg/m³)
        - mu: fluid viscosity (Pa·s)

    Returns:
    --------
    H_system : array, required head at each flow rate (m)
    """
    L = pipe_specs['L']
    D = pipe_specs['D']
    rho = pipe_specs['rho']
    mu = pipe_specs['mu']
    epsilon = pipe_specs['epsilon']

    # Calculate cross-sectional area
    A = np.pi * D**2 / 4

    H_system = np.zeros_like(Q_range)

    for i, Q in enumerate(Q_range):
        if Q == 0:
            H_system[i] = static_head
            continue

        # Calculate velocity
        V = Q / A

        # Reynolds number
        Re = Reynolds(V=V, D=D, rho=rho, mu=mu)

        # Friction factor
        eD = epsilon / D
        f = friction_factor(Re=Re, eD=eD)

        # Friction head loss (Darcy-Weisbach)
        h_friction = f * (L/D) * (V**2 / (2*9.81))

        # Total system head
        H_system[i] = static_head + h_friction

    return H_system

# Define system
pipe_specs = {
    'L': 200,  # m
    'D': 0.15,  # m
    'epsilon': 0.000045,  # m (steel)
    'rho': 1000,  # kg/m³
    'mu': 0.001  # Pa·s
}
static_head = 30  # m

# Generate system curve
Q_range = np.linspace(0, 0.1, 50)  # m³/s
H_system = system_curve(Q_range, static_head, pipe_specs)

# Example pump curve (quadratic fit)
# H = H0 - A*Q - B*Q²
H0 = 80  # m, shutoff head
A = 200  # head loss coefficient
B = 3000  # quadratic coefficient
H_pump = H0 - A*Q_range - B*Q_range**2

# Find operating point (intersection)
idx = np.argmin(np.abs(H_pump - H_system))
Q_op = Q_range[idx]
H_op = H_system[idx]

print(f"Operating Point:")
print(f"  Flow rate: {Q_op*3600:.1f} m³/h ({Q_op:.4f} m³/s)")
print(f"  Head: {H_op:.1f} m")

# Verification: At Q=0.05 m³/s
# V = Q/A = 0.05/(π*0.15²/4) = 2.83 m/s
# Re = 1000*2.83*0.15/0.001 = 424,500 (turbulent)
# f ≈ 0.0175 (from Moody chart)
# h_f = 0.0175*(200/0.15)*(2.83²/(2*9.81)) = 9.5 m
# H_total = 30 + 9.5 = 39.5 m ✓
```

### Example 2: Parallel Pump Configuration

```python
from fluids.pump import affinity_law_volume, affinity_law_head

def parallel_pumps(Q_total, n_pumps, single_pump_curve):
    """
    Calculate operating point for parallel pump configuration.

    For pumps in parallel:
    - Flow rates add: Q_total = n * Q_single
    - Head remains the same: H_total = H_single
    """
    # Single pump flow rate
    Q_single = Q_total / n_pumps

    # Head from single pump curve
    H = single_pump_curve(Q_single)

    return Q_single, H

# Single pump curve: H = 60 - 500*Q² (simplified)
def pump_curve(Q):
    return 60 - 500*Q**2

# System requires 150 m³/h at 40 m head
Q_required = 150/3600  # m³/s
n_pumps = 2

Q_single, H_operating = parallel_pumps(Q_required, n_pumps, pump_curve)

print(f"Parallel Pump Configuration ({n_pumps} pumps):")
print(f"  Total flow: {Q_required*3600:.1f} m³/h")
print(f"  Flow per pump: {Q_single*3600:.1f} m³/h")
print(f"  Operating head: {H_operating:.1f} m")

# Verification:
# Each pump delivers 75 m³/h (0.0208 m³/s)
# H = 60 - 500*(0.0208)² = 60 - 0.22 = 59.8 m ✓
```

### Example 3: Friction Factor Validation

```python
from fluids.friction import friction_factor, friction_factor_laminar

# Test Case 1: Laminar Flow (Poiseuille)
# Analytical solution: f = 64/Re
Re_laminar = 1000
f_calculated = friction_factor(Re=Re_laminar, eD=0)
f_analytical = 64/Re_laminar

print("Test 1: Laminar Flow")
print(f"  Re = {Re_laminar}")
print(f"  f (calculated) = {f_calculated:.6f}")
print(f"  f (analytical) = {f_analytical:.6f}")
print(f"  Error = {abs(f_calculated - f_analytical):.9f}")
assert abs(f_calculated - f_analytical) < 1e-9, "Laminar flow test failed"
print("  ✓ PASSED\n")

# Test Case 2: Turbulent Flow - Smooth Pipe
# From Moody diagram: Re=1e5, smooth pipe → f ≈ 0.0183
Re_turbulent = 1e5
f_smooth = friction_factor(Re=Re_turbulent, eD=0)

print("Test 2: Turbulent Flow (Smooth Pipe)")
print(f"  Re = {Re_turbulent:.0f}")
print(f"  f (calculated) = {f_smooth:.6f}")
print(f"  f (Moody chart) ≈ 0.0183")
print(f"  Error = {abs(f_smooth - 0.0183):.6f}")
assert abs(f_smooth - 0.0183) < 0.0005, "Smooth pipe test failed"
print("  ✓ PASSED\n")

# Test Case 3: Turbulent Flow - Rough Pipe
# Commercial steel: ε/D = 0.0002, Re = 1e6 → f ≈ 0.0144
Re_rough = 1e6
eD_rough = 0.0002
f_rough = friction_factor(Re=Re_rough, eD=eD_rough)

print("Test 3: Turbulent Flow (Rough Pipe)")
print(f"  Re = {Re_rough:.0f}")
print(f"  ε/D = {eD_rough}")
print(f"  f (calculated) = {f_rough:.6f}")
print(f"  f (Moody chart) ≈ 0.0144")
print(f"  Error = {abs(f_rough - 0.0144):.6f}")
assert abs(f_rough - 0.0144) < 0.0005, "Rough pipe test failed"
print("  ✓ PASSED\n")

# Test Case 4: From Crane TP-410 (Standard reference)
# Example 4-9: Water in 6-inch Schedule 40 pipe
D = 0.1541  # m (6-inch Schedule 40 ID)
V = 3.05  # m/s
rho = 998  # kg/m³
mu = 0.001  # Pa·s
epsilon = 0.000045  # m (steel)

Re_crane = Reynolds(V=V, D=D, rho=rho, mu=mu)
eD_crane = epsilon/D
f_crane = friction_factor(Re=Re_crane, eD=eD_crane)

# Crane TP-410 gives f ≈ 0.0172
print("Test 4: Crane TP-410 Example 4-9")
print(f"  Re = {Re_crane:.0f}")
print(f"  ε/D = {eD_crane:.6f}")
print(f"  f (calculated) = {f_crane:.6f}")
print(f"  f (Crane TP-410) ≈ 0.0172")
print(f"  Error = {abs(f_crane - 0.0172):.6f}")
assert abs(f_crane - 0.0172) < 0.0003, "Crane TP-410 test failed"
print("  ✓ PASSED\n")

print("All friction factor tests passed! ✓")
```

### Example 4: Compressible Flow in Pipe

```python
from fluids.compressible import isothermal_gas

# Gas pipeline calculation
# Problem: Natural gas (methane) pipeline
P1 = 5e6  # Pa, inlet pressure
T = 288.15  # K, temperature (15°C)
L = 50000  # m, pipeline length (50 km)
D = 0.5  # m, diameter
m = 10  # kg/s, mass flow rate
MW = 16.04  # g/mol, molecular weight (CH4)
k = 1.31  # heat capacity ratio

# Calculate outlet pressure using isothermal flow
# This accounts for friction and compressibility
from fluids.compressible import isothermal_gas
P2 = isothermal_gas(rho=None, P1=P1, P2=None, L=L, D=D, m=m,
                    T=T, Z=1, fd=0.012)  # Assuming f=0.012

print("Gas Pipeline Calculation:")
print(f"  Inlet pressure: {P1/1e6:.2f} MPa")
print(f"  Outlet pressure: {P2/1e6:.2f} MPa")
print(f"  Pressure drop: {(P1-P2)/1e6:.2f} MPa")
print(f"  Length: {L/1000:.0f} km")
print(f"  Mass flow: {m:.1f} kg/s")
```

## Best Practices

1. **Always check Reynolds number** before selecting friction factor correlation
2. **Validate results** against known test cases or handbook values
3. **Use consistent units** - fluids uses SI units throughout
4. **Consider safety factors** when sizing pumps (typically 10-15% margin)
5. **Check physical limits** - verify Mach number < 0.3 for incompressible assumption
6. **Use appropriate roughness values** from literature (Moody, Crane TP-410)

## Common Pipe Roughness Values

| Material | Roughness ε (m) |
|----------|----------------|
| Drawn tubing | 0.0000015 |
| Commercial steel | 0.000045 |
| Galvanized iron | 0.00015 |
| Cast iron | 0.00026 |
| Concrete | 0.0003 to 0.003 |
| Riveted steel | 0.0009 to 0.009 |

## Troubleshooting

### Issue: Friction factor doesn't converge
- **Cause**: Invalid Reynolds number or ε/D ratio
- **Solution**: Check that Re > 0 and 0 ≤ ε/D < 0.5

### Issue: Pump curves don't intersect system curve
- **Cause**: Pump undersized or oversized
- **Solution**: Adjust pump speed using affinity laws or select different pump

### Issue: Compressible flow results unrealistic
- **Cause**: Mach number > 0.3 but using incompressible equations
- **Solution**: Use compressible flow functions from fluids.compressible

## References

- Crane Technical Paper 410 (TP-410): "Flow of Fluids Through Valves, Fittings, and Pipe"
- Moody, L.F. (1944): "Friction factors for pipe flow"
- Colebrook, C.F. (1939): "Turbulent flow in pipes"
- Karassik's Pump Handbook (4th Edition)
- GPSA Engineering Data Book (14th Edition)

## Further Reading

- Official documentation: https://fluids.readthedocs.io/
- Source code: https://github.com/CalebBell/fluids
- Chemical Engineering Design Library: https://chemicals.readthedocs.io/
