# Fluids Package Reference

## Quick Reference Guide

This reference provides detailed information on the fluids library functions, correlations, and literature sources.

---

## Core Functions (fluids.core)

### Reynolds Number

```python
Reynolds(V, D, rho, mu)
```

**Description**: Calculates Reynolds number for pipe flow

**Parameters**:
- `V` (float): Velocity (m/s)
- `D` (float): Diameter (m)
- `rho` (float): Density (kg/m³)
- `mu` (float): Dynamic viscosity (Pa·s)

**Returns**: Reynolds number (dimensionless)

**Formula**: Re = ρVD/μ

**Flow Regimes**:
- Re < 2300: Laminar
- 2300 < Re < 4000: Transition
- Re > 4000: Turbulent

**Example**:
```python
Re = Reynolds(V=2.5, D=0.05, rho=1000, mu=0.001)
# Result: 125000
```

---

### Velocity from Flow Rate

```python
v_from_Q(Q, D)
```

**Description**: Calculate velocity from volumetric flow rate

**Parameters**:
- `Q` (float): Volumetric flow rate (m³/s)
- `D` (float): Pipe diameter (m)

**Returns**: Velocity (m/s)

**Formula**: V = Q / (πD²/4)

---

## Friction Functions (fluids.friction)

### Friction Factor (General)

```python
friction_factor(Re, eD=0.0)
```

**Description**: Calculate Darcy friction factor using the most appropriate correlation

**Parameters**:
- `Re` (float): Reynolds number
- `eD` (float): Relative roughness (ε/D)

**Returns**: Darcy friction factor (dimensionless)

**Correlations Used**:
- Laminar (Re < 2300): f = 64/Re
- Turbulent: Colebrook-White equation (implicit solution)

**Formula (Turbulent)**:
```
1/√f = -2.0 log₁₀(ε/(3.7D) + 2.51/(Re√f))
```

**Example**:
```python
f = friction_factor(Re=100000, eD=0.0001)
# Result: ~0.0178
```

---

### Colebrook-White Equation

```python
friction_factor_Colebrook(Re, eD)
```

**Description**: Implicit Colebrook-White equation for turbulent friction factor

**Valid Range**:
- Re > 4000
- 0 ≤ ε/D ≤ 0.05

**Accuracy**: Most accurate correlation for turbulent flow

**Reference**: Colebrook, C. F. (1939). "Turbulent flow in pipes, with particular reference to the transition region between the smooth and rough pipe laws". Journal of the Institution of Civil Engineers.

---

### Moody Correlation

```python
friction_factor_Moody(Re, eD)
```

**Description**: Explicit approximation to Colebrook-White equation

**Valid Range**:
- Re > 4000
- 0 ≤ ε/D ≤ 0.01

**Accuracy**: Within 5% of Colebrook-White

**Formula**:
```
f = 0.0055[1 + (20000ε/D + 10⁶/Re)^(1/3)]
```

**Reference**: Moody, L. F. (1944). "Friction factors for pipe flow". Transactions of the ASME.

---

### Swamee-Jain Correlation

```python
friction_factor_Swamee_Jain(Re, eD)
```

**Description**: Explicit approximation, widely used

**Valid Range**:
- 5000 < Re < 10⁸
- 0.00001 < ε/D < 0.01

**Accuracy**: Within 1% of Colebrook-White

**Formula**:
```
f = 0.25 / [log₁₀(ε/(3.7D) + 5.74/Re^0.9)]²
```

**Reference**: Swamee, P. K.; Jain, A. K. (1976). "Explicit equations for pipe-flow problems". Journal of the Hydraulics Division.

---

### Laminar Flow Friction Factor

```python
friction_factor_laminar(Re)
```

**Description**: Friction factor for laminar flow in circular pipes

**Valid Range**: Re < 2300

**Formula**: f = 64/Re

**Derivation**: From Hagen-Poiseuille equation

---

### Head Loss from Pressure Drop

```python
head_from_P(dP, rho)
```

**Description**: Convert pressure drop to head loss

**Parameters**:
- `dP` (float): Pressure drop (Pa)
- `rho` (float): Density (kg/m³)

**Returns**: Head loss (m)

**Formula**: h = ΔP / (ρg)

---

## Pump Functions (fluids.pump)

### Affinity Law - Volume

```python
affinity_law_volume(Q1, N1, N2, D1=None, D2=None)
```

**Description**: Calculate new flow rate from pump speed or diameter change

**Parameters**:
- `Q1` (float): Original flow rate (m³/s)
- `N1` (float): Original speed (rpm)
- `N2` (float): New speed (rpm)
- `D1` (float, optional): Original impeller diameter (m)
- `D2` (float, optional): New impeller diameter (m)

**Returns**: New flow rate Q2 (m³/s)

**Formulas**:
- Constant diameter: Q₂/Q₁ = N₂/N₁
- Constant speed: Q₂/Q₁ = D₂/D₁

---

### Affinity Law - Head

```python
affinity_law_head(H1, N1, N2, D1=None, D2=None)
```

**Description**: Calculate new head from pump speed or diameter change

**Parameters**: Same as affinity_law_volume

**Returns**: New head H2 (m)

**Formulas**:
- Constant diameter: H₂/H₁ = (N₂/N₁)²
- Constant speed: H₂/H₁ = (D₂/D₁)²

---

### Affinity Law - Power

```python
affinity_law_power(P1, N1, N2, D1=None, D2=None)
```

**Description**: Calculate new power from pump speed or diameter change

**Parameters**: Same as affinity_law_volume

**Returns**: New power P2 (W)

**Formulas**:
- Constant diameter: P₂/P₁ = (N₂/N₁)³
- Constant speed: P₂/P₁ = (D₂/D₁)³

---

### Specific Speed

```python
specific_speed(Q, H, N)
```

**Description**: Calculate pump specific speed (dimensionless)

**Parameters**:
- `Q` (float): Flow rate at BEP (m³/s)
- `H` (float): Head at BEP (m)
- `N` (float): Rotational speed (rpm)

**Returns**: Specific speed Ns (dimensionless)

**Formula**: Ns = N√Q / H^0.75

**Pump Type Selection**:
| Ns Range | Pump Type | Characteristics |
|----------|-----------|-----------------|
| < 0.5 | Centrifugal (radial) | High head, low flow |
| 0.5 - 1.0 | Francis (mixed) | Medium head, medium flow |
| > 1.0 | Propeller (axial) | Low head, high flow |

**Reference**: Karassik, I. J. et al. (2008). "Pump Handbook", 4th Edition.

---

### NPSH Available

```python
NPSH_available(P, Pvp, rho, z)
```

**Description**: Calculate Net Positive Suction Head Available

**Parameters**:
- `P` (float): Absolute pressure at pump inlet (Pa)
- `Pvp` (float): Vapor pressure of fluid (Pa)
- `rho` (float): Fluid density (kg/m³)
- `z` (float): Elevation of pump above source (m)

**Returns**: NPSH available (m)

**Formula**: NPSHa = (P - Pvp)/(ρg) - z

**Critical Requirement**: NPSHa > NPSHr (required) + safety margin (0.5-1.0 m)

---

## Compressible Flow Functions (fluids.compressible)

### Mach Number

```python
Mach(V, c)
```

**Description**: Calculate Mach number

**Parameters**:
- `V` (float): Flow velocity (m/s)
- `c` (float): Speed of sound (m/s)

**Returns**: Mach number (dimensionless)

**Formula**: Ma = V/c

**Speed of sound in ideal gas**: c = √(kRT/M)
- k: heat capacity ratio
- R: gas constant (8314 J/kmol·K)
- T: temperature (K)
- M: molecular weight (kg/kmol)

---

### Critical Pressure (Choked Flow)

```python
P_critical_flow(P, k)
```

**Description**: Calculate critical pressure ratio for choked flow

**Parameters**:
- `P` (float): Upstream pressure (Pa)
- `k` (float): Heat capacity ratio (γ = cp/cv)

**Returns**: Critical pressure (Pa)

**Formula**: P_crit/P = (2/(k+1))^(k/(k-1))

**For air (k=1.4)**: P_crit/P ≈ 0.528

**Choked Flow Condition**: Flow is choked when P_downstream < P_critical

---

### Isentropic Flow Relations

```python
isentropic_flow(Ma, k)
```

**Description**: Calculate isentropic flow properties

**Parameters**:
- `Ma` (float): Mach number
- `k` (float): Heat capacity ratio

**Returns**: Dictionary with:
- T/T₀: Temperature ratio
- P/P₀: Pressure ratio
- ρ/ρ₀: Density ratio

**Formulas**:
```
T/T₀ = 1 / (1 + (k-1)/2 × Ma²)
P/P₀ = (T/T₀)^(k/(k-1))
ρ/ρ₀ = (T/T₀)^(1/(k-1))
```

---

## Fitting Loss Coefficients (fluids.fittings)

### K Factor (Loss Coefficient)

```python
K_gate_valve(D1, D2, angle)
```

**Description**: Loss coefficient for gate valve

**Formula**: ΔP = K × (ρV²/2)

**Common K Values**:

| Fitting | K (fully open) |
|---------|----------------|
| Gate valve | 0.15 |
| Globe valve | 10.0 |
| Ball valve | 0.05 |
| 90° elbow (standard) | 0.30 |
| 90° elbow (long radius) | 0.20 |
| 45° elbow | 0.16 |
| Tee (straight through) | 0.10 |
| Tee (branch flow) | 0.60 |
| Pipe entrance (sharp) | 0.50 |
| Pipe entrance (rounded) | 0.04 |
| Pipe exit | 1.00 |

**Reference**: Crane Technical Paper No. 410 (TP-410), "Flow of Fluids Through Valves, Fittings, and Pipe"

---

## Common Correlations

### Pipe Roughness Values

**Absolute Roughness (ε)**:

| Material | Roughness (mm) | Roughness (m) |
|----------|----------------|----------------|
| Drawn tubing (brass, copper) | 0.0015 | 1.5 × 10⁻⁶ |
| Commercial steel, wrought iron | 0.045 | 4.5 × 10⁻⁵ |
| Galvanized iron | 0.15 | 1.5 × 10⁻⁴ |
| Cast iron (new) | 0.26 | 2.6 × 10⁻⁴ |
| Cast iron (old, corroded) | 1.5 | 1.5 × 10⁻³ |
| Concrete (smooth) | 0.30 | 3.0 × 10⁻⁴ |
| Concrete (rough) | 3.0 | 3.0 × 10⁻³ |
| Riveted steel | 0.9 to 9.0 | 9 × 10⁻⁴ to 9 × 10⁻³ |

**Reference**: Moody, L. F. (1944). "Friction factors for pipe flow"

---

### Fluid Properties (Water at 20°C)

- Density (ρ): 998.2 kg/m³
- Dynamic viscosity (μ): 1.002 × 10⁻³ Pa·s
- Kinematic viscosity (ν): 1.004 × 10⁻⁶ m²/s
- Vapor pressure: 2.34 kPa

### Fluid Properties (Air at 20°C, 1 atm)

- Density (ρ): 1.204 kg/m³
- Dynamic viscosity (μ): 1.82 × 10⁻⁵ Pa·s
- Kinematic viscosity (ν): 1.51 × 10⁻⁵ m²/s
- Heat capacity ratio (k): 1.40

---

## Dimensionless Numbers

### Reynolds Number (Re)
- **Definition**: Ratio of inertial to viscous forces
- **Formula**: Re = ρVD/μ = VD/ν
- **Significance**: Determines flow regime

### Froude Number (Fr)
- **Definition**: Ratio of inertial to gravitational forces
- **Formula**: Fr = V/√(gD)
- **Significance**: Important in open channel flow

### Weber Number (We)
- **Definition**: Ratio of inertial to surface tension forces
- **Formula**: We = ρV²L/σ
- **Significance**: Important in droplet/bubble formation

### Mach Number (Ma)
- **Definition**: Ratio of flow velocity to speed of sound
- **Formula**: Ma = V/c
- **Significance**: Determines compressibility effects

---

## Literature Sources

### Primary References

1. **Crane Technical Paper 410 (TP-410)**
   - "Flow of Fluids Through Valves, Fittings, and Pipe"
   - Industry standard for K factors and pressure drop calculations
   - Available: https://www.flowoffluids.com/

2. **Moody, L. F. (1944)**
   - "Friction factors for pipe flow"
   - Transactions of the ASME, Vol. 66, pp. 671-684
   - Original Moody diagram publication

3. **Colebrook, C. F. (1939)**
   - "Turbulent flow in pipes, with particular reference to the transition region"
   - Journal of the Institution of Civil Engineers, Vol. 11, pp. 133-156
   - Fundamental turbulent flow correlation

4. **Swamee, P. K.; Jain, A. K. (1976)**
   - "Explicit equations for pipe-flow problems"
   - Journal of the Hydraulics Division, ASCE, Vol. 102, No. HY5, pp. 657-664

5. **Karassik's Pump Handbook (4th Edition)**
   - Karassik, I. J.; Messina, J. P.; Cooper, P.; Heald, C. C. (2008)
   - McGraw-Hill Professional
   - ISBN: 978-0071460446
   - Comprehensive pump design and selection guide

### Secondary References

6. **White, F. M. (2016)**
   - "Fluid Mechanics" (8th Edition)
   - McGraw-Hill Education
   - ISBN: 978-0073398273
   - Fundamental fluid mechanics textbook

7. **Munson, B. R. et al. (2013)**
   - "Fundamentals of Fluid Mechanics" (7th Edition)
   - Wiley
   - ISBN: 978-1118116135

8. **Perry's Chemical Engineers' Handbook (9th Edition)**
   - Green, D. W.; Perry, R. H. (2019)
   - McGraw-Hill Education
   - ISBN: 978-0071834087
   - Comprehensive chemical engineering reference

9. **GPSA Engineering Data Book (14th Edition)**
   - Gas Processors Suppliers Association (2012)
   - Industry standard for gas processing

---

## Online Resources

### Official Documentation

- **Fluids Documentation**: https://fluids.readthedocs.io/
  - Complete API reference
  - Examples and tutorials
  - Installation guide

- **Fluids GitHub Repository**: https://github.com/CalebBell/fluids
  - Source code
  - Issue tracking
  - Contributing guidelines

### Related Libraries

- **Chemicals Library**: https://chemicals.readthedocs.io/
  - Thermodynamic properties
  - Component data
  - Integration with fluids

- **ThermoPy**: https://github.com/guillemborrell/thermopy
  - Thermodynamic calculations
  - Property estimation

### Standards Organizations

- **ASME (American Society of Mechanical Engineers)**: https://www.asme.org/
  - Pump standards
  - Piping codes

- **HI (Hydraulic Institute)**: https://www.pumps.org/
  - Pump standards and testing
  - Best practices

- **API (American Petroleum Institute)**: https://www.api.org/
  - Pipeline standards
  - Equipment specifications

---

## Validation and Testing

The fluids library includes extensive validation against:

1. **Analytical Solutions**
   - Laminar flow (Hagen-Poiseuille)
   - Fully developed turbulent flow

2. **Standard Test Cases**
   - Crane TP-410 examples
   - Moody diagram values
   - ASME standards

3. **Published Data**
   - Experimental measurements
   - CFD validations
   - Industry benchmarks

4. **Cross-Validation**
   - Comparison with commercial software
   - Multiple correlation methods
   - Uncertainty quantification

---

## Units and Conversions

### Standard SI Units (Used by Fluids)

| Quantity | SI Unit | Symbol |
|----------|---------|--------|
| Length | meter | m |
| Mass | kilogram | kg |
| Time | second | s |
| Temperature | Kelvin | K |
| Pressure | Pascal | Pa |
| Density | kg/m³ | kg/m³ |
| Viscosity (dynamic) | Pascal-second | Pa·s |
| Viscosity (kinematic) | m²/s | m²/s |
| Velocity | m/s | m/s |
| Flow rate | m³/s | m³/s |

### Common Conversions

**Flow Rate**:
- 1 m³/s = 3600 m³/h
- 1 m³/h = 4.403 US gpm
- 1 L/s = 0.001 m³/s

**Pressure**:
- 1 bar = 100,000 Pa = 100 kPa
- 1 psi = 6894.76 Pa
- 1 atm = 101,325 Pa

**Viscosity (Dynamic)**:
- 1 Pa·s = 1000 cP (centipoise)
- 1 cP = 0.001 Pa·s

**Viscosity (Kinematic)**:
- 1 m²/s = 10⁶ cSt (centistokes)
- 1 cSt = 10⁻⁶ m²/s

**Temperature**:
- K = °C + 273.15
- °F = 9/5 × °C + 32

---

## Error Handling

Common errors and solutions:

### ValueError: Reynolds number out of range
- **Cause**: Re < 0 or Re = infinity
- **Solution**: Check velocity, diameter, and viscosity values

### ValueError: Relative roughness invalid
- **Cause**: ε/D < 0 or ε/D > 0.5
- **Solution**: Verify roughness and diameter values

### RuntimeWarning: Friction factor iteration did not converge
- **Cause**: Extreme Re or ε/D values
- **Solution**: Use explicit correlation (Moody, Swamee-Jain)

### ValueError: NPSH available is negative
- **Cause**: Vapor pressure exceeds inlet pressure
- **Solution**: Increase inlet pressure or decrease temperature

---

## Performance Optimization

Tips for efficient calculations:

1. **Vectorization**: Use NumPy arrays for multiple calculations
2. **Explicit Correlations**: Use Swamee-Jain for large datasets
3. **Caching**: Store frequently used properties
4. **Profiling**: Use cProfile for bottleneck identification

Example:
```python
import numpy as np

# Efficient: Vectorized calculation
Re_array = np.array([1e4, 1e5, 1e6])
f_array = np.vectorize(friction_factor)(Re_array, 0.0001)

# Inefficient: Loop
f_list = [friction_factor(Re, 0.0001) for Re in Re_array]
```

---

## Contributing

To contribute to the fluids library:

1. Fork the repository: https://github.com/CalebBell/fluids
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

**Testing Requirements**:
- Unit tests for all functions
- Validation against published data
- Documentation with examples

---

## Citation

If using fluids in academic work:

```bibtex
@software{bell2024fluids,
  author = {Bell, Caleb},
  title = {fluids: Fluid dynamics component of Chemical Engineering Design Library},
  year = {2024},
  url = {https://github.com/CalebBell/fluids},
}
```

---

## Version History

- **v1.0**: Initial release with core friction and pump functions
- **v1.1**: Added compressible flow module
- **v1.2**: Enhanced fitting loss coefficients
- **v1.3**: Performance optimizations and NumPy vectorization

Check release notes: https://github.com/CalebBell/fluids/releases

---

## License

The fluids library is distributed under the MIT License.

See: https://github.com/CalebBell/fluids/blob/master/LICENSE.txt

---

## Support

For questions and support:

- **Documentation**: https://fluids.readthedocs.io/
- **Issues**: https://github.com/CalebBell/fluids/issues
- **Discussions**: https://github.com/CalebBell/fluids/discussions
- **Stack Overflow**: Tag questions with [fluids] and [python]

---

*Last Updated: 2025-11-07*
