# Thermodynamics Equations Reference

## First Law of Thermodynamics

### General Forms

**Closed System:**
```
ΔU = Q - W
```

Where:
- ΔU = change in internal energy (J)
- Q = heat transfer to system (J)
- W = work done by system (J)

**Differential form:**
```
dU = δQ - δW
```

**Open System (Control Volume) - General:**
```
dE_CV/dt = Q̇ - Ẇ + Σ ṁ_in·(h + c²/2 + gz)_in - Σ ṁ_out·(h + c²/2 + gz)_out
```

Where:
- E_CV = total energy in control volume (J)
- Q̇ = heat transfer rate (W)
- Ẇ = power (W)
- ṁ = mass flow rate (kg/s)
- h = specific enthalpy (J/kg)
- c = velocity (m/s)
- g = 9.81 m/s²
- z = elevation (m)

**Steady Flow Energy Equation (SFEE):**

For single inlet and outlet:
```
Q̇ - Ẇ = ṁ[(h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)]
```

**Per Unit Mass:**
```
q - w = (h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)
```

Where:
- q = Q̇/ṁ = specific heat transfer (J/kg)
- w = Ẇ/ṁ = specific work (J/kg)

### Pump-Specific Forms

**Pump Work (Shaft Work Input):**

```
Ẇshaft = ṁ[(h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)] + Q̇loss
```

**With Heat Loss:**
```
ẇshaft = (h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁) + q_loss
```

**Adiabatic Pump (Q̇ = 0):**
```
ẇshaft = (h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)
```

**Negligible Kinetic and Potential Energy:**
```
ẇshaft = h₂ - h₁ + q_loss
```

### Incompressible Liquid Forms

**Enthalpy Definition:**
```
h = u + P·v
```

For incompressible liquid (dv ≈ 0):
```
dh = du + v·dP = c_p·dT + v·dP
```

**Specific Work:**
```
w = (P₂ - P₁)/ρ + (c₂² - c₁²)/2 + g(z₂ - z₁)
```

**In Terms of Head:**

Head (m) = energy per unit weight:
```
H = w/g = (P₂ - P₁)/(ρg) + (c₂² - c₁²)/(2g) + (z₂ - z₁)
```

Components:
- **Pressure head:** H_p = (P₂ - P₁)/(ρg)
- **Velocity head:** H_v = (c₂² - c₁²)/(2g)
- **Elevation head:** H_z = z₂ - z₁

**With Losses:**
```
H_total = H_p + H_v + H_z + H_loss
```

**Power Equation:**
```
Ẇshaft = ṁ·g·H/η = ρ·g·Q·H/η
```

Where:
- Q = volumetric flow rate (m³/s)
- η = pump efficiency (dimensionless)

### Temperature Rise

**From Energy Balance:**

For adiabatic pump with incompressible liquid:
```
ΔT = (h₂ - h₁)/c_p - v·(P₂ - P₁)/c_p
```

**From Inefficiency:**

The mechanical energy not converted to fluid energy becomes heat:
```
ΔT = (1 - η)·w/c_p = (1 - η)·g·H/c_p
```

**With Heat Loss:**

If heat is lost to surroundings:
```
ΔT = (1 - η)·g·H/c_p - q_loss/c_p
```

**Dimensional Analysis:**

For water (c_p ≈ 4186 J/kg-K, g = 9.81 m/s²):
```
ΔT [°C] ≈ 0.00234·(1 - η)·H [m]
```

**Example:** H = 100 m, η = 0.80
```
ΔT = 0.00234·0.20·100 = 0.47°C
```

## Second Law of Thermodynamics

### Entropy Equations

**General Definition:**
```
dS ≥ δQ/T
```

Equality holds for reversible processes.

**Entropy Change of Pure Substance:**

**General:**
```
ds = δq_rev/T
```

**For Ideal Gas:**
```
s₂ - s₁ = c_p·ln(T₂/T₁) - R·ln(P₂/P₁)
```

Or:
```
s₂ - s₁ = c_v·ln(T₂/T₁) + R·ln(v₂/v₁)
```

**For Incompressible Liquid:**
```
s₂ - s₁ = c·ln(T₂/T₁)
```

Pressure effect usually negligible.

**For Phase Change (Constant T, P):**
```
Δs_fg = h_fg/T_sat
```

Where h_fg = latent heat of vaporization.

### Entropy Balance

**Control Volume:**
```
dS_CV/dt = Σ Q̇_k/T_k + Σ ṁ_in·s_in - Σ ṁ_out·s_out + Ṡ_gen
```

Where:
- S_CV = entropy in control volume (J/K)
- Q̇_k = heat transfer at boundary k (W)
- T_k = absolute temperature at boundary k (K)
- Ṡ_gen = entropy generation rate (W/K ≥ 0)

**Steady Flow:**
```
Ṡ_gen = Σ ṁ_out·s_out - Σ ṁ_in·s_in - Σ Q̇_k/T_k
```

**Single Stream:**
```
Ṡ_gen = ṁ(s₂ - s₁) - Q̇/T_boundary
```

**Adiabatic Process:**
```
Ṡ_gen = ṁ(s₂ - s₁) ≥ 0
```

For reversible adiabatic (isentropic): s₂ = s₁

### Irreversibility and Availability

**Irreversibility (Lost Work):**
```
İ = T₀·Ṡ_gen
```

Where T₀ = ambient temperature (K)

**Availability (Exergy) Function:**

**Non-flow (Closed System):**
```
φ = (u - u₀) + P₀(v - v₀) - T₀(s - s₀)
```

**Flow (Open System):**
```
ψ = (h - h₀) - T₀(s - s₀) + c²/2 + gz
```

Subscript 0 denotes dead state (ambient conditions).

**Exergy Balance:**
```
Ẋ_in - Ẋ_out - Ẋ_destroyed = dX_CV/dt
```

**Exergy Destruction:**
```
Ẋ_destroyed = T₀·Ṡ_gen
```

**Second Law Efficiency:**
```
η_II = Ẋ_out/Ẋ_in = 1 - Ẋ_destroyed/Ẋ_in
```

For pump:
```
η_II = (ψ₂ - ψ₁)/ẇshaft
```

## Efficiency Definitions

### Pump Efficiencies

**Overall Efficiency:**
```
η = η_h·η_vol·η_mech
```

**Hydraulic Efficiency:**
```
η_h = Water horsepower / (Shaft power + Leakage losses + Friction losses)
```

Or:
```
η_h = H_actual / H_theoretical
```

Typical range: 0.85-0.95

**Volumetric Efficiency:**
```
η_vol = Q_delivered / (Q_delivered + Q_leakage)
```

Typical range: 0.96-0.99

**Mechanical Efficiency:**
```
η_mech = (Power to fluid) / (Shaft power)
```

Accounts for bearing and seal friction.
Typical range: 0.95-0.98

**Isentropic Efficiency:**

For compressible fluids:
```
η_s = (h₂s - h₁) / (h₂ - h₁)
```

Where h₂s is enthalpy at P₂ with s₂s = s₁ (isentropic compression).

**Polytopic Efficiency:**

For multistage or continuously varying compression:
```
η_p = [(k-1)/k]·[ln(P₂/P₁)/ln(T₂/T₁)]
```

Where k = c_p/c_v (specific heat ratio)

### Cycle Efficiencies

**Thermal Efficiency:**
```
η_th = W_net/Q_in = (Q_in - Q_out)/Q_in
```

**Carnot Efficiency (Maximum Possible):**
```
η_Carnot = 1 - T_cold/T_hot
```

**Rankine Cycle Efficiency:**
```
η_Rankine = (w_turbine - w_pump)/q_in
```

Where:
- w_turbine = h₁ - h₂ (turbine work output)
- w_pump = h₄ - h₃ (pump work input)
- q_in = h₁ - h₄ (heat input in boiler)

**Coefficient of Performance (COP):**

For heat pumps and refrigeration:

**Heating:**
```
COP_heating = Q_hot/W_input
```

**Cooling:**
```
COP_cooling = Q_cold/W_input
```

**Carnot COP:**
```
COP_Carnot,heating = T_hot/(T_hot - T_cold)
COP_Carnot,cooling = T_cold/(T_hot - T_cold)
```

## Property Relations

### Fundamental Property Relations

**Maxwell Relations:**

From dU = TdS - PdV:
```
(∂T/∂V)_S = -(∂P/∂S)_V
```

From dH = TdS + VdP:
```
(∂T/∂P)_S = (∂V/∂S)_P
```

From dA = -SdT - PdV:
```
(∂S/∂V)_T = (∂P/∂T)_V
```

From dG = -SdT + VdP:
```
(∂S/∂P)_T = -(∂V/∂T)_P
```

### Specific Heat Relations

**Definitions:**
```
c_p = (∂h/∂T)_P  [Constant pressure]
c_v = (∂u/∂T)_V  [Constant volume]
```

**Relation:**
```
c_p - c_v = T·(∂P/∂T)_V·(∂V/∂T)_P
```

**For Ideal Gas:**
```
c_p - c_v = R
k = c_p/c_v
```

**For Incompressible Liquid:**
```
c_p ≈ c_v ≈ c
```

### Compressibility and Expansivity

**Isothermal Compressibility:**
```
β_T = -(1/V)·(∂V/∂P)_T
```

**Isentropic Compressibility:**
```
β_S = -(1/V)·(∂V/∂P)_S
```

**Volume Expansivity:**
```
α = (1/V)·(∂V/∂T)_P
```

**Relations:**
```
β_T/β_S = c_p/c_v = k
β_T = α²·T/(ρ·c_p)
```

### Joule-Thomson Coefficient

**Definition:**
```
μ_JT = (∂T/∂P)_h
```

Represents temperature change during throttling (constant enthalpy expansion).

**For Ideal Gas:** μ_JT = 0

**For Real Fluids:** μ_JT ≠ 0
- Positive: Cooling upon expansion (most gases at room temp)
- Negative: Heating upon expansion (H₂, He at room temp)

**Relation to Properties:**
```
μ_JT = (1/c_p)·[T·(∂V/∂T)_P - V]
```

## Cycle Equations

### Rankine Cycle

**State Points:**
1. Turbine inlet (high P, high T, superheated)
2. Turbine outlet (low P, two-phase)
3. Condenser outlet / Pump inlet (low P, saturated liquid)
4. Pump outlet / Boiler inlet (high P, compressed liquid)

**Work Terms:**

**Turbine work (per kg):**
```
w_turbine = h₁ - h₂
```

**Pump work (per kg):**
```
w_pump,s = h₄s - h₃ ≈ v₃·(P₄ - P₃)
w_pump = w_pump,s/η_pump
```

**Heat Transfer:**

**Boiler:**
```
q_in = h₁ - h₄
```

**Condenser:**
```
q_out = h₂ - h₃
```

**Cycle Efficiency:**
```
η = (w_turbine - w_pump)/q_in
```

**Back Work Ratio:**
```
BWR = w_pump/w_turbine
```

Typically 0.01-0.03 for Rankine cycles.

### Brayton Cycle (Gas Turbine)

**State Points:**
1. Compressor inlet (ambient)
2. Compressor outlet (high P)
3. Combustor outlet / Turbine inlet (high P, high T)
4. Turbine outlet (ambient P)

**Work Terms:**

**Compressor work:**
```
w_comp = h₂ - h₁ = c_p·(T₂ - T₁)
```

**Turbine work:**
```
w_turb = h₃ - h₄ = c_p·(T₃ - T₄)
```

**Heat Transfer:**

**Combustor:**
```
q_in = h₃ - h₂ = c_p·(T₃ - T₂)
```

**Cycle Efficiency:**
```
η = (w_turb - w_comp)/q_in = 1 - T₁/T₂ = 1 - (P₁/P₂)^((k-1)/k)
```

For ideal cycle with constant properties.

### Refrigeration Cycle (Vapor Compression)

**State Points:**
1. Compressor inlet (low P, saturated vapor)
2. Compressor outlet (high P, superheated vapor)
3. Condenser outlet (high P, saturated liquid)
4. Evaporator inlet (low P, two-phase)

**Work and Heat:**

**Compressor work:**
```
w_comp = h₂ - h₁
```

**Condenser heat rejection:**
```
q_out = h₂ - h₃
```

**Evaporator heat absorption:**
```
q_in = h₁ - h₄
```

**COP:**
```
COP = q_in/w_comp = (h₁ - h₄)/(h₂ - h₁)
```

**For Carnot Refrigeration Cycle:**
```
COP_Carnot = T_L/(T_H - T_L)
```

## Heat Transfer Equations

### Conduction

**Fourier's Law:**
```
q̇" = -k·dT/dx
```

Where:
- q̇" = heat flux (W/m²)
- k = thermal conductivity (W/m-K)
- dT/dx = temperature gradient (K/m)

**One-Dimensional Steady Conduction:**

**Plane wall:**
```
Q̇ = k·A·(T₁ - T₂)/L
```

**Cylindrical wall (pipe):**
```
Q̇ = 2π·k·L·(T_i - T_o)/ln(r_o/r_i)
```

**Thermal Resistance:**
```
R_th = ΔT/Q̇
```

For plane wall: R_th = L/(k·A)

### Convection

**Newton's Law of Cooling:**
```
Q̇ = h·A·(T_surface - T_fluid)
```

Where h = convection heat transfer coefficient (W/m²-K)

**Typical h Values:**
- Free convection, gas: 5-25 W/m²-K
- Free convection, liquid: 50-1000 W/m²-K
- Forced convection, gas: 25-250 W/m²-K
- Forced convection, liquid: 100-20,000 W/m²-K

**Dimensionless Numbers:**

**Nusselt Number:**
```
Nu = h·L/k
```

**Reynolds Number:**
```
Re = ρ·V·L/μ = V·L/ν
```

**Prandtl Number:**
```
Pr = c_p·μ/k = ν/α
```

**Grashof Number (Free Convection):**
```
Gr = g·β·ΔT·L³/ν²
```

**Correlations:**

**Forced convection in pipes (turbulent, Re > 2300):**
```
Nu = 0.023·Re^0.8·Pr^0.4
```

**Free convection, vertical plate:**
```
Nu = 0.59·(Gr·Pr)^0.25  [10⁴ < Gr·Pr < 10⁹]
```

### Radiation

**Stefan-Boltzmann Law:**
```
Q̇ = ε·σ·A·T⁴
```

Where:
- ε = emissivity (0 < ε < 1)
- σ = 5.67×10⁻⁸ W/m²-K⁴ (Stefan-Boltzmann constant)
- T = absolute temperature (K)

**Net Radiation Between Two Surfaces:**
```
Q̇₁₂ = σ·A·F₁₂·(T₁⁴ - T₂⁴)
```

Where F₁₂ is a view factor and emissivity function.

**Linearized Form (Small ΔT):**
```
Q̇ = h_r·A·(T₁ - T₂)
```

Where:
```
h_r = 4·ε·σ·T_avg³
```

### Overall Heat Transfer

**Combined Resistances:**
```
Q̇ = U·A·ΔT_overall
```

Where U = overall heat transfer coefficient:
```
1/U = 1/h_i + t_wall/k_wall + 1/h_o
```

**With Fouling:**
```
1/U = 1/h_i + R_fi + t_wall/k_wall + R_fo + 1/h_o
```

Where R_f = fouling resistance (m²-K/W)

## Fluid Properties

### Water and Steam (IAPWS-IF97)

Use CoolProp or steam tables for accurate properties.

**Approximate Values at 20°C, 1 atm:**
- Density: ρ = 998.2 kg/m³
- Specific heat: c_p = 4182 J/kg-K
- Viscosity: μ = 1.002×10⁻³ Pa·s
- Thermal conductivity: k = 0.598 W/m-K

**Saturation Properties (Example: 100°C, 1.013 bar):**
- h_f = 419.04 kJ/kg (saturated liquid)
- h_fg = 2257.0 kJ/kg (latent heat)
- h_g = 2676.0 kJ/kg (saturated vapor)
- s_f = 1.3069 kJ/kg-K
- s_fg = 6.048 kJ/kg-K
- s_g = 7.3549 kJ/kg-K

### Ideal Gas Properties

**Equation of State:**
```
P·v = R·T
or
P·V = n·R_u·T
```

Where:
- R = specific gas constant (J/kg-K)
- R_u = universal gas constant = 8.314 J/mol-K
- R = R_u/M (M = molecular weight)

**Common Gases:**

**Air:**
- M = 28.97 kg/kmol
- R = 287 J/kg-K
- c_p = 1005 J/kg-K (at 300K)
- c_v = 718 J/kg-K
- k = c_p/c_v = 1.4

**Nitrogen (N₂):**
- M = 28.01 kg/kmol
- R = 297 J/kg-K
- k = 1.4

**Carbon Dioxide (CO₂):**
- M = 44.01 kg/kmol
- R = 189 J/kg-K
- k = 1.3

### Real Fluids

**Compressibility Factor:**
```
Z = P·v/(R·T)
```

For ideal gas: Z = 1
For real gas: Z ≠ 1

**Reduced Properties:**
```
T_r = T/T_critical
P_r = P/P_critical
```

**Generalized Compressibility Charts:**

Z = f(T_r, P_r) - available in thermodynamics textbooks

**Van der Waals Equation:**
```
(P + a/v²)·(v - b) = R·T
```

Constants a, b depend on critical properties.

## Dimensional Analysis

### Common Unit Conversions

**Pressure:**
- 1 bar = 10⁵ Pa = 100 kPa
- 1 atm = 101.325 kPa = 1.01325 bar
- 1 psi = 6.895 kPa
- 1 m H₂O = 9.81 kPa

**Energy/Work:**
- 1 kJ = 1000 J
- 1 kW·h = 3600 kJ
- 1 Btu = 1.055 kJ
- 1 cal = 4.184 J

**Power:**
- 1 kW = 1000 W
- 1 hp = 745.7 W = 0.7457 kW
- 1 Btu/h = 0.293 W

**Temperature:**
- K = °C + 273.15
- °F = 1.8·°C + 32
- °R = °F + 459.67

**Specific Energy:**
- 1 kJ/kg = 1000 J/kg
- 1 Btu/lbm = 2.326 kJ/kg

### Dimensional Consistency

Always check dimensional consistency of equations.

**Example: Power**
```
Ẇ [W] = ṁ [kg/s] · g [m/s²] · H [m] / η [-]
     = kg/s · m/s² · m
     = kg·m²/s³
     = J/s
     = W  ✓
```

**Example: Temperature Rise**
```
ΔT [K] = g [m/s²] · H [m] / c_p [J/kg-K]
       = m/s² · m / (J/kg-K)
       = m²/s² / (J/kg-K)
       = (kg·m²/s²)/kg / (J/kg-K)
       = J/kg / (J/kg-K)
       = K  ✓
```

## Summary Tables

### Key Constants

| Constant | Symbol | Value | Units |
|----------|--------|-------|-------|
| Gravitational acceleration | g | 9.81 | m/s² |
| Universal gas constant | R_u | 8.314 | J/mol-K |
| Stefan-Boltzmann constant | σ | 5.67×10⁻⁸ | W/m²-K⁴ |
| Atmospheric pressure | P_atm | 101.325 | kPa |

### Water Properties (20°C, 1 atm)

| Property | Symbol | Value | Units |
|----------|--------|-------|-------|
| Density | ρ | 998.2 | kg/m³ |
| Specific heat | c_p | 4182 | J/kg-K |
| Dynamic viscosity | μ | 1.002×10⁻³ | Pa·s |
| Kinematic viscosity | ν | 1.004×10⁻⁶ | m²/s |
| Thermal conductivity | k | 0.598 | W/m-K |
| Vapor pressure | P_v | 2.339 | kPa |

### Air Properties (20°C, 1 atm)

| Property | Symbol | Value | Units |
|----------|--------|-------|-------|
| Density | ρ | 1.204 | kg/m³ |
| Specific heat (const P) | c_p | 1005 | J/kg-K |
| Specific heat (const V) | c_v | 718 | J/kg-K |
| Gas constant | R | 287 | J/kg-K |
| Specific heat ratio | k | 1.4 | - |
| Dynamic viscosity | μ | 1.825×10⁻⁵ | Pa·s |
| Thermal conductivity | k | 0.0257 | W/m-K |

## References

### Textbooks

1. **Moran, M.J., Shapiro, H.N., Boettner, D.D., and Bailey, M.B.**
   *"Fundamentals of Engineering Thermodynamics"* 9th Edition, Wiley (2018)
   - Comprehensive thermodynamics textbook
   - Chapters 4-6: Control volume analysis
   - Chapter 8: Vapor power systems

2. **Cengel, Y.A. and Boles, M.A.**
   *"Thermodynamics: An Engineering Approach"* 9th Edition, McGraw-Hill (2019)
   - Practical engineering focus
   - Chapter 5: Energy analysis of closed systems
   - Chapter 6: Energy analysis of control volumes
   - Chapter 7: Entropy

3. **Bejan, A.**
   *"Advanced Engineering Thermodynamics"* 4th Edition, Wiley (2016)
   - Advanced topics including exergy
   - Optimization and second law analysis
   - Chapter 5: Irreversibility and availability

### Standards

1. **ASME PTC 19.5 (2004)**
   *"Flow Measurement"*
   - Thermodynamic property calculations
   - Uncertainty analysis

2. **ISO 5198 (1987)**
   *"Centrifugal, Mixed Flow and Axial Pumps - Code for Hydraulic Performance Tests"*
   - Energy balance procedures
   - Efficiency calculations

3. **API 610 (2010)**
   *"Centrifugal Pumps for Petroleum, Petrochemical and Natural Gas Industries"*
   - Temperature rise calculations
   - Minimum flow requirements

### Property Databases

1. **NIST Chemistry WebBook**
   - Free online thermodynamic data
   - https://webbook.nist.gov/chemistry/fluid/

2. **CoolProp**
   - Open-source thermodynamic property library
   - http://www.coolprop.org/
   - Python, MATLAB, Excel interfaces

3. **REFPROP (NIST)**
   - Reference fluid thermodynamic properties
   - Most accurate property database (commercial)
