---
name: thermodynamics-workflow
description: "Analyze thermodynamic cycles and heat transfer in pump systems"
category: thinking
domain: thermodynamics
complexity: intermediate
dependencies:
  - CoolProp
  - numpy
---

# Thermodynamics Analysis Workflow

## Overview

This skill provides a systematic approach to thermodynamic analysis of pump systems, including energy balance, entropy analysis, heat transfer, and cycle calculations. The methodology follows classical thermodynamics principles applied specifically to pumping systems and fluid machinery.

## Thermodynamic Analysis Workflow

### 1. System Definition and Boundaries

**Define the Control Volume:**

The first step in any thermodynamic analysis is to clearly define the system boundaries.

**For Pump Systems:**
- **Control Volume:** The pump casing from suction flange to discharge flange
- **Boundaries:** Physical surfaces where fluid/energy crosses
- **Interfaces:**
  - Inlet: Suction flange (station 1)
  - Outlet: Discharge flange (station 2)
  - Shaft: Mechanical energy input
  - Walls: Heat transfer surface (if applicable)

**System Classification:**
- **Open System (Control Volume):** Most pump analyses
- **Steady Flow:** Constant properties with time
- **Transient Flow:** Start-up, shutdown, surge conditions

**Assumptions to Document:**
- [ ] Steady-state operation
- [ ] Uniform properties at inlet/outlet
- [ ] Negligible kinetic energy changes (if applicable)
- [ ] Adiabatic walls (or specify heat transfer)
- [ ] No chemical reactions
- [ ] Incompressible flow (for liquids) or compressible (for gases)

**Coordinate System:**
- Define datum for potential energy (z = 0)
- Specify flow direction (positive flow)
- Identify measurement planes

### 2. State Point Identification

**Thermodynamic State Definition:**

For a pure substance, two independent intensive properties define the state.

**State Point 1 (Inlet/Suction):**

Identify and measure:
```
Properties at inlet:
- Pressure: P₁ (Pa, bar, psi)
- Temperature: T₁ (K, °C)
- Velocity: c₁ (m/s)
- Elevation: z₁ (m)
- Phase: Liquid, vapor, two-phase
```

**State Point 2 (Outlet/Discharge):**

Identify and measure:
```
Properties at outlet:
- Pressure: P₂ (Pa, bar, psi)
- Temperature: T₂ (K, °C)
- Velocity: c₂ (m/s)
- Elevation: z₂ (m)
- Phase: Liquid, vapor, two-phase
```

**Additional State Points:**

For complex systems:
- **Interstage conditions** (multistage pumps)
- **Impeller exit** (before volute losses)
- **Intermediate heat exchangers**
- **Recirculation loops**

**Property Diagrams:**

Plot state points on:
- **P-h diagram:** Pressure vs. specific enthalpy
- **T-s diagram:** Temperature vs. specific entropy
- **P-v diagram:** Pressure vs. specific volume
- **Mollier diagram:** h-s chart

### 3. Property Evaluation

**Fluid Property Sources:**

**For Pure Substances:**
- **Water/Steam:** IAPWS-IF97 formulation (via CoolProp)
- **Refrigerants:** NIST REFPROP database
- **Hydrocarbons:** Peng-Robinson or Soave-Redlich-Kwong EOS
- **Cryogenic fluids:** NIST correlations

**For Incompressible Liquids:**
```
Specific volume: v ≈ constant
Specific heat: cp ≈ constant
Internal energy: du = cp·dT
Enthalpy: dh = cp·dT + v·dP
```

**For Ideal Gases:**
```
Equation of state: P·v = R·T
Specific heats: cp, cv = f(T)
Internal energy: du = cv·dT
Enthalpy: dh = cp·dT
Entropy: ds = cp·dT/T - R·dP/P
```

**Property Evaluation Methods:**

**Method 1: Using CoolProp (Recommended)**
```python
from CoolProp.CoolProp import PropsSI

# Water at 25°C, 2 bar
h = PropsSI('H', 'T', 298.15, 'P', 200000, 'Water')  # J/kg
s = PropsSI('S', 'T', 298.15, 'P', 200000, 'Water')  # J/kg-K
rho = PropsSI('D', 'T', 298.15, 'P', 200000, 'Water')  # kg/m³
```

**Method 2: Property Tables**
- Steam tables (superheated, saturated, compressed liquid)
- Interpolation between table values
- Linear or quadratic interpolation

**Method 3: Empirical Correlations**
- Density: ρ(T, P) correlations
- Viscosity: Andrade or Vogel equation
- Specific heat: polynomial fits

**Critical Properties to Evaluate:**

At each state point:
- **Specific enthalpy:** h (J/kg or kJ/kg)
- **Specific entropy:** s (J/kg-K or kJ/kg-K)
- **Specific volume:** v (m³/kg) or density ρ (kg/m³)
- **Specific internal energy:** u (J/kg or kJ/kg)
- **Quality:** x (for two-phase mixtures)
- **Viscosity:** μ (Pa·s) - for loss calculations
- **Thermal conductivity:** k (W/m-K) - for heat transfer

### 4. First Law Analysis (Energy Balance)

**General Energy Equation for Open Systems:**

For a control volume in steady flow:

```
Q̇ - Ẇ = ṁ[(h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)]
```

Where:
- Q̇ = heat transfer rate (W)
- Ẇ = work transfer rate (W)
- ṁ = mass flow rate (kg/s)
- h = specific enthalpy (J/kg)
- c = velocity (m/s)
- g = 9.81 m/s²
- z = elevation (m)

**For Pump Systems:**

Work is done ON the fluid (Ẇ < 0 by convention), so:

```
Ẇshaft = ṁ[(h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)] + Q̇loss
```

**Simplified Forms:**

**Incompressible Liquid (Constant Density):**

```
Ẇshaft/ṁ = (P₂ - P₁)/ρ + (c₂² - c₁²)/2 + g(z₂ - z₁) + q_loss
```

This can be written in terms of head:

```
H = (P₂ - P₁)/(ρg) + (c₂² - c₁²)/(2g) + (z₂ - z₁) + h_loss
```

**Adiabatic Pump (No Heat Loss):**

```
Ẇshaft = ṁ[(h₂ - h₁) + (c₂² - c₁²)/2 + g(z₂ - z₁)]
```

**Negligible Kinetic and Potential Energy Changes:**

```
Ẇshaft = ṁ(h₂ - h₁) + Q̇loss
```

**Energy Balance Components:**

1. **Enthalpy Change:**
   - Pressure rise: Δh_pressure = ∫v·dP
   - Temperature rise: Δh_thermal = cp·ΔT

2. **Kinetic Energy Change:**
   - Usually small for pumps (< 5% of total)
   - Important for high-velocity applications

3. **Potential Energy Change:**
   - Significant for vertical pumps
   - Negligible for horizontal installations

4. **Heat Transfer:**
   - Heat loss to environment: Q̇loss
   - Can be from temperature measurements
   - Typically 1-5% of shaft power

**Temperature Rise in Pumped Fluid:**

The temperature rise due to inefficiency:

```
ΔT = (1 - η)·Ẇshaft / (ṁ·cp)
```

Where η is the pump efficiency.

For incompressible liquids:

```
ΔT = (1 - η)·g·H / cp
```

**Example:** Water pump, H = 100 m, η = 0.80, cp = 4186 J/kg-K
```
ΔT = (1 - 0.80)·9.81·100 / 4186 = 0.47°C
```

### 5. Second Law Analysis (Entropy and Efficiency)

**Entropy Balance:**

For a control volume in steady flow:

```
ṁ(s₂ - s₁) = Q̇/T_boundary + Ṡgen
```

Where:
- s = specific entropy (J/kg-K)
- Ṡgen = entropy generation rate (W/K)
- T_boundary = boundary temperature (K)

**Entropy Generation:**

Measures irreversibility in the process:

```
Ṡgen = ṁ(s₂ - s₁) - Q̇/T_boundary
```

For adiabatic pump:

```
Ṡgen = ṁ(s₂ - s₁)
```

**Isentropic Efficiency:**

Compares actual process to ideal isentropic process:

```
η_isentropic = (h₂s - h₁) / (h₂ - h₁)
```

Where:
- h₂s = specific enthalpy at P₂ with s₂s = s₁ (isentropic)
- h₂ = actual specific enthalpy at outlet

**For Incompressible Liquids:**

The isentropic work is:

```
w_isentropic = ∫v·dP = v·(P₂ - P₁) = (P₂ - P₁)/ρ
```

**Exergy Analysis:**

Exergy (availability) represents maximum useful work:

```
ψ = (h - h₀) - T₀(s - s₀) + c²/2 + gz
```

Where subscript 0 denotes dead state (ambient conditions).

**Exergy destruction:**

```
Ẋdestroyed = T₀·Ṡgen
```

**Second Law Efficiency:**

```
η_II = Ẋuseful / Ẋinput = 1 - Ẋdestroyed / Ẋinput
```

This represents how effectively the pump uses available energy.

**Loss Breakdown by Entropy Generation:**

Different loss mechanisms contribute to entropy generation:

1. **Friction losses:** Fluid friction in impeller and volute
2. **Shock losses:** Flow separation and incidence
3. **Leakage losses:** Recirculation through clearances
4. **Disk friction:** Viscous shear on impeller surfaces
5. **Heat transfer:** Irreversible heat transfer to/from fluid

### 6. Cycle Calculations

**Thermodynamic Cycles Relevant to Pumps:**

While pumps are typically single components, they operate within larger thermodynamic cycles.

#### Rankine Cycle (Steam Power Plants)

**Components:**
1. Boiler (heat addition)
2. Turbine (power output)
3. Condenser (heat rejection)
4. **Feedwater pump** (pressure increase)

**Pump Analysis in Rankine Cycle:**

State points:
- State 3: Condenser outlet (saturated or subcooled liquid)
- State 4: Boiler inlet (high-pressure compressed liquid)

**Pump work (isentropic):**
```
w_pump,s = v₃·(P₄ - P₃) = h₄s - h₃
```

**Actual pump work:**
```
w_pump = w_pump,s / η_pump
```

**Impact on Cycle Efficiency:**
```
η_cycle = (w_turbine - w_pump) / q_in
```

Pump work is typically 1-3% of turbine work in Rankine cycles.

#### Brayton Cycle (Gas Turbines)

**Not directly applicable to liquid pumps, but relevant for:**
- Fuel pumps in gas turbines
- Liquid rocket propellant pumps
- High-pressure fluid systems

**Key Differences from Rankine:**
- Working fluid remains gas throughout
- Compression instead of pumping
- Higher temperature ratios

#### Refrigeration Cycle

**Relevant for:**
- Refrigerant pumps in heat pumps
- Chiller systems
- Cryogenic pumping

**Components:**
1. Evaporator
2. Compressor (or pump for liquid systems)
3. Condenser
4. Expansion valve

**Coefficient of Performance (COP):**
```
COP = Q_useful / W_input
```

Where pump work contributes to W_input.

## Applications to Pumps

### Temperature Rise in Pumped Fluid

**Physical Mechanism:**

Inefficiencies in pumps convert mechanical energy to thermal energy, causing temperature rise.

**Calculation Method:**

**From Energy Balance:**
```
ΔT = (h₂ - h₁) / cp - v·(P₂ - P₁) / cp
```

**Simplified for Liquids:**
```
ΔT = (1 - η)·v·(P₂ - P₁) / cp
```

Or in terms of head:
```
ΔT = (1 - η)·g·H / cp
```

**Factors Affecting Temperature Rise:**

1. **Pump efficiency:** Lower η → higher ΔT
2. **Pressure rise:** Higher ΔP → higher ΔT
3. **Specific heat:** Lower cp → higher ΔT
4. **Operating point:** Off-design → lower η → higher ΔT

**Practical Implications:**

- **Cold fluids:** May cause condensation in surroundings
- **Hot fluids:** May approach boiling point, causing cavitation
- **Viscous fluids:** Significant heating (η low, friction high)
- **Volatile fluids:** Risk of vaporization
- **Polymers:** May degrade at elevated temperatures

**Example Calculations:**

**Water pump:**
- H = 200 m, η = 0.75, cp = 4186 J/kg-K
- ΔT = 0.25·9.81·200/4186 = 1.17°C

**Hydraulic oil:**
- ΔP = 200 bar, η = 0.70, ρ = 850 kg/m³, cp = 2000 J/kg-K
- ΔT = 0.30·200e5/(850·2000) = 3.53°C

**High viscosity fluid:**
- H = 50 m, η = 0.50 (low due to viscosity), cp = 2500 J/kg-K
- ΔT = 0.50·9.81·50/2500 = 0.98°C

### Efficiency and Losses

**Overall Efficiency:**

```
η = η_h · η_vol · η_mech
```

**Component Efficiencies:**

1. **Hydraulic Efficiency (η_h):**
   - Accounts for fluid friction and shock losses
   - Typical: 0.85-0.95

2. **Volumetric Efficiency (η_vol):**
   - Accounts for internal leakage
   - Typical: 0.96-0.99

3. **Mechanical Efficiency (η_mech):**
   - Accounts for bearing and seal friction
   - Typical: 0.95-0.98

**Thermodynamic Loss Analysis:**

Each loss mechanism increases entropy:

**Friction Loss:**
```
Ṡ_friction = ṁ·g·h_friction / T_avg
```

**Leakage Loss:**
```
Ṡ_leakage = ṁ_leak·(s_discharge - s_suction)
```

**Heat Transfer Loss:**
```
Ṡ_heat = Q̇_loss / T_boundary - Q̇_loss / T_fluid
```

**Total Entropy Generation:**
```
Ṡ_total = Ṡ_friction + Ṡ_leakage + Ṡ_heat + Ṡ_mechanical
```

**Exergy Efficiency:**

```
η_exergy = 1 - T₀·Ṡ_total / Ẇshaft
```

This gives a more complete picture of thermodynamic performance.

### Heat Transfer Considerations

**Heat Transfer Modes in Pumps:**

1. **Conduction:**
   - Through pump casing walls
   - Through shaft to bearings
   - Through impeller material

2. **Convection:**
   - Forced convection with pumped fluid
   - Natural convection to ambient air
   - Internal recirculation loops

3. **Radiation:**
   - Usually negligible unless high temperatures
   - Can be significant for hot oil pumps (>150°C)

**Heat Loss Estimation:**

**From Casing Walls:**
```
Q̇_loss = h·A·(T_fluid - T_ambient)
```

Where:
- h = overall heat transfer coefficient (W/m²-K)
- A = surface area (m²)
- Typical h = 5-20 W/m²-K (natural convection, insulated)

**From Energy Balance:**
```
Q̇_loss = Ẇshaft - ṁ·(h₂ - h₁) - ṁ·(c₂² - c₁²)/2 - ṁ·g·(z₂ - z₁)
```

**Impact on Performance:**

- Reduces outlet temperature
- Affects fluid properties (viscosity, density)
- Can cause thermal gradients and distortion
- Important for cryogenic and high-temperature pumps

**Thermal Management:**

**Cooling Requirements:**
- External cooling jackets
- Flush plans for mechanical seals
- Heat exchangers in recirculation loops

**Insulation:**
- Prevents heat loss (hot fluids)
- Prevents condensation (cold fluids)
- Improves safety (hot surfaces)

### Multistage Pump Thermodynamics

**Staging Rationale:**

High head requirements divided among multiple impellers:

```
H_total = H_stage1 + H_stage2 + ... + H_stage_n
```

**Thermodynamic Analysis:**

**Between Stages:**
```
State 1 → Impeller 1 → State 2 → Diffuser 1 → State 3 →
Impeller 2 → State 4 → ... → Final discharge
```

**Temperature Rise per Stage:**
```
ΔT_stage = (1 - η_stage)·g·H_stage / cp
```

**Total Temperature Rise:**
```
ΔT_total = Σ ΔT_stage = (1 - η_overall)·g·H_total / cp
```

**Interstage Cooling:**

For high-head or high-temperature applications:

```
Q̇_cooling = ṁ·cp·(T_beforecooling - T_aftercooling)
```

Benefits:
- Reduces final discharge temperature
- Improves fluid properties (lower viscosity)
- Prevents vaporization
- Increases pump life

**Pressure Compounding Effects:**

As pressure increases through stages:
- Density may change (compressible fluids)
- Leakage flow increases (higher ΔP)
- Seal loads increase
- Structural requirements change

**Energy Distribution:**

Not all stages contribute equally:
- First stage: Highest NPSH requirement
- Middle stages: Maximum efficiency
- Last stage: Highest pressure loading

## Common Cycles Relevant to Pumps

### Rankine Cycle (Detailed)

**Steam Power Plant Configuration:**

```
Boiler → Turbine → Condenser → Feedwater Pump → Boiler
    ↓        ↓         ↓              ↓
   High P   Power   Low P          High P
   High T  Output   Low T        Moderate T
```

**State Point Analysis:**

**State 1:** Turbine exhaust entering condenser
- P₁ = 0.05-0.1 bar (vacuum conditions)
- T₁ = saturation temperature at P₁
- x₁ = 0.85-0.95 (wet steam)

**State 2:** Condenser outlet / Pump inlet
- P₂ = P₁ (pressure drop minimal)
- T₂ = T_sat(P₂) - subcooling (typically 2-5°C)
- x₂ = 0 (saturated or subcooled liquid)

**State 3:** Pump outlet / Boiler inlet
- P₃ = 100-300 bar (high pressure)
- T₃ ≈ T₂ + ΔT_pump (small temperature rise)
- Compressed liquid

**State 4:** Boiler outlet / Turbine inlet
- P₄ = P₃ - ΔP_boiler
- T₄ = 500-600°C (superheated steam)
- Superheated vapor

**Pump Work Calculation:**

**Isentropic work:**
```
w_pump,s = h₃s - h₂ = v₂·(P₃ - P₂)
```

For water at 40°C (v ≈ 0.001008 m³/kg):
```
w_pump,s = 0.001008·(250e5 - 0.08e5) ≈ 2.52 kJ/kg
```

**Actual work:**
```
w_pump = w_pump,s / η_pump = 2.52 / 0.85 ≈ 2.96 kJ/kg
```

**Temperature rise:**
```
ΔT = (h₃ - h₂) / cp ≈ 0.44 / 4.18 ≈ 0.11°C (negligible)
```

**Cycle Efficiency:**

```
η_Rankine = (w_turbine - w_pump) / q_in
```

Typical values:
- Simple cycle: 30-40%
- Reheat cycle: 38-45%
- Regenerative cycle: 40-50%

Pump work is small (~1-2% of turbine work) but critical for operation.

### Refrigeration Cycle

**Vapor Compression Cycle with Liquid Pump:**

Some industrial refrigeration systems use liquid pumps instead of or in addition to compressors.

**State Points:**

1. **Evaporator outlet:** Low P, low T, saturated vapor
2. **Compressor outlet:** High P, high T, superheated vapor
3. **Condenser outlet:** High P, moderate T, saturated liquid
4. **Pump outlet:** Very high P, moderate T, compressed liquid
5. **Expansion valve outlet:** Low P, low T, two-phase mixture

**Pump Application:**

- **Liquid overfeed systems:** Circulate liquid refrigerant
- **CO₂ transcritical systems:** Liquid CO₂ pumping
- **Cascade systems:** Pumping between pressure levels

**COP Calculation:**

```
COP = Q_evaporator / (W_compressor + W_pump)
```

Pump work usually small compared to compressor work.

### Closed-Loop Heat Transfer Cycles

**Applications:**
- Solar thermal systems
- Geothermal heat pumps
- Industrial heat recovery
- District heating/cooling

**Thermodynamic Analysis:**

**Heating Mode:**
```
Q̇_delivered = ṁ·cp·(T_supply - T_return)
W_pump = ṁ·g·H / η_pump
COP_heating = Q̇_delivered / W_pump
```

**Cooling Mode:**
```
Q̇_removed = ṁ·cp·(T_return - T_supply)
COP_cooling = Q̇_removed / W_pump
```

**Optimization:**
- Minimize flow rate (reduces pump power)
- Maximize ΔT (reduces flow for given load)
- Balance pipe friction vs. pump size

## Summary

Thermodynamic analysis of pump systems requires systematic application of:

1. **Control volume definition:** Clear boundaries and assumptions
2. **State point identification:** Complete property evaluation
3. **Property evaluation:** Accurate fluid properties from reliable sources
4. **First law:** Energy balance accounting for all forms of energy
5. **Second law:** Entropy generation and efficiency analysis
6. **Cycle analysis:** Understanding pump role in larger systems

Key outputs:
- Power requirements
- Temperature rise
- Efficiency breakdown
- Loss mechanisms
- Optimization opportunities

## References

See `equations-reference.md` for detailed equations and `examples.md` for worked examples.
