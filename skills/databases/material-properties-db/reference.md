# Material Properties Database - Technical Reference

Comprehensive reference for property correlation equations, data sources, and engineering standards.

## Table of Contents

1. [Viscosity Correlations](#viscosity-correlations)
2. [Density Correlations](#density-correlations)
3. [Vapor Pressure Correlations](#vapor-pressure-correlations)
4. [Transport Properties](#transport-properties)
5. [Dimensionless Numbers](#dimensionless-numbers)
6. [Data Sources and Standards](#data-sources-and-standards)
7. [Standard Conditions Definitions](#standard-conditions-definitions)
8. [Conversion Factors](#conversion-factors)
9. [Verification Data](#verification-data)

---

## Viscosity Correlations

### Dynamic Viscosity (μ)

Dynamic viscosity represents internal friction in a fluid. SI unit: Pa·s (Pascal-second)

**Common units:**
- Pa·s = N·s/m² = kg/(m·s)
- mPa·s = cP (centipoise) - common for liquids
- μPa·s - common for gases
- 1 Pa·s = 1000 mPa·s = 10⁶ μPa·s

### Andrade Equation (Liquids)

Simple exponential model for temperature dependence:

```
μ(T) = A · exp(B/T)
```

**Parameters:**
- μ = dynamic viscosity (Pa·s)
- T = absolute temperature (K)
- A = pre-exponential factor (Pa·s)
- B = activation energy parameter (K)

**Characteristics:**
- Two-parameter model
- Good for limited temperature ranges
- Typical accuracy: ±5-10%
- Fails near critical point

**Example: Water**
- A = 0.001002 Pa·s
- B = 1792 K
- Valid: 0-100°C
- Accuracy: ±5%

**Derivation basis:**
Based on Arrhenius-type activation energy for molecular motion:
```
μ ∝ exp(E_a / RT)
```
where E_a is activation energy for viscous flow.

### Vogel-Fulcher-Tammann (VFT) Equation

More accurate three-parameter model for liquids:

**Natural exponential form:**
```
μ(T) = A · exp(B / (T - C))
```

**Base-10 exponential form (common for water):**
```
μ(T) = A · 10^(B / (T - C))
```

**Parameters:**
- A = pre-exponential factor (Pa·s or mPa·s)
- B = pseudo-activation energy (K)
- C = ideal glass transition temperature (K)

**Characteristics:**
- Three-parameter model (better fit than Andrade)
- C typically 95-140 K for oils
- Excellent for oils and high-viscosity fluids
- Accuracy: ±2-5%

**Example: Water (base-10 form)**
- A = 0.02414 mPa·s
- B = 247.8 K
- C = 140 K
- Valid: 0-100°C
- Accuracy: ±2%

**Physical interpretation:**
- C represents temperature where viscosity would theoretically become infinite
- Related to glass transition in polymers and supercooled liquids

### Walther Equation (ASTM D341)

Industry standard for petroleum products:

```
log₁₀(log₁₀(ν + 0.7)) = A - B · log₁₀(T)
```

**Parameters:**
- ν = kinematic viscosity (cSt = mm²/s)
- T = absolute temperature (K)
- A, B = constants determined from two-point calibration
- 0.7 = empirical constant

**Two-point calibration:**
Given ν₁ at T₁ and ν₂ at T₂:

```
Z₁ = log₁₀(log₁₀(ν₁ + 0.7))
Z₂ = log₁₀(log₁₀(ν₂ + 0.7))

B = (Z₁ - Z₂) / (log₁₀(T₁) - log₁₀(T₂))
A = Z₁ + B · log₁₀(T₁)
```

**Characteristics:**
- ASTM D341 standard method
- Excellent for petroleum products
- Typical accuracy: ±2% over wide range
- Valid: -20°C to 150°C for oils

**Common calibration points:**
- ISO VG oils: 40°C and 100°C
- SAE grades: 100°C and 150°C

**Advantages:**
- Straight line on Walther chart (log-log axes)
- Accurate extrapolation
- Industry-accepted standard

### Sutherland's Law (Gases)

Temperature dependence for ideal gas viscosity:

```
μ(T) = μ₀ · (T/T₀)^(3/2) · (T₀ + S) / (T + S)
```

**Parameters:**
- μ₀ = reference viscosity at T₀ (Pa·s)
- T₀ = reference temperature (usually 273.15 K)
- S = Sutherland constant (K)
- T = absolute temperature (K)

**Sutherland constants for common gases:**

| Gas | μ₀ @ 273.15 K (μPa·s) | S (K) | Valid Range (K) |
|-----|----------------------|-------|-----------------|
| Air | 17.16 | 110.4 | 100-1900 |
| N₂ | 16.66 | 111 | 100-1900 |
| O₂ | 19.20 | 127 | 100-1900 |
| CO₂ | 13.73 | 240 | 200-1900 |
| H₂ | 8.41 | 72 | 100-1900 |
| He | 18.60 | 79.4 | 100-2000 |
| Ar | 21.05 | 142 | 100-2000 |
| CH₄ | 10.23 | 164 | 200-1500 |
| H₂O vapor | 8.85 | 961 | 380-1000 |

**Characteristics:**
- Based on kinetic theory of gases
- Valid for ideal gases (low to moderate pressure)
- Accuracy: ±2% within valid range
- Breaks down near critical point and at high pressure

**Theoretical basis:**
Sutherland improved Chapman-Enskog theory by accounting for intermolecular attractive forces. The S constant represents the effective temperature of attraction.

**Simplified form for small temperature changes:**
For |T - T₀| << T₀:
```
μ(T) ≈ μ₀ · (T/T₀)^0.76
```

### Kinematic Viscosity (ν)

Ratio of dynamic viscosity to density:

```
ν = μ / ρ
```

**SI Units:**
- ν in m²/s
- μ in Pa·s = kg/(m·s)
- ρ in kg/m³

**Common units:**
- cSt (centistokes) = mm²/s
- 1 m²/s = 10⁶ mm²/s = 10⁶ cSt
- 1 cSt = 10⁻⁶ m²/s

**Importance:**
- Used directly in Reynolds number: Re = v·D/ν
- ISO VG oil grades specified in cSt at 40°C
- Measured directly by capillary viscometers

### Viscosity Index (VI)

ASTM D2270 measure of viscosity-temperature relationship:

```
VI = [(L - U) / (L - H)] × 100
```

**Where:**
- U = kinematic viscosity at 40°C (cSt)
- L = viscosity at 40°C for VI=0 oil with same ν₁₀₀
- H = viscosity at 40°C for VI=100 oil with same ν₁₀₀
- ν₁₀₀ = kinematic viscosity at 100°C (cSt)

**L and H from ASTM D2270 tables:**
For ν₁₀₀ from 2 to 70 cSt, polynomial approximations available.

**Interpretation:**
- Higher VI = less viscosity change with temperature
- Mineral oils: VI = 90-110
- Synthetic oils: VI = 120-160
- Multigrade engine oils: VI = 140-180

**Modified VI for VI > 100:**
ASTM D2270 uses different equation for high-VI fluids.

---

## Density Correlations

### Linear Thermal Expansion (Liquids)

First-order approximation for incompressible liquids:

```
ρ(T) = ρ₀ · [1 - β(T - T₀)]
```

**Parameters:**
- ρ₀ = density at reference temperature T₀ (kg/m³)
- β = volumetric thermal expansion coefficient (K⁻¹)
- T, T₀ = temperatures (K or °C, difference is the same)

**Expansion coefficients (typical values at 20°C):**

| Fluid | β (×10⁻⁴ K⁻¹) |
|-------|---------------|
| Water | 2.07 |
| Seawater | 2.5 |
| Ethanol | 11.2 |
| Gasoline | 9.6 |
| Diesel | 8.5 |
| Hydraulic oil | 6.8 |
| Glycerin | 5.0 |
| Mercury | 1.82 |

**Validity:**
- Good for ΔT < 50°C
- Accuracy: ±0.5% for most liquids
- Fails near phase transitions

**Alternative form (volumetric):**
```
V(T) = V₀ · [1 + β(T - T₀)]
```

### Polynomial Fit (Water)

High-accuracy fit for water at atmospheric pressure:

```
ρ(T) = Σ aᵢ · Tⁱ  (i = 0 to n)
```

**5th order polynomial (0-100°C):**
```
ρ = a₀ + a₁T + a₂T² + a₃T³ + a₄T⁴ + a₅T⁵
```

**Coefficients (T in °C, ρ in kg/m³):**
- a₀ = 999.8395
- a₁ = 0.06798
- a₂ = -0.009106
- a₃ = 0.0001005
- a₄ = -0.0000011
- a₅ = 0.0000000065

**Accuracy:** ±0.01% vs IAPWS-95

**Simpler 3rd order (adequate for most engineering):**
```
ρ = 999.972 + 0.04886T - 0.00683T² + 0.0000397T³
```
Accuracy: ±0.1%

### IAPWS-95 Formulation (Water/Steam)

Full equation of state for water and steam:

**Helmholtz free energy formulation:**
```
f(ρ,T) = f⁰(ρ,T) + fʳ(ρ,T)
```

- f⁰ = ideal gas contribution
- fʳ = residual contribution (56 terms)

**Properties derived from f:**
```
P = ρ² ∂f/∂ρ |_T
s = -∂f/∂T |_ρ
u = f + Ts
h = u + P/ρ
```

**Validity:**
- Temperature: 0-1000°C
- Pressure: 0-1000 MPa
- All phases (liquid, vapor, supercritical)

**Accuracy:**
- Density: ±0.01% to ±0.1%
- Best available formulation for water

**Implementation:**
- Complex (56 terms in residual function)
- Use IAPWS-IF97 for industrial applications
- Use CoolProp library for practical calculations

### Ideal Gas Law (Gases)

```
ρ = P·M / (R·T)
```

**Parameters:**
- P = absolute pressure (Pa)
- M = molar mass (kg/mol)
- R = 8.314 J/(mol·K) - universal gas constant
- T = absolute temperature (K)

**Molar masses (kg/mol):**

| Gas | M |
|-----|---|
| Air (dry) | 0.02897 |
| N₂ | 0.02801 |
| O₂ | 0.03200 |
| CO₂ | 0.04401 |
| H₂ | 0.00201 |
| He | 0.00400 |
| Ar | 0.03995 |
| CH₄ | 0.01604 |
| H₂O vapor | 0.01802 |

**Validity:**
- Low to moderate pressure (< 10 bar typical)
- Away from critical point
- Non-condensing conditions

**Compressibility factor correction:**
For non-ideal gases:
```
ρ = P·M / (Z·R·T)
```
where Z = compressibility factor (Z = 1 for ideal gas)

### API Gravity (Petroleum)

American Petroleum Institute scale for oil density:

```
°API = (141.5 / SG₆₀°F) - 131.5
```

**Where:**
- SG₆₀°F = specific gravity at 60°F (15.56°C) vs water at 60°F

**Relationship to density:**
```
ρ₁₅°C = (141.5 / (°API + 131.5)) × 999.016 kg/m³
```

**Classification:**
- Light crude: °API > 31.1
- Medium crude: 22.3 < °API < 31.1
- Heavy crude: 10 < °API < 22.3
- Extra heavy: °API < 10

**Examples:**
- Water: 10 °API
- Light crude: 35-45 °API
- Heavy crude: 10-20 °API
- Gasoline: 60 °API

**Temperature correction (petroleum):**
```
ρ_T = ρ₁₅ - 0.65(T - 15)  [kg/m³, T in °C]
```
Approximate, varies by crude type.

---

## Vapor Pressure Correlations

### Antoine Equation

Most widely used vapor pressure correlation:

```
log₁₀(P) = A - B/(T + C)
```

**Parameters:**
- P = vapor pressure (units vary: mmHg, kPa, bar)
- T = temperature (°C or K, depending on coefficients)
- A, B, C = substance-specific constants

**Common fluids (T in °C, P in mmHg):**

| Substance | A | B | C | Range (°C) |
|-----------|---|---|---|------------|
| Water | 8.07131 | 1730.63 | 233.426 | 1-100 |
| Ethanol | 8.04494 | 1554.3 | 222.65 | 20-93 |
| Methanol | 8.89750 | 1474.08 | 229.13 | -10-65 |
| Acetone | 7.02447 | 1161.0 | 224.0 | -10-55 |
| Benzene | 6.90565 | 1211.033 | 220.790 | 8-80 |
| Toluene | 6.95464 | 1344.800 | 219.482 | 6-137 |
| n-Hexane | 6.87601 | 1171.530 | 224.366 | -26-69 |

**Unit conversion:**
- 1 mmHg = 0.133322 kPa
- 1 atm = 760 mmHg = 101.325 kPa
- 1 bar = 100 kPa

**Extended Antoine (more accurate):**
```
log₁₀(P) = A - B/(T + C) + D·T + E·T² + F·log₁₀(T)
```

**Characteristics:**
- Accuracy: ±1-5% within range
- Extrapolation unreliable
- Different coefficient sets for different ranges

**Data sources:**
- NIST Chemistry WebBook
- Yaws' Handbook
- Perry's Chemical Engineers' Handbook

### Clausius-Clapeyron Equation

Thermodynamic relationship for phase equilibrium:

```
dP/dT = ΔH_vap / (T · ΔV)
```

**For vapor-liquid equilibrium:**
```
dP/dT ≈ ΔH_vap · P / (R · T²)
```

**Integrated form:**
```
ln(P₂/P₁) = -ΔH_vap/R · (1/T₂ - 1/T₁)
```

**Parameters:**
- ΔH_vap = enthalpy of vaporization (J/mol)
- R = 8.314 J/(mol·K)
- T = absolute temperature (K)
- P = vapor pressure (Pa)

**Assumptions:**
- Ideal gas vapor phase
- ΔH_vap constant (good for small ΔT)
- Negligible liquid volume vs vapor

**Use cases:**
- Extrapolation from one known point
- Theoretical estimates
- Temperature correction of vapor pressure

**Heat of vaporization (typical values at normal boiling point):**

| Substance | ΔH_vap (kJ/mol) |
|-----------|-----------------|
| Water | 40.66 |
| Ethanol | 38.56 |
| Methanol | 35.21 |
| Acetone | 29.10 |
| Benzene | 30.72 |
| Ammonia | 23.33 |

### Boiling Point Elevation

Effect of pressure on boiling point:

**From Clausius-Clapeyron:**
```
T_b(P) = T_b,0 / (1 - (R·T_b,0/ΔH_vap)·ln(P/P₀))
```

**Approximate for water near 1 atm:**
```
ΔT_b ≈ 0.03 K per mbar pressure change
```

**Examples:**
- At altitude: P decreases → T_boil decreases
- At 2000 m (0.8 bar): Water boils at ~93°C
- In pressure cooker (2 bar): Water boils at ~121°C

---

## Transport Properties

### Thermal Conductivity (k)

Heat transfer by conduction:

**Fourier's law:**
```
q̇ = -k · A · dT/dx
```

**Typical values at 20°C (W/m·K):**

| Material | k |
|----------|---|
| Air | 0.026 |
| Water | 0.598 |
| Ethanol | 0.169 |
| Engine oil | 0.145 |
| Glycerin | 0.286 |
| Refrigerant R134a | 0.082 |

**Temperature dependence:**
- Liquids: k decreases slightly with T
- Gases: k increases with T (similar to viscosity)

### Specific Heat (Cp)

Heat capacity at constant pressure:

```
Q = m · Cp · ΔT
```

**Typical values at 20°C (kJ/kg·K):**

| Fluid | Cp |
|-------|-----|
| Water | 4.182 |
| Air | 1.005 |
| Ethanol | 2.44 |
| Engine oil | 1.88 |
| Glycerin | 2.43 |
| Ammonia | 4.70 |

**Temperature dependence:**
- Liquids: Cp increases slightly with T
- Water has unusually high Cp (important for cooling)

### Prandtl Number (Pr)

Ratio of momentum diffusivity to thermal diffusivity:

```
Pr = ν / α = μ·Cp / k
```

**Where:**
- ν = kinematic viscosity (m²/s)
- α = thermal diffusivity = k/(ρ·Cp) (m²/s)

**Typical values:**

| Fluid | Pr @ 20°C |
|-------|-----------|
| Water | 7.0 |
| Air | 0.71 |
| Engine oil | 1000-10000 |
| Liquid metals | 0.01-0.1 |

**Interpretation:**
- Pr << 1: Thermal diffusion dominates (liquid metals)
- Pr ≈ 1: Similar momentum and thermal diffusion (gases)
- Pr >> 1: Momentum diffusion dominates (oils)

**Importance in heat transfer:**
- Correlations for convection use Pr as key parameter
- Nusselt number correlations: Nu = f(Re, Pr)

---

## Dimensionless Numbers

### Reynolds Number (Re)

Ratio of inertial to viscous forces:

```
Re = ρ·v·D/μ = v·D/ν
```

**Parameters:**
- ρ = density (kg/m³)
- v = velocity (m/s)
- D = characteristic length (m) - pipe diameter
- μ = dynamic viscosity (Pa·s)
- ν = kinematic viscosity (m²/s)

**Flow regimes in pipes:**
- Laminar: Re < 2300
- Transitional: 2300 < Re < 4000
- Turbulent: Re > 4000

**Critical Reynolds numbers:**
- Pipe flow: Re_crit ≈ 2300
- Flat plate: Re_crit ≈ 5×10⁵
- Sphere: Re_crit ≈ 2×10⁵
- Pump impeller: Re > 10⁵ (always turbulent)

**Applications:**
- Friction factor determination (Moody chart)
- Flow regime identification
- Similarity scaling
- Pump performance prediction

### Friction Factor (f)

Darcy-Weisbach friction factor for pipe flow:

```
ΔP = f · (L/D) · (ρ·v²/2)
```

**Laminar flow (Re < 2300):**
```
f = 64/Re
```

**Turbulent flow smooth pipe (Re > 4000):**

**Blasius (Re < 10⁵):**
```
f = 0.316/Re^0.25
```

**Prandtl-von Karman (Re > 10⁵):**
```
1/√f = 2.0·log₁₀(Re·√f) - 0.8
```

**Colebrook-White (rough pipes):**
```
1/√f = -2.0·log₁₀(ε/(3.7D) + 2.51/(Re·√f))
```

**Where:**
- ε = absolute roughness (m)
- D = pipe diameter (m)
- ε/D = relative roughness

**Typical roughness (mm):**
- Drawn tubing: 0.0015
- Commercial steel: 0.045
- Galvanized iron: 0.15
- Cast iron: 0.26
- Concrete: 0.3-3.0

### Froude Number (Fr)

Ratio of inertial to gravitational forces:

```
Fr = v/√(g·L)
```

**Parameters:**
- v = velocity (m/s)
- g = 9.81 m/s²
- L = characteristic length (m)

**Open channel flow:**
- Subcritical: Fr < 1 (tranquil flow)
- Critical: Fr = 1 (minimum energy)
- Supercritical: Fr > 1 (rapid flow)

**Pump applications:**
- Free surface in sumps and intakes
- Vortex formation
- Air entrainment risk

### Weber Number (We)

Ratio of inertial to surface tension forces:

```
We = ρ·v²·L/σ
```

**Parameters:**
- σ = surface tension (N/m)

**Applications:**
- Droplet formation
- Bubble dynamics
- Atomization in sprays

**Typical surface tension (mN/m at 20°C):**
- Water-air: 72.8
- Ethanol-air: 22.3
- Mercury-air: 486
- Engine oil-air: 35

### Cavitation Number (σ)

Measure of cavitation tendency:

```
σ = (P - P_vap)/(0.5·ρ·v²)
```

**Parameters:**
- P = local static pressure (Pa)
- P_vap = vapor pressure at fluid temperature (Pa)
- v = reference velocity (m/s)

**Interpretation:**
- σ >> 1: No cavitation risk
- σ ≈ 1: Cavitation inception
- σ < 1: Cavitation likely

**Thoma cavitation parameter (pumps):**
```
σ = NPSH_a / H
```
- NPSH_a = net positive suction head available (m)
- H = pump head (m)

---

## Data Sources and Standards

### Primary Databases

#### NIST (National Institute of Standards and Technology)

**NIST Chemistry WebBook**
- URL: https://webbook.nist.gov/chemistry/
- Coverage: 40,000+ compounds
- Properties: Thermodynamic, transport, spectroscopic
- Accuracy: Research-grade, peer-reviewed
- Free access: Yes

**Data available:**
- Vapor pressure (Antoine coefficients)
- Enthalpy of vaporization
- Critical properties
- Phase transition data
- Viscosity (limited coverage)

**NIST REFPROP**
- Commercial software (~$300)
- 147 pure fluids, mixtures
- Equations of state (Helmholtz formulation)
- Highest accuracy available
- Transport properties included

#### IAPWS (International Association for Properties of Water and Steam)

**IAPWS-95**
- Definitive formulation for H₂O
- Helmholtz free energy EOS
- Valid: 0-1000°C, 0-1000 MPa
- Accuracy: 0.01-0.1%

**IAPWS-IF97 (Industrial)**
- Simplified for speed
- Valid: 0-800°C, 0-100 MPa
- Accuracy: 0.1-0.5%
- Used in power plant design

**IAPWS releases:**
- Viscosity (2008)
- Thermal conductivity (2011)
- Surface tension (1994)
- Ice properties (2006)

**Implementation:**
- Use CoolProp library (free)
- freesteam library (C++, free)
- XSteam library (Excel/MATLAB)

### Engineering Handbooks

#### Perry's Chemical Engineers' Handbook
- McGraw-Hill, 9th Edition (2019)
- Section 2: Physical and Chemical Data
- 2000+ pages of property data
- Correlation equations included
- Industry standard reference

**Key sections:**
- 2-95 to 2-141: Prediction and correlation of physical properties
- 2-309 to 2-497: Physical property data tables

#### Yaws' Critical Property Data
- Over 7,900 compounds
- Physical properties and correlations
- Antoine equation coefficients
- Viscosity correlations
- Electronic database available

#### CRC Handbook of Chemistry and Physics
- 100+ editions since 1913
- Fluid properties section
- Quick reference data
- Online version available

#### ASHRAE Handbooks
- American Society of Heating, Refrigerating, and Air-Conditioning Engineers
- Four volumes (updated on 4-year cycle):
  - Fundamentals (psychrometrics, fluid properties)
  - Refrigeration (refrigerants, brines)
  - HVAC Systems and Equipment
  - HVAC Applications

**Refrigerant data:**
- Pressure-enthalpy diagrams
- Saturation tables
- Superheated vapor tables
- Transport properties

### Standards Organizations

#### ASTM International

**Viscosity standards:**
- **ASTM D341**: Viscosity-temperature charts for liquid petroleum products
- **ASTM D445**: Kinematic viscosity of transparent and opaque liquids
- **ASTM D2270**: Calculating viscosity index from kinematic viscosity
- **ASTM D446**: Specifications for glass capillary kinematic viscometers
- **ASTM D7042**: Dynamic viscosity and density by stabinger viscometer

**Density standards:**
- **ASTM D1298**: Density, relative density, or API gravity by hydrometer
- **ASTM D4052**: Density by digital density meter
- **ASTM D5002**: Density and relative density by digital density meter

**Petroleum standards:**
- **ASTM D6751**: Biodiesel fuel blend stock specification
- **ASTM D4057**: Manual sampling of petroleum liquids

#### ISO (International Organization for Standardization)

**ISO 3448**: Industrial liquid lubricants - ISO viscosity classification
- Defines VG grades (VG 2 through VG 3200)
- Mid-point viscosity at 40°C in cSt
- ±10% tolerance

**ISO VG grades:**
| Grade | ν @ 40°C (cSt) | Min | Max |
|-------|----------------|-----|-----|
| VG 10 | 10 | 9.0 | 11.0 |
| VG 15 | 15 | 13.5 | 16.5 |
| VG 22 | 22 | 19.8 | 24.2 |
| VG 32 | 32 | 28.8 | 35.2 |
| VG 46 | 46 | 41.4 | 50.6 |
| VG 68 | 68 | 61.2 | 74.8 |
| VG 100 | 100 | 90.0 | 110.0 |

**Other ISO standards:**
- **ISO 12185**: Crude petroleum density determination
- **ISO 2909**: Petroleum measurement tables
- **ISO 4259**: Petroleum products precision of measurement

#### SAE (Society of Automotive Engineers)

**SAE J300**: Engine oil viscosity classification
- Single grade: SAE 20, 30, 40, 50, 60
- Multigrade: SAE 0W-40, 5W-30, 10W-40, 15W-40, 20W-50

**Viscosity requirements (SAE J300):**
- "W" grades: max viscosity at cold temperature (-40°C to -10°C)
- High temperature: min viscosity at 100°C (cSt)

#### API (American Petroleum Institute)

**API gravity scale:**
- Characterizes crude oil and petroleum products
- Higher °API = lighter oil
- Water = 10 °API

**API standards:**
- API 2540: Petroleum measurement tables
- API 11N: Care and use of subsurface pumps

### Software Libraries

#### CoolProp (Open Source)
- URL: http://www.coolprop.org/
- Language: C++ with Python, MATLAB, Excel wrappers
- Coverage: 122 pure fluids, mixtures
- License: MIT (free, open source)
- Accuracy: Comparable to REFPROP for most fluids
- See `coolprop-db` skill for usage

**Installation:**
```bash
pip install CoolProp
```

**Key features:**
- High-accuracy equations of state
- Transport properties (viscosity, thermal conductivity)
- All phases (liquid, vapor, supercritical, two-phase)
- IAPWS-95 for water

#### thermo (Python package)
- Pure Python implementation
- 20,000+ chemicals
- Focus on organic compounds
- Includes mixing rules

**Installation:**
```bash
pip install thermo
```

#### chemicals (Python package)
- Lightweight, pure Python
- Property correlations and data
- No complex EOS

---

## Standard Conditions Definitions

### Standard Temperature and Pressure (STP)

**IUPAC (Current standard):**
- Temperature: 0°C (273.15 K)
- Pressure: 100 kPa (1 bar)
- Molar volume: 22.711 L/mol (ideal gas)

**NIST/IUPAC (Old standard, pre-1982):**
- Temperature: 0°C (273.15 K)
- Pressure: 101.325 kPa (1 atm)
- Molar volume: 22.414 L/mol (ideal gas)

**Note:** Always verify which standard is being used!

### Normal Temperature and Pressure (NTP)

**Common definition:**
- Temperature: 20°C (293.15 K)
- Pressure: 101.325 kPa (1 atm)

**Alternative (ISO):**
- Temperature: 25°C (298.15 K)
- Pressure: 100 kPa (1 bar)

### Standard Ambient Temperature and Pressure (SATP)

**IUPAC (1982):**
- Temperature: 25°C (298.15 K)
- Pressure: 100 kPa (1 bar)
- Molar volume: 24.790 L/mol (ideal gas)

### Standard Conditions for Petroleum

**API/ASTM:**
- Temperature: 60°F (15.56°C)
- Pressure: 14.696 psia (101.325 kPa)

**ISO:**
- Temperature: 15°C (288.15 K)
- Pressure: 101.325 kPa (1 atm)

### Water Reference Conditions

**Density maximum:**
- Temperature: 3.98°C
- Pressure: 101.325 kPa
- Density: 999.972 kg/m³

**Steam tables reference:**
- Triple point: 0.01°C, 0.6117 kPa
- Often use 0°C or 25°C as reference for h, s

---

## Conversion Factors

### Pressure

| From | To | Multiply by |
|------|-----|-------------|
| Pa | kPa | 0.001 |
| kPa | bar | 0.01 |
| bar | psi | 14.504 |
| psi | kPa | 6.895 |
| atm | kPa | 101.325 |
| mmHg (Torr) | kPa | 0.133322 |
| inH₂O | Pa | 249.09 |

**Absolute vs gauge:**
- Absolute pressure = gauge pressure + atmospheric pressure
- P_abs (kPa) = P_gauge (kPa) + 101.325 kPa (at sea level)

### Viscosity

**Dynamic viscosity:**

| From | To | Multiply by |
|------|-----|-------------|
| Pa·s | mPa·s (cP) | 1000 |
| mPa·s | μPa·s | 1000 |
| lbf·s/ft² | Pa·s | 47.88 |
| poise (P) | Pa·s | 0.1 |
| cP | mPa·s | 1.0 |

**Kinematic viscosity:**

| From | To | Multiply by |
|------|-----|-------------|
| m²/s | cSt (mm²/s) | 10⁶ |
| cSt | m²/s | 10⁻⁶ |
| ft²/s | m²/s | 0.0929 |
| stokes (St) | m²/s | 10⁻⁴ |
| cSt | St | 0.01 |

### Density

| From | To | Multiply by |
|------|-----|-------------|
| kg/m³ | g/cm³ | 0.001 |
| g/cm³ | kg/L | 1.0 |
| lb/ft³ | kg/m³ | 16.018 |
| lb/gal (US) | kg/L | 0.1198 |
| SG (vs water) | kg/m³ | 1000 (approx) |

### Temperature

**Formulas:**
```
K = °C + 273.15
°C = (°F - 32) × 5/9
°F = °C × 9/5 + 32
°R = °F + 459.67
```

**Common values:**
- Absolute zero: 0 K = -273.15°C = -459.67°F
- Ice point: 273.15 K = 0°C = 32°F
- Steam point: 373.15 K = 100°C = 212°F

---

## Verification Data

### Water Properties Validation (IAPWS-95)

**At 25°C, 101.325 kPa:**

| Property | Value | Unit | Uncertainty |
|----------|-------|------|-------------|
| Density | 997.0 | kg/m³ | ±0.01% |
| Viscosity | 0.8903 | mPa·s | ±0.3% |
| Thermal conductivity | 0.6065 | W/m·K | ±1% |
| Specific heat | 4.1813 | kJ/kg·K | ±0.1% |
| Prandtl number | 6.14 | - | ±1% |
| Vapor pressure | 3.1688 | kPa | ±0.02% |

**At 100°C, 101.325 kPa (saturated liquid):**

| Property | Value | Unit |
|----------|-------|------|
| Density | 958.4 | kg/m³ |
| Viscosity | 0.282 | mPa·s |
| Vapor pressure | 101.325 | kPa |
| Latent heat | 2257 | kJ/kg |

### ISO VG 46 Oil Validation

**Typical mineral oil, VI ≈ 95:**

| Temperature | Viscosity (cSt) | Source |
|-------------|-----------------|--------|
| 0°C | 370 | Walther equation |
| 40°C | 46.0 | ISO 3448 definition |
| 100°C | 6.8 | Extrapolated |

**Density:**
- 15°C: 870 kg/m³
- 40°C: 854 kg/m³
- 100°C: 815 kg/m³

### Air Properties Validation (NIST)

**At 25°C, 101.325 kPa:**

| Property | Value | Unit |
|----------|-------|------|
| Density | 1.184 | kg/m³ |
| Viscosity | 18.46 | μPa·s |
| Thermal conductivity | 0.02624 | W/m·K |
| Specific heat | 1.005 | kJ/kg·K |
| Prandtl number | 0.707 | - |

**Sutherland's law validation:**
- μ(0°C) = 17.16 μPa·s (exact by definition)
- μ(25°C) = 18.46 μPa·s (error < 0.5%)
- μ(100°C) = 21.86 μPa·s (error < 1%)

---

## Best Practices for Property Evaluation

### 1. Select Appropriate Correlation

**For liquids:**
- Simple, narrow range: Andrade equation
- Wide range, oils: Walther equation (ASTM D341)
- High accuracy, water: IAPWS-IF97 or CoolProp
- Refrigerants: Always use CoolProp or REFPROP

**For gases:**
- Ideal gases: Sutherland's law
- High pressure: Use real EOS (CoolProp, REFPROP)

### 2. Verify Validity Range

Always check:
- Temperature range of correlation
- Pressure range (especially for gases)
- Phase region (liquid, vapor, two-phase)

### 3. Check Against Multiple Sources

For critical applications:
- Cross-reference at least two independent sources
- Compare against experimental data if available
- Use highest-accuracy source (NIST, IAPWS) when possible

### 4. Document Assumptions

Always record:
- Property source (NIST, CoolProp, correlation)
- Temperature and pressure conditions
- Units (SI strongly recommended)
- Date accessed (for online databases)
- Assumptions (purity, composition, phase)

### 5. Sensitivity Analysis

For design calculations:
- Vary properties by ±10% (or known uncertainty)
- Check impact on key results (NPSH, power, Re)
- Design for worst-case conditions

### 6. Measurement Validation

When possible:
- Measure critical properties (viscosity, density)
- Use calibrated instruments (viscometer, hydrometer)
- Compare measurements to predictions
- Account for impurities and aging (especially oils)

---

## References

### Key Papers

1. **Vogel, H.** (1921). "Das Temperaturabhängigkeitsgesetz der Viskosität von Flüssigkeiten." *Physikalische Zeitschrift*, 22, 645-646.

2. **Sutherland, W.** (1893). "The viscosity of gases and molecular force." *Philosophical Magazine*, Series 5, 36(223), 507-531.

3. **Antoine, C.** (1888). "Tensions des vapeurs; nouvelle relation entre les tensions et les températures." *Comptes Rendus des Séances de l'Académie des Sciences*, 107, 681-684, 778-780, 836-837.

4. **Wagner, W., & Pruß, A.** (2002). "The IAPWS formulation 1995 for the thermodynamic properties of ordinary water substance for general and scientific use." *Journal of Physical and Chemical Reference Data*, 31(2), 387-535.

5. **Huber, M. L., et al.** (2009). "New international formulation for the viscosity of H₂O." *Journal of Physical and Chemical Reference Data*, 38(2), 101-125.

### Online Resources

- **NIST Chemistry WebBook**: https://webbook.nist.gov/chemistry/
- **IAPWS Official Site**: http://www.iapws.org/
- **CoolProp Documentation**: http://www.coolprop.org/
- **ASTM Standards**: https://www.astm.org/
- **ISO Standards**: https://www.iso.org/

---

*This reference document provides comprehensive technical details for material property evaluation in engineering applications. Always verify critical data against primary sources.*
