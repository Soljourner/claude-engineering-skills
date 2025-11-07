# NIST REFPROP Database Reference

Comprehensive reference for NIST REFPROP thermophysical property database queries.

## Table of Contents

1. [Supported Fluids](#supported-fluids)
2. [Property Codes](#property-codes)
3. [Input-Output Specifications](#input-output-specifications)
4. [Equation of State Information](#equation-of-state-information)
5. [Mixture Calculations](#mixture-calculations)
6. [Comparison with CoolProp](#comparison-with-coolprop)
7. [Unit Systems](#unit-systems)
8. [API Function Reference](#api-function-reference)
9. [External Resources](#external-resources)

---

## Supported Fluids

NIST REFPROP Version 10.0 includes **147 pure fluids** with reference-quality equations of state.

### Water and Air
- **WATER** - H₂O (IAPWS-95 formulation, highest accuracy available)
- **AIR** - Dry air (Lemmon et al. pseudo-pure mixture)

### Noble Gases (6 fluids)
- **HELIUM** - He
- **NEON** - Ne
- **ARGON** - Ar
- **KRYPTON** - Kr
- **XENON** - Xe
- **RADON** - Rn

### Light Gases (7 fluids)
- **HYDROGEN** - H₂
- **PARAHYD** - Para-hydrogen
- **ORTHOHYD** - Ortho-hydrogen
- **DEUTERIUM** - D₂ (Heavy hydrogen)
- **NITROGEN** - N₂
- **OXYGEN** - O₂
- **FLUORINE** - F₂

### Carbon and Sulfur Compounds (6 fluids)
- **CO** - Carbon monoxide
- **CO2** - Carbon dioxide (R744)
- **COS** - Carbonyl sulfide
- **SO2** - Sulfur dioxide (R764)
- **H2S** - Hydrogen sulfide
- **SF6** - Sulfur hexafluoride

### Alkanes - Normal (18 fluids)
- **METHANE** - CH₄ (C1)
- **ETHANE** - C₂H₆ (C2)
- **PROPANE** - C₃H₈ (C3, R290)
- **BUTANE** - C₄H₁₀ (n-C4, R600)
- **PENTANE** - C₅H₁₂ (n-C5)
- **HEXANE** - C₆H₁₄ (n-C6)
- **HEPTANE** - C₇H₁₆ (n-C7)
- **OCTANE** - C₈H₁₈ (n-C8)
- **NONANE** - C₉H₂₀ (n-C9)
- **DECANE** - C₁₀H₂₂ (n-C10)
- **UNDECANE** - C₁₁H₂₄ (n-C11)
- **DODECANE** - C₁₂H₂₆ (n-C12)
- **C13**, **C14**, **C15**, **C16**, **C20**, **C22**

### Alkanes - Branched (8 fluids)
- **ISOBUTAN** - Isobutane (i-C4, R600a)
- **IPENTANE** - Isopentane (i-C5)
- **IHEXANE** - 2-Methylpentane
- **IOCTANE** - Isooctane (2,2,4-Trimethylpentane)
- **NEOPENTN** - Neopentane
- **C12_2METYL** - 2-Methylundecane
- **C12_5METYL** - 5-Methylundecane
- **C12_35DIME** - 3,5-Dimethyldecane

### Alkenes (9 fluids)
- **ETHYLENE** - Ethene (C₂H₄)
- **PROPENE** - Propylene (C₃H₆, R1270)
- **PROPYLEN** - Propylene (alias)
- **1BUTENE** - 1-Butene
- **IBUTENE** - Isobutene
- **C2BUTENE** - cis-2-Butene
- **T2BUTENE** - trans-2-Butene
- **C3CC6** - 1-Hexene
- **C4CC6** - 1-Heptene

### Cyclic Hydrocarbons (7 fluids)
- **CYCLOPEN** - Cyclopentane (C₅H₁₀)
- **CYCLOHEX** - Cyclohexane (C₆H₁₂)
- **CYCLOPRO** - Cyclopropane (C₃H₆)
- **BENZENE** - Benzene (C₆H₆)
- **TOLUENE** - Toluene (C₇H₈)
- **MXYLENE** - m-Xylene (C₈H₁₀)
- **OXYLENE** - o-Xylene (C₈H₁₀)
- **PXYLENE** - p-Xylene (C₈H₁₀)
- **EBENZENE** - Ethylbenzene (C₈H₁₀)

### Alcohols (6 fluids)
- **METHANOL** - CH₃OH
- **ETHANOL** - C₂H₅OH
- **PROPANOL** - 1-Propanol (n-Propanol)
- **IPROPANOL** - 2-Propanol (Isopropanol)
- **BUTANOL** - 1-Butanol
- **IBUTANOL** - 2-Methyl-1-propanol

### CFCs - Chlorofluorocarbons (9 fluids)
- **R11** - CFCl₃ (Trichlorofluoromethane)
- **R12** - CF₂Cl₂ (Dichlorodifluoromethane)
- **R13** - CClF₃ (Chlorotrifluoromethane)
- **R14** - CF₄ (Tetrafluoromethane)
- **R113** - C₂Cl₃F₃
- **R114** - C₂Cl₂F₄
- **R115** - C₂ClF₅
- **R116** - C₂F₆ (Hexafluoroethane)
- **RC318** - C₄F₈ (Octafluorocyclobutane)

### HCFCs - Hydrochlorofluorocarbons (6 fluids)
- **R21** - CHCl₂F
- **R22** - CHClF₂ (Chlorodifluoromethane)
- **R123** - CHCl₂CF₃
- **R124** - CHClFCF₃
- **R141B** - CH₃CCl₂F
- **R142B** - CH₃CClF₂

### HFCs - Hydrofluorocarbons (20 fluids)
- **R23** - CHF₃ (Trifluoromethane)
- **R32** - CH₂F₂ (Difluoromethane)
- **R41** - CH₃F (Fluoromethane)
- **R125** - CHF₂CF₃ (Pentafluoroethane)
- **R134A** - CH₂FCF₃ (1,1,1,2-Tetrafluoroethane) - **Most common**
- **R143A** - CH₃CF₃
- **R152A** - CH₃CHF₂
- **R161** - CH₃CH₂F
- **R227EA** - CF₃CHFCF₃
- **R236EA** - CHF₂CHFCF₃
- **R236FA** - CF₃CH₂CF₃
- **R245CA** - CH₂FCF₂CHF₂
- **R245FA** - CHF₂CH₂CF₃
- **R365MFC** - CH₃CF₂CH₂CF₃
- **R4310MEE** - CF₃CHFCHFCF₂CF₃

### HFOs - Hydrofluoroolefins (Low GWP) (8 fluids)
- **R1234YF** - CF₃CF=CH₂ (GWP < 1, automotive AC)
- **R1234ZE** - CF₃CH=CHF (trans-1,3,3,3-Tetrafluoropropene)
- **R1233ZDE** - CF₃CH=CHCl (trans-1-Chloro-3,3,3-trifluoropropene)
- **R1336MZZZ** - (Z)-1,1,1,4,4,4-Hexafluoro-2-butene
- **R1243ZF** - CF₃CH=CH₂
- **R1216** - C₃F₆ (Hexafluoropropylene)

### Natural Refrigerants (5 fluids)
- **R717** - NH₃ (Ammonia)
- **R744** - CO₂ (Carbon dioxide, transcritical systems)
- **R290** - C₃H₈ (Propane)
- **R600** - C₄H₁₀ (n-Butane)
- **R600A** - i-C₄H₁₀ (Isobutane)
- **R1270** - C₃H₆ (Propylene)

### Predefined Refrigerant Blends (20+ mixtures)

**Zeotropic Blends (Temperature Glide):**
- **R404A** - R125/143a/134a (44/52/4%)
- **R407A** - R32/125/134a (20/40/40%)
- **R407C** - R32/125/134a (23/25/52%)
- **R407F** - R32/125/134a (30/30/40%)
- **R410A** - R32/125 (50/50%) - **Most common residential AC**
- **R507A** - R125/143a (50/50%)
- **R448A** - R32/125/1234yf/134a/1234ze(E) (26/26/20/21/7%)
- **R449A** - R32/125/1234yf/134a (24.3/24.7/25.3/25.7%)
- **R450A** - R134a/1234ze(E) (42/58%)
- **R513A** - R1234yf/134a (56/44%)

**Azeotropic Blends:**
- **R500** - R12/152a (73.8/26.2%)
- **R502** - R22/115 (48.8/51.2%)
- **R503** - R23/13 (40.1/59.9%)
- **R504** - R32/115 (48.2/51.8%)
- **R508B** - R23/116 (46/54%)

### Ammonia and Derivatives (3 fluids)
- **AMMONIA** - NH₃ (R717)
- **NH3** - Ammonia (alias)
- **HYDRAZIN** - N₂H₄ (Hydrazine)

### Siloxanes (ORC Applications) (9 fluids)
- **D4** - Octamethylcyclotetrasiloxane
- **D5** - Decamethylcyclopentasiloxane
- **D6** - Dodecamethylcyclohexasiloxane
- **MD2M** - Decamethyltetrasiloxane
- **MD3M** - Dodecamethylpentasiloxane
- **MD4M** - Tetradecamethylhexasiloxane
- **MDM** - Octamethyltrisiloxane
- **MM** - Hexamethyldisiloxane
- **MDM** - Octamethyltrisiloxane

### Other Organic Compounds (10 fluids)
- **ACETONE** - CH₃COCH₃
- **ACETYLENE** - C₂H₂
- **DME** - CH₃OCH₃ (Dimethyl ether, RE170)
- **DMC** - Dimethyl carbonate
- **DEE** - C₄H₁₀O (Diethyl ether)
- **ETHER** - Diethyl ether (alias)
- **ETHYLENE** - C₂H₄
- **ETHYLENEOXIDE** - C₂H₄O
- **NEOPENTAN** - C₅H₁₂ (Neopentane)
- **METHYL** - Methyl (various derivatives)

### Heavy Fluids and Other (5 fluids)
- **HEAVYWAT** - D₂O (Heavy water)
- **THERMINOL66** - Heat transfer fluid
- **T2BUTENE** - trans-2-Butene
- **NF3** - Nitrogen trifluoride
- **HCL** - Hydrogen chloride

---

## Property Codes

### Output Property Codes

These codes are used in the output specification string (3rd parameter):

| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `T` | Temperature | K | Absolute temperature |
| `P` | Pressure | kPa | Absolute pressure |
| `D` | Density | kg/m³ | Mass density |
| `E` | Internal energy | J/kg | Specific internal energy |
| `H` | Enthalpy | J/kg | Specific enthalpy |
| `S` | Entropy | J/kg/K | Specific entropy |
| `CV` | Isochoric heat capacity | J/kg/K | Cp at constant volume |
| `CP` | Isobaric heat capacity | J/kg/K | Cp at constant pressure |
| `W` | Speed of sound | m/s | Sonic velocity |
| `Z` | Compressibility factor | - | Z = PV/(RT) |
| `JT` | Joule-Thomson coefficient | K/kPa | μ = (∂T/∂P)ₕ |
| `Q` | Quality | - | Vapor mass fraction (0-1) |
| `PHASE` | Phase index | - | Phase identifier |

### Transport Property Codes

| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `VIS` | Dynamic viscosity | Pa·s | η, μ |
| `TCX` | Thermal conductivity | W/m/K | λ, k |
| `KV` | Kinematic viscosity | m²/s | ν = μ/ρ |
| `STN` | Surface tension | N/m | σ |
| `PRANDTL` | Prandtl number | - | Pr = μCp/k |

### Derivative Properties

| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `DPDT` | dP/dT at const ρ | kPa/K | Isochoric pressure derivative |
| `DPDL` | dP/dρ at const T | kPa·m³/kg | Isothermal pressure derivative |
| `D2PDT2` | d²P/dT² at const ρ | kPa/K² | Second derivative |
| `DHDT_D` | dH/dT at const ρ | J/kg/K | Enthalpy temperature derivative |
| `DHDT_P` | dH/dT at const P | J/kg/K | Same as Cp |
| `DHDP_T` | dH/dP at const T | J/kg/kPa | Enthalpy pressure derivative |
| `CHEMPOT` | Chemical potential | J/mol | μ |
| `FUGACITY` | Fugacity | kPa | f |
| `KAPPA` | Isothermal compressibility | 1/kPa | κ = -1/V(∂V/∂P)ₜ |
| `BETA` | Volume expansivity | 1/K | β = 1/V(∂V/∂T)ₚ |

### Critical and Reference Point Properties

| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `TCRIT` | Critical temperature | K | Tс |
| `PCRIT` | Critical pressure | kPa | Pс |
| `DCRIT` | Critical density | kg/m³ | ρс |
| `TTRP` | Triple point temperature | K | Tₜ |
| `PTRP` | Triple point pressure | kPa | Pₜ |
| `TNBP` | Normal boiling point | K | At 101.325 kPa |
| `TMIN` | Minimum temperature | K | Lower validity limit |
| `TMAX` | Maximum temperature | K | Upper validity limit |
| `PMAX` | Maximum pressure | kPa | Upper validity limit |

### Molecular Properties

| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `MM` | Molar mass | kg/mol | M, molecular weight |
| `ACF` | Acentric factor | - | ω (Pitzer) |
| `DIPOLE` | Dipole moment | Debye | μ |
| `RGAS` | Gas constant | J/mol/K | R = 8.314... J/mol/K |

### Environmental Properties

| Code | Property | SI Unit | Description |
|------|----------|---------|-------------|
| `ODP` | Ozone depletion potential | - | Relative to R11 |
| `GWP` | Global warming potential | - | 100-year ITH |
| `FH` | ASHRAE flammability | - | 1=Low, 2=Moderate, 3=High |
| `HH` | ASHRAE health hazard | - | A=Low, B=High |
| `NFL` | ASHRAE safety class | - | e.g., A1, A2L, B2 |

---

## Input-Output Specifications

### Input Property Pairs

The 2nd parameter to `REFPROPdll()` specifies the input property pair:

| Code | Properties | Units | Valid Regions |
|------|------------|-------|---------------|
| `TP` | Temperature, Pressure | K, kPa | Single-phase only |
| `TH` | Temperature, Enthalpy | K, J/kg | All regions |
| `TS` | Temperature, Entropy | K, J/kg/K | All regions |
| `TE` | Temperature, Internal energy | K, J/kg | All regions |
| `TD` | Temperature, Density | K, kg/m³ | All regions |
| `PH` | Pressure, Enthalpy | kPa, J/kg | All regions |
| `PS` | Pressure, Entropy | kPa, J/kg/K | All regions |
| `PE` | Pressure, Internal energy | kPa, J/kg | All regions |
| `PD` | Pressure, Density | kPa, kg/m³ | All regions |
| `TQ` | Temperature, Quality | K, - | Two-phase (0≤Q≤1) |
| `PQ` | Pressure, Quality | kPa, - | Two-phase (0≤Q≤1) |
| `DH` | Density, Enthalpy | kg/m³, J/kg | All regions |
| `DS` | Density, Entropy | kg/m³, J/kg/K | All regions |
| `DE` | Density, Internal energy | kg/m³, J/kg | All regions |
| `HS` | Enthalpy, Entropy | J/kg, J/kg/K | All regions |

### Multiple Output Request

Request multiple properties in one call by separating with semicolons:

```python
# Request density, enthalpy, entropy, Cp, viscosity, conductivity
output = "D;H;S;CP;VIS;TCX"
r = RP.REFPROPdll(fluid, "TP", output, MOLAR_BASE_SI, 0, 0, T, P, z)

density = r.Output[0]
enthalpy = r.Output[1]
entropy = r.Output[2]
cp = r.Output[3]
viscosity = r.Output[4]
conductivity = r.Output[5]
```

---

## Equation of State Information

### Helmholtz Energy Formulation

REFPROP uses explicit Helmholtz energy equations of state for most pure fluids:

```
a(δ, τ) = a⁰(δ, τ) + aʳ(δ, τ)
```

Where:
- `a = A/(RT)` - Dimensionless Helmholtz energy
- `δ = ρ/ρc` - Reduced density
- `τ = Tc/T` - Inverse reduced temperature
- `a⁰` - Ideal gas contribution
- `aʳ` - Residual (real fluid) contribution

All thermodynamic properties are derived from partial derivatives of this function with respect to δ and τ.

### Reference-Quality Equations

REFPROP contains the most accurate equations available:

#### Water - IAPWS-95
- **Source**: Wagner & Pruß (2002)
- **Accuracy**: ±0.001% for density, ±0.01% for heat capacity
- **Range**: 251.2 K to 1273 K, up to 1000 MPa
- **Publication**: J. Phys. Chem. Ref. Data, 31, 387-535

#### Refrigerants
Most refrigerants use equations from:
- **Tillner-Roth equations**: R134a (±0.01% density)
- **Lemmon-Span short equations**: R32, R125, R143a
- **Extended corresponding states**: Older refrigerants

#### Hydrocarbons
- **Natural gas**: GERG-2008 equation (Kunz & Wagner)
- **Accuracy**: ±0.1% for density of natural gas mixtures
- **Pure alkanes**: Span, Lemmon, Thol equations

#### HFOs (Low-GWP Refrigerants)
- **R1234yf**: Richter et al. (2011) - ±0.1% density
- **R1234ze(E)**: Thol & Lemmon (2016)
- **R1233zd(E)**: Mondejar et al. (2015)

### Equation Types by Fluid Class

| Fluid Class | Equation Type | Typical Accuracy (Density) |
|-------------|---------------|----------------------------|
| Water | IAPWS-95 | ±0.001% |
| Air | Lemmon et al. | ±0.1% |
| Hydrocarbons | Span, Lemmon equations | ±0.05-0.2% |
| HFCs | Tillner-Roth, Lemmon | ±0.01-0.1% |
| HFOs | Modern multi-parameter | ±0.1-0.2% |
| Natural gases | GERG-2008 | ±0.1% |
| Siloxanes | Colonna et al. | ±0.2-0.5% |

### Mixing Rules for Mixtures

REFPROP uses advanced mixing rules:

- **GERG-2008**: For natural gas mixtures (up to 21 components)
- **Kunz-Wagner**: High-accuracy for LNG
- **Lemmon-Jacobsen**: General hydrocarbons
- **Tillner-Roth & Friend**: Ammonia-water mixtures
- **Extended Corresponding States (ECS)**: Default for other mixtures

### Transport Properties

Transport properties use separate correlations:

- **Viscosity**: Extended corresponding states, fluid-specific correlations
- **Thermal conductivity**: Enhanced ECS model, critical enhancement
- **Surface tension**: Correlation based on molecular structure

**Typical Uncertainties:**
- Viscosity: ±1-5%
- Thermal conductivity: ±2-5%
- Surface tension: ±0.5-2%

---

## Mixture Calculations

### Defining Mixtures

```python
# Binary mixture
fluid = "METHANE;ETHANE"
z = [0.9, 0.1]  # Mole fractions (must sum to 1.0)

# Ternary mixture
fluid = "NITROGEN;METHANE;ETHANE"
z = [0.05, 0.90, 0.05]

# Complex mixture (natural gas)
fluid = "NITROGEN;METHANE;ETHANE;PROPANE;NBUTANE;ISOBUTAN;CO2"
z = [0.01, 0.92, 0.04, 0.02, 0.005, 0.005, 0.01]
```

### Predefined Refrigerant Blends

REFPROP includes predefined blends with optimized mixing parameters:

```python
# Use predefined blend name
fluid = "R410A.MIX"  # or just "R410A"
z = [1.0]  # Composition already defined in .mix file

# Available blend files: R404A.MIX, R407C.MIX, R410A.MIX, R507A.MIX, etc.
```

### Temperature Glide

Zeotropic mixtures exhibit temperature glide (temperature changes during phase change):

```python
# R407C at 10 bar
fluid = "R407C.MIX"
P = 1000.0  # kPa

# Bubble point (Q=0)
r = RP.REFPROPdll(fluid, "PQ", "T", MOLAR_BASE_SI, 0, 0, P, 0.0, [1.0])
T_bubble = r.Output[0]

# Dew point (Q=1)
r = RP.REFPROPdll(fluid, "PQ", "T", MOLAR_BASE_SI, 0, 0, P, 1.0, [1.0])
T_dew = r.Output[0]

glide = T_dew - T_bubble  # Temperature glide in K
```

### VLE Calculations

Vapor-liquid equilibrium for mixtures:

```python
# Use TPFLSH for VLE calculation
r = RP.TPFLSHdll(T, P, z)

# Returns both vapor and liquid compositions
x_liquid = r.x  # Liquid mole fractions
y_vapor = r.y   # Vapor mole fractions
Q = r.q        # Overall quality
```

---

## Comparison with CoolProp

### Accuracy Comparison

| Property | REFPROP | CoolProp | Notes |
|----------|---------|----------|-------|
| Density (liquid) | ±0.01% | ±0.1% | REFPROP 10× more accurate |
| Vapor pressure | ±0.02% | ±0.1% | REFPROP better near critical point |
| Heat capacity | ±0.5% | ±1% | Both good for most applications |
| Viscosity | ±1-2% | ±2-5% | REFPROP has better correlations |
| Thermal conductivity | ±2% | ±3-5% | REFPROP slightly better |
| Speed of sound | ±0.1% | ±0.5% | REFPROP significantly better |

### Fluid Coverage

| Category | REFPROP | CoolProp | Advantage |
|----------|---------|----------|-----------|
| Pure fluids | 147 | ~122 | REFPROP: +25 fluids |
| Refrigerants | 60+ | 40+ | REFPROP: Latest HFOs |
| Hydrocarbons | 50+ | 30+ | REFPROP: More isomers |
| Siloxanes | 9 | 9 | Similar |
| Mixtures | Unlimited | Limited | REFPROP: Advanced mixing |
| Predefined blends | 20+ | 10+ | REFPROP: More blends |

### Mixture Capabilities

| Feature | REFPROP | CoolProp |
|---------|---------|----------|
| Binary mixtures | ✓ Excellent | ✓ Good |
| Multi-component (3-20) | ✓ Excellent | ✗ Limited |
| GERG-2008 (natural gas) | ✓ Yes | ✓ Yes (limited) |
| Temperature glide | ✓ Full support | ~ Approximate |
| VLE calculations | ✓ Advanced | ~ Basic |
| Custom mixing rules | ✓ Multiple options | ✗ No |

### Validation

- **CoolProp is validated against REFPROP**
- REFPROP serves as the reference standard
- Deviations reported in CoolProp documentation are relative to REFPROP
- For most engineering applications (±1% accuracy), CoolProp is sufficient
- For research, critical design, or publication, REFPROP is preferred

### Cost-Benefit Analysis

| Factor | REFPROP | CoolProp |
|--------|---------|----------|
| Cost | ~$300 USD | Free |
| Accuracy | Highest (reference) | Excellent (0.1-1%) |
| Support | NIST official | Community |
| Updates | Regular (5+ years) | Continuous |
| Mixtures | Superior | Limited |
| Licensing | Commercial | Open-source (MIT) |
| Use case | Critical/research | General engineering |

**Recommendation:**
- **Use REFPROP** for: Critical design, mixtures, highest accuracy, publication work, regulatory compliance
- **Use CoolProp** for: General engineering, pure fluids, budget constraints, open-source requirements

---

## Unit Systems

REFPROP supports multiple unit systems via the `iUnits` parameter:

### Unit System Codes

```python
# Get unit system enumerations
MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum  # 21
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum    # 20
SI_WITH_C = RP.GETENUMdll(0, "SI WITH C").iEnum          # 2
ENGLISH = RP.GETENUMdll(0, "ENGLISH").iEnum              # 1
```

### Unit System Details

| System | Code | T | P | ρ | h | s |
|--------|------|---|---|---|---|---|
| SI (mass) | 20 | K | kPa | kg/m³ | J/kg | J/kg/K |
| SI (molar) | 21 | K | kPa | mol/L | J/mol | J/mol/K |
| SI with °C | 2 | °C | kPa | kg/m³ | kJ/kg | kJ/kg/K |
| English | 1 | °F | psia | lbm/ft³ | BTU/lbm | BTU/lbm/°R |

### Recommended Practice

**Use MOLAR_BASE_SI or MASS_BASE_SI** for consistency with scientific literature and avoid unit conversion errors.

```python
# Recommended initialization
MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
r = RP.REFPROPdll(fluid, "TP", "D;H;S", MOLAR_BASE_SI, 0, 0, T, P, z)
```

---

## API Function Reference

### REFPROPdll - Primary Function

```python
r = RP.REFPROPdll(hFld, hIn, hOut, iUnits, iMass, iFlag, a, b, z)
```

**Parameters:**
- `hFld` (str): Fluid name or mixture (e.g., "WATER", "METHANE;ETHANE")
- `hIn` (str): Input property pair (e.g., "TP", "PH", "PS")
- `hOut` (str): Output properties (e.g., "D;H;S;CP;VIS")
- `iUnits` (int): Unit system (use `MOLAR_BASE_SI`)
- `iMass` (int): 0=molar, 1=mass basis (usually 0)
- `iFlag` (int): Flag (usually 0)
- `a` (float): First input value
- `b` (float): Second input value
- `z` (list): Composition array (mole fractions)

**Returns:**
- `r.Output[i]` - Array of output values
- `r.ierr` - Error code (0=success, >0=error, <0=warning)
- `r.herr` - Error message string

### TPFLSHdll - Flash Calculation

```python
r = RP.TPFLSHdll(T, P, z)
```

Temperature-pressure flash for phase equilibrium.

**Returns:**
- `r.D` - Total density
- `r.Dl` - Liquid density
- `r.Dv` - Vapor density
- `r.q` - Quality (vapor fraction)
- `r.x` - Liquid composition
- `r.y` - Vapor composition
- `r.e`, `r.h`, `r.s`, `r.cv`, `r.cp`, `r.w` - Properties

### SATTdll - Saturation at Temperature

```python
r = RP.SATTdll(T, z, kph=1)
```

**kph:** 1=bubble point, 2=dew point, 3=both

**Returns:**
- `r.P` - Saturation pressure
- `r.Dl` - Liquid density
- `r.Dv` - Vapor density
- `r.x` - Liquid composition
- `r.y` - Vapor composition

### SATPdll - Saturation at Pressure

```python
r = RP.SATPdll(P, z, kph=1)
```

Similar to SATTdll but at specified pressure.

### CRITPdll - Critical Point

```python
r = RP.CRITPdll(z)
```

**Returns:**
- `r.Tc` - Critical temperature
- `r.Pc` - Critical pressure
- `r.Dc` - Critical density

### INFOdll - Fluid Information

```python
r = RP.INFOdll(icomp)
```

Get information about component `icomp`.

**Returns:**
- `r.wmm` - Molar mass
- `r.Ttrp` - Triple point temperature
- `r.Tnbpt` - Normal boiling point
- `r.Tc` - Critical temperature
- `r.Pc` - Critical pressure
- `r.Dc` - Critical density
- `r.Zc` - Critical compressibility
- `r.acf` - Acentric factor

---

## External Resources

### Official NIST Resources

1. **NIST REFPROP Product Page**
   - URL: https://www.nist.gov/srd/refprop
   - Purchase, documentation, updates

2. **REFPROP Download (License Holders)**
   - URL: https://www.nist.gov/srd/refprop-download
   - Software downloads for registered users

3. **REFPROP FAQ**
   - URL: https://pages.nist.gov/REFPROP-docs/
   - Frequently asked questions

4. **REFPROP Documentation (PDF)**
   - Included with installation
   - Location: `<REFPROP>\REFPROP.pdf`
   - Comprehensive 300+ page manual

### Python Wrapper Resources

5. **ctREFPROP GitHub Repository**
   - URL: https://github.com/usnistgov/REFPROP-wrappers/tree/master/wrappers/python
   - Python wrapper source and examples

6. **ctREFPROP Documentation**
   - URL: https://pages.nist.gov/REFPROP-wrappers/python/
   - API documentation and tutorials

7. **REFPROP-wrappers (All Languages)**
   - URL: https://github.com/usnistgov/REFPROP-wrappers
   - Wrappers for Python, MATLAB, C++, C#, etc.

### Support and Community

8. **REFPROP Google Group**
   - URL: https://groups.google.com/g/refprop
   - User community forum

9. **REFPROP Issue Tracker**
   - URL: https://github.com/usnistgov/REFPROP-issues
   - Report bugs and request features

10. **NIST Technical Support**
    - Email: refprop@nist.gov
    - Official support for license holders

### Key Publications

11. **REFPROP Version 10.0 Reference**
    - Lemmon, E.W., Bell, I.H., Huber, M.L., McLinden, M.O. (2018)
    - "NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties-REFPROP, Version 10.0"
    - National Institute of Standards and Technology, Boulder, CO

12. **GERG-2008 Equation (Natural Gas)**
    - Kunz, O. and Wagner, W. (2012)
    - "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures"
    - J. Chem. Eng. Data, 57(11), 3032-3091

13. **IAPWS-95 (Water)**
    - Wagner, W. and Pruß, A. (2002)
    - "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance"
    - J. Phys. Chem. Ref. Data, 31, 387-535

### Comparison and Validation

14. **CoolProp Validation vs REFPROP**
    - URL: http://www.coolprop.org/validation/index.html
    - Shows deviations between CoolProp and REFPROP

### Related NIST Databases

15. **NIST Chemistry WebBook**
    - URL: https://webbook.nist.gov/chemistry/
    - Free thermochemical data (less comprehensive than REFPROP)

16. **NIST ThermoData Engine (TDE)**
    - URL: https://www.nist.gov/mml/acmd/trc/thermodata-engine
    - Advanced thermodynamic database (enterprise solution)

### Academic and Training

17. **REFPROP Training Materials**
    - Check NIST website for workshops and tutorials
    - Sometimes offered at ASHRAE and AIChE conferences

18. **Example Scripts and Notebooks**
    - URL: https://github.com/usnistgov/REFPROP-wrappers/tree/master/wrappers/python/examples
    - Jupyter notebooks and example Python scripts

### Licensing and Purchase

19. **NIST SRD (Standard Reference Data) Program**
    - URL: https://www.nist.gov/srd
    - Overview of all NIST databases

20. **Academic Licensing**
    - Contact NIST for educational institution pricing
    - Site licenses available

---

## Quick Command Reference

```python
# Initialization
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os

os.environ['RPPREFIX'] = r'C:\Program Files (x86)\REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum

# Pure fluid query
r = RP.REFPROPdll("WATER", "TP", "D;H;S;CP", MOLAR_BASE_SI, 0, 0, 300.0, 101.325, [1.0])

# Mixture query
r = RP.REFPROPdll("METHANE;ETHANE", "TP", "D;H;S", MOLAR_BASE_SI, 0, 0, 200.0, 5000.0, [0.9, 0.1])

# Saturation properties
r = RP.SATTdll(300.0, [1.0], kph=1)  # Bubble point

# Critical point
r = RP.CRITPdll([1.0])

# Error checking
if r.ierr > 0:
    print(f"Error: {r.herr}")
elif r.ierr < 0:
    print(f"Warning: {r.herr}")
```

---

*This reference document is maintained as part of the Claude Engineering Skills Library. REFPROP is a product of NIST and requires a commercial license. For the most up-to-date information, consult the official NIST REFPROP documentation at https://www.nist.gov/srd/refprop*
