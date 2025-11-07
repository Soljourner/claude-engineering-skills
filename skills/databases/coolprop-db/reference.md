# CoolProp Database Reference

Comprehensive reference for CoolProp thermophysical property database queries.

## Table of Contents

1. [Supported Fluids](#supported-fluids)
2. [Property Parameters](#property-parameters)
3. [Input Parameter Pairs](#input-parameter-pairs)
4. [Equation of State Information](#equation-of-state-information)
5. [Unit Conversions](#unit-conversions)
6. [Phase Identification](#phase-identification)
7. [API Function Reference](#api-function-reference)
8. [External Resources](#external-resources)

---

## Supported Fluids

### Pure Fluids (100+ Fluids)

CoolProp includes high-accuracy equations of state for the following pure fluids:

#### Water and Air
- `Water` - H₂O (IAPWS-95 formulation)
- `Air` - Dry air (Lemmon et al. formulation, pseudo-pure)

#### Noble Gases
- `Helium` - He
- `Neon` - Ne
- `Argon` - Ar
- `Krypton` - Kr
- `Xenon` - Xe

#### Light Gases
- `Hydrogen` - H₂
- `Deuterium` - D₂ (Heavy hydrogen)
- `ParaHydrogen` - Para-H₂
- `OrthoHydrogen` - Ortho-H₂
- `Nitrogen` - N₂
- `Oxygen` - O₂

#### Carbon Compounds
- `CarbonMonoxide` - CO
- `CarbonDioxide` / `CO2` - CO₂
- `CarbonylSulfide` - COS

#### Hydrocarbons (Alkanes)
- `Methane` - CH₄
- `Ethane` - C₂H₆
- `Propane` - C₃H₈
- `n-Butane` - C₄H₁₀
- `IsoButane` - i-C₄H₁₀
- `n-Pentane` - C₅H₁₂
- `Isopentane` - i-C₅H₁₂
- `n-Hexane` - C₆H₁₄
- `Isohexane` - i-C₆H₁₄
- `n-Heptane` - C₇H₁₆
- `n-Octane` - C₈H₁₈
- `n-Nonane` - C₉H₂₀
- `n-Decane` - C₁₀H₂₂
- `n-Undecane` - C₁₁H₂₄
- `n-Dodecane` - C₁₂H₂₆

#### Alkenes
- `Ethylene` - C₂H₄
- `Propylene` / `Propene` - C₃H₆
- `1-Butene` - C₄H₈
- `IsoButene` - i-C₄H₈
- `cis-2-Butene` - C₄H₈
- `trans-2-Butene` - C₄H₈

#### Cyclic Hydrocarbons
- `CycloHexane` - C₆H₁₂
- `CycloPentane` - C₅H₁₀
- `CycloPropane` - C₃H₆
- `Benzene` - C₆H₆
- `Toluene` - C₇H₈
- `m-Xylene` - C₈H₁₀
- `o-Xylene` - C₈H₁₀
- `p-Xylene` - C₈H₁₀
- `EthylBenzene` - C₈H₁₀

#### Alcohols
- `Methanol` - CH₃OH
- `Ethanol` - C₂H₅OH
- `1-Propanol` / `n-Propanol` - C₃H₇OH
- `2-Propanol` / `Isopropanol` - C₃H₇OH
- `1-Butanol` - C₄H₉OH

#### Refrigerants (HFCs - Hydrofluorocarbons)
- `R11` - CFCl₃ (Trichlorofluoromethane)
- `R12` - CF₂Cl₂ (Dichlorodifluoromethane)
- `R13` - CClF₃ (Chlorotrifluoromethane)
- `R14` - CF₄ (Tetrafluoromethane)
- `R22` - CHClF₂ (Chlorodifluoromethane, HCFC)
- `R23` - CHF₃ (Trifluoromethane)
- `R32` - CH₂F₂ (Difluoromethane)
- `R41` - CH₃F (Fluoromethane)
- `R113` - C₂Cl₃F₃
- `R114` - C₂Cl₂F₄
- `R115` - C₂ClF₅
- `R116` - C₂F₆ (Hexafluoroethane)
- `R123` - CHCl₂CF₃
- `R124` - CHClFCF₃
- `R125` - CHF₂CF₃ (Pentafluoroethane)
- `R134a` - CH₂FCF₃ (1,1,1,2-Tetrafluoroethane) **Most common automotive/residential refrigerant**
- `R141b` - CH₃CCl₂F
- `R142b` - CH₃CClF₂
- `R143a` - CH₃CF₃
- `R152A` - CH₃CHF₂
- `R161` - CH₃CH₂F
- `R21` - CHCl₂F
- `R227EA` - CF₃CHFCF₃
- `R236EA` - CHF₂CHFCF₃
- `R236FA` - CF₃CH₂CF₃
- `R245ca` - CH₂FCF₂CHF₂
- `R245fa` - CHF₂CH₂CF₃
- `RC318` - C₄F₈ (Octafluorocyclobutane)

#### Refrigerant Blends (Pseudo-pure)
- `R404A` - R125/143a/134a (44/52/4%) - Commercial refrigeration
- `R407C` - R32/125/134a (23/25/52%) - Air conditioning
- `R410A` - R32/125 (50/50%) - Residential AC, heat pumps
- `R507A` - R125/143a (50/50%) - Low-temperature refrigeration
- `R508A` - R23/116 (39/61%) - Ultra-low temperature
- `R508B` - R23/116 (46/54%) - Ultra-low temperature

#### Natural Refrigerants
- `R717` - Ammonia (NH₃) - Industrial refrigeration
- `R744` - Carbon dioxide (CO₂) - Transcritical systems
- `R290` - Propane - Hydrocarbon refrigerant
- `R600` - n-Butane - Hydrocarbon refrigerant
- `R600a` - Isobutane - Hydrocarbon refrigerant
- `R1270` - Propylene/Propene - Hydrocarbon refrigerant

#### Sulfur and Halogen Compounds
- `SulfurDioxide` / `SO2` - SO₂
- `SulfurHexafluoride` - SF₆
- `HydrogenSulfide` / `H2S` - H₂S

#### Siloxanes
- `D4` - Octamethylcyclotetrasiloxane
- `D5` - Decamethylcyclopentasiloxane
- `D6` - Dodecamethylcyclohexasiloxane
- `MD2M` - Decamethyltetrasiloxane
- `MD3M` - Dodecamethylpentasiloxane
- `MD4M` - Tetradecamethylhexasiloxane
- `MDM` - Octamethyltrisiloxane
- `MM` - Hexamethyldisiloxane

#### Other Organic Compounds
- `Acetone` - CH₃COCH₃
- `DiethylEther` - C₄H₁₀O
- `DimethylEther` / `RE170` - CH₃OCH₃
- `Ethylene` - C₂H₄
- `EthyleneOxide` - C₂H₄O

#### Ammonia and Derivatives
- `Ammonia` - NH₃

---

## Property Parameters

### Complete List of Property Codes

All property codes for use with `PropsSI()` and `Props1SI()`.

#### Thermodynamic Properties (Two-input functions)

| Parameter | Description | SI Unit | Symbol |
|-----------|-------------|---------|--------|
| `T` | Temperature | K | T |
| `P` | Pressure | Pa | P |
| `D` | Mass density | kg/m³ | ρ |
| `Dmass` | Mass density | kg/m³ | ρ |
| `Dmolar` | Molar density | mol/m³ | ρₘ |
| `H` | Mass specific enthalpy | J/kg | h |
| `Hmass` | Mass specific enthalpy | J/kg | h |
| `Hmolar` | Molar specific enthalpy | J/mol | hₘ |
| `S` | Mass specific entropy | J/kg/K | s |
| `Smass` | Mass specific entropy | J/kg/K | s |
| `Smolar` | Molar specific entropy | J/mol/K | sₘ |
| `U` | Mass specific internal energy | J/kg | u |
| `Umass` | Mass specific internal energy | J/kg | u |
| `Umolar` | Molar specific internal energy | J/mol | uₘ |
| `G` | Mass specific Gibbs energy | J/kg | g |
| `Gmass` | Mass specific Gibbs energy | J/kg | g |
| `Gmolar` | Molar specific Gibbs energy | J/mol | gₘ |
| `Helmholtzmass` | Mass specific Helmholtz energy | J/kg | a |
| `Helmholtzmolar` | Molar specific Helmholtz energy | J/mol | aₘ |
| `Q` | Vapor quality (mass fraction) | - | x |
| `Phase` | Phase index | - | - |

#### Transport Properties

| Parameter | Description | SI Unit | Symbol |
|-----------|-------------|---------|--------|
| `V` | Viscosity (dynamic) | Pa·s | μ |
| `viscosity` | Viscosity (dynamic) | Pa·s | μ |
| `L` | Thermal conductivity | W/m/K | k, λ |
| `conductivity` | Thermal conductivity | W/m/K | k, λ |
| `I` | Surface tension | N/m | σ |
| `surface_tension` | Surface tension | N/m | σ |
| `Prandtl` | Prandtl number | - | Pr |

#### Specific Heats

| Parameter | Description | SI Unit | Symbol |
|-----------|-------------|---------|--------|
| `C` | Mass specific constant pressure specific heat | J/kg/K | cₚ |
| `Cpmass` | Mass specific constant pressure specific heat | J/kg/K | cₚ |
| `Cpmolar` | Molar specific constant pressure specific heat | J/mol/K | cₚₘ |
| `O` | Mass specific constant volume specific heat | J/kg/K | cᵥ |
| `Cvmass` | Mass specific constant volume specific heat | J/kg/K | cᵥ |
| `Cvmolar` | Molar specific constant volume specific heat | J/mol/K | cᵥₘ |
| `Cp0mass` | Ideal gas mass specific constant pressure specific heat | J/kg/K | cₚ⁰ |
| `Cp0molar` | Ideal gas molar specific constant pressure specific heat | J/mol/K | cₚₘ⁰ |

#### Derivative Properties

| Parameter | Description | SI Unit | Symbol |
|-----------|-------------|---------|--------|
| `isobaric_expansion_coefficient` | Isobaric expansion coefficient | 1/K | β |
| `isothermal_compressibility` | Isothermal compressibility | 1/Pa | κₜ |
| `fundamental_derivative_of_gas_dynamics` | Fundamental derivative | - | Γ |
| `speed_of_sound` | Speed of sound | m/s | w |

#### Saturation Properties (Use with T or P and Q)

| Parameter | Description | SI Unit | Notes |
|-----------|-------------|---------|-------|
| `T` with `Q` | Saturation temperature | K | Given pressure |
| `P` with `Q` | Saturation pressure | Pa | Given temperature |

#### Single-Input Properties (Use with Props1SI)

| Parameter | Description | SI Unit | Symbol |
|-----------|-------------|---------|--------|
| `Tcrit` | Critical temperature | K | Tс |
| `Pcrit` | Critical pressure | Pa | Pс |
| `Dcrit` | Critical density | kg/m³ | ρс |
| `Ttriple` | Triple point temperature | K | Tₜ |
| `Ptriple` | Triple point pressure | Pa | Pₜ |
| `Tmin` | Minimum temperature | K | - |
| `Tmax` | Maximum temperature | K | - |
| `Pmin` | Minimum pressure | Pa | - |
| `Pmax` | Maximum pressure | Pa | - |
| `T_reducing` | Reducing temperature | K | - |
| `P_reducing` | Reducing pressure | Pa | - |
| `rhomass_reducing` | Reducing density | kg/m³ | - |
| `M` | Molar mass | kg/mol | M |
| `molemass` | Molar mass | kg/mol | M |
| `acentric` | Acentric factor | - | ω |
| `accentric_factor` | Acentric factor | - | ω |
| `dipole_moment` | Dipole moment | Debye | - |
| `GWP20` | 20-year Global Warming Potential | - | - |
| `GWP100` | 100-year Global Warming Potential | - | - |
| `GWP500` | 500-year Global Warming Potential | - | - |
| `ODP` | Ozone Depletion Potential | - | - |
| `FH` | Flame Hazard | - | - |
| `HH` | Health Hazard | - | - |
| `PH` | Physical Hazard | - | - |

---

## Input Parameter Pairs

Valid combinations of input parameters for `PropsSI()`.

### Valid Input Pairs by Phase Region

| Input Pair | Subcooled Liquid | Two-Phase | Superheated Vapor | Supercritical | Notes |
|------------|------------------|-----------|-------------------|---------------|-------|
| `(T, P)` | ✓ | ✗ | ✓ | ✓ | Not valid in two-phase region |
| `(P, Q)` | ✓* | ✓ | ✓* | ✗ | Q=0 (liquid) or Q=1 (vapor) at boundaries |
| `(T, Q)` | ✓* | ✓ | ✓* | ✗ | Q=0 (liquid) or Q=1 (vapor) at boundaries |
| `(P, H)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |
| `(P, S)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |
| `(P, D)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |
| `(P, U)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |
| `(H, S)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere, useful for isentropic |
| `(D, H)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |
| `(D, S)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |
| `(D, U)` | ✓ | ✓ | ✓ | ✓ | Valid everywhere |

**Legend:**
- ✓ = Valid
- ✗ = Invalid
- ✓* = Valid only at saturation boundaries (Q=0 or Q=1)

### Recommended Input Pairs

| Application | Recommended Pair | Reason |
|-------------|------------------|--------|
| Known temperature and pressure | `(T, P)` | Direct physical measurements |
| Saturation properties | `(T, Q)` or `(P, Q)` | Explicitly specifies phase |
| Compressor/pump outlet | `(P, H)` or `(P, S)` | Pressure set by system, enthalpy/entropy from energy balance |
| Isentropic processes | `(P, S)` or `(H, S)` | Entropy remains constant |
| Heat exchangers | `(P, H)` | Pressure known, enthalpy from energy balance |
| Two-phase calculations | `(P, Q)` or `(T, Q)` | Quality explicitly defined |
| General robust queries | `(P, H)` or `(P, S)` | Valid in all phase regions |

---

## Equation of State Information

CoolProp uses highly accurate multi-parameter equations of state.

### Helmholtz Energy Formulation

Most pure fluids use explicit Helmholtz energy equations of state:

```
a(δ, τ) = a⁰(δ, τ) + aʳ(δ, τ)
```

Where:
- `a` = Helmholtz energy / (RT)
- `δ` = ρ/ρc (reduced density)
- `τ` = Tc/T (inverse reduced temperature)
- `a⁰` = Ideal gas contribution
- `aʳ` = Residual (real gas) contribution

All thermodynamic properties are derived from derivatives of this function.

### Equation Sources by Fluid Type

#### Water
- **IAPWS-95**: International Association for the Properties of Water and Steam
- **Accuracy**: ±0.001% for density, ±0.01% for speed of sound
- **Reference**: Wagner, W. and Pruß, A. (2002). J. Phys. Chem. Ref. Data, 31, 387-535

#### Refrigerants
Most refrigerants use equations from:
- **Tillner-Roth**: R134a, R32, etc.
- **Lemmon-Span**: Short equations for light refrigerants
- **NIST REFPROP**: High-accuracy multiparameter equations

#### Hydrocarbons
- **Natural Gas**: Kunz-Wagner (GERG-2008) for mixtures
- **Pure Components**: Span, Lemmon, and other NIST researchers

#### Industrial Gases
- **Nitrogen, Oxygen, Argon**: Span et al. reference equations
- **Carbon Dioxide**: Span-Wagner equation
- **Ammonia**: Tillner-Roth equation

### Accuracy

Typical uncertainties in the equations of state:

| Property | Liquid Phase | Vapor Phase | Critical Region |
|----------|--------------|-------------|-----------------|
| Density | ±0.01-0.1% | ±0.05-0.2% | ±0.5-1% |
| Vapor Pressure | ±0.02-0.1% | - | ±0.5% |
| Speed of Sound | ±0.02-0.5% | ±0.1-1% | ±1-5% |
| Heat Capacity | ±0.5-2% | ±0.5-2% | ±5-10% |

Transport properties (viscosity, thermal conductivity) typically have:
- **Viscosity**: ±1-5% uncertainty
- **Thermal Conductivity**: ±2-5% uncertainty
- **Surface Tension**: ±0.5-2% uncertainty

---

## Unit Conversions

CoolProp uses **SI units exclusively**. Common conversions:

### Temperature
```
T[K] = T[°C] + 273.15
T[K] = (T[°F] + 459.67) × 5/9
T[K] = T[°R] × 5/9
```

### Pressure
```
P[Pa] = P[bar] × 1e5
P[Pa] = P[kPa] × 1000
P[Pa] = P[MPa] × 1e6
P[Pa] = P[psi] × 6894.76
P[Pa] = P[atm] × 101325
P[Pa] = P[mmHg] × 133.322
```

### Enthalpy / Internal Energy
```
h[J/kg] = h[kJ/kg] × 1000
h[J/kg] = h[BTU/lb] × 2326
```

### Entropy
```
s[J/kg/K] = s[kJ/kg/K] × 1000
s[J/kg/K] = s[BTU/lb/°R] × 4186.8
```

### Viscosity
```
μ[Pa·s] = μ[cP] / 1000  (1 cP = 1 mPa·s)
μ[Pa·s] = μ[lbm/ft/s] × 1.488164
```

### Density
```
ρ[kg/m³] = ρ[g/cm³] × 1000
ρ[kg/m³] = ρ[lb/ft³] × 16.0185
```

---

## Phase Identification

### Phase Codes

The `Phase` parameter returns an integer code:

| Code | Phase | Description |
|------|-------|-------------|
| 0 | Liquid | Subcooled liquid or compressed liquid |
| 3 | Supercritical liquid | T > Tc, ρ > ρc |
| 4 | Supercritical gas | T > Tc, ρ < ρc |
| 5 | Supercritical | T > Tc and P > Pc |
| 6 | Gas | Superheated vapor or low-pressure gas |
| 8 | Not imposed | Two-phase region |

### Quality (Vapor Fraction)

Quality `Q` is defined as:

```
Q = m_vapor / m_total
```

- **Q = 0**: Saturated liquid (bubble point)
- **0 < Q < 1**: Two-phase mixture
- **Q = 1**: Saturated vapor (dew point)
- **Q < 0 or Q > 1**: Subcooled liquid or superheated vapor (undefined)

### Phase Diagram Regions

```
        Supercritical
             │
    Liquid   │   Vapor
       ╲     │     ╱
        ╲    │    ╱
         ╲   │   ╱
          ╲  │  ╱
           ╲ │ ╱
            ╲│╱
        Critical Point
             │
    ─────────┼─────────  Saturation Curve
             │
        Two-Phase
```

---

## API Function Reference

### PropsSI - Primary Function

```python
PropsSI(output, input1_name, input1_value, input2_name, input2_value, fluid)
```

**Parameters:**
- `output` (str): Property to calculate (e.g., 'T', 'P', 'D', 'H', 'S')
- `input1_name` (str): First input property name
- `input1_value` (float): First input property value (SI units)
- `input2_name` (str): Second input property name
- `input2_value` (float): Second input property value (SI units)
- `fluid` (str): Fluid name

**Returns:** float - Calculated property value in SI units

**Example:**
```python
density = PropsSI('D', 'T', 298.15, 'P', 101325, 'Water')
```

### Props1SI - Single Input Function

```python
Props1SI(output, fluid)
```

**Parameters:**
- `output` (str): Property to retrieve (e.g., 'Tcrit', 'Pcrit', 'M')
- `fluid` (str): Fluid name

**Returns:** float - Property value in SI units

**Example:**
```python
T_crit = Props1SI('Tcrit', 'Water')
```

### PhaseSI - Phase Determination

```python
PhaseSI(input1_name, input1_value, input2_name, input2_value, fluid)
```

**Returns:** int - Phase code

**Example:**
```python
from CoolProp.CoolProp import PhaseSI
phase = PhaseSI('T', 300, 'P', 101325, 'Water')
```

### get_global_param_string - Global Parameters

```python
get_global_param_string(param)
```

Retrieve global information about CoolProp.

**Common Parameters:**
- `'version'` - CoolProp version
- `'gitrevision'` - Git commit hash
- `'fluids_list'` - Comma-separated list of all fluids

**Example:**
```python
from CoolProp.CoolProp import get_global_param_string
version = get_global_param_string('version')
fluids = get_global_param_string('fluids_list').split(',')
```

### get_fluid_param_string - Fluid-Specific Information

```python
get_fluid_param_string(fluid, param)
```

**Common Parameters:**
- `'CAS'` - CAS registry number
- `'formula'` - Chemical formula
- `'ASHRAE34'` - ASHRAE safety classification
- `'REFPROPname'` - REFPROP fluid name

**Example:**
```python
from CoolProp.CoolProp import get_fluid_param_string
cas = get_fluid_param_string('R134a', 'CAS')
formula = get_fluid_param_string('R134a', 'formula')
```

---

## External Resources

### Official Documentation

1. **CoolProp Main Website**
   - URL: http://www.coolprop.org/
   - Description: Official documentation homepage

2. **Python Documentation**
   - URL: http://www.coolprop.org/coolprop/HighLevelAPI.html
   - Description: High-level API documentation for Python

3. **Fluid Properties Database**
   - URL: http://www.coolprop.org/fluid_properties/PurePseudoPure.html
   - Description: Complete list of fluids with properties

4. **Validation Tables**
   - URL: http://www.coolprop.org/validation/index.html
   - Description: Validation data against NIST REFPROP

5. **Development Documentation**
   - URL: http://www.coolprop.org/dev/index.html
   - Description: Developer documentation and internals

### Source Code and Support

6. **GitHub Repository**
   - URL: https://github.com/CoolProp/CoolProp
   - Description: Source code, issues, pull requests

7. **Issue Tracker**
   - URL: https://github.com/CoolProp/CoolProp/issues
   - Description: Bug reports and feature requests

8. **Forum/Discussions**
   - URL: https://github.com/CoolProp/CoolProp/discussions
   - Description: Community questions and discussions

### Key Publications

9. **Primary Reference Paper**
   - Bell, I. H., Wronski, J., Quoilin, S., & Lemort, V. (2014)
   - "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp"
   - *Industrial & Engineering Chemistry Research*, 53(6), 2498-2508
   - DOI: 10.1021/ie4033999

10. **NIST REFPROP Reference**
    - Lemmon, E. W., Bell, I. H., Huber, M. L., & McLinden, M. O. (2018)
    - NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties (REFPROP), Version 10.0
    - National Institute of Standards and Technology

### Related Standards

11. **IAPWS (Water Properties)**
    - URL: http://www.iapws.org/
    - Description: International Association for the Properties of Water and Steam

12. **ASHRAE Refrigerants**
    - URL: https://www.ashrae.org/
    - Description: American Society of Heating, Refrigerating and Air-Conditioning Engineers

13. **ISO Standards**
    - ISO 17584: Refrigerant properties

### Tutorials and Examples

14. **Tutorial Videos**
    - URL: http://www.coolprop.org/coolprop/examples.html
    - Description: Example code in multiple languages

15. **Jupyter Notebooks**
    - URL: https://github.com/CoolProp/CoolProp-museum
    - Description: Interactive examples and demonstrations

### Python Package

16. **PyPI Package**
    - URL: https://pypi.org/project/CoolProp/
    - Description: Python package index page
    - Install: `pip install CoolProp`

17. **Conda Package**
    - URL: https://anaconda.org/conda-forge/coolprop
    - Description: Conda-forge package
    - Install: `conda install -c conda-forge coolprop`

### Comparison Tools

18. **Online Calculator**
    - URL: http://www.coolprop.org/coolprop/Fluid Properties.html
    - Description: Web-based property calculator

19. **REFPROP Comparison**
    - URL: http://www.coolprop.org/fluid_properties/Validation.html
    - Description: Validation against NIST REFPROP data

### Academic and Research

20. **Citing CoolProp**
    ```bibtex
    @article{CoolProp2014,
        title = {Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp},
        author = {Bell, Ian H. and Wronski, Jorrit and Quoilin, Sylvain and Lemort, Vincent},
        journal = {Industrial \& Engineering Chemistry Research},
        volume = {53},
        number = {6},
        pages = {2498--2508},
        year = {2014},
        doi = {10.1021/ie4033999}
    }
    ```

### Additional Language Bindings

21. **MATLAB/Octave**
    - URL: http://www.coolprop.org/coolprop/wrappers/MATLAB/index.html

22. **C++**
    - URL: http://www.coolprop.org/coolprop/wrappers/CPlusPlus/index.html

23. **Excel Add-in**
    - URL: http://www.coolprop.org/coolprop/wrappers/Excel/index.html

24. **LabVIEW**
    - URL: http://www.coolprop.org/coolprop/wrappers/LabVIEW/index.html

---

## Quick Reference Command Summary

```python
# Import
from CoolProp.CoolProp import PropsSI, Props1SI

# Two-input property query
value = PropsSI('Output', 'Input1', value1, 'Input2', value2, 'Fluid')

# Single-input property query
value = Props1SI('Property', 'Fluid')

# Examples
T = PropsSI('T', 'P', 101325, 'Q', 0, 'Water')  # Saturation temp at 1 atm
rho = PropsSI('D', 'T', 300, 'P', 101325, 'Air')  # Air density
T_crit = Props1SI('Tcrit', 'CO2')  # Critical temperature of CO2
```

---

*This reference document is maintained as part of the Claude Engineering Skills Library. For the most up-to-date information, always consult the official CoolProp documentation at http://www.coolprop.org/*
