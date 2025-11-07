---
name: nist-refprop
description: "Query high-accuracy thermodynamic properties from NIST REFPROP database (commercial)"
category: databases
domain: fluids
complexity: intermediate
dependencies:
  - ctREFPROP
---

# NIST REFPROP Database Skill

Query high-accuracy thermodynamic and transport properties using the NIST REFPROP (Reference Fluid Thermodynamic and Transport Properties) database - the gold standard for thermophysical property calculations.

## Overview

NIST REFPROP is a commercial-grade thermophysical property database developed by the National Institute of Standards and Technology (NIST). It provides:

- **Pure Fluids**: 147+ pure fluids with the highest available accuracy
- **Mixtures**: Predefined and custom mixtures with advanced mixing rules
- **Refrigerant Blends**: All common refrigerant blends with temperature glide
- **Transport Properties**: Viscosity, thermal conductivity, surface tension
- **Thermodynamic Properties**: Enthalpy, entropy, density, specific heat
- **Phase Equilibrium**: VLE, LLE, VLLE calculations for mixtures
- **Highest Accuracy**: Reference-quality data validated against experimental measurements

REFPROP is the industry and research standard for:
- Critical engineering design requiring highest accuracy
- Validation of other property libraries (CoolProp is validated against REFPROP)
- Mixture calculations with complex phase behavior
- Refrigerant blend applications with temperature glide
- Research and publication-quality data

## Licensing

**COMMERCIAL LICENSE REQUIRED**

NIST REFPROP is **NOT** free or open-source software. You must purchase a license from NIST:

- **Standard License**: ~$300 USD (one-time purchase)
- **Academic Pricing**: Available for educational institutions
- **Purchase**: https://www.nist.gov/srd/refprop

**License includes:**
- REFPROP software (Windows GUI)
- Fluid database files
- Excel add-in
- FORTRAN/DLL libraries
- Python wrapper compatibility

**Important Notes:**
- License is per-user, perpetual (one-time fee)
- Updates and new fluid files provided for ~5 years
- Required for commercial applications
- Cannot be redistributed

## Installation

### Step 1: Purchase and Install REFPROP

1. Purchase REFPROP license from NIST: https://www.nist.gov/srd/refprop
2. Download and install REFPROP (typically to `C:\Program Files (x86)\REFPROP\`)
3. Note the installation directory path

### Step 2: Install Python Wrapper

```bash
pip install ctREFPROP
```

The `ctREFPROP` package provides a Python interface to REFPROP using ctypes.

### Step 3: Set RPPREFIX Environment Variable

The Python wrapper needs to know where REFPROP is installed:

**Windows:**
```bash
# PowerShell
$env:RPPREFIX = "C:\Program Files (x86)\REFPROP\"

# Or set permanently via System Properties > Environment Variables
```

**Linux/Mac:**
```bash
export RPPREFIX="/path/to/REFPROP/"
```

**Python (alternative):**
```python
import os
os.environ['RPPREFIX'] = r'C:\Program Files (x86)\REFPROP'
```

### Step 4: Verify Installation

```python
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# Initialize REFPROP
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])

# Test query
print(RP.RPVersion())
```

## Available Fluids

REFPROP includes **147+ pure fluids** organized by class:

### Common Industrial Fluids
- **Water** (H₂O) - IAPWS-95 formulation
- **Air** - Dry air mixture
- **Carbon Dioxide** (CO₂)
- **Nitrogen** (N₂)
- **Oxygen** (O₂)
- **Argon** (Ar)
- **Helium** (He)
- **Hydrogen** (H₂)

### Hydrocarbons (40+ fluids)
- **Alkanes**: Methane, Ethane, Propane, Butane, Pentane, Hexane, Heptane, Octane, Nonane, Decane (and isomers)
- **Alkenes**: Ethylene, Propylene, Butene
- **Aromatics**: Benzene, Toluene, Xylenes
- **Cyclic**: Cyclohexane, Cyclopentane

### Refrigerants (60+ fluids)
- **CFCs**: R11, R12, R13, R113, R114, R115
- **HCFCs**: R22, R123, R124, R141b, R142b
- **HFCs**: R23, R32, R125, R134a, R143a, R152a, R227ea, R236fa, R245fa, R365mfc
- **HFOs**: R1234yf, R1234ze(E), R1233zd(E), R1336mzz(Z)
- **Natural**: R717 (Ammonia), R744 (CO₂), R290 (Propane), R600a (Isobutane)
- **Blends**: R404A, R407C, R410A, R507A, R407F, R448A, R449A, R450A, R513A

### Cryogenic Fluids
- Helium, Neon, Argon, Krypton, Xenon
- Hydrogen, Deuterium, Parahydrogen
- Nitrogen, Oxygen, Fluorine
- Methane, Ethane (LNG components)

### Specialty Fluids
- **Siloxanes**: D4, D5, D6, MD2M, MD3M, MD4M, MDM, MM (ORC applications)
- **Alcohols**: Methanol, Ethanol, Propanol, Butanol
- **Other Organics**: Acetone, Benzene, Toluene
- **Inorganics**: Ammonia, Sulfur dioxide, Hydrogen sulfide

See `reference.md` for complete fluid list.

## Property Queries

### Core Functions

REFPROP uses different function calls than CoolProp:

```python
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os

# Initialize
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])

# Define fluid
MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum
fluid = "WATER"

# Single component
r = RP.REFPROPdll(fluid, "TP", "D;H;S;CP;CV;W;VIS;TCX", MOLAR_BASE_SI, 0, 0, 298.15, 101.325, [1.0])

# Access results
density = r.Output[0]      # kg/m³
enthalpy = r.Output[1]     # J/kg
entropy = r.Output[2]      # J/kg/K
cp = r.Output[3]           # J/kg/K
cv = r.Output[4]           # J/kg/K
sound_speed = r.Output[5]  # m/s
viscosity = r.Output[6]    # Pa·s
conductivity = r.Output[7] # W/m/K
```

### Input Specifications

| Code | Property | Unit |
|------|----------|------|
| `TP` | Temperature, Pressure | K, kPa |
| `TH` | Temperature, Enthalpy | K, J/kg |
| `TS` | Temperature, Entropy | K, J/kg/K |
| `TD` | Temperature, Density | K, kg/m³ |
| `PH` | Pressure, Enthalpy | kPa, J/kg |
| `PS` | Pressure, Entropy | kPa, J/kg/K |
| `PD` | Pressure, Density | kPa, kg/m³ |
| `TQ` | Temperature, Quality | K, - |
| `PQ` | Pressure, Quality | kPa, - |

### Output Properties

| Code | Property | Unit |
|------|----------|------|
| `D` | Density | kg/m³ |
| `H` | Enthalpy | J/kg |
| `S` | Entropy | J/kg/K |
| `U` | Internal energy | J/kg |
| `CP` | Heat capacity (const P) | J/kg/K |
| `CV` | Heat capacity (const V) | J/kg/K |
| `W` | Speed of sound | m/s |
| `VIS` | Viscosity | Pa·s |
| `TCX` | Thermal conductivity | W/m/K |
| `KV` | Kinematic viscosity | m²/s |
| `PRANDTL` | Prandtl number | - |
| `STN` | Surface tension | N/m |
| `P` | Pressure | kPa |
| `T` | Temperature | K |
| `Q` | Quality (vapor fraction) | - |

Multiple properties can be requested in one call, separated by semicolons: `"D;H;S;CP;VIS;TCX"`

## Mixture Calculations

One of REFPROP's key advantages is sophisticated mixture handling:

```python
# Define mixture by components and mole fractions
fluid = "METHANE;ETHANE;PROPANE"
z = [0.9, 0.07, 0.03]  # Mole fractions (must sum to 1.0)

# Query mixture properties
r = RP.REFPROPdll(fluid, "TP", "D;H;S", MOLAR_BASE_SI, 0, 0, 250.0, 5000.0, z)
```

**Mixture Features:**
- Up to 20 components in a mixture
- Advanced mixing rules (GERG-2008, Kunz-Wagner, etc.)
- Temperature glide calculations for zeotropic blends
- Vapor-liquid equilibrium (VLE) calculations
- Liquid-liquid equilibrium (LLE) for certain systems

## Advantages Over CoolProp

While CoolProp is excellent and free, REFPROP offers:

### 1. Higher Accuracy
- REFPROP is the **reference standard** that other libraries validate against
- More precise equations of state, especially near critical point
- Better transport property correlations
- Typical accuracy: 0.01% for density vs 0.1% for CoolProp

### 2. More Fluids
- 147+ pure fluids vs ~120 in CoolProp
- All latest refrigerants (HFOs: R1234yf, R1234ze, R1233zd, etc.)
- More industrial and specialty chemicals
- Frequently updated with new fluids

### 3. Better Mixture Handling
- Superior mixing rules (GERG-2008, advanced models)
- Handles complex non-ideal mixtures
- Temperature glide for zeotropic refrigerant blends
- VLE/LLE calculations
- Up to 20 components vs limited mixture support in CoolProp

### 4. Extended Property Range
- Valid over wider temperature and pressure ranges
- Better extrapolation behavior
- More reliable in extreme conditions

### 5. Additional Capabilities
- Humidity calculations (psychrometrics)
- Flash calculations (PT, PH, PS, etc.)
- Saturation properties with Jacobian
- Ancillary equations
- Pseudo-pure fluid approximations

### 6. Official Support
- Maintained by NIST (National Institute of Standards and Technology)
- Technical support available
- Regular updates with latest research
- Publication-quality results

### When to Use REFPROP vs CoolProp

**Use REFPROP when:**
- Highest accuracy is required (±0.01%)
- Working with complex mixtures
- Need latest refrigerant data (HFOs)
- Publication or critical design work
- Budget allows commercial license
- Working with zeotropic blends (temperature glide)

**Use CoolProp when:**
- Budget constraints (free/open-source)
- Basic property needs (±0.1% accuracy acceptable)
- Pure fluids only
- Open-source license required
- Quick prototyping

**Validation Note:** CoolProp is validated against REFPROP, confirming REFPROP as the reference standard.

## Common Workflow Example

```python
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os

# Setup
os.environ['RPPREFIX'] = r'C:\Program Files (x86)\REFPROP'
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])

MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum

# Query water at 25°C, 1 bar
fluid = "WATER"
T = 298.15  # K
P = 100.0   # kPa

r = RP.REFPROPdll(fluid, "TP", "D;H;S;CP;VIS;TCX", MOLAR_BASE_SI, 0, 0, T, P, [1.0])

print(f"Water at {T} K, {P} kPa:")
print(f"  Density: {r.Output[0]:.4f} kg/m³")
print(f"  Enthalpy: {r.Output[1]:.2f} J/kg")
print(f"  Entropy: {r.Output[2]:.2f} J/kg/K")
print(f"  Cp: {r.Output[3]:.2f} J/kg/K")
print(f"  Viscosity: {r.Output[4]*1000:.4f} mPa·s")
print(f"  Conductivity: {r.Output[5]:.4f} W/m/K")
```

## Error Handling

REFPROP returns error codes and messages:

```python
r = RP.REFPROPdll(fluid, "TP", "D", MOLAR_BASE_SI, 0, 0, T, P, z)

if r.ierr > 0:
    print(f"Error {r.ierr}: {r.herr}")
elif r.ierr < 0:
    print(f"Warning {r.ierr}: {r.herr}")
else:
    print(f"Success: {r.Output[0]}")
```

## Units

REFPROP uses different unit conventions than CoolProp:

**Default Units (with MOLAR BASE SI):**
- Temperature: K (Kelvin)
- Pressure: kPa (kilopascals) - **Note: CoolProp uses Pa**
- Density: kg/m³
- Enthalpy: J/kg
- Entropy: J/kg/K
- Viscosity: Pa·s
- Thermal conductivity: W/m/K

**Unit Systems Available:**
- SI with molar units
- SI with mass units (default shown above)
- English units
- Custom units

## Best Practices

1. **Always set RPPREFIX** before importing ctREFPROP
2. **Check ierr** for errors after each call
3. **Use exact fluid names** (case-sensitive, from FLUIDS directory)
4. **Pre-calculate compositions** for mixtures (must sum to 1.0)
5. **Cache RP instance** - don't reinitialize repeatedly
6. **Use appropriate input pairs** for each phase region
7. **Validate results** against known values initially
8. **Read error messages** - REFPROP provides detailed diagnostics

## References

### Official Resources

1. **NIST REFPROP Website**
   - URL: https://www.nist.gov/srd/refprop
   - Purchase and download

2. **REFPROP Documentation**
   - Included with installation
   - `REFPROP.pdf` in installation directory

3. **ctREFPROP Documentation**
   - URL: https://github.com/usnistgov/REFPROP-wrappers/tree/master/wrappers/python
   - GitHub: https://github.com/usnistgov/REFPROP-wrappers

### Key Publications

- Lemmon, E.W., Bell, I.H., Huber, M.L., McLinden, M.O. (2018). "NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties-REFPROP, Version 10.0", National Institute of Standards and Technology.

### Support

- **NIST Support**: refprop@nist.gov
- **Issue Tracker**: https://github.com/usnistgov/REFPROP-issues
- **Forum**: https://groups.google.com/g/refprop

## Quick Reference

```python
# Initialize
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])
MOLAR_BASE_SI = RP.GETENUMdll(0, "MOLAR BASE SI").iEnum

# Pure fluid
r = RP.REFPROPdll("WATER", "TP", "D;H;S", MOLAR_BASE_SI, 0, 0, 300.0, 101.325, [1.0])

# Mixture
r = RP.REFPROPdll("METHANE;ETHANE", "TP", "D;H;S", MOLAR_BASE_SI, 0, 0, 200.0, 5000.0, [0.9, 0.1])

# Check for errors
if r.ierr != 0:
    print(f"Error/Warning: {r.herr}")
```

---

*NIST REFPROP is the gold standard for thermophysical properties. While it requires a commercial license, the investment is justified for applications requiring the highest accuracy, mixture calculations, or publication-quality results. For open-source alternatives, see the `coolprop-db` skill.*
