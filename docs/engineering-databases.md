# Engineering Databases Reference

This document provides comprehensive information about engineering databases for thermodynamic, fluid, and material properties. Each database includes access methods, authentication, query examples, and practical Python code.

## Table of Contents
- [NIST REFPROP](#nist-refprop)
- [CoolProp](#coolprop)
- [NASA Earthdata](#nasa-earthdata)
- [Material Properties Databases](#material-properties-databases)
- [Pump Performance Databases](#pump-performance-databases)
- [Fluid Viscosity Databases](#fluid-viscosity-databases)
- [Cavitation and Vapor Pressure Data](#cavitation-and-vapor-pressure-data)
- [Turbulence Model Parameters](#turbulence-model-parameters)
- [Hydraulic Components Database](#hydraulic-components-database)
- [Thermodynamic Tables](#thermodynamic-tables)

---

## NIST REFPROP

### Description
NIST Reference Fluid Thermodynamic and Transport Properties Database (REFPROP) is the industry-standard database for accurate thermodynamic and transport properties of pure fluids and mixtures. Provides high-accuracy equations of state for 147+ fluids.

### Access Methods
- **Desktop Application**: REFPROP GUI (Windows/Mac/Linux)
- **Python Package**: `ctREFPROP` (wrapper for shared library)
- **MATLAB**: Built-in interface
- **Excel Add-in**: Available with purchase

### Authentication & Licensing
- **License Required**: Yes - Commercial license ($325-$425 USD)
- **Purchase**: https://www.nist.gov/srd/refprop
- **Academic Discount**: Available for educational institutions
- **Installation**: Requires license file and compiled libraries

### Installation

```bash
# Purchase and download REFPROP from NIST
# Install the main REFPROP package (includes shared libraries)

# Install Python wrapper
pip install ctREFPROP

# Set environment variable to REFPROP installation
export RPPREFIX=/opt/refprop  # Linux/Mac
# or set RPPREFIX in Windows environment variables
```

### Python Usage Examples

```python
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import os

# Initialize REFPROP
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
RP.SETPATHdll(os.environ['RPPREFIX'])

# Set units (SI, MASS BASIS)
MASS_BASE_SI = RP.GETENUMdll(0, "MASS BASE SI").iEnum

# Example 1: Water properties at 20°C and 1 bar
fluid = "Water"
T = 20 + 273.15  # K
P = 1.0  # bar = 100 kPa

# Call REFPROP
result = RP.REFPROPdll(fluid, "TP", "D;H;S;CP;CV;W;VIS;TCX;PRANDTL",
                       MASS_BASE_SI, 0, 0, T, P * 1e5, [1.0])

print(f"Water Properties at {T-273.15}°C, {P} bar:")
print(f"  Density: {result.Output[0]:.3f} kg/m³")
print(f"  Enthalpy: {result.Output[1]:.2f} J/kg")
print(f"  Entropy: {result.Output[2]:.2f} J/(kg·K)")
print(f"  Cp: {result.Output[3]:.1f} J/(kg·K)")
print(f"  Cv: {result.Output[4]:.1f} J/(kg·K)")
print(f"  Speed of sound: {result.Output[5]:.1f} m/s")
print(f"  Viscosity: {result.Output[6]*1e6:.3f} μPa·s")
print(f"  Thermal conductivity: {result.Output[7]:.4f} W/(m·K)")
print(f"  Prandtl number: {result.Output[8]:.3f}")

# Example 2: Refrigerant R134a saturation properties
fluid = "R134a"
T_sat = 25 + 273.15  # Saturation temperature (K)

# Get saturation pressure and properties
result = RP.REFPROPdll(fluid, "T", "P;DL;DV;HL;HV;SIGMA",
                       MASS_BASE_SI, 0, 0, T_sat, 0, [1.0])

P_sat = result.Output[0] / 1e5  # Convert to bar
rho_liquid = result.Output[1]
rho_vapor = result.Output[2]
h_liquid = result.Output[3]
h_vapor = result.Output[4]
surface_tension = result.Output[5]

print(f"\nR134a Saturation Properties at {T_sat-273.15}°C:")
print(f"  Saturation pressure: {P_sat:.3f} bar")
print(f"  Liquid density: {rho_liquid:.2f} kg/m³")
print(f"  Vapor density: {rho_vapor:.3f} kg/m³")
print(f"  Liquid enthalpy: {h_liquid:.1f} J/kg")
print(f"  Vapor enthalpy: {h_vapor:.1f} J/kg")
print(f"  Latent heat: {h_vapor - h_liquid:.1f} J/kg")
print(f"  Surface tension: {surface_tension*1e3:.3f} mN/m")

# Example 3: Mixture properties (60% Methane, 40% Ethane by mole)
mixture = "Methane;Ethane"
composition = [0.6, 0.4]  # Mole fractions
T = 300  # K
P = 50e5  # Pa

result = RP.REFPROPdll(mixture, "TP", "D;H;CP;VIS;TCX",
                       MASS_BASE_SI, 0, 0, T, P, composition)

print(f"\nMixture (60% CH4, 40% C2H6) at {T}K, {P/1e5} bar:")
print(f"  Density: {result.Output[0]:.3f} kg/m³")
print(f"  Enthalpy: {result.Output[1]:.2f} J/kg")
print(f"  Cp: {result.Output[2]:.1f} J/(kg·K)")
print(f"  Viscosity: {result.Output[3]*1e6:.3f} μPa·s")
print(f"  Thermal conductivity: {result.Output[4]:.4f} W/(m·K)")

# Example 4: Property along isentropic expansion
fluid = "Nitrogen"
T1 = 300  # K
P1 = 100e5  # Pa
P2 = 10e5  # Pa

# Get initial entropy
result1 = RP.REFPROPdll(fluid, "TP", "S;H;D", MASS_BASE_SI, 0, 0, T1, P1, [1.0])
s1 = result1.Output[0]
h1 = result1.Output[1]

# Get final state at constant entropy
result2 = RP.REFPROPdll(fluid, "PS", "T;H;D", MASS_BASE_SI, 0, 0, P2, s1, [1.0])
T2 = result2.Output[0]
h2 = result2.Output[1]

print(f"\nIsentropic Expansion of Nitrogen:")
print(f"  Initial: T={T1}K, P={P1/1e5} bar, h={h1:.1f} J/kg")
print(f"  Final: T={T2:.1f}K, P={P2/1e5} bar, h={h2:.1f} J/kg")
print(f"  Isentropic work: {h1-h2:.1f} J/kg")
```

### Data Format
- **Input**: Temperature (K), Pressure (Pa), or other state variables
- **Output**: SI units (configurable)
- **Fluids**: Pure substances and predefined mixtures
- **Accuracy**: Typically 0.1% for density, 1% for transport properties

### Rate Limits
- No rate limits (local library)
- Performance: ~0.1-1 ms per property call

### Best Practices
1. Cache property calls for repeated calculations
2. Use appropriate fluid files (verify spelling)
3. Check quality flag in output for two-phase regions
4. Validate results against REFPROP GUI for critical applications
5. Use IAPWS-95 formulation for water (most accurate)

---

## CoolProp

### Description
Open-source, cross-platform thermodynamic and transport property database supporting 122 pure fluids, pseudo-pure fluids (air), humid air, and incompressible fluids. Free alternative to REFPROP with good accuracy.

### Access Methods
- **Python Package**: `CoolProp`
- **C++/C**: Native library
- **MATLAB/Octave**: Wrapper available
- **Excel Add-in**: Available
- **Web API**: Online calculator available

### Authentication & Licensing
- **License**: MIT License (Free, Open Source)
- **No authentication required**
- **Commercial use**: Permitted

### Installation

```bash
# Install via pip
pip install CoolProp

# Install from conda
conda install -c conda-forge coolprop

# Install development version
pip install git+https://github.com/CoolProp/CoolProp
```

### Python Usage Examples

```python
import CoolProp.CoolProp as CP
from CoolProp.HumidAirProp import HAPropsSI
import numpy as np
import matplotlib.pyplot as plt

# Example 1: Basic property call for water
T = 20 + 273.15  # Temperature (K)
P = 101325  # Pressure (Pa)

rho = CP.PropsSI('D', 'T', T, 'P', P, 'Water')  # Density
cp = CP.PropsSI('C', 'T', T, 'P', P, 'Water')  # Specific heat
mu = CP.PropsSI('V', 'T', T, 'P', P, 'Water')  # Viscosity
k = CP.PropsSI('L', 'T', T, 'P', P, 'Water')  # Thermal conductivity
Pr = CP.PropsSI('PRANDTL', 'T', T, 'P', P, 'Water')  # Prandtl number

print(f"Water at {T-273.15}°C, {P/1000} kPa:")
print(f"  Density: {rho:.2f} kg/m³")
print(f"  Specific heat: {cp:.1f} J/(kg·K)")
print(f"  Dynamic viscosity: {mu*1000:.3f} mPa·s")
print(f"  Thermal conductivity: {k:.3f} W/(m·K)")
print(f"  Prandtl number: {Pr:.2f}")

# Example 2: Saturation properties
fluid = "R134a"
T_sat = 25 + 273.15  # K

P_sat = CP.PropsSI('P', 'T', T_sat, 'Q', 0, fluid)
h_liquid = CP.PropsSI('H', 'T', T_sat, 'Q', 0, fluid)
h_vapor = CP.PropsSI('H', 'T', T_sat, 'Q', 1, fluid)
s_liquid = CP.PropsSI('S', 'T', T_sat, 'Q', 0, fluid)
s_vapor = CP.PropsSI('S', 'T', T_sat, 'Q', 1, fluid)

print(f"\n{fluid} Saturation at {T_sat-273.15}°C:")
print(f"  Pressure: {P_sat/1000:.2f} kPa")
print(f"  h_fg: {(h_vapor-h_liquid)/1000:.2f} kJ/kg")
print(f"  s_fg: {(s_vapor-s_liquid)/1000:.3f} kJ/(kg·K)")

# Example 3: Vapor pressure calculation
T_range = np.linspace(0, 100, 50) + 273.15
P_vapor = [CP.PropsSI('P', 'T', T, 'Q', 0, 'Water') for T in T_range]

plt.figure(figsize=(10, 6))
plt.semilogy(T_range - 273.15, np.array(P_vapor) / 1000, 'b-', linewidth=2)
plt.xlabel('Temperature (°C)', fontsize=12)
plt.ylabel('Vapor Pressure (kPa)', fontsize=12)
plt.title('Water Vapor Pressure Curve', fontsize=14, fontweight='bold')
plt.grid(True, which='both', alpha=0.3)
plt.savefig('vapor_pressure_curve.png', dpi=300, bbox_inches='tight')
print("\n✓ Vapor pressure curve saved")

# Example 4: Refrigeration cycle calculation
# Evaporator at -10°C, Condenser at 40°C
T_evap = -10 + 273.15
T_cond = 40 + 273.15
fluid = "R134a"

# State 1: Saturated vapor leaving evaporator
P1 = CP.PropsSI('P', 'T', T_evap, 'Q', 1, fluid)
h1 = CP.PropsSI('H', 'T', T_evap, 'Q', 1, fluid)
s1 = CP.PropsSI('S', 'T', T_evap, 'Q', 1, fluid)

# State 2: After isentropic compression
P2 = CP.PropsSI('P', 'T', T_cond, 'Q', 0, fluid)
h2s = CP.PropsSI('H', 'P', P2, 'S', s1, fluid)  # Isentropic
h2 = h1 + (h2s - h1) / 0.85  # Actual (85% efficiency)

# State 3: Saturated liquid leaving condenser
h3 = CP.PropsSI('H', 'T', T_cond, 'Q', 0, fluid)

# State 4: After throttling (h4 = h3)
h4 = h3

# Performance metrics
Q_evap = h1 - h4  # Cooling capacity per unit mass
W_comp = h2 - h1  # Compressor work per unit mass
COP = Q_evap / W_comp

print(f"\nRefrigeration Cycle ({fluid}):")
print(f"  Evaporator temp: {T_evap-273.15}°C ({P1/1e5:.2f} bar)")
print(f"  Condenser temp: {T_cond-273.15}°C ({P2/1e5:.2f} bar)")
print(f"  Cooling capacity: {Q_evap/1000:.2f} kJ/kg")
print(f"  Compressor work: {W_comp/1000:.2f} kJ/kg")
print(f"  COP: {COP:.2f}")

# Example 5: Humid air properties
T_db = 25 + 273.15  # Dry bulb temperature (K)
P = 101325  # Pressure (Pa)
RH = 0.60  # Relative humidity (60%)

# Calculate humid air properties
h_ha = HAPropsSI('H', 'T', T_db, 'P', P, 'R', RH)  # Enthalpy
W = HAPropsSI('W', 'T', T_db, 'P', P, 'R', RH)  # Humidity ratio
T_dp = HAPropsSI('D', 'T', T_db, 'P', P, 'R', RH)  # Dew point
T_wb = HAPropsSI('B', 'T', T_db, 'P', P, 'R', RH)  # Wet bulb
V = HAPropsSI('V', 'T', T_db, 'P', P, 'R', RH)  # Specific volume

print(f"\nHumid Air at {T_db-273.15}°C, {RH*100}% RH:")
print(f"  Enthalpy: {h_ha/1000:.2f} kJ/kg_da")
print(f"  Humidity ratio: {W*1000:.2f} g/kg_da")
print(f"  Dew point: {T_dp-273.15:.1f}°C")
print(f"  Wet bulb: {T_wb-273.15:.1f}°C")
print(f"  Specific volume: {V:.3f} m³/kg_da")

# Example 6: Getting fluid information
fluids = CP.get_global_param_string("fluids_list").split(',')
print(f"\nAvailable fluids: {len(fluids)}")
print(f"First 10: {fluids[:10]}")

# Get critical point
T_crit = CP.PropsSI('Tcrit', 'Water')
P_crit = CP.PropsSI('pcrit', 'Water')
print(f"\nWater critical point:")
print(f"  T_crit: {T_crit-273.15:.2f}°C")
print(f"  P_crit: {P_crit/1e6:.2f} MPa")
```

### Available Fluids
**Common fluids**: Water, Air, Nitrogen, Oxygen, CO2, Ammonia, Propane, Methane
**Refrigerants**: R134a, R410A, R32, R404A, R407C, R22, R717 (NH3)
**Hydrocarbons**: Methane, Ethane, Propane, n-Butane, Isobutane
**Cryogens**: Helium, Hydrogen, Neon, Argon
**Incompressibles**: Water-glycol mixtures, thermal oils

### Data Format
- **Input/Output**: SI units by default
- **Temperature**: K
- **Pressure**: Pa
- **Specific properties**: Per kg
- **Quality**: 0 (liquid) to 1 (vapor)

### Rate Limits
- No rate limits (local library)
- Performance: 0.01-0.1 ms per call

### Accuracy
- Pure fluids: 0.1-1% for most properties
- Humid air: 0.5% for psychrometric properties
- Transport properties: 1-5%
- Less accurate than REFPROP but sufficient for most engineering

### Best Practices
1. Use `PropsSI()` for general fluid properties
2. Use `HAPropsSI()` specifically for humid air
3. Cache calls when computing property tables
4. Check for two-phase region with quality (Q)
5. Use try-except for property calls near critical point
6. Validate against experimental data for critical applications

### Common Input/Output Keys
```
# Inputs
T: Temperature (K)
P: Pressure (Pa)
D: Density (kg/m³)
H: Enthalpy (J/kg)
S: Entropy (J/kg/K)
Q: Quality (0-1)

# Outputs
C: Specific heat (J/kg/K)
L: Thermal conductivity (W/m/K)
V: Dynamic viscosity (Pa·s)
A: Speed of sound (m/s)
PRANDTL: Prandtl number
CONDUCTIVITY: Thermal conductivity
```

---

## NASA Earthdata

### Description
NASA's Earth Observing System Data and Information System (EOSDIS) provides atmospheric, oceanographic, and aerospace data. Useful for atmospheric properties, aerospace design, and environmental analysis.

### Access Methods
- **Web Portal**: https://earthdata.nasa.gov/
- **API**: REST API with JSON responses
- **Python Package**: `earthdata` (unofficial wrappers available)
- **Direct Download**: FTP/HTTPS data access
- **OPeNDAP**: Remote data access protocol

### Authentication & Licensing
- **License**: Free (US taxpayer-funded)
- **Registration Required**: Yes - Create Earthdata Login
- **API Token**: Required for programmatic access
- **Rate Limits**: Reasonable use policy

### Setup & Registration

```bash
# 1. Register at https://urs.earthdata.nasa.gov/users/new
# 2. Create .netrc file for authentication

cat > ~/.netrc << EOF
machine urs.earthdata.nasa.gov
login YOUR_USERNAME
password YOUR_PASSWORD
EOF

chmod 0600 ~/.netrc

# 3. Install Python packages
pip install requests
pip install h5py  # For HDF5 data files
pip install netCDF4  # For NetCDF data files
pip install pydap  # For OPeNDAP access
```

### Python Usage Examples

```python
import requests
from datetime import datetime
import json
import numpy as np
import matplotlib.pyplot as plt

# Example 1: Search for atmospheric data products
def search_earthdata(short_name, start_date, end_date, bbox=None):
    """
    Search NASA Earthdata CMR (Common Metadata Repository)

    short_name: Dataset short name (e.g., 'AIRS3STD', 'MERRA2')
    start_date/end_date: 'YYYY-MM-DD' format
    bbox: [west, south, east, north] in degrees
    """
    url = "https://cmr.earthdata.nasa.gov/search/granules.json"

    params = {
        'short_name': short_name,
        'temporal': f"{start_date}T00:00:00Z,{end_date}T23:59:59Z",
        'page_size': 100
    }

    if bbox:
        params['bounding_box'] = ','.join(map(str, bbox))

    response = requests.get(url, params=params)
    data = response.json()

    return data['feed']['entry']

# Search for AIRS atmospheric data (temperature, humidity profiles)
results = search_earthdata(
    short_name='AIRS3STD',
    start_date='2024-01-01',
    end_date='2024-01-02',
    bbox=[-180, -90, 180, 90]
)

print(f"Found {len(results)} granules")
print(f"Example granule: {results[0]['title']}")

# Example 2: Download atmospheric profile data
def download_file(url, filename, auth):
    """Download file from NASA Earthdata"""
    with requests.get(url, auth=auth, stream=True) as r:
        r.raise_for_status()
        with open(filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return filename

# Earthdata credentials (use .netrc or provide directly)
from requests.auth import HTTPBasicAuth
auth = HTTPBasicAuth('YOUR_USERNAME', 'YOUR_PASSWORD')

# Example 3: Access MERRA-2 reanalysis data (atmospheric properties)
# MERRA-2 provides global atmospheric conditions at various altitudes

def get_merra2_atmospheric_profile(lat, lon, date):
    """
    Get atmospheric profile from MERRA-2
    (Simplified example - actual implementation requires OPeNDAP)
    """
    # This is a conceptual example
    # Real implementation uses OPeNDAP or downloads HDF files

    # Standard atmosphere approximation for demonstration
    altitudes = np.linspace(0, 50000, 50)  # meters

    # US Standard Atmosphere 1976
    T0 = 288.15  # K at sea level
    P0 = 101325  # Pa at sea level
    L = 0.0065  # K/m temperature lapse rate
    g = 9.80665  # m/s²
    M = 0.0289644  # kg/mol
    R = 8.31447  # J/(mol·K)

    temperatures = np.zeros_like(altitudes)
    pressures = np.zeros_like(altitudes)
    densities = np.zeros_like(altitudes)

    for i, h in enumerate(altitudes):
        if h <= 11000:  # Troposphere
            T = T0 - L * h
            P = P0 * (T / T0) ** (g * M / (R * L))
        elif h <= 25000:  # Lower stratosphere
            T = 216.65  # K (isothermal)
            P = 22632 * np.exp(-g * M * (h - 11000) / (R * T))
        else:
            T = 216.65 + 0.0028 * (h - 25000)
            P = 2488.5 * (T / 216.65) ** (-g * M / (R * 0.0028))

        rho = P * M / (R * T)

        temperatures[i] = T
        pressures[i] = P
        densities[i] = rho

    return altitudes, temperatures, pressures, densities

# Get atmospheric profile
alt, temp, pres, dens = get_merra2_atmospheric_profile(40.0, -105.0, '2024-01-01')

print("\nAtmospheric Profile (Standard Atmosphere):")
print(f"Altitude (km)  Temp (K)  Pressure (kPa)  Density (kg/m³)")
for i in [0, 10, 20, 30, 40]:
    print(f"{alt[i]/1000:8.1f}     {temp[i]:7.2f}    {pres[i]/1000:8.2f}      {dens[i]:8.4f}")

# Plot atmospheric profile
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

ax1.plot(temp - 273.15, alt / 1000, 'r-', linewidth=2)
ax1.set_xlabel('Temperature (°C)', fontsize=12)
ax1.set_ylabel('Altitude (km)', fontsize=12)
ax1.set_title('Temperature Profile', fontweight='bold')
ax1.grid(True, alpha=0.3)

ax2.semilogx(pres / 1000, alt / 1000, 'b-', linewidth=2)
ax2.set_xlabel('Pressure (kPa)', fontsize=12)
ax2.set_ylabel('Altitude (km)', fontsize=12)
ax2.set_title('Pressure Profile', fontweight='bold')
ax2.grid(True, alpha=0.3)

ax3.semilogx(dens, alt / 1000, 'g-', linewidth=2)
ax3.set_xlabel('Density (kg/m³)', fontsize=12)
ax3.set_ylabel('Altitude (km)', fontsize=12)
ax3.set_title('Density Profile', fontweight='bold')
ax3.grid(True, alpha=0.3)

# Speed of sound
a = np.sqrt(1.4 * 287 * temp)  # m/s
ax4.plot(a, alt / 1000, 'm-', linewidth=2)
ax4.set_xlabel('Speed of Sound (m/s)', fontsize=12)
ax4.set_ylabel('Altitude (km)', fontsize=12)
ax4.set_title('Speed of Sound Profile', fontweight='bold')
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('atmospheric_profile.png', dpi=300, bbox_inches='tight')
print("\n✓ Atmospheric profile saved")

# Example 4: Calculate dynamic pressure for aerospace applications
def calculate_flight_conditions(altitude, velocity):
    """Calculate flight conditions at given altitude and velocity"""
    # Get atmospheric properties at altitude
    idx = np.argmin(np.abs(alt - altitude))
    T = temp[idx]
    P = pres[idx]
    rho = dens[idx]

    # Calculate derived quantities
    mu = 1.458e-6 * T**1.5 / (T + 110.4)  # Sutherland's formula
    a = np.sqrt(1.4 * 287 * T)  # Speed of sound
    M = velocity / a  # Mach number
    q = 0.5 * rho * velocity**2  # Dynamic pressure
    Re_per_m = rho * velocity / mu  # Reynolds number per meter

    return {
        'altitude': altitude,
        'velocity': velocity,
        'temperature': T - 273.15,
        'pressure': P / 1000,
        'density': rho,
        'viscosity': mu,
        'mach': M,
        'dynamic_pressure': q / 1000,
        'reynolds_per_m': Re_per_m
    }

# Calculate for typical cruise conditions
flight = calculate_flight_conditions(altitude=11000, velocity=250)  # 11km, 250 m/s

print(f"\nFlight Conditions at {flight['altitude']/1000}km, {flight['velocity']} m/s:")
print(f"  Temperature: {flight['temperature']:.1f}°C")
print(f"  Pressure: {flight['pressure']:.2f} kPa")
print(f"  Density: {flight['density']:.4f} kg/m³")
print(f"  Mach number: {flight['mach']:.3f}")
print(f"  Dynamic pressure: {flight['dynamic_pressure']:.2f} kPa")
print(f"  Reynolds number/m: {flight['reynolds_per_m']:.2e}")
```

### Available Datasets
- **AIRS**: Atmospheric Infrared Sounder (temperature/humidity profiles)
- **MERRA-2**: Modern-Era Retrospective analysis (global atmospheric reanalysis)
- **MODIS**: Moderate Resolution Imaging Spectroradiometer
- **GOES**: Geostationary weather satellites
- **COSMIC**: GPS radio occultation atmospheric profiles

### Data Format
- **HDF5**: Hierarchical Data Format
- **NetCDF**: Network Common Data Form
- **GeoTIFF**: Georeferenced raster data
- **JSON**: Metadata responses

### Rate Limits
- No hard limits but reasonable use policy
- Bulk downloads: Contact EOSDIS for large requests
- Recommended: < 1000 requests/hour

### Best Practices
1. Use OPeNDAP for subsetting data before download
2. Cache downloaded data locally
3. Use bulk download for large time series
4. Check data quality flags
5. Cite datasets in publications

---

## Material Properties Databases

### Description
Comprehensive databases for structural and engineering materials including metals, alloys, composites, ceramics, and polymers. Includes mechanical, thermal, and electrical properties.

### Access Methods
- **NIST Materials Data**: https://www.nist.gov/mml/materials-data
- **MatWeb**: http://www.matweb.com/ (Free with registration)
- **ASM Handbooks**: Subscription required
- **Python Package**: `matminer` (materials data mining)
- **CES EduPack**: Educational license available

### Installation

```bash
# Install matminer for materials informatics
pip install matminer

# Install materials project API client
pip install mp-api

# Optional: ASM handbooks digital access (subscription required)
```

### Python Usage Examples

```python
from matminer.datasets import load_dataset
import pandas as pd
import numpy as np

# Example 1: Load common materials database
# Note: First run will download dataset
df_steels = load_dataset("steel_strength")
print(f"Steel database: {len(df_steels)} entries")
print(df_steels.head())

# Example 2: Manual material property database
class MaterialsDatabase:
    """Simple materials property database"""

    def __init__(self):
        self.materials = {
            # Common structural steels
            'A36': {
                'name': 'ASTM A36 Structural Steel',
                'type': 'Carbon Steel',
                'density': 7850,  # kg/m³
                'yield_strength': 250e6,  # Pa
                'ultimate_strength': 400e6,  # Pa
                'elastic_modulus': 200e9,  # Pa
                'poisson_ratio': 0.26,
                'thermal_conductivity': 51.9,  # W/(m·K)
                'specific_heat': 486,  # J/(kg·K)
                'thermal_expansion': 12e-6,  # 1/K
                'elongation': 0.20  # 20%
            },
            '304_SS': {
                'name': '304 Stainless Steel',
                'type': 'Austenitic Stainless',
                'density': 8000,
                'yield_strength': 215e6,
                'ultimate_strength': 505e6,
                'elastic_modulus': 193e9,
                'poisson_ratio': 0.29,
                'thermal_conductivity': 16.2,
                'specific_heat': 500,
                'thermal_expansion': 17.3e-6,
                'elongation': 0.40
            },
            '6061-T6': {
                'name': '6061-T6 Aluminum Alloy',
                'type': 'Aluminum Alloy',
                'density': 2700,
                'yield_strength': 276e6,
                'ultimate_strength': 310e6,
                'elastic_modulus': 68.9e9,
                'poisson_ratio': 0.33,
                'thermal_conductivity': 167,
                'specific_heat': 896,
                'thermal_expansion': 23.6e-6,
                'elongation': 0.12
            },
            'Ti-6Al-4V': {
                'name': 'Titanium Ti-6Al-4V (Grade 5)',
                'type': 'Titanium Alloy',
                'density': 4430,
                'yield_strength': 880e6,
                'ultimate_strength': 950e6,
                'elastic_modulus': 113.8e9,
                'poisson_ratio': 0.342,
                'thermal_conductivity': 6.7,
                'specific_heat': 526,
                'thermal_expansion': 8.6e-6,
                'elongation': 0.14
            },
            'Inconel_718': {
                'name': 'Inconel 718 (High-Temp Nickel Alloy)',
                'type': 'Nickel Superalloy',
                'density': 8190,
                'yield_strength': 1035e6,
                'ultimate_strength': 1275e6,
                'elastic_modulus': 200e9,
                'poisson_ratio': 0.29,
                'thermal_conductivity': 11.4,
                'specific_heat': 435,
                'thermal_expansion': 13e-6,
                'elongation': 0.12
            },
            'Copper_C11000': {
                'name': 'Copper C11000 (ETP)',
                'type': 'Pure Copper',
                'density': 8940,
                'yield_strength': 69e6,
                'ultimate_strength': 220e6,
                'elastic_modulus': 117e9,
                'poisson_ratio': 0.34,
                'thermal_conductivity': 391,
                'specific_heat': 385,
                'thermal_expansion': 17e-6,
                'elongation': 0.45
            },
            'Brass_C36000': {
                'name': 'Brass C36000 (Free-Cutting)',
                'type': 'Copper Alloy',
                'density': 8500,
                'yield_strength': 124e6,
                'ultimate_strength': 338e6,
                'elastic_modulus': 97e9,
                'poisson_ratio': 0.33,
                'thermal_conductivity': 115,
                'specific_heat': 380,
                'thermal_expansion': 20.5e-6,
                'elongation': 0.53
            }
        }

    def get_property(self, material, property_name):
        """Get specific property of a material"""
        if material not in self.materials:
            raise ValueError(f"Material '{material}' not found in database")

        mat = self.materials[material]
        if property_name not in mat:
            raise ValueError(f"Property '{property_name}' not available for {material}")

        return mat[property_name]

    def get_all_properties(self, material):
        """Get all properties of a material"""
        if material not in self.materials:
            raise ValueError(f"Material '{material}' not found")
        return self.materials[material]

    def compare_materials(self, materials, property_name):
        """Compare specific property across multiple materials"""
        comparison = {}
        for mat in materials:
            comparison[mat] = self.get_property(mat, property_name)
        return comparison

    def strength_to_weight_ratio(self, material):
        """Calculate strength-to-weight ratio"""
        rho = self.get_property(material, 'density')
        sigma_y = self.get_property(material, 'yield_strength')
        return sigma_y / rho  # Specific strength (Pa·m³/kg = m²/s²)

# Initialize database
mat_db = MaterialsDatabase()

# Example 3: Get material properties
material = 'Ti-6Al-4V'
props = mat_db.get_all_properties(material)

print(f"\n{props['name']}:")
print(f"  Density: {props['density']} kg/m³")
print(f"  Yield strength: {props['yield_strength']/1e6:.0f} MPa")
print(f"  Ultimate strength: {props['ultimate_strength']/1e6:.0f} MPa")
print(f"  Elastic modulus: {props['elastic_modulus']/1e9:.1f} GPa")
print(f"  Thermal conductivity: {props['thermal_conductivity']:.1f} W/(m·K)")
print(f"  Specific heat: {props['specific_heat']} J/(kg·K)")

# Example 4: Material selection based on strength-to-weight ratio
materials = ['A36', '304_SS', '6061-T6', 'Ti-6Al-4V', 'Inconel_718']
specific_strengths = {mat: mat_db.strength_to_weight_ratio(mat)
                     for mat in materials}

print("\nSpecific Strength Comparison (Yield/Density):")
for mat, spec_str in sorted(specific_strengths.items(),
                            key=lambda x: x[1], reverse=True):
    print(f"  {mat}: {spec_str/1000:.1f} kJ/kg = {np.sqrt(spec_str):.0f} m/s")

# Example 5: Thermal stress calculation
def thermal_stress(material, delta_T, db):
    """Calculate thermal stress if constrained"""
    E = db.get_property(material, 'elastic_modulus')
    alpha = db.get_property(material, 'thermal_expansion')
    sigma_thermal = E * alpha * delta_T
    return sigma_thermal

delta_T = 100  # K temperature change
print(f"\nThermal stress for ΔT = {delta_T}K (if fully constrained):")
for mat in ['A36', '304_SS', '6061-T6', 'Ti-6Al-4V']:
    sigma = thermal_stress(mat, delta_T, mat_db)
    sigma_y = mat_db.get_property(mat, 'yield_strength')
    safety_factor = sigma_y / sigma
    print(f"  {mat}: {sigma/1e6:.1f} MPa (SF = {safety_factor:.2f})")

# Example 6: Beam deflection calculation
def beam_deflection(material, length, width, height, load, db):
    """Calculate beam deflection under center load (simply supported)"""
    E = db.get_property(material, 'elastic_modulus')
    I = width * height**3 / 12  # Moment of inertia
    delta = load * length**3 / (48 * E * I)
    return delta

# Compare beam deflection for different materials
L = 1.0  # m
W = 0.05  # m
H = 0.1  # m
F = 1000  # N

print(f"\nBeam deflection (L={L}m, {W*1000}x{H*1000}mm, F={F}N):")
for mat in ['A36', '6061-T6', 'Ti-6Al-4V']:
    delta = beam_deflection(mat, L, W, H, F, mat_db)
    print(f"  {mat}: {delta*1000:.3f} mm")
```

### Data Sources & APIs

```python
# Materials Project API (requires API key from materialsproject.org)
from mp_api.client import MPRester

# Get free API key from https://materialsproject.org/api
API_KEY = "YOUR_API_KEY_HERE"

# Example: Query Materials Project
with MPRester(API_KEY) as mpr:
    # Search for silicon
    docs = mpr.materials.summary.search(
        formula="Si",
        fields=["material_id", "formula_pretty", "band_gap", "density"]
    )

    for doc in docs[:5]:
        print(f"{doc.material_id}: {doc.formula_pretty}")
        print(f"  Band gap: {doc.band_gap} eV")
        print(f"  Density: {doc.density:.2f} g/cm³")
```

### Best Practices
1. Verify material grade and heat treatment condition
2. Check temperature dependence of properties
3. Use safety factors for design calculations
4. Cross-reference multiple sources for critical applications
5. Consider manufacturing tolerances and variability

---

## Pump Performance Databases

### Description
Manufacturer pump performance data and empirical correlations for centrifugal, axial, and positive displacement pumps. Includes performance curves, NPSH requirements, and efficiency data.

### Access Methods
- **Manufacturer Catalogs**: PDF/Web (Grundfos, KSB, Flowserve, etc.)
- **Pump Selection Software**: Vendor-specific (free registration)
- **Empirical Correlations**: Published in handbooks
- **Custom Database**: Build from manufacturer data

### Python Implementation

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

class PumpDatabase:
    """Pump performance database with empirical correlations"""

    def __init__(self):
        self.pumps = {}
        self._load_standard_pumps()

    def _load_standard_pumps(self):
        """Load typical pump performance data"""

        # Example: End-suction centrifugal pump (based on Grundfos NB series)
        self.pumps['NB_50-200'] = {
            'type': 'end_suction_centrifugal',
            'manufacturer': 'Generic',
            'impeller_diameter': 0.200,  # m
            'speed': 1450,  # rpm
            'Q_BEP': 80 / 3600,  # m³/s (80 m³/h)
            'H_BEP': 32,  # m
            'eta_BEP': 0.72,
            'P_rated': 11000,  # W
            'NPSH_req_BEP': 3.5,  # m
            # Performance curve coefficients: H = a + b*Q + c*Q²
            'head_coeffs': [40.0, -50.0, -2500.0],
            # Efficiency curve parameters
            'eta_coeffs': [0.72, 80/3600],  # max_eta, Q_at_max_eta
            # NPSH curve: NPSH = a + b*Q + c*Q²
            'NPSH_coeffs': [2.5, 0, 5000]
        }

        # Multistage centrifugal pump
        self.pumps['CR_64-3'] = {
            'type': 'multistage_centrifugal',
            'manufacturer': 'Generic',
            'stages': 3,
            'impeller_diameter': 0.170,
            'speed': 2900,  # rpm
            'Q_BEP': 60 / 3600,
            'H_BEP': 60,
            'eta_BEP': 0.75,
            'P_rated': 15000,
            'NPSH_req_BEP': 4.0,
            'head_coeffs': [72.0, -80.0, -4000.0],
            'eta_coeffs': [0.75, 60/3600],
            'NPSH_coeffs': [3.0, 0, 6000]
        }

        # Submersible pump
        self.pumps['SP_95-4'] = {
            'type': 'submersible',
            'manufacturer': 'Generic',
            'stages': 4,
            'impeller_diameter': 0.180,
            'speed': 2900,
            'Q_BEP': 95 / 3600,
            'H_BEP': 85,
            'eta_BEP': 0.78,
            'P_rated': 28000,
            'NPSH_req_BEP': 5.5,
            'head_coeffs': [100.0, -100.0, -3500.0],
            'eta_coeffs': [0.78, 95/3600],
            'NPSH_coeffs': [4.5, 0, 5500]
        }

    def head_curve(self, pump_id, Q):
        """Calculate head at given flow rate"""
        pump = self.pumps[pump_id]
        a, b, c = pump['head_coeffs']
        H = a + b * Q + c * Q**2
        return max(0, H)  # Non-negative head

    def efficiency_curve(self, pump_id, Q):
        """Calculate efficiency at given flow rate"""
        pump = self.pumps[pump_id]
        eta_max, Q_BEP = pump['eta_coeffs']
        # Gaussian-like efficiency curve
        eta = eta_max * np.exp(-3 * ((Q - Q_BEP) / Q_BEP)**2)
        return np.clip(eta, 0, 1)

    def power_curve(self, pump_id, Q, rho=998, g=9.81):
        """Calculate shaft power at given flow rate"""
        H = self.head_curve(pump_id, Q)
        eta = self.efficiency_curve(pump_id, Q)

        if eta < 0.01:  # Avoid division by very small efficiency
            return np.nan

        P = rho * g * Q * H / eta
        return P

    def NPSH_required(self, pump_id, Q):
        """Calculate NPSH required at given flow rate"""
        pump = self.pumps[pump_id]
        a, b, c = pump['NPSH_coeffs']
        NPSH = a + b * Q + c * Q**2
        return max(0, NPSH)

    def specific_speed(self, pump_id):
        """Calculate specific speed (Ns)"""
        pump = self.pumps[pump_id]
        n = pump['speed']
        Q = pump['Q_BEP']
        H = pump['H_BEP']

        # For multistage, use head per stage
        if 'stages' in pump:
            H = H / pump['stages']

        Ns = n * Q**0.5 / H**0.75
        return Ns

    def plot_performance(self, pump_id, save_path=None):
        """Plot pump performance curves"""
        pump = self.pumps[pump_id]
        Q_max = 1.5 * pump['Q_BEP']
        Q_range = np.linspace(0, Q_max, 100)

        # Calculate curves
        H = np.array([self.head_curve(pump_id, q) for q in Q_range])
        eta = np.array([self.efficiency_curve(pump_id, q) for q in Q_range])
        P = np.array([self.power_curve(pump_id, q) for q in Q_range])
        NPSH = np.array([self.NPSH_required(pump_id, q) for q in Q_range])

        # Create plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

        # Head curve
        ax1.plot(Q_range * 3600, H, 'b-', linewidth=2)
        ax1.plot(pump['Q_BEP'] * 3600, pump['H_BEP'], 'ro', markersize=10, label='BEP')
        ax1.set_xlabel('Flow Rate (m³/h)', fontsize=11)
        ax1.set_ylabel('Head (m)', fontsize=11)
        ax1.set_title(f"{pump_id} - Head Curve", fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.legend()

        # Efficiency curve
        ax2.plot(Q_range * 3600, eta * 100, 'g-', linewidth=2)
        ax2.plot(pump['Q_BEP'] * 3600, pump['eta_BEP'] * 100, 'ro', markersize=10, label='BEP')
        ax2.set_xlabel('Flow Rate (m³/h)', fontsize=11)
        ax2.set_ylabel('Efficiency (%)', fontsize=11)
        ax2.set_title('Efficiency Curve', fontweight='bold')
        ax2.set_ylim(0, 100)
        ax2.grid(True, alpha=0.3)
        ax2.legend()

        # Power curve
        ax3.plot(Q_range * 3600, P / 1000, 'r-', linewidth=2)
        ax3.plot(pump['Q_BEP'] * 3600, pump['P_rated'] / 1000, 'ro', markersize=10, label='Rated')
        ax3.set_xlabel('Flow Rate (m³/h)', fontsize=11)
        ax3.set_ylabel('Shaft Power (kW)', fontsize=11)
        ax3.set_title('Power Curve', fontweight='bold')
        ax3.grid(True, alpha=0.3)
        ax3.legend()

        # NPSH required
        ax4.plot(Q_range * 3600, NPSH, 'm-', linewidth=2)
        ax4.plot(pump['Q_BEP'] * 3600, pump['NPSH_req_BEP'], 'ro', markersize=10, label='BEP')
        ax4.set_xlabel('Flow Rate (m³/h)', fontsize=11)
        ax4.set_ylabel('NPSH Required (m)', fontsize=11)
        ax4.set_title('NPSH Curve', fontweight='bold')
        ax4.grid(True, alpha=0.3)
        ax4.legend()

        plt.suptitle(f"Pump Performance: {pump_id} ({pump['type']})",
                     fontsize=14, fontweight='bold', y=1.00)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')

        return fig

# Example usage
pump_db = PumpDatabase()

# List available pumps
print("Available pumps in database:")
for pump_id, pump in pump_db.pumps.items():
    print(f"\n{pump_id} ({pump['type']}):")
    print(f"  BEP: {pump['Q_BEP']*3600:.1f} m³/h @ {pump['H_BEP']:.1f}m")
    print(f"  Efficiency: {pump['eta_BEP']*100:.1f}%")
    print(f"  Speed: {pump['speed']} rpm")
    print(f"  Specific speed: {pump_db.specific_speed(pump_id):.1f}")

# Example: Get performance at operating point
pump_id = 'NB_50-200'
Q_operating = 75 / 3600  # m³/s

H = pump_db.head_curve(pump_id, Q_operating)
eta = pump_db.efficiency_curve(pump_id, Q_operating)
P = pump_db.power_curve(pump_id, Q_operating)
NPSH_req = pump_db.NPSH_required(pump_id, Q_operating)

print(f"\n{pump_id} at {Q_operating*3600:.1f} m³/h:")
print(f"  Head: {H:.2f} m")
print(f"  Efficiency: {eta*100:.1f}%")
print(f"  Power: {P/1000:.2f} kW")
print(f"  NPSH required: {NPSH_req:.2f} m")

# Plot performance curves
pump_db.plot_performance('NB_50-200', 'pump_performance.png')
print("\n✓ Pump performance curves saved")

# Example: Pump selection based on duty point
def select_pump(Q_required, H_required, db):
    """Select pump based on required duty point"""
    candidates = []

    for pump_id in db.pumps.keys():
        # Check if pump can deliver required head at required flow
        H_available = db.head_curve(pump_id, Q_required)

        if H_available >= H_required:
            eta = db.efficiency_curve(pump_id, Q_required)
            P = db.power_curve(pump_id, Q_required)

            # Calculate how close to BEP
            Q_BEP = db.pumps[pump_id]['Q_BEP']
            BEP_ratio = Q_required / Q_BEP

            candidates.append({
                'pump_id': pump_id,
                'head': H_available,
                'efficiency': eta,
                'power': P,
                'BEP_ratio': BEP_ratio,
                'score': eta * (1 - abs(1 - BEP_ratio))  # Prefer operation near BEP
            })

    # Sort by score (efficiency and proximity to BEP)
    candidates.sort(key=lambda x: x['score'], reverse=True)

    return candidates

# Select pump for duty point
Q_req = 70 / 3600  # m³/s
H_req = 30  # m

candidates = select_pump(Q_req, H_req, pump_db)

print(f"\nPump selection for {Q_req*3600:.0f} m³/h @ {H_req}m:")
for i, cand in enumerate(candidates, 1):
    print(f"{i}. {cand['pump_id']}:")
    print(f"   Head: {cand['head']:.1f}m, Efficiency: {cand['efficiency']*100:.1f}%")
    print(f"   Power: {cand['power']/1000:.1f} kW, BEP ratio: {cand['BEP_ratio']:.2f}")
```

### Data Format
- Flow rate: m³/s or m³/h
- Head: m (meters of fluid)
- Efficiency: 0-1 or percentage
- Power: W or kW
- Speed: rpm
- NPSH: m

### Best Practices
1. Always operate pumps within 70-120% of BEP flow for good efficiency
2. Verify NPSH available > NPSH required + 0.5m safety margin
3. Consider system curve when selecting pumps
4. Account for wear and fouling in long-term operation
5. Check manufacturer data for specific impeller trims

---

## Fluid Viscosity Databases

### Description
Comprehensive viscosity data for liquids and gases as functions of temperature and pressure. Critical for Reynolds number calculations, pressure drop analysis, and heat transfer.

### Access Methods
- **NIST Chemistry WebBook**: https://webbook.nist.gov/chemistry/fluid/
- **CoolProp**: Includes viscosity for many fluids
- **DIPPR Database**: Subscription required
- **Empirical Correlations**: Sutherland, Vogel, Andrade equations

### Python Implementation

```python
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

class ViscosityDatabase:
    """Viscosity database with empirical correlations"""

    @staticmethod
    def sutherland_gas(T, T0=273.15, mu0=1.716e-5, S=110.4):
        """
        Sutherland's formula for gas viscosity (air by default)

        T: Temperature (K)
        T0: Reference temperature (K)
        mu0: Reference viscosity (Pa·s)
        S: Sutherland constant (K)
        """
        mu = mu0 * (T / T0)**1.5 * (T0 + S) / (T + S)
        return mu

    @staticmethod
    def vogel_liquid(T, A, B, C):
        """
        Vogel equation for liquid viscosity

        log10(μ) = A + B / (T - C)

        Common parameters:
        Water: A=-3.7188, B=578.919, C=-137.546 (T in K, μ in Pa·s)
        """
        log_mu = A + B / (T - C)
        mu = 10**log_mu
        return mu

    @staticmethod
    def andrade_liquid(T, A, B):
        """
        Andrade equation for liquid viscosity

        μ = A * exp(B/T)
        """
        mu = A * np.exp(B / T)
        return mu

    def water_viscosity(self, T):
        """Water viscosity using empirical correlation"""
        # Valid for 0-100°C at atmospheric pressure
        # Result in Pa·s
        if T < 273.15 or T > 373.15:
            print(f"Warning: Temperature {T}K outside valid range (273-373K)")

        # Vogel equation fitted to IAPWS data
        A = -3.7188
        B = 578.919
        C = -137.546
        return self.vogel_liquid(T, A, B, C)

    def air_viscosity(self, T):
        """Air viscosity using Sutherland's formula"""
        return self.sutherland_gas(T)

    def oil_viscosity(self, T, grade='SAE_30'):
        """
        Engine oil viscosity approximation

        Grades: SAE_10W, SAE_30, SAE_50
        """
        # Simplified Vogel parameters for different grades
        params = {
            'SAE_10W': {'A': -2.5, 'B': 1200, 'C': -120},
            'SAE_30': {'A': -2.0, 'B': 1500, 'C': -100},
            'SAE_50': {'A': -1.5, 'B': 1800, 'C': -90}
        }

        if grade not in params:
            raise ValueError(f"Oil grade {grade} not supported")

        p = params[grade]
        return self.vogel_liquid(T, p['A'], p['B'], p['C'])

# Example usage
visc_db = ViscosityDatabase()

# Example 1: Water viscosity vs temperature
T_range = np.linspace(0, 100, 50) + 273.15
mu_water = [visc_db.water_viscosity(T) for T in T_range]
mu_water_cp = [CP.PropsSI('V', 'T', T, 'P', 101325, 'Water') for T in T_range]

print("Water Viscosity:")
print("T(°C)  μ_correlation(mPa·s)  μ_CoolProp(mPa·s)")
for i in [0, 10, 20, 30, 40]:
    idx = int(i / 2)
    print(f"{T_range[idx]-273.15:5.0f}      {mu_water[idx]*1000:8.4f}           "
          f"{mu_water_cp[idx]*1000:8.4f}")

# Example 2: Air viscosity vs temperature
mu_air = [visc_db.air_viscosity(T) for T in T_range]

print("\nAir Viscosity (Sutherland):")
for T, mu in zip([273.15, 293.15, 373.15, 473.15],
                 [visc_db.air_viscosity(T) for T in [273.15, 293.15, 373.15, 473.15]]):
    print(f"  {T-273.15:5.0f}°C: {mu*1e6:6.3f} μPa·s")

# Example 3: Engine oil viscosity
T_oil = np.linspace(-20, 100, 50) + 273.15
oils = ['SAE_10W', 'SAE_30', 'SAE_50']

plt.figure(figsize=(12, 6))

for oil in oils:
    mu_oil = [visc_db.oil_viscosity(T, oil) for T in T_oil]
    plt.semilogy(T_oil - 273.15, np.array(mu_oil) * 1000, linewidth=2, label=oil)

plt.xlabel('Temperature (°C)', fontsize=12)
plt.ylabel('Dynamic Viscosity (mPa·s)', fontsize=12)
plt.title('Engine Oil Viscosity vs Temperature', fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(True, which='both', alpha=0.3)
plt.savefig('oil_viscosity.png', dpi=300, bbox_inches='tight')
print("\n✓ Oil viscosity plot saved")

# Example 4: Kinematic viscosity calculation
def kinematic_viscosity(T, fluid='Water'):
    """Calculate kinematic viscosity"""
    mu = CP.PropsSI('V', 'T', T, 'P', 101325, fluid)
    rho = CP.PropsSI('D', 'T', T, 'P', 101325, fluid)
    nu = mu / rho
    return nu

print("\nKinematic Viscosity (ν = μ/ρ):")
for T in [273.15, 293.15, 313.15, 333.15]:
    nu = kinematic_viscosity(T, 'Water')
    print(f"  Water at {T-273.15:.0f}°C: {nu*1e6:.3f} mm²/s (cSt)")

# Example 5: Reynolds number calculation
def reynolds_number(V, D, T, fluid='Water'):
    """Calculate Reynolds number"""
    nu = kinematic_viscosity(T, fluid)
    Re = V * D / nu
    return Re

V = 2.0  # m/s
D = 0.1  # m
T = 293.15  # K

Re = reynolds_number(V, D, T, 'Water')
print(f"\nReynolds Number:")
print(f"  Velocity: {V} m/s, Diameter: {D} m, Temp: {T-273.15}°C")
print(f"  Re = {Re:.0f} ({'Turbulent' if Re > 4000 else 'Transitional' if Re > 2300 else 'Laminar'})")

# Example 6: Comparison with CoolProp for various fluids
fluids_to_test = ['Water', 'Ethanol', 'Methanol', 'Ammonia', 'R134a']
T_test = 293.15  # 20°C

print(f"\nViscosity comparison at {T_test-273.15}°C:")
print("Fluid          μ (mPa·s)    ν (mm²/s)")
for fluid in fluids_to_test:
    try:
        mu = CP.PropsSI('V', 'T', T_test, 'P', 101325, fluid)
        nu = kinematic_viscosity(T_test, fluid)
        print(f"{fluid:12}   {mu*1000:8.4f}     {nu*1e6:8.4f}")
    except:
        print(f"{fluid:12}   {'N/A':>8}     {'N/A':>8}")
```

### Common Viscosity Ranges (20°C, 1 atm)

| Fluid | Dynamic Viscosity (mPa·s) | Kinematic Viscosity (mm²/s) |
|-------|---------------------------|------------------------------|
| Water | 1.002 | 1.004 |
| Air | 0.0181 | 15.1 |
| Ethanol | 1.095 | 1.39 |
| Glycerin | 1412 | 1121 |
| SAE 30 Oil | ~200 | ~220 |
| Gasoline | 0.4-0.6 | 0.5-0.7 |

### Best Practices
1. Always specify temperature when reporting viscosity
2. Use dynamic viscosity (μ) for stress calculations
3. Use kinematic viscosity (ν) for Reynolds number
4. Check if viscosity model is valid for pressure range
5. Consider non-Newtonian behavior for polymers, slurries

---

## Cavitation and Vapor Pressure Data

### Description
Vapor pressure data for cavitation analysis in pumps, valves, and hydraulic systems. Critical for NPSH calculations and avoiding cavitation damage.

### Access Methods
- **CoolProp**: `PropsSI('P', 'T', T, 'Q', 0, fluid)`
- **NIST Chemistry WebBook**: Vapor pressure tables
- **Antoine Equation**: Empirical correlation

### Python Implementation

```python
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

class CavitationDatabase:
    """Cavitation and vapor pressure database"""

    @staticmethod
    def antoine_equation(T, A, B, C):
        """
        Antoine equation for vapor pressure

        log10(P) = A - B / (T + C)

        P in mmHg, T in °C (original form)
        Returns P in Pa
        """
        T_celsius = T - 273.15
        log_P_mmHg = A - B / (T_celsius + C)
        P_mmHg = 10**log_P_mmHg
        P_Pa = P_mmHg * 133.322  # Convert mmHg to Pa
        return P_Pa

    @staticmethod
    def clausius_clapeyron(T, T_ref, P_ref, L_vap, R=8.314, M=0.018):
        """
        Clausius-Clapeyron equation for vapor pressure

        ln(P/P_ref) = -(L_vap*M/R) * (1/T - 1/T_ref)

        L_vap: Latent heat of vaporization (J/mol)
        M: Molar mass (kg/mol)
        """
        ln_P_ratio = -(L_vap * M / R) * (1/T - 1/T_ref)
        P = P_ref * np.exp(ln_P_ratio)
        return P

    def vapor_pressure(self, T, fluid='Water'):
        """Get vapor pressure using CoolProp"""
        try:
            P_vap = CP.PropsSI('P', 'T', T, 'Q', 0, fluid)
            return P_vap
        except:
            print(f"Error getting vapor pressure for {fluid} at {T}K")
            return None

    def NPSH_available(self, P_suction, P_vapor, rho, v_suction, z_suction=0, g=9.81):
        """
        Calculate NPSH available

        NPSH_a = (P_suction - P_vapor)/(ρ*g) + z_suction - v²/(2g)
        """
        NPSH_a = (P_suction - P_vapor) / (rho * g) + z_suction - v_suction**2 / (2 * g)
        return NPSH_a

    def cavitation_number(self, P_ref, P_vapor, rho, V):
        """
        Calculate cavitation number (sigma)

        σ = (P_ref - P_vapor) / (0.5 * ρ * V²)
        """
        sigma = (P_ref - P_vapor) / (0.5 * rho * V**2)
        return sigma

    def incipient_cavitation_pressure(self, P_vapor, rho, V, sigma_i=1.0):
        """
        Calculate pressure for incipient cavitation

        P_cav = P_vapor + σ_i * 0.5 * ρ * V²
        """
        P_cav = P_vapor + sigma_i * 0.5 * rho * V**2
        return P_cav

# Example usage
cav_db = CavitationDatabase()

# Example 1: Water vapor pressure vs temperature
T_range = np.linspace(0, 100, 50) + 273.15
P_vapor_water = [cav_db.vapor_pressure(T, 'Water') for T in T_range]

print("Water Vapor Pressure:")
print("T(°C)   P_vapor(kPa)")
for i in [0, 10, 20, 40]:
    idx = int(i / 2)
    print(f"{T_range[idx]-273.15:5.0f}      {P_vapor_water[idx]/1000:7.2f}")

# Example 2: NPSH available calculation
T_water = 60 + 273.15  # Hot water
P_suction = 101325  # Atmospheric pressure at pump inlet
rho = CP.PropsSI('D', 'T', T_water, 'P', P_suction, 'Water')
P_vapor = cav_db.vapor_pressure(T_water, 'Water')
v_suction = 2.0  # m/s in suction pipe

NPSH_a = cav_db.NPSH_available(P_suction, P_vapor, rho, v_suction)

print(f"\nNPSH Available for water at {T_water-273.15}°C:")
print(f"  Suction pressure: {P_suction/1000:.1f} kPa")
print(f"  Vapor pressure: {P_vapor/1000:.2f} kPa")
print(f"  Density: {rho:.1f} kg/m³")
print(f"  Suction velocity: {v_suction} m/s")
print(f"  NPSH_a: {NPSH_a:.2f} m")

# Example 3: Maximum water temperature for given NPSH
def max_temperature_for_NPSH(NPSH_required, P_suction, v_suction, z_suction=0):
    """Find maximum water temperature that provides required NPSH"""
    g = 9.81

    # Search for temperature
    for T in np.linspace(273.15, 373.15, 100):
        try:
            rho = CP.PropsSI('D', 'T', T, 'P', P_suction, 'Water')
            P_vapor = CP.PropsSI('P', 'T', T, 'Q', 0, 'Water')

            NPSH_a = (P_suction - P_vapor) / (rho * g) + z_suction - v_suction**2 / (2 * g)

            if NPSH_a < NPSH_required:
                return T - 273.15  # Return previous temperature
        except:
            continue

    return 100  # Maximum

NPSH_req = 4.0  # m
T_max = max_temperature_for_NPSH(NPSH_req, 101325, 2.0)
print(f"\nMaximum water temperature for NPSH_req = {NPSH_req}m:")
print(f"  T_max = {T_max:.1f}°C")

# Example 4: Cavitation number for valve
V_valve = 10  # m/s through valve
P_upstream = 500000  # Pa
T = 293.15  # K
rho = CP.PropsSI('D', 'T', T, 'P', P_upstream, 'Water')
P_vapor = cav_db.vapor_pressure(T, 'Water')

sigma = cav_db.cavitation_number(P_upstream, P_vapor, rho, V_valve)

print(f"\nCavitation Number for valve:")
print(f"  Upstream pressure: {P_upstream/1000:.0f} kPa")
print(f"  Velocity: {V_valve} m/s")
print(f"  Temperature: {T-273.15}°C")
print(f"  Cavitation number σ: {sigma:.2f}")
print(f"  Risk: {'Low' if sigma > 2 else 'Medium' if sigma > 1 else 'High'}")

# Example 5: Plot vapor pressure for multiple fluids
fluids = ['Water', 'Ethanol', 'Ammonia', 'R134a']
T_plot = np.linspace(-20, 100, 60) + 273.15

plt.figure(figsize=(12, 7))

for fluid in fluids:
    P_vap = []
    T_valid = []

    for T in T_plot:
        try:
            P = CP.PropsSI('P', 'T', T, 'Q', 0, fluid)
            if P > 0:  # Valid data point
                P_vap.append(P / 1000)  # Convert to kPa
                T_valid.append(T - 273.15)
        except:
            continue

    plt.semilogy(T_valid, P_vap, linewidth=2, marker='o', markersize=3, label=fluid)

# Add atmospheric pressure reference line
plt.axhline(y=101.325, color='black', linestyle='--', linewidth=1, label='Atmospheric pressure')

plt.xlabel('Temperature (°C)', fontsize=12)
plt.ylabel('Vapor Pressure (kPa)', fontsize=12)
plt.title('Vapor Pressure vs Temperature', fontsize=14, fontweight='bold')
plt.legend(fontsize=11)
plt.grid(True, which='both', alpha=0.3)
plt.savefig('vapor_pressure_comparison.png', dpi=300, bbox_inches='tight')
print("\n✓ Vapor pressure comparison saved")

# Example 6: Boiling point at different pressures
pressures = [50000, 101325, 200000, 500000]  # Pa

print("\nWater Boiling Point at Different Pressures:")
print("Pressure (kPa)  Boiling Point (°C)")
for P in pressures:
    try:
        T_boil = CP.PropsSI('T', 'P', P, 'Q', 0, 'Water')
        print(f"   {P/1000:7.1f}           {T_boil-273.15:8.2f}")
    except:
        print(f"   {P/1000:7.1f}           Error")
```

### Best Practices
1. Always check NPSH_available > NPSH_required + 0.5m margin
2. Account for elevation changes in suction line
3. Consider vapor pressure increase with temperature
4. Use conservative safety factors for critical applications
5. Monitor for cavitation noise and vibration in operation
6. Reduce suction velocities to minimize velocity head loss

---

## Turbulence Model Parameters

### Description
Parameters and coefficients for turbulence models used in CFD simulations including k-ε, k-ω, SST, and LES models.

### Access Methods
- **OpenFOAM**: Built-in model coefficients
- **ANSYS Fluent**: Model constants database
- **Literature**: Research papers and textbooks
- **CFD Online**: Community database

### Python Implementation

```python
import numpy as np

class TurbulenceModelDatabase:
    """Turbulence model parameters and coefficients"""

    def __init__(self):
        self.models = {
            'k_epsilon_standard': {
                'name': 'Standard k-ε Model',
                'reference': 'Launder & Spalding (1974)',
                'C_mu': 0.09,
                'C_1epsilon': 1.44,
                'C_2epsilon': 1.92,
                'sigma_k': 1.0,
                'sigma_epsilon': 1.3,
                'description': 'General purpose RANS model, good for free shear flows',
                'wall_function_required': True
            },
            'k_epsilon_realizable': {
                'name': 'Realizable k-ε Model',
                'reference': 'Shih et al. (1995)',
                'C_1epsilon': 1.44,
                'C_2epsilon': 1.9,
                'sigma_k': 1.0,
                'sigma_epsilon': 1.2,
                'C_mu': 'variable',  # Function of strain rate
                'description': 'Better for flows with strong streamline curvature and rotation',
                'wall_function_required': True
            },
            'k_epsilon_RNG': {
                'name': 'RNG k-ε Model',
                'reference': 'Yakhot & Orszag (1986)',
                'C_mu': 0.0845,
                'C_1epsilon': 1.42,
                'C_2epsilon': 1.68,
                'sigma_k': 0.7179,
                'sigma_epsilon': 0.7179,
                'eta_0': 4.38,
                'beta': 0.012,
                'description': 'Better for separated flows and swirl',
                'wall_function_required': True
            },
            'k_omega_standard': {
                'name': 'Standard k-ω Model',
                'reference': 'Wilcox (1988)',
                'beta_star': 0.09,
                'alpha': 5/9,
                'beta': 0.075,
                'sigma_k': 2.0,
                'sigma_omega': 2.0,
                'description': 'Good for wall-bounded flows, sensitive to freestream values',
                'wall_function_required': False
            },
            'k_omega_SST': {
                'name': 'k-ω SST Model',
                'reference': 'Menter (1994)',
                'beta_star': 0.09,
                # Inner region (k-ω)
                'alpha_1': 0.553,
                'beta_1': 0.075,
                'sigma_k1': 1.176,
                'sigma_omega1': 2.0,
                # Outer region (k-ε transformed to k-ω)
                'alpha_2': 0.44,
                'beta_2': 0.0828,
                'sigma_k2': 1.0,
                'sigma_omega2': 1.168,
                'a1': 0.31,  # Limiter for eddy viscosity
                'description': 'Industry standard, combines k-ω near wall with k-ε in free stream',
                'wall_function_required': False,
                'recommended_y_plus': '< 1 for accurate boundary layer'
            },
            'Spalart_Allmaras': {
                'name': 'Spalart-Allmaras Model',
                'reference': 'Spalart & Allmaras (1994)',
                'c_b1': 0.1355,
                'c_b2': 0.622,
                'sigma': 2/3,
                'c_v1': 7.1,
                'c_w2': 0.3,
                'c_w3': 2.0,
                'kappa': 0.41,  # von Karman constant
                'description': 'Single equation model, good for aerospace applications',
                'wall_function_required': False
            },
            'LES_Smagorinsky': {
                'name': 'Smagorinsky SGS Model',
                'reference': 'Smagorinsky (1963)',
                'C_s': 0.1,  # Smagorinsky constant (0.1-0.2 typical)
                'description': 'Basic SGS model for LES, simple but dissipative',
                'grid_requirement': 'Δx ~ Kolmogorov scale, expensive'
            },
            'DES': {
                'name': 'Detached Eddy Simulation',
                'reference': 'Spalart et al. (1997)',
                'C_DES': 0.65,
                'description': 'Hybrid RANS-LES, RANS near wall and LES in separated regions',
                'grid_requirement': 'Refined in regions where LES is desired'
            }
        }

    def get_model(self, model_key):
        """Get turbulence model parameters"""
        if model_key not in self.models:
            raise ValueError(f"Model '{model_key}' not found")
        return self.models[model_key]

    def calculate_turbulent_intensity(self, Re_D, intensity_type='pipe'):
        """
        Calculate turbulence intensity for boundary conditions

        intensity_type: 'pipe', 'low', 'medium', 'high'
        """
        if intensity_type == 'pipe':
            # Empirical correlation for pipe flow
            I = 0.16 * Re_D**(-1/8)
        elif intensity_type == 'low':
            I = 0.01  # 1%
        elif intensity_type == 'medium':
            I = 0.05  # 5%
        elif intensity_type == 'high':
            I = 0.10  # 10%
        else:
            raise ValueError("Unknown intensity type")

        return I

    def calculate_k_epsilon_inlet(self, U_inf, I, L_ref):
        """
        Calculate k and epsilon inlet conditions

        U_inf: Reference velocity (m/s)
        I: Turbulence intensity (0-1)
        L_ref: Reference length scale (m) - typically 0.07 * hydraulic diameter
        """
        k = 1.5 * (U_inf * I)**2
        C_mu = 0.09
        epsilon = C_mu**0.75 * k**1.5 / L_ref

        return k, epsilon

    def calculate_k_omega_inlet(self, U_inf, I, L_ref):
        """
        Calculate k and omega inlet conditions

        omega = epsilon / (C_mu * k)
        """
        k, epsilon = self.calculate_k_epsilon_inlet(U_inf, I, L_ref)
        C_mu = 0.09
        omega = epsilon / (C_mu * k)

        return k, omega

    def calculate_wall_yplus(self, rho, u_tau, y, mu):
        """
        Calculate y+ for wall-adjacent cell

        y+ = ρ * u_τ * y / μ

        u_τ: Friction velocity (m/s)
        y: Wall distance (m)
        """
        y_plus = rho * u_tau * y / mu
        return y_plus

    def estimate_first_cell_height(self, U_inf, L, nu, Re_L, y_plus_target=1):
        """
        Estimate first cell height for target y+

        For y+ < 1: Resolve viscous sublayer
        For y+ = 30-300: Use wall functions
        """
        # Estimate wall shear stress using flat plate correlation
        C_f = 0.026 / Re_L**0.143  # Turbulent flat plate
        tau_w = 0.5 * C_f * (U_inf**2 / nu) * nu  # ρ cancels
        u_tau = np.sqrt(tau_w * nu)  # Friction velocity

        # Calculate required wall distance
        y = y_plus_target * nu / u_tau

        return y, u_tau

    def print_model_summary(self):
        """Print summary of all available models"""
        print("Available Turbulence Models:")
        print("=" * 80)

        for key, model in self.models.items():
            print(f"\n{model['name']}")
            print(f"  Reference: {model['reference']}")
            print(f"  Description: {model['description']}")
            print(f"  Key: '{key}'")

# Example usage
turb_db = TurbulenceModelDatabase()

# Print all models
turb_db.print_model_summary()

# Example 1: Get k-ω SST parameters
sst_model = turb_db.get_model('k_omega_SST')
print(f"\n\nk-ω SST Model Parameters:")
for key, value in sst_model.items():
    if isinstance(value, float):
        print(f"  {key}: {value}")

# Example 2: Calculate inlet turbulence conditions
U_inlet = 5.0  # m/s
D_h = 0.1  # Hydraulic diameter (m)
nu = 1e-6  # Kinematic viscosity (m²/s)
Re_D = U_inlet * D_h / nu

I = turb_db.calculate_turbulent_intensity(Re_D, 'pipe')
L_ref = 0.07 * D_h  # Turbulence length scale

k, epsilon = turb_db.calculate_k_epsilon_inlet(U_inlet, I, L_ref)
k_omega, omega = turb_db.calculate_k_omega_inlet(U_inlet, I, L_ref)

print(f"\nInlet Turbulence Conditions:")
print(f"  Velocity: {U_inlet} m/s")
print(f"  Reynolds number: {Re_D:.0f}")
print(f"  Turbulence intensity: {I*100:.2f}%")
print(f"\n  k-ε model:")
print(f"    k = {k:.6f} m²/s²")
print(f"    ε = {epsilon:.6f} m²/s³")
print(f"\n  k-ω model:")
print(f"    k = {k_omega:.6f} m²/s²")
print(f"    ω = {omega:.4f} 1/s")

# Example 3: Calculate required mesh resolution for y+ = 1
L = 1.0  # Characteristic length (m)
Re_L = U_inlet * L / nu
y_1, u_tau = turb_db.estimate_first_cell_height(U_inlet, L, nu, Re_L, y_plus_target=1)

print(f"\nMesh Requirements for y+ = 1:")
print(f"  First cell height: {y_1*1e6:.2f} μm = {y_1*1e3:.4f} mm")
print(f"  Friction velocity: {u_tau:.4f} m/s")

# For y+ = 30 (wall functions)
y_30, _ = turb_db.estimate_first_cell_height(U_inlet, L, nu, Re_L, y_plus_target=30)
print(f"\nMesh Requirements for y+ = 30 (wall functions):")
print(f"  First cell height: {y_30*1e3:.3f} mm")

# Example 4: Model selection guide
print("\n\nTurbulence Model Selection Guide:")
print("=" * 80)
print("\nk-ε Standard:")
print("  ✓ Free shear flows, jets, wakes")
print("  ✓ External aerodynamics (with wall functions)")
print("  ✗ Flows with adverse pressure gradients")
print("  ✗ Rotating flows, swirl")

print("\nk-ω SST:")
print("  ✓ Boundary layers, adverse pressure gradients")
print("  ✓ Separated flows")
print("  ✓ Internal flows (pipes, pumps, valves)")
print("  ✓ General purpose - RECOMMENDED for most applications")
print("  ⚠ Requires fine mesh near walls (y+ < 1)")

print("\nSpalart-Allmaras:")
print("  ✓ Aerospace applications (airfoils, wings)")
print("  ✓ Attached boundary layers")
print("  ✓ Computationally efficient (1 equation)")
print("  ✗ Free shear flows")

print("\nLES:")
print("  ✓ Unsteady flows with large-scale structures")
print("  ✓ Separated flows, vortex shedding")
print("  ✓ Mixing, combustion")
print("  ⚠ VERY expensive computationally")
print("  ⚠ Requires very fine mesh (millions of cells)")
```

### Best Practices
1. Use k-ω SST as default for most engineering applications
2. Ensure y+ < 1 for k-ω SST (resolve boundary layer)
3. Use wall functions with k-ε if y+ = 30-300
4. Set appropriate turbulence intensity at inlets (1-10%)
5. Validate turbulence model choice with experimental data
6. Check mesh independence of results

---

## Hydraulic Components Database

### Description
Pressure loss coefficients (K-factors) for pipes, fittings, valves, and other hydraulic components.

### Access Methods
- **Crane TP-410**: Industry-standard handbook
- **ASHRAE Handbook**: HVAC applications
- **fluids Python package**: `fittings` module

### Python Implementation

```python
import numpy as np
from fluids import fittings

class HydraulicComponentsDatabase:
    """Database of K-factors for hydraulic components"""

    def __init__(self):
        # K-factors for various components (based on Crane TP-410)
        self.fittings = {
            # Elbows
            '90_elbow_standard': {'K': 30, 'description': '90° standard elbow'},
            '90_elbow_long_radius': {'K': 20, 'description': '90° long radius elbow'},
            '45_elbow': {'K': 16, 'description': '45° elbow'},

            # Tees
            'tee_branch_flow': {'K': 60, 'description': 'Tee, flow through branch'},
            'tee_line_flow': {'K': 20, 'description': 'Tee, flow through line'},

            # Valves (fully open)
            'gate_valve': {'K': 8, 'description': 'Gate valve, fully open'},
            'globe_valve': {'K': 340, 'description': 'Globe valve, fully open'},
            'ball_valve': {'K': 3, 'description': 'Ball valve, fully open'},
            'butterfly_valve': {'K': 45, 'description': 'Butterfly valve, fully open'},
            'check_valve_swing': {'K': 50, 'description': 'Swing check valve'},

            # Entrances/Exits
            'entrance_sharp': {'K': 0.5, 'description': 'Sharp-edged entrance'},
            'entrance_rounded': {'K': 0.04, 'description': 'Well-rounded entrance (r/d=0.15)'},
            'exit': {'K': 1.0, 'description': 'Exit to reservoir'},

            # Expansions/Contractions
            'sudden_expansion': {'K': 'variable', 'description': 'Sudden expansion (use area ratio)'},
            'sudden_contraction': {'K': 'variable', 'description': 'Sudden contraction (use area ratio)'},
        }

    def get_K_factor(self, component):
        """Get K-factor for component"""
        if component not in self.fittings:
            raise ValueError(f"Component '{component}' not found")

        K = self.fittings[component]['K']
        if K == 'variable':
            raise ValueError(f"{component} requires additional parameters")

        return K

    def pressure_loss(self, K, V, rho):
        """
        Calculate pressure loss through fitting

        ΔP = K * (ρ * V² / 2)

        K: Loss coefficient
        V: Velocity (m/s)
        rho: Density (kg/m³)
        """
        delta_P = K * 0.5 * rho * V**2
        return delta_P

    def head_loss(self, K, V, g=9.81):
        """
        Calculate head loss through fitting

        h_L = K * V² / (2g)
        """
        h_L = K * V**2 / (2 * g)
        return h_L

    def expansion_loss(self, D1, D2, V1):
        """
        Sudden expansion loss

        K = (1 - (D1/D2)²)²
        """
        K = (1 - (D1/D2)**2)**2
        h_L = K * V1**2 / (2 * 9.81)
        return h_L, K

    def contraction_loss(self, D1, D2, V2):
        """
        Sudden contraction loss

        K_c = 0.5 * (1 - (D2/D1)²)
        Based on downstream velocity V2
        """
        K = 0.5 * (1 - (D2/D1)**2)
        h_L = K * V2**2 / (2 * 9.81)
        return h_L, K

    def darcy_friction_loss(self, f, L, D, V, g=9.81):
        """
        Darcy-Weisbach friction loss

        h_f = f * (L/D) * V²/(2g)
        """
        h_f = f * (L / D) * V**2 / (2 * g)
        return h_f

    def equivalent_length(self, K, f):
        """
        Convert K-factor to equivalent pipe length

        L_eq = K * D / f
        """
        # This returns L_eq/D ratio
        return K / f

    def list_components(self):
        """List all available components"""
        print("Available Hydraulic Components:")
        print("=" * 80)
        for key, data in self.fittings.items():
            K_str = f"K = {data['K']}" if isinstance(data['K'], (int, float)) else "K = variable"
            print(f"  {key:25} : {data['description']:40} ({K_str})")

# Example usage
hyd_db = HydraulicComponentsDatabase()

# List components
hyd_db.list_components()

# Example 1: Calculate pressure loss through valve
V = 3.0  # m/s
rho = 998  # kg/m³
component = 'gate_valve'

K = hyd_db.get_K_factor(component)
delta_P = hyd_db.pressure_loss(K, V, rho)
h_L = hyd_db.head_loss(K, V)

print(f"\nPressure Loss Through {component}:")
print(f"  Velocity: {V} m/s")
print(f"  K-factor: {K}")
print(f"  ΔP: {delta_P/1000:.2f} kPa")
print(f"  Head loss: {h_L:.2f} m")

# Example 2: Total system loss calculation
def calculate_system_loss(components_list, V, D, L_pipe, f, rho=998, g=9.81):
    """
    Calculate total head loss in piping system

    components_list: List of component names
    V: Velocity (m/s)
    D: Pipe diameter (m)
    L_pipe: Total straight pipe length (m)
    f: Darcy friction factor
    """
    # Pipe friction loss
    h_friction = hyd_db.darcy_friction_loss(f, L_pipe, D, V)

    # Minor losses from fittings
    h_minor = 0
    K_total = 0

    for component in components_list:
        try:
            K = hyd_db.get_K_factor(component)
            h_L = hyd_db.head_loss(K, V)
            h_minor += h_L
            K_total += K
        except:
            print(f"Skipping {component}")

    h_total = h_friction + h_minor

    # Convert to equivalent length
    L_eq = hyd_db.equivalent_length(K_total, f) * D

    return {
        'h_friction': h_friction,
        'h_minor': h_minor,
        'h_total': h_total,
        'K_total': K_total,
        'L_equivalent': L_eq
    }

# Example piping system
components = [
    'entrance_rounded',
    '90_elbow_long_radius',
    '90_elbow_long_radius',
    'gate_valve',
    'tee_line_flow',
    'exit'
]

V = 2.5  # m/s
D = 0.1  # m
L_pipe = 50  # m
f = 0.02  # Friction factor

losses = calculate_system_loss(components, V, D, L_pipe, f)

print(f"\nPiping System Analysis:")
print(f"  Pipe: {D*1000}mm diameter, {L_pipe}m length")
print(f"  Velocity: {V} m/s")
print(f"  Friction factor: {f}")
print(f"\n  Friction loss: {losses['h_friction']:.2f} m")
print(f"  Minor losses: {losses['h_minor']:.2f} m")
print(f"  TOTAL HEAD LOSS: {losses['h_total']:.2f} m")
print(f"\n  Total K-factor: {losses['K_total']:.1f}")
print(f"  Equivalent length: {losses['L_equivalent']:.1f} m")
print(f"  Total equivalent: {L_pipe + losses['L_equivalent']:.1f} m")

# Example 3: Using fluids package
from fluids import friction_factor, Reynolds

# Calculate friction factor
Re = Reynolds(V=V, D=D, nu=1e-6)
roughness = 0.046e-3  # m (commercial steel)
f_calculated = friction_factor(Re=Re, eD=roughness/D)

print(f"\nFriction Factor Calculation:")
print(f"  Reynolds number: {Re:.0f}")
print(f"  Relative roughness: {roughness/D:.6f}")
print(f"  Friction factor (Colebrook): {f_calculated:.5f}")

# Example 4: Valve sizing - Cv coefficient
def calculate_Cv(Q_gpm, delta_P_psi, SG=1.0):
    """
    Calculate valve flow coefficient (Cv)

    Q_gpm: Flow rate (US gallons/minute)
    delta_P_psi: Pressure drop (psi)
    SG: Specific gravity
    """
    Cv = Q_gpm * np.sqrt(SG / delta_P_psi)
    return Cv

Q_gpm = 100  # gpm
delta_P_psi = 10  # psi

Cv = calculate_Cv(Q_gpm, delta_P_psi)
print(f"\nValve Sizing:")
print(f"  Flow rate: {Q_gpm} gpm")
print(f"  Pressure drop: {delta_P_psi} psi")
print(f"  Required Cv: {Cv:.1f}")
```

### Best Practices
1. Use K-factors from consistent source (Crane TP-410 standard)
2. Account for all fittings and components
3. Verify valve positions (open/closed/throttled)
4. Consider velocity changes in expansions/contractions
5. Use equivalent length method for hand calculations
6. Verify Reynolds number range for friction factor correlations

---

## Thermodynamic Tables

### Description
Steam tables, refrigerant properties, and thermodynamic property tables for common working fluids.

### Access Methods
- **IAPWS**: International standard for water/steam (via CoolProp)
- **ASHRAE Handbooks**: Refrigerant tables
- **CoolProp**: Comprehensive property library
- **pyromat**: Thermodynamic properties library

### Installation & Usage

```bash
# Install pyromat for thermodynamic properties
pip install pyromat

# CoolProp already installed from previous examples
```

### Python Implementation

```python
import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt

# Try to import pyromat (optional)
try:
    import pyromat as pm
    PYROMAT_AVAILABLE = True
except ImportError:
    PYROMAT_AVAILABLE = False
    print("pyromat not available, using CoolProp only")

class ThermodynamicTables:
    """Generate thermodynamic property tables"""

    @staticmethod
    def steam_table_sat_temp(T_range):
        """
        Generate saturated steam table (temperature basis)

        Returns: DataFrame with P_sat, v_f, v_g, h_f, h_fg, h_g, s_f, s_fg, s_g
        """
        import pandas as pd

        data = []

        for T in T_range:
            try:
                # Saturation pressure
                P_sat = CP.PropsSI('P', 'T', T, 'Q', 0, 'Water')

                # Liquid properties (Q=0)
                rho_f = CP.PropsSI('D', 'T', T, 'Q', 0, 'Water')
                h_f = CP.PropsSI('H', 'T', T, 'Q', 0, 'Water')
                s_f = CP.PropsSI('S', 'T', T, 'Q', 0, 'Water')

                # Vapor properties (Q=1)
                rho_g = CP.PropsSI('D', 'T', T, 'Q', 1, 'Water')
                h_g = CP.PropsSI('H', 'T', T, 'Q', 1, 'Water')
                s_g = CP.PropsSI('S', 'T', T, 'Q', 1, 'Water')

                # Specific volumes
                v_f = 1 / rho_f
                v_g = 1 / rho_g

                # Latent heat and entropy
                h_fg = h_g - h_f
                s_fg = s_g - s_f

                data.append({
                    'T(°C)': T - 273.15,
                    'P(kPa)': P_sat / 1000,
                    'v_f(m³/kg)': v_f,
                    'v_g(m³/kg)': v_g,
                    'h_f(kJ/kg)': h_f / 1000,
                    'h_fg(kJ/kg)': h_fg / 1000,
                    'h_g(kJ/kg)': h_g / 1000,
                    's_f(kJ/kg·K)': s_f / 1000,
                    's_fg(kJ/kg·K)': s_fg / 1000,
                    's_g(kJ/kg·K)': s_g / 1000
                })
            except:
                continue

        return pd.DataFrame(data)

    @staticmethod
    def superheated_steam_table(P, T_range):
        """
        Generate superheated steam table at constant pressure
        """
        import pandas as pd

        # Get saturation temperature at this pressure
        T_sat = CP.PropsSI('T', 'P', P, 'Q', 0, 'Water')

        data = []

        for T in T_range:
            if T <= T_sat:  # Skip if below saturation
                continue

            try:
                rho = CP.PropsSI('D', 'T', T, 'P', P, 'Water')
                h = CP.PropsSI('H', 'T', T, 'P', P, 'Water')
                s = CP.PropsSI('S', 'T', T, 'P', P, 'Water')
                cp = CP.PropsSI('C', 'T', T, 'P', P, 'Water')

                v = 1 / rho

                data.append({
                    'T(°C)': T - 273.15,
                    'v(m³/kg)': v,
                    'h(kJ/kg)': h / 1000,
                    's(kJ/kg·K)': s / 1000,
                    'cp(kJ/kg·K)': cp / 1000
                })
            except:
                continue

        return pd.DataFrame(data)

# Example 1: Generate saturated steam table
import pandas as pd
pd.set_option('display.float_format', '{:.4f}'.format)

T_sat_range = np.linspace(0, 200, 21) + 273.15  # 0-200°C
steam_table = ThermodynamicTables.steam_table_sat_temp(T_sat_range)

print("Saturated Steam Table (Water/Steam):")
print("=" * 100)
print(steam_table.to_string(index=False))

# Example 2: Superheated steam at 1 MPa
P_superheat = 1e6  # Pa (1 MPa)
T_superheat_range = np.linspace(200, 500, 16) + 273.15
superheat_table = ThermodynamicTables.superheated_steam_table(P_superheat, T_superheat_range)

print(f"\n\nSuperheated Steam at P = {P_superheat/1e6} MPa:")
print("=" * 80)
print(superheat_table.to_string(index=False))

# Example 3: Rankine cycle analysis using steam tables
def rankine_cycle_analysis(P_high, P_low, T_superheat=None, eta_turbine=0.85, eta_pump=0.75):
    """
    Analyze simple Rankine cycle

    P_high: Boiler pressure (Pa)
    P_low: Condenser pressure (Pa)
    T_superheat: Superheater exit temperature (K), None for saturated
    """
    print(f"\nRankine Cycle Analysis")
    print("=" * 80)

    # State 1: Saturated liquid leaving condenser
    T1 = CP.PropsSI('T', 'P', P_low, 'Q', 0, 'Water')
    h1 = CP.PropsSI('H', 'P', P_low, 'Q', 0, 'Water')
    s1 = CP.PropsSI('S', 'P', P_low, 'Q', 0, 'Water')

    print(f"\nState 1 (Condenser exit - saturated liquid):")
    print(f"  T = {T1-273.15:.2f}°C, P = {P_low/1000:.1f} kPa")
    print(f"  h = {h1/1000:.2f} kJ/kg, s = {s1/1000:.4f} kJ/(kg·K)")

    # State 2: After pump (isentropic)
    v1 = 1 / CP.PropsSI('D', 'P', P_low, 'Q', 0, 'Water')
    w_pump_s = v1 * (P_high - P_low)
    w_pump_actual = w_pump_s / eta_pump
    h2 = h1 + w_pump_actual

    T2 = CP.PropsSI('T', 'P', P_high, 'H', h2, 'Water')
    s2 = CP.PropsSI('S', 'P', P_high, 'H', h2, 'Water')

    print(f"\nState 2 (Pump exit):")
    print(f"  T = {T2-273.15:.2f}°C, P = {P_high/1000:.1f} kPa")
    print(f"  h = {h2/1000:.2f} kJ/kg, s = {s2/1000:.4f} kJ/(kg·K)")
    print(f"  Pump work = {w_pump_actual/1000:.2f} kJ/kg")

    # State 3: After boiler/superheater
    if T_superheat is None:
        # Saturated vapor
        h3 = CP.PropsSI('H', 'P', P_high, 'Q', 1, 'Water')
        s3 = CP.PropsSI('S', 'P', P_high, 'Q', 1, 'Water')
        T3 = CP.PropsSI('T', 'P', P_high, 'Q', 1, 'Water')
    else:
        # Superheated
        h3 = CP.PropsSI('H', 'T', T_superheat, 'P', P_high, 'Water')
        s3 = CP.PropsSI('S', 'T', T_superheat, 'P', P_high, 'Water')
        T3 = T_superheat

    q_in = h3 - h2

    print(f"\nState 3 (Boiler exit):")
    print(f"  T = {T3-273.15:.2f}°C, P = {P_high/1000:.1f} kPa")
    print(f"  h = {h3/1000:.2f} kJ/kg, s = {s3/1000:.4f} kJ/(kg·K)")
    print(f"  Heat input = {q_in/1000:.2f} kJ/kg")

    # State 4: After turbine (isentropic, then actual)
    h4s = CP.PropsSI('H', 'P', P_low, 'S', s3, 'Water')
    w_turbine_s = h3 - h4s
    w_turbine_actual = eta_turbine * w_turbine_s
    h4 = h3 - w_turbine_actual

    T4 = CP.PropsSI('T', 'P', P_low, 'H', h4, 'Water')
    s4 = CP.PropsSI('S', 'P', P_low, 'H', h4, 'Water')
    x4 = CP.PropsSI('Q', 'P', P_low, 'H', h4, 'Water')  # Quality

    print(f"\nState 4 (Turbine exit):")
    print(f"  T = {T4-273.15:.2f}°C, P = {P_low/1000:.1f} kPa")
    print(f"  h = {h4/1000:.2f} kJ/kg, s = {s4/1000:.4f} kJ/(kg·K)")
    print(f"  Quality x = {x4:.4f} ({'two-phase' if 0 < x4 < 1 else 'error'})")
    print(f"  Turbine work = {w_turbine_actual/1000:.2f} kJ/kg")

    # Condenser heat rejection
    q_out = h4 - h1

    print(f"\nCondenser heat rejection = {q_out/1000:.2f} kJ/kg")

    # Cycle performance
    w_net = w_turbine_actual - w_pump_actual
    eta_thermal = w_net / q_in

    # Carnot efficiency for comparison
    T_H = T3
    T_L = T1
    eta_carnot = 1 - T_L / T_H

    print(f"\n{'='*80}")
    print(f"Cycle Performance:")
    print(f"  Net work = {w_net/1000:.2f} kJ/kg")
    print(f"  Thermal efficiency = {eta_thermal*100:.2f}%")
    print(f"  Carnot efficiency = {eta_carnot*100:.2f}% (theoretical maximum)")
    print(f"  Back work ratio = {w_pump_actual/w_turbine_actual*100:.2f}%")

    # T-s diagram data
    return {
        'states': [(T1, s1), (T2, s2), (T3, s3), (T4, s4)],
        'eta': eta_thermal,
        'w_net': w_net
    }

# Run Rankine cycle analysis
P_boiler = 10e6  # 10 MPa
P_condenser = 10e3  # 10 kPa
T_superheat = 500 + 273.15  # 500°C

cycle_data = rankine_cycle_analysis(P_boiler, P_condenser, T_superheat)

print("\n✓ Thermodynamic analysis complete")
```

### Best Practices
1. Always specify phase (saturated/superheated/compressed)
2. Use quality (x) to determine two-phase region (0 < x < 1)
3. Verify units consistency (kJ vs J, kPa vs Pa)
4. Check critical point limits for properties
5. Use IAPWS-97 formulation for high accuracy (CoolProp default)

---

## Summary

This comprehensive engineering databases reference provides:

1. **Thermodynamic Properties**: REFPROP, CoolProp for accurate fluid properties
2. **Atmospheric Data**: NASA Earthdata for aerospace applications
3. **Material Properties**: Structural and thermal data for engineering materials
4. **Pump Performance**: Selection and analysis tools
5. **Fluid Viscosity**: Temperature-dependent correlations
6. **Cavitation Data**: NPSH and vapor pressure calculations
7. **Turbulence Models**: CFD simulation parameters
8. **Hydraulic Components**: K-factors and pressure losses
9. **Steam Tables**: Thermodynamic cycle analysis

### Quick Reference

| Database | Best For | License | Accuracy |
|----------|----------|---------|----------|
| NIST REFPROP | High-accuracy thermodynamics | Commercial | Excellent |
| CoolProp | General engineering | Free | Good |
| NASA Earthdata | Aerospace/atmosphere | Free | Very Good |
| Material DBs | Structural design | Mixed | Good |
| Pump DBs | Selection/analysis | Free | Moderate |
| Viscosity | Fluid mechanics | Free | Good |
| Turbulence | CFD | Free | Good |
| Hydraulic | System design | Free | Good |
| Steam Tables | Power cycles | Free | Excellent |

### Installation Summary

```bash
# Essential packages
pip install CoolProp fluids pint matplotlib numpy scipy pandas

# Optional advanced packages
pip install ctREFPROP  # Requires REFPROP license
pip install mp-api  # Materials Project
pip install matminer  # Materials data mining
pip install pyromat  # Thermodynamic properties
```

Use these databases and tools to perform accurate, efficient engineering calculations!
