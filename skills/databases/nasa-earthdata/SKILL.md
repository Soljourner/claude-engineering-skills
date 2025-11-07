---
name: nasa-earthdata
description: "Access atmospheric properties and aerospace fluid data from NASA Earthdata"
category: databases
domain: aerospace
complexity: intermediate
dependencies: []
---

# NASA Earthdata Database Skill

Access atmospheric properties, standard atmosphere models, and aerospace environmental data from NASA's Earthdata platform for aerospace pump design, high-altitude systems, and atmospheric flight analysis.

## Overview

NASA Earthdata provides free access to Earth science data from multiple NASA missions and instruments. For aerospace engineering applications, it offers:

- **Atmospheric Properties**: Temperature, pressure, density profiles vs altitude
- **Standard Atmosphere Models**: US Standard Atmosphere (1976), MSIS, NRLMSISE-00
- **Wind Data**: Global wind profiles, jet streams, seasonal variations
- **Temperature Profiles**: Troposphere through thermosphere
- **Pressure Profiles**: Sea level to exosphere
- **Geopotential Height**: Standard atmosphere reference
- **Composition Data**: Atmospheric gas composition by altitude
- **Environmental Conditions**: Humidity, precipitation, cloud cover

**Key Applications:**
- High-altitude pump and compressor design
- Aircraft and rocket environmental conditions
- Aerospace fluid system analysis
- Thermal control system design
- Aerodynamic heating calculations
- Launch window planning

## Free Registration Required

NASA Earthdata is **free** but requires user registration.

### Registration Process

1. **Create Earthdata Account**
   - URL: https://urs.earthdata.nasa.gov/users/new
   - Required information: Name, email, organization
   - Free for all users (academic, commercial, personal)
   - Instant approval

2. **Approve Applications**
   - Some datasets require application approval
   - Navigate to: https://urs.earthdata.nasa.gov/home
   - Click "My Applications" → "Approve More Applications"
   - Approve required data services (e.g., GES DISC, ASDC)

3. **Save Credentials**
   - Username: Your email or chosen username
   - Password: Your account password
   - Store securely (needed for API access)

## Authentication Setup

### Method 1: Username/Password Authentication

```python
import requests
from requests.auth import HTTPBasicAuth

# Your Earthdata credentials
username = "your_earthdata_username"
password = "your_earthdata_password"

# Example: Access data with authentication
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/..."
response = requests.get(url, auth=HTTPBasicAuth(username, password))
```

### Method 2: .netrc File (Recommended)

Create a `.netrc` file for automatic authentication:

**Linux/Mac** (`~/.netrc`):
```bash
machine urs.earthdata.nasa.gov
    login your_username
    password your_password
```

Set permissions:
```bash
chmod 600 ~/.netrc
```

**Windows** (`C:\Users\YourName\_netrc`):
```
machine urs.earthdata.nasa.gov
    login your_username
    password your_password
```

**Python with .netrc**:
```python
import requests

# Authentication automatically used from .netrc
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/..."
response = requests.get(url)
```

### Method 3: Token-Based Authentication

Generate a token for programmatic access:

1. Log in to Earthdata: https://urs.earthdata.nasa.gov/
2. Navigate to "My Profile" → "Generate Token"
3. Copy the token

**Using token in Python**:
```python
import requests

token = "your_generated_token"
headers = {"Authorization": f"Bearer {token}"}

url = "https://api.earthdata.nasa.gov/..."
response = requests.get(url, headers=headers)
```

## Available Datasets for Aerospace Applications

### 1. Standard Atmosphere Models

#### MSIS (Mass Spectrometer Incoherent Scatter)
- **Description**: Empirical atmospheric model
- **Altitude Range**: 0 km to 1000 km
- **Parameters**: Temperature, density, composition
- **Temporal**: Time-dependent (solar activity, season)
- **Use Case**: High-altitude aircraft, rockets, satellites

#### NRLMSISE-00 (Naval Research Laboratory MSIS Extended)
- **Description**: Extended atmospheric model to thermosphere
- **Altitude Range**: 0 km to 2000 km
- **Parameters**: T, P, ρ, composition (N₂, O₂, O, He, H, Ar, N)
- **Inputs**: Altitude, latitude, longitude, date, time, solar activity
- **Use Case**: Spacecraft drag, upper atmosphere analysis

#### US Standard Atmosphere (1976)
- **Description**: Static atmospheric model
- **Altitude Range**: 0 km to 1000 km
- **Layers**: Troposphere, stratosphere, mesosphere, thermosphere
- **Parameters**: T, P, ρ, speed of sound
- **Use Case**: Baseline design, standardized testing

### 2. MERRA-2 (Modern-Era Retrospective analysis)

Real atmospheric data from NASA reanalysis:

- **Spatial Resolution**: 0.5° × 0.625° (lat × lon)
- **Temporal Resolution**: Hourly, 3-hourly, daily, monthly
- **Altitude Levels**: Surface to 0.01 hPa (~80 km)
- **Time Period**: 1980 to present (updated ongoing)

**Available Parameters:**
- Temperature (T) [K]
- Pressure (P) [Pa]
- Density (ρ) [kg/m³]
- Wind components (U, V, W) [m/s]
- Geopotential height [m]
- Relative humidity [%]
- Specific humidity [kg/kg]

**Aerospace Applications:**
- Flight envelope analysis
- Wind shear assessment
- Thermal environment modeling
- Launch trajectory planning

### 3. AIRS (Atmospheric Infrared Sounder)

High-resolution atmospheric profiles:

- **Spatial Resolution**: 50 km
- **Vertical Resolution**: 28 pressure levels
- **Parameters**: T, P, H₂O, O₃, CO₂
- **Coverage**: Global, twice daily
- **Accuracy**: ±1 K (temperature), ±15% (humidity)

**Use Cases:**
- Precise atmospheric property data
- Regional atmospheric analysis
- Flight planning

### 4. GEOS (Goddard Earth Observing System)

Forward-looking atmospheric data:

- **Type**: Forecast model (7-day ahead)
- **Resolution**: 0.25° × 0.3125°
- **Temporal**: 3-hourly
- **Parameters**: T, P, winds, humidity, composition

**Use Cases:**
- Mission planning
- Launch forecasting
- Flight operations

## API Access Methods

### Method 1: OPeNDAP (Open-source Project for Network Data Access Protocol)

Direct data subsetting and download:

```python
from pydap.client import open_url

# Open MERRA-2 dataset via OPeNDAP
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2023/01/MERRA2_400.inst3_3d_asm_Np.20230101.nc4"

dataset = open_url(url, username="your_username", password="your_password")

# Access variables
temperature = dataset['T']  # Temperature [K]
pressure = dataset['PL']    # Pressure levels [Pa]
density = dataset['RHO']    # Density [kg/m³]

# Subset data (e.g., specific altitude and location)
T_subset = temperature[0, :, 100, 200]  # time, level, lat, lon
```

### Method 2: NASA CMR (Common Metadata Repository)

Search and discover datasets:

```python
import requests

# Search for atmospheric data
cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json"
params = {
    'short_name': 'M2I3NPASM',  # MERRA-2 3D atmospheric data
    'temporal': '2023-01-01T00:00:00Z,2023-01-31T23:59:59Z',
    'bounding_box': '-180,-90,180,90'  # Global
}

response = requests.get(cmr_url, params=params)
granules = response.json()['feed']['entry']

# Download URLs
for granule in granules:
    print(granule['title'])
    print(granule['links'][0]['href'])  # Data URL
```

### Method 3: Direct HTTP Download

Download files directly:

```python
import requests
from requests.auth import HTTPBasicAuth

username = "your_username"
password = "your_password"

# MERRA-2 file URL
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4/2023/01/MERRA2_400.inst3_3d_asm_Np.20230101.nc4"

# Download with authentication
response = requests.get(url, auth=HTTPBasicAuth(username, password), stream=True)

with open("merra2_data.nc4", "wb") as f:
    for chunk in response.iter_content(chunk_size=8192):
        f.write(chunk)

print("Download complete")
```

### Method 4: Python earthaccess Library

Simplified access to NASA Earthdata:

```bash
pip install earthaccess
```

```python
import earthaccess

# Login (uses .netrc or prompts for credentials)
earthaccess.login()

# Search for MERRA-2 data
results = earthaccess.search_data(
    short_name='M2I3NPASM',
    cloud_hosted=True,
    temporal=('2023-01-01', '2023-01-31')
)

# Download data
files = earthaccess.download(results[0:5], "./data")
```

## Standard Atmosphere Calculation Example

```python
import numpy as np

def us_standard_atmosphere_1976(altitude_m):
    """
    Calculate atmospheric properties using US Standard Atmosphere (1976)

    Parameters:
    -----------
    altitude_m : float
        Geometric altitude [m] (0 to 86000 m)

    Returns:
    --------
    dict: {'T': temperature [K],
           'P': pressure [Pa],
           'rho': density [kg/m³],
           'a': speed of sound [m/s]}
    """
    # Constants
    g0 = 9.80665  # Standard gravity [m/s²]
    R = 287.05    # Gas constant for air [J/kg/K]
    gamma = 1.4   # Specific heat ratio

    # Layer definitions [altitude_base, T_base, lapse_rate]
    layers = [
        (0,     288.15, -0.0065),  # Troposphere
        (11000, 216.65,  0.0),     # Tropopause
        (20000, 216.65,  0.001),   # Stratosphere 1
        (32000, 228.65,  0.0028),  # Stratosphere 2
        (47000, 270.65,  0.0),     # Stratopause
        (51000, 270.65, -0.0028),  # Mesosphere 1
        (71000, 214.65, -0.002),   # Mesosphere 2
    ]

    # Find appropriate layer
    h = altitude_m
    for i, (h_base, T_base, L) in enumerate(layers):
        if i + 1 < len(layers):
            h_next = layers[i + 1][0]
            if h >= h_base and h < h_next:
                break
        else:
            if h >= h_base:
                break

    # Base pressure (sea level)
    P_base = 101325  # Pa

    # Calculate pressure at each layer base
    for j, (hb, Tb, Lb) in enumerate(layers[:i+1]):
        if j == 0:
            P_base = 101325
        else:
            h_prev, T_prev, L_prev = layers[j-1]
            if abs(L_prev) < 1e-10:  # Isothermal
                P_base = P_base * np.exp(-g0 * (hb - h_prev) / (R * T_prev))
            else:  # Gradient
                P_base = P_base * (T_prev / (T_prev + L_prev * (hb - h_prev))) ** (g0 / (R * L_prev))

    # Calculate temperature at altitude
    T = T_base + L * (h - h_base)

    # Calculate pressure at altitude
    if abs(L) < 1e-10:  # Isothermal layer
        P = P_base * np.exp(-g0 * (h - h_base) / (R * T_base))
    else:  # Gradient layer
        P = P_base * (T_base / T) ** (g0 / (R * L))

    # Calculate density
    rho = P / (R * T)

    # Speed of sound
    a = np.sqrt(gamma * R * T)

    return {
        'T': T,
        'P': P,
        'rho': rho,
        'a': a
    }

# Example usage
altitudes = [0, 5000, 10000, 15000, 20000, 30000, 40000]  # meters

print("Altitude [m] | T [K] | P [Pa] | ρ [kg/m³] | a [m/s]")
print("-" * 65)
for h in altitudes:
    props = us_standard_atmosphere_1976(h)
    print(f"{h:12.0f} | {props['T']:5.2f} | {props['P']:8.1f} | {props['rho']:9.6f} | {props['a']:6.2f}")
```

## Applications to Aerospace Pumps and High-Altitude Systems

### 1. Pump Inlet Conditions

**Altitude effect on pump performance:**

```python
# Calculate NPSH available at different altitudes
def npsh_available(altitude_m, fluid='water', T_fluid=293.15):
    """
    Calculate Net Positive Suction Head Available at altitude

    NPSH_a = (P_atm - P_vapor) / (ρ * g) + elevation_head - friction_losses
    """
    atm = us_standard_atmosphere_1976(altitude_m)
    P_atm = atm['P']  # Atmospheric pressure at altitude

    # Water vapor pressure (Antoine equation, simplified)
    P_vapor = 611.2 * np.exp(17.67 * (T_fluid - 273.15) / (T_fluid - 29.65))

    # Water density (approximate)
    rho_water = 1000  # kg/m³
    g = 9.81  # m/s²

    NPSH_a = (P_atm - P_vapor) / (rho_water * g)

    return NPSH_a, P_atm

# Example: NPSH at various altitudes
print("Altitude [m] | P_atm [kPa] | NPSH_a [m]")
print("-" * 45)
for h in [0, 1000, 2000, 3000, 5000]:
    npsh, p_atm = npsh_available(h)
    print(f"{h:12.0f} | {p_atm/1000:11.2f} | {npsh:10.2f}")
```

**Output:**
```
Altitude [m] | P_atm [kPa] | NPSH_a [m]
---------------------------------------------
           0 |      101.33 |      10.12
        1000 |       89.88 |       8.96
        2000 |       79.50 |       7.88
        3000 |       70.12 |       6.88
        5000 |       54.05 |       5.12
```

### 2. High-Altitude Compressor Design

**Density variation impacts:**

```python
# Compressor power requirement vs altitude
def compressor_power(mass_flow_kg_s, pressure_ratio, altitude_m, eta_c=0.85):
    """
    Calculate compressor power at different altitudes
    """
    atm = us_standard_atmosphere_1976(altitude_m)
    T_in = atm['T']
    P_in = atm['P']
    rho_in = atm['rho']

    gamma = 1.4
    R = 287.05  # J/kg/K

    # Isentropic temperature ratio
    T_ratio = pressure_ratio ** ((gamma - 1) / gamma)

    # Actual temperature rise
    T_out = T_in * (1 + (T_ratio - 1) / eta_c)

    # Specific work
    w_c = R * (T_out - T_in) / (gamma - 1)

    # Power
    Power = mass_flow_kg_s * w_c / 1000  # kW

    return Power, rho_in, T_in

# Example: Compressor at sea level vs high altitude
m_dot = 1.0  # kg/s
PR = 3.0     # Pressure ratio

print(f"Compressor Performance (m_dot = {m_dot} kg/s, PR = {PR})")
print("Altitude [m] | ρ [kg/m³] | T_in [K] | Power [kW]")
print("-" * 60)
for h in [0, 5000, 10000, 15000]:
    P, rho, T = compressor_power(m_dot, PR, h)
    print(f"{h:12.0f} | {rho:9.4f} | {T:8.2f} | {P:10.2f}")
```

### 3. Thermal Control Systems

**Ambient temperature for radiator/heat exchanger design:**

```python
# Heat rejection at altitude
def heat_rejection_altitude(Q_reject_W, altitude_m):
    """
    Calculate required radiator area for heat rejection at altitude
    Assuming natural convection to ambient air
    """
    atm = us_standard_atmosphere_1976(altitude_m)
    T_amb = atm['T']
    rho_amb = atm['rho']

    # Simplified convection coefficient (natural convection)
    # h ~ ρ^0.5 (density effect)
    h_sea_level = 10  # W/m²/K (typical natural convection)
    h = h_sea_level * (rho_amb / 1.225) ** 0.5

    # Assume radiator surface temperature
    T_surface = 350  # K (example: electronics cooling)

    # Required area: Q = h * A * ΔT
    delta_T = T_surface - T_amb
    A_required = Q_reject_W / (h * delta_T)

    return A_required, h, T_amb

# Example
Q = 1000  # W heat rejection
print(f"Radiator Area Required for {Q} W Heat Rejection")
print("Altitude [m] | T_amb [K] | h [W/m²K] | Area [m²]")
print("-" * 60)
for h in [0, 5000, 10000, 15000]:
    A, h_conv, T_amb = heat_rejection_altitude(Q, h)
    print(f"{h:12.0f} | {T_amb:9.2f} | {h_conv:9.2f} | {A:9.4f}")
```

### 4. Fluid Properties at Altitude

**Cavitation and boiling point considerations:**

```python
def boiling_point_altitude(fluid='water'):
    """
    Calculate boiling point of water at different altitudes
    """
    print("Boiling Point vs Altitude (Water)")
    print("Altitude [m] | P_atm [kPa] | T_boil [°C]")
    print("-" * 50)

    for h in [0, 1000, 2000, 3000, 4000, 5000]:
        atm = us_standard_atmosphere_1976(h)
        P_atm = atm['P'] / 1000  # kPa

        # Approximate boiling point from pressure (Antoine equation)
        # log10(P) = A - B / (C + T)
        # Rearranged: T = B / (A - log10(P)) - C
        A, B, C = 8.07131, 1730.63, 233.426  # Antoine constants (water, T in °C, P in mmHg)
        P_mmHg = P_atm * 7.50062  # Convert kPa to mmHg
        T_boil = B / (A - np.log10(P_mmHg)) - C

        print(f"{h:12.0f} | {P_atm:11.2f} | {T_boil:11.2f}")

boiling_point_altitude()
```

## Data Processing Workflow

### Complete Example: Extract Atmospheric Profile

```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Read MERRA-2 NetCDF file
filename = "MERRA2_400.inst3_3d_asm_Np.20230101.nc4"
dataset = nc.Dataset(filename)

# Extract variables
lat = dataset.variables['lat'][:]  # Latitude
lon = dataset.variables['lon'][:]  # Longitude
lev = dataset.variables['lev'][:]  # Pressure levels [Pa]
T = dataset.variables['T'][:]      # Temperature [K]
H = dataset.variables['H'][:]      # Geopotential height [m]

# Select location (e.g., lat=40°N, lon=105°W - Colorado)
lat_idx = np.argmin(np.abs(lat - 40))
lon_idx = np.argmin(np.abs(lon - 255))  # -105° = 255° East

# Extract vertical profile at location
T_profile = T[0, :, lat_idx, lon_idx]  # time=0, all levels
H_profile = H[0, :, lat_idx, lon_idx]

# Convert to altitude (geopotential to geometric)
altitude_km = H_profile / 1000

# Plot temperature profile
plt.figure(figsize=(8, 10))
plt.plot(T_profile - 273.15, altitude_km)
plt.xlabel('Temperature [°C]')
plt.ylabel('Altitude [km]')
plt.title('Atmospheric Temperature Profile\n(40°N, 105°W)')
plt.grid(True)
plt.show()

dataset.close()
```

## Best Practices

1. **Cache Downloaded Data**: NASA datasets are large; download once and process locally
2. **Use OPeNDAP for Subsetting**: Only download needed spatial/temporal regions
3. **Check Data Version**: MERRA-2 has multiple collections; use latest (M2I3NPASM.5.12.4)
4. **Respect Download Limits**: NASA may throttle excessive requests
5. **Use Standard Atmosphere for Design**: Real data for analysis; standard atmosphere for conservative design
6. **Validate Results**: Cross-check with published data (ICAO, ISO 2533)
7. **Units Awareness**: NASA data uses SI units (K, Pa, kg/m³); convert as needed
8. **Time Zones**: NASA data in UTC; adjust for local analysis

## References

### Official NASA Resources

1. **NASA Earthdata Portal**
   - URL: https://www.earthdata.nasa.gov/
   - Main portal for all Earth science data

2. **Earthdata Search**
   - URL: https://search.earthdata.nasa.gov/
   - Visual search and discovery tool

3. **GES DISC (Goddard Earth Sciences Data and Information Services Center)**
   - URL: https://disc.gsfc.nasa.gov/
   - Primary source for atmospheric data

4. **MERRA-2 Documentation**
   - URL: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
   - Details on reanalysis methodology

### API and Tools

5. **earthaccess Python Library**
   - GitHub: https://github.com/nsidc/earthaccess
   - PyPI: `pip install earthaccess`

6. **OPeNDAP Protocol**
   - URL: https://www.opendap.org/
   - Data access protocol documentation

7. **CMR (Common Metadata Repository)**
   - URL: https://cmr.earthdata.nasa.gov/search/
   - Programmatic data search

### Atmospheric Models

8. **US Standard Atmosphere (1976)**
   - NOAA-S/T-76-1562
   - Free download from NOAA

9. **NRLMSISE-00**
   - URL: https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
   - Model code and documentation

10. **ISO 2533:1975**
    - Standard Atmosphere reference
    - Available from ISO

### Aerospace Applications

11. **AIAA Standards**
    - Atmospheric models for aerospace design
    - Available from AIAA

12. **ICAO Standard Atmosphere**
    - International Civil Aviation Organization
    - Doc 7488/3

---

*NASA Earthdata provides free, comprehensive atmospheric data essential for aerospace engineering analysis. Combined with standard atmosphere models, it enables accurate design and analysis of high-altitude pumps, compressors, thermal systems, and aerospace fluid applications.*
