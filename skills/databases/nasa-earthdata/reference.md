# NASA Earthdata Reference

Comprehensive reference for accessing atmospheric and aerospace environmental data from NASA Earthdata.

## Table of Contents

1. [Dataset Catalog](#dataset-catalog)
2. [Atmospheric Data Collections](#atmospheric-data-collections)
3. [API Documentation](#api-documentation)
4. [Query Examples](#query-examples)
5. [Data Formats](#data-formats)
6. [Standard Atmosphere Models](#standard-atmosphere-models)
7. [Coordinate Systems and Units](#coordinate-systems-and-units)
8. [External Resources](#external-resources)

---

## Dataset Catalog

### MERRA-2 (Modern-Era Retrospective analysis for Research and Applications, Version 2)

**Overview:**
- **Description**: NASA's latest atmospheric reanalysis dataset
- **Temporal Coverage**: 1980-present (ongoing updates)
- **Spatial Resolution**: 0.5° latitude × 0.625° longitude (~50 km)
- **Vertical Levels**: 72 hybrid-eta levels (surface to 0.01 hPa, ~80 km)
- **Temporal Resolution**: Hourly, 3-hourly, daily, monthly
- **Update Frequency**: 2-3 week lag

**Data Collections:**

| Collection | Description | Variables | Temporal Resolution |
|------------|-------------|-----------|---------------------|
| M2I3NPASM | 3D Atmospheric State | T, U, V, H, RH, QV | 3-hourly |
| M2T1NXSLV | Single-Level Diagnostics | Surface T, P, winds | Hourly |
| M2I1NXASM | 2D Atmospheric State | Column integrals | Hourly |
| M2T3NVASM | 3D Analyzed State | Full atmosphere | 3-hourly |
| M2I3NVAER | 3D Aerosol | Dust, sea salt, smoke | 3-hourly |

**Key Variables for Aerospace:**

| Variable | Symbol | Units | Description |
|----------|--------|-------|-------------|
| T | T | K | Air temperature |
| U | u | m/s | Eastward wind component |
| V | v | m/s | Northward wind component |
| OMEGA | ω | Pa/s | Vertical pressure velocity |
| H | H | m | Geopotential height |
| RH | RH | % | Relative humidity |
| QV | q | kg/kg | Specific humidity |
| PL | p | Pa | Pressure level |
| PS | Ps | Pa | Surface pressure |

**Access URL:**
```
https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/[COLLECTION]/[YEAR]/[MONTH]/
```

**File Naming Convention:**
```
MERRA2_[stream].[collection].[YYYYMMDD].nc4

Example:
MERRA2_400.inst3_3d_asm_Np.20230101.nc4
```

---

### AIRS (Atmospheric Infrared Sounder)

**Overview:**
- **Description**: High-resolution atmospheric sounding from satellite
- **Temporal Coverage**: 2002-present
- **Spatial Resolution**: 50 km at nadir
- **Vertical Levels**: 28 standard pressure levels
- **Temporal Resolution**: Twice daily (ascending/descending passes)

**Pressure Levels [hPa]:**
```
1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100,
70, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5, 1, 0.5
```

**Key Variables:**

| Variable | Units | Description | Accuracy |
|----------|-------|-------------|----------|
| Temperature | K | Air temperature profile | ±1 K |
| Water Vapor | kg/kg | Specific humidity | ±15% |
| Ozone | mol/mol | O₃ mixing ratio | ±10% |
| CO₂ | ppmv | Carbon dioxide | ±2 ppmv |
| Geopotential Height | m | Altitude of pressure surfaces | ±30 m |

**Access URL:**
```
https://acdisc.gesdisc.eosdis.nasa.gov/data/Aqua_AIRS_Level3/AIRX3STD.006/
```

---

### GEOS-FP (Goddard Earth Observing System - Forward Processing)

**Overview:**
- **Description**: Near-real-time atmospheric analysis and forecast
- **Temporal Coverage**: Current - 7 days ahead
- **Spatial Resolution**: 0.25° × 0.3125° (~25 km)
- **Vertical Levels**: 72 levels (surface to 0.01 hPa)
- **Temporal Resolution**: 3-hourly analysis, 3-hourly forecasts
- **Update Frequency**: Daily

**Forecast Products:**

| Product | Description | Forecast Range |
|---------|-------------|----------------|
| assim | Assimilated analysis | Current time |
| asm | Analysis product | Current time |
| fcst | Forecast | +6, +12, +18, +24 hours |
| tavg3 | 3-hour time average | Forecast period |

**Use Cases:**
- Mission planning
- Launch window analysis
- Flight operations
- Near-term atmospheric forecasting

**Access URL:**
```
https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/
```

---

### CALIPSO (Cloud-Aerosol Lidar and Infrared Pathfinder Satellite Observation)

**Overview:**
- **Description**: Vertical structure of clouds and aerosols
- **Temporal Coverage**: 2006-present
- **Vertical Resolution**: 30-60 m (below 20 km)
- **Spatial Resolution**: 333 m horizontal, 60 m vertical
- **Measurement**: Lidar backscatter at 532 nm and 1064 nm

**Aerospace Applications:**
- Cloud top heights
- Aerosol layers
- Volcanic ash detection
- Flight hazard assessment

---

### Standard Atmosphere Datasets

#### MSISE-00 (Mass Spectrometer and Incoherent Scatter Extended Model)

**Overview:**
- **Altitude Range**: 0-1000 km
- **Variables**: T, ρ, composition (N₂, O₂, O, He, H, Ar, N, O*, H₂)
- **Inputs**: Altitude, latitude, longitude, day of year, UT, F10.7, Ap
- **Update**: Empirical model based on satellite and ground-based data

**Model Inputs:**

| Input | Description | Range | Units |
|-------|-------------|-------|-------|
| alt | Altitude | 0-1000 | km |
| lat | Latitude | -90 to 90 | degrees |
| lon | Longitude | -180 to 180 | degrees |
| doy | Day of year | 1-366 | - |
| sec | Seconds in day | 0-86400 | s |
| F107 | Solar radio flux | ~70-300 | sfu |
| F107a | 81-day avg F10.7 | ~70-300 | sfu |
| Ap | Geomagnetic index | 0-400 | - |

**Python Implementation:**
```python
from pyglow import pyglow

# Create atmospheric point
pt = pyglow.Point(
    dn=datetime.datetime(2023, 1, 1, 12, 0, 0),  # Date and time
    lat=40,    # Latitude [deg]
    lon=-105,  # Longitude [deg]
    alt=100    # Altitude [km]
)

# Run MSISE-00
pt.run_msise()

# Access results
temperature = pt.Tn_msis  # Neutral temperature [K]
density = pt.nn['N2']     # N₂ number density [m⁻³]
total_density = pt.rho    # Total mass density [kg/m³]
```

#### US Standard Atmosphere (1976)

**Layers:**

| Layer | Altitude Range [km] | Lapse Rate [K/km] | Base Temp [K] |
|-------|---------------------|-------------------|---------------|
| Troposphere | 0 - 11 | -6.5 | 288.15 |
| Tropopause | 11 - 20 | 0 | 216.65 |
| Stratosphere 1 | 20 - 32 | +1.0 | 216.65 |
| Stratosphere 2 | 32 - 47 | +2.8 | 228.65 |
| Stratopause | 47 - 51 | 0 | 270.65 |
| Mesosphere 1 | 51 - 71 | -2.8 | 270.65 |
| Mesosphere 2 | 71 - 86 | -2.0 | 214.65 |

**Reference Properties (Sea Level):**
- Temperature: 288.15 K (15°C)
- Pressure: 101325 Pa (1 atm)
- Density: 1.225 kg/m³
- Speed of sound: 340.29 m/s

---

## Atmospheric Data Collections

### Data Centers

#### GES DISC (Goddard Earth Sciences Data and Information Services Center)

**Primary for atmospheric data**

**Collections:**
- MERRA-2 (reanalysis)
- AIRS (satellite sounding)
- GEOS-5 (model)
- GPM (precipitation)
- TRMM (tropical rainfall)

**Base URL:**
```
https://disc.gsfc.nasa.gov/
```

#### ASDC (Atmospheric Science Data Center)

**Primary for radiation and clouds**

**Collections:**
- CERES (radiation budget)
- CALIPSO (lidar)
- MISR (multi-angle imaging)
- TES (tropospheric emission)

**Base URL:**
```
https://asdc.larc.nasa.gov/
```

#### LAADS DAAC (Level-1 and Atmosphere Archive & Distribution System)

**Primary for aerosol and cloud properties**

**Collections:**
- MODIS (aerosol, cloud products)
- VIIRS (visible/infrared imaging)
- MISR (aerosol)

**Base URL:**
```
https://ladsweb.modaps.eosdis.nasa.gov/
```

---

## API Documentation

### 1. OPeNDAP Access

**OPeNDAP URL Structure:**
```
https://[data_center]/opendap/[collection]/[version]/[year]/[month]/[filename].nc4
```

**Example:**
```
https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2023/01/MERRA2_400.inst3_3d_asm_Np.20230101.nc4
```

**Python Access:**

```python
from pydap.client import open_url
from pydap.cas.urs import setup_session

# Setup authentication
username = "your_earthdata_username"
password = "your_earthdata_password"
session = setup_session(username, password, check_url=url)

# Open dataset
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2023/01/MERRA2_400.inst3_3d_asm_Np.20230101.nc4"
dataset = open_url(url, session=session)

# List variables
print(dataset.keys())

# Access specific variable
temperature = dataset['T']  # 4D array: [time, lev, lat, lon]

# Get metadata
print(temperature.attributes)

# Subset data (OPeNDAP does server-side subsetting)
# Syntax: variable[time_start:time_end:time_stride, lev_start:lev_end, lat_start:lat_end, lon_start:lon_end]
T_subset = temperature[0, 20:30, 100:110, 200:210]  # Specific region

# Convert to numpy array
T_array = T_subset[:]
```

**OPeNDAP Subsetting Advantages:**
- Server-side data reduction (download only what you need)
- Faster access for small regions
- Reduced bandwidth and storage

---

### 2. NASA CMR (Common Metadata Repository) API

**Search for Data Granules:**

**CMR Search Endpoint:**
```
https://cmr.earthdata.nasa.gov/search/granules.json
```

**Query Parameters:**

| Parameter | Description | Example |
|-----------|-------------|---------|
| `short_name` | Dataset short name | `M2I3NPASM` |
| `version` | Dataset version | `5.12.4` |
| `temporal` | Time range | `2023-01-01T00:00:00Z,2023-01-31T23:59:59Z` |
| `bounding_box` | Spatial extent | `west,south,east,north` |
| `page_size` | Results per page | `100` |

**Python Example:**

```python
import requests

def search_earthdata(short_name, temporal_range, bbox=None):
    """
    Search NASA CMR for data granules

    Parameters:
    -----------
    short_name : str
        Dataset short name (e.g., 'M2I3NPASM')
    temporal_range : tuple
        (start_date, end_date) in ISO format
    bbox : tuple, optional
        (west, south, east, north) bounding box

    Returns:
    --------
    list : List of data granules with download URLs
    """
    cmr_url = "https://cmr.earthdata.nasa.gov/search/granules.json"

    params = {
        'short_name': short_name,
        'temporal': f"{temporal_range[0]},{temporal_range[1]}",
        'page_size': 100
    }

    if bbox:
        params['bounding_box'] = f"{bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"

    response = requests.get(cmr_url, params=params)
    granules = response.json()['feed']['entry']

    results = []
    for granule in granules:
        result = {
            'title': granule['title'],
            'time': granule['time_start'],
            'links': [link['href'] for link in granule['links'] if link.get('rel') == 'http://esipfed.org/ns/fedsearch/1.1/data#']
        }
        results.append(result)

    return results

# Example usage
granules = search_earthdata(
    short_name='M2I3NPASM',
    temporal_range=('2023-01-01T00:00:00Z', '2023-01-31T23:59:59Z'),
    bbox=(-110, 35, -100, 45)  # Colorado region
)

for g in granules[:5]:  # First 5 results
    print(f"Date: {g['time']}")
    print(f"Title: {g['title']}")
    print(f"URL: {g['links'][0]}\n")
```

---

### 3. Direct HTTP/HTTPS Download

**Authentication Required:**

```python
import requests
from requests.auth import HTTPBasicAuth

def download_earthdata_file(url, output_filename, username, password):
    """
    Download file from NASA Earthdata with authentication

    Parameters:
    -----------
    url : str
        Direct URL to data file
    output_filename : str
        Local filename to save
    username : str
        Earthdata username
    password : str
        Earthdata password
    """
    # Create session to handle redirects
    session = requests.Session()
    session.auth = (username, password)

    # Request file
    response = session.get(url, stream=True)

    # Check if authentication succeeded
    if response.status_code == 401:
        raise ValueError("Authentication failed. Check username and password.")

    # Download with progress
    total_size = int(response.headers.get('content-length', 0))
    block_size = 8192
    downloaded = 0

    with open(output_filename, 'wb') as f:
        for chunk in response.iter_content(chunk_size=block_size):
            if chunk:
                f.write(chunk)
                downloaded += len(chunk)
                progress = (downloaded / total_size) * 100 if total_size > 0 else 0
                print(f"\rProgress: {progress:.1f}%", end='')

    print(f"\nDownload complete: {output_filename}")

# Example usage
url = "https://goldsmr4.gesdisc.eosdis.nasa.gov/data/MERRA2/M2I3NPASM.5.12.4/2023/01/MERRA2_400.inst3_3d_asm_Np.20230101.nc4"
download_earthdata_file(url, "merra2_data.nc4", "your_username", "your_password")
```

---

### 4. earthaccess Library (Recommended)

**Installation:**
```bash
pip install earthaccess
```

**Features:**
- Simplified authentication
- Unified search interface
- Cloud-optimized access (AWS S3)
- Automatic session management

**Python Example:**

```python
import earthaccess
import xarray as xr

# Login (interactive prompt or uses .netrc)
earthaccess.login()

# Search for data
results = earthaccess.search_data(
    short_name='M2I3NPASM',
    cloud_hosted=True,
    temporal=('2023-01-01', '2023-01-31'),
    bounding_box=(-110, 35, -100, 45)
)

print(f"Found {len(results)} granules")

# Download files
files = earthaccess.download(
    results[0:5],  # Download first 5 files
    local_path='./data'
)

# Or open directly in cloud (if on AWS)
datasets = earthaccess.open(results[0:5])

# Open with xarray
ds = xr.open_mfdataset(files, engine='netcdf4')

# Access variables
temperature = ds['T']
pressure = ds['lev']
```

---

## Query Examples

### Example 1: Atmospheric Profile at Specific Location

**Extract temperature and pressure profile:**

```python
import netCDF4 as nc
import numpy as np

def get_atmospheric_profile(filename, lat_target, lon_target):
    """
    Extract atmospheric profile at specified location

    Parameters:
    -----------
    filename : str
        Path to MERRA-2 NetCDF file
    lat_target : float
        Target latitude [-90, 90]
    lon_target : float
        Target longitude [-180, 180]

    Returns:
    --------
    dict : {'altitude': altitude [km],
            'temperature': T [K],
            'pressure': P [Pa],
            'density': rho [kg/m³]}
    """
    # Open dataset
    ds = nc.Dataset(filename)

    # Get coordinates
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    lev = ds.variables['lev'][:]  # Pressure levels [Pa]

    # Find nearest grid point
    lat_idx = np.argmin(np.abs(lat - lat_target))
    lon_idx = np.argmin(np.abs(lon - (lon_target % 360)))  # Convert to 0-360

    # Extract vertical profile (time=0)
    T_profile = ds.variables['T'][0, :, lat_idx, lon_idx]  # Temperature [K]
    H_profile = ds.variables['H'][0, :, lat_idx, lon_idx]  # Geopotential height [m]
    QV_profile = ds.variables['QV'][0, :, lat_idx, lon_idx]  # Specific humidity [kg/kg]

    # Calculate density from ideal gas law
    R_specific = 287.05  # J/kg/K (dry air)
    rho_profile = lev / (R_specific * T_profile)

    # Convert geopotential height to altitude
    altitude_km = H_profile / 1000

    ds.close()

    return {
        'altitude': altitude_km,
        'temperature': T_profile,
        'pressure': lev,
        'density': rho_profile,
        'specific_humidity': QV_profile
    }

# Example usage
profile = get_atmospheric_profile(
    'MERRA2_400.inst3_3d_asm_Np.20230101.nc4',
    lat_target=40.0,   # 40°N
    lon_target=-105.0  # 105°W (Colorado)
)

# Print results
print("Altitude [km] | Temp [K] | Pressure [Pa] | Density [kg/m³]")
print("-" * 65)
for i in range(0, len(profile['altitude']), 5):  # Every 5 levels
    print(f"{profile['altitude'][i]:13.2f} | {profile['temperature'][i]:8.2f} | "
          f"{profile['pressure'][i]:13.1f} | {profile['density'][i]:15.6f}")
```

### Example 2: Wind Profile for Flight Planning

**Extract wind components at cruise altitude:**

```python
def get_wind_profile(filename, lat_target, lon_target, pressure_level_pa=25000):
    """
    Extract wind data at cruise altitude (~250 hPa ≈ 10 km)

    Returns:
    --------
    dict : {'U': eastward wind [m/s],
            'V': northward wind [m/s],
            'speed': wind speed [m/s],
            'direction': wind direction [degrees from N]}
    """
    ds = nc.Dataset(filename)

    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    lev = ds.variables['lev'][:]

    # Find indices
    lat_idx = np.argmin(np.abs(lat - lat_target))
    lon_idx = np.argmin(np.abs(lon - (lon_target % 360)))
    lev_idx = np.argmin(np.abs(lev - pressure_level_pa))

    # Extract wind components
    U = ds.variables['U'][0, lev_idx, lat_idx, lon_idx]  # Eastward [m/s]
    V = ds.variables['V'][0, lev_idx, lat_idx, lon_idx]  # Northward [m/s]

    # Calculate magnitude and direction
    speed = np.sqrt(U**2 + V**2)
    direction = np.degrees(np.arctan2(U, V)) % 360  # Direction FROM

    ds.close()

    return {
        'U': U,
        'V': V,
        'speed': speed,
        'direction': direction,
        'pressure': lev[lev_idx],
        'lat': lat[lat_idx],
        'lon': lon[lon_idx]
    }

# Example usage
wind = get_wind_profile(
    'MERRA2_400.inst3_3d_asm_Np.20230101.nc4',
    lat_target=40.0,
    lon_target=-105.0,
    pressure_level_pa=25000  # ~250 hPa = ~34,000 ft
)

print(f"Wind at {wind['pressure']/100:.0f} hPa ({wind['lat']:.1f}°N, {wind['lon']:.1f}°E):")
print(f"  U (eastward): {wind['U']:6.2f} m/s")
print(f"  V (northward): {wind['V']:6.2f} m/s")
print(f"  Speed: {wind['speed']:6.2f} m/s ({wind['speed']*1.944:.1f} knots)")
print(f"  Direction: {wind['direction']:6.1f}° (from)")
```

### Example 3: Temporal Variation Analysis

**Analyze daily temperature variation at altitude:**

```python
def analyze_daily_variation(file_list, lat_target, lon_target, altitude_km):
    """
    Analyze temperature variation over multiple time steps

    Parameters:
    -----------
    file_list : list
        List of MERRA-2 files (e.g., hourly or 3-hourly)
    altitude_km : float
        Target altitude [km]

    Returns:
    --------
    dict : {'time': time array,
            'temperature': T array [K]}
    """
    temperatures = []
    times = []

    for filename in file_list:
        ds = nc.Dataset(filename)

        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
        H = ds.variables['H'][:]  # Geopotential height [m]

        # Find location
        lat_idx = np.argmin(np.abs(lat - lat_target))
        lon_idx = np.argmin(np.abs(lon - (lon_target % 360)))

        # Find altitude level
        H_profile = H[0, :, lat_idx, lon_idx] / 1000  # Convert to km
        alt_idx = np.argmin(np.abs(H_profile - altitude_km))

        # Extract temperature at this level
        T = ds.variables['T'][0, alt_idx, lat_idx, lon_idx]

        # Get time
        time = ds.variables['time'][0]

        temperatures.append(T)
        times.append(time)

        ds.close()

    return {
        'time': np.array(times),
        'temperature': np.array(temperatures)
    }
```

### Example 4: Spatial Map at Constant Altitude

**Create temperature map at flight level:**

```python
def create_altitude_map(filename, altitude_km):
    """
    Create 2D temperature map at constant altitude

    Returns:
    --------
    dict : {'lat': latitude array,
            'lon': longitude array,
            'temperature': 2D temperature array [K]}
    """
    ds = nc.Dataset(filename)

    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    H = ds.variables['H'][0, :, :, :]  # [lev, lat, lon]

    # For each grid point, find closest altitude level
    T_map = np.zeros((len(lat), len(lon)))

    for i in range(len(lat)):
        for j in range(len(lon)):
            H_column = H[:, i, j] / 1000  # km
            alt_idx = np.argmin(np.abs(H_column - altitude_km))
            T_map[i, j] = ds.variables['T'][0, alt_idx, i, j]

    ds.close()

    return {
        'lat': lat,
        'lon': lon,
        'temperature': T_map
    }

# Plot with matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

data = create_altitude_map('MERRA2_400.inst3_3d_asm_Np.20230101.nc4', altitude_km=10)

plt.figure(figsize=(12, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
cs = ax.contourf(data['lon'], data['lat'], data['temperature'],
                 levels=20, transform=ccrs.PlateCarree(), cmap='RdYlBu_r')
ax.coastlines()
plt.colorbar(cs, label='Temperature [K]')
plt.title('Temperature at 10 km Altitude')
plt.show()
```

---

## Data Formats

### NetCDF-4 Format

**Structure:**

NASA Earthdata primarily uses **NetCDF-4** (Network Common Data Form) format.

**File Extension:** `.nc4`

**Key Features:**
- Self-describing (metadata embedded)
- Binary format (efficient storage)
- Multi-dimensional arrays
- Compression support
- CF (Climate and Forecast) conventions

**Dimensions:**

| Dimension | Description | Typical Size |
|-----------|-------------|--------------|
| `time` | Time coordinate | 1-8 (per file) |
| `lev` | Vertical level (pressure) | 42-72 |
| `lat` | Latitude | 361 (0.5°) or 721 (0.25°) |
| `lon` | Longitude | 576 (0.625°) or 1152 (0.3125°) |

**Variable Structure:**
```python
# Variable dimensions: [time, lev, lat, lon]
# Example: T[0, 20, 100, 200]
#   - time index 0 (first time step)
#   - level index 20 (~300 hPa)
#   - lat index 100 (~10°N)
#   - lon index 200 (~125°E)
```

**Metadata (Attributes):**
```python
import netCDF4 as nc

ds = nc.Dataset('file.nc4')
var = ds.variables['T']

# Variable attributes
print(var.long_name)    # "air_temperature"
print(var.units)        # "K"
print(var.standard_name) # CF standard name
print(var.missing_value) # Fill value for missing data
print(var.valid_range)   # [min, max] valid values

# Global attributes
print(ds.Title)         # Dataset title
print(ds.Institution)   # "NASA GSFC"
print(ds.History)       # Processing history
```

**Reading NetCDF Files:**

```python
import netCDF4 as nc
import numpy as np

# Open file
ds = nc.Dataset('MERRA2_400.inst3_3d_asm_Np.20230101.nc4', 'r')

# List all variables
print("Variables:", ds.variables.keys())

# List dimensions
print("Dimensions:", ds.dimensions.keys())

# Get dimension sizes
print(f"Time steps: {len(ds.dimensions['time'])}")
print(f"Levels: {len(ds.dimensions['lev'])}")
print(f"Latitude: {len(ds.dimensions['lat'])}")
print(f"Longitude: {len(ds.dimensions['lon'])}")

# Read entire variable (caution: can be large!)
temperature = ds.variables['T'][:]  # Returns numpy array

# Read subset (recommended)
T_subset = ds.variables['T'][0, 20:30, 100:110, 200:210]

# Close file
ds.close()
```

**Using xarray (Recommended):**

```python
import xarray as xr

# Open dataset (lazily loaded)
ds = xr.open_dataset('MERRA2_400.inst3_3d_asm_Np.20230101.nc4')

# View dataset info
print(ds)

# Access variables with metadata
temperature = ds['T']
print(temperature.attrs)  # Metadata
print(temperature.dims)   # Dimensions

# Subset using coordinate values (not indices)
T_subset = ds['T'].sel(
    time='2023-01-01T00:00:00',
    lev=25000,           # 250 hPa
    lat=40,              # 40°N
    lon=255,             # 105°W (255°E)
    method='nearest'     # Nearest neighbor
)

# Arithmetic operations
T_celsius = temperature - 273.15

# Aggregations
T_mean = temperature.mean(dim='time')
T_std = temperature.std(dim='time')

# Close
ds.close()
```

---

### HDF-EOS Format

Some NASA datasets use **HDF-EOS** (Hierarchical Data Format - Earth Observing System).

**File Extension:** `.hdf` or `.he5`

**Reading HDF-EOS:**

```python
from pyhdf.SD import SD, SDC

# Open HDF file
hdf = SD('file.hdf', SDC.READ)

# List datasets
print(hdf.datasets())

# Read dataset
data = hdf.select('Temperature')
temperature = data.get()

# Get attributes
attrs = data.attributes()

hdf.end()
```

---

### ASCII/CSV Format (Derived Products)

Some processed data available as ASCII:

**Example: Radiosonde Profile**
```
# Altitude [km], Temperature [K], Pressure [Pa], Density [kg/m³]
0.0, 288.15, 101325, 1.225
1.0, 281.65, 89875, 1.112
2.0, 275.15, 79495, 1.007
5.0, 255.65, 54048, 0.736
10.0, 223.25, 26436, 0.412
```

---

## Standard Atmosphere Models

### US Standard Atmosphere (1976)

**Reference Properties:**

| Altitude [km] | Temp [K] | Pressure [Pa] | Density [kg/m³] | Speed of Sound [m/s] |
|---------------|----------|---------------|-----------------|----------------------|
| 0 | 288.15 | 101325 | 1.2250 | 340.29 |
| 5 | 255.68 | 54048 | 0.7361 | 320.53 |
| 10 | 223.25 | 26436 | 0.4127 | 299.53 |
| 15 | 216.65 | 12044 | 0.1937 | 295.07 |
| 20 | 216.65 | 5474 | 0.0880 | 295.07 |
| 25 | 221.55 | 2549 | 0.0401 | 298.45 |
| 30 | 226.51 | 1172 | 0.0180 | 301.71 |
| 40 | 250.35 | 277 | 0.0039 | 317.19 |
| 50 | 270.65 | 76 | 0.0010 | 329.80 |

**Equations:**

**Temperature:**
```
T(h) = T_base + L * (h - h_base)
```
Where L is lapse rate [K/m]

**Pressure (gradient layer):**
```
P(h) = P_base * (T_base / T) ^ (g₀ / (R * L))
```

**Pressure (isothermal layer):**
```
P(h) = P_base * exp(-g₀ * (h - h_base) / (R * T))
```

**Density:**
```
ρ = P / (R * T)
```

**Speed of Sound:**
```
a = √(γ * R * T)
```
Where γ = 1.4 for air

---

## Coordinate Systems and Units

### Coordinate Systems

**Geographic Coordinates:**
- **Latitude**: -90° (South Pole) to +90° (North Pole)
- **Longitude**: -180° to +180° or 0° to 360°
- **Note**: MERRA-2 uses 0-360° longitude

**Vertical Coordinates:**

| System | Description | Range | Units |
|--------|-------------|-------|-------|
| Pressure levels | Standard pressure surfaces | 1000 - 0.01 hPa | Pa or hPa |
| Geopotential height | Gravity-adjusted altitude | 0 - 80 km | m |
| Geometric altitude | True altitude MSL | 0 - 1000 km | km |
| Model levels | Hybrid sigma-pressure | 1 - 72 | - |

**Conversion: Geopotential to Geometric Altitude:**
```python
def geopotential_to_geometric(H):
    """
    Convert geopotential height to geometric altitude

    H : geopotential height [m]
    Returns: geometric altitude [m]
    """
    R_earth = 6356766  # Earth radius [m]
    h = R_earth * H / (R_earth - H)
    return h
```

### Units

**MERRA-2 Standard Units:**

| Variable | Units | SI Equivalent |
|----------|-------|---------------|
| Temperature | K | Kelvin |
| Pressure | Pa | Pascal |
| Density | kg/m³ | kg/m³ |
| Wind | m/s | m/s |
| Height | m | meters |
| Humidity (specific) | kg/kg | kg water / kg air |
| Humidity (relative) | % | Percent |
| Time | minutes since epoch | UTC |

**Common Conversions:**

```python
# Temperature
T_C = T_K - 273.15          # Kelvin to Celsius
T_F = T_C * 9/5 + 32        # Celsius to Fahrenheit

# Pressure
P_hPa = P_Pa / 100          # Pa to hPa (millibar)
P_kPa = P_Pa / 1000         # Pa to kPa
P_psi = P_Pa / 6894.76      # Pa to psi
P_atm = P_Pa / 101325       # Pa to atm

# Wind speed
v_kt = v_ms * 1.94384       # m/s to knots
v_mph = v_ms * 2.23694      # m/s to mph
v_kmh = v_ms * 3.6          # m/s to km/h

# Altitude
alt_ft = alt_m * 3.28084    # meters to feet
```

---

## External Resources

### NASA Official Resources

1. **NASA Earthdata Homepage**
   - URL: https://www.earthdata.nasa.gov/
   - Portal for all Earth science data

2. **Earthdata Search**
   - URL: https://search.earthdata.nasa.gov/
   - Visual data discovery tool

3. **Earthdata Login**
   - URL: https://urs.earthdata.nasa.gov/
   - Account creation and management

4. **GES DISC (Goddard Earth Sciences Data and Information Services Center)**
   - URL: https://disc.gsfc.nasa.gov/
   - Primary atmospheric data archive

5. **MERRA-2 Documentation**
   - URL: https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2/
   - Comprehensive reanalysis documentation

6. **AIRS Mission**
   - URL: https://airs.jpl.nasa.gov/
   - Atmospheric sounding mission

### API and Tools

7. **earthaccess Library**
   - GitHub: https://github.com/nsidc/earthaccess
   - Docs: https://earthaccess.readthedocs.io/
   - Install: `pip install earthaccess`

8. **OPeNDAP Documentation**
   - URL: https://www.opendap.org/
   - Data access protocol

9. **CMR API Documentation**
   - URL: https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html
   - Search API reference

10. **pydap Library**
    - GitHub: https://github.com/pydap/pydap
    - OPeNDAP client for Python

### Atmospheric Models

11. **MSISE-00 Model**
    - URL: https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
    - Online calculator and code

12. **pyglow Library**
    - GitHub: https://github.com/timduly4/pyglow
    - Python interface to atmospheric models
    - Install: `pip install pyglow`

13. **US Standard Atmosphere (1976)**
    - NOAA-S/T-76-1562
    - PDF: https://ntrs.nasa.gov/citations/19770009539

14. **ISO 2533:1975**
    - Standard Atmosphere
    - Available from ISO

### Data Processing

15. **NetCDF Documentation**
    - URL: https://www.unidata.ucar.edu/software/netcdf/
    - NetCDF file format reference

16. **xarray Documentation**
    - URL: https://xarray.pydata.org/
    - Python library for multidimensional arrays

17. **CF Conventions**
    - URL: https://cfconventions.org/
    - Climate and Forecast metadata conventions

### Tutorials and Examples

18. **NASA Earthdata Webinars**
    - URL: https://www.earthdata.nasa.gov/learn/webinars
    - Training materials and tutorials

19. **GES DISC Recipes**
    - URL: https://disc.gsfc.nasa.gov/information/recipes
    - Code examples for data access

20. **MERRA-2 Python Examples**
    - URL: https://github.com/GEOS-ESM/MERRA2-Tools
    - Official processing tools

### Aerospace Standards

21. **ICAO Standard Atmosphere**
    - Doc 7488/3
    - International Civil Aviation Organization

22. **AIAA Atmospheric Models**
    - AIAA G-003C-2010
    - Guide to reference and standard atmosphere models

23. **ISO 2533:1975**
    - Standard atmosphere for aerospace use

---

## Quick Reference Command Summary

```python
# --- Authentication ---
import earthaccess
earthaccess.login()  # Interactive or uses .netrc

# --- Search Data ---
results = earthaccess.search_data(
    short_name='M2I3NPASM',  # MERRA-2 3D atmospheric
    temporal=('2023-01-01', '2023-01-31'),
    bounding_box=(-110, 35, -100, 45)
)

# --- Download ---
files = earthaccess.download(results[0:5], './data')

# --- Read NetCDF ---
import xarray as xr
ds = xr.open_dataset(files[0])

# --- Extract Profile ---
profile = ds['T'].sel(
    time='2023-01-01T00:00:00',
    lat=40,
    lon=255,
    method='nearest'
)

# --- Standard Atmosphere ---
from your_module import us_standard_atmosphere_1976
atm = us_standard_atmosphere_1976(altitude_m=10000)
print(f"T: {atm['T']} K, P: {atm['P']} Pa, ρ: {atm['rho']} kg/m³")
```

---

*This reference provides comprehensive information for accessing NASA Earthdata atmospheric datasets for aerospace engineering applications. For the most current information, always consult official NASA Earthdata documentation.*
