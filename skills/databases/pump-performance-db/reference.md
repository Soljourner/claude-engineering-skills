# Pump Performance Database - Reference

Quick reference for accessing manufacturer pump databases, data formats, and custom pump integration.

## Manufacturer Selection Tools

### Major Manufacturers

| Manufacturer | Selection Tool | URL | Data Access |
|--------------|---------------|-----|-------------|
| **Grundfos** | Product Center / WebCAPS | https://product-selection.grundfos.com/ | Web tool, Desktop app, PDF exports |
| **KSB** | EasySelect | https://www.ksb.com/en-us/products-and-solutions/tools-services/easyselect | Web tool, PDF datasheets |
| **Flowserve** | Pump Configurator | https://www.flowserve.com/ | Contact rep, Custom quotes |
| **Goulds/ITT** | Pump Expert | https://www.gouldspumps.com/ | Downloadable software |
| **Xylem** | Various by brand | https://www.xylem.com/ | Multiple tools per brand |
| **Bell & Gossett** | ESP-PLUS / System Syzer | https://www.bellgossett.com/ | Web tool, Desktop app |
| **Armstrong** | Pump Manager | https://www.armstrongfluidtechnology.com/ | Desktop software |
| **Wilo** | Wilo-Select | https://wilo.com/en/Services/Tools/Wilo-Select/ | Web tool, Desktop app |
| **Sulzer** | Blue Box | https://www.sulzer.com/en/products/selection-tools | Web configurator |
| **Pentair** | AquaSuite | https://www.pentair.com/ | Desktop software |
| **Ebara** | E-COMS | https://www.ebarapump.com/ | Desktop software |
| **Ruhrpumpen** | Quick Selector | https://www.ruhrpumpen.com/ | Web tool |
| **Weir** | Various by brand | https://www.global.weir/ | Contact rep |

### Grundfos Resources

**Product Center (Desktop Application)**
- Download: https://product-selection.grundfos.com/download
- Features: Full catalog access, offline use, detailed curves
- Export: PDF, DXF, Excel
- Installation: Windows only, ~2GB

**WebCAPS (Online)**
- URL: https://webcaps.grundfos.com/
- Features: Quick selection, system design
- Export: PDF reports, curve images
- Account: Registration required for full features

**Data Formats**:
- PDF Datasheets: Standard format with embedded curves
- GPA files: Proprietary Grundfos format (requires Product Center)
- Excel exports: Limited availability for some models

### KSB Resources

**EasySelect Online**
- URL: https://www.ksb.com/easyselect
- Features: Pump selection, system design, affinity laws
- Export: PDF datasheets, technical drawings
- Coverage: Most standard pump series

**KSB FluidFuture**
- Desktop application for process pumps
- Advanced hydraulic calculations
- System simulation capabilities

**Data Access**:
- PDF datasheets: Standard format
- DXF/DWG: Technical drawings
- Limited direct curve data export

### Goulds/ITT Resources

**Goulds Pump Expert**
- Download: Available through website
- Features: Comprehensive catalog, curve generation
- Export: PDF, technical specifications
- Updates: Regular software updates

**Online Product Catalog**
- URL: https://www.gouldspumps.com/en-US/Products/
- Browse by type, application
- PDF literature downloads

### Flowserve Resources

**Pump Selection**
- Contact: Local sales representative
- Process: Custom selection for most models
- Documentation: Detailed submittal packages
- Note: Limited self-service selection tools

**Product Literature**
- URL: https://www.flowserve.com/en/products/pumps
- Technical bulletins
- Installation/maintenance manuals
- Performance curves (limited)

## Data Format Specifications

### JSON Format Standard

```json
{
  "metadata": {
    "format_version": "1.0",
    "created_date": "2025-11-07",
    "source": "Manufacturer datasheet",
    "verified": true
  },
  "pump": {
    "manufacturer": "Grundfos",
    "model": "CR 10-5",
    "series": "CR",
    "type": "vertical_multistage_centrifugal",
    "speed_rpm": 2900,
    "frequency_hz": 50,
    "stages": 5,
    "impeller_diameter_mm": 116,
    "impeller_trim": "full"
  },
  "performance": {
    "flow_units": "m3/h",
    "head_units": "m",
    "power_units": "kW",
    "efficiency_units": "percent",
    "npsh_units": "m",
    "curves": [
      {
        "flow": [0, 2, 4, 6, 8, 10, 12],
        "head": [48, 47, 45, 42, 38, 32, 24],
        "efficiency": [0, 45, 62, 70, 72, 68, 55],
        "power": [0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7],
        "npsh_required": [0.5, 0.6, 0.8, 1.2, 1.8, 2.6, 3.6]
      }
    ]
  },
  "specifications": {
    "physical": {
      "weight_kg": 45,
      "length_mm": 580,
      "width_mm": 180,
      "height_mm": 290
    },
    "connections": {
      "suction": "DN 32",
      "discharge": "DN 25",
      "standard": "DIN"
    },
    "operating_limits": {
      "max_pressure_bar": 16,
      "max_temperature_c": 120,
      "min_temperature_c": -10,
      "ph_range": [4, 10]
    },
    "materials": {
      "casing": "Cast Iron EN-GJL-250",
      "impeller": "Stainless Steel 304",
      "shaft": "Stainless Steel 316",
      "seal": "EPDM mechanical seal"
    },
    "motor": {
      "power_kw": 2.2,
      "voltage_v": 400,
      "phases": 3,
      "enclosure": "IE3",
      "insulation_class": "F"
    }
  },
  "application": {
    "suitable_fluids": ["water", "glycol", "light_chemicals"],
    "typical_applications": ["HVAC", "pressure_boosting", "industrial_water"]
  }
}
```

### CSV Format for Curve Data

**Filename Convention**: `[Manufacturer]_[Model]_[Speed]rpm_[ImpellerDia]mm.csv`

Example: `Grundfos_CR10-5_2900rpm_116mm.csv`

```csv
Flow_m3h,Head_m,Efficiency_pct,Power_kW,NPSH_m,Notes
0.0,48.0,0.0,0.80,0.5,Shutoff
2.0,47.0,45.0,1.00,0.6,
4.0,45.0,62.0,1.20,0.8,
6.0,42.0,70.0,1.40,1.2,Near BEP
8.0,38.0,72.0,1.50,1.8,BEP
10.0,32.0,68.0,1.60,2.6,
12.0,24.0,55.0,1.70,3.6,Max flow
```

**Column Definitions**:
- `Flow_m3h`: Volumetric flow rate (m³/h)
- `Head_m`: Total dynamic head (m)
- `Efficiency_pct`: Pump efficiency (%)
- `Power_kW`: Shaft power or motor power (specify in notes)
- `NPSH_m`: Net Positive Suction Head Required (m)
- `Notes`: Optional annotations (BEP, operating limits, etc.)

### SQLite Database Schema

```sql
-- Complete database schema for pump performance data

CREATE TABLE manufacturers (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE,
    website TEXT,
    country TEXT,
    contact_info TEXT,
    notes TEXT
);

CREATE TABLE pump_series (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    manufacturer_id INTEGER,
    series_name TEXT NOT NULL,
    pump_type TEXT,
    description TEXT,
    FOREIGN KEY (manufacturer_id) REFERENCES manufacturers(id)
);

CREATE TABLE pump_models (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    series_id INTEGER,
    model_number TEXT NOT NULL,
    speed_rpm INTEGER,
    frequency_hz INTEGER,
    stages INTEGER,
    impeller_diameter_mm REAL,
    impeller_trim TEXT,
    max_pressure_bar REAL,
    max_temperature_c REAL,
    min_temperature_c REAL,
    weight_kg REAL,
    connection_suction TEXT,
    connection_discharge TEXT,
    date_added TEXT,
    date_modified TEXT,
    data_source TEXT,
    verified BOOLEAN,
    notes TEXT,
    FOREIGN KEY (series_id) REFERENCES pump_series(id)
);

CREATE TABLE performance_curves (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pump_id INTEGER,
    curve_type TEXT,  -- 'head', 'efficiency', 'power', 'npsh'
    flow_m3h REAL,
    value REAL,
    units TEXT,
    FOREIGN KEY (pump_id) REFERENCES pump_models(id)
);

CREATE TABLE materials (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pump_id INTEGER,
    component TEXT,
    material TEXT,
    material_standard TEXT,
    FOREIGN KEY (pump_id) REFERENCES pump_models(id)
);

CREATE TABLE applications (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pump_id INTEGER,
    application_type TEXT,
    suitable BOOLEAN,
    notes TEXT,
    FOREIGN KEY (pump_id) REFERENCES pump_models(id)
);

-- Indexes for performance
CREATE INDEX idx_pump_models_model ON pump_models(model_number);
CREATE INDEX idx_performance_pump_flow ON performance_curves(pump_id, flow_m3h);
CREATE INDEX idx_materials_pump ON materials(pump_id);
```

## Instructions for Adding Custom Pumps

### Method 1: Manual Data Entry

**Step 1: Gather Data**
- Obtain manufacturer datasheet (PDF or hard copy)
- Record pump specifications
- Extract performance curve data

**Step 2: Digitize Curves** (if needed)

Use WebPlotDigitizer:
1. Go to https://automeris.io/WebPlotDigitizer/
2. Load datasheet image (screenshot or scan)
3. Calibrate axes (set known points)
4. Select curve type (e.g., spline)
5. Click points along the curve
6. Export as CSV

**Step 3: Create JSON File**

```python
# Template for manual pump entry
pump_data = {
    "metadata": {
        "format_version": "1.0",
        "created_date": "2025-11-07",
        "source": "Grundfos Datasheet CR_50Hz_2018.pdf",
        "verified": True,
        "notes": "Digitized from PDF using WebPlotDigitizer"
    },
    "pump": {
        "manufacturer": "Grundfos",
        "model": "CR 10-5",
        # ... complete specifications
    },
    "performance": {
        # ... curve data
    }
}

# Save to file
import json
with open('grundfos_cr10-5.json', 'w') as f:
    json.dump(pump_data, f, indent=2)
```

**Step 4: Import to Database**

```python
from query_examples import PumpDatabase, parse_json_pump_data

db = PumpDatabase('my_pump_database.db')
pump = parse_json_pump_data('grundfos_cr10-5.json')
pump_id = db.add_pump(pump)
print(f"Added pump with ID: {pump_id}")
```

### Method 2: Batch Import from CSV

**Prepare CSV File**:

File: `pump_batch_import.csv`

```csv
Manufacturer,Model,Type,Speed,Stages,ImpellerDia,MaxPressure,MaxTemp,CurveFile
Grundfos,CR 10-5,Vertical Multistage,2900,5,116,16,120,curves/grundfos_cr10-5.csv
KSB,Etanorm 080-250,End Suction,1450,1,252,16,140,curves/ksb_etanorm_080-250.csv
```

**Import Script**:

```python
import csv
from query_examples import PumpDatabase, PumpData, PumpSpecifications, PumpPerformance

def batch_import_pumps(csv_file, database):
    """Import multiple pumps from CSV"""
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse curve data from separate file
            curve_data = parse_csv_curve_data(row['CurveFile'])

            # Create pump specifications
            specs = PumpSpecifications(
                model=row['Model'],
                manufacturer=row['Manufacturer'],
                pump_type=row['Type'],
                speed=int(row['Speed']),
                stages=int(row['Stages']),
                impeller_diameter=float(row['ImpellerDia']),
                max_pressure=float(row['MaxPressure']),
                max_temperature=float(row['MaxTemp']),
                connection_size='',
                materials={}
            )

            pump = PumpData(specifications=specs, performance=curve_data)
            pump_id = database.add_pump(pump)
            print(f"Added: {row['Manufacturer']} {row['Model']} (ID: {pump_id})")

# Usage
db = PumpDatabase('my_pump_database.db')
batch_import_pumps('pump_batch_import.csv', db)
```

### Method 3: Direct Database Entry

```python
from query_examples import PumpDatabase

db = PumpDatabase('my_pump_database.db')
cursor = db.conn.cursor()

# Add manufacturer
cursor.execute("""
    INSERT INTO manufacturers (name, website)
    VALUES (?, ?)
""", ('Custom Pumps Inc.', 'https://custompumps.example.com'))
mfg_id = cursor.lastrowid

# Add pump model
cursor.execute("""
    INSERT INTO pump_models
    (manufacturer_id, model_number, pump_type, speed_rpm, stages)
    VALUES (?, ?, ?, ?, ?)
""", (mfg_id, 'CP-100', 'Centrifugal', 1750, 1))
pump_id = cursor.lastrowid

# Add performance curve points
curve_data = [
    (0, 100, 0, 5.0, 1.0),
    (50, 95, 60, 7.5, 1.5),
    (100, 85, 75, 10.0, 2.5),
    (150, 70, 70, 12.5, 4.0),
    (200, 50, 55, 15.0, 6.5),
]

for flow, head, eff, power, npsh in curve_data:
    cursor.execute("""
        INSERT INTO performance_curves
        (pump_id, flow_m3h, head_m, efficiency_pct, power_kw, npsh_m)
        VALUES (?, ?, ?, ?, ?, ?)
    """, (pump_id, flow, head, eff, power, npsh))

db.conn.commit()
print(f"Added custom pump with ID: {pump_id}")
```

### Method 4: Web Scraping (Advanced)

**Important**: Review manufacturer Terms of Service before scraping.

```python
import requests
from bs4 import BeautifulSoup
import time

def scrape_pump_catalog(base_url, max_pages=10):
    """
    Example web scraping template
    NOTE: Customize for specific manufacturer website structure
    """
    pumps = []

    for page in range(1, max_pages + 1):
        url = f"{base_url}?page={page}"

        # Respect rate limiting
        time.sleep(2)

        response = requests.get(url)
        soup = BeautifulSoup(response.content, 'html.parser')

        # Example: Find pump listings
        # (Actual selectors depend on website structure)
        pump_listings = soup.find_all('div', class_='pump-item')

        for listing in pump_listings:
            pump_data = {
                'model': listing.find('span', class_='model').text,
                'type': listing.find('span', class_='type').text,
                'datasheet_url': listing.find('a', class_='datasheet')['href']
            }
            pumps.append(pump_data)

    return pumps

# Use with caution and respect robots.txt
# pumps = scrape_pump_catalog('https://example.com/pumps')
```

## Data Validation Checklist

When adding custom pumps, verify:

- [ ] Flow and head have positive correlation at BEP
- [ ] Efficiency curve is realistic (0-100%, peak < 90% for most pumps)
- [ ] Power increases with flow (generally)
- [ ] NPSH increases with flow
- [ ] Shutoff head > BEP head
- [ ] All units are consistent
- [ ] Data points cover operating range
- [ ] Manufacturer and model clearly identified
- [ ] Source documented for traceability
- [ ] Materials compatible with application

## Common Data Issues and Solutions

### Issue: Missing NPSH Data

**Solution**: Use conservative estimates or manufacturer's standard curves
```python
# Approximate NPSH based on specific speed
import numpy as np

def estimate_npsh(flow, head, speed_rpm):
    """Rough NPSH estimate (use with caution)"""
    specific_speed = speed_rpm * np.sqrt(flow) / head**0.75
    npsh_factor = (specific_speed / 1000) ** 2
    return flow * npsh_factor / 100
```

### Issue: PDF Curves Unreadable

**Solutions**:
1. Contact manufacturer for tabular data
2. Use high-resolution scan and image enhancement
3. Request data from distributor
4. Use manufacturer's selection software

### Issue: Inconsistent Units

**Solution**: Convert to standard units
```python
def convert_units(value, from_unit, to_unit):
    """Unit conversion for pump data"""
    conversions = {
        ('gpm', 'm3/h'): 0.227124,
        ('ft', 'm'): 0.3048,
        ('hp', 'kW'): 0.745699,
        ('psi', 'bar'): 0.0689476,
        ('F', 'C'): lambda x: (x - 32) * 5/9
    }

    if (from_unit, to_unit) in conversions:
        factor = conversions[(from_unit, to_unit)]
        return value * factor if not callable(factor) else factor(value)
    return value
```

## Recommended Directory Structure

```
pump_database/
├── database/
│   └── pumps.db                      # SQLite database
├── json/
│   ├── grundfos/
│   │   ├── cr_series/
│   │   │   ├── cr10-5_50hz.json
│   │   │   └── cr10-5_60hz.json
│   │   └── metadata.json
│   ├── ksb/
│   └── goulds/
├── datasheets/
│   ├── grundfos/
│   │   └── CR_50Hz_2018.pdf
│   └── ksb/
├── curves_csv/
│   ├── grundfos_cr10-5.csv
│   └── ksb_etanorm_080-250.csv
├── scripts/
│   ├── import_pumps.py
│   ├── export_pumps.py
│   └── validate_data.py
└── README.md
```

## Additional Resources

### Tools for Curve Digitization
- **WebPlotDigitizer**: https://automeris.io/WebPlotDigitizer/ (Free, browser-based)
- **Engauge Digitizer**: http://markummitchell.github.io/engauge-digitizer/ (Free, desktop)
- **GraphClick**: https://www.arizona-software.ch/graphclick/ (Mac only)

### Pump Selection Software
- **AFT Fathom/Arrow**: System simulation with pump database
- **PIPE-FLO**: Piping network analysis
- **EPANET**: Water distribution modeling (free)

### Standards and References
- **Hydraulic Institute**: https://www.pumps.org/
  - ANSI/HI 9.6.3: Pump Curves
  - ANSI/HI 1.3: Rotodynamic Pumps
- **ISO 9906**: Rotodynamic pumps - Hydraulic performance acceptance tests

### Community Resources
- **Eng-Tips Forums**: https://www.eng-tips.com/threadminder.cfm?pid=1
- **Pump Engineers Network**: LinkedIn groups
- **r/engineering**: Reddit community

## Legal Disclaimer

When using manufacturer data:
- Respect copyright and intellectual property
- Review Terms of Service for each manufacturer's website
- Use data for legitimate engineering purposes
- Verify critical selections with manufacturer
- Do not redistribute proprietary data without permission
- Consult manufacturer for warranty and liability information

**This skill is for educational and engineering analysis purposes. Always consult manufacturer documentation and representatives for final pump selection and procurement.**
