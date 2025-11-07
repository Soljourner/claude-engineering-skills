---
name: pump-performance-db
description: "Access manufacturer pump curves and specifications from Grundfos, KSB, and other databases"
category: databases
domain: mechanical
complexity: basic
dependencies: []
---

# Pump Performance Database

Access manufacturer pump curves and specifications from Grundfos, KSB, and other major pump manufacturers. This skill provides methods for querying pump performance data, parsing manufacturer databases, and creating custom pump libraries.

## Overview of Pump Databases

Pump performance databases contain critical information for pump selection and analysis:

- **Performance Curves**: Flow rate (Q) vs. Head (H) relationships
- **Efficiency Curves**: Efficiency vs. flow rate
- **NPSH Requirements**: Net Positive Suction Head Required vs. flow rate
- **Power Curves**: Power consumption vs. flow rate
- **Pump Curves**: Multiple impeller diameters and speeds
- **Physical Specifications**: Size, weight, connection types
- **Operating Limits**: Temperature, pressure, fluid compatibility
- **Selection Data**: Model numbers, materials of construction

## Manufacturer Resources

### Grundfos Product Center

**Website**: https://product-selection.grundfos.com/

**Features**:
- Comprehensive product catalog
- WebCAPS online selection tool
- Downloadable product data
- Performance curve exports
- 3D CAD models
- Pump configuration tools

**Data Access**:
- Web-based selection interface
- PDF datasheets with curves
- Excel export capabilities
- API access (limited, requires partnership)
- Desktop software: Grundfos Product Center (offline)

**Typical Data Format**:
- Curve data often in proprietary formats
- PDF datasheets require parsing/digitizing
- Some Excel exports with tabular data

### KSB EasySelect

**Website**: https://www.ksb.com/en-us/products-and-solutions/tools-services/easyselect

**Features**:
- Extensive pump catalog
- Selection based on duty point
- Multiple pump types (centrifugal, positive displacement)
- Performance curve visualization
- Technical documentation

**Data Access**:
- Web-based selection tool
- Downloadable PDF datasheets
- Limited direct data export
- Desktop application available

### Flowserve Pump Selector

**Website**: https://www.flowserve.com/

**Features**:
- Industrial pump selection
- Custom engineered pumps
- Extensive API and process pump catalog
- Material selection guidance

**Data Access**:
- Contact local representative for selection
- PDF technical literature
- Custom quotes for specific applications
- Limited online selection tools

### ITT/Goulds Selection Tools

**Website**: https://www.gouldspumps.com/

**Features**:
- Industrial and commercial pump selection
- i-ALERT condition monitoring integration
- Performance curve data
- Application-specific pump series

**Data Access**:
- Goulds Pump Expert software (downloadable)
- Online product catalog
- PDF literature and curves
- Technical support for custom selections

### Other Manufacturer Resources

- **Xylem/Bell & Gossett**: System Syzer web tool
- **Armstrong**: Pump Manager software
- **Wilo**: Wilo-Select online selection
- **Sulzer**: Blue Box online configurator
- **Pentair**: AquaSuite selection software
- **Ebara**: E-COMS selection software

## Web Scraping Considerations

### Legal and Ethical Considerations

**Before scraping manufacturer websites**:
- Review Terms of Service (ToS)
- Check for robots.txt restrictions
- Respect rate limiting
- Consider requesting API access
- Evaluate if data use complies with licensing

### Technical Approaches

**PDF Parsing**:
```python
# Extract curve data from PDF datasheets
import pdfplumber
import re
import numpy as np

# Digitize performance curves from images
from scipy.interpolate import interp1d
```

**Web Scraping Tools**:
- BeautifulSoup for HTML parsing
- Selenium for JavaScript-heavy sites
- Requests for API-like interfaces
- Tabula for PDF table extraction

**Challenges**:
- Dynamic content loading
- Authentication requirements
- CAPTCHA systems
- Frequent website structure changes
- Data in image format requiring OCR/digitization

## API Access

### Commercial APIs

Most major manufacturers do not provide public APIs. API access typically requires:

1. **Partnership Agreements**: OEM or distributor relationships
2. **Enterprise Licensing**: Large-scale procurement agreements
3. **Software Integration**: Building applications with manufacturer SDKs

### Alternative Data Sources

**Third-Party Aggregators**:
- Pump selection software (e.g., Selecore, PumpBase)
- Industry databases (requires subscription)
- Engineering software integrations (AFT, PIPE-FLO)

**Open Data Projects**:
- Research institutions
- University pump databases
- Community-contributed libraries

## Data Format and Parsing

### Typical Pump Data Structure

```python
{
    "model": "CR 10-5",
    "manufacturer": "Grundfos",
    "type": "Vertical Multistage Centrifugal",
    "impeller_diameter": 116,  # mm
    "speed": 2900,  # rpm
    "stages": 5,
    "performance_curve": {
        "flow": [0, 2, 4, 6, 8, 10, 12],  # m³/h
        "head": [48, 47, 45, 42, 38, 32, 24],  # m
        "efficiency": [0, 45, 62, 70, 72, 68, 55],  # %
        "power": [0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7]  # kW
    },
    "npsh_required": {
        "flow": [0, 2, 4, 6, 8, 10, 12],
        "npsh": [0.5, 0.6, 0.8, 1.2, 1.8, 2.6, 3.6]  # m
    },
    "specifications": {
        "min_flow": 0,  # m³/h
        "max_flow": 12,  # m³/h
        "max_head": 48,  # m
        "max_pressure": 16,  # bar
        "max_temperature": 120,  # °C
        "connection_size": "DN 32/25",
        "materials": {
            "casing": "Cast Iron",
            "impeller": "Stainless Steel",
            "shaft": "Stainless Steel"
        }
    }
}
```

### Parsing Tabular Data

**From Excel/CSV**:
```python
import pandas as pd

# Read manufacturer data export
df = pd.read_excel('pump_data.xlsx', sheet_name='Performance')

# Extract curve data
flow = df['Flow_m3h'].values
head = df['Head_m'].values
efficiency = df['Efficiency_pct'].values
```

**From PDF Tables**:
```python
import tabula

# Extract tables from PDF
tables = tabula.read_pdf('datasheet.pdf', pages='all')
performance_data = tables[0]  # First table
```

### Curve Digitization

**Manual Digitization**:
- Use WebPlotDigitizer (https://automeris.io/WebPlotDigitizer/)
- Export coordinates to CSV
- Import into Python/Excel

**Automated Approaches**:
```python
# Image processing for curve extraction
import cv2
from skimage import morphology
import matplotlib.pyplot as plt

# Read curve image
img = cv2.imread('pump_curve.png', 0)

# Threshold and extract curve points
# (Implementation depends on image quality)
```

## Custom Pump Database Creation

### Database Structure

**SQLite Database Example**:

```sql
CREATE TABLE manufacturers (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    website TEXT
);

CREATE TABLE pump_models (
    id INTEGER PRIMARY KEY,
    manufacturer_id INTEGER,
    model_number TEXT NOT NULL,
    pump_type TEXT,
    speed_rpm INTEGER,
    stages INTEGER,
    impeller_diameter_mm REAL,
    FOREIGN KEY (manufacturer_id) REFERENCES manufacturers(id)
);

CREATE TABLE performance_curves (
    id INTEGER PRIMARY KEY,
    pump_id INTEGER,
    flow_m3h REAL,
    head_m REAL,
    efficiency_pct REAL,
    power_kw REAL,
    npsh_m REAL,
    FOREIGN KEY (pump_id) REFERENCES pump_models(id)
);

CREATE TABLE specifications (
    id INTEGER PRIMARY KEY,
    pump_id INTEGER,
    spec_name TEXT,
    spec_value TEXT,
    spec_unit TEXT,
    FOREIGN KEY (pump_id) REFERENCES pump_models(id)
);
```

### JSON Database Format

```json
{
    "database_version": "1.0",
    "last_updated": "2025-11-07",
    "manufacturers": [
        {
            "id": "grundfos",
            "name": "Grundfos",
            "pumps": [
                {
                    "model": "CR 10-5",
                    "type": "vertical_multistage",
                    "performance": { ... },
                    "specifications": { ... }
                }
            ]
        }
    ]
}
```

### Adding Pump Data Programmatically

```python
import sqlite3
import json

class PumpDatabase:
    def __init__(self, db_path='pump_database.db'):
        self.conn = sqlite3.connect(db_path)
        self.create_tables()

    def add_pump(self, manufacturer, model, pump_type, curve_data, specs):
        """Add a new pump to the database"""
        cursor = self.conn.cursor()

        # Get or create manufacturer
        cursor.execute(
            "INSERT OR IGNORE INTO manufacturers (name) VALUES (?)",
            (manufacturer,)
        )
        cursor.execute(
            "SELECT id FROM manufacturers WHERE name = ?",
            (manufacturer,)
        )
        mfg_id = cursor.fetchone()[0]

        # Insert pump model
        cursor.execute("""
            INSERT INTO pump_models
            (manufacturer_id, model_number, pump_type, speed_rpm, stages)
            VALUES (?, ?, ?, ?, ?)
        """, (mfg_id, model, pump_type,
              specs.get('speed'), specs.get('stages')))

        pump_id = cursor.lastrowid

        # Insert performance curve points
        for i, flow in enumerate(curve_data['flow']):
            cursor.execute("""
                INSERT INTO performance_curves
                (pump_id, flow_m3h, head_m, efficiency_pct, power_kw)
                VALUES (?, ?, ?, ?, ?)
            """, (pump_id, flow, curve_data['head'][i],
                  curve_data['efficiency'][i], curve_data['power'][i]))

        self.conn.commit()
```

## Adding Proprietary Pump Data

### Data Collection Methods

1. **Manufacturer Datasheets**:
   - Download PDF datasheets
   - Extract curve data (digitize if necessary)
   - Record specifications

2. **Factory Test Data**:
   - Witness test reports
   - Performance verification data
   - Actual vs. guaranteed performance

3. **Field Measurements**:
   - Installed pump performance
   - In-service efficiency
   - Degradation tracking

4. **Custom/Modified Pumps**:
   - Trimmed impellers
   - VFD operation at different speeds
   - Special configurations

### Data Organization

```python
# Directory structure for proprietary data
proprietary_pumps/
├── manufacturer_name/
│   ├── model_series/
│   │   ├── model_variant_1.json
│   │   ├── model_variant_2.json
│   │   └── datasheets/
│   │       ├── model_variant_1.pdf
│   │       └── model_variant_2.pdf
│   └── metadata.json
└── custom_pumps/
    ├── project_name/
    │   ├── pump_1.json
    │   └── test_data.csv
    └── metadata.json
```

### Version Control for Pump Data

```python
# Include metadata for tracking
{
    "model": "CR 10-5",
    "data_version": "2.1",
    "last_updated": "2025-11-07",
    "source": "Grundfos Product Center 2025",
    "verified": true,
    "verification_date": "2025-10-15",
    "notes": "Updated efficiency curve based on factory test data"
}
```

## Best Practices

1. **Data Validation**:
   - Verify curve data makes physical sense
   - Check for discontinuities or errors
   - Compare against multiple sources when possible

2. **Units Management**:
   - Store data in consistent units (SI preferred)
   - Document unit conventions
   - Provide conversion utilities

3. **Interpolation**:
   - Use appropriate interpolation for curve queries
   - Avoid extrapolation beyond curve limits
   - Mark interpolated vs. actual data points

4. **Updates and Maintenance**:
   - Track manufacturer updates
   - Version control for database changes
   - Archive deprecated models

5. **Documentation**:
   - Source all data properly
   - Record assumptions and limitations
   - Maintain audit trail for critical selections

## Example Workflows

### Pump Selection Workflow

1. **Define Requirements**: Flow rate, head, fluid properties
2. **Query Database**: Search for suitable pumps
3. **Filter Results**: By type, size, efficiency
4. **Evaluate Performance**: Check operating point on curve
5. **Verify NPSH**: Ensure adequate suction conditions
6. **Review Specifications**: Materials, temperature, pressure ratings
7. **Generate Report**: Document selection basis

### Database Maintenance Workflow

1. **Monitor Manufacturer Updates**: Check for new models/data
2. **Download New Data**: Acquire latest datasheets
3. **Parse and Validate**: Extract and verify curve data
4. **Update Database**: Add new models, update existing
5. **Version Control**: Commit changes with notes
6. **Notify Users**: Document significant changes

## Limitations and Considerations

- **Manufacturer data copyright**: Respect intellectual property
- **Data accuracy**: Verify critical selections with manufacturer
- **Tolerances**: Published curves have manufacturing tolerances
- **Operating conditions**: Curves are for specific fluids (usually water)
- **Application-specific**: Consult manufacturer for special applications
- **Warranty implications**: Unofficial data may void warranties

## Related Skills

- `pump-curves`: Working with performance curves
- `pump-selection`: Systematic pump selection process
- `hydraulic-analysis`: System curve and operating point analysis
- `npsh-analysis`: Net Positive Suction Head calculations
