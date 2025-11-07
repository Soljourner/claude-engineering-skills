"""
Pump Performance Database - Query Examples

This module demonstrates querying and working with pump performance databases,
including searching by requirements, filtering by type, and parsing manufacturer data.
"""

import json
import sqlite3
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import numpy as np
from scipy.interpolate import interp1d


# ============================================================================
# Example Pump Data Structures
# ============================================================================

@dataclass
class PumpPerformance:
    """Structure for pump performance curve data"""
    flow: List[float]  # m³/h
    head: List[float]  # m
    efficiency: List[float]  # %
    power: List[float]  # kW
    npsh_required: List[float]  # m


@dataclass
class PumpSpecifications:
    """Structure for pump specifications"""
    model: str
    manufacturer: str
    pump_type: str
    speed: int  # rpm
    stages: int
    impeller_diameter: float  # mm
    max_pressure: float  # bar
    max_temperature: float  # °C
    connection_size: str
    materials: Dict[str, str]


@dataclass
class PumpData:
    """Complete pump data including performance and specifications"""
    specifications: PumpSpecifications
    performance: PumpPerformance

    def get_performance_at_flow(self, flow_rate: float) -> Dict[str, float]:
        """
        Interpolate pump performance at a specific flow rate

        Args:
            flow_rate: Flow rate in m³/h

        Returns:
            Dictionary with head, efficiency, power, and NPSH at the flow rate
        """
        if flow_rate < min(self.performance.flow) or flow_rate > max(self.performance.flow):
            raise ValueError(f"Flow rate {flow_rate} outside pump curve range")

        # Create interpolation functions
        head_interp = interp1d(self.performance.flow, self.performance.head, kind='cubic')
        eff_interp = interp1d(self.performance.flow, self.performance.efficiency, kind='cubic')
        power_interp = interp1d(self.performance.flow, self.performance.power, kind='cubic')
        npsh_interp = interp1d(self.performance.flow, self.performance.npsh_required, kind='cubic')

        return {
            'flow': flow_rate,
            'head': float(head_interp(flow_rate)),
            'efficiency': float(eff_interp(flow_rate)),
            'power': float(power_interp(flow_rate)),
            'npsh_required': float(npsh_interp(flow_rate))
        }

    def find_bep(self) -> Tuple[float, float, float]:
        """
        Find Best Efficiency Point (BEP)

        Returns:
            Tuple of (flow, head, efficiency) at BEP
        """
        max_eff_idx = np.argmax(self.performance.efficiency)
        return (
            self.performance.flow[max_eff_idx],
            self.performance.head[max_eff_idx],
            self.performance.efficiency[max_eff_idx]
        )


# Example pump data instances
EXAMPLE_PUMP_1 = PumpData(
    specifications=PumpSpecifications(
        model="CR 10-5",
        manufacturer="Grundfos",
        pump_type="Vertical Multistage Centrifugal",
        speed=2900,
        stages=5,
        impeller_diameter=116,
        max_pressure=16,
        max_temperature=120,
        connection_size="DN 32/25",
        materials={
            "casing": "Cast Iron",
            "impeller": "Stainless Steel",
            "shaft": "Stainless Steel"
        }
    ),
    performance=PumpPerformance(
        flow=[0, 2, 4, 6, 8, 10, 12],
        head=[48, 47, 45, 42, 38, 32, 24],
        efficiency=[0, 45, 62, 70, 72, 68, 55],
        power=[0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7],
        npsh_required=[0.5, 0.6, 0.8, 1.2, 1.8, 2.6, 3.6]
    )
)

EXAMPLE_PUMP_2 = PumpData(
    specifications=PumpSpecifications(
        model="Etanorm 080-250",
        manufacturer="KSB",
        pump_type="End Suction Centrifugal",
        speed=1450,
        stages=1,
        impeller_diameter=252,
        max_pressure=16,
        max_temperature=140,
        connection_size="DN 80/65",
        materials={
            "casing": "Cast Iron",
            "impeller": "Bronze",
            "shaft": "Stainless Steel"
        }
    ),
    performance=PumpPerformance(
        flow=[0, 50, 100, 150, 200, 250, 300],
        head=[35, 34, 32, 28, 23, 16, 8],
        efficiency=[0, 48, 68, 78, 80, 75, 60],
        power=[5.0, 6.5, 8.5, 11.0, 13.5, 15.5, 17.0],
        npsh_required=[1.5, 1.8, 2.3, 3.2, 4.5, 6.5, 9.0]
    )
)

EXAMPLE_PUMP_3 = PumpData(
    specifications=PumpSpecifications(
        model="3196 MT 3x4-10",
        manufacturer="Goulds",
        pump_type="Process Centrifugal (ANSI)",
        speed=3560,
        stages=1,
        impeller_diameter=254,
        max_pressure=20,
        max_temperature=200,
        connection_size="3\" x 4\"",
        materials={
            "casing": "Ductile Iron",
            "impeller": "316 Stainless Steel",
            "shaft": "316 Stainless Steel"
        }
    ),
    performance=PumpPerformance(
        flow=[0, 100, 200, 300, 400, 500, 600],
        head=[92, 90, 86, 80, 72, 60, 45],
        efficiency=[0, 42, 64, 76, 82, 80, 70],
        power=[15, 20, 28, 38, 48, 58, 68],
        npsh_required=[2.0, 2.5, 3.5, 5.0, 7.5, 11.0, 15.5]
    )
)


# ============================================================================
# Database Class for Pump Data Management
# ============================================================================

class PumpDatabase:
    """SQLite-based pump performance database"""

    def __init__(self, db_path: str = 'pump_database.db'):
        """Initialize database connection and create tables if needed"""
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Enable column access by name
        self.create_tables()

    def create_tables(self):
        """Create database tables if they don't exist"""
        cursor = self.conn.cursor()

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS manufacturers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL UNIQUE,
                website TEXT
            )
        """)

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS pump_models (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                manufacturer_id INTEGER,
                model_number TEXT NOT NULL,
                pump_type TEXT,
                speed_rpm INTEGER,
                stages INTEGER,
                impeller_diameter_mm REAL,
                max_pressure_bar REAL,
                max_temperature_c REAL,
                connection_size TEXT,
                FOREIGN KEY (manufacturer_id) REFERENCES manufacturers(id)
            )
        """)

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS performance_curves (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pump_id INTEGER,
                flow_m3h REAL,
                head_m REAL,
                efficiency_pct REAL,
                power_kw REAL,
                npsh_m REAL,
                FOREIGN KEY (pump_id) REFERENCES pump_models(id)
            )
        """)

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS materials (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pump_id INTEGER,
                component TEXT,
                material TEXT,
                FOREIGN KEY (pump_id) REFERENCES pump_models(id)
            )
        """)

        self.conn.commit()

    def add_pump(self, pump_data: PumpData) -> int:
        """
        Add a pump to the database

        Args:
            pump_data: PumpData object containing specifications and performance

        Returns:
            pump_id: Database ID of the added pump
        """
        cursor = self.conn.cursor()
        specs = pump_data.specifications
        perf = pump_data.performance

        # Get or create manufacturer
        cursor.execute(
            "INSERT OR IGNORE INTO manufacturers (name) VALUES (?)",
            (specs.manufacturer,)
        )
        cursor.execute(
            "SELECT id FROM manufacturers WHERE name = ?",
            (specs.manufacturer,)
        )
        mfg_id = cursor.fetchone()[0]

        # Insert pump model
        cursor.execute("""
            INSERT INTO pump_models
            (manufacturer_id, model_number, pump_type, speed_rpm, stages,
             impeller_diameter_mm, max_pressure_bar, max_temperature_c, connection_size)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (mfg_id, specs.model, specs.pump_type, specs.speed, specs.stages,
              specs.impeller_diameter, specs.max_pressure, specs.max_temperature,
              specs.connection_size))

        pump_id = cursor.lastrowid

        # Insert performance curve points
        for i in range(len(perf.flow)):
            cursor.execute("""
                INSERT INTO performance_curves
                (pump_id, flow_m3h, head_m, efficiency_pct, power_kw, npsh_m)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (pump_id, perf.flow[i], perf.head[i], perf.efficiency[i],
                  perf.power[i], perf.npsh_required[i]))

        # Insert materials
        for component, material in specs.materials.items():
            cursor.execute("""
                INSERT INTO materials (pump_id, component, material)
                VALUES (?, ?, ?)
            """, (pump_id, component, material))

        self.conn.commit()
        return pump_id

    def get_pump(self, pump_id: int) -> Optional[PumpData]:
        """
        Retrieve a pump from the database by ID

        Args:
            pump_id: Database ID of the pump

        Returns:
            PumpData object or None if not found
        """
        cursor = self.conn.cursor()

        # Get pump specifications
        cursor.execute("""
            SELECT pm.*, m.name as manufacturer_name
            FROM pump_models pm
            JOIN manufacturers m ON pm.manufacturer_id = m.id
            WHERE pm.id = ?
        """, (pump_id,))

        pump_row = cursor.fetchone()
        if not pump_row:
            return None

        # Get performance curve
        cursor.execute("""
            SELECT flow_m3h, head_m, efficiency_pct, power_kw, npsh_m
            FROM performance_curves
            WHERE pump_id = ?
            ORDER BY flow_m3h
        """, (pump_id,))

        curve_data = cursor.fetchall()

        # Get materials
        cursor.execute("""
            SELECT component, material
            FROM materials
            WHERE pump_id = ?
        """, (pump_id,))

        materials = {row['component']: row['material'] for row in cursor.fetchall()}

        # Build PumpData object
        specs = PumpSpecifications(
            model=pump_row['model_number'],
            manufacturer=pump_row['manufacturer_name'],
            pump_type=pump_row['pump_type'],
            speed=pump_row['speed_rpm'],
            stages=pump_row['stages'],
            impeller_diameter=pump_row['impeller_diameter_mm'],
            max_pressure=pump_row['max_pressure_bar'],
            max_temperature=pump_row['max_temperature_c'],
            connection_size=pump_row['connection_size'],
            materials=materials
        )

        performance = PumpPerformance(
            flow=[row['flow_m3h'] for row in curve_data],
            head=[row['head_m'] for row in curve_data],
            efficiency=[row['efficiency_pct'] for row in curve_data],
            power=[row['power_kw'] for row in curve_data],
            npsh_required=[row['npsh_m'] for row in curve_data]
        )

        return PumpData(specifications=specs, performance=performance)


# ============================================================================
# Search by Requirements (Q, H)
# ============================================================================

def search_pumps_by_duty_point(
    database: PumpDatabase,
    required_flow: float,
    required_head: float,
    head_tolerance: float = 0.1,
    npsh_available: Optional[float] = None
) -> List[Tuple[int, PumpData, Dict[str, float]]]:
    """
    Search for pumps that can meet specified flow and head requirements

    Args:
        database: PumpDatabase instance
        required_flow: Required flow rate in m³/h
        required_head: Required head in m
        head_tolerance: Acceptable head deviation (default 10%)
        npsh_available: Available NPSH in m (optional)

    Returns:
        List of tuples: (pump_id, PumpData, operating_point)
    """
    cursor = database.conn.cursor()

    # Get all pumps
    cursor.execute("SELECT id FROM pump_models")
    pump_ids = [row[0] for row in cursor.fetchall()]

    suitable_pumps = []

    for pump_id in pump_ids:
        pump = database.get_pump(pump_id)
        if not pump:
            continue

        # Check if flow is within pump curve range
        flow_min = min(pump.performance.flow)
        flow_max = max(pump.performance.flow)

        if required_flow < flow_min or required_flow > flow_max:
            continue

        try:
            # Get performance at required flow
            op_point = pump.get_performance_at_flow(required_flow)

            # Check if head meets requirements (with tolerance)
            head_min = required_head * (1 - head_tolerance)
            head_max = required_head * (1 + head_tolerance)

            if head_min <= op_point['head'] <= head_max:
                # Check NPSH if specified
                if npsh_available is not None:
                    if op_point['npsh_required'] > npsh_available:
                        continue  # Insufficient NPSH

                suitable_pumps.append((pump_id, pump, op_point))

        except ValueError:
            continue

    # Sort by efficiency (highest first)
    suitable_pumps.sort(key=lambda x: x[2]['efficiency'], reverse=True)

    return suitable_pumps


def search_pumps_simple(
    required_flow: float,
    required_head: float,
    pump_list: List[PumpData]
) -> List[Tuple[PumpData, Dict[str, float]]]:
    """
    Simple search through a list of pump objects (no database)

    Args:
        required_flow: Required flow rate in m³/h
        required_head: Required head in m
        required_head_tolerance: Acceptable deviation (default 10%)
        pump_list: List of PumpData objects to search

    Returns:
        List of tuples: (PumpData, operating_point)
    """
    suitable_pumps = []

    for pump in pump_list:
        flow_min = min(pump.performance.flow)
        flow_max = max(pump.performance.flow)

        if flow_min <= required_flow <= flow_max:
            try:
                op_point = pump.get_performance_at_flow(required_flow)

                # Check if pump can provide required head
                if op_point['head'] >= required_head:
                    suitable_pumps.append((pump, op_point))

            except ValueError:
                continue

    # Sort by efficiency
    suitable_pumps.sort(key=lambda x: x[1]['efficiency'], reverse=True)

    return suitable_pumps


# ============================================================================
# Filter by Pump Type
# ============================================================================

def filter_by_pump_type(
    database: PumpDatabase,
    pump_type: str
) -> List[PumpData]:
    """
    Filter pumps by type

    Args:
        database: PumpDatabase instance
        pump_type: Pump type to filter (e.g., "Centrifugal", "Multistage")

    Returns:
        List of PumpData objects matching the type
    """
    cursor = database.conn.cursor()

    cursor.execute("""
        SELECT id FROM pump_models
        WHERE pump_type LIKE ?
    """, (f'%{pump_type}%',))

    pump_ids = [row[0] for row in cursor.fetchall()]

    return [database.get_pump(pump_id) for pump_id in pump_ids
            if database.get_pump(pump_id) is not None]


def filter_pumps_advanced(
    database: PumpDatabase,
    pump_type: Optional[str] = None,
    manufacturer: Optional[str] = None,
    min_flow: Optional[float] = None,
    max_flow: Optional[float] = None,
    min_head: Optional[float] = None,
    max_head: Optional[float] = None,
    max_temperature: Optional[float] = None,
    material: Optional[str] = None
) -> List[PumpData]:
    """
    Advanced filtering with multiple criteria

    Args:
        database: PumpDatabase instance
        pump_type: Pump type substring to match (optional)
        manufacturer: Manufacturer name (optional)
        min_flow: Minimum flow capacity in m³/h (optional)
        max_flow: Maximum flow capacity in m³/h (optional)
        min_head: Minimum head capacity in m (optional)
        max_head: Maximum head capacity in m (optional)
        max_temperature: Required temperature rating in °C (optional)
        material: Material requirement (optional)

    Returns:
        List of PumpData objects matching all criteria
    """
    cursor = database.conn.cursor()

    # Build dynamic query
    query = """
        SELECT DISTINCT pm.id
        FROM pump_models pm
        JOIN manufacturers m ON pm.manufacturer_id = m.id
        LEFT JOIN performance_curves pc ON pm.id = pc.pump_id
        LEFT JOIN materials mat ON pm.id = mat.pump_id
        WHERE 1=1
    """
    params = []

    if pump_type:
        query += " AND pm.pump_type LIKE ?"
        params.append(f'%{pump_type}%')

    if manufacturer:
        query += " AND m.name LIKE ?"
        params.append(f'%{manufacturer}%')

    if max_temperature:
        query += " AND pm.max_temperature_c >= ?"
        params.append(max_temperature)

    if material:
        query += " AND mat.material LIKE ?"
        params.append(f'%{material}%')

    cursor.execute(query, params)
    pump_ids = [row[0] for row in cursor.fetchall()]

    # Filter by flow and head ranges (requires checking curve data)
    filtered_pumps = []

    for pump_id in pump_ids:
        pump = database.get_pump(pump_id)
        if not pump:
            continue

        pump_flow_min = min(pump.performance.flow)
        pump_flow_max = max(pump.performance.flow)
        pump_head_min = min(pump.performance.head)
        pump_head_max = max(pump.performance.head)

        if min_flow and pump_flow_max < min_flow:
            continue
        if max_flow and pump_flow_min > max_flow:
            continue
        if min_head and pump_head_max < min_head:
            continue
        if max_head and pump_head_min > max_head:
            continue

        filtered_pumps.append(pump)

    return filtered_pumps


# ============================================================================
# Parse Manufacturer Data
# ============================================================================

def parse_json_pump_data(json_file: str) -> PumpData:
    """
    Parse pump data from JSON file

    Args:
        json_file: Path to JSON file containing pump data

    Returns:
        PumpData object
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    specs = PumpSpecifications(
        model=data['model'],
        manufacturer=data['manufacturer'],
        pump_type=data['type'],
        speed=data['speed'],
        stages=data.get('stages', 1),
        impeller_diameter=data.get('impeller_diameter', 0),
        max_pressure=data['specifications']['max_pressure'],
        max_temperature=data['specifications']['max_temperature'],
        connection_size=data['specifications'].get('connection_size', ''),
        materials=data['specifications'].get('materials', {})
    )

    perf_data = data['performance_curve']
    npsh_data = data.get('npsh_required', {})

    performance = PumpPerformance(
        flow=perf_data['flow'],
        head=perf_data['head'],
        efficiency=perf_data['efficiency'],
        power=perf_data['power'],
        npsh_required=npsh_data.get('npsh', [0] * len(perf_data['flow']))
    )

    return PumpData(specifications=specs, performance=performance)


def parse_csv_curve_data(csv_file: str) -> PumpPerformance:
    """
    Parse performance curve data from CSV file

    Expected CSV format:
    Flow(m3/h),Head(m),Efficiency(%),Power(kW),NPSH(m)
    0,48,0,0.8,0.5
    2,47,45,1.0,0.6
    ...

    Args:
        csv_file: Path to CSV file

    Returns:
        PumpPerformance object
    """
    import csv

    flow, head, efficiency, power, npsh = [], [], [], [], []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            flow.append(float(row['Flow(m3/h)']))
            head.append(float(row['Head(m)']))
            efficiency.append(float(row['Efficiency(%)']))
            power.append(float(row['Power(kW)']))
            npsh.append(float(row.get('NPSH(m)', 0)))

    return PumpPerformance(
        flow=flow,
        head=head,
        efficiency=efficiency,
        power=power,
        npsh_required=npsh
    )


def export_pump_to_json(pump: PumpData, output_file: str):
    """
    Export pump data to JSON file

    Args:
        pump: PumpData object
        output_file: Path to output JSON file
    """
    specs = pump.specifications
    perf = pump.performance

    data = {
        "model": specs.model,
        "manufacturer": specs.manufacturer,
        "type": specs.pump_type,
        "speed": specs.speed,
        "stages": specs.stages,
        "impeller_diameter": specs.impeller_diameter,
        "performance_curve": {
            "flow": perf.flow,
            "head": perf.head,
            "efficiency": perf.efficiency,
            "power": perf.power
        },
        "npsh_required": {
            "flow": perf.flow,
            "npsh": perf.npsh_required
        },
        "specifications": {
            "max_pressure": specs.max_pressure,
            "max_temperature": specs.max_temperature,
            "connection_size": specs.connection_size,
            "materials": specs.materials
        }
    }

    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)


# ============================================================================
# Example Usage
# ============================================================================

def main():
    """Example usage of pump database functions"""

    # Create database and add example pumps
    print("Creating pump database...")
    db = PumpDatabase('example_pump_database.db')

    print("\nAdding example pumps to database...")
    pump_id_1 = db.add_pump(EXAMPLE_PUMP_1)
    pump_id_2 = db.add_pump(EXAMPLE_PUMP_2)
    pump_id_3 = db.add_pump(EXAMPLE_PUMP_3)

    print(f"Added pumps with IDs: {pump_id_1}, {pump_id_2}, {pump_id_3}")

    # Search for pumps by duty point
    print("\n" + "="*60)
    print("Example 1: Search by Duty Point")
    print("="*60)
    required_flow = 8.0  # m³/h
    required_head = 38.0  # m
    npsh_available = 3.0  # m

    print(f"\nSearching for pumps that can deliver:")
    print(f"  Flow: {required_flow} m³/h")
    print(f"  Head: {required_head} m")
    print(f"  NPSH Available: {npsh_available} m")

    results = search_pumps_by_duty_point(
        db, required_flow, required_head,
        head_tolerance=0.1, npsh_available=npsh_available
    )

    print(f"\nFound {len(results)} suitable pump(s):")
    for pump_id, pump, op_point in results:
        specs = pump.specifications
        print(f"\n  {specs.manufacturer} {specs.model}")
        print(f"    Type: {specs.pump_type}")
        print(f"    Operating Point:")
        print(f"      Flow: {op_point['flow']:.1f} m³/h")
        print(f"      Head: {op_point['head']:.1f} m")
        print(f"      Efficiency: {op_point['efficiency']:.1f}%")
        print(f"      Power: {op_point['power']:.1f} kW")
        print(f"      NPSH Required: {op_point['npsh_required']:.2f} m")

    # Filter by pump type
    print("\n" + "="*60)
    print("Example 2: Filter by Pump Type")
    print("="*60)
    pump_type = "Centrifugal"
    print(f"\nFiltering for pump type: {pump_type}")

    centrifugal_pumps = filter_by_pump_type(db, pump_type)
    print(f"\nFound {len(centrifugal_pumps)} centrifugal pump(s):")
    for pump in centrifugal_pumps:
        specs = pump.specifications
        bep_flow, bep_head, bep_eff = pump.find_bep()
        print(f"\n  {specs.manufacturer} {specs.model}")
        print(f"    Type: {specs.pump_type}")
        print(f"    BEP: {bep_flow:.0f} m³/h @ {bep_head:.0f} m ({bep_eff:.0f}%)")

    # Advanced filtering
    print("\n" + "="*60)
    print("Example 3: Advanced Filtering")
    print("="*60)

    filtered_pumps = filter_pumps_advanced(
        db,
        pump_type="Centrifugal",
        min_flow=50,
        max_flow=500,
        min_head=20,
        material="Stainless Steel"
    )

    print(f"\nFiltered pumps (Centrifugal, 50-500 m³/h, >20m head, SS impeller):")
    print(f"Found {len(filtered_pumps)} pump(s):")
    for pump in filtered_pumps:
        specs = pump.specifications
        print(f"\n  {specs.manufacturer} {specs.model}")
        print(f"    Flow range: {min(pump.performance.flow)}-{max(pump.performance.flow)} m³/h")
        print(f"    Head range: {min(pump.performance.head)}-{max(pump.performance.head)} m")
        print(f"    Impeller: {specs.materials.get('impeller', 'N/A')}")

    # Export example
    print("\n" + "="*60)
    print("Example 4: Export Pump Data")
    print("="*60)

    export_file = "pump_export_example.json"
    export_pump_to_json(EXAMPLE_PUMP_1, export_file)
    print(f"\nExported pump data to: {export_file}")

    print("\n" + "="*60)
    print("Examples completed successfully!")
    print("="*60)


if __name__ == "__main__":
    main()
