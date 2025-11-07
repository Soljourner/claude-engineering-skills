---
name: pump-selection-helper
description: "Decision tree for selecting pump type based on flow, head, and fluid properties"
category: helpers
domain: mechanical
complexity: basic
dependencies: []
---

# Pump Selection Helper

A practical decision-tree tool for selecting the appropriate pump type based on operating conditions, flow requirements, head, and fluid properties.

## Pump Type Overview

### Centrifugal Pumps (Dynamic)

#### Radial Flow (Centrifugal)
- **Flow range**: 10 - 100,000 gpm (0.6 - 6,300 L/s)
- **Head range**: 50 - 5,000 ft (15 - 1,500 m)
- **Specific speed (Ns)**: 500 - 4,000 (US units)
- **Applications**: General purpose, high head, moderate to high flow
- **Advantages**: Simple, reliable, low maintenance, handles solids
- **Limitations**: Poor efficiency at low flow, not suitable for high viscosity

#### Mixed Flow
- **Flow range**: 500 - 20,000 gpm (30 - 1,260 L/s)
- **Head range**: 20 - 200 ft (6 - 60 m)
- **Specific speed (Ns)**: 4,000 - 9,000 (US units)
- **Applications**: Irrigation, flood control, water supply
- **Advantages**: Good efficiency, handles moderate flow and head
- **Limitations**: Limited head capability

#### Axial Flow (Propeller)
- **Flow range**: 2,000 - 100,000 gpm (125 - 6,300 L/s)
- **Head range**: 5 - 50 ft (1.5 - 15 m)
- **Specific speed (Ns)**: 9,000 - 15,000 (US units)
- **Applications**: Circulation, cooling water, drainage
- **Advantages**: Very high flow, compact
- **Limitations**: Low head only, sensitive to flow variations

### Positive Displacement Pumps

#### Gear Pumps
- **Flow range**: 1 - 1,500 gpm (0.06 - 95 L/s)
- **Pressure range**: Up to 3,000 psi (200 bar)
- **Viscosity range**: 1 - 1,000,000 cP
- **Applications**: Lubrication oils, fuel transfer, hydraulics
- **Advantages**: Self-priming, handles viscous fluids, constant flow
- **Limitations**: Cannot handle abrasives, pulsating flow

#### Piston/Plunger Pumps
- **Flow range**: 0.1 - 5,000 gpm (0.006 - 315 L/s)
- **Pressure range**: Up to 50,000 psi (3,400 bar)
- **Viscosity range**: 1 - 100,000 cP
- **Applications**: High-pressure cleaning, oil/gas, chemical injection
- **Advantages**: Very high pressure, accurate metering
- **Limitations**: Pulsating flow, high maintenance, expensive

#### Diaphragm Pumps
- **Flow range**: 0.1 - 800 gpm (0.006 - 50 L/s)
- **Pressure range**: Up to 1,000 psi (70 bar)
- **Applications**: Corrosive chemicals, slurries, hazardous fluids
- **Advantages**: Seal-less, handles abrasives and solids
- **Limitations**: Limited pressure, pulsating flow

#### Screw Pumps (Progressive Cavity)
- **Flow range**: 1 - 2,000 gpm (0.06 - 125 L/s)
- **Pressure range**: Up to 1,500 psi (100 bar)
- **Viscosity range**: 1 - 1,000,000 cP
- **Applications**: Viscous fluids, slurries, shear-sensitive fluids
- **Advantages**: Non-pulsating, handles high viscosity, gentle pumping
- **Limitations**: Wear on rotor/stator, limited to moderate pressure

### Specialty Pumps

#### Turbine Pumps (Vertical)
- **Flow range**: 50 - 10,000 gpm (3 - 630 L/s)
- **Head range**: 50 - 1,000 ft (15 - 300 m)
- **Applications**: Deep wells, booster stations, cooling towers
- **Advantages**: Space-efficient, handles high head
- **Limitations**: Complex installation, difficult maintenance

#### Jet Pumps
- **Flow range**: 5 - 100 gpm (0.3 - 6.3 L/s)
- **Head range**: 50 - 300 ft (15 - 90 m)
- **Applications**: Shallow/deep wells, remote locations
- **Advantages**: No moving parts in fluid, simple
- **Limitations**: Low efficiency (25-35%)

#### Airlift Pumps
- **Flow range**: 10 - 5,000 gpm (0.6 - 315 L/s)
- **Applications**: Wells, wastewater, sand/gravel
- **Advantages**: Simple, handles solids and corrosives
- **Limitations**: Very low efficiency, requires air compressor

## Selection Criteria

### 1. Flow Rate Requirements

- **Low flow** (< 10 gpm / 0.6 L/s): Positive displacement preferred
- **Medium flow** (10 - 1,000 gpm / 0.6 - 63 L/s): Centrifugal or PD
- **High flow** (> 1,000 gpm / 63 L/s): Centrifugal (mixed or axial flow)

### 2. Head Requirements

- **Low head** (< 50 ft / 15 m): Axial flow centrifugal or PD
- **Medium head** (50 - 500 ft / 15 - 150 m): Radial centrifugal or PD
- **High head** (> 500 ft / 150 m): High-speed centrifugal or piston pumps

### 3. Specific Speed (Ns)

Specific speed determines the pump impeller type:

**US Units**: Ns = N × √Q / H^0.75
**SI Units**: Ns = N × √Q / H^0.75

Where:
- N = rotational speed (rpm)
- Q = flow rate (gpm or m³/h)
- H = head (ft or m)

**Classification**:
- Ns < 2,000: Radial flow (high head, low flow)
- Ns = 2,000 - 5,000: Francis vane (medium head/flow)
- Ns = 5,000 - 10,000: Mixed flow
- Ns > 10,000: Axial flow (low head, high flow)

### 4. Fluid Viscosity Effects

**Low viscosity** (< 100 cP):
- Centrifugal pumps work well
- No significant correction needed

**Medium viscosity** (100 - 1,000 cP):
- Centrifugal efficiency drops
- Consider positive displacement
- Apply viscosity corrections

**High viscosity** (> 1,000 cP):
- Positive displacement required
- Gear, screw, or piston pumps
- Centrifugal pumps ineffective

### 5. NPSH Requirements

**Net Positive Suction Head** prevents cavitation:

- **NPSHa** (Available) = Atmospheric pressure + Static head - Vapor pressure - Friction losses
- **NPSHr** (Required) = From pump curve (manufacturer data)
- **Requirement**: NPSHa > NPSHr + Safety margin (3-5 ft)

**High NPSHr concerns**:
- Use double suction impeller
- Lower pump speed
- Use inducer or booster pump
- Positive displacement (self-priming)

### 6. Efficiency Considerations

**Best Efficiency Point (BEP)**:
- Centrifugal: Operate within 70-120% of BEP flow
- Peak efficiency: 60-85% for centrifugal
- PD pumps: 70-90% (less flow-dependent)

**Energy cost calculation**:
Annual cost = (BHP × 0.746 × Operating hours × kWh rate) / Efficiency

### 7. Cost Factors

**Initial Cost**:
- Centrifugal: $$ (lowest)
- Gear/Screw: $$$ (moderate)
- Piston/Plunger: $$$$ (highest)

**Operating Cost**:
- Energy consumption
- Maintenance frequency
- Spare parts availability

**Life Cycle Cost** = Initial + Installation + Energy + Maintenance + Downtime

## Decision Tree for Pump Selection

```
START
  |
  ├─ Is fluid viscosity > 1,000 cP?
  |    YES → POSITIVE DISPLACEMENT
  |           ├─ High pressure (> 1,000 psi)? → PISTON/PLUNGER
  |           ├─ Shear-sensitive? → SCREW PUMP
  |           ├─ Abrasive/corrosive? → DIAPHRAGM
  |           └─ General viscous? → GEAR PUMP
  |
  NO ↓
  |
  ├─ Is constant flow required despite pressure changes?
  |    YES → POSITIVE DISPLACEMENT (Gear, Piston, or Screw)
  |
  NO ↓
  |
  ├─ Calculate Specific Speed: Ns = N × √Q / H^0.75
  |
  ├─ Ns < 500? (Very high head, low flow)
  |    YES → TURBINE or HIGH-SPEED CENTRIFUGAL
  |
  ├─ Ns = 500 - 4,000? (High head, moderate flow)
  |    YES → RADIAL CENTRIFUGAL
  |           ├─ Deep well? → VERTICAL TURBINE
  |           └─ Surface? → HORIZONTAL CENTRIFUGAL
  |
  ├─ Ns = 4,000 - 9,000? (Moderate head, high flow)
  |    YES → MIXED FLOW CENTRIFUGAL
  |
  ├─ Ns > 9,000? (Low head, very high flow)
  |    YES → AXIAL FLOW (PROPELLER)
  |
  └─ Special Conditions?
       ├─ Self-priming required? → PD or JET PUMP
       ├─ No electricity available? → ENGINE-DRIVEN
       ├─ Solids > 10% by volume? → DIAPHRAGM or SCREW
       └─ Metering accuracy critical? → PISTON or DIAPHRAGM
```

## Application-Specific Recommendations

### Water Supply
- **Municipal**: Horizontal split-case centrifugal (high reliability)
- **Wells**: Vertical turbine or submersible
- **Booster**: Multistage centrifugal
- **Typical**: Q = 100-5,000 gpm, H = 50-500 ft

### HVAC/Cooling
- **Chilled water**: End-suction or inline centrifugal
- **Condenser water**: Horizontal split-case
- **Typical**: Q = 50-2,000 gpm, H = 30-150 ft

### Chemical Processing
- **Corrosive**: Lined centrifugal or diaphragm
- **Viscous**: Gear or screw pumps
- **Metering**: Diaphragm or piston
- **Typical**: Q = 1-500 gpm, P = 50-500 psi

### Oil & Gas
- **Transfer**: Centrifugal or screw
- **Injection**: High-pressure piston
- **Crude oil**: Screw pumps (viscous)
- **Typical**: Q = 10-1,000 gpm, P = 100-5,000 psi

### Wastewater
- **Raw sewage**: Submersible non-clog centrifugal
- **Sludge**: Progressive cavity (screw)
- **Chemical feed**: Diaphragm metering
- **Typical**: Q = 50-5,000 gpm, H = 20-200 ft

### Agriculture/Irrigation
- **Surface water**: Horizontal centrifugal
- **Wells**: Vertical turbine
- **Drip irrigation**: Centrifugal with filtration
- **Typical**: Q = 100-5,000 gpm, H = 50-300 ft

### Mining/Slurry
- **Heavy slurry**: Horizontal slurry pump (rubber-lined)
- **Abrasive**: Hard-metal or ceramic-lined
- **Dewatering**: Submersible or horizontal centrifugal
- **Typical**: Q = 100-10,000 gpm, H = 50-500 ft

### Food & Beverage
- **Sanitary**: Centrifugal (polished, 3A certified)
- **Viscous products**: Lobe or screw pumps
- **CIP/Cleaning**: Centrifugal
- **Typical**: Q = 10-500 gpm, P = 50-150 psi

## Usage Guide

### Using selector.py

Run the interactive selector:
```bash
python selector.py
```

Or use programmatically:
```python
from selector import select_pump, calculate_specific_speed

# Example 1: Water supply
result = select_pump(
    flow_rate=500,      # gpm
    head=200,           # ft
    viscosity=1,        # cP
    fluid_type="water",
    speed=1750          # rpm
)
print(result)

# Example 2: High viscosity
result = select_pump(
    flow_rate=50,
    head=100,
    viscosity=5000,
    fluid_type="oil",
    temp=100            # °F
)
print(result)
```

### Quick Selection Guidelines

1. **Start with flow and head** - These are primary factors
2. **Check fluid properties** - Viscosity, corrosiveness, abrasiveness
3. **Calculate specific speed** - Determines centrifugal type
4. **Verify NPSH** - Ensure adequate suction conditions
5. **Consider operating range** - Pump should operate near BEP
6. **Evaluate life cycle cost** - Not just initial cost
7. **Check maintenance access** - Space for service
8. **Review manufacturer curves** - Verify actual performance

### Common Mistakes to Avoid

- ❌ Oversizing pumps (reduces efficiency, increases cost)
- ❌ Ignoring viscosity effects on centrifugal pumps
- ❌ Insufficient NPSH margin (causes cavitation)
- ❌ Operating far from BEP (premature wear)
- ❌ Selecting based on initial cost only
- ❌ Not considering future expansion needs
- ❌ Ignoring system curve changes

## References

See `reference.md` for:
- Detailed pump selection charts
- Specific speed ranges from industry standards
- Manufacturer selection guides
- Performance curve examples
- NPSH calculation methods
