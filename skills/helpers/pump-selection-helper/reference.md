# Pump Selection Reference Guide

Comprehensive reference data for pump selection including charts, standards, and manufacturer guidelines.

## Table of Contents

1. [Pump Selection Charts](#pump-selection-charts)
2. [Specific Speed Ranges](#specific-speed-ranges)
3. [Viscosity Correction Factors](#viscosity-correction-factors)
4. [NPSH Calculations](#npsh-calculations)
5. [Application Guidelines](#application-guidelines)
6. [Manufacturer Resources](#manufacturer-resources)
7. [Industry Standards](#industry-standards)

---

## Pump Selection Charts

### Flow vs Head Operating Ranges

```
HEAD (ft)
10,000 |                    * Piston/Plunger
       |                   *
 5,000 |                  *
       |                 *  Multi-stage
 2,000 |                *   Centrifugal
       |               *
 1,000 |              *
       |         ****
   500 |     ****      Single-stage
       | ****          Centrifugal
   200 |**
       |**
   100 |***
       | ****
    50 |  *****  Mixed Flow
       |     *******
    20 |         ********
       |              *****  Axial Flow
    10 |                ******
       +-----|-----|-----|-----|-----|-----|-----|-----
            1    10   100  1K  10K 100K  1M  FLOW (gpm)

Legend:
* = Operating range boundary
```

### Pump Type Selection by Specific Speed

```
Ns Range    | Pump Type           | H/Q Ratio | Best Application
------------|---------------------|-----------|------------------
< 500       | Turbine/High-speed  | Very High | Deep wells, high head
500-2,000   | Radial centrifugal  | High      | General high head
2,000-4,000 | Francis/Radial      | Medium    | Balanced head/flow
4,000-9,000 | Mixed flow          | Low       | Irrigation, drainage
9,000+      | Axial flow          | Very Low  | Circulation, flood

Ns = N × √Q / H^0.75
where: N = rpm, Q = gpm, H = ft
```

---

## Specific Speed Ranges

### Detailed Specific Speed Classifications

#### US Units (rpm, gpm, ft)

| Ns Range | Impeller Type | Characteristics | Typical Applications |
|----------|---------------|-----------------|---------------------|
| 500-1,000 | Narrow radial | Very high head, narrow vanes | Boiler feed, high-rise buildings |
| 1,000-1,500 | Radial | High head, efficient | Chemical plants, refineries |
| 1,500-2,500 | Radial | General purpose | HVAC, water supply |
| 2,500-4,000 | Francis vane | Balanced design | Municipal water, irrigation |
| 4,000-7,000 | Mixed flow | Higher flow capacity | Large water systems |
| 7,000-9,000 | Mixed flow | High flow, moderate head | Flood control, cooling |
| 9,000-15,000 | Axial flow | Very high flow, low head | Circulation, drainage |

#### SI Units (rpm, m³/h, m)

| Ns Range | Impeller Type | Notes |
|----------|---------------|-------|
| 10-30 | Narrow radial | Multistage applications |
| 30-50 | Radial | Standard industrial |
| 50-80 | Radial/Francis | Wide application range |
| 80-160 | Mixed flow | Large volume pumps |
| 160-300 | Axial flow | Very high capacity |

### Conversion Between Units

**US to SI**: Ns(SI) = Ns(US) × 0.0196

**SI to US**: Ns(US) = Ns(SI) × 51.0

---

## Viscosity Correction Factors

### Centrifugal Pump Performance Correction

For fluids with viscosity > 1 cP, centrifugal pump performance must be corrected:

#### Correction Factor Table

| Viscosity (SSU) | Flow Factor (CQ) | Head Factor (CH) | Efficiency Factor (Cη) |
|-----------------|------------------|------------------|------------------------|
| 100 (1 cP) | 1.00 | 1.00 | 1.00 |
| 500 (10 cP) | 0.98 | 0.98 | 0.90 |
| 1,000 (60 cP) | 0.95 | 0.96 | 0.75 |
| 2,000 (200 cP) | 0.90 | 0.92 | 0.55 |
| 5,000 (1,000 cP) | 0.80 | 0.85 | 0.35 |
| 10,000 (2,200 cP) | 0.68 | 0.75 | 0.20 |

**Corrected Values**:
- Q_actual = Q_water × CQ
- H_actual = H_water × CH
- η_actual = η_water × Cη

#### Hydraulic Institute Method

For more accurate corrections, use HI 9.6.7 Correction Charts:

1. Calculate capacity coefficient: CQ = 0.95 × (ν/100)^-0.1
2. Calculate head coefficient: CH = 0.97 × (ν/100)^-0.05
3. Calculate efficiency: Cη = 1 - 0.015 × log(ν)

Where ν = kinematic viscosity in SSU

### Viscosity Conversion

**Centipoise (cP) to SSU**:
- SSU = 4.632 × cP (for cP < 70)
- SSU = 4.664 × cP - 0.164 (for cP > 70)

**SSU to Centistokes (cSt)**:
- cSt = 0.226 × SSU - 195/SSU (for 32 < SSU < 100)
- cSt = 0.220 × SSU - 135/SSU (for SSU > 100)

---

## NPSH Calculations

### Net Positive Suction Head Available (NPSHa)

**Formula**:
```
NPSHa = (Pa / γ) + Zs - (Pvap / γ) - hL - hf

Where:
Pa   = Absolute pressure at liquid surface (psi or kPa)
Pvap = Vapor pressure of liquid at pumping temp (psi or kPa)
γ    = Specific weight of liquid (lb/ft³ or kg/m³)
Zs   = Static suction head (+) or lift (-) in ft or m
hL   = Suction line losses (ft or m)
hf   = Fitting losses (ft or m)
```

**Simplified (for water at sea level)**:
```
NPSHa = 34 ft + Static head - Vapor pressure head - Friction losses
```

### Vapor Pressure of Water

| Temperature (°F) | Vapor Pressure (psia) | Head (ft) |
|------------------|----------------------|-----------|
| 60 | 0.26 | 0.6 |
| 80 | 0.51 | 1.2 |
| 100 | 0.95 | 2.2 |
| 140 | 2.89 | 6.7 |
| 180 | 7.51 | 17.3 |
| 212 | 14.70 | 34.0 |

### NPSH Margin Requirements

**Minimum safety margin**: NPSHa > NPSHr + 3 to 5 ft

**Application-specific margins**:
- Clean water, steady service: +3 ft
- Hydrocarbons, light fluids: +5 ft
- Hot water (> 200°F): +10 ft
- Variable suction conditions: +10 ft
- Slurry service: +15 ft

### Suction Specific Speed (Nss)

Indicates suction performance capability:

```
Nss = N × √Q / (NPSHr)^0.75

Where:
N = rpm
Q = gpm (at impeller eye)
NPSHr = required NPSH in ft
```

**Guidelines**:
- Nss < 8,500: Good suction performance
- Nss = 8,500-11,000: Average
- Nss > 11,000: Poor (cavitation risk)
- Nss > 15,000: High risk, avoid

---

## Application Guidelines

### Municipal Water Supply

**Typical Requirements**:
- Flow: 100 - 10,000 gpm
- Head: 50 - 500 ft
- Efficiency target: > 75%
- Reliability: 99.5%+

**Recommended Pumps**:
1. **Horizontal split-case** (primary choice)
   - Flows: 500-10,000 gpm
   - Easy maintenance (remove top half)
   - Double-suction: low NPSH

2. **Vertical turbine** (for wells)
   - Space-efficient
   - Multi-stage for high head

3. **Submersible** (deep wells)
   - No priming needed
   - Protected from weather

**Design Considerations**:
- Operate at 80-110% of BEP
- Provide standby pump (N+1 redundancy)
- VFD for energy savings
- Soft starters to reduce water hammer

### HVAC/Chilled Water

**Typical Requirements**:
- Flow: 50 - 2,000 gpm
- Head: 30 - 150 ft
- Temperature: 40-180°F
- Constant or variable flow

**Recommended Pumps**:
1. **End-suction centrifugal** (< 500 gpm)
   - Close-coupled or frame-mounted
   - Bronze fitted for corrosion resistance

2. **Inline centrifugal** (space-limited)
   - Vertical orientation
   - Easy pipe mounting

3. **Horizontal split-case** (> 500 gpm)
   - High efficiency
   - Low NPSH

**Design Considerations**:
- Size for design flow + 10-15%
- VFD for variable flow systems
- Isolation valves and gauges
- System curve analysis critical

### Chemical Processing

**Typical Requirements**:
- Flow: 1 - 500 gpm
- Pressure: 50 - 1,000 psi
- Corrosive, toxic, or hazardous fluids
- Temperature: -50 to 400°F

**Recommended Pumps**:
1. **Sealless magnetic drive** (no leaks)
   - Corrosive service
   - Zero emissions

2. **Canned motor** (hermetically sealed)
   - High reliability
   - No shaft seal

3. **Diaphragm metering** (precise dosing)
   - Accurate flow control
   - Handles slurries

4. **PTFE-lined centrifugal** (strong acids)
   - Excellent chemical resistance
   - Lower cost than exotic alloys

**Design Considerations**:
- Material compatibility (corrosion charts)
- Secondary containment
- API 685 for process pumps
- Failure modes analysis

### Oil & Gas Production

**Typical Requirements**:
- Flow: 10 - 2,000 gpm
- Pressure: 500 - 10,000 psi
- Viscous, abrasive, multiphase
- Remote/hazardous locations

**Recommended Pumps**:
1. **Screw pumps** (crude oil transfer)
   - Handles gas entrainment
   - High viscosity capability
   - Gentle handling

2. **Plunger pumps** (injection)
   - Very high pressure
   - Chemical/water injection
   - API 674 rated

3. **Multiphase** (wellhead)
   - Handles gas/liquid mix
   - Reduces separation needs

4. **Centrifugal** (water flood)
   - High volume
   - Reliable operation

**Design Considerations**:
- API 610/674 compliance
- Explosion-proof motors
- Winterization for cold climates
- Remote monitoring/control

### Wastewater/Sewage

**Typical Requirements**:
- Flow: 50 - 10,000 gpm
- Head: 20 - 200 ft
- Solids up to 3" diameter
- Corrosive, abrasive environment

**Recommended Pumps**:
1. **Submersible non-clog** (lift stations)
   - N-impeller or recessed vortex
   - Built-in check valve
   - No priming

2. **Horizontal non-clog** (dry pit)
   - Easy maintenance access
   - Hard iron construction

3. **Progressive cavity** (sludge)
   - High solids content (up to 50%)
   - Gentle, non-pulsating

4. **Grinder pumps** (residential)
   - Macerate solids
   - Smaller diameter piping

**Design Considerations**:
- Minimum 3" solids passage
- Hardened wearing surfaces
- Easily replaceable wear parts
- Clog detection/alarms

### Mining/Slurry

**Typical Requirements**:
- Flow: 100 - 20,000 gpm
- Head: 50 - 1,000 ft
- Solids: 10-70% by weight
- Highly abrasive

**Recommended Pumps**:
1. **Horizontal slurry** (primary)
   - Rubber or metal-lined
   - Recessed impellers
   - Heavy-duty construction

2. **Vertical sump** (pit dewatering)
   - Submerged operation
   - No priming needed

3. **Dredge pumps** (sand/gravel)
   - Large solids passage
   - Centrifugal or cutterhead

**Design Considerations**:
- Slurry velocity: 6-12 ft/s (avoid settling)
- Wear allowances (3-6mm)
- Replaceable liners
- Flushing/cleanout provisions

---

## Manufacturer Resources

### Major Pump Manufacturers

#### Centrifugal Pumps

**Flowserve**
- Website: flowserve.com
- Product lines: Durco, IDP, Worthington
- Selection software: FlowSelect
- Applications: Chemical, oil & gas, power

**Sulzer**
- Website: sulzer.com
- Product lines: CPT, HSB, SMV
- Selection software: Sulzer Select
- Applications: Water, wastewater, industry

**KSB**
- Website: ksb.com
- Product lines: Etanorm, Megabloc, Amarex
- Selection software: KSB FlowManager
- Applications: Water supply, mining, energy

**Grundfos**
- Website: grundfos.com
- Product lines: CR, NB, TP, SEG
- Selection software: Grundfos Product Center (online)
- Applications: HVAC, water supply, wastewater

**Xylem (Goulds)**
- Website: xylem.com/goulds
- Product lines: 3196, 3298, 3400
- Selection software: Xylem Pro
- Applications: Industrial, municipal, building

**ITT Goulds Pumps**
- Website: gouldspumps.com
- Product lines: 3196, AF, SSH
- Standards: ANSI, API 610
- Applications: Chemical, pulp & paper, mining

#### Positive Displacement Pumps

**Viking Pump (IDEX)**
- Website: vikingpump.com
- Product lines: Gear, lobe, vane pumps
- Applications: Viscous fluids, food, chemical

**Moyno (NOV)**
- Website: moyno.com
- Product lines: Progressive cavity
- Applications: Sludge, slurry, viscous

**Wilden (PSG Dover)**
- Website: wildenpump.com
- Product lines: Air-operated diaphragm
- Applications: Chemical, mining, paint

**Cat Pumps**
- Website: catpumps.com
- Product lines: Plunger, piston
- Applications: High pressure cleaning, oil & gas

**Blackmer (PSG Dover)**
- Website: blackmer.com
- Product lines: Sliding vane, diaphragm
- Applications: Truck loading, transfer

### Online Selection Tools

#### Free Web-Based Tools

1. **Pump Selector Pro** (pumpselector.com)
   - Multi-manufacturer comparison
   - Centrifugal and PD pumps
   - Basic specifications

2. **Engineering ToolBox** (engineeringtoolbox.com)
   - Pump calculations
   - Fluid properties
   - Piping design

3. **Pump Fundamentals** (pumpfundamentals.com)
   - Educational resource
   - Selection guides
   - Calculations

4. **Hydraulic Institute** (pumps.org)
   - Industry standards
   - Best practices
   - Training resources

#### Manufacturer-Specific Tools

**Grundfos Product Center** (product-selection.grundfos.com)
- Comprehensive online tool
- All Grundfos products
- Export to CAD/BIM
- Performance curves

**KSB FlowManager** (flowmanager.ksb.com)
- Web-based selection
- KSB product range
- Technical documentation
- Energy analysis

**Flowserve FlowSelect**
- Desktop/online tool
- Chemical compatibility
- Material selection
- Lifecycle cost

### Mobile Apps

1. **Grundfos GO Balance** (iOS/Android)
   - HVAC balancing
   - Pump selection
   - Product lookup

2. **KSB Service Tool** (iOS/Android)
   - Product information
   - Spare parts
   - Service support

3. **Xylem Pro** (iOS/Android)
   - Pump selection
   - Performance curves
   - Troubleshooting

---

## Industry Standards

### Pump Design and Performance

**API 610** - Centrifugal Pumps for Petroleum, Petrochemical and Natural Gas Industries
- Most widely used in oil & gas
- Rigorous design requirements
- OH (overhung), BB (between bearings), VS (vertical)
- Reliability and interchangeability focus

**API 676** - Positive Displacement Pumps - Rotary
- Gear, screw, vane, lobe pumps
- Similar rigor to API 610
- Oil & gas applications

**API 674** - Positive Displacement Pumps - Reciprocating
- Piston, plunger, diaphragm
- High-pressure service
- Instrumentation requirements

**ANSI/HI 9.6.3** - Rotodynamic (Centrifugal and Vertical) Pumps - Guidelines for Allowable Operating Regions
- Operating envelope limits
- Minimum flow requirements
- Protection against damage

**ANSI/HI 9.6.7** - Effects of Liquid Viscosity on Rotodynamic Pump Performance
- Viscosity correction methods
- Performance prediction
- Testing procedures

**ANSI/HI 1.3** - Rotodynamic (Centrifugal) Pumps for Design and Application
- Comprehensive design standard
- Dimensional standards
- Installation guidelines

**ISO 5199** - Technical Specifications for Centrifugal Pumps - Class II
- International standard
- Similar to ANSI B73.1
- End-suction pumps

**ISO 9906** - Rotodynamic Pumps - Hydraulic Performance Acceptance Tests
- Test procedures
- Tolerance grades (1, 2, 3)
- Acceptance criteria

### Testing and Performance

**HI 1.6** - Centrifugal Pump Tests
- Standard test procedures
- Equipment requirements
- Data analysis

**HI 14.6** - Rotodynamic Pumps for Hydraulic Performance Acceptance Tests
- Shop test requirements
- Witness test procedures
- Documentation

**ASME PTC 8.2** - Centrifugal Pumps
- Performance test code
- Precision testing
- Uncertainty analysis

### Materials and Metallurgy

**ASTM Standards**:
- A48: Gray iron castings
- A536: Ductile iron
- A743: Stainless steel castings
- B584: Copper alloy castings
- A276: Stainless steel bar/shapes

**NACE MR0175/ISO 15156** - Petroleum and Natural Gas Industries - Materials for Use in H₂S-Containing Environments
- Sour service requirements
- Material restrictions
- Testing requirements

### Installation and Maintenance

**HI 9.6.4** - Rotodynamic Pumps Guideline for NPSH Margin
- NPSH requirements
- Margin recommendations
- Application-specific guidance

**HI 9.6.6** - Rotodynamic Pumps Guideline for PIping of Pumps
- Suction and discharge piping
- Support requirements
- Flow considerations

**HI 9.8** - Pump Intake Design
- Approach flow conditions
- Anti-vortex devices
- Physical modeling

### Energy Efficiency

**HI 40.6** - Energy Efficiency of Pumping Equipment
- Energy Index (EI)
- Efficiency ratings
- Variable speed drives

**European Union ErP Directive** (2009/125/EC)
- Minimum efficiency index (MEI)
- Motor efficiency requirements
- Implementation timeline

---

## Quick Reference Tables

### Power Calculation Formulas

**Water Horsepower (WHP)**:
```
WHP = (Q × H) / 3960

Where: Q = gpm, H = ft
```

**Brake Horsepower (BHP)**:
```
BHP = WHP / Pump Efficiency

Or: BHP = (Q × H × SG) / (3960 × η)
```

**Motor Size**:
```
Motor HP = BHP × Service Factor (typically 1.15)
```

**SI Units**:
```
Power (kW) = (Q × H × ρ × g) / (3600000 × η)

Where:
Q = m³/h
H = m
ρ = kg/m³
g = 9.81 m/s²
η = efficiency (decimal)
```

### Affinity Laws

For same impeller diameter:

| Parameter | Relationship |
|-----------|--------------|
| Flow | Q₂/Q₁ = N₂/N₁ |
| Head | H₂/H₁ = (N₂/N₁)² |
| Power | P₂/P₁ = (N₂/N₁)³ |

For same speed:

| Parameter | Relationship |
|-----------|--------------|
| Flow | Q₂/Q₁ = (D₂/D₁)³ |
| Head | H₂/H₁ = (D₂/D₁)² |
| Power | P₂/P₁ = (D₂/D₁)⁵ |

### Unit Conversions

| From | To | Multiply by |
|------|-----|-------------|
| gpm | L/s | 0.0631 |
| gpm | m³/h | 0.227 |
| ft | m | 0.3048 |
| psi | bar | 0.0689 |
| psi | kPa | 6.895 |
| HP | kW | 0.746 |
| cP | cSt | 1/SG |
| °F | °C | (°F-32)/1.8 |

### Typical Pump Speeds

| Motor Frequency | Synchronous Speeds (rpm) |
|-----------------|-------------------------|
| 60 Hz | 3600, 1800, 1200, 900 |
| 50 Hz | 3000, 1500, 1000, 750 |

Actual speeds (accounting for slip):
- 2-pole: 3450/2850 rpm
- 4-pole: 1750/1450 rpm
- 6-pole: 1150/960 rpm
- 8-pole: 850/720 rpm

---

## Additional Resources

### Books

1. **"Centrifugal Pumps"** by Johann Friedrich Gülich
   - Comprehensive technical reference
   - Design and theory

2. **"Pump Handbook"** by Igor Karassik et al.
   - Industry bible
   - All pump types covered

3. **"Pumping Machinery Theory and Practice"** by Hassan M. Badr & Wael H. Ahmed
   - Academic textbook
   - Fundamentals and applications

4. **"Practical Centrifugal Pumps"** by Paresh Girdhar & Octo Moniz
   - Design, operation, maintenance
   - Troubleshooting guide

### Training Organizations

- **Hydraulic Institute** (pumps.org) - Certification programs
- **ASME** (asme.org) - Professional development
- **FSA** (fluidsealing.com) - Mechanical seal training
- **Pump School** (mypumpschool.com) - Online courses

### Technical Papers and Journals

- **Pumps & Systems Magazine** (pump-zone.com)
- **World Pumps** (worldpumps.com)
- **ASME Journal of Fluids Engineering**
- **Hydraulic Institute Technical Papers**

---

*Last updated: 2025-01-07*
*For the most current information, always consult manufacturer data sheets and current industry standards.*
