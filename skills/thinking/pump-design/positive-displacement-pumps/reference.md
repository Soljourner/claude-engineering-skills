# Positive Displacement Pumps - Reference Guide

## Design Equations

### General PD Pump Relationships

**Theoretical Flow Rate:**
```
Q_theoretical = V_d × N / 60
```
- Q_theoretical: volumetric flow rate (m³/s or L/min)
- V_d: displacement per revolution (m³ or cm³)
- N: pump speed (rpm)

**Actual Flow Rate:**
```
Q_actual = η_v × Q_theoretical
Q_actual = Q_theoretical - Q_slip
```
- η_v: volumetric efficiency (0-1)
- Q_slip: slip flow due to leakage

**Slip Flow (Viscosity Dependent):**
```
Q_slip = C × ΔP / μ
```
- C: slip coefficient (m³/(Pa·s)), depends on clearances
- ΔP: differential pressure (Pa)
- μ: dynamic viscosity (Pa·s)

**Volumetric Efficiency:**
```
η_v = Q_actual / Q_theoretical
η_v = 1 - (Q_slip / Q_theoretical)
```

### Power Calculations

**Hydraulic Power:**
```
P_hydraulic = Q × ΔP
```
- P: power (W)
- Q: flow rate (m³/s)
- ΔP: pressure rise (Pa)

**Brake Power (Shaft Power):**
```
P_brake = P_hydraulic / η_overall
η_overall = η_v × η_m
```
- η_m: mechanical efficiency (accounts for friction)

**Motor Power:**
```
P_motor = P_brake × SF
```
- SF: safety factor (typically 1.15-1.25)

**Torque:**
```
T = P_brake × 60 / (2π × N)
```
- T: torque (N·m)
- N: speed (rpm)

### Gear Pump Equations

**External Gear Pump Displacement:**
```
V_d = 2 × π × b × (D_o² - D_i²) / 4
```
- b: gear width (face width)
- D_o: outer diameter (tip diameter)
- D_i: inner diameter (root diameter)

**Simplified (with module):**
```
V_d = 2 × π × b × m × z × h
```
- m: module
- z: number of teeth
- h: tooth height

**Internal Gear Pump Displacement:**
```
V_d = π × b × (D_o² - D_i²) / 4
```
(Approximately half that of external gear with same dimensions)

**Gear Pump Slip Coefficient:**
```
C = k × (h³ × L_clearance) / (12 × L_seal)
```
- h: radial clearance
- L_clearance: length of clearance path
- L_seal: sealing length
- k: geometry factor (typically 1-2)

### Piston Pump Equations

**Single Piston Displacement:**
```
V_piston = (π × d² / 4) × L
```
- d: piston diameter
- L: stroke length

**Single-Acting Pump:**
```
V_d = n × V_piston
```
- n: number of pistons

**Double-Acting Pump:**
```
V_d = n × [V_piston + (π × (d² - d_rod²) / 4) × L]
```
- d_rod: piston rod diameter

**Piston Velocity:**
```
v = ω × (L/2) × sin(ωt)
```
- ω: angular velocity (rad/s)
- t: time (s)

**Instantaneous Flow (Single Piston):**
```
Q(t) = A × v = (π × d² / 4) × ω × (L/2) × sin(ωt)
```
(Only positive values during discharge stroke)

**Multi-Piston Total Flow:**
```
Q_total(t) = Σ Q_i(t - t_i)
```
- t_i: phase offset for piston i

### Screw Pump Equations

**Single Screw (Progressive Cavity):**
```
V_d = 4 × e × D × p
```
- e: eccentricity
- D: rotor diameter
- p: pitch

**Alternative form:**
```
V_d = 4 × A_cavity × p
```
- A_cavity: cross-sectional area of cavity

**Twin Screw:**
```
V_d = π × L_effective × (D² - d²) × tan(α)
```
- L_effective: effective screw length
- D: outer diameter
- d: inner diameter
- α: helix angle

**Pressure Capability (Progressive Cavity):**
```
ΔP_max = n_stages × ΔP_per_stage
```
- n_stages: number of stages (pitch lengths)
- ΔP_per_stage: typically 0.5-1 bar for standard elastomers

### Diaphragm Pump Equations

**AODD (Air Operated Double Diaphragm):**
```
Q_theoretical = 2 × V_chamber × f
```
- V_chamber: chamber volume
- f: cycle frequency (Hz)

**Air Consumption:**
```
Q_air = Q_fluid × (P_discharge + P_atm) / (P_air - P_atm) × k
```
- k: inefficiency factor (typically 1.5-2.0)

**Mechanical Diaphragm:**
```
V_d = (π × d² / 4) × s
```
- d: diaphragm diameter
- s: diaphragm stroke

### Pulsation Analysis

**Pulsation Index:**
```
PI = (Q_max - Q_min) / Q_avg × 100%
```

**Fundamental Pulsation Frequency:**
```
f_pulsation = N × n / 60
```
- N: pump speed (rpm)
- n: number of pumping chambers

**For Reciprocating Pumps:**

**Simplex (1 piston):**
```
PI ≈ 100-150%
f = N/60 Hz
```

**Duplex (2 pistons):**
```
PI ≈ 50-60%
f = 2N/60 Hz
```

**Triplex (3 pistons):**
```
PI ≈ 10-20%
f = 3N/60 Hz
```

**Quintuplex (5 pistons):**
```
PI ≈ 5-8%
f = 5N/60 Hz
```

### Pulsation Dampener Sizing

**Dampener Volume:**
```
V_dampener = (Q_theoretical × C_d) / (η_p × f)
```
- C_d: dampening coefficient (5-10)
- η_p: pulsation reduction efficiency (0.90-0.95)
- f: pulsation frequency (Hz)

**Gas Pre-Charge Pressure:**
```
P_precharge = k × P_operating
```
- k: typically 0.6-0.8

**Pressure Ratio (for bladder accumulators):**
```
V_gas_1 / V_gas_2 = (P_2 / P_1)^(1/n)
```
- n: polytropic exponent (1.0 for isothermal, 1.4 for adiabatic)
- Subscript 1: min pressure state
- Subscript 2: max pressure state

### NPSH Calculations

**NPSH Available:**
```
NPSH_a = P_atm/ρg + h_static - h_friction - P_vapor/ρg
```
- P_atm: atmospheric pressure
- h_static: static head (positive if fluid above pump)
- h_friction: friction losses in suction line
- P_vapor: vapor pressure of fluid
- ρ: density
- g: gravitational acceleration

**NPSH Required (typical):**
- Gear pumps: 0.5-2 m
- Piston pumps: 1-3 m
- Screw pumps: 0.5-1.5 m
- Diaphragm pumps: 0.5-2 m

**Safety Margin:**
```
NPSH_a ≥ NPSH_r + 0.5 m
```

## Efficiency Correlations

### Volumetric Efficiency

**Gear Pumps:**
```
η_v = 1 - (k_1 × ΔP) / (μ × N)
```
- k_1: constant depending on clearances (typically 1e-8 to 1e-9)

**Typical ranges:**
- Water (low viscosity): 70-85%
- Hydraulic oil: 85-95%
- Heavy oils: 90-95%

**Piston Pumps:**
```
η_v = 1 - (k_2 × ΔP) / (μ × N)
```
- k_2: valve leakage constant

**Typical ranges:**
- Low pressure (<100 bar): 92-98%
- High pressure (>200 bar): 88-95%

**Screw Pumps:**
```
η_v = 1 - (k_3 × ΔP) / (μ × N)
```

**Progressive cavity typical:**
- Thin fluids: 85-92%
- Viscous fluids: 90-95%
- With solids: 80-90%

**Twin/triple screw typical:**
- 85-95% across wide range

### Mechanical Efficiency

**Bearing and seal friction:**
```
η_m = 1 - (P_friction / P_hydraulic)
```

**Typical values:**
- Gear pumps: 85-92%
- Piston pumps: 88-93%
- Screw pumps: 85-90%
- Diaphragm pumps: 60-80% (AODD), 80-90% (mechanical)

### Overall Efficiency

```
η_overall = η_v × η_m
```

**Typical overall efficiencies:**
- Gear pumps: 70-85%
- Piston pumps: 80-92%
- Screw pumps: 70-85%
- AODD pumps: 30-60%
- Mechanical diaphragm: 70-80%

### Efficiency vs. Operating Conditions

**Effect of Speed:**
```
η_v(N) = η_v0 - k_s × (N - N_0)
```
- Higher speed → slightly lower volumetric efficiency (more turbulence)
- Lower speed → lower mechanical efficiency (higher relative friction)
- Optimum typically at 60-80% of maximum speed

**Effect of Pressure:**
```
η_v(ΔP) = η_v0 × [1 - (ΔP / ΔP_crit)^m]
```
- m: empirical exponent (typically 0.5-1.0)
- ΔP_crit: pressure at which efficiency drops significantly

**Effect of Viscosity:**
```
η_v(μ) = η_v_min + (η_v_max - η_v_min) × [1 - exp(-μ / μ_0)]
```
- μ_0: characteristic viscosity (typically 0.01-0.05 Pa·s)
- Higher viscosity → better sealing → higher volumetric efficiency

## Application Guidelines

### Selection Matrix

#### By Viscosity Range

**Low Viscosity (< 10 cP)**
- Piston pumps: Excellent
- Gear pumps: Fair to Good
- Screw pumps: Good
- AODD: Good
- Centrifugal: Preferred (if pressure allows)

**Medium Viscosity (10-500 cP)**
- Piston pumps: Excellent
- Gear pumps: Excellent
- Screw pumps: Excellent
- AODD: Good
- Centrifugal: Poor

**High Viscosity (500-10,000 cP)**
- Piston pumps: Good (with heated cylinder)
- Gear pumps: Excellent
- Screw pumps: Excellent
- AODD: Fair
- Centrifugal: Not suitable

**Very High Viscosity (> 10,000 cP)**
- Gear pumps: Excellent (at low speed)
- Screw pumps: Excellent
- Others: Not recommended

#### By Pressure Range

**Low Pressure (< 10 bar)**
- All PD pumps: Suitable
- Centrifugal: Often more economical
- Consider: Initial cost, efficiency

**Medium Pressure (10-100 bar)**
- Gear pumps: Excellent
- Piston pumps: Excellent
- Screw pumps: Good
- Centrifugal: Possible but inefficient

**High Pressure (100-400 bar)**
- Piston pumps (plunger): Preferred
- Gear pumps: Limited to ~250 bar
- Screw pumps: Not typical

**Very High Pressure (> 400 bar)**
- Plunger pumps: Only practical option
- Multiple stages may be required

#### By Flow Rate

**Very Low Flow (< 1 L/min)**
- Metering pumps (diaphragm): Excellent
- Small gear pumps: Good
- Peristaltic: Good for sterile applications

**Low Flow (1-50 L/min)**
- Gear pumps: Excellent
- Piston pumps: Excellent
- Screw pumps: Good

**Medium Flow (50-500 L/min)**
- Gear pumps: Good (multiple units or large size)
- Piston pumps: Good
- Screw pumps: Excellent
- Lobe pumps: Excellent

**High Flow (> 500 L/min)**
- Screw pumps: Good
- Lobe pumps: Good
- Centrifugal: Often preferred (if suitable)

#### By Application Type

**Chemical Injection / Metering**
- Preferred: Diaphragm metering pump or plunger pump
- Accuracy: ±1-2%
- Features: Adjustable stroke, pulsation dampener

**High-Pressure Cleaning**
- Preferred: Triplex plunger pump
- Pressure: 100-400 bar typical
- Features: Ceramic plungers, high-pressure seals

**Food Processing**
- Preferred: Lobe pump or sanitary progressive cavity
- Requirements: 3-A sanitary design, CIP capability
- Features: FDA-approved elastomers

**Viscous Oils / Polymers**
- Preferred: Gear pump or twin-screw
- Features: Jacketed (heated), low-shear

**Slurries / Solids**
- Preferred: AODD or progressive cavity
- Can handle: Up to 50% solids by volume
- Features: Large port sizes, wear-resistant materials

**Paint / Coatings**
- Preferred: Diaphragm or piston pump
- Requirements: Pulsation dampening, precise control
- Features: Easy color change, cleanability

**Wastewater / Sludge**
- Preferred: Progressive cavity or AODD
- Can handle: Solids, rags, debris
- Features: Heavy-duty construction, easy maintenance

### Material Selection

#### Wetted Materials

**For water and mild chemicals:**
- Cast iron (economics)
- Bronze or brass (corrosion resistance)
- 304 Stainless steel

**For corrosive chemicals:**
- 316 Stainless steel (standard)
- Alloy 20 (severe service)
- Hastelloy (highly corrosive)
- PTFE/PVDF lining (cost-effective alternative)

**For abrasive slurries:**
- Hardened steel
- Chrome-plated
- Tungsten carbide (severe abrasion)
- Ceramic (ultra-abrasive)

**For food / pharmaceutical:**
- 316L Stainless steel (electropolished)
- FDA-approved elastomers
- 3-A sanitary certification

#### Seal Materials

**General purpose:**
- Nitrile (NBR): Hydraulic oils, petroleum
- EPDM: Water, steam, mild chemicals
- Viton (FKM): Wide chemical resistance

**Special applications:**
- PTFE: Universal chemical resistance
- Kalrez: Extreme chemicals and temperature
- UHMWPE: Abrasive services

#### Elastomer Selection (for cavity/lobe pumps)

**Temperature range:**
- Nitrile: -40°C to +100°C
- EPDM: -40°C to +150°C
- Viton: -20°C to +200°C

**Chemical compatibility:**
- Consult manufacturer charts
- Consider swelling, hardness change
- Test samples in actual fluid if critical

### Installation Guidelines

#### Suction Line Design

**Minimize NPSH losses:**
```
h_friction = f × (L/D) × (v²/2g)
```

**Best practices:**
- Short, direct route
- Minimum bends (use large radius)
- One size larger than pump inlet
- Upward slope to pump (no high points)
- Strainer with adequate area (low velocity)

**Velocity limits:**
- Suction line: < 1.5 m/s (< 0.6 m/s for high viscosity)
- Discharge line: < 3 m/s

#### Piping Recommendations

**Pressure pulsation considerations:**
- Support piping independently (not on pump)
- Flexible connectors at pump inlet/outlet
- Secure piping to prevent vibration
- Pulsation dampener close to pump

**For viscous fluids:**
- Heat tracing if needed
- Insulation to maintain temperature
- Large bore to minimize friction loss

#### Relief Valve Sizing

**Always required for PD pumps:**
```
Relief valve setting = 1.1 × P_design
```

**Capacity:**
```
Q_relief ≥ 1.1 × Q_max
```

**Location:**
- Between pump and first isolation valve
- Discharge back to source or surge tank

#### Foundation and Mounting

**Rigid mounting for high-pressure pumps:**
- Concrete foundation or heavy steel base
- Grouted baseplates
- Alignment within 0.05 mm

**Flexible mounting for AODD pumps:**
- Rubber isolators
- Flexible hose connections

### Maintenance Planning

#### Inspection Intervals

**Daily:**
- Visual check for leaks
- Listen for unusual noise
- Check operating pressure

**Weekly:**
- Check lubrication level
- Monitor performance (flow, pressure, power)
- Inspect for vibration

**Monthly:**
- Check seal/packing condition
- Verify relief valve operation
- Review trend data

**Quarterly:**
- Detailed inspection
- Check alignment
- Test instrumentation
- Change filters

**Annually:**
- Major overhaul (as needed)
- Replace wear parts
- Recalibrate instruments

#### Performance Monitoring

**Key indicators of wear:**

**Reduced volumetric efficiency:**
```
η_v_current = Q_actual_current / Q_theoretical
If: η_v_current < 0.85 × η_v_rated → Plan overhaul
```

**Increased power consumption:**
```
If: P_current > 1.15 × P_normal → Investigate
```

**Possible causes:**
- Worn internal clearances
- Valve seat wear (piston pumps)
- Seal degradation
- Bearing wear

#### Spare Parts Inventory

**Critical spares (always stock):**
- Seal/packing sets
- Valves (for piston pumps)
- Diaphragms (for diaphragm pumps)
- Bearings

**Major components (as needed):**
- Complete rotor assemblies
- Stators (for progressive cavity)
- Gear sets

### Troubleshooting Guide

#### Insufficient Flow

**Possible causes:**
1. Speed too low → Check motor speed
2. Internal wear → Measure slip, plan overhaul
3. Air in suction → Check NPSH, seal
4. Viscosity too high → Heat fluid or reduce speed
5. Relief valve leaking → Test and replace

#### Excessive Noise/Vibration

**Possible causes:**
1. Cavitation → Check NPSH, reduce speed
2. Air entrainment → Check suction line
3. Misalignment → Re-align coupling
4. Bearing wear → Replace bearings
5. Resonance → Check mounting, add dampener

#### Overheating

**Possible causes:**
1. Excessive pressure → Check relief valve
2. High viscosity → Reduce speed or heat fluid
3. Insufficient cooling → Check jacket flow
4. Internal recirculation → Check for wear

#### Seal/Packing Leaks

**Possible causes:**
1. Normal wear → Replace seals
2. Improper installation → Reinstall correctly
3. Shaft wear → Repair or replace shaft
4. Excessive pressure → Check system

#### Pressure Pulsation

**Solutions:**
1. Install pulsation dampener
2. Increase number of pistons
3. Add flexible discharge line
4. Check for air in system

## Design Checklist

### Initial Selection

- [ ] Define required flow rate (at operating pressure)
- [ ] Define pressure requirements (min, max, normal)
- [ ] Characterize fluid (density, viscosity, temperature)
- [ ] Identify any solids or abrasives
- [ ] Determine required accuracy/consistency
- [ ] Consider space and installation constraints
- [ ] Establish budget (capital and operating)

### Pump Sizing

- [ ] Select pump type based on application
- [ ] Account for volumetric efficiency
- [ ] Calculate required displacement
- [ ] Select operating speed
- [ ] Verify within manufacturer's range
- [ ] Check turndown ratio if variable flow needed
- [ ] Size motor with adequate safety factor

### System Design

- [ ] Calculate NPSH available
- [ ] Verify NPSH_a > NPSH_r + margin
- [ ] Size suction piping (low velocity)
- [ ] Size discharge piping
- [ ] Select relief valve and size properly
- [ ] Evaluate need for pulsation dampener
- [ ] Design foundation/mounting

### Material Selection

- [ ] Verify chemical compatibility
- [ ] Check temperature limits
- [ ] Consider abrasion resistance
- [ ] Select appropriate seals/elastomers
- [ ] Ensure food-grade if required
- [ ] Specify surface finish if sanitary

### Instrumentation

- [ ] Flow measurement (if needed)
- [ ] Pressure gauges (suction and discharge)
- [ ] Temperature sensors (if heated)
- [ ] Vibration monitoring (critical services)
- [ ] Level switches (source and destination)

### Documentation

- [ ] P&ID showing all components
- [ ] Equipment data sheets
- [ ] Installation drawings
- [ ] Operating procedures
- [ ] Maintenance schedule
- [ ] Spare parts list

## Key References

### Standards

- **API 674**: Positive Displacement Pumps - Reciprocating
- **API 676**: Positive Displacement Pumps - Rotary
- **ISO 9906**: Rotodynamic and positive displacement pumps
- **ASME B73.1**: Horizontal end suction centrifugal pumps
- **3-A Sanitary Standards**: For food/pharma applications

### Typical Performance Ranges

```
Pump Type        | Flow Range     | Pressure Range | Viscosity Range
-----------------+----------------+----------------+------------------
Gear (external)  | 1-500 L/min   | 10-250 bar    | 1-10,000 cP
Gear (internal)  | 1-200 L/min   | 5-150 bar     | 1-10,000 cP
Piston (plunger) | 0.1-500 L/min | 50-1500 bar   | 0.5-1000 cP
Diaphragm (AODD) | 5-500 L/min   | 0-8 bar       | 0.5-5000 cP
Diaphragm (mech) | 0.01-100 L/min| 1-400 bar     | 0.5-2000 cP
Prog. cavity     | 1-100 m³/h    | 0-48 bar      | 1-100,000 cP
Twin screw       | 5-500 m³/h    | 5-100 bar     | 1-50,000 cP
Lobe             | 5-500 m³/h    | 1-15 bar      | 1-10,000 cP
```

### Conversion Factors

```
Pressure:
  1 bar = 100 kPa = 14.5 psi = 10.2 m H₂O

Viscosity:
  1 cP = 1 mPa·s
  1 cSt × density (g/cm³) = cP

Flow:
  1 m³/h = 16.67 L/min = 4.4 US gpm

Power:
  1 hp = 0.746 kW
```
