---
name: positive-displacement-pumps
description: "Design and analyze gear, piston, and screw pumps with volumetric efficiency"
category: thinking
domain: mechanical
complexity: intermediate
dependencies: [numpy]
---

# Positive Displacement Pumps

Positive displacement (PD) pumps move fluid by trapping a fixed volume and forcing it into a discharge pipe. Unlike centrifugal pumps, they deliver nearly constant flow regardless of discharge pressure.

## Types of Positive Displacement Pumps

### Gear Pumps

**External Gear Pumps**
- Two meshing gears rotate in opposite directions
- Fluid trapped between gear teeth and casing
- Simple, reliable, good for clean, viscous fluids
- Flow proportional to speed
- Typical efficiency: 80-95%

**Internal Gear Pumps**
- One gear inside another
- Crescent-shaped seal between gears
- Smoother flow, less pulsation
- Good for viscous fluids
- Compact design

### Piston Pumps (Reciprocating)

**Single-Acting Piston**
- Fluid displaced on one stroke only
- High pulsation
- Simple construction

**Double-Acting Piston**
- Fluid displaced on both strokes
- Reduced pulsation
- Higher efficiency

**Multi-Piston (Triplex, Quintuplex)**
- Multiple pistons offset in phase
- Smoother flow
- Common in high-pressure applications
- Typical efficiency: 85-95%

### Diaphragm Pumps

**Air-Operated Double Diaphragm (AODD)**
- Two flexible diaphragms
- Air pressure drives operation
- Self-priming, can run dry
- Excellent for slurries and solids
- Lower efficiency (~30-70%)

**Mechanically Driven**
- Diaphragm actuated by mechanical linkage
- Higher efficiency than AODD
- Good for metering applications

### Screw Pumps

**Single Screw (Progressive Cavity)**
- Rotor rotates within stator
- Continuous, non-pulsating flow
- Excellent for viscous, shear-sensitive fluids
- Self-priming

**Twin/Triple Screw**
- Two or three intermeshing screws
- Low pulsation
- Good for high-pressure applications
- Typical efficiency: 75-90%

### Lobe Pumps

- Two or more lobes rotate in opposite directions
- Gentle handling of product
- Common in food, pharmaceutical industries
- Easy to clean (sanitary designs)
- Typical efficiency: 50-80%

## Key Characteristics

### Constant Flow Behavior

**Ideal Behavior:**
- Flow rate independent of discharge pressure
- Flow proportional to speed only
- Q = N × V_d

Where:
- Q = volumetric flow rate
- N = pump speed (rpm)
- V_d = displacement per revolution

**Real Behavior:**
- Flow decreases slightly with pressure (slip)
- Efficiency varies with operating conditions

### Volumetric Efficiency

Volumetric efficiency accounts for internal leakage (slip):

η_v = Q_actual / Q_theoretical

η_v = (Q_theoretical - Q_slip) / Q_theoretical

Factors affecting volumetric efficiency:
- Clearances and wear
- Fluid viscosity (higher = better sealing)
- Differential pressure (higher = more leakage)
- Operating speed

### Slip and Leakage

**Slip Flow:**
Q_slip = C × ΔP / μ

Where:
- C = slip coefficient (depends on clearances)
- ΔP = differential pressure
- μ = dynamic viscosity

**Implications:**
- Viscous fluids: less slip, higher efficiency
- High pressures: more slip, lower efficiency
- Worn pumps: increased clearances, more slip

### Pulsation

**Causes:**
- Discrete volume displacement
- Reciprocating motion
- Gear tooth engagement/disengagement

**Pulsation Index:**
PI = (Q_max - Q_min) / Q_avg × 100%

**Typical Pulsation Levels:**
- Single piston: Very high (100%+)
- Duplex piston: High (~50-60%)
- Triplex piston: Moderate (~10-20%)
- Gear pumps: Low to moderate (~5-15%)
- Screw pumps: Very low (<5%)

### Self-Priming Capability

Most PD pumps are self-priming:
- Create vacuum on suction side
- Can lift fluid from below pump
- Can evacuate air from suction line

**Limitations:**
- Maximum suction lift ~8m (limited by atmospheric pressure)
- Requires reasonable seal condition
- May need priming for high-viscosity fluids

## Design Calculations

### Displacement Per Revolution

**Gear Pump:**
V_d = 2 × π × b × (D_o² - D_i²) / 4

Where:
- b = gear width
- D_o = outer diameter
- D_i = inner (root) diameter

**Piston Pump:**
V_d = (π × d² / 4) × L × n

Where:
- d = piston diameter
- L = stroke length
- n = number of pistons (single-acting) or 2n (double-acting)

**Screw Pump:**
V_d = 4 × A_c × p

Where:
- A_c = cavity area
- p = pitch

### Theoretical Flow Rate

Q_theoretical = V_d × N / 60

Where:
- Q_theoretical in m³/s or L/min
- V_d in m³ or L
- N in rpm

### Actual Flow Rate (Accounting for Slip)

Q_actual = η_v × Q_theoretical

Q_actual = Q_theoretical - Q_slip

**Slip as function of pressure and viscosity:**
Q_actual = Q_theoretical - (C × ΔP / μ)

**For design:**
- Specify desired flow at operating pressure
- Account for expected volumetric efficiency
- Select pump with adequate theoretical capacity

### Power Requirements

**Hydraulic Power:**
P_hydraulic = Q × ΔP

Where:
- P in Watts
- Q in m³/s
- ΔP in Pa

**Brake Power (Shaft Power):**
P_brake = P_hydraulic / η_overall

η_overall = η_v × η_m

Where:
- η_v = volumetric efficiency
- η_m = mechanical efficiency (bearings, seals)

**Typical Overall Efficiencies:**
- Gear pumps: 70-85%
- Piston pumps: 80-90%
- Screw pumps: 70-85%
- Diaphragm pumps: 30-70%

**Motor Power (with safety factor):**
P_motor = P_brake × SF

SF typically 1.15-1.25

### NPSH Requirements

PD pumps generally require lower NPSH than centrifugal pumps:

NPSH_required = P_atm/ρg - h_suction - h_friction - P_vapor/ρg - safety_margin

**Typical NPSH_required:**
- 0.5-2 m for most PD pumps
- Higher for high-speed pumps
- Lower for slow-speed pumps

## When to Use PD vs Centrifugal Pumps

### Choose Positive Displacement When:

1. **High-Pressure Applications**
   - ΔP > 10-20 bar
   - PD pumps maintain efficiency at high pressure
   - Centrifugal pumps become impractical

2. **Viscous Fluids**
   - μ > 100 cP
   - PD efficiency improves with viscosity
   - Centrifugal efficiency drops dramatically

3. **Constant Flow Required**
   - Metering and dosing
   - Flow independent of pressure variations
   - Predictable delivery

4. **Low Flow, High Pressure**
   - Centrifugal pumps inefficient at low flow
   - PD pumps excel in this range

5. **Self-Priming Required**
   - Suction lift needed
   - Air entrainment possible
   - Dry-run capability

6. **Shear-Sensitive Fluids**
   - Food products, polymers
   - Use lobe or progressive cavity pumps
   - Gentle handling

### Choose Centrifugal When:

1. **High Flow, Low Pressure**
   - Q > 100 m³/h, ΔP < 10 bar
   - More economical
   - Simpler maintenance

2. **Low Viscosity Fluids**
   - μ < 50 cP (water-like)
   - Centrifugal pumps efficient
   - Less expensive

3. **Continuous, Smooth Flow**
   - No pulsation acceptable
   - Variable flow needed
   - Throttling control

4. **Particulate Handling**
   - Large solids
   - PD pumps can jam
   - Centrifugal more forgiving

5. **Lower Initial Cost**
   - Simple installation
   - Standard motors
   - Lower maintenance

## Pulsation Dampening

Pulsation can cause:
- Vibration and noise
- Inaccurate flow measurement
- Pressure spikes
- System fatigue

### Dampening Methods

**1. Pulsation Dampener (Accumulator)**

**Gas-Charged Bladder Type:**
- Bladder separates gas and fluid
- Gas compresses during pressure peaks
- Gas expands during pressure valleys
- Smooth flow output

**Sizing:**
V_dampener = (Q_theoretical × C_d) / (η_p × f)

Where:
- C_d = dampening coefficient (typically 5-10)
- η_p = pulsation reduction efficiency (0.9-0.95)
- f = pump frequency (Hz)

**Gas Pre-Charge Pressure:**
P_precharge = 0.6 × P_operating (typical)

**2. Multiple Pistons**

**Triplex Pump (3 pistons at 120°):**
- Pulsation reduced ~90%
- Common in high-pressure applications

**Quintuplex Pump (5 pistons at 72°):**
- Pulsation reduced ~95%
- Smoother than triplex

**3. Air Chambers**

Simple expansion chamber on discharge:
- Gas cushion absorbs pulsation
- Requires regular air charging
- Lower cost than bladder type

**4. Flexible Discharge Line**

- Hose instead of rigid pipe (short section)
- Elasticity absorbs pulses
- Simple, low cost
- Limited effectiveness

**5. Flow Stabilizer**

- Restrictor orifice
- Creates back pressure
- Dampens pressure fluctuations
- Energy loss

### Design Considerations

**Critical Frequencies:**
Avoid resonance with system natural frequency:

f_pump = N × n / 60

Where:
- N = pump speed (rpm)
- n = number of pistons or pumping chambers

**Dampener Location:**
- Close to pump discharge
- Before flow meter (if smooth flow required)
- Consider accessibility for maintenance

**Maintenance:**
- Check bladder integrity
- Verify gas pre-charge pressure
- Inspect for leaks

## Application Selection Guide

### High-Pressure Chemical Injection
**Recommended:** Plunger pump (triplex)
- High pressure capability (up to 1000+ bar)
- Good metering accuracy
- Pulsation dampener required

### Viscous Oil Transfer
**Recommended:** Gear pump or screw pump
- Handles high viscosity well
- Self-priming
- Relatively smooth flow

### Slurry and Solids
**Recommended:** Diaphragm pump or progressive cavity
- Can handle solids
- Won't jam easily
- Gentle action

### Food and Pharmaceutical
**Recommended:** Lobe pump or sanitary diaphragm
- Hygienic design
- Easy to clean
- Gentle product handling

### Metering and Dosing
**Recommended:** Diaphragm metering pump or plunger pump
- Excellent accuracy
- Adjustable stroke
- Handles chemicals

### General Water Transfer
**Recommended:** Centrifugal pump
- Lower cost
- Lower maintenance
- Adequate for low viscosity

## Design Example Workflow

1. **Define Requirements:**
   - Flow rate (actual, at operating pressure)
   - Differential pressure
   - Fluid properties (density, viscosity)
   - Operating conditions

2. **Select Pump Type:**
   - Based on application guide above
   - Consider fluid compatibility
   - Evaluate cost constraints

3. **Calculate Theoretical Capacity:**
   - Account for volumetric efficiency
   - Q_theoretical = Q_actual / η_v

4. **Size Displacement:**
   - V_d = Q_theoretical × 60 / N
   - Choose operating speed N

5. **Calculate Power:**
   - P_hydraulic = Q × ΔP
   - P_brake = P_hydraulic / η_overall
   - P_motor = P_brake × SF

6. **Check NPSH:**
   - NPSH_available > NPSH_required + margin
   - Size suction piping appropriately

7. **Evaluate Pulsation:**
   - Calculate pulsation index
   - Determine if dampening needed
   - Size dampener if required

8. **Verify Operating Range:**
   - Minimum/maximum speed
   - Pressure limitations
   - Viscosity range

## Performance Monitoring

**Key Parameters to Track:**
- Flow rate vs. speed (check for increased slip)
- Discharge pressure
- Power consumption (increased = wear)
- Vibration levels
- Temperature

**Indicators of Wear:**
- Reduced flow at same speed
- Increased power consumption
- Increased noise/vibration
- Reduced volumetric efficiency

**Maintenance Planning:**
- Replace seals/gaskets per schedule
- Monitor clearances in gear pumps
- Check valve seats in piston pumps
- Inspect diaphragms regularly

## Summary

Positive displacement pumps are essential for:
- High-pressure applications
- Viscous fluid handling
- Metering and constant flow
- Self-priming requirements

Key design considerations:
- Account for volumetric efficiency (slip)
- Size for actual flow needed
- Consider pulsation dampening
- Match pump type to application

Trade-offs vs. centrifugal:
- Higher pressure capability
- Better viscosity handling
- Pulsating flow
- Higher initial cost
- More maintenance
