---
name: pump-system-integration
description: "Design complete pump systems including piping, controls, and parallel/series configurations"
category: thinking
domain: mechanical
complexity: advanced
dependencies: [numpy, networkx]
---

# Pump System Integration

## Overview

Comprehensive pump system design integrating pumps with piping networks, controls, and multiple pump configurations for industrial applications.

## System Design Workflow

### 1. System Requirements Analysis

**Define Operating Parameters:**
- Required flow rate range (Q_min, Q_normal, Q_max)
- Total head requirements (static + dynamic)
- Fluid properties (density, viscosity, temperature)
- Operating conditions (continuous, intermittent, variable)
- Redundancy requirements (N+1, duty/standby)
- Environmental constraints

**System Classification:**
- Open vs. closed loop systems
- Constant vs. variable volume systems
- Batch vs. continuous process
- Critical vs. non-critical service

**Performance Metrics:**
- Efficiency targets
- Energy consumption limits
- Lifecycle cost constraints
- Reliability requirements (MTBF)

### 2. System Curve Development

**Static Head Components:**
- Elevation difference: H_static = Z_discharge - Z_suction
- Tank pressure differences: ΔP_tank
- Process pressure requirements

**Dynamic Head Components:**
- Friction losses (Darcy-Weisbach): h_f = f × (L/D) × (v²/2g)
- Minor losses (fittings, valves): h_m = K × (v²/2g)
- Velocity head changes
- Equipment pressure drops

**System Curve Equation:**
```
H_system = H_static + K_system × Q²
where K_system accounts for all friction losses
```

**Variable System Curves:**
- Multiple operating scenarios
- Seasonal variations
- Worst-case design points
- Safety margins (10-15% typical)

### 3. Pump Selection and Matching

**Duty Point Selection:**
- Locate operating point: intersection of pump curve and system curve
- Verify pump efficiency at duty point (within 85-90% of BEP)
- Check NPSH available vs. required (margin ≥ 0.5-1.0 m)
- Ensure operation within preferred operating region (POR)

**Pump Sizing Criteria:**
- Normal flow: 80-110% of BEP flow
- Avoid operation below 70% or above 120% of BEP
- Consider future capacity requirements
- Evaluate turndown ratio for variable flow

**Multiple Pump Options:**
- Single large pump vs. multiple smaller pumps
- Fixed speed vs. variable speed
- Trade-offs: capital cost, efficiency, redundancy, control

### 4. Piping Layout and Sizing

**Velocity Limits:**
- Suction piping: 0.6-1.5 m/s (minimize losses, prevent cavitation)
- Discharge piping: 1.5-3.0 m/s (balance cost vs. friction)
- High-velocity concerns: erosion, noise, water hammer

**Piping Best Practices:**
- Eccentric reducers on suction (flat side up to prevent air pockets)
- Straight pipe runs: 5-10D upstream, 2-3D downstream
- Minimize elbows and fittings near pump
- Support piping to prevent stress on pump nozzles
- Isolation valves for maintenance
- Drain and vent provisions

**Pipe Sizing:**
```
D = √(4Q / πv)
Select standard pipe size (DN, NPS)
Calculate actual velocity and pressure drop
```

**Pressure Drop Calculation:**
- Friction factor (Moody diagram or Colebrook equation)
- Reynolds number: Re = ρvD/μ
- Relative roughness: ε/D
- Total equivalent length: L_eq = L_pipe + ΣL_fittings

### 5. Parallel vs Series Pump Configurations

**Parallel Pumps Configuration:**

*Applications:*
- Variable flow requirements
- Redundancy (N+1 reliability)
- High flow, moderate head systems
- Turndown capability

*Combined Characteristics:*
```
At same head: Q_total = Q_1 + Q_2 + ... + Q_n
Flow distributes based on individual pump curves
```

*Combined H-Q Curve:*
- For identical pumps: double flow at each head
- For different pumps: add flows at constant heads
- Operating point shifts right (higher flow)

*Flow Distribution:*
- Pumps operate at same discharge pressure
- Flow divides according to individual resistances
- Check valve prevents backflow through idle pumps
- Balance piping resistances to equalize flow

*Control Sequencing:*
1. Lead/lag selection (rotate for wear equalization)
2. Staging based on flow demand
3. Minimum runtime between starts
4. Automatic switchover on failure

*Parallel Pump Issues:*
- Unstable operation if pumps have drooping curves
- One pump may dominate if curves differ significantly
- System resistance must be sufficient (avoid runout)
- Need check valves to prevent reverse flow

**Series Pumps Configuration:**

*Applications:*
- High head, moderate flow requirements
- Boosting pressure in long pipelines
- Multi-stage processes
- Overcoming elevation changes

*Head Addition:*
```
At same flow: H_total = H_1 + H_2 + ... + H_n
Heads add at each flow rate
```

*Combined H-Q Curve:*
- For identical pumps: double head at each flow
- For different pumps: add heads at constant flows
- Operating point shifts up (higher head)

*Staging Considerations:*
- Pump order based on head contribution
- Inter-stage pressure limitations
- NPSH available for downstream pumps
- Temperature rise considerations

*Series Pump Issues:*
- Downstream pumps must handle increased pressure
- Check shaft power requirements at all stages
- Inter-stage piping may need pressure class upgrade
- Fluid temperature increases with each stage

### 6. Control Strategies

**Variable Frequency Drive (VFD) Control:**

*Affinity Laws for Speed Variation:*
```
Q_2/Q_1 = N_2/N_1
H_2/H_1 = (N_2/N_1)²
P_2/P_1 = (N_2/N_1)³
```

*Advantages:*
- Energy savings at reduced flow (cubic relationship)
- Smooth flow control
- Soft start reduces mechanical stress
- Eliminates throttling losses

*Disadvantages:*
- Higher capital cost
- Potential harmonic issues
- Minimum speed limitations (cooling, lubrication)
- Not suitable for all applications (steep system curves)

*VFD Energy Savings:*
```
Power savings at reduced flow:
P_ratio = (Q_ratio)³ for VFD
P_ratio = varies for throttling (less efficient)
```

**Throttling Control (Discharge Valve):**

*Method:*
- Increase system resistance by closing valve
- Shifts system curve upward
- Operating point moves left on pump curve

*Characteristics:*
- Simple and low cost
- Energy inefficient (creates artificial resistance)
- Suitable for occasional flow adjustment
- Better for flat system curves

**Bypass Control:**
- Recirculates excess flow back to suction
- Maintains minimum flow through pump
- Prevents overheating at low flows
- Very energy inefficient

**On/Off Control:**
- Simplest method
- Suitable for batch processes
- Pressure or level switches
- Consider cycle frequency limits

**Advanced Control:**
- PID control loops (pressure, flow, level)
- Cascade control systems
- Predictive control algorithms
- System optimization (minimize total energy)

### 7. Transient Analysis (Water Hammer)

**Water Hammer Causes:**
- Rapid valve closure
- Pump startup/shutdown
- Check valve slam
- Air pocket collapse

**Joukowsky Equation (Pressure Surge):**
```
ΔP = ρ × a × Δv
where:
a = wave speed = √(K/ρ) / √(1 + (K/E)(D/t))
K = bulk modulus of fluid
E = elastic modulus of pipe
Δv = velocity change
```

**Critical Closure Time:**
```
T_critical = 2L/a
If closure time < T_critical: maximum pressure surge
If closure time > T_critical: reduced pressure surge
```

**Mitigation Strategies:**
- Slow valve closure (increase closure time)
- Surge tanks or accumulators
- Air chambers (cushioning effect)
- Surge relief valves
- Pump bypass systems
- Controlled pump coastdown

**Analysis Requirements:**
- Transient simulation software (AFT Impulse, WANDA, Hammer)
- Pipe material properties
- Valve closure characteristics
- Pump inertia and performance curves
- Allowable pressure limits (pipe rating ×1.5)

## Parallel Pumps Analysis

### Combined H-Q Curves

**For Identical Pumps:**
```python
# At each head value:
Q_combined = n × Q_single
where n = number of pumps running

# Operating point:
Solve: H_pump(Q_combined/n) = H_system(Q_combined)
```

**For Different Pumps:**
```python
# Create combined curve:
for H in head_range:
    Q_total = Q_pump1(H) + Q_pump2(H) + ...

# Find intersection with system curve
```

**Graphical Method:**
1. Plot individual pump curves
2. At each head, add flows horizontally
3. Plot system curve
4. Intersection = operating point
5. Individual flow = total flow × (individual resistance ratio)

### Flow Distribution

**Equal Resistance Distribution:**
- Identical pumps, symmetric piping: equal flow
- Q_each = Q_total / n

**Unequal Resistance:**
```
Flow distributes inversely with resistance:
Q_1/Q_2 = √(R_2/R_1)

For n pumps with different heads at same flow:
All pumps operate at same discharge pressure P_d
Q_i determined by individual pump curve at P_d
```

**Balancing Techniques:**
- Balance valves in each pump discharge
- Symmetric piping layout (equal length, fittings)
- Flow meters for monitoring
- Pressure taps for verification

### Control Sequencing

**Demand-Based Staging:**
```
Flow ranges for n-pump operation:
1 pump: Q_min to Q_1max
2 pumps: Q_1max to Q_2max
3 pumps: Q_2max to Q_3max

Start additional pump when:
- Current flow > upper limit
- Current pump > 90% capacity
- Pressure drops below setpoint

Stop pump when:
- Flow < lower limit with n pumps
- Redundant capacity available
- Pressure exceeds setpoint
```

**Lead/Lag Rotation:**
- Equalize runtime across pumps
- Alternate lead pump daily/weekly
- Automatic rotation algorithm
- Manual override capability

**Sequencing Logic:**
```
IF system_flow > Q_setpoint_high AND pumps_running < pumps_available:
    START next_pump
    DELAY minimum_time

IF system_flow < Q_setpoint_low AND pumps_running > 1:
    STOP lagging_pump
    DELAY minimum_time

ROTATE lead_pump EVERY rotation_period
```

## Series Pumps Analysis

### Head Addition

**Identical Pumps in Series:**
```python
# At each flow value:
H_combined = n × H_single
where n = number of pumps in series

# Operating point:
Solve: n × H_pump(Q) = H_system(Q)
```

**Different Pumps in Series:**
```python
# Create combined curve:
for Q in flow_range:
    H_total = H_pump1(Q) + H_pump2(Q) + ...

# Find intersection with system curve
```

**Graphical Method:**
1. Plot individual pump curves
2. At each flow, add heads vertically
3. Plot system curve
4. Intersection = operating point
5. Each pump delivers same flow, different heads

### Staging Considerations

**Pump Order:**
- Generally high-head pumps downstream
- Consider NPSH requirements
- Minimize inter-stage pressure
- Temperature rise accumulation

**Inter-stage Pressure:**
```
P_interstage = P_suction + H_1 × ρ × g
Must not exceed:
- Downstream pump casing rating
- Piping pressure class
- Seal pressure limits
```

**NPSH Cascade:**
```
NPSH_available for pump n = P_discharge(n-1) - P_vapor + H_static - H_losses
Verify: NPSH_a > NPSH_r + margin for all stages
```

**Temperature Rise:**
```
ΔT = (H × g × (1 - η)) / (c_p × η)
where:
H = head per pump
η = pump efficiency
c_p = specific heat

Total rise = Σ(ΔT_i) for all pumps
Check: T_final < fluid limits, seal limits
```

**Staging Strategies:**
- All pumps always on (simple, reliable)
- Sequential staging (variable head applications)
- Bypass first stage at low demand
- VFD on first stage, fixed speed on boosters

## Piping Network Analysis

### Network Modeling

**Node-Based Analysis:**
- Nodes: junctions, tanks, pumps
- Links: pipes, valves, fittings
- Continuity at each node: ΣQ = 0
- Energy balance around loops: ΣH = 0

**Hardy-Cross Method:**
```
Iterative solution for complex networks:
1. Assume flow distribution
2. Calculate head loss in each pipe: h = K × Q²
3. Calculate loop corrections: ΔQ = -Σh / (2ΣK|Q|)
4. Update flows: Q_new = Q_old + ΔQ
5. Repeat until convergence
```

**Matrix Methods:**
```
[A] × [H] = [Q]
where:
[A] = network connectivity matrix
[H] = nodal heads
[Q] = nodal demands

Solve using linear algebra
```

### Common Network Configurations

**Branching System:**
- Tree structure, single path to each point
- Simple analysis: sequential calculation
- No loops, no redundancy

**Looped System:**
- Multiple paths between points
- Hardy-Cross or matrix solution required
- Redundancy, better pressure distribution

**Distributed Demand:**
- Withdrawal along pipe length
- Approximate with multiple point loads
- Use equivalent length methods

### Pressure Analysis

**Minimum Pressure Requirements:**
- Process equipment inlet pressures
- End-user pressure requirements
- Elevation effects
- Future expansion margin

**Maximum Pressure Limitations:**
- Pipe pressure ratings (ANSI class)
- Equipment pressure limits
- Valve and fitting ratings
- Safety factors (typically 1.5-2.0)

**Critical Points:**
- Highest elevation (low pressure risk)
- Furthest from pump (highest resistance)
- High demand nodes
- Branch takeoffs

## System Optimization

### Energy Optimization

**Minimize Total Energy Consumption:**
```
E_total = Σ(P_pump × t) + E_losses
where:
P_pump = (ρ × g × Q × H) / η
E_losses = pumping, piping, control losses

Optimize:
- Pump efficiency (operate near BEP)
- Pipe sizing (balance capital vs. energy)
- Control strategy (VFD vs. throttling)
- System curve reduction (improve layout)
```

**Lifecycle Cost Analysis:**
```
LCC = C_capital + C_energy + C_maintenance + C_downtime

C_energy = (P_pump × hours × $/kWh × years)
Typically dominant for continuously operating systems
```

**Optimization Variables:**
- Pump selection (size, type, speed)
- Pipe diameter (velocity vs. friction)
- Number of pumps (parallel/series)
- Control method (VFD, throttle, on/off)

### Reliability Optimization

**Redundancy Configurations:**
- N+1: One standby for N duty pumps
- N+2: Two standbys (critical applications)
- 2×100%: Two pumps, each capable of full duty
- 3×50%: Three pumps, any two handle full load

**Reliability Calculations:**
```
Series reliability: R_series = R_1 × R_2 × ... × R_n
Parallel reliability: R_parallel = 1 - (1-R_1)(1-R_2)...(1-R_n)

Availability = MTBF / (MTBF + MTTR)
```

**Maintenance Strategies:**
- Predictive maintenance (vibration, performance)
- Preventive maintenance schedules
- Condition monitoring
- Spare parts inventory

### Design Optimization Process

1. **Define Objective Function:**
   - Minimize energy cost
   - Minimize capital cost
   - Minimize lifecycle cost
   - Maximize reliability

2. **Identify Constraints:**
   - Flow rate limits
   - Pressure limits
   - Budget constraints
   - Space limitations
   - Standardization requirements

3. **Select Design Variables:**
   - Pump size and number
   - Pipe diameters
   - Control strategy
   - Operating schedules

4. **Optimization Methods:**
   - Parametric studies
   - Gradient-based optimization
   - Genetic algorithms
   - Multi-objective optimization

5. **Sensitivity Analysis:**
   - Vary key parameters
   - Identify critical factors
   - Assess robustness
   - Evaluate risk

## Practical Guidelines

### System Design Checklist

- [ ] System requirements clearly defined
- [ ] System curve developed for all operating scenarios
- [ ] Pump selection optimized for efficiency
- [ ] NPSH margin verified (min 0.5-1.0 m)
- [ ] Piping velocities within limits
- [ ] Parallel/series configuration analyzed
- [ ] Control strategy selected and sized
- [ ] Water hammer analysis completed
- [ ] Pressure ratings verified throughout system
- [ ] Energy consumption calculated
- [ ] Lifecycle cost evaluated
- [ ] Redundancy requirements met
- [ ] Maintenance accessibility verified
- [ ] Instrumentation and monitoring specified

### Common Design Mistakes

1. **Undersized NPSH margin:** Use adequate safety factor
2. **Poor pump selection:** Operating far from BEP
3. **Excessive pipe velocities:** Erosion and noise issues
4. **Inadequate transient analysis:** Water hammer damage
5. **Improper parallel pump application:** Flow instability
6. **Ignoring system curve changes:** Future conditions
7. **Over-reliance on throttling:** Energy waste
8. **Insufficient redundancy:** Reliability issues
9. **Poor piping layout:** Cavitation, air entrainment
10. **Neglecting lifecycle costs:** Focus only on capital cost

### Performance Verification

**Acceptance Testing:**
- Flow rate verification (±5% tolerance)
- Head measurement at duty point
- Power consumption check
- Vibration and noise levels
- NPSH margin confirmation

**Monitoring Points:**
- Suction and discharge pressure
- Flow rate
- Power consumption
- Vibration levels
- Bearing temperature
- Seal condition

**Performance Trending:**
- Track efficiency over time
- Detect degradation early
- Schedule maintenance proactively
- Optimize operating conditions

## Examples

### Example 1: Parallel Pump System Design
See `network-model.py` for complete implementation with:
- Two pumps in parallel serving variable demand
- Combined H-Q curve calculation
- Operating point determination
- Flow distribution analysis
- Control sequencing logic

### Example 2: Series Pump System Design
See `network-model.py` for complete implementation with:
- Two-stage series configuration
- Head addition calculation
- Inter-stage pressure verification
- NPSH cascade analysis
- Temperature rise calculation

### Example 3: Piping Network Analysis
See `network-model.py` for complete implementation with:
- Hardy-Cross network solver
- Multiple pump/demand configuration
- Pressure distribution calculation
- Optimization for pipe sizing

## References

- See `reference.md` for detailed standards and codes
- Hydraulic Institute Standards (HI 9.6.3 - Pump System Design)
- ASME B31.1/B31.3 Piping Codes
- AFT Software - Transient Analysis
- Crane TP-410 - Flow of Fluids

## Notes

- Always verify pump operation within POR (70-120% of BEP)
- Consider future expansion in initial design
- Energy costs typically dominate lifecycle costs
- Reliability requirements drive redundancy decisions
- System optimization requires balanced approach: capital, energy, maintenance
