# Centrifugal Pump Design - Technical Reference

## Euler Turbine Equation Derivation

### Fundamental Principle

The Euler turbine equation describes energy transfer between a rotating impeller and fluid. It is derived from the angular momentum theorem applied to steady flow through a control volume.

### Derivation

**Step 1: Angular Momentum Balance**

For steady flow through a rotating impeller, the torque equals the rate of change of angular momentum:

```
T = ṁ · (r₂·c_u₂ - r₁·c_u₁)
```

Where:
- T = torque on the shaft (N·m)
- ṁ = mass flow rate (kg/s)
- r = radius (m)
- c_u = tangential (whirl) component of absolute velocity (m/s)
- Subscripts 1, 2 denote inlet and outlet

**Step 2: Power Transfer**

Power transmitted to the fluid:

```
P = T·ω = ṁ·ω·(r₂·c_u₂ - r₁·c_u₁)
```

Where ω = angular velocity (rad/s)

**Step 3: Blade Velocity**

The blade velocity at any radius:

```
u = ω·r
```

Therefore:

```
P = ṁ·(u₂·c_u₂ - u₁·c_u₁)
```

**Step 4: Head Calculation**

Power can also be expressed as:

```
P = ṁ·g·H
```

Where:
- g = gravitational acceleration (m/s²)
- H = head (m)

**Step 5: Euler Equation**

Equating the two expressions for power:

```
ṁ·g·H = ṁ·(u₂·c_u₂ - u₁·c_u₁)
```

Simplifying:

```
H = (u₂·c_u₂ - u₁·c_u₁) / g
```

**This is the Euler turbine equation.**

### Special Cases

**Case 1: Radial Entry (No Pre-rotation)**

If the flow enters radially with no tangential component:

```
c_u₁ = 0
```

Then:

```
H = u₂·c_u₂ / g
```

**Case 2: Backward-Curved Blades**

From velocity triangle geometry:

```
c_u₂ = u₂ - c_m₂/tan(β₂)
```

Substituting:

```
H = u₂·[u₂ - c_m₂/tan(β₂)] / g = u₂²/g - u₂·c_m₂/(g·tan(β₂))
```

This shows that head consists of:
1. Centrifugal component: u₂²/g
2. Blade angle component: -u₂·c_m₂/(g·tan(β₂))

**Case 3: Radial Blades**

If β₂ = 90°, then tan(β₂) → ∞ and:

```
H = u₂²/g
```

Maximum theoretical head, but unstable operation.

### Theoretical vs. Actual Head

**Slip Factor**

Real pumps have finite number of blades, causing "slip" - the fluid doesn't follow blade angles exactly.

Wiesner slip factor:

```
σ = 1 - √(sin(β₂))/Z^0.7
```

Or simplified Stodola:

```
σ = 1 - π·sin(β₂)/Z
```

Where Z = number of blades

**Actual Head:**

```
H_actual = σ·H_theoretical
```

Typical slip factors: σ = 0.85-0.95

## Specific Speed

### Definition and Physical Meaning

Specific speed (Ns) is a dimensionless parameter characterizing pump geometry and type. It represents the speed at which a geometrically similar pump would operate to deliver 1 m³/s at 1 m head.

### Formulations

**European Convention (nq):**

```
nq = N·√Q / H^(3/4)
```

Units:
- N in rpm
- Q in m³/s
- H in m (per stage)

**US Convention (Ns):**

```
Ns = N·√(Q_gpm) / H_ft^(3/4)
```

Units:
- N in rpm
- Q in gpm
- H in ft

**Conversion:**

```
Ns(US) = 51.6·nq(European)
```

Or:

```
Ns(US) = 2733·n_s(dimensionless)
```

### Derivation from Dimensional Analysis

**Step 1: Relevant Parameters**

For a pump operating at peak efficiency:
- Q = flow rate
- H = head
- N = speed
- D = impeller diameter
- ρ = fluid density
- μ = fluid viscosity
- g = gravity

**Step 2: Dimensionless Groups**

From Buckingham Pi theorem, three main groups:

```
Φ = Flow coefficient = Q/(N·D³)
Ψ = Head coefficient = g·H/(N²·D²)
Re = Reynolds number = ρ·N·D²/μ
```

**Step 3: Specific Speed Derivation**

At peak efficiency, pumps of the same family have constant Φ/Ψ^(3/4).

Eliminating diameter:

```
Ns = N·√Q / H^(3/4) = constant
```

This is independent of pump size but characterizes pump type.

### Specific Speed Ranges

**Table 1: Pump Types by Specific Speed**

| nq (European) | Ns (US) | Pump Type | Characteristics |
|---------------|---------|-----------|----------------|
| 10-25 | 500-1300 | Slow, radial | Narrow impeller, high head |
| 25-40 | 1300-2000 | Normal, radial | Medium width, good efficiency |
| 40-60 | 2000-3000 | Francis-vane | Wide impeller, very efficient |
| 60-100 | 3000-5000 | Mixed flow | Diagonal flow, high capacity |
| 100-200 | 5000-10000 | Mixed/axial | Very wide, high flow |
| 200-400 | 10000-20000 | Axial flow | Propeller type, low head |

**Physical Interpretation:**

- **Low Ns:** High head, low flow → narrow, radial impellers
- **Medium Ns:** Balanced head and flow → Francis-vane impellers
- **High Ns:** Low head, high flow → wide, axial impellers

### Suction Specific Speed

For cavitation analysis:

```
S = N·√Q / NPSH^(3/4)
```

Typical limits:
- S < 8,500 (US units): Good suction performance
- S > 11,000 (US units): Cavitation risk

## Design Correlations from Literature

### Stepanoff Correlations (1957)

**1. Outlet Diameter:**

```
D₂ = 84.6·√(H/N)  [D₂ in mm, H in m, N in rpm]
```

Or:

```
D₂ = 9.37·√(H/N)  [D₂ in inches, H in ft, N in rpm]
```

**Derivation:**

Based on optimal peripheral velocity:

```
u₂ = K·√(2·g·H)
```

Where K ≈ 1.0-1.1 (typically 1.05)

**2. Outlet Width:**

```
b₂ = K_b·D₂·(Q/(N·D₂³))^(0.65)
```

Where K_b = 2.5-3.5 (smaller for higher specific speeds)

**3. Eye Diameter:**

```
D_eye = 2·√(Q/(π·c_m))
```

Where c_m = 2-4 m/s (axial velocity at eye)

### Gülich Correlations (2014)

**1. Optimal Efficiency:**

```
η_opt = η_size·η_vol·η_Re·η_roughness
```

Where:

**Size effect:**
```
η_size = 1 - 0.03·(D_ref/D₂)^(0.3)
```

D_ref = 400 mm reference diameter

**Volumetric efficiency:**
```
η_vol = 1 - 0.068·(D₂/D_eye)^(-2) / nq^(0.5)
```

**Reynolds number effect:**
```
η_Re = 1 - (Re_ref/Re)^m
```

Where m = 0.1-0.15

**2. Head Coefficient:**

```
Ψ = g·H/(u₂²) ≈ 0.45-0.65  (for backward-curved blades)
```

Increases with decreasing specific speed.

**3. Flow Coefficient:**

```
Φ = Q/(π/4·D₂²·u₂) ≈ 0.05-0.15
```

Increases with increasing specific speed.

### Pfleiderer Formulas

**1. Number of Blades:**

```
Z = 6.5·(D₂ + D₁)/(D₂ - D₁)·sin((β₁ + β₂)/2)
```

Round to nearest integer: typically Z = 5-9

**2. Blade Thickness Effect:**

Effective outlet angle considering blade blockage:

```
β₂_effective = arctan(tan(β₂)·(1 - t·Z/(π·D₂·sin(β₂))))
```

Where t = blade thickness

**3. Slip Factor:**

```
σ = 1 - 0.63·π/Z
```

For backward-curved blades with β₂ ≈ 25°

### Karassik Rules of Thumb

**1. Outlet-to-Inlet Diameter Ratio:**

```
D₂/D₁ = 1.5-2.5  (radial pumps)
D₂/D₁ = 1.2-1.6  (mixed flow pumps)
```

**2. Width Ratio:**

```
b₂/D₂ = 0.012·(1000/nq)^0.65
```

**3. Blade Angles:**

```
β₁ = 15-28° (inlet)
β₂ = 17-35° (outlet, backward-curved)
```

Optimal β₂ ≈ 22-27° for most applications

**4. Number of Blades:**

- Z = 5-6: High specific speed (nq > 60)
- Z = 6-7: Medium specific speed (nq = 30-60)
- Z = 7-9: Low specific speed (nq < 30)

### Pump Affinity Laws

For geometrically similar pumps:

**Speed change (same diameter):**

```
Q₂/Q₁ = N₂/N₁
H₂/H₁ = (N₂/N₁)²
P₂/P₁ = (N₂/N₁)³
```

**Diameter change (same speed):**

```
Q₂/Q₁ = (D₂/D₁)³
H₂/H₁ = (D₂/D₁)²
P₂/P₁ = (D₂/D₁)⁵
```

**Combined:**

```
Q₂/Q₁ = (N₂/N₁)·(D₂/D₁)³
H₂/H₁ = (N₂/N₁)²·(D₂/D₁)²
P₂/P₁ = (N₂/N₁)³·(D₂/D₁)⁵
```

## Loss Analysis

### Major Loss Components

**1. Hydraulic Losses (h_h):**

Friction in impeller and volute:

```
h_h = k_h·u₂²/(2·g)
```

Where k_h = 0.05-0.15 (depends on surface roughness and Reynolds number)

**2. Shock Losses (h_shock):**

At off-design conditions:

```
h_shock = k_shock·(c - c_design)²/(2·g)
```

**3. Leakage Losses:**

Through wear rings and balancing holes:

```
Q_leak = C_d·A_gap·√(2·g·ΔP/ρ)
```

Where:
- C_d = discharge coefficient ≈ 0.7
- A_gap = clearance area
- ΔP = pressure difference

**4. Disk Friction (P_df):**

Power lost to disk friction:

```
P_df = k·ρ·ω³·r₅
```

Where k depends on clearance and surface condition.

**5. Volumetric Losses:**

```
η_vol = Q/(Q + Q_leak)
```

Typically η_vol = 0.96-0.99 for well-designed pumps.

### Efficiency Estimation

**Overall efficiency:**

```
η = η_h·η_vol·η_mech
```

**Component values:**

- Hydraulic: η_h = 0.85-0.94 (depends on Ns, size, finish)
- Volumetric: η_vol = 0.96-0.99 (depends on clearances)
- Mechanical: η_mech = 0.95-0.98 (depends on bearings, seals)

**Empirical correlation (Anderson):**

```
η = 1 - 0.095/nq^0.2  [for nq in European units]
```

Valid for nq = 15-80, medium to large pumps.

## NPSH and Cavitation

### Net Positive Suction Head

**Definition:**

```
NPSH = (P_suction - P_vapor)/(ρ·g) + c_suction²/(2·g)
```

Where:
- P_suction = pressure at pump inlet
- P_vapor = vapor pressure at fluid temperature
- c_suction = velocity at suction

**Required NPSH:**

Empirical (Hydraulic Institute):

```
NPSH_req = (c₁²)/(2·g)·(1 + σ_cavitation)
```

Where σ_cavitation ≈ 1.5-2.5

**Thoma Cavitation Parameter:**

```
σ = NPSH_req/H
```

Typical values:
- σ = 0.05-0.12 (radial, Ns < 2000)
- σ = 0.10-0.20 (mixed flow, Ns = 2000-5000)
- σ = 0.20-0.40 (axial, Ns > 5000)

### Suction Specific Speed

```
S = N·√Q / NPSH_req^(3/4)
```

**Design guidelines:**

- S < 8,500 (US): Excellent suction characteristics
- S = 8,500-11,000: Good, acceptable
- S > 11,000: Poor, cavitation likely

**In European units (nss):**

```
nss = 0.0195·S(US)
```

Limit: nss < 215 (European units)

## Material Selection and Structural Design

### Material Properties

**Common impeller materials:**

| Material | Yield Strength | Density | Applications |
|----------|---------------|---------|--------------|
| Cast Iron | 250 MPa | 7200 kg/m³ | Clean water |
| Bronze | 300 MPa | 8800 kg/m³ | Seawater |
| 316 SS | 290 MPa | 8000 kg/m³ | Chemicals |
| Duplex SS | 450 MPa | 7800 kg/m³ | Corrosive + high stress |

### Rotational Stress

**Centrifugal stress at rim:**

```
σ_r = ρ_mat·ω²·r²/3
```

For solid disk rotating at ω.

**Maximum stress criterion:**

```
σ_max < σ_yield/SF
```

Where SF = safety factor = 3-5

**Critical speed check:**

```
N_critical = (60/2π)·√(k/m_eff)
```

Operating speed should be < 0.7·N_critical or > 1.3·N_critical

## References and Further Reading

### Essential Textbooks

1. **Stepanoff, A.J. (1957)**
   *"Centrifugal and Axial Flow Pumps: Theory, Design, and Application"*
   - 2nd Edition, John Wiley & Sons
   - Classic text, comprehensive design methodology
   - Empirical correlations widely used
   - Chapter 5: Design calculations
   - Chapter 6: Specific speed applications

2. **Gülich, J.F. (2014)**
   *"Centrifugal Pumps"*
   - 3rd Edition, Springer
   - Most comprehensive modern reference
   - 1200+ pages, detailed theory and practice
   - Chapter 3: Pump hydraulics and velocity triangles
   - Chapter 8: Design methodology
   - Chapter 10: Efficiency and losses
   - Extensive correlations and data

3. **Karassik, I.J., et al. (2008)**
   *"Pump Handbook"*
   - 4th Edition, McGraw-Hill
   - Industry standard reference
   - Section 2.1: Centrifugal pump theory
   - Section 2.2: Design and construction
   - Section 2.3: Performance characteristics

4. **Lobanoff, V.S. and Ross, R.R. (1992)**
   *"Centrifugal Pumps: Design and Application"*
   - 2nd Edition, Gulf Publishing
   - Practical design focus
   - Chapter 3: Hydraulic design
   - Chapter 4: Mechanical design

5. **Brennen, C.E. (1994)**
   *"Hydrodynamics of Pumps"*
   - Oxford University Press
   - Advanced fluid dynamics
   - Detailed cavitation analysis
   - CFD validation

### Technical Standards

1. **API 610 (2010)**
   *"Centrifugal Pumps for Petroleum, Petrochemical and Natural Gas Industries"*
   - Design standards
   - Material specifications
   - Testing requirements

2. **ISO 9906 (2012)**
   *"Rotodynamic Pumps - Hydraulic Performance Acceptance Tests"*
   - Performance testing procedures
   - Tolerance grades

3. **Hydraulic Institute Standards (2014)**
   *"ANSI/HI 1.3 - Rotodynamic Centrifugal Pumps for Design and Application"*
   - US industry standards
   - NPSH requirements
   - Operating envelope

### Key Journal Papers

1. **Wiesner, F.J. (1967)**
   "A Review of Slip Factors for Centrifugal Impellers"
   *ASME Journal of Engineering for Power*, Vol. 89, pp. 558-572
   - Comprehensive slip factor analysis
   - Various correlations compared

2. **Pfleiderer, C. (1961)**
   "Die Kreiselpumpen für Flüssigkeiten und Gase"
   *Springer-Verlag*
   - Classic German text
   - Blade number formula
   - Secondary flow analysis

3. **Anderson, H.H. (1980)**
   *"Centrifugal Pumps"*
   - Trade & Technical Press
   - Practical design procedures
   - Efficiency correlations

4. **Lazarkiewicz, S. and Troskolanski, A.T. (1965)**
   *"Impeller Pumps"*
   - Pergamon Press
   - Detailed velocity triangle analysis
   - Loss coefficient data

### Online Resources

1. **Pump Fundamentals**
   - PumpScout.com/articles
   - Engineering ToolBox
   - Comprehensive reference data

2. **Goulds Pumps Application Engineering Manual**
   - Free download from ITT
   - Practical design examples
   - System design guidance

3. **Sulzer Centrifugal Pump Handbook**
   - Comprehensive design guide
   - Selection charts
   - Troubleshooting

### CFD Validation References

1. **Gonzalez, J. et al. (2002)**
   "Numerical Simulation of the Dynamic Effects Due to Impeller-Volute Interaction"
   *ASME Journal of Fluids Engineering*, Vol. 124

2. **Majidi, K. (2005)**
   "Numerical Study of Unsteady Flow in a Centrifugal Pump"
   *ASME Journal of Turbomachinery*, Vol. 127

3. **Barrio, R. et al. (2010)**
   "Performance Characteristics and Internal Flow Patterns in a Reverse-Running Pump-Turbine"
   *Proceedings of IMechE Part C*

### Historical References

1. **Euler, L. (1754)**
   "Théorie plus complète des machines qui sont mises en mouvement par la réaction de l'eau"
   - Original turbine equation
   - Fundamental energy transfer

2. **Osborne Reynolds (1875)**
   "On the Efficiency of Centrifugal Pumps"
   - Early efficiency analysis
   - Slip phenomenon

## Nomenclature

### Geometric Parameters

- D₁, D₂ = impeller inlet/outlet diameter (m)
- D_eye = eye (suction) diameter (m)
- b₁, b₂ = impeller inlet/outlet width (m)
- β₁, β₂ = blade inlet/outlet angle (degrees)
- Z = number of blades
- t = blade thickness (m)

### Velocities

- u = blade peripheral velocity (m/s)
- c = absolute velocity (m/s)
- w = relative velocity (m/s)
- c_m = meridional (through-flow) velocity component (m/s)
- c_u = tangential (whirl) velocity component (m/s)
- α = absolute flow angle (degrees)
- β = relative flow angle (degrees)

### Performance Parameters

- Q = flow rate (m³/s)
- H = head (m)
- N = rotational speed (rpm)
- ω = angular velocity (rad/s)
- P = power (W or kW)
- T = torque (N·m)
- η = efficiency (dimensionless)
- nq, Ns = specific speed (various conventions)
- S = suction specific speed

### Fluid Properties

- ρ = density (kg/m³)
- μ = dynamic viscosity (Pa·s)
- ν = kinematic viscosity (m²/s)
- g = gravitational acceleration (9.81 m/s²)

### Dimensionless Groups

- Re = Reynolds number
- Φ = flow coefficient
- Ψ = head coefficient
- σ = cavitation parameter or slip factor

## Appendix: Unit Conversions

### Flow Rate

- 1 m³/s = 3600 m³/h
- 1 m³/s = 15850 gpm (US gallons per minute)
- 1 m³/h = 4.403 gpm

### Head

- 1 m = 3.281 ft
- 1 m = 0.0981 bar
- 1 m = 1.422 psi

### Power

- 1 kW = 1.341 hp (horsepower)
- 1 hp = 0.746 kW

### Specific Speed Conversion

```
Ns(US) = 51.6·nq(European)
```

Where:
- Ns(US): rpm, gpm, ft
- nq(European): rpm, m³/s, m

### Angular Velocity

```
ω (rad/s) = 2π·N/60
```

Where N is in rpm.
