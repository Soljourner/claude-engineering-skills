# Fluid Dynamics Equations Reference

A comprehensive reference for fundamental and advanced equations in fluid dynamics, computational fluid dynamics, and turbulence modeling.

## Table of Contents

1. [Fundamental Equations](#fundamental-equations)
2. [Navier-Stokes Equations](#navier-stokes-equations)
3. [Turbulence Models](#turbulence-models)
4. [Dimensionless Numbers](#dimensionless-numbers)
5. [Boundary Layer Equations](#boundary-layer-equations)
6. [Special Cases](#special-cases)

---

## Fundamental Equations

### Conservation of Mass (Continuity)

**General form (compressible)**:
```
∂ρ/∂t + ∇·(ρV) = 0

Expanded:
∂ρ/∂t + ∂(ρu)/∂x + ∂(ρv)/∂y + ∂(ρw)/∂z = 0
```

**Incompressible flow** (ρ = constant):
```
∇·V = 0

Expanded:
∂u/∂x + ∂v/∂y + ∂w/∂z = 0
```

**Cylindrical coordinates** (r, θ, z):
```
∂u_r/∂r + u_r/r + (1/r)∂u_θ/∂θ + ∂u_z/∂z = 0
```

**Spherical coordinates** (r, θ, φ):
```
(1/r²)∂(r²u_r)/∂r + (1/(r sin θ))∂(u_θ sin θ)/∂θ + (1/(r sin θ))∂u_φ/∂φ = 0
```

---

## Navier-Stokes Equations

### Incompressible Navier-Stokes

**Vector form**:
```
ρ(∂V/∂t + V·∇V) = -∇p + μ∇²V + ρg

Where:
- V·∇V = u∂V/∂x + v∂V/∂y + w∂V/∂z (convective acceleration)
- ∇²V = ∂²V/∂x² + ∂²V/∂y² + ∂²V/∂z² (viscous term)
```

**Cartesian components**:

**X-momentum**:
```
ρ(∂u/∂t + u∂u/∂x + v∂u/∂y + w∂u/∂z) = -∂p/∂x + μ(∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²) + ρg_x
```

**Y-momentum**:
```
ρ(∂v/∂t + u∂v/∂x + v∂v/∂y + w∂v/∂z) = -∂p/∂y + μ(∂²v/∂x² + ∂²v/∂y² + ∂²v/∂z²) + ρg_y
```

**Z-momentum**:
```
ρ(∂w/∂t + u∂w/∂x + v∂w/∂y + w∂w/∂z) = -∂p/∂z + μ(∂²w/∂x² + ∂²w/∂y² + ∂²w/∂z²) + ρg_z
```

### Non-dimensional Form

**Using characteristic scales** (L, V₀):
```
Re(∂u*/∂t* + u*·∇*u*) = -∇*p* + ∇*²u* + (1/Fr²)g*

Where:
- Re = ρV₀L/μ (Reynolds number)
- Fr = V₀/√(gL) (Froude number)
- * denotes non-dimensional quantities
```

### Compressible Navier-Stokes

**Momentum (with stress tensor)**:
```
ρ(∂V/∂t + V·∇V) = -∇p + ∇·τ + ρg

Stress tensor:
τᵢⱼ = μ[(∂u_i/∂x_j + ∂u_j/∂x_i) - (2/3)δᵢⱼ(∇·V)] + λ(∇·V)δᵢⱼ

Stokes hypothesis:
λ = -2μ/3
```

**Total energy equation**:
```
∂(ρE)/∂t + ∇·[(ρE + p)V] = ∇·(k∇T) + ∇·(τ·V) + ρV·g

Where:
E = e + V²/2 (total energy = internal + kinetic)
```

**Alternative form (temperature)**:
```
ρc_p(∂T/∂t + V·∇T) = ∂p/∂t + V·∇p + ∇·(k∇T) + Φ

Viscous dissipation:
Φ = τᵢⱼ(∂u_i/∂x_j)
```

### Euler Equations (Inviscid)

**Conservation form**:
```
∂U/∂t + ∂F/∂x + ∂G/∂y + ∂H/∂z = 0

Where:
U = [ρ, ρu, ρv, ρw, ρE]ᵀ

F = [ρu, ρu² + p, ρuv, ρuw, u(ρE + p)]ᵀ
```

---

## Turbulence Models

### Reynolds-Averaged Navier-Stokes (RANS)

**Reynolds decomposition**:
```
u = ū + u'
p = p̄ + p'

Where:
- ¯ denotes time average
- ' denotes fluctuation
```

**RANS equations**:
```
ρ(∂ūᵢ/∂t + ūⱼ∂ūᵢ/∂xⱼ) = -∂p̄/∂xᵢ + ∂/∂xⱼ[μ(∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ) - ρu'ᵢu'ⱼ]

Reynolds stress tensor:
τᵢⱼᴿ = -ρu'ᵢu'ⱼ
```

### k-ε Models

#### Standard k-ε

**Turbulent kinetic energy (k)**:
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₖ)∂k/∂xⱼ] + Pₖ - ρε

Where:
k = (1/2)u'ᵢu'ᵢ (turbulent kinetic energy)
```

**Dissipation rate (ε)**:
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₑ)∂ε/∂xⱼ] + C₁ₑ(ε/k)Pₖ - C₂ₑρ(ε²/k)

Where:
ε = ν(∂u'ᵢ/∂xⱼ)(∂u'ᵢ/∂xⱼ) (dissipation rate)
```

**Production term**:
```
Pₖ = μₜ(∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ)(∂ūᵢ/∂xⱼ) - (2/3)ρk(∂ūₖ/∂xₖ)

Simplified:
Pₖ = μₜS² (S = strain rate magnitude)
```

**Turbulent viscosity**:
```
μₜ = ρC_μk²/ε

Effective viscosity:
μₑff = μ + μₜ
```

**Standard constants**:
```
C_μ = 0.09
C₁ₑ = 1.44
C₂ₑ = 1.92
σₖ = 1.0
σₑ = 1.3
```

#### Realizable k-ε

**Modified dissipation equation**:
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₑ)∂ε/∂xⱼ] + ρC₁Sε - ρC₂(ε²/(k + √(νε)))

Where:
C₁ = max[0.43, η/(η + 5)]
η = S(k/ε)
S = √(2SᵢⱼSᵢⱼ)
```

**Variable C_μ**:
```
C_μ = 1/(A₀ + Aₛ(kU*/ε))

U* = √(SᵢⱼSᵢⱼ + Ω̃ᵢⱼΩ̃ᵢⱼ)
Ω̃ᵢⱼ = Ωᵢⱼ - 2εᵢⱼₖωₖ (rotation rate tensor)
```

**Constants**:
```
C₂ = 1.9
σₖ = 1.0
σₑ = 1.2
```

#### RNG k-ε

**Additional term in ε equation**:
```
Rₑ = -C_μρη³(1 - η/η₀)/(1 + βη³) × ε²/k

Where:
η = Sk/ε
η₀ = 4.38
β = 0.012
```

### k-ω Models

#### Standard k-ω (Wilcox)

**Turbulent kinetic energy (k)**:
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσₖ)∂k/∂xⱼ] + Pₖ - β*ρkω
```

**Specific dissipation rate (ω)**:
```
∂(ρω)/∂t + ∂(ρωūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσ_ω)∂ω/∂xⱼ] + α(ω/k)Pₖ - βρω²

Where:
ω = ε/(C_μk) (specific dissipation rate, 1/s)
```

**Turbulent viscosity**:
```
μₜ = ρk/ω
```

**Constants**:
```
α = 5/9
β = 3/40
β* = 9/100
σₖ = 2
σ_ω = 2
```

#### k-ω SST (Menter)

**Blending between k-ω (near wall) and k-ε (far field)**:

**Blending function**:
```
F₁ = tanh(arg₁⁴)

arg₁ = min[max(√k/(β*ωy), 500ν/(y²ω)), 4ρσ_ω₂k/(CDₖ_ω y²)]

Where:
CDₖ_ω = max(2ρσ_ω₂(1/ω)∂k/∂xⱼ∂ω/∂xⱼ, 10⁻¹⁰)
```

**k equation**:
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσₖ)∂k/∂xⱼ] + P̃ₖ - β*ρkω

P̃ₖ = min(Pₖ, 10β*ρkω) (production limiter)
```

**ω equation**:
```
∂(ρω)/∂t + ∂(ρωūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσ_ω)∂ω/∂xⱼ] + α(ω/k)Pₖ - βρω² + 2(1-F₁)ρσ_ω₂(1/ω)∂k/∂xⱼ∂ω/∂xⱼ
```

**Turbulent viscosity with strain rate limiter**:
```
μₜ = ρa₁k/max(a₁ω, SF₂)

F₂ = tanh(arg₂²)
arg₂ = max(2√k/(β*ωy), 500ν/(y²ω))
S = √(2SᵢⱼSᵢⱼ)
a₁ = 0.31
```

**Constants (blended)**:
```
φ = F₁φ₁ + (1 - F₁)φ₂

Inner constants (k-ω):
α₁ = 5/9, β₁ = 0.075, σₖ₁ = 2, σ_ω₁ = 2

Outer constants (k-ε transformed):
α₂ = 0.44, β₂ = 0.0828, σₖ₂ = 1, σ_ω₂ = 1/0.856

β* = 0.09
```

### Spalart-Allmaras (One-Equation Model)

**Transport equation for ν̃**:
```
∂ν̃/∂t + ūⱼ∂ν̃/∂xⱼ = C_b₁S̃ν̃ - C_w₁f_w(ν̃/d)² + (1/σ)[∂/∂xⱼ((ν + ν̃)∂ν̃/∂xⱼ) + C_b₂(∂ν̃/∂xⱼ)²]

Where:
ν̃ = modified turbulent viscosity
d = distance to nearest wall
```

**Turbulent viscosity**:
```
νₜ = ν̃f_v₁

f_v₁ = χ³/(χ³ + C³_v₁)
χ = ν̃/ν
```

**Modified vorticity**:
```
S̃ = S + (ν̃/κ²d²)f_v₂

f_v₂ = 1 - χ/(1 + χf_v₁)
S = |Ωᵢⱼ| = vorticity magnitude
```

**Destruction function**:
```
f_w = g[(1 + C⁶_w₃)/(g⁶ + C⁶_w₃)]^(1/6)

g = r + C_w₂(r⁶ - r)
r = min(ν̃/(S̃κ²d²), 10)
```

**Constants**:
```
σ = 2/3
C_b₁ = 0.1355
C_b₂ = 0.622
C_v₁ = 7.1
C_w₁ = C_b₁/κ² + (1 + C_b₂)/σ
C_w₂ = 0.3
C_w₃ = 2
κ = 0.41 (von Kármán constant)
```

### Reynolds Stress Model (RSM)

**Transport equations for Reynolds stresses**:
```
∂(ρu'ᵢu'ⱼ)/∂t + ūₖ∂(ρu'ᵢu'ⱼ)/∂xₖ = Pᵢⱼ + Φᵢⱼ - εᵢⱼ + Dᵢⱼ

Production:
Pᵢⱼ = -ρ(u'ᵢu'ₖ∂ūⱼ/∂xₖ + u'ⱼu'ₖ∂ūᵢ/∂xₖ)

Pressure-strain:
Φᵢⱼ = p'(∂u'ᵢ/∂xⱼ + ∂u'ⱼ/∂xᵢ)

Dissipation:
εᵢⱼ = 2μ(∂u'ᵢ/∂xₖ)(∂u'ⱼ/∂xₖ)

Diffusion:
Dᵢⱼ = ∂/∂xₖ[μ∂(u'ᵢu'ⱼ)/∂xₖ + ρu'ᵢu'ⱼu'ₖ + p'(δₖⱼu'ᵢ + δᵢₖu'ⱼ)]
```

**ε equation** (same as k-ε):
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₑ)∂ε/∂xⱼ] + C₁ₑ(ε/k)Pₖ - C₂ₑρ(ε²/k)
```

**Turbulent viscosity**:
```
μₜ = ρC_μk²/ε
k = (1/2)u'ᵢu'ᵢ
```

### Large Eddy Simulation (LES)

**Filtered Navier-Stokes**:
```
∂ū̄ᵢ/∂t + ∂(ū̄ᵢū̄ⱼ)/∂xⱼ = -(1/ρ)∂p̄̄/∂xᵢ + ν∂²ū̄ᵢ/∂xⱼ² - ∂τᵢⱼˢᵍˢ/∂xⱼ

Subgrid-scale stress:
τᵢⱼˢᵍˢ = uᵢuⱼ - ū̄ᵢū̄ⱼ
```

#### Smagorinsky Model

**Subgrid viscosity**:
```
νₜ = (Cₛ Δ)²|S̄̄|

|S̄̄| = √(2S̄̄ᵢⱼS̄̄ᵢⱼ)
Δ = (ΔxΔyΔz)^(1/3) (filter width)
Cₛ = 0.1 - 0.2 (Smagorinsky constant)
```

**Subgrid stress**:
```
τᵢⱼˢᵍˢ - (1/3)τₖₖˢᵍˢδᵢⱼ = -2νₜS̄̄ᵢⱼ
```

#### WALE Model

**Better near-wall behavior**:
```
νₜ = (C_w Δ)² [(S̄̄ᵈᵢⱼS̄̄ᵈᵢⱼ)^(3/2)] / [(S̄̄ᵢⱼS̄̄ᵢⱼ)^(5/2) + (S̄̄ᵈᵢⱼS̄̄ᵈᵢⱼ)^(5/4)]

Traceless symmetric part:
S̄̄ᵈᵢⱼ = (1/2)(ḡ²ᵢⱼ + ḡ²ⱼᵢ) - (1/3)δᵢⱼḡ²ₖₖ
ḡᵢⱼ = ∂ū̄ᵢ/∂xⱼ

C_w = 0.325
```

---

## Dimensionless Numbers

### Reynolds Number

**Definition**:
```
Re = ρVL/μ = VL/ν = (inertial forces)/(viscous forces)

Where:
- ρ = density
- V = characteristic velocity
- L = characteristic length
- μ = dynamic viscosity
- ν = kinematic viscosity
```

**Critical values**:
```
Pipes:         Re_crit ≈ 2300
Flat plate:    Re_x,crit ≈ 5×10⁵
Sphere:        Re_crit ≈ 2×10⁵ (drag crisis)
```

### Mach Number

**Definition**:
```
Ma = V/c = V/√(γRT)

Where:
- c = speed of sound
- γ = specific heat ratio (1.4 for air)
- R = specific gas constant
- T = temperature
```

**Flow regimes**:
```
Ma < 0.3:      Incompressible
0.3 < Ma < 0.8: Subsonic compressible
0.8 < Ma < 1.2: Transonic
1.2 < Ma < 5:   Supersonic
Ma > 5:         Hypersonic
```

### Froude Number

**Definition**:
```
Fr = V/√(gL) = (inertial forces)/(gravitational forces)

For open channel:
Fr = V/√(gh) (h = depth)
```

**Flow regimes**:
```
Fr < 1:  Subcritical (controlled by downstream)
Fr = 1:  Critical flow
Fr > 1:  Supercritical (controlled by upstream)
```

### Strouhal Number

**Definition**:
```
St = fL/V = (oscillation frequency)/(flow frequency)

Where:
- f = oscillation frequency
- L = characteristic length
- V = velocity
```

**Typical values**:
```
Cylinder vortex shedding:  St ≈ 0.2 (Re > 1000)
Bluff body wakes:          St ≈ 0.1 - 0.3
```

### Prandtl Number

**Definition**:
```
Pr = ν/α = μc_p/k = (momentum diffusivity)/(thermal diffusivity)

Where:
- α = k/(ρc_p) = thermal diffusivity
- c_p = specific heat
- k = thermal conductivity
```

**Typical values**:
```
Liquid metals:  Pr ≈ 0.01
Gases:         Pr ≈ 0.7
Water:         Pr ≈ 7
Oils:          Pr ≈ 100 - 1000
```

### Schmidt Number

**Definition**:
```
Sc = ν/D = μ/(ρD) = (momentum diffusivity)/(mass diffusivity)

Where:
- D = mass diffusion coefficient
```

### Grashof Number

**Definition**:
```
Gr = gβΔTL³/ν² = (buoyancy forces)/(viscous forces)

Where:
- β = thermal expansion coefficient
- ΔT = temperature difference
```

**Natural convection regime**:
```
Ra = Gr × Pr (Rayleigh number)

Ra < 10⁸:  Laminar natural convection
Ra > 10⁹:  Turbulent natural convection
```

### Nusselt Number

**Definition**:
```
Nu = hL/k = (convective heat transfer)/(conductive heat transfer)

Where:
- h = convective heat transfer coefficient
- L = characteristic length
- k = thermal conductivity
```

**Correlations**:
```
Flat plate (laminar):    Nu_x = 0.332Re_x^(1/2)Pr^(1/3)
Flat plate (turbulent):  Nu_x = 0.0296Re_x^(4/5)Pr^(1/3)
Cylinder:                Nu = CRe^m Pr^(1/3)
```

### Weber Number

**Definition**:
```
We = ρV²L/σ = (inertial forces)/(surface tension forces)

Where:
- σ = surface tension
```

### Capillary Number

**Definition**:
```
Ca = μV/σ = (viscous forces)/(surface tension forces)
```

### Eckert Number

**Definition**:
```
Ec = V²/(c_p ΔT) = (kinetic energy)/(enthalpy difference)
```

---

## Boundary Layer Equations

### Blasius Solution (Flat Plate, Laminar)

**Similarity variable**:
```
η = y√(U_∞/(νx))

Stream function:
ψ = √(U_∞νx) f(η)
```

**Blasius equation**:
```
f''' + (1/2)ff'' = 0

Boundary conditions:
f(0) = 0, f'(0) = 0, f'(∞) = 1
```

**Velocity profile**:
```
u/U_∞ = f'(η)
v/U_∞ = (1/2)√(ν/(U_∞x))[ηf'(η) - f(η)]
```

**Boundary layer thickness**:
```
δ = 5.0x/√Re_x
δ* = 1.721x/√Re_x (displacement thickness)
θ = 0.664x/√Re_x (momentum thickness)
```

**Skin friction**:
```
C_f = τ_w/(1/2 ρU²_∞) = 0.664/√Re_x

Average (0 to x):
C̄_f = 1.328/√Re_x
```

### Prandtl Boundary Layer Equations

**Assumptions**:
- δ << L
- u >> v
- ∂/∂x << ∂/∂y

**Simplified equations**:
```
∂u/∂x + ∂v/∂y = 0

u∂u/∂x + v∂u/∂y = U_∞dU_∞/dx + ν∂²u/∂y²

Thermal:
u∂T/∂x + v∂T/∂y = α∂²T/∂y²
```

### von Kármán Momentum Integral

**Integral form**:
```
dθ/dx + (2 + H)θ/U_∞ dU_∞/dx = C_f/2

Where:
θ = ∫₀^∞ (u/U_∞)(1 - u/U_∞) dy (momentum thickness)
δ* = ∫₀^∞ (1 - u/U_∞) dy (displacement thickness)
H = δ*/θ (shape factor)
```

**Shape factor**:
```
Laminar:   H ≈ 2.6
Turbulent: H ≈ 1.3 - 1.4
Separated: H > 2.4
```

### Turbulent Boundary Layer

**Power law velocity profile**:
```
u/U_∞ = (y/δ)^(1/n)

n ≈ 7 (typical)
```

**Logarithmic law of the wall**:
```
u⁺ = (1/κ)ln(y⁺) + B

u⁺ = u/u_τ (velocity in wall units)
y⁺ = ρu_τy/μ (distance in wall units)
u_τ = √(τ_w/ρ) (friction velocity)
κ = 0.41 (von Kármán constant)
B = 5.0 (constant)
```

**Three-layer structure**:
```
Viscous sublayer (y⁺ < 5):
u⁺ = y⁺

Buffer layer (5 < y⁺ < 30):
Transition region

Log layer (30 < y⁺ < 0.2Re_θ):
u⁺ = (1/κ)ln(y⁺) + B
```

**Skin friction (Prandtl-Schlichting)**:
```
C_f = 0.0592Re_x^(-1/5)  (5×10⁵ < Re_x < 10⁷)

C_f = 0.370(log₁₀Re_x)^(-2.584)  (Re_x > 10⁷)
```

---

## Special Cases

### Stokes Flow (Creeping Flow)

**Re << 1, neglect inertia**:
```
0 = -∇p + μ∇²V + ρg
∇·V = 0
```

**Drag on sphere** (Stokes drag):
```
F_D = 6πμRV

C_D = 24/Re
```

### Potential Flow (Inviscid, Irrotational)

**Velocity potential**:
```
V = ∇φ
∇²φ = 0 (Laplace equation)
```

**Stream function** (2D, incompressible):
```
u = ∂ψ/∂y
v = -∂ψ/∂x
∇²ψ = -ω_z (vorticity)
```

**Complex potential** (2D):
```
W(z) = φ + iψ
dW/dz = u - iv
```

### Bernoulli Equation

**Along streamline** (steady, inviscid, incompressible):
```
p + (1/2)ρV² + ρgh = constant

Pressure coefficient:
C_p = (p - p_∞)/(1/2 ρV²_∞)
```

**Compressible** (isentropic):
```
(γ/(γ-1))(p/ρ) + V²/2 = constant

Total pressure:
p₀/p = (1 + ((γ-1)/2)Ma²)^(γ/(γ-1))
```

### Rankine-Hugoniot Relations (Shock Waves)

**Normal shock**:
```
ρ₂/ρ₁ = ((γ+1)Ma₁²) / (2 + (γ-1)Ma₁²)

p₂/p₁ = 1 + (2γ/(γ+1))(Ma₁² - 1)

Ma₂² = ((γ-1)Ma₁² + 2) / (2γMa₁² - (γ-1))
```

### Prandtl-Meyer Expansion

**Expansion fan** (supersonic):
```
ν(Ma) = √((γ+1)/(γ-1)) arctan(√(((γ-1)/(γ+1))(Ma² - 1))) - arctan(√(Ma² - 1))

Flow turning:
θ = ν(Ma₂) - ν(Ma₁)
```

---

## Summary Tables

### Turbulence Model Selection Guide

| Flow Type | Recommended Model | Alternative |
|-----------|------------------|-------------|
| External aerodynamics | k-ω SST | Spalart-Allmaras |
| Internal flows | k-ε Realizable | k-ω SST |
| Separated flows | k-ω SST | LES |
| Rotating flows | RSM | k-ω SST |
| Natural convection | k-ε RNG | k-ω SST |
| Transitional | γ-Re_θ SST | k-ω SST |
| Highly unsteady | LES | URANS |

### Typical Dimensionless Numbers

| Number | Symbol | Typical Range |
|--------|--------|---------------|
| Reynolds | Re | 10⁰ - 10⁸ |
| Mach | Ma | 0 - 5+ |
| Froude | Fr | 0.1 - 10 |
| Prandtl (air) | Pr | ~0.7 |
| Prandtl (water) | Pr | ~7 |
| Grashof | Gr | 10³ - 10¹² |
| Nusselt | Nu | 1 - 1000+ |

---

## References

1. Pope, S.B. "Turbulent Flows" (2000)
2. Wilcox, D.C. "Turbulence Modeling for CFD" (2006)
3. White, F.M. "Viscous Fluid Flow" (2006)
4. Anderson, J.D. "Computational Fluid Dynamics" (1995)
5. Schlichting, H. "Boundary Layer Theory" (2000)
