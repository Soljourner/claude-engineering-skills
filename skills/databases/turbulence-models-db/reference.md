# Turbulence Models Technical Reference

Detailed mathematical formulations, constants, validation test cases, and academic references for turbulence models.

## Table of Contents

1. [Fundamental Equations](#fundamental-equations)
2. [RANS Model Equations](#rans-model-equations)
3. [Model Constants](#model-constants)
4. [Validation Test Cases](#validation-test-cases)
5. [Academic References](#academic-references)

---

## Fundamental Equations

### Reynolds-Averaged Navier-Stokes (RANS)

The time-averaged continuity and momentum equations:

**Continuity:**
```
∂ρ/∂t + ∂(ρūᵢ)/∂xᵢ = 0
```

**Momentum:**
```
∂(ρūᵢ)/∂t + ∂(ρūᵢūⱼ)/∂xⱼ = -∂p̄/∂xᵢ + ∂/∂xⱼ[μ(∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ) - ρu'ᵢu'ⱼ] + ρgᵢ
```

Where:
- ūᵢ = mean velocity component
- u'ᵢ = fluctuating velocity component
- ρu'ᵢu'ⱼ = Reynolds stress tensor (requires modeling)

### Reynolds Stress Tensor

The Reynolds stress represents the effect of turbulent fluctuations on the mean flow:

```
τᵢⱼᴿᴬᴺˢ = -ρu'ᵢu'ⱼ
```

### Boussinesq Hypothesis

Most two-equation models use the Boussinesq approximation:

```
-ρu'ᵢu'ⱼ = μₜ(∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ) - (2/3)ρkδᵢⱼ
```

Where:
- μₜ = turbulent (eddy) viscosity
- k = turbulent kinetic energy = (1/2)u'ᵢu'ᵢ
- δᵢⱼ = Kronecker delta

### Turbulent Kinetic Energy

```
k = (1/2)(u'² + v'² + w'²)
```

### Turbulent Dissipation Rate

```
ε = ν ∂u'ᵢ/∂xⱼ ∂u'ᵢ/∂xⱼ
```

Where ν is the kinematic viscosity.

### Specific Dissipation Rate

```
ω = ε/(β*k)
```

Where β* = 0.09 for k-ω models.

---

## RANS Model Equations

### Standard k-ε Model

#### Turbulent Viscosity
```
μₜ = ρCμk²/ε
```

#### Turbulent Kinetic Energy Equation
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₖ)∂k/∂xⱼ] + Pₖ - ρε
```

#### Dissipation Rate Equation
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σε)∂ε/∂xⱼ] + C₁ε(ε/k)Pₖ - C₂ερε(ε/k)
```

#### Production Term
```
Pₖ = μₜS²
```

Where:
```
S = √(2SᵢⱼSᵢⱼ)
Sᵢⱼ = (1/2)(∂ūᵢ/∂xⱼ + ∂ūⱼ/∂xᵢ)
```

#### Model Constants
- Cμ = 0.09
- C₁ε = 1.44
- C₂ε = 1.92
- σₖ = 1.0
- σε = 1.3

---

### RNG k-ε Model

#### Turbulent Viscosity
```
μₜ = ρCμk²/ε
```

Where Cμ = 0.0845 (different from standard k-ε)

#### k Equation
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₖ)∂k/∂xⱼ] + Pₖ - ρε
```

#### ε Equation
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σε)∂ε/∂xⱼ] + C₁ε(ε/k)Pₖ - C₂ε*ρ(ε²/k)
```

#### Modified C₂ε
```
C₂ε* = C₂ε + Cμη³(1 - η/η₀)/(1 + βη³)
```

Where:
```
η = Sk/ε
S = √(2SᵢⱼSᵢⱼ)
```

#### Model Constants
- Cμ = 0.0845
- C₁ε = 1.42
- C₂ε = 1.68
- σₖ = 0.7179
- σε = 0.7179
- η₀ = 4.38
- β = 0.012

---

### Realizable k-ε Model

#### Turbulent Viscosity
```
μₜ = ρCμk²/ε
```

Where Cμ is now a function of the flow:

```
Cμ = 1/(A₀ + Aₛ(kU*/ε))
```

Where:
```
U* = √(SᵢⱼSᵢⱼ + Ω̃ᵢⱼΩ̃ᵢⱼ)
Ω̃ᵢⱼ = Ωᵢⱼ - 2εᵢⱼₖωₖ (in rotating frame)
Ωᵢⱼ = (1/2)(∂ūᵢ/∂xⱼ - ∂ūⱼ/∂xᵢ)
```

Constants:
- A₀ = 4.04
- Aₛ = √6 cos(φ)
- φ = (1/3)cos⁻¹(√6W)
- W = (SᵢⱼSⱼₖSₖᵢ)/S³
- S = √(2SᵢⱼSᵢⱼ)

#### k Equation
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σₖ)∂k/∂xⱼ] + Pₖ - ρε
```

#### ε Equation
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σε)∂ε/∂xⱼ] + ρC₁Sε - ρC₂(ε²/(k + √(νε)))
```

Where:
```
C₁ = max[0.43, η/(η + 5)]
η = Sk/ε
S = √(2SᵢⱼSᵢⱼ)
```

#### Model Constants
- C₁ε = 1.44
- C₂ = 1.9
- σₖ = 1.0
- σε = 1.2

---

### Standard k-ω Model (Wilcox 1988)

#### Turbulent Viscosity
```
μₜ = ρk/ω
```

#### k Equation
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσₖ)∂k/∂xⱼ] + Pₖ - β*ρkω
```

#### ω Equation
```
∂(ρω)/∂t + ∂(ρωūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσω)∂ω/∂xⱼ] + α(ω/k)Pₖ - βρω²
```

#### Production Term
```
Pₖ = μₜS²
S = √(2SᵢⱼSᵢⱼ)
```

#### Model Constants (1988 version)
- α = 5/9
- β = 0.075
- β* = 0.09
- σₖ = 2.0
- σω = 2.0

#### Model Constants (2006 version)
- α = 13/25
- β = β₀fᵦ
- β₀ = 0.0708
- β* = 0.09
- σₖ = 0.6
- σω = 0.5

---

### k-ω SST Model (Menter 1994)

The SST model blends the k-ω model near walls with the k-ε model in the freestream using a blending function.

#### Turbulent Viscosity
```
μₜ = ρa₁k/max(a₁ω, SF₂)
```

Where:
- a₁ = 0.31
- S = √(2SᵢⱼSᵢⱼ)
- F₂ = blending function

#### k Equation
```
∂(ρk)/∂t + ∂(ρkūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσₖ)∂k/∂xⱼ] + P̃ₖ - β*ρkω
```

Where:
```
P̃ₖ = min(Pₖ, 10β*ρkω)  (production limiter)
Pₖ = μₜS²
```

#### ω Equation
```
∂(ρω)/∂t + ∂(ρωūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜσω)∂ω/∂xⱼ] + α(ω/k)Pₖ - βρω² + 2(1-F₁)ρσω₂(1/ω)(∂k/∂xⱼ)(∂ω/∂xⱼ)
```

#### Blending Functions

**F₁:**
```
F₁ = tanh(arg₁⁴)
```

Where:
```
arg₁ = min[max(√k/(β*ωy), 500ν/(y²ω)), 4ρσω₂k/(CDₖωy²)]
CDₖω = max(2ρσω₂(1/ω)(∂k/∂xⱼ)(∂ω/∂xⱼ), 10⁻¹⁰)
```

**F₂:**
```
F₂ = tanh(arg₂²)
```

Where:
```
arg₂ = max(2√k/(β*ωy), 500ν/(y²ω))
```

In the above, y is the distance to the nearest wall.

#### Blended Constants
Any constant φ is computed as:
```
φ = F₁φ₁ + (1 - F₁)φ₂
```

#### Model Constants

**Set 1 (k-ω inner):**
- α₁ = 5/9
- β₁ = 0.075
- β* = 0.09
- σₖ₁ = 2.0
- σω₁ = 2.0

**Set 2 (k-ε outer, transformed):**
- α₂ = 0.44
- β₂ = 0.0828
- β* = 0.09
- σₖ₂ = 1.0
- σω₂ = 1.168

**Additional:**
- a₁ = 0.31 (SST limiter)
- κ = 0.41 (von Karman constant)

---

### Spalart-Allmaras Model

One-equation model for modified turbulent viscosity ν̃.

#### Transport Equation
```
∂(ρν̃)/∂t + ∂(ρν̃ūᵢ)/∂xᵢ = Cᵦ₁ρS̃ν̃ + (1/σ)∂/∂xⱼ[(μ + ρν̃)∂ν̃/∂xⱼ] + Cᵦ₂(ρ/σ)(∂ν̃/∂xⱼ)² - [Cw₁fw - (Cᵦ₁/κ²)fₜ₂]ρ(ν̃/d)²
```

#### Turbulent Viscosity
```
μₜ = ρν̃fᵥ₁
```

Where:
```
fᵥ₁ = χ³/(χ³ + Cᵥ₁³)
χ = ν̃/ν
```

#### Modified Vorticity
```
S̃ = S + (ν̃/κ²d²)fᵥ₂
```

Where:
```
S = √(2ΩᵢⱼΩᵢⱼ)
Ωᵢⱼ = (1/2)(∂ūᵢ/∂xⱼ - ∂ūⱼ/∂xᵢ)
fᵥ₂ = 1 - χ/(1 + χfᵥ₁)
```

#### Wall Destruction Function
```
fw = g[(1 + Cw₃⁶)/(g⁶ + Cw₃⁶)]^(1/6)
g = r + Cw₂(r⁶ - r)
r = min[ν̃/(S̃κ²d²), 10]
```

#### Trip Function (for transition)
```
fₜ₂ = Cₜ₃exp(-Cₜ₄χ²)
```

#### Model Constants
- Cᵦ₁ = 0.1355
- Cᵦ₂ = 0.622
- σ = 2/3
- κ = 0.41 (von Karman constant)
- Cw₁ = Cᵦ₁/κ² + (1 + Cᵦ₂)/σ = 3.239
- Cw₂ = 0.3
- Cw₃ = 2.0
- Cᵥ₁ = 7.1
- Cₜ₃ = 1.2 (for trip)
- Cₜ₄ = 0.5 (for trip)

---

### Reynolds Stress Model (RSM)

Solves transport equations for all six Reynolds stress components plus dissipation rate.

#### Reynolds Stress Transport Equation
```
∂(ρu'ᵢu'ⱼ)/∂t + ∂(ρūₖu'ᵢu'ⱼ)/∂xₖ = Pᵢⱼ + Dᵢⱼ + Φᵢⱼ - εᵢⱼ
```

Where:
- Pᵢⱼ = Production
- Dᵢⱼ = Diffusion
- Φᵢⱼ = Pressure-strain
- εᵢⱼ = Dissipation

#### Production
```
Pᵢⱼ = -ρ(u'ᵢu'ₖ ∂ūⱼ/∂xₖ + u'ⱼu'ₖ ∂ūᵢ/∂xₖ)
```

#### Diffusion (modeled)
```
Dᵢⱼ = ∂/∂xₖ[(μ + ρCₛk/ε)∂u'ᵢu'ⱼ/∂xₖ]
```

#### Pressure-Strain (Linear model)
```
Φᵢⱼ = -C₁ρ(ε/k)[u'ᵢu'ⱼ - (2/3)kδᵢⱼ] - C₂[Pᵢⱼ - (2/3)Pₖδᵢⱼ]
```

Where Pₖ = (1/2)Pᵢᵢ

#### Dissipation
```
εᵢⱼ = (2/3)ρεδᵢⱼ
```

#### Dissipation Rate Equation (same as k-ε)
```
∂(ρε)/∂t + ∂(ρεūᵢ)/∂xᵢ = ∂/∂xⱼ[(μ + μₜ/σε)∂ε/∂xⱼ] + C₁ε(ε/k)Pₖ - C₂ερ(ε²/k)
```

#### Model Constants (Linear RSM)
- C₁ = 1.8
- C₂ = 0.6
- Cₛ = 0.22
- C₁ε = 1.44
- C₂ε = 1.92
- σε = 1.3

---

## LES Equations

### Filtered Navier-Stokes

#### Continuity
```
∂ρ̄/∂t + ∂(ρ̄ũᵢ)/∂xᵢ = 0
```

#### Momentum
```
∂(ρ̄ũᵢ)/∂t + ∂(ρ̄ũᵢũⱼ)/∂xⱼ = -∂p̄/∂xᵢ + ∂/∂xⱼ[μ(∂ũᵢ/∂xⱼ + ∂ũⱼ/∂xᵢ) - τᵢⱼˢᵍˢ]
```

Where:
- ũᵢ = filtered velocity
- τᵢⱼˢᵍˢ = subgrid-scale stress tensor (requires modeling)

### Subgrid-Scale Stress

```
τᵢⱼˢᵍˢ = ρ̄(uᵢuⱼ̃ - ũᵢũⱼ)
```

---

### Smagorinsky Model

#### Subgrid Viscosity
```
μₛₘₛ = ρ̄(CₛΔ)²|S̃|
```

Where:
```
|S̃| = √(2S̃ᵢⱼS̃ᵢⱼ)
S̃ᵢⱼ = (1/2)(∂ũᵢ/∂xⱼ + ∂ũⱼ/∂xᵢ)
Δ = (ΔxΔyΔz)^(1/3) (filter width)
```

#### Subgrid-Scale Stress
```
τᵢⱼˢᵍˢ - (1/3)τₖₖˢᵍˢδᵢⱼ = -2μₛₘₛS̃ᵢⱼ
```

#### Model Constant
- Cₛ = 0.1-0.2 (typically 0.17 for homogeneous turbulence)
- Cₛ ≈ 0.1 for wall-bounded flows

---

### Dynamic Smagorinsky Model

The Smagorinsky constant is computed dynamically:

```
Cₛ² = -⟨LᵢⱼMᵢⱼ⟩/(2⟨MᵢⱼMᵢⱼ⟩)
```

Where:
```
Lᵢⱼ = ũᵢũⱼ̂ - ũ̂ᵢũ̂ⱼ (Leonard stress)
Mᵢⱼ = Δ̂²|S̃̂|S̃ᵢⱼ̂ - Δ²|S̃|S̃ᵢⱼ^
```

The hat (^) denotes test filtering at scale Δ̂ = 2Δ.
⟨·⟩ denotes spatial averaging (homogeneous directions or local).

---

### WALE Model

#### Subgrid Viscosity
```
μₛₘₛ = ρ̄(CwΔ)²[(SᵢⱼᵈSᵢⱼᵈ)^(3/2)]/[(S̃ᵢⱼS̃ᵢⱼ)^(5/2) + (SᵢⱼᵈSᵢⱼᵈ)^(5/4)]
```

Where:
```
Sᵢⱼᵈ = (1/2)(g̃ᵢⱼ² + g̃ⱼᵢ²) - (1/3)δᵢⱼg̃ₖₖ²
g̃ᵢⱼ = ∂ũᵢ/∂xⱼ
```

#### Model Constant
- Cw = 0.325

**Advantage:** Returns correct y³ scaling near walls (μₛₘₛ ~ y³ as y → 0).

---

## Wall Functions

### Standard Wall Functions

Based on the law of the wall:

**Viscous sublayer (y⁺ < 5):**
```
u⁺ = y⁺
```

**Log-law region (y⁺ > 30):**
```
u⁺ = (1/κ)ln(y⁺) + B
```

Where:
- u⁺ = u/uτ
- y⁺ = ρuτy/μ
- uτ = √(τw/ρ) (friction velocity)
- κ = 0.41 (von Karman constant)
- B = 5.2 (smooth wall constant)

**Turbulent kinetic energy at wall:**
```
k = uτ²/√Cμ
```

**Dissipation rate at wall:**
```
ε = uτ³/(κy)
```

**Specific dissipation rate at wall:**
```
ω = uτ/(√Cμ κy)
```

### Scalable Wall Functions

Modified to prevent deterioration for y⁺ < 11.25:

```
y*⁺ = max(y⁺, 11.25)
```

Then use y*⁺ in place of y⁺ in standard wall functions.

### Enhanced Wall Treatment

Two-layer approach:

**Fully turbulent region (y⁺ > y⁺*):**
- Standard wall functions

**Near-wall region (y⁺ < y⁺*):**
- One-equation model or blended formulation
- y⁺* ≈ 11.25

Blending function:
```
φ = φwall·λ + φlog·(1 - λ)
```

Where λ varies from 1 (at wall) to 0 (at y⁺*).

---

## Validation Test Cases

### 1. Flat Plate Boundary Layer

**Description:** Zero-pressure-gradient turbulent boundary layer over flat plate

**Reynolds Number:** Re_x = 10⁶ - 10⁷

**Measured Quantities:**
- Skin friction coefficient: Cf
- Velocity profile: u⁺ vs y⁺
- Boundary layer thickness: δ, δ*, θ

**Validation Metrics:**
- Cf = 0.027/Re_θ^(1/7) (empirical correlation)
- Log-law slope in overlap region
- Wake parameter

**Best Models:**
- k-ω SST (excellent)
- Spalart-Allmaras (excellent)
- Realizable k-ε with EWT (good)
- Standard k-ε with wall functions (acceptable)

**Reference Data:** Coles (1968), Wieghardt & Tillmann

---

### 2. Channel Flow (DNS)

**Description:** Fully developed turbulent flow between parallel plates

**Reynolds Number:** Re_τ = uτδ/ν = 180, 395, 590, 1000

**Measured Quantities:**
- Mean velocity profile
- Reynolds stresses: u'u', v'v', w'w', u'v'
- Turbulent kinetic energy budget

**Validation Metrics:**
- Log-law region verification
- Peak Reynolds stress location and magnitude
- Turbulence production and dissipation profiles

**Best Models:**
- k-ω SST (very good)
- RSM (good for Reynolds stresses)
- k-ε variants (fair)

**Reference Data:**
- Moser, Kim & Mansour (1999) - Re_τ = 180, 395, 590
- Hoyas & Jiménez (2006) - Re_τ = 950, 2003

---

### 3. Backward-Facing Step

**Description:** Flow separation and reattachment over sudden expansion

**Geometry:** Expansion ratio ER = 1.2 - 2.0

**Reynolds Number:** Re_h = 5,000 - 50,000 (based on step height)

**Measured Quantities:**
- Reattachment length: Xr/h
- Separation bubble size
- Pressure distribution
- Velocity profiles

**Validation Metrics:**
- Reattachment length (critical metric)
- Skin friction distribution
- Turbulence intensity in recirculation zone

**Best Models:**
- k-ω SST (very good) - Xr within 5-10%
- Realizable k-ε (fair) - tends to overpredict Xr
- LES (excellent) - Xr within 2-5%
- Standard k-ε (poor) - significantly overpredicts Xr

**Reference Data:**
- Driver & Seegmiller (1985) - Re_h = 37,000, ER = 1.125
- Jovic & Driver (1994) - Re_h = 5,000-45,000

**Typical Results:**
- Experimental: Xr/h ≈ 6-7 (for Re_h = 37,000, ER = 1.125)
- k-ω SST: Xr/h ≈ 6.5
- Realizable k-ε: Xr/h ≈ 7.5-8.5
- Standard k-ε: Xr/h ≈ 9-11

---

### 4. Flow Over Circular Cylinder

**Description:** Separated flow with vortex shedding

**Reynolds Number:**
- Sub-critical: Re = 10⁴ - 10⁵
- Critical: Re = 3×10⁵
- Super-critical: Re > 3×10⁵

**Measured Quantities:**
- Drag coefficient: Cd
- Separation angle: θsep
- Strouhal number: St = fD/U
- Pressure distribution: Cp

**Validation Metrics:**
- Cd (critical parameter)
- Base pressure coefficient
- Separation point

**Best Models:**
- URANS k-ω SST (for unsteady simulation)
- DES/DDES (good)
- LES (excellent)
- Steady RANS (poor - cannot capture vortex shedding)

**Reference Data:**
- Roshko (1961)
- Achenbach (1968)
- Cantwell & Coles (1983)

**Typical Results (Re = 3.6×10⁴):**
- Experimental: Cd ≈ 1.0-1.2, St ≈ 0.20
- URANS k-ω SST: Cd ≈ 1.1, St ≈ 0.19-0.21
- LES: Cd ≈ 1.05, St ≈ 0.20

---

### 5. Axisymmetric Jet

**Description:** Round turbulent jet in co-flow or quiescent environment

**Reynolds Number:** Re_D = 10⁴ - 10⁶

**Measured Quantities:**
- Jet half-width growth rate: r1/2(x)
- Centerline velocity decay: Uc(x)
- Velocity profiles (self-similar region)
- Turbulence intensity

**Validation Metrics:**
- Spreading rate: dr1/2/dx
- Decay rate: dUc/dx
- Virtual origin location

**Best Models:**
- Realizable k-ε (excellent)
- RNG k-ε (very good)
- k-ω SST (good, but slightly slow spreading)
- Standard k-ε (overpredicts spreading)

**Reference Data:**
- Hussein, Capp & George (1994)
- Panchapakesan & Lumley (1993)

**Typical Results:**
- Experimental: dr1/2/dx ≈ 0.094
- Realizable k-ε: dr1/2/dx ≈ 0.095
- Standard k-ε: dr1/2/dx ≈ 0.110
- k-ω SST: dr1/2/dx ≈ 0.085

---

### 6. Flow in a Pipe (Fully Developed)

**Description:** Turbulent pipe flow (circumferentially averaged)

**Reynolds Number:** Re_D = 10⁴ - 10⁶

**Measured Quantities:**
- Velocity profile: u(r)
- Friction factor: f = Δp/(ρU²/2)(D/L)
- Turbulent stresses

**Validation Metrics:**
- Compare with Blasius correlation: f = 0.3164/Re_D^0.25
- Log-law region
- Centerline velocity

**Best Models:**
- All RANS models acceptable for mean velocity
- k-ω SST (best for friction factor)
- RSM (best for Reynolds stresses)

**Reference Data:**
- Laufer (1954)
- Zagarola & Smits (1998)

---

### 7. NASA Hump

**Description:** Wall-mounted smooth contour with separation

**Geometry:** 2D hump (bump) on flat plate

**Reynolds Number:** Re_c = 9.36×10⁵ (based on hump chord)

**Measured Quantities:**
- Pressure coefficient: Cp
- Separation and reattachment locations
- Velocity profiles
- Skin friction coefficient

**Validation Metrics:**
- Separation point: x/c ≈ 0.665
- Reattachment point: x/c ≈ 1.11
- Cp distribution

**Best Models:**
- k-ω SST (very good)
- DES (excellent)
- LES (benchmark)
- Standard k-ε (poor separation prediction)

**Reference Data:**
- NASA Langley Research Center (Greenblatt et al. 2006)
- 2004 CFD Validation Workshop

---

### 8. NACA 0012 Airfoil

**Description:** Symmetric airfoil at various angles of attack

**Reynolds Number:** Re_c = 3×10⁶ - 9×10⁶

**Angle of Attack:** α = 0° - 15°

**Measured Quantities:**
- Lift coefficient: Cl
- Drag coefficient: Cd
- Moment coefficient: Cm
- Pressure distribution: Cp
- Stall angle

**Validation Metrics:**
- Cl vs α (linear region and stall)
- Cd at various α
- Stall prediction

**Best Models:**
- k-ω SST (industry standard, excellent)
- Spalart-Allmaras (excellent)
- DES (for stall/post-stall)
- k-ε variants (poor for stall)

**Reference Data:**
- Abbott & Von Doenhoff (1959)
- Gregory & O'Reilly (1970)
- Ladson (1988) - NASA

**Typical Results (Re = 6×10⁶, α = 10°):**
- Experimental: Cl ≈ 1.08, Cd ≈ 0.010
- k-ω SST: Cl ≈ 1.06-1.09, Cd ≈ 0.009-0.011

---

### 9. Mixing Layer

**Description:** Plane mixing layer between two streams

**Velocity Ratio:** r = U2/U1 = 0.3 - 0.6

**Reynolds Number:** Re_θ = 30 - 300 (based on momentum thickness)

**Measured Quantities:**
- Growth rate: dδ/dx
- Velocity profiles (self-similar)
- Turbulence intensity
- Reynolds stresses

**Validation Metrics:**
- Spreading rate
- Self-similarity of profiles

**Best Models:**
- Realizable k-ε (very good)
- Standard k-ε (good)
- k-ω SST (acceptable, slow growth)

**Reference Data:**
- Bell & Mehta (1990)
- Champagne, Pao & Wygnanski (1976)

---

### 10. Swirling Flow (Turbulent Swirl)

**Description:** Pipe flow with swirl component

**Swirl Number:** S = 0.3 - 1.5

**Reynolds Number:** Re_D = 10⁴ - 10⁶

**Measured Quantities:**
- Axial and tangential velocity profiles
- Pressure drop
- Swirl decay rate
- Recirculation zones

**Validation Metrics:**
- Tangential velocity distribution
- Vortex breakdown prediction (high swirl)

**Best Models:**
- RNG k-ε (very good)
- RSM (excellent for complex swirl)
- k-ω SST (good)
- Standard k-ε (poor)

**Reference Data:**
- Reader-Harris (1994)
- Kitoh (1991)

---

## Summary: Model Performance by Test Case

| Test Case | k-ω SST | Realizable k-ε | RNG k-ε | Spalart-Allmaras | RSM | Standard k-ε |
|-----------|---------|----------------|---------|------------------|-----|--------------|
| Flat Plate | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ |
| Channel Flow | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |
| Backward Step | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| Cylinder | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ | ⭐⭐ | ⭐⭐ |
| Jet | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |
| Pipe Flow | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ |
| NASA Hump | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| NACA 0012 | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐ |
| Mixing Layer | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ |
| Swirl | ⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ |

⭐⭐⭐⭐⭐ = Excellent
⭐⭐⭐⭐ = Very Good
⭐⭐⭐ = Good/Acceptable
⭐⭐ = Fair/Poor
⭐ = Not Recommended

---

## Academic References

### Fundamental Turbulence Theory

1. **Pope, S.B.** (2000). *Turbulent Flows*. Cambridge University Press.
   - Comprehensive theoretical foundation
   - Statistical description of turbulence
   - Closure modeling approaches

2. **Wilcox, D.C.** (2006). *Turbulence Modeling for CFD* (3rd ed.). DCW Industries.
   - Practical guide to turbulence models
   - Detailed model derivations
   - Validation and application

3. **Tennekes, H., & Lumley, J.L.** (1972). *A First Course in Turbulence*. MIT Press.
   - Classic introduction to turbulence
   - Physical understanding of turbulent flows

4. **Hinze, J.O.** (1975). *Turbulence* (2nd ed.). McGraw-Hill.
   - Fundamental turbulence concepts
   - Experimental observations

5. **Durbin, P.A., & Pettersson Reif, B.A.** (2011). *Statistical Theory and Modeling for Turbulent Flows* (2nd ed.). Wiley.
   - Advanced turbulence modeling
   - Modern RANS approaches

### k-ε Models

6. **Launder, B.E., & Spalding, D.B.** (1974). "The numerical computation of turbulent flows." *Computer Methods in Applied Mechanics and Engineering*, 3(2), 269-289.
   - Original standard k-ε model

7. **Yakhot, V., & Orszag, S.A.** (1986). "Renormalization group analysis of turbulence. I. Basic theory." *Journal of Scientific Computing*, 1(1), 3-51.
   - RNG k-ε model derivation

8. **Shih, T.-H., Liou, W.W., Shabbir, A., Yang, Z., & Zhu, J.** (1995). "A new k-ε eddy viscosity model for high Reynolds number turbulent flows." *Computers & Fluids*, 24(3), 227-238.
   - Realizable k-ε model

### k-ω Models

9. **Wilcox, D.C.** (1988). "Reassessment of the scale-determining equation for advanced turbulence models." *AIAA Journal*, 26(11), 1299-1310.
   - Original k-ω model

10. **Wilcox, D.C.** (2008). "Formulation of the k-ω turbulence model revisited." *AIAA Journal*, 46(11), 2823-2838.
    - Updated k-ω (2006) model

11. **Menter, F.R.** (1994). "Two-equation eddy-viscosity turbulence models for engineering applications." *AIAA Journal*, 32(8), 1598-1605.
    - k-ω SST model (landmark paper)

12. **Menter, F.R., Kuntz, M., & Langtry, R.** (2003). "Ten years of industrial experience with the SST turbulence model." *Turbulence, Heat and Mass Transfer*, 4, 625-632.
    - SST model applications and validation

### Spalart-Allmaras Model

13. **Spalart, P.R., & Allmaras, S.R.** (1994). "A one-equation turbulence model for aerodynamic flows." *La Recherche Aérospatiale*, 1, 5-21.
    - Original S-A model

14. **Spalart, P.R.** (2000). "Strategies for turbulence modelling and simulations." *International Journal of Heat and Fluid Flow*, 21(3), 252-263.
    - S-A model variants and applications

### Reynolds Stress Models

15. **Launder, B.E., Reece, G.J., & Rodi, W.** (1975). "Progress in the development of a Reynolds-stress turbulence closure." *Journal of Fluid Mechanics*, 68(3), 537-566.
    - Linear RSM (LRR model)

16. **Speziale, C.G., Sarkar, S., & Gatski, T.B.** (1991). "Modelling the pressure-strain correlation of turbulence: an invariant dynamical systems approach." *Journal of Fluid Mechanics*, 227, 245-272.
    - SSG Reynolds stress model

### Large Eddy Simulation

17. **Smagorinsky, J.** (1963). "General circulation experiments with the primitive equations: I. The basic experiment." *Monthly Weather Review*, 91(3), 99-164.
    - Smagorinsky SGS model

18. **Germano, M., Piomelli, U., Moin, P., & Cabot, W.H.** (1991). "A dynamic subgrid-scale eddy viscosity model." *Physics of Fluids A*, 3(7), 1760-1765.
    - Dynamic SGS model

19. **Nicoud, F., & Ducros, F.** (1999). "Subgrid-scale stress modelling based on the square of the velocity gradient tensor." *Flow, Turbulence and Combustion*, 62(3), 183-200.
    - WALE SGS model

20. **Piomelli, U., & Chasnov, J.R.** (1996). "Large-eddy simulations: theory and applications." In *Turbulence and Transition Modelling* (pp. 269-336). Springer.
    - Comprehensive LES review

### Hybrid RANS-LES

21. **Spalart, P.R., Jou, W.-H., Strelets, M., & Allmaras, S.R.** (1997). "Comments on the feasibility of LES for wings, and on a hybrid RANS/LES approach." *Advances in DNS/LES*, 1, 4-8.
    - Detached Eddy Simulation (DES)

22. **Spalart, P.R., Deck, S., Shur, M.L., Squires, K.D., Strelets, M.K., & Travin, A.** (2006). "A new version of detached-eddy simulation, resistant to ambiguous grid densities." *Theoretical and Computational Fluid Dynamics*, 20(3), 181-195.
    - Delayed DES (DDES)

23. **Menter, F.R., & Egorov, Y.** (2010). "The scale-adaptive simulation method for unsteady turbulent flow predictions. Part 1: Theory and model description." *Flow, Turbulence and Combustion*, 85(1), 113-138.
    - Scale-Adaptive Simulation (SAS)

### Wall Functions and Near-Wall Modeling

24. **Launder, B.E., & Spalding, D.B.** (1974). "The numerical computation of turbulent flows." (Same as reference 6)
    - Standard wall functions

25. **Grotjans, H., & Menter, F.R.** (1998). "Wall functions for general application CFD codes." *ECCOMAS 98 Proceedings*, 4, 1112-1117.
    - Automatic wall treatment

26. **Kader, B.A.** (1981). "Temperature and concentration profiles in fully turbulent boundary layers." *International Journal of Heat and Mass Transfer*, 24(9), 1541-1544.
    - Enhanced wall treatment theory

### Validation and Benchmarking

27. **Driver, D.M., & Seegmiller, H.L.** (1985). "Features of a reattaching turbulent shear layer in divergent channel flow." *AIAA Journal*, 23(2), 163-171.
    - Backward-facing step experiments

28. **Moser, R.D., Kim, J., & Mansour, N.N.** (1999). "Direct numerical simulation of turbulent channel flow up to Re_τ = 590." *Physics of Fluids*, 11(4), 943-945.
    - DNS channel flow benchmark

29. **Hussein, H.J., Capp, S.P., & George, W.K.** (1994). "Velocity measurements in a high-Reynolds-number, momentum-conserving, axisymmetric, turbulent jet." *Journal of Fluid Mechanics*, 258, 31-75.
    - Round jet experiments

30. **Greenblatt, D., Paschal, K.B., Yao, C.-S., Harris, J., Schaeffler, N.W., & Washburn, A.E.** (2006). "Experimental investigation of separation control part 1: Baseline and steady suction." *AIAA Journal*, 44(12), 2820-2830.
    - NASA wall-mounted hump

### CFD Textbooks with Turbulence Modeling

31. **Versteeg, H.K., & Malalasekera, W.** (2007). *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* (2nd ed.). Pearson.
    - Practical CFD with turbulence models

32. **Ferziger, J.H., & Perić, M.** (2002). *Computational Methods for Fluid Dynamics* (3rd ed.). Springer.
    - Numerical methods and turbulence

33. **Blazek, J.** (2015). *Computational Fluid Dynamics: Principles and Applications* (3rd ed.). Butterworth-Heinemann.
    - Industrial CFD applications

### Standards and Best Practice Guidelines

34. **ERCOFTAC** (European Research Community on Flow, Turbulence and Combustion). *Best Practice Guidelines*.
    - Available at: http://www.ercoftac.org/
    - Quality and trust in CFD

35. **AIAA** (American Institute of Aeronautics and Astronautics). "Guide for the Verification and Validation of Computational Fluid Dynamics Simulations." AIAA G-077-1998.
    - V&V standards

36. **NASA** Turbulence Modeling Resource: https://turbmodels.larc.nasa.gov/
    - Model verification
    - Validation test cases
    - Code implementations

---

## Online Resources

### Databases and Repositories

1. **NASA Turbulence Modeling Resource**
   - https://turbmodels.larc.nasa.gov/
   - Verification cases, validation data, model implementations

2. **ERCOFTAC Database**
   - http://www.ercoftac.org/
   - Classic experimental test cases

3. **CFD Online**
   - https://www.cfd-online.com/Wiki/Turbulence_models
   - Community wiki and forums

### Software Documentation

4. **ANSYS Fluent Theory Guide**
   - Detailed turbulence model implementations
   - Wall treatment descriptions

5. **OpenFOAM User Guide**
   - Open-source turbulence model implementations
   - Source code available

6. **STAR-CCM+ Documentation**
   - Commercial CFD turbulence models

---

## Model Implementation Notes

### Numerical Considerations

1. **Under-Relaxation Factors:**
   - k, ε, ω: 0.5-0.8
   - Turbulent viscosity: 0.8-1.0
   - Start lower, increase as convergence progresses

2. **Initialization:**
   - Turbulent intensity: I = 1-5%
   - Turbulent length scale: l ≈ 0.07L
   - Or turbulent viscosity ratio: μₜ/μ = 1-10

3. **Boundary Conditions:**
   - Inlet: Specify I and l (or k and ε/ω)
   - Wall: Automatic (no-slip + wall treatment)
   - Outlet: Zero gradient
   - Symmetry: Zero gradient

4. **Convergence Criteria:**
   - Residuals: < 10⁻⁴ (momentum), < 10⁻⁴ (turbulence)
   - Monitor: Forces, mass flow, separation point
   - Ensure all monitored quantities stabilize

---

## Quick Reference: When to Use Which Model

### Use k-ω SST when:
- External aerodynamics (airfoils, vehicles, aircraft)
- Accurate separation prediction needed
- Heat transfer is important
- Adverse pressure gradients present
- You can afford y+ ≈ 1 mesh

### Use Realizable k-ε when:
- Internal flows (pipes, ducts, HVAC)
- Jets, wakes, mixing flows
- Large domains (computational efficiency matters)
- Wall functions acceptable (y+ = 30-100)
- Buoyancy-driven flows

### Use Spalart-Allmaras when:
- Aerospace applications
- External aerodynamics (single-element airfoils)
- Need lowest computational cost
- Attached or mildly separated flows

### Use RNG k-ε when:
- Swirling or rotating flows
- Flows with streamline curvature
- Need improvement over standard k-ε
- Can use wall functions

### Use RSM when:
- Highly anisotropic turbulence
- Strong swirl, rotation, or curvature
- Complex 3D flows
- Have computational resources
- RANS models fail

### Use LES when:
- Acoustics (noise prediction)
- Vortex-dominated flows
- Combustion/mixing
- Have massive computational resources
- Time-accurate structures needed

### Use DES/DDES when:
- Massively separated flows
- External aerodynamics at high angles of attack
- Can't afford full LES
- Need both attached boundary layers and separated regions

---

This reference should be used in conjunction with SKILL.md for practical turbulence model selection and setup in CFD simulations.
