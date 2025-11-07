# OpenFOAM Boundary Conditions Examples

Comprehensive collection of boundary condition examples for common CFD scenarios in OpenFOAM.

## Table of Contents
- [Velocity Boundary Conditions](#velocity-boundary-conditions)
- [Pressure Boundary Conditions](#pressure-boundary-conditions)
- [Turbulence Boundary Conditions](#turbulence-boundary-conditions)
- [Temperature Boundary Conditions](#temperature-boundary-conditions)
- [Advanced Boundary Conditions](#advanced-boundary-conditions)

---

## Velocity Boundary Conditions

### Fixed Inlet Velocity

```cpp
inlet
{
    type            fixedValue;
    value           uniform (1 0 0);  // 1 m/s in x-direction
}
```

### Parabolic Inlet Profile (Poiseuille Flow)

```cpp
inlet
{
    type            codedFixedValue;
    value           uniform (0 0 0);
    name            parabolicVelocity;
    code
    #{
        const fvPatch& boundaryPatch = patch();
        const vectorField& Cf = boundaryPatch.Cf();
        vectorField& field = *this;

        scalar H = 0.1;      // Channel height
        scalar Umax = 1.5;   // Maximum velocity

        forAll(Cf, faceI)
        {
            scalar y = Cf[faceI].y();
            scalar u = Umax * 4 * y * (H - y) / (H * H);
            field[faceI] = vector(u, 0, 0);
        }
    #};
}
```

### Power Law Velocity Profile (Turbulent Pipe)

```cpp
inlet
{
    type            codedFixedValue;
    value           uniform (0 0 0);
    name            powerLawProfile;
    code
    #{
        const fvPatch& boundaryPatch = patch();
        const vectorField& Cf = boundaryPatch.Cf();
        vectorField& field = *this;

        scalar R = 0.05;     // Pipe radius
        scalar Umax = 2.0;   // Centerline velocity
        scalar n = 7.0;      // Power law exponent (1/7 for turbulent)

        forAll(Cf, faceI)
        {
            scalar r = sqrt(sqr(Cf[faceI].y()) + sqr(Cf[faceI].z()));
            scalar u = Umax * pow(1.0 - r/R, 1.0/n);
            field[faceI] = vector(u, 0, 0);
        }
    #};
}
```

### Time-Varying Inlet (Pulsatile Flow)

```cpp
inlet
{
    type            uniformFixedValue;
    uniformValue    table
    (
        (0      (0 0 0))
        (0.5    (1 0 0))
        (1.0    (2 0 0))
        (1.5    (1 0 0))
        (2.0    (0 0 0))
    );
}
```

### Outlet - Zero Gradient

```cpp
outlet
{
    type            zeroGradient;
}
```

### Outlet - Pressure Inlet/Outlet

```cpp
outlet
{
    type            pressureInletOutletVelocity;
    value           uniform (0 0 0);
}
```

### Wall - No Slip

```cpp
wall
{
    type            noSlip;
}
```

### Wall - Slip

```cpp
wall
{
    type            slip;
}
```

### Wall - Moving Wall

```cpp
movingWall
{
    type            movingWallVelocity;
    value           uniform (1 0 0);  // Wall velocity
}
```

### Wall - Rotating Wall

```cpp
rotatingWall
{
    type            rotatingWallVelocity;
    origin          (0 0 0);
    axis            (0 0 1);
    omega           10.47;  // 100 RPM = 10.47 rad/s
}
```

### Symmetry Plane

```cpp
symmetry
{
    type            symmetryPlane;
}
```

### Wedge (Axisymmetric)

```cpp
wedge
{
    type            wedge;
}
```

---

## Pressure Boundary Conditions

### Fixed Value (Outlet)

```cpp
outlet
{
    type            fixedValue;
    value           uniform 0;
}
```

### Zero Gradient (Inlet/Wall)

```cpp
inlet
{
    type            zeroGradient;
}

wall
{
    type            zeroGradient;
}
```

### Total Pressure (Inlet)

```cpp
inlet
{
    type            totalPressure;
    p0              uniform 101325;  // Total pressure (Pa)
    U               U;
    phi             phi;
    rho             none;
    psi             none;
    gamma           1;
    value           uniform 101325;
}
```

### Fan Pressure Jump

```cpp
fan
{
    type            fan;
    patchType       cyclic;
    jump            uniform 100;  // Pressure jump (Pa)
    value           uniform 0;
    jumpTable       constant 100;
}
```

### Pressure Inlet/Outlet Velocity

```cpp
outlet
{
    type            inletOutlet;
    inletValue      uniform 0;
    value           uniform 0;
}
```

---

## Turbulence Boundary Conditions

### Turbulent Kinetic Energy (k)

**Inlet - From Intensity**
```cpp
inlet
{
    type            turbulentIntensityKineticEnergyInlet;
    intensity       0.05;  // 5% turbulence intensity
    value           uniform 0.0015;
}
```

**Inlet - Fixed Value**
```cpp
inlet
{
    type            fixedValue;
    value           uniform 0.0015;
}
```

**Outlet - Zero Gradient**
```cpp
outlet
{
    type            zeroGradient;
}
```

**Wall - Wall Function**
```cpp
wall
{
    type            kqRWallFunction;
    value           uniform 0.0015;
}
```

**Wall - Low-Re (Resolved)**
```cpp
wall
{
    type            fixedValue;
    value           uniform 1e-10;
}
```

### Dissipation Rate (epsilon)

**Inlet - From Mixing Length**
```cpp
inlet
{
    type            turbulentMixingLengthDissipationRateInlet;
    mixingLength    0.005;  // 5 mm
    value           uniform 0.01;
}
```

**Inlet - Fixed Value**
```cpp
inlet
{
    type            fixedValue;
    value           uniform 0.01;
}
```

**Wall - Wall Function**
```cpp
wall
{
    type            epsilonWallFunction;
    value           uniform 0.01;
}
```

**Wall - Low-Re (Resolved)**
```cpp
wall
{
    type            epsilonLowReWallFunction;
    value           uniform 0.01;
}
```

### Specific Dissipation Rate (omega)

**Inlet - From Mixing Length**
```cpp
inlet
{
    type            turbulentMixingLengthFrequencyInlet;
    mixingLength    0.005;
    value           uniform 1;
}
```

**Inlet - Fixed Value**
```cpp
inlet
{
    type            fixedValue;
    value           uniform 1;
}
```

**Wall - Wall Function**
```cpp
wall
{
    type            omegaWallFunction;
    value           uniform 1;
}
```

### Turbulent Viscosity (nut)

**Inlet - Calculated**
```cpp
inlet
{
    type            calculated;
    value           uniform 0;
}
```

**Wall - Wall Function**
```cpp
wall
{
    type            nutkWallFunction;
    value           uniform 0;
}
```

**Wall - Rough Wall**
```cpp
wall
{
    type            nutkRoughWallFunction;
    Ks              uniform 1e-4;  // Roughness height (m)
    Cs              uniform 0.5;   // Roughness constant
    value           uniform 0;
}
```

---

## Temperature Boundary Conditions

### Fixed Temperature

```cpp
hotWall
{
    type            fixedValue;
    value           uniform 373;  // 100°C in Kelvin
}
```

### Fixed Heat Flux

```cpp
heatedWall
{
    type            fixedGradient;
    gradient        uniform 1000;  // Temperature gradient (K/m)
}
```

### Convective Heat Transfer

```cpp
wall
{
    type            externalWallHeatFluxTemperature;
    mode            coefficient;
    h               uniform 10;      // Heat transfer coefficient (W/m²K)
    Ta              constant 300;    // Ambient temperature (K)
    value           uniform 300;
    kappaMethod     fluidThermo;
}
```

### Adiabatic Wall

```cpp
adiabaticWall
{
    type            zeroGradient;
}
```

### Conjugate Heat Transfer Interface

```cpp
interface
{
    type            compressible::turbulentTemperatureCoupledBaffleMixed;
    Tnbr            T;
    kappaMethod     fluidThermo;
    value           uniform 300;
}
```

---

## Advanced Boundary Conditions

### Inlet from CSV File

```cpp
inlet
{
    type            timeVaryingMappedFixedValue;
    offset          (0 0 0);
    setAverage      off;
}
```

With `constant/boundaryData/inlet/points` and `time/U` files.

### Mapped Boundary (Periodic)

```cpp
inlet
{
    type            mapped;
    sampleMode      nearestCell;
    sampleRegion    region0;
    samplePatch     outlet;
    offset          (0 0 0);
    value           uniform (1 0 0);
}
```

### Atmospheric Boundary Layer

```cpp
ground
{
    type            atmBoundaryLayerInletVelocity;
    flowDir         (1 0 0);
    zDir            (0 0 1);
    Uref            10;       // Reference velocity (m/s)
    Zref            10;       // Reference height (m)
    z0              uniform 0.1;  // Roughness length (m)
    d               uniform 0;    // Zero plane displacement
    value           uniform (10 0 0);
}
```

### Wave Boundary (Multiphase)

```cpp
inlet
{
    type            waveVelocity;
    value           uniform (0 0 0);
}
```

### Free Stream

```cpp
farfield
{
    type            freestreamVelocity;
    freestreamValue uniform (10 0 0);
}
```

### Supersonic Inlet

```cpp
inlet
{
    type            supersonicFreestream;
    pInf            101325;
    TInf            300;
    UInf            (680 0 0);  // Mach 2 at 300K
    gamma           1.4;
}
```

### Outflow

```cpp
outlet
{
    type            outletPhaseMeanVelocity;
    Umean           1.5;  // Mean velocity
    alpha           alpha.water;
    value           uniform (1.5 0 0);
}
```

---

## Calculating Turbulence Parameters

### Formula Reference

```
Turbulence intensity:    I = u'/U = 0.16 * Re^(-1/8)
Turbulent kinetic energy: k = 3/2 * (U * I)²
Dissipation rate:        ε = C_μ^(3/4) * k^(3/2) / L
Specific dissipation:    ω = k^(1/2) / (C_μ^(1/4) * L)
Turbulent viscosity:     νt = C_μ * k² / ε

Where:
  U = mean velocity
  Re = Reynolds number
  L = turbulent length scale ≈ 0.07 * characteristic_length
  C_μ = 0.09
```

### Example Calculation

For a pipe with D=0.1m, U=1m/s, ν=1.5e-5 m²/s (air):

```python
Re = U * D / ν = 1 * 0.1 / 1.5e-5 = 6667
I = 0.16 * Re^(-1/8) = 0.16 * 6667^(-1/8) = 0.0516 (5.16%)
k = 3/2 * (U * I)² = 1.5 * (1 * 0.0516)² = 0.00399 m²/s²
L = 0.07 * D = 0.007 m
ε = 0.09^(3/4) * k^(3/2) / L = 0.164 * 0.00399^1.5 / 0.007 = 0.0412 m²/s³
ω = k^0.5 / (0.09^0.25 * L) = 0.0632 / (0.548 * 0.007) = 16.5 1/s
```

Use the provided `calculate_turbulence.py` script for automatic calculation.

---

## Common Boundary Condition Combinations

### Internal Flow (Pipe/Channel)

**Velocity-Pressure**
- Inlet: Fixed velocity, zero gradient pressure
- Outlet: Zero gradient velocity, fixed pressure
- Walls: No-slip velocity, zero gradient pressure

### External Flow (Aerodynamics)

**Velocity-Pressure**
- Inlet: Fixed velocity, zero gradient pressure
- Outlet: Zero gradient velocity, zero gradient pressure
- Walls: No-slip velocity, zero gradient pressure
- Farfield: Freestream conditions

### Natural Convection

**Velocity-Pressure-Temperature**
- All walls: No-slip velocity, zero gradient pressure, fixed temperature
- Hot wall: Fixed high temperature
- Cold wall: Fixed low temperature

### Rotating Machinery

**Velocity-Pressure**
- Inlet: Fixed velocity, zero gradient pressure
- Outlet: Zero gradient velocity, fixed pressure
- Rotating walls: Moving wall velocity, zero gradient pressure
- Stationary walls: No-slip velocity, zero gradient pressure

---

## Tips and Best Practices

1. **Consistency**: Ensure boundary conditions are physically consistent
2. **Convergence**: Use gradual changes for better convergence
3. **Turbulence**: Match k, epsilon/omega values to flow conditions
4. **Walls**: Use wall functions for y+ > 30, resolve for y+ < 1
5. **Pressure**: At least one pressure boundary must be specified
6. **Validation**: Compare with experimental or analytical solutions

---

## References

- OpenFOAM User Guide: Boundary Conditions
- OpenFOAM Programmer's Guide: Boundary Condition Implementation
- CFD Online: OpenFOAM Boundary Conditions
