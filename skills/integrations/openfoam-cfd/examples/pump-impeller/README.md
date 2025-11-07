# Pump Impeller Simulation

## Description

This case demonstrates a transient simulation of flow through a simplified pump impeller using the k-omega SST turbulence model and the pimpleFoam solver. The impeller rotates at 1000 RPM using the Multiple Reference Frame (MRF) or sliding mesh approach.

**Note**: This is a simplified 2D representation. For actual pump design, a full 3D geometry with proper blade profiles is recommended.

## Geometry

- **Inner radius**: 0.05 m (impeller hub)
- **Outer radius**: 0.15 m (casing)
- **Axial height**: 0.2 m
- **Depth**: 0.01 m (2D simulation, single cell in z-direction)
- **Mesh**: 40 x 40 x 1 cells
- **Rotation speed**: 1000 RPM (104.72 rad/s)

## Flow Conditions

- **Inlet velocity**: 2 m/s (axial)
- **Kinematic viscosity**: 1e-6 m²/s (water at 20°C)
- **Fluid density**: 1000 kg/m³
- **Reynolds number**: ~200,000 (based on blade chord)
- **Turbulence model**: k-omega SST

## Boundary Conditions

### Velocity (U)
- **inlet**: Fixed value (0, 0, 2) m/s (axial inlet)
- **outlet**: Zero gradient
- **impeller**: Moving wall velocity (rotates with impeller)
- **casing**: No slip (stationary)
- **frontAndBack**: Empty (2D)

### Pressure (p)
- **inlet**: Zero gradient
- **outlet**: Fixed value (0 Pa)
- **impeller**: Zero gradient
- **casing**: Zero gradient
- **frontAndBack**: Empty (2D)

### Turbulence (k, omega)
- **inlet**: Calculated from turbulence intensity and length scale
- **outlet**: Zero gradient
- **impeller**: Wall functions
- **casing**: Wall functions
- **frontAndBack**: Empty (2D)

## Rotating Mesh Setup

The impeller rotation is handled using the `solidBodyMotionFvMesh` with rotatingMotion:
- **Origin**: (0, 0, 0)
- **Axis**: (0, 0, 1) - z-axis
- **Angular velocity**: 104.72 rad/s (1000 RPM)

## How to Run

```bash
# Navigate to case directory
cd pump-impeller

# Generate mesh
blockMesh

# Check mesh quality
checkMesh

# Run the solver (transient)
pimpleFoam

# View results
paraFoam
```

## Expected Results

The simulation should show:
1. Axial flow entering at the inlet
2. Flow acceleration due to rotation
3. Pressure increase across the impeller
4. Turbulent wake formation
5. Transient vortex structures

## Post-Processing

### In ParaView:
- Velocity vectors showing swirl
- Pressure contours showing head rise
- Vorticity magnitude
- Turbulent kinetic energy
- Streamlines colored by velocity magnitude

### Performance Metrics:
```bash
# View forces on impeller
cat postProcessing/forces/0/forces.dat

# Calculate torque and power
# Torque = Force × Radius
# Power = Torque × Angular velocity
```

## Typical Runtime

- Mesh generation: < 1 second
- Simulation: 10-60 minutes (depending on hardware)
- Time steps: ~10,000 for 1 second simulation

## Important Notes

### MRF vs Sliding Mesh

**This example uses sliding mesh (actual rotation)**. For steady-state analysis, consider using MRF (Multiple Reference Frame) instead:

1. Use `simpleFoam` with MRF
2. Add `MRFProperties` dictionary
3. Define rotating zone
4. Much faster convergence for steady operation

### 3D Considerations

For realistic pump simulations:
- Create 3D blade geometry (STL file)
- Use `snappyHexMesh` for meshing
- Include volute or diffuser
- Perform grid independence study
- Validate against experimental data

## Modifications

To change operating conditions:

**Rotation Speed:**
```cpp
// constant/dynamicMeshDict
omega  constant 209.44;  // 2000 RPM
```

**Flow Rate:**
```cpp
// 0/U inlet
value  uniform (0 0 5);  // 5 m/s
```

**Turbulence Model:**
```cpp
// constant/turbulenceProperties
RASModel  kEpsilon;  // or realizableKE
```

**Mesh Resolution:**
```cpp
// system/blockMeshDict
hex (...) (80 80 1) ...  // Double resolution
```

## Troubleshooting

**Large Courant Number:**
- Reduce time step in controlDict
- Use adjustable time stepping (already enabled)

**Divergence:**
- Reduce omega (rotation speed)
- Improve initial conditions
- Increase mesh resolution near walls

**Mesh Motion Errors:**
- Check dynamicMeshDict syntax
- Ensure impeller patch is correctly defined
- Verify cell zone exists

## Validation

Compare results with:
- Pump performance curves (head vs. flow rate)
- Manufacturer data
- CFD benchmarks
- Experimental measurements (PIV, LDA)

## References

- OpenFOAM rotating machinery tutorials
- Pump design handbooks
- Turbomachinery CFD literature
- ANSYS/Fluent validation cases for comparison
