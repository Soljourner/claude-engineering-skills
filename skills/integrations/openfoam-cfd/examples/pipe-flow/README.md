# Turbulent Pipe Flow Example

## Description

This case demonstrates a turbulent flow simulation through a rectangular pipe using the k-epsilon turbulence model and the simpleFoam solver for steady-state incompressible flow.

## Geometry

- **Length**: 3 meters
- **Height**: 0.1 meters
- **Depth**: 0.02 meters (2D simulation, single cell in z-direction)
- **Mesh**: 150 x 50 x 1 cells

## Flow Conditions

- **Inlet velocity**: 1 m/s in x-direction
- **Kinematic viscosity**: 1e-5 m²/s (air at 20°C)
- **Reynolds number**: ~10,000 (based on height)
- **Turbulence intensity**: 5%
- **Turbulent length scale**: 0.005 m

## Boundary Conditions

### Velocity (U)
- **inlet**: Fixed value (1, 0, 0) m/s
- **outlet**: Zero gradient
- **wall**: No slip
- **frontAndBack**: Empty (2D)

### Pressure (p)
- **inlet**: Zero gradient
- **outlet**: Fixed value (0 Pa)
- **wall**: Zero gradient
- **frontAndBack**: Empty (2D)

### Turbulence (k, epsilon)
- **inlet**: Calculated from turbulence intensity and length scale
- **outlet**: Zero gradient
- **wall**: Wall functions
- **frontAndBack**: Empty (2D)

## How to Run

```bash
# Navigate to case directory
cd pipe-flow

# Generate mesh
blockMesh

# Check mesh quality
checkMesh

# Run the solver
simpleFoam

# View results
paraFoam
```

## Expected Results

The flow should:
1. Develop from the inlet profile
2. Show turbulent boundary layers near walls
3. Reach fully developed flow downstream
4. Converge residuals to < 1e-5

## Post-Processing

View the results in ParaView:
- Velocity magnitude contours
- Pressure distribution
- Turbulent kinetic energy
- Wall shear stress

## Monitoring

Check convergence:
```bash
# Monitor residuals
tail -f log.simpleFoam

# Plot residuals (if gnuplot installed)
foamLog log.simpleFoam
gnuplot -persist postProcessing/residuals/0/residuals.gnuplot
```

## Typical Runtime

- Mesh generation: < 1 second
- Simulation: 1-5 minutes (depending on hardware)
- Convergence: ~500-1000 iterations

## Modifications

To change flow conditions, edit:
- `0/U` - Change inlet velocity
- `constant/transportProperties` - Change viscosity
- `0/k` and `0/epsilon` - Adjust turbulence parameters
- `system/blockMeshDict` - Modify geometry and mesh resolution
