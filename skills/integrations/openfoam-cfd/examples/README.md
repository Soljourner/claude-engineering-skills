# OpenFOAM Examples

This directory contains complete, working OpenFOAM case examples and utility scripts.

## Directory Structure

```
examples/
├── pipe-flow/              # Turbulent pipe flow case (simpleFoam)
├── pump-impeller/          # Rotating impeller case (pimpleFoam)
├── scripts/                # Python utilities for case generation
│   ├── generate_case.py           # Automated case generator
│   └── calculate_turbulence.py    # Turbulence parameter calculator
└── boundary-conditions-examples.md # Comprehensive BC reference
```

## Quick Start

### Example 1: Pipe Flow

Run a complete turbulent pipe flow simulation:

```bash
cd pipe-flow
blockMesh
checkMesh
simpleFoam
paraFoam
```

**Description**: Steady-state turbulent flow through a rectangular pipe using k-epsilon model.
**Runtime**: ~2-5 minutes
**Key Learning**: Basic OpenFOAM workflow, k-epsilon turbulence

See [pipe-flow/README.md](pipe-flow/README.md) for details.

### Example 2: Pump Impeller

Run a rotating machinery simulation:

```bash
cd pump-impeller
blockMesh
checkMesh
pimpleFoam
paraFoam
```

**Description**: Transient flow through a simplified pump impeller with rotation.
**Runtime**: ~10-60 minutes
**Key Learning**: Transient solver, rotating meshes, k-omega SST

See [pump-impeller/README.md](pump-impeller/README.md) for details.

## Utility Scripts

### Case Generator

Automatically generate OpenFOAM case directories:

```bash
python3 scripts/generate_case.py myNewCase \
    --velocity 2 0 0 \
    --nu 1.5e-5 \
    --turbulence kOmegaSST
```

**Features**:
- Creates complete directory structure (0/, constant/, system/)
- Generates all required dictionaries
- Supports k-epsilon and k-omega SST models
- Configurable flow parameters

**Usage**:
```bash
./scripts/generate_case.py --help
```

### Turbulence Calculator

Calculate turbulence parameters for boundary conditions:

```bash
python3 scripts/calculate_turbulence.py \
    --velocity 1.0 \
    --length 0.1 \
    --nu 1.5e-5
```

**Calculates**:
- Turbulence intensity (I)
- Turbulent kinetic energy (k)
- Dissipation rate (epsilon)
- Specific dissipation rate (omega)
- Mixing length
- Turbulent viscosity

**Output**: Ready-to-use OpenFOAM boundary condition snippets

**Usage**:
```bash
./scripts/calculate_turbulence.py --help
```

## Boundary Conditions Reference

The file [boundary-conditions-examples.md](boundary-conditions-examples.md) contains comprehensive examples of:

- **Velocity BCs**: Fixed value, profiles, time-varying, rotating walls
- **Pressure BCs**: Fixed value, total pressure, fan jump
- **Turbulence BCs**: k, epsilon, omega for inlets, outlets, walls
- **Temperature BCs**: Fixed value, heat flux, conjugate heat transfer
- **Advanced BCs**: Mapped, atmospheric, waves, supersonic

## Example Modifications

### Change Flow Velocity

Edit the inlet boundary condition in `0/U`:

```cpp
inlet
{
    type            fixedValue;
    value           uniform (2 0 0);  // Change from (1 0 0) to (2 0 0)
}
```

Recalculate turbulence parameters:

```bash
python3 ../scripts/calculate_turbulence.py --velocity 2.0 --length 0.1 --nu 1.5e-5
```

Update `0/k` and `0/epsilon` with new values.

### Change Fluid Properties

Edit `constant/transportProperties`:

```cpp
nu              [0 2 -1 0 0 0 0] 1e-06;  // Change to water
```

### Increase Mesh Resolution

Edit `system/blockMeshDict`:

```cpp
hex (0 1 2 3 4 5 6 7) (300 100 1) simpleGrading (1 1 1)  // Double from (150 50 1)
```

Regenerate mesh:

```bash
blockMesh
checkMesh
```

### Change Turbulence Model

Edit `constant/turbulenceProperties`:

```cpp
RASModel        kOmegaSST;  // Change from kEpsilon
```

Add `0/omega` field and adjust boundary conditions accordingly.

## Creating Custom Cases

### Method 1: Use the Generator Script

```bash
python3 scripts/generate_case.py myCase \
    --velocity 5 0 0 \
    --nu 1e-6 \
    --turbulence kOmegaSST \
    --dir ~/OpenFOAM/run
```

Then:
1. Create/edit `system/blockMeshDict` for your geometry
2. Run `blockMesh`
3. Run solver
4. Post-process

### Method 2: Copy and Modify Existing Case

```bash
cp -r pipe-flow myNewCase
cd myNewCase
# Edit dictionaries as needed
blockMesh
simpleFoam
```

### Method 3: Start from Tutorials

```bash
cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily myCase
cd myCase
# Modify for your application
```

## Common Workflow

1. **Create Case Structure**
   ```bash
   mkdir -p myCase/{0,constant,system}
   ```

2. **Define Geometry** (system/blockMeshDict or STL)

3. **Set Boundary Conditions** (0/U, 0/p, etc.)

4. **Define Physics** (constant/transportProperties, turbulenceProperties)

5. **Configure Numerics** (system/fvSchemes, fvSolution)

6. **Set Simulation Control** (system/controlDict)

7. **Generate Mesh**
   ```bash
   blockMesh  # or snappyHexMesh
   checkMesh
   ```

8. **Initialize** (optional)
   ```bash
   potentialFoam  # Better initial velocity field
   ```

9. **Run Solver**
   ```bash
   simpleFoam > log.simpleFoam 2>&1 &
   tail -f log.simpleFoam
   ```

10. **Post-Process**
    ```bash
    paraFoam
    ```

## Parallelization

To run in parallel:

1. **Create** `system/decomposeParDict`:
   ```cpp
   numberOfSubdomains 4;
   method          scotch;
   ```

2. **Decompose**:
   ```bash
   decomposePar
   ```

3. **Run**:
   ```bash
   mpirun -np 4 simpleFoam -parallel
   ```

4. **Reconstruct**:
   ```bash
   reconstructPar -latestTime
   ```

## Troubleshooting

### Mesh Quality Issues

```bash
checkMesh
```

Look for:
- Non-orthogonality > 70
- Skewness > 4
- Negative volumes

Fix by refining mesh or adjusting `nNonOrthogonalCorrectors` in `system/fvSolution`.

### Divergence

- Reduce relaxation factors in `system/fvSolution`
- Use more diffusive schemes in `system/fvSchemes` (e.g., upwind)
- Improve initial conditions
- Reduce time step (for transient)

### Slow Convergence

- Check boundary conditions
- Use better initial conditions (`potentialFoam`)
- Adjust solver tolerances
- Use GAMG for pressure

### Poor Results

- Perform mesh independence study
- Check y+ values (should be 30-300 for wall functions, <1 for resolved)
- Validate against experimental data
- Check residuals (should be < 1e-4 for engineering accuracy)

## Validation

Always validate your results:

1. **Grid Independence**: Run with 2x and 4x mesh resolution
2. **Time Step Independence**: Run with 0.5x and 0.25x time step (transient)
3. **Analytical Comparison**: Compare with known solutions where available
4. **Experimental Data**: Validate against measurements
5. **Code Verification**: Use method of manufactured solutions

## Tips and Best Practices

1. **Start Simple**: Begin with coarse mesh and simple schemes
2. **Monitor Residuals**: Should decrease steadily
3. **Check Courant Number**: Keep Co < 1 for transient (adjustable time stepping recommended)
4. **Use Wall Functions**: Unless y+ < 1 is required
5. **Version Control**: Use git for your cases
6. **Document**: Keep notes on settings and results
7. **Backup**: Save converged solutions
8. **Visualize Early**: Check mesh and initial conditions in ParaView before running

## Resources

- Main skill documentation: [../SKILL.md](../SKILL.md)
- Detailed reference: [../reference.md](../reference.md)
- Boundary conditions: [boundary-conditions-examples.md](boundary-conditions-examples.md)

## Getting Help

- OpenFOAM Forums: https://www.cfd-online.com/Forums/openfoam/
- Documentation: https://www.openfoam.org/
- Issue Tracker: https://bugs.openfoam.org/

## Contributing

To add new examples:

1. Create a complete, working case
2. Include a detailed README.md
3. Test on a clean OpenFOAM installation
4. Document expected runtime and results
5. Include post-processing instructions

## License

These examples are provided for educational purposes. OpenFOAM is licensed under GPLv3.
