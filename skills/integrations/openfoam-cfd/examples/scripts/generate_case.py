#!/usr/bin/env python3
"""
OpenFOAM Case Generator
Automates the creation of OpenFOAM case directories with templated configurations.
"""

import os
import shutil
import argparse
from pathlib import Path


class OpenFOAMCaseGenerator:
    """Generate OpenFOAM case directories with standard configurations."""

    def __init__(self, case_name, case_dir="."):
        self.case_name = case_name
        self.case_dir = Path(case_dir) / case_name

    def create_directory_structure(self):
        """Create the basic OpenFOAM case directory structure."""
        dirs = [
            self.case_dir / "0",
            self.case_dir / "constant" / "polyMesh",
            self.case_dir / "system",
        ]

        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)

        print(f"Created case directory: {self.case_dir}")

    def write_control_dict(self, solver="simpleFoam", end_time=1000,
                          write_interval=100, delta_t=1):
        """Write system/controlDict file."""
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

application     {solver};

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time};
deltaT          {delta_t};

writeControl    timeStep;
writeInterval   {write_interval};
purgeWrite      2;

writeFormat     ascii;
writePrecision  6;
writeCompression off;

timeFormat      general;
timePrecision   6;

runTimeModifiable true;

// ************************************************************************* //
"""

        filepath = self.case_dir / "system" / "controlDict"
        filepath.write_text(content)
        print(f"Created: {filepath}")

    def write_fv_schemes(self, steady_state=True):
        """Write system/fvSchemes file."""
        ddt_scheme = "steadyState" if steady_state else "Euler"

        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{{
    default         {ddt_scheme};
}}

gradSchemes
{{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}}

divSchemes
{{
    default         none;
    div(phi,U)      bounded Gauss linearUpwind grad(U);
    div(phi,k)      bounded Gauss upwind;
    div(phi,epsilon) bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}

laplacianSchemes
{{
    default         Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         corrected;
}}

wallDist
{{
    method          meshWave;
}}

// ************************************************************************* //
"""

        filepath = self.case_dir / "system" / "fvSchemes"
        filepath.write_text(content)
        print(f"Created: {filepath}")

    def write_fv_solution(self, algorithm="SIMPLE"):
        """Write system/fvSolution file."""
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{{
    p
    {{
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.01;
        smoother        GaussSeidel;
        nPreSweeps      0;
        nPostSweeps     2;
        cacheAgglomeration true;
        nCellsInCoarsestLevel 10;
        agglomerator    faceAreaPair;
        mergeLevels     1;
    }}

    U
    {{
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-07;
        relTol          0.1;
        nSweeps         1;
    }}

    "(k|epsilon|omega)"
    {{
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-08;
        relTol          0.1;
        nSweeps         1;
    }}
}}

{algorithm}
{{
    nNonOrthogonalCorrectors 0;
    consistent      yes;

    residualControl
    {{
        p               1e-5;
        U               1e-5;
        "(k|epsilon|omega)" 1e-5;
    }}
}}

relaxationFactors
{{
    fields
    {{
        p               0.3;
    }}
    equations
    {{
        U               0.7;
        k               0.7;
        epsilon         0.7;
        omega           0.7;
    }}
}}

// ************************************************************************* //
"""

        filepath = self.case_dir / "system" / "fvSolution"
        filepath.write_text(content)
        print(f"Created: {filepath}")

    def write_transport_properties(self, nu=1e-5):
        """Write constant/transportProperties file."""
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

transportModel  Newtonian;

nu              [0 2 -1 0 0 0 0] {nu};

// ************************************************************************* //
"""

        filepath = self.case_dir / "constant" / "transportProperties"
        filepath.write_text(content)
        print(f"Created: {filepath}")

    def write_turbulence_properties(self, model="kEpsilon"):
        """Write constant/turbulenceProperties file."""
        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  RAS;

RAS
{{
    RASModel        {model};

    turbulence      on;

    printCoeffs     on;
}}

// ************************************************************************* //
"""

        filepath = self.case_dir / "constant" / "turbulenceProperties"
        filepath.write_text(content)
        print(f"Created: {filepath}")

    def write_boundary_field(self, field_name, field_class, dimensions,
                            internal_value, boundary_conditions):
        """Write a boundary field file in the 0 directory."""
        bc_text = ""
        for patch_name, bc_dict in boundary_conditions.items():
            bc_text += f"    {patch_name}\n    {{\n"
            for key, value in bc_dict.items():
                bc_text += f"        {key}            {value};\n"
            bc_text += "    }\n\n"

        content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  11                                    |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {field_class};
    location    "0";
    object      {field_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      {dimensions};

internalField   uniform {internal_value};

boundaryField
{{
{bc_text}}}

// ************************************************************************* //
"""

        filepath = self.case_dir / "0" / field_name
        filepath.write_text(content)
        print(f"Created: {filepath}")

    def generate_standard_case(self, velocity=(1, 0, 0), nu=1e-5,
                               turbulence_model="kEpsilon"):
        """Generate a standard incompressible turbulent flow case."""
        print(f"\nGenerating OpenFOAM case: {self.case_name}")
        print("=" * 60)

        # Create directories
        self.create_directory_structure()

        # Write system files
        self.write_control_dict()
        self.write_fv_schemes(steady_state=True)
        self.write_fv_solution(algorithm="SIMPLE")

        # Write constant files
        self.write_transport_properties(nu=nu)
        self.write_turbulence_properties(model=turbulence_model)

        # Write 0 directory fields
        u_bcs = {
            "inlet": {
                "type": "fixedValue",
                "value": f"uniform ({velocity[0]} {velocity[1]} {velocity[2]})"
            },
            "outlet": {
                "type": "zeroGradient"
            },
            "walls": {
                "type": "noSlip"
            }
        }
        self.write_boundary_field("U", "volVectorField", "[0 1 -1 0 0 0 0]",
                                 "(0 0 0)", u_bcs)

        p_bcs = {
            "inlet": {
                "type": "zeroGradient"
            },
            "outlet": {
                "type": "fixedValue",
                "value": "uniform 0"
            },
            "walls": {
                "type": "zeroGradient"
            }
        }
        self.write_boundary_field("p", "volScalarField", "[0 2 -2 0 0 0 0]",
                                 "0", p_bcs)

        k_bcs = {
            "inlet": {
                "type": "turbulentIntensityKineticEnergyInlet",
                "intensity": "0.05",
                "value": "uniform 0.0015"
            },
            "outlet": {
                "type": "zeroGradient"
            },
            "walls": {
                "type": "kqRWallFunction",
                "value": "uniform 0.0015"
            }
        }
        self.write_boundary_field("k", "volScalarField", "[0 2 -2 0 0 0 0]",
                                 "0.0015", k_bcs)

        if turbulence_model in ["kEpsilon", "realizableKE"]:
            epsilon_bcs = {
                "inlet": {
                    "type": "turbulentMixingLengthDissipationRateInlet",
                    "mixingLength": "0.005",
                    "value": "uniform 0.01"
                },
                "outlet": {
                    "type": "zeroGradient"
                },
                "walls": {
                    "type": "epsilonWallFunction",
                    "value": "uniform 0.01"
                }
            }
            self.write_boundary_field("epsilon", "volScalarField",
                                     "[0 2 -3 0 0 0 0]", "0.01", epsilon_bcs)

        if turbulence_model in ["kOmegaSST", "kOmega"]:
            omega_bcs = {
                "inlet": {
                    "type": "turbulentMixingLengthFrequencyInlet",
                    "mixingLength": "0.005",
                    "value": "uniform 1"
                },
                "outlet": {
                    "type": "zeroGradient"
                },
                "walls": {
                    "type": "omegaWallFunction",
                    "value": "uniform 1"
                }
            }
            self.write_boundary_field("omega", "volScalarField",
                                     "[0 0 -1 0 0 0 0]", "1", omega_bcs)

        nut_bcs = {
            "inlet": {
                "type": "calculated",
                "value": "uniform 0"
            },
            "outlet": {
                "type": "calculated",
                "value": "uniform 0"
            },
            "walls": {
                "type": "nutkWallFunction",
                "value": "uniform 0"
            }
        }
        self.write_boundary_field("nut", "volScalarField", "[0 2 -1 0 0 0 0]",
                                 "0", nut_bcs)

        print("=" * 60)
        print(f"Case generation complete!\n")
        print(f"Next steps:")
        print(f"  1. cd {self.case_dir}")
        print(f"  2. Create/edit system/blockMeshDict for your geometry")
        print(f"  3. Run: blockMesh")
        print(f"  4. Run: simpleFoam")
        print(f"  5. Post-process: paraFoam")


def main():
    parser = argparse.ArgumentParser(
        description="Generate OpenFOAM case directories automatically"
    )
    parser.add_argument("case_name", help="Name of the case directory")
    parser.add_argument("--dir", default=".", help="Parent directory for case")
    parser.add_argument("--velocity", nargs=3, type=float, default=[1, 0, 0],
                       help="Inlet velocity vector (Ux Uy Uz)")
    parser.add_argument("--nu", type=float, default=1e-5,
                       help="Kinematic viscosity (m²/s)")
    parser.add_argument("--turbulence", default="kEpsilon",
                       choices=["kEpsilon", "kOmegaSST", "realizableKE"],
                       help="Turbulence model")

    args = parser.parse_args()

    generator = OpenFOAMCaseGenerator(args.case_name, args.dir)
    generator.generate_standard_case(
        velocity=tuple(args.velocity),
        nu=args.nu,
        turbulence_model=args.turbulence
    )


if __name__ == "__main__":
    main()
