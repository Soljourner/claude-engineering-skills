#!/usr/bin/env python3
"""
Turbulence Parameter Calculator for OpenFOAM
Calculates k, epsilon, omega, and other turbulence parameters from flow conditions.
"""

import math
import argparse


class TurbulenceCalculator:
    """Calculate turbulence parameters for OpenFOAM boundary conditions."""

    def __init__(self, velocity, length_scale, nu=1e-5, intensity=None):
        """
        Initialize calculator.

        Args:
            velocity: Mean velocity magnitude (m/s)
            length_scale: Characteristic length scale (m)
            nu: Kinematic viscosity (m²/s)
            intensity: Turbulence intensity (0-1), if None will be calculated
        """
        self.U = velocity
        self.L = length_scale
        self.nu = nu
        self.Cmu = 0.09  # Standard k-epsilon constant

        # Calculate Reynolds number
        self.Re = self.U * self.L / self.nu

        # Calculate or use provided turbulence intensity
        if intensity is None:
            self.I = self.calculate_turbulence_intensity()
        else:
            self.I = intensity

    def calculate_turbulence_intensity(self):
        """
        Calculate turbulence intensity from Reynolds number.
        Uses empirical correlation: I = 0.16 * Re^(-1/8)
        """
        if self.Re < 1:
            return 0.01  # Minimum for numerical stability

        I = 0.16 * (self.Re ** (-1.0/8.0))
        return I

    def calculate_k(self):
        """
        Calculate turbulent kinetic energy.
        k = 3/2 * (U * I)²
        """
        k = 1.5 * (self.U * self.I) ** 2
        return k

    def calculate_epsilon(self, length_scale=None):
        """
        Calculate turbulent dissipation rate.
        ε = Cμ^(3/4) * k^(3/2) / L
        """
        if length_scale is None:
            # Use default turbulent length scale (7% of characteristic length)
            L_turb = 0.07 * self.L
        else:
            L_turb = length_scale

        k = self.calculate_k()
        epsilon = (self.Cmu ** 0.75) * (k ** 1.5) / L_turb
        return epsilon

    def calculate_omega(self, length_scale=None):
        """
        Calculate specific dissipation rate.
        ω = k^(1/2) / (Cμ^(1/4) * L)
        """
        if length_scale is None:
            L_turb = 0.07 * self.L
        else:
            L_turb = length_scale

        k = self.calculate_k()
        omega = (k ** 0.5) / ((self.Cmu ** 0.25) * L_turb)
        return omega

    def calculate_nut(self):
        """
        Calculate turbulent viscosity.
        νt = Cμ * k² / ε
        """
        k = self.calculate_k()
        epsilon = self.calculate_epsilon()
        nut = self.Cmu * (k ** 2) / epsilon
        return nut

    def calculate_mixing_length(self):
        """
        Calculate turbulent mixing length.
        L = k^(3/2) / ε
        """
        k = self.calculate_k()
        epsilon = self.calculate_epsilon()
        L_mix = (k ** 1.5) / epsilon
        return L_mix

    def print_summary(self):
        """Print a formatted summary of all turbulence parameters."""
        print("\n" + "=" * 70)
        print("TURBULENCE PARAMETER CALCULATOR FOR OPENFOAM")
        print("=" * 70)
        print("\nInput Parameters:")
        print(f"  Velocity (U):              {self.U:12.6f} m/s")
        print(f"  Length Scale (L):          {self.L:12.6f} m")
        print(f"  Kinematic Viscosity (ν):   {self.nu:12.2e} m²/s")
        print(f"  Reynolds Number (Re):      {self.Re:12.2f}")

        print("\nTurbulence Parameters:")
        print(f"  Turbulence Intensity (I):  {self.I:12.6f} ({self.I*100:.2f}%)")
        print(f"  Turb. Kinetic Energy (k):  {self.calculate_k():12.6e} m²/s²")
        print(f"  Dissipation Rate (ε):      {self.calculate_epsilon():12.6e} m²/s³")
        print(f"  Spec. Dissipation (ω):     {self.calculate_omega():12.6e} 1/s")
        print(f"  Turbulent Viscosity (νt):  {self.calculate_nut():12.6e} m²/s")
        print(f"  Mixing Length (L_mix):     {self.calculate_mixing_length():12.6e} m")
        print(f"  νt/ν Ratio:                {self.calculate_nut()/self.nu:12.2f}")

        print("\n" + "=" * 70)
        print("OPENFOAM BOUNDARY CONDITIONS")
        print("=" * 70)

        k = self.calculate_k()
        epsilon = self.calculate_epsilon()
        omega = self.calculate_omega()
        L_turb = 0.07 * self.L

        print("\nFor k-epsilon model:")
        print("-" * 70)
        print("0/k:")
        print(f"""
    inlet
    {{
        type            turbulentIntensityKineticEnergyInlet;
        intensity       {self.I:.6f};
        value           uniform {k:.6e};
    }}
""")

        print("0/epsilon:")
        print(f"""
    inlet
    {{
        type            turbulentMixingLengthDissipationRateInlet;
        mixingLength    {L_turb:.6f};
        value           uniform {epsilon:.6e};
    }}
""")

        print("\nFor k-omega SST model:")
        print("-" * 70)
        print("0/k:")
        print(f"""
    inlet
    {{
        type            turbulentIntensityKineticEnergyInlet;
        intensity       {self.I:.6f};
        value           uniform {k:.6e};
    }}
""")

        print("0/omega:")
        print(f"""
    inlet
    {{
        type            turbulentMixingLengthFrequencyInlet;
        mixingLength    {L_turb:.6f};
        value           uniform {omega:.6e};
    }}
""")

        print("\nAlternative (fixed values):")
        print("-" * 70)
        print(f"internalField   uniform {k:.6e};  // for 0/k")
        print(f"internalField   uniform {epsilon:.6e};  // for 0/epsilon")
        print(f"internalField   uniform {omega:.6e};  // for 0/omega")

        print("\n" + "=" * 70)

    def generate_boundary_condition_snippet(self, model="kEpsilon"):
        """Generate OpenFOAM boundary condition code snippet."""
        k = self.calculate_k()
        L_turb = 0.07 * self.L

        if model == "kEpsilon":
            epsilon = self.calculate_epsilon()
            return f"""
// k-epsilon boundary conditions
// 0/k
inlet
{{
    type            turbulentIntensityKineticEnergyInlet;
    intensity       {self.I:.6f};
    value           uniform {k:.6e};
}}

// 0/epsilon
inlet
{{
    type            turbulentMixingLengthDissipationRateInlet;
    mixingLength    {L_turb:.6f};
    value           uniform {epsilon:.6e};
}}
"""
        elif model == "kOmegaSST":
            omega = self.calculate_omega()
            return f"""
// k-omega SST boundary conditions
// 0/k
inlet
{{
    type            turbulentIntensityKineticEnergyInlet;
    intensity       {self.I:.6f};
    value           uniform {k:.6e};
}}

// 0/omega
inlet
{{
    type            turbulentMixingLengthFrequencyInlet;
    mixingLength    {L_turb:.6f};
    value           uniform {omega:.6e};
}}
"""


def main():
    parser = argparse.ArgumentParser(
        description="Calculate turbulence parameters for OpenFOAM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Pipe flow: D=0.1m, U=1m/s, air
  %(prog)s --velocity 1 --length 0.1 --nu 1.5e-5

  # Channel flow: H=0.05m, U=2m/s, water
  %(prog)s --velocity 2 --length 0.05 --nu 1e-6

  # With specified turbulence intensity
  %(prog)s --velocity 10 --length 1 --intensity 0.05
        """
    )

    parser.add_argument("--velocity", "-U", type=float, required=True,
                       help="Mean velocity magnitude (m/s)")
    parser.add_argument("--length", "-L", type=float, required=True,
                       help="Characteristic length scale (m)")
    parser.add_argument("--nu", type=float, default=1.5e-5,
                       help="Kinematic viscosity (m²/s), default: 1.5e-5 (air)")
    parser.add_argument("--intensity", "-I", type=float, default=None,
                       help="Turbulence intensity (0-1), if not provided will be calculated")
    parser.add_argument("--model", choices=["kEpsilon", "kOmegaSST"],
                       default="kEpsilon", help="Turbulence model")
    parser.add_argument("--brief", action="store_true",
                       help="Show brief output (values only)")

    args = parser.parse_args()

    calc = TurbulenceCalculator(
        velocity=args.velocity,
        length_scale=args.length,
        nu=args.nu,
        intensity=args.intensity
    )

    if args.brief:
        print(f"k={calc.calculate_k():.6e}")
        if args.model == "kEpsilon":
            print(f"epsilon={calc.calculate_epsilon():.6e}")
        else:
            print(f"omega={calc.calculate_omega():.6e}")
    else:
        calc.print_summary()


if __name__ == "__main__":
    main()
