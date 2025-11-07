#!/usr/bin/env python3
"""
Pump Selection Helper - Interactive pump type selector
Helps engineers choose the right pump based on operating conditions
"""

import math
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


@dataclass
class PumpRecommendation:
    """Data class for pump recommendations"""
    pump_type: str
    suitability_score: float  # 0-100
    reasons: List[str]
    warnings: List[str]
    typical_efficiency: str
    cost_range: str
    maintenance: str


class PumpSelector:
    """Main pump selection logic"""

    def __init__(self):
        self.recommendations = []

    @staticmethod
    def calculate_specific_speed(flow_gpm: float, head_ft: float, speed_rpm: float = 1750) -> float:
        """
        Calculate specific speed (US units)
        Ns = N * sqrt(Q) / H^0.75

        Args:
            flow_gpm: Flow rate in gallons per minute
            head_ft: Total head in feet
            speed_rpm: Rotational speed in RPM (default 1750)

        Returns:
            Specific speed (dimensionless)
        """
        if head_ft <= 0:
            raise ValueError("Head must be positive")

        ns = speed_rpm * math.sqrt(flow_gpm) / (head_ft ** 0.75)
        return ns

    @staticmethod
    def calculate_specific_speed_si(flow_m3h: float, head_m: float, speed_rpm: float = 1750) -> float:
        """
        Calculate specific speed (SI units)

        Args:
            flow_m3h: Flow rate in m³/h
            head_m: Total head in meters
            speed_rpm: Rotational speed in RPM

        Returns:
            Specific speed (SI)
        """
        if head_m <= 0:
            raise ValueError("Head must be positive")

        ns = speed_rpm * math.sqrt(flow_m3h) / (head_m ** 0.75)
        return ns

    @staticmethod
    def calculate_brake_horsepower(flow_gpm: float, head_ft: float,
                                   efficiency: float = 0.75, sg: float = 1.0) -> float:
        """
        Calculate brake horsepower
        BHP = (Q * H * SG) / (3960 * Eff)

        Args:
            flow_gpm: Flow rate in GPM
            head_ft: Total head in feet
            efficiency: Pump efficiency (0-1)
            sg: Specific gravity

        Returns:
            Brake horsepower
        """
        bhp = (flow_gpm * head_ft * sg) / (3960 * efficiency)
        return bhp

    @staticmethod
    def check_npsh(atm_pressure_psi: float, static_head_ft: float,
                   vapor_pressure_psi: float, friction_loss_ft: float) -> float:
        """
        Calculate available NPSH
        NPSHa = (Patm / 2.31 / SG) + Static head - (Pvapor / 2.31 / SG) - Friction losses

        Args:
            atm_pressure_psi: Atmospheric pressure (typically 14.7 psi at sea level)
            static_head_ft: Static suction head (positive) or lift (negative)
            vapor_pressure_psi: Vapor pressure of fluid at operating temp
            friction_loss_ft: Friction losses in suction line

        Returns:
            Available NPSH in feet
        """
        npsh_a = (atm_pressure_psi * 2.31) + static_head_ft - (vapor_pressure_psi * 2.31) - friction_loss_ft
        return npsh_a

    def evaluate_centrifugal(self, flow_gpm: float, head_ft: float,
                            viscosity_cp: float, ns: float) -> Optional[PumpRecommendation]:
        """Evaluate suitability of centrifugal pumps"""

        reasons = []
        warnings = []
        score = 100

        # Check viscosity
        if viscosity_cp > 500:
            warnings.append(f"High viscosity ({viscosity_cp} cP) - efficiency will be reduced")
            score -= 30
        if viscosity_cp > 1000:
            warnings.append("Viscosity too high for centrifugal - consider positive displacement")
            return None

        # Determine centrifugal type based on specific speed
        if ns < 500:
            pump_subtype = "High-speed radial/turbine"
            reasons.append(f"Very low specific speed ({ns:.0f}) indicates high head, low flow")
        elif ns < 2000:
            pump_subtype = "Radial flow (volute or diffuser)"
            reasons.append(f"Low specific speed ({ns:.0f}) - good for high head applications")
        elif ns < 4000:
            pump_subtype = "Francis vane (radial to mixed)"
            reasons.append(f"Medium specific speed ({ns:.0f}) - balanced head/flow design")
        elif ns < 9000:
            pump_subtype = "Mixed flow"
            reasons.append(f"High specific speed ({ns:.0f}) - efficient for moderate head, high flow")
        else:
            pump_subtype = "Axial flow (propeller)"
            reasons.append(f"Very high specific speed ({ns:.0f}) - designed for high flow, low head")
            if head_ft > 50:
                warnings.append("Head may be too high for axial flow pump")
                score -= 20

        # Flow range check
        if flow_gpm < 10:
            warnings.append("Low flow - centrifugal may be inefficient")
            score -= 20
        elif flow_gpm > 100000:
            warnings.append("Very high flow - verify manufacturer capabilities")
            score -= 10

        # Head range check
        if head_ft < 10:
            warnings.append("Very low head - consider axial flow or other options")
            score -= 10

        reasons.append("Simple, reliable design with low maintenance")
        reasons.append("Can handle some solids and variations in flow")

        return PumpRecommendation(
            pump_type=f"Centrifugal - {pump_subtype}",
            suitability_score=max(score, 0),
            reasons=reasons,
            warnings=warnings,
            typical_efficiency="60-85%",
            cost_range="$ - $$",
            maintenance="Low - quarterly inspection, annual bearing service"
        )

    def evaluate_gear_pump(self, flow_gpm: float, head_ft: float,
                          viscosity_cp: float) -> Optional[PumpRecommendation]:
        """Evaluate suitability of gear pumps"""

        reasons = []
        warnings = []
        score = 100

        # Flow range check
        if flow_gpm > 1500:
            return None  # Out of typical range
        elif flow_gpm > 500:
            warnings.append("High flow for gear pump - verify manufacturer capabilities")
            score -= 20

        # Pressure calculation (approximate)
        pressure_psi = head_ft * 0.433  # Convert head to psi (water)

        if pressure_psi > 3000:
            warnings.append("Pressure exceeds typical gear pump rating")
            return None

        # Viscosity advantage
        if viscosity_cp > 100:
            reasons.append(f"Good for viscous fluids ({viscosity_cp} cP)")
            score += 10
        if viscosity_cp > 1000:
            reasons.append("Excellent choice for high viscosity")
            score += 20

        reasons.append("Self-priming capability")
        reasons.append("Constant flow regardless of pressure variations")
        reasons.append("Compact design")

        warnings.append("Cannot handle abrasives")
        warnings.append("Pulsating flow - may need dampener")

        if viscosity_cp < 10:
            warnings.append("Low viscosity may reduce volumetric efficiency")
            score -= 10

        return PumpRecommendation(
            pump_type="Gear Pump (positive displacement)",
            suitability_score=min(score, 100),
            reasons=reasons,
            warnings=warnings,
            typical_efficiency="70-85%",
            cost_range="$$",
            maintenance="Moderate - regular gear inspection and replacement"
        )

    def evaluate_piston_pump(self, flow_gpm: float, head_ft: float,
                            viscosity_cp: float) -> Optional[PumpRecommendation]:
        """Evaluate suitability of piston/plunger pumps"""

        reasons = []
        warnings = []
        score = 100

        pressure_psi = head_ft * 0.433

        # Best for high pressure
        if pressure_psi < 500:
            score -= 30
            warnings.append("Piston pumps are typically oversized for low pressure")
        else:
            reasons.append(f"Excellent for high pressure ({pressure_psi:.0f} psi)")
            score += 20

        if pressure_psi > 5000:
            reasons.append("One of few pump types capable of very high pressure")
            score += 30

        # Flow range
        if flow_gpm > 5000:
            return None
        elif flow_gpm > 1000:
            warnings.append("High flow - may require multiple pumps")
            score -= 10

        # Viscosity
        if viscosity_cp > 100000:
            warnings.append("Very high viscosity - verify with manufacturer")
            score -= 10

        reasons.append("Accurate metering and flow control")
        reasons.append("Can achieve very high pressures")

        warnings.append("Pulsating flow - may need pulsation dampener")
        warnings.append("Higher maintenance than centrifugal")
        warnings.append("Higher initial cost")

        return PumpRecommendation(
            pump_type="Piston/Plunger Pump (positive displacement)",
            suitability_score=min(score, 100),
            reasons=reasons,
            warnings=warnings,
            typical_efficiency="75-90%",
            cost_range="$$$ - $$$$",
            maintenance="High - regular packing/seal replacement, valve service"
        )

    def evaluate_diaphragm_pump(self, flow_gpm: float, head_ft: float,
                               is_corrosive: bool = False,
                               is_abrasive: bool = False) -> Optional[PumpRecommendation]:
        """Evaluate suitability of diaphragm pumps"""

        reasons = []
        warnings = []
        score = 80  # Start lower as specialty pump

        pressure_psi = head_ft * 0.433

        # Flow range
        if flow_gpm > 800:
            return None

        # Pressure range
        if pressure_psi > 1000:
            warnings.append("Pressure exceeds typical diaphragm pump rating")
            return None

        # Special features
        if is_corrosive:
            reasons.append("Seal-less design - excellent for corrosive fluids")
            score += 20

        if is_abrasive:
            reasons.append("Can handle abrasive slurries")
            score += 15

        reasons.append("No seals - zero leakage of hazardous fluids")
        reasons.append("Can run dry without damage")
        reasons.append("Self-priming")

        warnings.append("Pulsating flow")
        warnings.append("Diaphragm requires periodic replacement")

        return PumpRecommendation(
            pump_type="Diaphragm Pump (positive displacement)",
            suitability_score=min(score, 100),
            reasons=reasons,
            warnings=warnings,
            typical_efficiency="70-80%",
            cost_range="$$ - $$$",
            maintenance="Moderate - diaphragm and valve replacement"
        )

    def evaluate_screw_pump(self, flow_gpm: float, head_ft: float,
                           viscosity_cp: float,
                           is_shear_sensitive: bool = False) -> Optional[PumpRecommendation]:
        """Evaluate suitability of progressive cavity (screw) pumps"""

        reasons = []
        warnings = []
        score = 90

        pressure_psi = head_ft * 0.433

        # Flow range
        if flow_gpm > 2000:
            return None

        # Pressure range
        if pressure_psi > 1500:
            warnings.append("Pressure high for screw pump")
            score -= 20

        # Viscosity - major advantage
        if viscosity_cp > 1000:
            reasons.append(f"Excellent for high viscosity ({viscosity_cp} cP)")
            score += 10
        if viscosity_cp > 10000:
            reasons.append("One of best choices for very viscous fluids")
            score += 20

        # Shear sensitive
        if is_shear_sensitive:
            reasons.append("Low shear - gentle pumping action")
            score += 15

        reasons.append("Smooth, non-pulsating flow")
        reasons.append("Self-priming")
        reasons.append("Can handle solids up to 50mm")

        warnings.append("Rotor/stator wear - requires replacement")
        warnings.append("Speed limited to prevent stator heating")

        return PumpRecommendation(
            pump_type="Progressive Cavity/Screw Pump (positive displacement)",
            suitability_score=min(score, 100),
            reasons=reasons,
            warnings=warnings,
            typical_efficiency="75-85%",
            cost_range="$$ - $$$",
            maintenance="Moderate - rotor/stator replacement every 1-3 years"
        )

    def select_pump(self, flow_gpm: float, head_ft: float,
                   viscosity_cp: float = 1.0,
                   fluid_type: str = "water",
                   speed_rpm: float = 1750,
                   is_corrosive: bool = False,
                   is_abrasive: bool = False,
                   is_shear_sensitive: bool = False,
                   requires_constant_flow: bool = False) -> List[PumpRecommendation]:
        """
        Main pump selection function

        Args:
            flow_gpm: Flow rate in gallons per minute
            head_ft: Total head in feet
            viscosity_cp: Fluid viscosity in centipoise (default 1.0 for water)
            fluid_type: Type of fluid (water, oil, chemical, slurry, etc.)
            speed_rpm: Expected pump speed in RPM
            is_corrosive: Whether fluid is corrosive
            is_abrasive: Whether fluid contains abrasives
            is_shear_sensitive: Whether fluid is shear-sensitive
            requires_constant_flow: Whether constant flow is required despite pressure changes

        Returns:
            List of PumpRecommendation objects, sorted by suitability score
        """

        recommendations = []

        # Calculate specific speed for centrifugal evaluation
        ns = self.calculate_specific_speed(flow_gpm, head_ft, speed_rpm)

        # Evaluate centrifugal (unless viscosity too high)
        if viscosity_cp < 1000:
            centrifugal = self.evaluate_centrifugal(flow_gpm, head_ft, viscosity_cp, ns)
            if centrifugal:
                # Penalize if constant flow required
                if requires_constant_flow:
                    centrifugal.warnings.append("Flow varies with pressure - not ideal for constant flow requirement")
                    centrifugal.suitability_score -= 30
                recommendations.append(centrifugal)

        # Evaluate gear pump
        gear = self.evaluate_gear_pump(flow_gpm, head_ft, viscosity_cp)
        if gear:
            if requires_constant_flow:
                gear.reasons.append("Provides constant flow regardless of pressure")
                gear.suitability_score += 10
            recommendations.append(gear)

        # Evaluate piston pump
        piston = self.evaluate_piston_pump(flow_gpm, head_ft, viscosity_cp)
        if piston:
            if requires_constant_flow:
                piston.reasons.append("Positive displacement provides consistent flow")
                piston.suitability_score += 5
            recommendations.append(piston)

        # Evaluate diaphragm pump
        diaphragm = self.evaluate_diaphragm_pump(flow_gpm, head_ft, is_corrosive, is_abrasive)
        if diaphragm:
            recommendations.append(diaphragm)

        # Evaluate screw pump
        screw = self.evaluate_screw_pump(flow_gpm, head_ft, viscosity_cp, is_shear_sensitive)
        if screw:
            if requires_constant_flow:
                screw.reasons.append("Non-pulsating constant flow")
                screw.suitability_score += 10
            recommendations.append(screw)

        # Sort by suitability score
        recommendations.sort(key=lambda x: x.suitability_score, reverse=True)

        return recommendations


def format_recommendations(recommendations: List[PumpRecommendation],
                          flow_gpm: float, head_ft: float,
                          viscosity_cp: float, ns: float) -> str:
    """Format recommendations for display"""

    output = []
    output.append("=" * 80)
    output.append("PUMP SELECTION RESULTS")
    output.append("=" * 80)
    output.append(f"\nInput Parameters:")
    output.append(f"  Flow rate:      {flow_gpm:.1f} gpm")
    output.append(f"  Total head:     {head_ft:.1f} ft")
    output.append(f"  Viscosity:      {viscosity_cp:.1f} cP")
    output.append(f"  Specific speed: {ns:.0f}")
    output.append(f"\n{'─' * 80}\n")

    if not recommendations:
        output.append("No suitable pumps found for the given parameters.")
        return "\n".join(output)

    for i, rec in enumerate(recommendations, 1):
        output.append(f"RECOMMENDATION #{i}: {rec.pump_type}")
        output.append(f"Suitability Score: {rec.suitability_score:.0f}/100")
        output.append(f"")

        output.append("Reasons for selection:")
        for reason in rec.reasons:
            output.append(f"  ✓ {reason}")

        if rec.warnings:
            output.append("\nWarnings/Considerations:")
            for warning in rec.warnings:
                output.append(f"  ⚠ {warning}")

        output.append(f"\nPerformance:")
        output.append(f"  Typical efficiency: {rec.typical_efficiency}")
        output.append(f"  Cost range:         {rec.cost_range}")
        output.append(f"  Maintenance:        {rec.maintenance}")

        if i < len(recommendations):
            output.append(f"\n{'─' * 80}\n")

    output.append(f"\n{'=' * 80}")

    return "\n".join(output)


def interactive_mode():
    """Interactive pump selection"""

    print("\n" + "=" * 80)
    print("PUMP SELECTION HELPER - Interactive Mode")
    print("=" * 80 + "\n")

    try:
        # Get basic parameters
        flow_gpm = float(input("Enter flow rate (gpm): "))
        head_ft = float(input("Enter total head (ft): "))
        viscosity_cp = float(input("Enter fluid viscosity (cP, default=1 for water): ") or "1")
        speed_rpm = float(input("Enter pump speed (rpm, default=1750): ") or "1750")

        # Get fluid properties
        print("\nFluid properties (y/n):")
        is_corrosive = input("  Is fluid corrosive? (y/n): ").lower().startswith('y')
        is_abrasive = input("  Does fluid contain abrasives? (y/n): ").lower().startswith('y')
        is_shear_sensitive = input("  Is fluid shear-sensitive? (y/n): ").lower().startswith('y')
        requires_constant_flow = input("  Requires constant flow? (y/n): ").lower().startswith('y')

        # Calculate specific speed
        selector = PumpSelector()
        ns = selector.calculate_specific_speed(flow_gpm, head_ft, speed_rpm)

        # Get recommendations
        recommendations = selector.select_pump(
            flow_gpm=flow_gpm,
            head_ft=head_ft,
            viscosity_cp=viscosity_cp,
            speed_rpm=speed_rpm,
            is_corrosive=is_corrosive,
            is_abrasive=is_abrasive,
            is_shear_sensitive=is_shear_sensitive,
            requires_constant_flow=requires_constant_flow
        )

        # Display results
        print("\n" + format_recommendations(recommendations, flow_gpm, head_ft, viscosity_cp, ns))

        # Calculate power requirement
        if recommendations:
            best = recommendations[0]
            eff = 0.75  # Default efficiency
            bhp = selector.calculate_brake_horsepower(flow_gpm, head_ft, eff)
            print(f"\nEstimated brake horsepower (at 75% efficiency): {bhp:.2f} HP")
            print(f"Estimated motor size: {math.ceil(bhp * 1.15)} HP (with 15% service factor)")

    except ValueError as e:
        print(f"\nError: Invalid input - {e}")
    except KeyboardInterrupt:
        print("\n\nSelection cancelled.")


def main():
    """Main entry point"""

    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "--example":
        # Run examples
        print("\nEXAMPLE 1: Water supply pump")
        print("-" * 40)
        selector = PumpSelector()
        recs = selector.select_pump(flow_gpm=500, head_ft=200, viscosity_cp=1)
        ns = selector.calculate_specific_speed(500, 200)
        print(format_recommendations(recs, 500, 200, 1, ns))

        print("\n\nEXAMPLE 2: High viscosity oil transfer")
        print("-" * 40)
        recs = selector.select_pump(flow_gpm=50, head_ft=100, viscosity_cp=5000)
        ns = selector.calculate_specific_speed(50, 100)
        print(format_recommendations(recs, 50, 100, 5000, ns))

        print("\n\nEXAMPLE 3: Corrosive chemical metering")
        print("-" * 40)
        recs = selector.select_pump(flow_gpm=10, head_ft=150, viscosity_cp=10,
                                    is_corrosive=True, requires_constant_flow=True)
        ns = selector.calculate_specific_speed(10, 150)
        print(format_recommendations(recs, 10, 150, 10, ns))
    else:
        # Interactive mode
        interactive_mode()


if __name__ == "__main__":
    main()
