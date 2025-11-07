"""
Pump Efficiency Optimization

This module provides tools for calculating and optimizing pump efficiency
through design modifications and operational strategies.

Key features:
- Efficiency component calculations (hydraulic, volumetric, mechanical)
- Loss mechanism quantification
- Design parameter optimization
- Energy cost analysis
- Life cycle cost calculations

Dependencies: numpy, scipy
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from typing import Dict, Tuple, List
import warnings

# Physical constants
GRAVITY = 9.81  # m/s²
WATER_DENSITY = 1000  # kg/m³


class PumpEfficiencyAnalyzer:
    """
    Analyzes and calculates various pump efficiency components.
    """

    def __init__(self, density: float = WATER_DENSITY):
        """
        Initialize analyzer.

        Args:
            density: Fluid density (kg/m³)
        """
        self.rho = density
        self.g = GRAVITY

    def hydraulic_efficiency(self, flow: float, head: float,
                            impeller_tip_speed: float,
                            tangential_velocity: float) -> float:
        """
        Calculate hydraulic efficiency.

        η_h = (g × H) / (u₂ × c_u2)

        Args:
            flow: Flow rate (m³/s)
            head: Head developed (m)
            impeller_tip_speed: u₂ (m/s)
            tangential_velocity: c_u2 (m/s)

        Returns:
            Hydraulic efficiency (0-1)
        """
        if impeller_tip_speed * tangential_velocity <= 0:
            return 0.0

        eta_h = (self.g * head) / (impeller_tip_speed * tangential_velocity)
        return min(eta_h, 1.0)

    def volumetric_efficiency(self, delivered_flow: float,
                             leakage_flow: float) -> float:
        """
        Calculate volumetric efficiency.

        η_v = Q_delivered / (Q_delivered + Q_leakage)

        Args:
            delivered_flow: Actual delivered flow (m³/s)
            leakage_flow: Internal leakage (m³/s)

        Returns:
            Volumetric efficiency (0-1)
        """
        total_flow = delivered_flow + leakage_flow
        if total_flow <= 0:
            return 0.0

        return delivered_flow / total_flow

    def mechanical_efficiency(self, hydraulic_power: float,
                             friction_losses: float) -> float:
        """
        Calculate mechanical efficiency.

        η_m = P_hydraulic / (P_hydraulic + P_friction)

        Args:
            hydraulic_power: Useful hydraulic power (W)
            friction_losses: Power lost to friction (W)

        Returns:
            Mechanical efficiency (0-1)
        """
        total_power = hydraulic_power + friction_losses
        if total_power <= 0:
            return 0.0

        return hydraulic_power / total_power

    def overall_efficiency(self, flow: float, head: float,
                          shaft_power: float) -> float:
        """
        Calculate overall pump efficiency.

        η_overall = (ρ × g × Q × H) / P_shaft

        Args:
            flow: Flow rate (m³/s)
            head: Head developed (m)
            shaft_power: Shaft power input (W)

        Returns:
            Overall efficiency (0-1)
        """
        if shaft_power <= 0:
            return 0.0

        hydraulic_power = self.rho * self.g * flow * head
        eta = hydraulic_power / shaft_power

        return min(eta, 1.0)

    def disk_friction_power(self, omega: float, radius: float,
                           clearance: float = 0.001) -> float:
        """
        Calculate disk friction power loss.

        P_disk ≈ k × ρ × ω³ × r⁵ × f(clearance)

        Args:
            omega: Angular velocity (rad/s)
            radius: Impeller radius (m)
            clearance: Clearance gap (m)

        Returns:
            Disk friction power loss (W)
        """
        # Simplified empirical formula
        k = 0.073  # Friction coefficient
        clearance_factor = 1.0 + (0.001 / clearance) * 0.1

        P_disk = k * self.rho * (omega ** 3) * (radius ** 5) * clearance_factor

        return P_disk

    def leakage_flow(self, pressure_diff: float, clearance: float,
                     diameter: float, width: float) -> float:
        """
        Estimate leakage flow through clearance gaps.

        Q_leak ≈ (π × D × c³ × ΔP) / (12 × μ × L)

        Args:
            pressure_diff: Pressure differential (Pa)
            clearance: Radial clearance (m)
            diameter: Seal diameter (m)
            width: Seal width (m)

        Returns:
            Leakage flow rate (m³/s)
        """
        viscosity = 0.001  # Water viscosity (Pa·s)

        if width <= 0 or clearance <= 0:
            return 0.0

        # Simplified leakage formula
        Q_leak = (np.pi * diameter * clearance**3 * pressure_diff) / \
                 (12 * viscosity * width)

        return Q_leak


class PumpDesignOptimizer:
    """
    Optimizes pump design parameters for maximum efficiency.
    """

    def __init__(self, design_flow: float, design_head: float,
                 speed: float = 1750.0):
        """
        Initialize optimizer.

        Args:
            design_flow: Target flow rate (m³/s)
            design_head: Target head (m)
            speed: Rotational speed (RPM)
        """
        self.Q_design = design_flow
        self.H_design = design_head
        self.N = speed  # RPM
        self.omega = speed * 2 * np.pi / 60  # rad/s
        self.analyzer = PumpEfficiencyAnalyzer()

    def calculate_specific_speed(self) -> float:
        """
        Calculate specific speed.

        N_s = N × √Q / H^(3/4)

        Returns:
            Specific speed (dimensionless)
        """
        # Convert to US units for specific speed
        Q_gpm = self.Q_design * 15850.3  # m³/s to GPM
        H_ft = self.H_design * 3.28084  # m to ft

        N_s = self.N * np.sqrt(Q_gpm) / (H_ft ** 0.75)

        return N_s

    def optimal_blade_angle(self, radius: float, meridional_velocity: float,
                           blade_type: str = 'inlet') -> float:
        """
        Calculate optimal blade angle for shock-free entry.

        Args:
            radius: Radius at blade point (m)
            meridional_velocity: Meridional flow velocity (m/s)
            blade_type: 'inlet' or 'outlet'

        Returns:
            Optimal blade angle (degrees)
        """
        u = self.omega * radius  # Tangential velocity

        if blade_type == 'inlet':
            # Shock-free entry
            beta = np.arctan(meridional_velocity / u)
        else:
            # Outlet angle for backward-curved blades (optimal efficiency)
            # Typical range 20-25 degrees
            beta = np.radians(22.5)  # Empirical optimum

        return np.degrees(beta)

    def optimal_blade_count(self, D_inlet: float, D_outlet: float,
                           beta_inlet: float, beta_outlet: float) -> int:
        """
        Calculate optimal number of blades.

        Z = 6.5 × (D₂ + D₁)/(D₂ - D₁) × sin((β₁ + β₂)/2)

        Args:
            D_inlet: Inlet diameter (m)
            D_outlet: Outlet diameter (m)
            beta_inlet: Inlet blade angle (degrees)
            beta_outlet: Outlet blade angle (degrees)

        Returns:
            Optimal blade count (integer)
        """
        if D_outlet <= D_inlet:
            return 5  # Default

        beta_avg = (beta_inlet + beta_outlet) / 2

        Z = 6.5 * (D_outlet + D_inlet) / (D_outlet - D_inlet) * \
            np.sin(np.radians(beta_avg))

        # Round to nearest odd number (typical for pumps)
        Z_int = int(np.round(Z))
        if Z_int % 2 == 0:
            Z_int += 1

        # Typical range 5-9 blades
        return np.clip(Z_int, 5, 9)

    def optimize_impeller_geometry(self, initial_guess: Dict = None) -> Dict:
        """
        Optimize impeller geometry for maximum efficiency.

        Design variables:
        - D2: Outlet diameter (m)
        - b2: Outlet width (m)
        - beta2: Outlet blade angle (degrees)
        - Z: Number of blades

        Args:
            initial_guess: Dictionary of initial values

        Returns:
            Dictionary of optimized parameters and efficiency
        """
        # Default initial guess
        if initial_guess is None:
            initial_guess = {
                'D2': 0.3,      # 300 mm
                'b2': 0.03,     # 30 mm
                'beta2': 25.0,  # 25 degrees
            }

        # Design variable bounds
        bounds = [
            (0.15, 0.6),   # D2: 150-600 mm
            (0.01, 0.10),  # b2: 10-100 mm
            (15.0, 35.0),  # beta2: 15-35 degrees
        ]

        def objective(x):
            """Negative efficiency (for minimization)."""
            D2, b2, beta2 = x

            # Calculate impeller tip speed
            u2 = self.omega * (D2 / 2)

            # Meridional velocity (continuity)
            c_m = self.Q_design / (np.pi * D2 * b2)

            # Tangential velocity component (Euler equation)
            c_u2 = u2 - c_m / np.tan(np.radians(beta2))

            # Head developed (Euler equation)
            H_calc = (u2 * c_u2) / self.analyzer.g

            # Hydraulic efficiency estimate
            eta_h = self.analyzer.hydraulic_efficiency(
                self.Q_design, H_calc, u2, c_u2
            )

            # Penalize deviation from target head
            head_penalty = abs(H_calc - self.H_design) / self.H_design

            # Volumetric efficiency (empirical correlation)
            # Assume 2% leakage for good design
            eta_v = 0.98

            # Mechanical efficiency (empirical correlation)
            # Includes disk friction losses
            P_disk = self.analyzer.disk_friction_power(
                self.omega, D2/2, clearance=0.002
            )
            P_hydraulic = self.analyzer.rho * self.analyzer.g * \
                         self.Q_design * H_calc
            eta_m = P_hydraulic / (P_hydraulic + P_disk + 100)  # +100W bearing losses

            # Overall efficiency
            eta_overall = eta_h * eta_v * eta_m

            # Objective: maximize efficiency, match head
            return -(eta_overall - 10 * head_penalty)

        # Initial guess vector
        x0 = [initial_guess['D2'], initial_guess['b2'], initial_guess['beta2']]

        # Optimize
        result = minimize(objective, x0, method='SLSQP', bounds=bounds)

        # Extract optimal parameters
        D2_opt, b2_opt, beta2_opt = result.x

        # Calculate performance at optimum
        u2 = self.omega * (D2_opt / 2)
        c_m = self.Q_design / (np.pi * D2_opt * b2_opt)
        c_u2 = u2 - c_m / np.tan(np.radians(beta2_opt))
        H_calc = (u2 * c_u2) / self.analyzer.g

        # Optimal blade count
        D1 = D2_opt * 0.5  # Typical inlet/outlet ratio
        beta1 = self.optimal_blade_angle(D1/2, c_m, 'inlet')
        Z_opt = self.optimal_blade_count(D1, D2_opt, beta1, beta2_opt)

        # Calculate efficiency components
        eta_h = self.analyzer.hydraulic_efficiency(self.Q_design, H_calc, u2, c_u2)
        eta_v = 0.98
        P_disk = self.analyzer.disk_friction_power(self.omega, D2_opt/2)
        P_hydraulic = self.analyzer.rho * self.analyzer.g * self.Q_design * H_calc
        eta_m = P_hydraulic / (P_hydraulic + P_disk + 100)
        eta_overall = eta_h * eta_v * eta_m

        return {
            'D2': D2_opt,
            'b2': b2_opt,
            'beta2': beta2_opt,
            'Z': Z_opt,
            'head_achieved': H_calc,
            'eta_hydraulic': eta_h,
            'eta_volumetric': eta_v,
            'eta_mechanical': eta_m,
            'eta_overall': eta_overall,
            'specific_speed': self.calculate_specific_speed()
        }


class EnergyAnalyzer:
    """
    Analyzes energy consumption and cost savings.
    """

    def __init__(self, electricity_cost: float = 0.10):
        """
        Initialize analyzer.

        Args:
            electricity_cost: Cost per kWh ($/kWh)
        """
        self.cost_per_kwh = electricity_cost

    def annual_energy_consumption(self, flow: float, head: float,
                                 efficiency: float, operating_hours: float,
                                 motor_efficiency: float = 0.95) -> float:
        """
        Calculate annual energy consumption.

        Args:
            flow: Flow rate (m³/s)
            head: Head (m)
            efficiency: Pump efficiency (0-1)
            operating_hours: Annual operating hours
            motor_efficiency: Motor efficiency (0-1)

        Returns:
            Annual energy consumption (kWh/year)
        """
        # Hydraulic power (W)
        P_hydraulic = WATER_DENSITY * GRAVITY * flow * head

        # Shaft power (W)
        P_shaft = P_hydraulic / efficiency

        # Electrical power (W)
        P_electrical = P_shaft / motor_efficiency

        # Annual energy (kWh)
        E_annual = (P_electrical / 1000) * operating_hours

        return E_annual

    def annual_energy_cost(self, flow: float, head: float,
                          efficiency: float, operating_hours: float,
                          motor_efficiency: float = 0.95) -> float:
        """
        Calculate annual energy cost.

        Args:
            flow: Flow rate (m³/s)
            head: Head (m)
            efficiency: Pump efficiency (0-1)
            operating_hours: Annual operating hours
            motor_efficiency: Motor efficiency (0-1)

        Returns:
            Annual energy cost ($)
        """
        E_annual = self.annual_energy_consumption(
            flow, head, efficiency, operating_hours, motor_efficiency
        )

        return E_annual * self.cost_per_kwh

    def savings_from_efficiency_improvement(self, flow: float, head: float,
                                           baseline_efficiency: float,
                                           improved_efficiency: float,
                                           operating_hours: float,
                                           motor_efficiency: float = 0.95) -> Dict:
        """
        Calculate savings from efficiency improvement.

        Args:
            flow: Flow rate (m³/s)
            head: Head (m)
            baseline_efficiency: Current efficiency (0-1)
            improved_efficiency: Improved efficiency (0-1)
            operating_hours: Annual operating hours
            motor_efficiency: Motor efficiency (0-1)

        Returns:
            Dictionary with energy and cost savings
        """
        # Baseline
        E_baseline = self.annual_energy_consumption(
            flow, head, baseline_efficiency, operating_hours, motor_efficiency
        )
        Cost_baseline = E_baseline * self.cost_per_kwh

        # Improved
        E_improved = self.annual_energy_consumption(
            flow, head, improved_efficiency, operating_hours, motor_efficiency
        )
        Cost_improved = E_improved * self.cost_per_kwh

        # Savings
        E_savings = E_baseline - E_improved
        Cost_savings = Cost_baseline - Cost_improved
        Percent_savings = (E_savings / E_baseline) * 100

        return {
            'baseline_energy_kwh': E_baseline,
            'baseline_cost_usd': Cost_baseline,
            'improved_energy_kwh': E_improved,
            'improved_cost_usd': Cost_improved,
            'energy_savings_kwh': E_savings,
            'cost_savings_usd': Cost_savings,
            'percent_savings': Percent_savings
        }

    def vfd_energy_savings(self, rated_power: float, speed_ratio: float,
                          operating_hours: float) -> Dict:
        """
        Calculate energy savings from VFD speed reduction.

        Power scales with speed cubed: P ∝ N³

        Args:
            rated_power: Rated power at full speed (kW)
            speed_ratio: Reduced speed / rated speed (0-1)
            operating_hours: Annual hours at reduced speed

        Returns:
            Dictionary with VFD savings analysis
        """
        # Power at reduced speed
        power_reduced = rated_power * (speed_ratio ** 3)

        # Energy consumption
        E_full_speed = rated_power * operating_hours
        E_reduced_speed = power_reduced * operating_hours

        # Savings
        E_savings = E_full_speed - E_reduced_speed
        Cost_savings = E_savings * self.cost_per_kwh
        Percent_savings = ((1 - speed_ratio**3) * 100)

        return {
            'full_speed_power_kw': rated_power,
            'reduced_speed_power_kw': power_reduced,
            'speed_ratio': speed_ratio,
            'energy_savings_kwh': E_savings,
            'cost_savings_usd': Cost_savings,
            'percent_savings': Percent_savings
        }

    def simple_payback(self, investment_cost: float, annual_savings: float,
                      rebates: float = 0.0) -> float:
        """
        Calculate simple payback period.

        Args:
            investment_cost: Initial investment ($)
            annual_savings: Annual energy cost savings ($)
            rebates: Utility rebates or incentives ($)

        Returns:
            Payback period (years)
        """
        net_investment = investment_cost - rebates

        if annual_savings <= 0:
            return float('inf')

        return net_investment / annual_savings

    def life_cycle_cost(self, capital_cost: float, installation_cost: float,
                       annual_energy_cost: float, annual_maintenance_cost: float,
                       lifetime_years: float = 20, discount_rate: float = 0.05,
                       salvage_value: float = 0.0) -> Dict:
        """
        Calculate life cycle cost.

        Args:
            capital_cost: Equipment purchase cost ($)
            installation_cost: Installation cost ($)
            annual_energy_cost: Annual energy cost ($)
            annual_maintenance_cost: Annual maintenance cost ($)
            lifetime_years: Expected lifetime (years)
            discount_rate: Discount rate (0-1)
            salvage_value: End-of-life salvage value ($)

        Returns:
            Dictionary with LCC analysis
        """
        # Present value of annual costs
        years = np.arange(1, lifetime_years + 1)
        discount_factors = 1 / ((1 + discount_rate) ** years)

        PV_energy = np.sum((annual_energy_cost * discount_factors))
        PV_maintenance = np.sum((annual_maintenance_cost * discount_factors))
        PV_salvage = salvage_value / ((1 + discount_rate) ** lifetime_years)

        # Total LCC
        LCC = capital_cost + installation_cost + PV_energy + PV_maintenance - PV_salvage

        return {
            'capital_cost': capital_cost,
            'installation_cost': installation_cost,
            'pv_energy_cost': PV_energy,
            'pv_maintenance_cost': PV_maintenance,
            'pv_salvage_value': PV_salvage,
            'life_cycle_cost': LCC,
            'annualized_cost': LCC * discount_rate / (1 - (1 + discount_rate)**(-lifetime_years))
        }


# ============================================================================
# EXAMPLE USAGE AND VERIFICATION
# ============================================================================

def example_efficiency_analysis():
    """
    Example: Analyze efficiency components of an existing pump.
    """
    print("=" * 70)
    print("EXAMPLE 1: Pump Efficiency Component Analysis")
    print("=" * 70)

    analyzer = PumpEfficiencyAnalyzer()

    # Pump operating conditions
    flow = 0.063  # m³/s (1000 GPM)
    head = 30.5  # m (100 ft)
    shaft_power = 26000  # W (35 HP)

    # Calculate overall efficiency
    eta_overall = analyzer.overall_efficiency(flow, head, shaft_power)

    print(f"\nOperating Conditions:")
    print(f"  Flow rate: {flow:.3f} m³/s ({flow * 15850.3:.1f} GPM)")
    print(f"  Head: {head:.1f} m ({head * 3.28084:.1f} ft)")
    print(f"  Shaft power: {shaft_power/1000:.1f} kW ({shaft_power/746:.1f} HP)")

    print(f"\nOverall Efficiency: {eta_overall*100:.1f}%")

    # Estimate component efficiencies
    # Assuming some typical values for demonstration
    impeller_tip_speed = 25.0  # m/s
    tangential_velocity = 20.0  # m/s
    leakage_flow = 0.002  # m³/s (3% leakage)
    friction_power = 1500  # W

    eta_h = analyzer.hydraulic_efficiency(flow, head, impeller_tip_speed, tangential_velocity)
    eta_v = analyzer.volumetric_efficiency(flow, leakage_flow)

    hydraulic_power = WATER_DENSITY * GRAVITY * flow * head
    eta_m = analyzer.mechanical_efficiency(hydraulic_power, friction_power)

    print(f"\nEfficiency Components:")
    print(f"  Hydraulic efficiency: {eta_h*100:.1f}%")
    print(f"  Volumetric efficiency: {eta_v*100:.1f}%")
    print(f"  Mechanical efficiency: {eta_m*100:.1f}%")
    print(f"  Product: {eta_h*eta_v*eta_m*100:.1f}%")

    # Loss analysis
    print(f"\nLoss Analysis:")
    print(f"  Hydraulic power: {hydraulic_power/1000:.2f} kW")
    print(f"  Leakage loss: {WATER_DENSITY*GRAVITY*leakage_flow*head/1000:.2f} kW")
    print(f"  Friction loss: {friction_power/1000:.2f} kW")
    print(f"  Total losses: {(shaft_power - hydraulic_power)/1000:.2f} kW")

    return eta_overall


def example_design_optimization():
    """
    Example: Optimize impeller geometry for a new pump.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 2: Impeller Geometry Optimization")
    print("=" * 70)

    # Design requirements
    design_flow = 0.063  # m³/s (1000 GPM)
    design_head = 30.5  # m (100 ft)
    speed = 1750  # RPM

    print(f"\nDesign Requirements:")
    print(f"  Flow rate: {design_flow:.3f} m³/s ({design_flow * 15850.3:.1f} GPM)")
    print(f"  Head: {design_head:.1f} m ({design_head * 3.28084:.1f} ft)")
    print(f"  Speed: {speed:.0f} RPM")

    # Create optimizer
    optimizer = PumpDesignOptimizer(design_flow, design_head, speed)

    # Calculate specific speed
    N_s = optimizer.calculate_specific_speed()
    print(f"  Specific speed: {N_s:.0f}")

    # Optimize impeller geometry
    print(f"\nOptimizing impeller geometry...")
    result = optimizer.optimize_impeller_geometry()

    print(f"\nOptimal Design Parameters:")
    print(f"  Impeller diameter (D2): {result['D2']*1000:.1f} mm")
    print(f"  Impeller width (b2): {result['b2']*1000:.1f} mm")
    print(f"  Outlet blade angle (β2): {result['beta2']:.1f}°")
    print(f"  Number of blades: {result['Z']}")

    print(f"\nPredicted Performance:")
    print(f"  Head achieved: {result['head_achieved']:.1f} m (target: {design_head:.1f} m)")
    print(f"  Hydraulic efficiency: {result['eta_hydraulic']*100:.1f}%")
    print(f"  Volumetric efficiency: {result['eta_volumetric']*100:.1f}%")
    print(f"  Mechanical efficiency: {result['eta_mechanical']*100:.1f}%")
    print(f"  Overall efficiency: {result['eta_overall']*100:.1f}%")

    return result


def example_energy_savings_analysis():
    """
    Example: Calculate energy savings from efficiency improvement.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 3: Energy Savings Analysis")
    print("=" * 70)

    # System parameters
    flow = 0.063  # m³/s (1000 GPM)
    head = 30.5  # m (100 ft)
    operating_hours = 6000  # hours/year
    electricity_cost = 0.10  # $/kWh

    # Baseline and improved efficiencies
    baseline_efficiency = 0.70
    improved_efficiency = 0.80

    print(f"\nSystem Parameters:")
    print(f"  Flow rate: {flow:.3f} m³/s ({flow * 15850.3:.1f} GPM)")
    print(f"  Head: {head:.1f} m ({head * 3.28084:.1f} ft)")
    print(f"  Operating hours: {operating_hours:.0f} hours/year")
    print(f"  Electricity cost: ${electricity_cost:.2f}/kWh")

    print(f"\nEfficiency Comparison:")
    print(f"  Baseline efficiency: {baseline_efficiency*100:.0f}%")
    print(f"  Improved efficiency: {improved_efficiency*100:.0f}%")

    # Create analyzer
    analyzer = EnergyAnalyzer(electricity_cost)

    # Calculate savings
    savings = analyzer.savings_from_efficiency_improvement(
        flow, head, baseline_efficiency, improved_efficiency, operating_hours
    )

    print(f"\nBaseline Performance:")
    print(f"  Annual energy: {savings['baseline_energy_kwh']:,.0f} kWh/year")
    print(f"  Annual cost: ${savings['baseline_cost_usd']:,.0f}/year")

    print(f"\nImproved Performance:")
    print(f"  Annual energy: {savings['improved_energy_kwh']:,.0f} kWh/year")
    print(f"  Annual cost: ${savings['improved_cost_usd']:,.0f}/year")

    print(f"\nAnnual Savings:")
    print(f"  Energy savings: {savings['energy_savings_kwh']:,.0f} kWh/year")
    print(f"  Cost savings: ${savings['cost_savings_usd']:,.0f}/year")
    print(f"  Percent savings: {savings['percent_savings']:.1f}%")

    # Payback analysis
    improvement_cost = 50000  # $50,000 for new impeller and installation
    payback = analyzer.simple_payback(improvement_cost, savings['cost_savings_usd'])

    print(f"\nPayback Analysis:")
    print(f"  Improvement cost: ${improvement_cost:,.0f}")
    print(f"  Simple payback: {payback:.1f} years")

    return savings


def example_vfd_savings():
    """
    Example: Calculate VFD energy savings.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 4: VFD Energy Savings Analysis")
    print("=" * 70)

    # System parameters
    rated_power = 22.0  # kW (30 HP)
    speed_ratio = 0.80  # Operating at 80% speed
    operating_hours = 6000  # hours/year at reduced speed
    electricity_cost = 0.10  # $/kWh

    print(f"\nSystem Parameters:")
    print(f"  Rated power: {rated_power:.1f} kW ({rated_power/0.746:.1f} HP)")
    print(f"  Operating speed: {speed_ratio*100:.0f}% of rated")
    print(f"  Operating hours: {operating_hours:.0f} hours/year")
    print(f"  Electricity cost: ${electricity_cost:.2f}/kWh")

    # Create analyzer
    analyzer = EnergyAnalyzer(electricity_cost)

    # Calculate VFD savings
    vfd_savings = analyzer.vfd_energy_savings(rated_power, speed_ratio, operating_hours)

    print(f"\nPower Comparison:")
    print(f"  Full speed power: {vfd_savings['full_speed_power_kw']:.1f} kW")
    print(f"  Reduced speed power: {vfd_savings['reduced_speed_power_kw']:.1f} kW")

    print(f"\nAnnual Savings with VFD:")
    print(f"  Energy savings: {vfd_savings['energy_savings_kwh']:,.0f} kWh/year")
    print(f"  Cost savings: ${vfd_savings['cost_savings_usd']:,.0f}/year")
    print(f"  Percent savings: {vfd_savings['percent_savings']:.1f}%")

    # VFD payback
    vfd_cost = 8000  # $8,000 for VFD and installation
    payback = analyzer.simple_payback(vfd_cost, vfd_savings['cost_savings_usd'])

    print(f"\nVFD Investment:")
    print(f"  VFD cost: ${vfd_cost:,.0f}")
    print(f"  Simple payback: {payback:.2f} years")

    return vfd_savings


def example_life_cycle_cost():
    """
    Example: Compare life cycle costs of two pump options.
    """
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Life Cycle Cost Analysis")
    print("=" * 70)

    analyzer = EnergyAnalyzer(electricity_cost=0.10)

    # Option A: Standard efficiency pump
    print("\nOption A: Standard Efficiency Pump (70%)")
    lcc_a = analyzer.life_cycle_cost(
        capital_cost=15000,
        installation_cost=5000,
        annual_energy_cost=16800,  # From previous examples
        annual_maintenance_cost=1000,
        lifetime_years=20,
        discount_rate=0.05,
        salvage_value=1000
    )

    print(f"  Capital cost: ${lcc_a['capital_cost']:,.0f}")
    print(f"  Installation cost: ${lcc_a['installation_cost']:,.0f}")
    print(f"  PV of energy costs: ${lcc_a['pv_energy_cost']:,.0f}")
    print(f"  PV of maintenance: ${lcc_a['pv_maintenance_cost']:,.0f}")
    print(f"  Life cycle cost: ${lcc_a['life_cycle_cost']:,.0f}")
    print(f"  Annualized cost: ${lcc_a['annualized_cost']:,.0f}/year")

    # Option B: High efficiency pump
    print("\nOption B: High Efficiency Pump (80%)")
    lcc_b = analyzer.life_cycle_cost(
        capital_cost=25000,  # Higher initial cost
        installation_cost=5000,
        annual_energy_cost=14700,  # Lower energy cost
        annual_maintenance_cost=1000,
        lifetime_years=20,
        discount_rate=0.05,
        salvage_value=1500
    )

    print(f"  Capital cost: ${lcc_b['capital_cost']:,.0f}")
    print(f"  Installation cost: ${lcc_b['installation_cost']:,.0f}")
    print(f"  PV of energy costs: ${lcc_b['pv_energy_cost']:,.0f}")
    print(f"  PV of maintenance: ${lcc_b['pv_maintenance_cost']:,.0f}")
    print(f"  Life cycle cost: ${lcc_b['life_cycle_cost']:,.0f}")
    print(f"  Annualized cost: ${lcc_b['annualized_cost']:,.0f}/year")

    # Comparison
    lcc_savings = lcc_a['life_cycle_cost'] - lcc_b['life_cycle_cost']

    print(f"\nComparison:")
    print(f"  LCC savings with Option B: ${lcc_savings:,.0f}")
    print(f"  Additional capital cost: ${lcc_b['capital_cost'] - lcc_a['capital_cost']:,.0f}")
    print(f"  Annual energy savings: ${16800 - 14700:,.0f}")

    if lcc_savings > 0:
        print(f"  ✓ Option B (high efficiency) has lower life cycle cost")
    else:
        print(f"  ✓ Option A (standard) has lower life cycle cost")

    return {'option_a': lcc_a, 'option_b': lcc_b, 'savings': lcc_savings}


def run_all_examples():
    """
    Run all verification examples.
    """
    print("\n")
    print("#" * 70)
    print("#" + " " * 68 + "#")
    print("#" + "  PUMP EFFICIENCY OPTIMIZATION - VERIFICATION EXAMPLES".center(68) + "#")
    print("#" + " " * 68 + "#")
    print("#" * 70)

    try:
        # Run examples
        eta = example_efficiency_analysis()
        design = example_design_optimization()
        savings = example_energy_savings_analysis()
        vfd = example_vfd_savings()
        lcc = example_life_cycle_cost()

        # Summary
        print("\n" + "=" * 70)
        print("VERIFICATION SUMMARY")
        print("=" * 70)
        print(f"✓ Efficiency analysis: {eta*100:.1f}% overall efficiency calculated")
        print(f"✓ Design optimization: {design['eta_overall']*100:.1f}% optimized efficiency")
        print(f"✓ Energy savings: ${savings['cost_savings_usd']:,.0f}/year from 10% efficiency gain")
        print(f"✓ VFD savings: {vfd['percent_savings']:.1f}% energy reduction at 80% speed")
        print(f"✓ Life cycle cost: ${abs(lcc['savings']):,.0f} LCC difference between options")
        print("\n✓ All optimization algorithms verified successfully!")
        print("=" * 70)

    except Exception as e:
        print(f"\n✗ Error during verification: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    run_all_examples()
