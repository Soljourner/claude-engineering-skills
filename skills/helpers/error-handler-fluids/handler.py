"""
Error Handler for Fluid Calculations
=====================================

Robust error handling, validation, and debugging tools for production-ready
fluid mechanics calculations.

Features:
- Input validation functions
- Physical bounds checkers
- Error handling decorators
- Safe iterative solvers
- Logging and diagnostics

Units: SI unless otherwise specified
"""

import math
import warnings
import logging
from functools import wraps
import inspect


# ============================================================================
# LOGGING CONFIGURATION
# ============================================================================

# Create logger for fluid calculations
logger = logging.getLogger('fluid_error_handler')
logger.setLevel(logging.INFO)

# Create console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.WARNING)
console_formatter = logging.Formatter('%(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)
logger.addHandler(console_handler)


def configure_logging(level=logging.INFO, log_file=None):
    """
    Configure logging for fluid calculations.

    Args:
        level: Logging level (logging.DEBUG, INFO, WARNING, ERROR)
        log_file: Optional file path for logging
    """
    logger.setLevel(level)

    # Add file handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)


# ============================================================================
# INPUT VALIDATION FUNCTIONS
# ============================================================================

def validate_temperature(T, T_min=-273.15, T_max=1000, unit='C', param_name='Temperature'):
    """
    Validate temperature is within physical and practical limits.

    Args:
        T: Temperature value
        T_min: Minimum allowed temperature (default: absolute zero)
        T_max: Maximum allowed temperature
        unit: 'C' for Celsius, 'K' for Kelvin
        param_name: Parameter name for error messages

    Raises:
        ValueError: If temperature is invalid

    Returns:
        T: Validated temperature
    """
    # Check for NaN/Inf
    if math.isnan(T):
        raise ValueError(f"{param_name} is NaN")
    if math.isinf(T):
        raise ValueError(f"{param_name} is infinite")

    # Check absolute zero
    if unit == 'C':
        if T <= -273.15:
            raise ValueError(f"{param_name} {T}°C below absolute zero (-273.15°C)")
    elif unit == 'K':
        if T <= 0:
            raise ValueError(f"{param_name} {T}K below absolute zero")
    else:
        raise ValueError(f"Unknown temperature unit: {unit}")

    # Check validity range
    if T < T_min:
        raise ValueError(
            f"{param_name} {T}{unit} below minimum allowed value {T_min}{unit}"
        )
    if T > T_max:
        raise ValueError(
            f"{param_name} {T}{unit} exceeds maximum allowed value {T_max}{unit}"
        )

    logger.debug(f"{param_name} validated: {T}{unit}")
    return T


def validate_pressure(P, P_min=0, P_max=None, allow_vacuum=False, param_name='Pressure'):
    """
    Validate pressure is positive (or allow vacuum if specified).

    Args:
        P: Pressure in Pa
        P_min: Minimum allowed pressure (default: 0 for absolute)
        P_max: Maximum allowed pressure
        allow_vacuum: Allow negative values (for gauge pressure)
        param_name: Parameter name for error messages

    Raises:
        ValueError: If pressure is invalid

    Returns:
        P: Validated pressure
    """
    # Check for NaN/Inf
    if math.isnan(P):
        raise ValueError(f"{param_name} is NaN")
    if math.isinf(P):
        raise ValueError(f"{param_name} is infinite")

    # Check positivity (unless vacuum allowed)
    if not allow_vacuum and P < 0:
        raise ValueError(f"Absolute {param_name} {P} Pa cannot be negative")

    # Check minimum
    if P_min is not None and P < P_min:
        raise ValueError(f"{param_name} {P} Pa below minimum {P_min} Pa")

    # Check maximum
    if P_max is not None and P > P_max:
        raise ValueError(f"{param_name} {P} Pa exceeds maximum {P_max} Pa")

    logger.debug(f"{param_name} validated: {P} Pa")
    return P


def validate_positive(value, param_name='Parameter', allow_zero=False):
    """
    Validate that a parameter is positive (or non-negative if zero allowed).

    Args:
        value: Value to check
        param_name: Parameter name for error messages
        allow_zero: Allow value to be exactly zero

    Raises:
        ValueError: If value is not positive

    Returns:
        value: Validated value
    """
    # Check for NaN/Inf
    if math.isnan(value):
        raise ValueError(f"{param_name} is NaN")
    if math.isinf(value):
        raise ValueError(f"{param_name} is infinite")

    # Check positivity
    if allow_zero:
        if value < 0:
            raise ValueError(f"{param_name} must be non-negative, got {value}")
    else:
        if value <= 0:
            raise ValueError(f"{param_name} must be positive, got {value}")

    logger.debug(f"{param_name} validated: {value}")
    return value


def validate_range(value, min_val=None, max_val=None, param_name='Parameter'):
    """
    Validate that a parameter is within specified range.

    Args:
        value: Value to check
        min_val: Minimum allowed value (None = no minimum)
        max_val: Maximum allowed value (None = no maximum)
        param_name: Parameter name for error messages

    Raises:
        ValueError: If value is outside range

    Returns:
        value: Validated value
    """
    # Check for NaN/Inf
    if math.isnan(value):
        raise ValueError(f"{param_name} is NaN")
    if math.isinf(value):
        raise ValueError(f"{param_name} is infinite")

    # Check minimum
    if min_val is not None and value < min_val:
        raise ValueError(f"{param_name} = {value} below minimum {min_val}")

    # Check maximum
    if max_val is not None and value > max_val:
        raise ValueError(f"{param_name} = {value} exceeds maximum {max_val}")

    logger.debug(f"{param_name} validated: {value}")
    return value


def validate_reynolds_number(Re, warn_transition=True):
    """
    Validate Reynolds number and optionally warn about transition region.

    Args:
        Re: Reynolds number
        warn_transition: Warn if in transition region (2300-4000)

    Raises:
        ValueError: If Reynolds number is invalid

    Returns:
        Re: Validated Reynolds number
    """
    # Check for NaN/Inf
    if math.isnan(Re):
        raise ValueError("Reynolds number is NaN")
    if math.isinf(Re):
        raise ValueError("Reynolds number is infinite")

    # Must be non-negative
    if Re < 0:
        raise ValueError(f"Reynolds number {Re} cannot be negative")

    # Warn about transition region
    if warn_transition and 2300 < Re < 4000:
        warnings.warn(
            f"Reynolds number {Re:.0f} in transition region (2300-4000). "
            "Flow regime uncertain."
        )

    logger.debug(f"Reynolds number validated: {Re:.0f}")
    return Re


# ============================================================================
# PHYSICAL BOUNDS CHECKERS
# ============================================================================

def check_cavitation_risk(P, P_vapor, location='system'):
    """
    Check if pressure is above vapor pressure to prevent cavitation.

    Args:
        P: Local pressure in Pa
        P_vapor: Vapor pressure at fluid temperature in Pa
        location: Description of location for warning message

    Raises:
        ValueError: If cavitation conditions detected

    Returns:
        margin: Pressure margin above vapor pressure
    """
    margin = P - P_vapor

    if margin < 0:
        raise ValueError(
            f"Cavitation risk at {location}: "
            f"Pressure {P/1000:.2f} kPa below vapor pressure {P_vapor/1000:.2f} kPa"
        )

    # Warn if margin is small
    if margin < 0.1 * P_vapor:
        warnings.warn(
            f"Low cavitation margin at {location}: "
            f"{margin/1000:.2f} kPa above vapor pressure"
        )

    logger.debug(f"Cavitation check passed at {location}: margin = {margin/1000:.2f} kPa")
    return margin


def check_compressibility(M, M_max=0.3):
    """
    Check if Mach number is low enough for incompressible flow assumption.

    Args:
        M: Mach number
        M_max: Maximum Mach number for incompressible flow (default: 0.3)

    Returns:
        compressible: True if compressibility effects are significant
    """
    if M > M_max:
        warnings.warn(
            f"Mach number {M:.2f} exceeds {M_max}. "
            "Compressibility effects may be significant."
        )
        return True

    logger.debug(f"Compressibility check passed: M = {M:.3f}")
    return False


def check_numerical_stability(value, min_val=1e-10, max_val=1e10, param_name='Value'):
    """
    Check if value is in numerically stable range.

    Args:
        value: Value to check
        min_val: Minimum for underflow warning
        max_val: Maximum for overflow warning
        param_name: Parameter name for messages

    Raises:
        ValueError: If value is NaN or Inf
    """
    if math.isnan(value):
        raise ValueError(f"{param_name} is NaN")

    if math.isinf(value):
        raise ValueError(f"{param_name} is infinite")

    if abs(value) < min_val:
        warnings.warn(
            f"{param_name} = {value:.2e} may cause underflow (< {min_val:.2e})"
        )

    if abs(value) > max_val:
        warnings.warn(
            f"{param_name} = {value:.2e} may cause overflow (> {max_val:.2e})"
        )


# ============================================================================
# ERROR HANDLING DECORATORS
# ============================================================================

def validate_positive_args(*arg_names):
    """
    Decorator to validate that specified arguments are positive.

    Usage:
        @validate_positive_args('velocity', 'diameter')
        def reynolds_number(velocity, diameter, viscosity):
            return velocity * diameter / viscosity
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get function signature
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # Check specified arguments
            for arg_name in arg_names:
                if arg_name in bound_args.arguments:
                    value = bound_args.arguments[arg_name]
                    if value <= 0:
                        raise ValueError(
                            f"Argument '{arg_name}' must be positive, got {value}"
                        )

            return func(*args, **kwargs)
        return wrapper
    return decorator


def validate_range_decorator(param_name, min_val=None, max_val=None):
    """
    Decorator to validate parameter is within specified range.

    Usage:
        @validate_range_decorator('temperature', min_val=0, max_val=100)
        def water_density(temperature):
            return 1000 - 0.1 * temperature
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Get function signature
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # Check parameter
            if param_name in bound_args.arguments:
                value = bound_args.arguments[param_name]

                if min_val is not None and value < min_val:
                    raise ValueError(
                        f"{param_name} = {value} below minimum {min_val}"
                    )

                if max_val is not None and value > max_val:
                    raise ValueError(
                        f"{param_name} = {value} exceeds maximum {max_val}"
                    )

            return func(*args, **kwargs)
        return wrapper
    return decorator


def log_errors(log_level=logging.ERROR):
    """
    Decorator to log errors and re-raise.

    Usage:
        @log_errors()
        def risky_calculation(x):
            return 1 / x
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logger.log(
                    log_level,
                    f"Error in {func.__name__}: {e}",
                    exc_info=True
                )
                raise
        return wrapper
    return decorator


def safe_divide(numerator, denominator, default=None, epsilon=1e-15):
    """
    Safely divide two numbers, checking for zero denominator.

    Args:
        numerator: Numerator value
        denominator: Denominator value
        default: Default value if division fails (None = raise error)
        epsilon: Threshold for considering denominator as zero

    Returns:
        result: numerator / denominator, or default if fails

    Raises:
        ValueError: If denominator is zero and no default provided
    """
    if abs(denominator) < epsilon:
        if default is not None:
            warnings.warn(f"Division by near-zero value: returning default {default}")
            return default
        else:
            raise ValueError(f"Division by zero: denominator = {denominator}")

    return numerator / denominator


# ============================================================================
# ROBUST CALCULATION FUNCTIONS
# ============================================================================

@validate_positive_args('velocity', 'length')
def safe_reynolds_number(velocity, length, kinematic_viscosity=None,
                        density=None, dynamic_viscosity=None):
    """
    Calculate Reynolds number with comprehensive error handling.

    Args:
        velocity: Flow velocity in m/s
        length: Characteristic length in m
        kinematic_viscosity: Kinematic viscosity in m²/s
        OR
        density: Density in kg/m³
        dynamic_viscosity: Dynamic viscosity in Pa·s

    Returns:
        Re: Reynolds number

    Raises:
        ValueError: If inputs are invalid
    """
    try:
        # Method 1: Using kinematic viscosity
        if kinematic_viscosity is not None:
            validate_positive(kinematic_viscosity, 'kinematic_viscosity')
            Re = safe_divide(velocity * length, kinematic_viscosity)

        # Method 2: Using density and dynamic viscosity
        elif density is not None and dynamic_viscosity is not None:
            validate_positive(density, 'density')
            validate_positive(dynamic_viscosity, 'dynamic_viscosity')
            Re = safe_divide(density * velocity * length, dynamic_viscosity)

        else:
            raise ValueError(
                "Must provide either kinematic_viscosity or "
                "(density and dynamic_viscosity)"
            )

        # Validate result
        validate_reynolds_number(Re)

        logger.info(f"Reynolds number calculated: {Re:.0f}")
        return Re

    except Exception as e:
        logger.error(f"Failed to calculate Reynolds number: {e}")
        raise


def safe_friction_factor(Re, roughness=0, diameter=1, method='auto', max_iter=100):
    """
    Calculate friction factor with multiple fallback methods.

    Args:
        Re: Reynolds number
        roughness: Absolute roughness in m (default: 0 for smooth)
        diameter: Pipe diameter in m
        method: 'auto', 'colebrook', 'swamee-jain', 'blasius', or 'laminar'
        max_iter: Maximum iterations for Colebrook

    Returns:
        f: Darcy friction factor

    Raises:
        ValueError: If calculation fails with all methods
    """
    try:
        # Validate inputs
        validate_reynolds_number(Re, warn_transition=True)
        validate_positive(diameter, 'diameter')
        validate_positive(roughness, 'roughness', allow_zero=True)

        if roughness > diameter:
            raise ValueError(f"Roughness {roughness} m exceeds diameter {diameter} m")

        # Laminar flow
        if Re < 2300 and method in ['auto', 'laminar']:
            f = 64 / Re
            logger.info(f"Laminar friction factor: {f:.5f}")
            return f

        # Transition region
        if 2300 <= Re < 4000 and method == 'auto':
            # Linear interpolation
            f_lam = 64 / 2300
            if roughness == 0:
                f_turb = 0.316 / 4000**0.25
            else:
                rel_roughness = roughness / diameter
                f_turb = 0.25 / (math.log10(rel_roughness/3.7 + 5.74/4000**0.9))**2

            f = f_lam + (f_turb - f_lam) * (Re - 2300) / (4000 - 2300)
            logger.info(f"Transition friction factor: {f:.5f}")
            return f

        # Turbulent flow - try methods in order
        if method == 'auto' or method == 'colebrook':
            try:
                f = friction_factor_colebrook_safe(Re, roughness, diameter, max_iter)
                logger.info(f"Colebrook friction factor: {f:.5f}")
                return f
            except Exception as e:
                if method == 'colebrook':
                    raise
                logger.warning(f"Colebrook failed: {e}, trying Swamee-Jain")

        if method == 'auto' or method == 'swamee-jain':
            try:
                rel_roughness = roughness / diameter
                f = 0.25 / (math.log10(rel_roughness/3.7 + 5.74/Re**0.9))**2
                logger.info(f"Swamee-Jain friction factor: {f:.5f}")
                return f
            except Exception as e:
                if method == 'swamee-jain':
                    raise
                logger.warning(f"Swamee-Jain failed: {e}, trying Blasius")

        if method == 'auto' or method == 'blasius':
            if Re < 100000 or method == 'blasius':
                f = 0.316 / Re**0.25
                logger.info(f"Blasius friction factor: {f:.5f}")
                return f

        # All methods failed
        raise ValueError("Could not calculate friction factor with any method")

    except Exception as e:
        logger.error(f"Failed to calculate friction factor: {e}")
        raise


def friction_factor_colebrook_safe(Re, roughness, diameter, max_iter=100, tol=1e-6):
    """
    Solve Colebrook equation with convergence monitoring.

    Args:
        Re: Reynolds number
        roughness: Absolute roughness in m
        diameter: Pipe diameter in m
        max_iter: Maximum iterations
        tol: Convergence tolerance

    Returns:
        f: Darcy friction factor

    Raises:
        RuntimeError: If iteration fails to converge or diverges
    """
    rel_roughness = roughness / diameter

    # Initial guess (Swamee-Jain)
    try:
        f = 0.25 / (math.log10(rel_roughness/3.7 + 5.74/Re**0.9))**2
    except:
        f = 0.02  # Fallback guess

    convergence_history = []

    for iteration in range(max_iter):
        f_old = f

        # Colebrook-White equation: 1/√f = -2 log₁₀(ε/3.7D + 2.51/(Re√f))
        try:
            term = rel_roughness/3.7 + 2.51/(Re * math.sqrt(f))

            if term <= 0:
                raise ValueError("Invalid argument to logarithm")

            f = 1 / (-2 * math.log10(term))**2

            # Check for NaN/Inf
            check_numerical_stability(f, param_name='friction factor')

        except Exception as e:
            raise RuntimeError(f"Colebrook calculation error at iteration {iteration}: {e}")

        # Check convergence
        error = abs(f - f_old)
        convergence_history.append(error)

        if error < tol:
            logger.debug(f"Colebrook converged in {iteration} iterations")
            return f

        # Check for divergence
        if iteration > 10 and error > max(convergence_history[-5:]):
            raise RuntimeError(f"Colebrook iteration diverging at iteration {iteration}")

    raise RuntimeError(f"Colebrook failed to converge after {max_iter} iterations")


@validate_range_decorator('temperature', min_val=0, max_val=100)
def safe_water_density(temperature):
    """
    Calculate water density with validation.

    Args:
        temperature: Temperature in °C (valid: 0-100°C)

    Returns:
        rho: Density in kg/m³
    """
    # Simplified polynomial correlation
    rho = (999.83952
           + 16.945176 * temperature
           - 7.9870401e-3 * temperature**2
           - 46.170461e-6 * temperature**3
           + 105.56302e-9 * temperature**4
           - 280.54253e-12 * temperature**5) / (1 + 16.879850e-3 * temperature)

    check_numerical_stability(rho, param_name='water density')

    logger.debug(f"Water density at {temperature}°C: {rho:.2f} kg/m³")
    return rho


# ============================================================================
# DIAGNOSTICS
# ============================================================================

def diagnose_calculation(func, *args, **kwargs):
    """
    Diagnose why a calculation is failing.

    Args:
        func: Function to diagnose
        *args: Positional arguments
        **kwargs: Keyword arguments

    Returns:
        result: Calculation result if successful, None otherwise
    """
    print("=" * 70)
    print("CALCULATION DIAGNOSTICS")
    print("=" * 70)

    print(f"\nFunction: {func.__name__}")
    print(f"\nPositional arguments:")
    for i, arg in enumerate(args):
        print(f"  arg[{i}]: {type(arg).__name__} = {arg}")

    print(f"\nKeyword arguments:")
    for key, value in kwargs.items():
        print(f"  {key}: {type(value).__name__} = {value}")

    # Check for numerical issues
    print("\nNumerical checks:")
    all_args = list(args) + list(kwargs.values())
    for i, arg in enumerate(all_args):
        if isinstance(arg, (int, float)):
            issues = []
            if math.isnan(arg):
                issues.append("NaN")
            if math.isinf(arg):
                issues.append("Inf")
            if arg == 0:
                issues.append("Zero")
            if arg < 0:
                issues.append("Negative")

            if issues:
                print(f"  arg[{i}]: {', '.join(issues)}")

    # Attempt calculation
    print("\nAttempting calculation:")
    try:
        result = func(*args, **kwargs)
        print(f"  ✓ Success! Result: {result}")

        # Check result
        if isinstance(result, (int, float)):
            if math.isnan(result):
                print("  ⚠ Warning: Result is NaN")
            elif math.isinf(result):
                print("  ⚠ Warning: Result is Inf")

        print("=" * 70)
        return result

    except Exception as e:
        print(f"  ✗ Failed with {type(e).__name__}: {e}")
        print("\nTraceback:")
        import traceback
        traceback.print_exc()
        print("=" * 70)
        return None


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("Error Handler for Fluid Calculations - Examples")
    print("=" * 70)

    # Configure logging for examples
    configure_logging(level=logging.INFO)

    # Example 1: Safe Reynolds number
    print("\n1. Safe Reynolds number calculation:")
    try:
        Re = safe_reynolds_number(velocity=2.5, length=0.1, kinematic_viscosity=1e-6)
        print(f"   Reynolds number: {Re:.0f}")
    except Exception as e:
        print(f"   Error: {e}")

    # Example 2: Invalid Reynolds number (negative velocity)
    print("\n2. Invalid Reynolds number (should fail):")
    try:
        Re = safe_reynolds_number(velocity=-2.5, length=0.1, kinematic_viscosity=1e-6)
        print(f"   Reynolds number: {Re:.0f}")
    except Exception as e:
        print(f"   Expected error: {e}")

    # Example 3: Friction factor with fallback
    print("\n3. Friction factor with automatic method selection:")
    try:
        f = safe_friction_factor(Re=50000, roughness=0.045e-3, diameter=0.1)
        print(f"   Friction factor: {f:.5f}")
    except Exception as e:
        print(f"   Error: {e}")

    # Example 4: Temperature validation
    print("\n4. Temperature validation:")
    try:
        T = validate_temperature(25, T_min=0, T_max=100)
        print(f"   Valid temperature: {T}°C")
    except Exception as e:
        print(f"   Error: {e}")

    # Example 5: Temperature out of range
    print("\n5. Temperature out of range (should fail):")
    try:
        T = validate_temperature(150, T_min=0, T_max=100)
        print(f"   Temperature: {T}°C")
    except Exception as e:
        print(f"   Expected error: {e}")

    # Example 6: Safe division
    print("\n6. Safe division with zero denominator:")
    try:
        result = safe_divide(10, 1e-20, default=0)
        print(f"   Result: {result}")
    except Exception as e:
        print(f"   Error: {e}")

    # Example 7: Cavitation check
    print("\n7. Cavitation risk check:")
    try:
        P = 50000  # 50 kPa
        P_vapor = 2340  # Water at 20°C
        margin = check_cavitation_risk(P, P_vapor, location='pump inlet')
        print(f"   Cavitation margin: {margin/1000:.2f} kPa")
    except Exception as e:
        print(f"   Error: {e}")

    print("\n" + "=" * 70)
    print("Examples completed!")
