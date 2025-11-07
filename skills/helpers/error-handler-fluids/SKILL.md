---
name: error-handler-fluids
description: "Handle common numerical errors in fluid calculations with validation checks"
category: helpers
domain: fluids
complexity: basic
dependencies: []
---

# Error Handler for Fluid Calculations

Robust error handling and validation framework for production-ready fluid mechanics calculations. Prevents common numerical errors, validates physical bounds, and provides graceful degradation strategies.

## Overview

Fluid calculations are prone to numerical instabilities and physical constraint violations. This helper provides comprehensive error handling, input validation, and debugging tools to create production-ready code.

## Common Errors in Fluid Calculations

### 1. Division by Zero Errors

**Mach Number at Zero Velocity**
```python
# Dangerous - can divide by zero
M = V / a  # If a (speed of sound) = 0, this fails
```

**Causes:**
- Zero or near-zero speed of sound at extreme temperatures
- Zero velocity in stagnation point calculations
- Undefined fluid properties at phase boundaries

**Prevention:**
- Check denominators before division
- Use epsilon values for near-zero comparisons
- Validate property calculations

### 2. Invalid Reynolds Numbers

**Negative Reynolds Numbers**
```python
# Can occur with invalid inputs
Re = rho * V * D / mu  # If any parameter is negative
```

**Causes:**
- Negative density (unphysical)
- Negative viscosity (calculation error)
- Negative dimensions (input error)
- Invalid property correlations

**Prevention:**
- Validate all inputs are positive
- Check property correlation validity ranges
- Ensure dimensional consistency

### 3. Temperature Out of Range

**Property Correlation Failures**
```python
# Most correlations have limited validity
rho = water_density(T)  # Fails if T outside 0-100°C range
```

**Causes:**
- Input temperatures outside correlation validity
- Supercritical conditions
- Phase change regions
- Cryogenic or high-temperature conditions

**Prevention:**
- Check temperature bounds before property calculations
- Validate against correlation validity ranges
- Provide appropriate warnings

### 4. Pressure Below Vapor Pressure

**Cavitation Risk**
```python
# Physical impossibility in liquid systems
if P < P_vapor(T):
    # Liquid will vaporize - calculation invalid
```

**Causes:**
- Low pressure in pump suction
- High velocity in pipe restrictions
- Elevated temperatures reducing vapor pressure margin
- Altitude effects on atmospheric pressure

**Prevention:**
- Calculate and check vapor pressure
- Implement Net Positive Suction Head (NPSH) checks
- Validate system pressure throughout flow path

### 5. Non-Converging Iterative Solvers

**Friction Factor Iterations**
```python
# Colebrook equation may not converge
for i in range(max_iter):
    f_new = solve_colebrook(Re, roughness, f_old)
    if not converged:
        # What to do here?
```

**Causes:**
- Poor initial guess
- Extreme parameter values
- Stiff equations
- Numerical instability

**Prevention:**
- Set maximum iterations
- Check convergence criteria
- Use stable initial guesses
- Implement fallback methods

### 6. Domain Errors in Mathematical Functions

**Logarithm and Square Root Errors**
```python
# Common in pressure drop calculations
log(negative_value)  # Math domain error
sqrt(negative_value)  # Complex result in real calculation
```

**Causes:**
- Negative arguments to log/sqrt functions
- Complex number results in real calculations
- Numerical precision issues

**Prevention:**
- Validate function arguments
- Use absolute values where appropriate
- Check for numerical underflow/overflow

### 7. Unit Inconsistency Errors

**Mixed Unit Systems**
```python
# Subtle bugs from unit mixing
P_psi = 150
rho_SI = 1000  # kg/m³
# Calculation uses mixed units - wrong result!
```

**Causes:**
- Mixing imperial and SI units
- Inconsistent temperature scales (°C vs K)
- Pressure units (absolute vs gauge)

**Prevention:**
- Standardize to single unit system
- Use unit-aware libraries (pint)
- Document units in variable names

## Validation Checks

### Physical Limits Validation

#### Temperature Bounds
```python
def validate_temperature(T, T_min=-273.15, T_max=1000, unit='C'):
    """
    Validate temperature is within physical and practical limits.

    Physical minimum: absolute zero (-273.15°C)
    Practical maximum: depends on application
    """
    if unit == 'C':
        if T <= -273.15:
            raise ValueError(f"Temperature {T}°C below absolute zero")
    elif unit == 'K':
        if T <= 0:
            raise ValueError(f"Temperature {T}K below absolute zero")

    if T < T_min or T > T_max:
        raise ValueError(f"Temperature {T}{unit} outside valid range [{T_min}, {T_max}]{unit}")
```

#### Pressure Bounds
```python
def validate_pressure(P, P_min=0, P_max=None, allow_vacuum=False):
    """
    Validate pressure is positive (or allow vacuum if specified).

    Absolute pressure must be positive.
    Gauge pressure can be negative (vacuum).
    """
    if not allow_vacuum and P < 0:
        raise ValueError(f"Absolute pressure {P} Pa cannot be negative")

    if P_min is not None and P < P_min:
        raise ValueError(f"Pressure {P} Pa below minimum {P_min} Pa")

    if P_max is not None and P > P_max:
        raise ValueError(f"Pressure {P} Pa exceeds maximum {P_max} Pa")
```

#### Flow Rate Validation
```python
def validate_flow_rate(Q, Q_min=0, Q_max=None):
    """
    Validate volumetric flow rate is non-negative.
    """
    if Q < 0:
        raise ValueError(f"Flow rate {Q} m³/s cannot be negative")

    if Q_max is not None and Q > Q_max:
        raise ValueError(f"Flow rate {Q} m³/s exceeds maximum {Q_max} m³/s")
```

#### Dimensionless Number Validation
```python
def validate_reynolds_number(Re, warn_transition=True):
    """
    Validate Reynolds number and warn about transition region.
    """
    if Re < 0:
        raise ValueError(f"Reynolds number {Re} cannot be negative")

    if warn_transition and 2300 < Re < 4000:
        warnings.warn(f"Reynolds number {Re} in transition region (2300-4000)")

    return Re
```

### Unit Consistency Checks

```python
def check_dimensional_consistency(velocity, length, viscosity):
    """
    Verify dimensions are consistent for Reynolds number calculation.

    Re = V * L / nu
    [dimensionless] = [m/s] * [m] / [m²/s]
    """
    # Using pint for dimensional analysis
    import pint
    ureg = pint.UnitRegistry()

    try:
        V = ureg.Quantity(velocity)
        L = ureg.Quantity(length)
        nu = ureg.Quantity(viscosity)

        Re = (V * L / nu).to_base_units()

        if not Re.dimensionless:
            raise ValueError("Reynolds number calculation not dimensionless!")

    except pint.DimensionalityError as e:
        raise ValueError(f"Dimensional inconsistency: {e}")
```

### Numerical Stability Checks

```python
def check_numerical_stability(value, min_val=1e-10, max_val=1e10):
    """
    Check if value is in numerically stable range.

    Prevents overflow, underflow, and loss of precision.
    """
    if abs(value) < min_val:
        warnings.warn(f"Value {value} may cause underflow")

    if abs(value) > max_val:
        warnings.warn(f"Value {value} may cause overflow")

    if math.isnan(value):
        raise ValueError("Calculation resulted in NaN")

    if math.isinf(value):
        raise ValueError("Calculation resulted in infinity")
```

## Error Handling Strategies

### 1. Try-Except Blocks

**Basic Error Handling**
```python
def safe_reynolds_number(V, D, nu):
    """
    Calculate Reynolds number with error handling.
    """
    try:
        # Validate inputs
        if V < 0 or D <= 0 or nu <= 0:
            raise ValueError("Invalid inputs for Reynolds number")

        # Check for division by zero
        if nu < 1e-15:
            raise ValueError("Kinematic viscosity too small (division by zero risk)")

        Re = V * D / nu

        # Check result
        if math.isnan(Re) or math.isinf(Re):
            raise ValueError("Invalid Reynolds number result")

        return Re

    except ValueError as e:
        print(f"Error calculating Reynolds number: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None
```

### 2. Input Validation

**Pre-calculation Validation**
```python
def validate_pipe_flow_inputs(V, D, rho, mu, roughness):
    """
    Comprehensive input validation for pipe flow calculations.
    """
    errors = []

    # Velocity
    if V < 0:
        errors.append(f"Velocity {V} m/s cannot be negative")
    elif V > 100:
        errors.append(f"Velocity {V} m/s unusually high (check units)")

    # Diameter
    if D <= 0:
        errors.append(f"Diameter {D} m must be positive")
    elif D > 10:
        errors.append(f"Diameter {D} m unusually large (check units)")

    # Density
    if rho <= 0:
        errors.append(f"Density {rho} kg/m³ must be positive")
    elif rho < 0.1 or rho > 20000:
        errors.append(f"Density {rho} kg/m³ outside typical range")

    # Viscosity
    if mu <= 0:
        errors.append(f"Viscosity {mu} Pa·s must be positive")
    elif mu < 1e-6 or mu > 10:
        errors.append(f"Viscosity {mu} Pa·s outside typical range")

    # Roughness
    if roughness < 0:
        errors.append(f"Roughness {roughness} m cannot be negative")
    elif roughness > D:
        errors.append(f"Roughness {roughness} m cannot exceed diameter {D} m")

    if errors:
        raise ValueError("Input validation failed:\n" + "\n".join(errors))

    return True
```

### 3. Graceful Degradation

**Fallback Methods**
```python
def friction_factor_robust(Re, roughness, diameter, method='auto'):
    """
    Calculate friction factor with multiple fallback methods.
    """
    # Try exact Colebrook iteration first
    if method == 'auto' or method == 'colebrook':
        try:
            return friction_factor_colebrook(Re, roughness, diameter, max_iter=50)
        except (ValueError, RuntimeError) as e:
            warnings.warn(f"Colebrook method failed: {e}, trying Swamee-Jain")

    # Fallback to Swamee-Jain approximation
    if method == 'auto' or method == 'swamee-jain':
        try:
            return friction_factor_swamee_jain(Re, roughness, diameter)
        except (ValueError, RuntimeError) as e:
            warnings.warn(f"Swamee-Jain failed: {e}, trying simple correlation")

    # Final fallback to simple turbulent correlation
    if Re >= 4000:
        return 0.316 / Re**0.25
    else:
        return 64 / Re
```

### 4. Bounded Iterations

**Safe Iterative Solvers**
```python
def solve_colebrook_safe(Re, rel_roughness, max_iter=100, tol=1e-6):
    """
    Solve Colebrook equation with convergence monitoring.
    """
    # Initial guess
    f = 0.02

    convergence_history = []

    for iteration in range(max_iter):
        f_old = f

        # Colebrook-White equation
        try:
            term = rel_roughness/3.7 + 2.51/(Re * math.sqrt(f))
            if term <= 0:
                raise ValueError("Invalid argument to logarithm")

            f = 1 / (-2 * math.log10(term))**2

        except (ValueError, ZeroDivisionError) as e:
            raise RuntimeError(f"Colebrook calculation error at iteration {iteration}: {e}")

        # Check convergence
        error = abs(f - f_old)
        convergence_history.append(error)

        if error < tol:
            return f, iteration, convergence_history

        # Check for divergence
        if iteration > 10 and error > convergence_history[iteration-5]:
            raise RuntimeError("Iteration diverging")

    raise RuntimeError(f"Failed to converge after {max_iter} iterations")
```

## Error Handling Decorators

### Input Validation Decorator

```python
def validate_positive(*arg_names):
    """
    Decorator to validate that specified arguments are positive.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Get function signature
            import inspect
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            # Check specified arguments
            for arg_name in arg_names:
                if arg_name in bound_args.arguments:
                    value = bound_args.arguments[arg_name]
                    if value <= 0:
                        raise ValueError(f"{arg_name} must be positive, got {value}")

            return func(*args, **kwargs)
        return wrapper
    return decorator


@validate_positive('velocity', 'diameter', 'viscosity')
def reynolds_number(velocity, diameter, viscosity):
    """Reynolds number with automatic validation."""
    return velocity * diameter / viscosity
```

### Physical Bounds Decorator

```python
def validate_range(param_name, min_val=None, max_val=None):
    """
    Decorator to validate parameter is within specified range.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            import inspect
            sig = inspect.signature(func)
            bound_args = sig.bind(*args, **kwargs)
            bound_args.apply_defaults()

            if param_name in bound_args.arguments:
                value = bound_args.arguments[param_name]

                if min_val is not None and value < min_val:
                    raise ValueError(f"{param_name} = {value} below minimum {min_val}")

                if max_val is not None and value > max_val:
                    raise ValueError(f"{param_name} = {value} exceeds maximum {max_val}")

            return func(*args, **kwargs)
        return wrapper
    return decorator


@validate_range('temperature', min_val=0, max_val=100)
def water_density(temperature):
    """Water density with automatic temperature validation."""
    return 1000 - 0.1 * temperature  # Simplified correlation
```

### Error Logging Decorator

```python
def log_errors(logger=None):
    """
    Decorator to log errors and optionally re-raise.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                import logging
                log = logger or logging.getLogger(func.__module__)
                log.error(f"Error in {func.__name__}: {e}", exc_info=True)
                raise
        return wrapper
    return decorator
```

## Logging and Debugging

### Calculation Logging

```python
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('fluid_calculations.log'),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger('fluid_calcs')

def calculate_pressure_drop(Q, D, L, rho, mu, roughness):
    """
    Calculate pressure drop with comprehensive logging.
    """
    logger.info(f"Starting pressure drop calculation")
    logger.debug(f"Inputs: Q={Q}, D={D}, L={L}, rho={rho}, mu={mu}, roughness={roughness}")

    try:
        # Calculate Reynolds number
        V = 4 * Q / (math.pi * D**2)
        Re = rho * V * D / mu
        logger.debug(f"Velocity: {V:.3f} m/s, Reynolds: {Re:.0f}")

        # Calculate friction factor
        f = friction_factor(Re, roughness, D)
        logger.debug(f"Friction factor: {f:.5f}")

        # Calculate pressure drop
        dP = f * (L/D) * (rho * V**2 / 2)
        logger.info(f"Pressure drop: {dP:.2f} Pa")

        return dP

    except Exception as e:
        logger.error(f"Calculation failed: {e}", exc_info=True)
        raise
```

### Debug Mode

```python
class FluidCalculator:
    """
    Fluid calculator with debug mode for detailed output.
    """
    def __init__(self, debug=False):
        self.debug = debug

    def _log(self, message):
        """Internal logging."""
        if self.debug:
            print(f"[DEBUG] {message}")

    def calculate_reynolds(self, V, D, nu):
        """Calculate Reynolds number with debug output."""
        self._log(f"Input: V={V}, D={D}, nu={nu}")

        # Validation
        if nu < 1e-15:
            self._log("Warning: viscosity very small")

        Re = V * D / nu
        self._log(f"Result: Re={Re}")

        # Flow regime
        if Re < 2300:
            regime = "laminar"
        elif Re < 4000:
            regime = "transition"
        else:
            regime = "turbulent"

        self._log(f"Flow regime: {regime}")

        return Re
```

### Diagnostic Tools

```python
def diagnose_calculation_failure(func, *args, **kwargs):
    """
    Diagnose why a calculation is failing.
    """
    print("=" * 60)
    print("CALCULATION DIAGNOSTICS")
    print("=" * 60)

    print(f"\nFunction: {func.__name__}")
    print(f"Arguments: {args}")
    print(f"Keyword arguments: {kwargs}")

    # Check argument types
    print("\nArgument types:")
    for i, arg in enumerate(args):
        print(f"  arg[{i}]: {type(arg).__name__} = {arg}")

    # Check for NaN/Inf
    print("\nNumerical checks:")
    for i, arg in enumerate(args):
        if isinstance(arg, (int, float)):
            if math.isnan(arg):
                print(f"  arg[{i}] is NaN!")
            elif math.isinf(arg):
                print(f"  arg[{i}] is Inf!")
            elif arg == 0:
                print(f"  arg[{i}] is zero (potential division issue)")

    # Try to execute
    print("\nAttempting calculation:")
    try:
        result = func(*args, **kwargs)
        print(f"  Success! Result: {result}")
        return result
    except Exception as e:
        print(f"  Failed with error: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
        return None
```

## Best Practices Summary

1. **Always validate inputs** before calculations
2. **Check denominators** before division operations
3. **Verify physical constraints** (positive values, valid ranges)
4. **Use try-except blocks** for error-prone operations
5. **Implement fallback methods** for iterative solvers
6. **Log calculations** for debugging and auditing
7. **Test edge cases** (zero, very large/small values)
8. **Document assumptions** and validity ranges
9. **Use decorators** for consistent validation
10. **Provide clear error messages** with context

## When to Use This Helper

**Essential for:**
- Production code in safety-critical applications
- Automated calculation pipelines
- User-facing applications
- Long-running simulations
- Systems with varying input quality

**May be overkill for:**
- Quick prototypes
- Well-controlled research calculations
- Educational examples
- Single-use scripts

## Additional Resources

See `handler.py` for complete implementation examples and `best-practices.md` for detailed guidance on production-ready fluid calculations.
