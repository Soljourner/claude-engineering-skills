# Best Practices for Error Handling in Fluid Calculations

Comprehensive guide for writing robust, production-ready fluid mechanics code with proper error handling, validation, and debugging capabilities.

## Table of Contents

1. [Input Validation Strategy](#input-validation-strategy)
2. [Physical Constraints](#physical-constraints)
3. [Numerical Stability](#numerical-stability)
4. [Error Recovery](#error-recovery)
5. [Logging and Monitoring](#logging-and-monitoring)
6. [Testing Edge Cases](#testing-edge-cases)
7. [Documentation Standards](#documentation-standards)
8. [Performance Considerations](#performance-considerations)

---

## Input Validation Strategy

### Always Validate at Entry Points

Validate all inputs at the beginning of functions before any calculations:

```python
def calculate_pressure_drop(Q, D, L, rho, mu, roughness):
    """
    Calculate pressure drop in pipe.

    VALIDATION FIRST - before any calculations
    """
    # Validate all inputs
    validate_positive(Q, 'flow_rate')
    validate_positive(D, 'diameter')
    validate_positive(L, 'length')
    validate_positive(rho, 'density')
    validate_positive(mu, 'viscosity')
    validate_positive(roughness, 'roughness', allow_zero=True)

    # Now safe to calculate
    V = 4 * Q / (math.pi * D**2)
    Re = rho * V * D / mu
    f = friction_factor(Re, roughness, D)
    dP = f * (L/D) * (rho * V**2 / 2)

    return dP
```

### Validate Ranges, Not Just Positivity

Check that inputs are in reasonable ranges, not just positive:

```python
def validate_pipe_inputs(D, roughness):
    """Comprehensive pipe geometry validation."""

    # Check positive
    if D <= 0:
        raise ValueError(f"Diameter {D} must be positive")

    # Check reasonable range
    if D < 0.001:  # 1 mm
        warnings.warn(f"Diameter {D} m is very small (< 1 mm)")
    elif D > 10:  # 10 m
        warnings.warn(f"Diameter {D} m is very large (> 10 m)")

    # Check roughness relative to diameter
    if roughness > 0.1 * D:
        warnings.warn(
            f"Roughness/diameter ratio {roughness/D:.3f} is very high (> 0.1)"
        )
```

### Use Type Hints and Runtime Checking

Combine type hints with runtime validation:

```python
def reynolds_number(velocity: float, diameter: float, viscosity: float) -> float:
    """
    Calculate Reynolds number with type hints and validation.

    Args:
        velocity: Flow velocity in m/s (must be positive)
        diameter: Characteristic length in m (must be positive)
        viscosity: Kinematic viscosity in m²/s (must be positive)

    Returns:
        Reynolds number (dimensionless, non-negative)
    """
    # Runtime validation even with type hints
    for name, value in [('velocity', velocity),
                         ('diameter', diameter),
                         ('viscosity', viscosity)]:
        if not isinstance(value, (int, float)):
            raise TypeError(f"{name} must be numeric, got {type(value)}")
        if value <= 0:
            raise ValueError(f"{name} must be positive, got {value}")

    return velocity * diameter / viscosity
```

---

## Physical Constraints

### Check Thermodynamic Consistency

Ensure properties are physically realizable:

```python
def validate_state_point(T, P, fluid='water'):
    """
    Validate that temperature and pressure define valid state.

    Checks:
    - Temperature above absolute zero
    - Pressure positive
    - State not in two-phase region (for liquid calculations)
    """
    # Temperature check
    validate_temperature(T, T_min=-273.15, T_max=1000)

    # Pressure check
    validate_pressure(P, P_min=0)

    # Check against vapor pressure
    if fluid == 'water':
        P_sat = saturation_pressure(T)

        # For liquid calculations, P must be > P_sat
        if P < P_sat:
            raise ValueError(
                f"Pressure {P/1000:.2f} kPa below saturation pressure "
                f"{P_sat/1000:.2f} kPa at {T}°C - two-phase region"
            )
```

### Cavitation Prevention

Always check NPSH (Net Positive Suction Head) for pumps:

```python
def check_npsh_available(P_suction, P_vapor, rho, V_suction):
    """
    Calculate and validate NPSH available.

    NPSH_available = (P_suction - P_vapor)/(rho*g) + V²/(2g)

    Must exceed NPSH_required (from pump curve) to prevent cavitation.
    """
    g = 9.81  # m/s²

    # Static head contribution
    h_static = (P_suction - P_vapor) / (rho * g)

    # Velocity head contribution
    h_velocity = V_suction**2 / (2 * g)

    NPSH_a = h_static + h_velocity

    if NPSH_a < 0:
        raise ValueError(
            f"NPSH available is negative ({NPSH_a:.2f} m) - "
            "severe cavitation risk"
        )

    if NPSH_a < 2:
        warnings.warn(
            f"NPSH available is low ({NPSH_a:.2f} m) - "
            "check against pump NPSH required"
        )

    return NPSH_a
```

### Compressibility Checks

Validate incompressible flow assumption:

```python
def validate_incompressible_assumption(V, a, location='flow'):
    """
    Check if incompressible flow assumption is valid.

    Mach number M = V/a should be < 0.3 for incompressible flow.
    """
    M = V / a

    if M > 0.3:
        warnings.warn(
            f"Mach number {M:.2f} at {location} exceeds 0.3. "
            "Compressibility effects significant. "
            "Incompressible flow assumption may be invalid."
        )
        return False

    return True
```

---

## Numerical Stability

### Safe Division Operations

Always check denominators:

```python
def safe_divide(numerator, denominator, epsilon=1e-15, param_name='denominator'):
    """
    Safely divide with zero checking.
    """
    if abs(denominator) < epsilon:
        raise ValueError(
            f"Division by near-zero value: {param_name} = {denominator:.2e}"
        )

    result = numerator / denominator

    # Check result
    if math.isnan(result):
        raise ValueError("Division resulted in NaN")
    if math.isinf(result):
        raise ValueError("Division resulted in infinity")

    return result


# Usage in calculations
def mach_number(V, a):
    """Calculate Mach number with safe division."""
    return safe_divide(V, a, param_name='speed_of_sound')
```

### Prevent Overflow/Underflow

Check intermediate calculations:

```python
def pressure_drop_pipe(f, L, D, rho, V):
    """
    Calculate pressure drop with overflow checking.

    ΔP = f * (L/D) * (ρV²/2)
    """
    # Check each term
    term1 = f * L / D
    if term1 > 1e6:
        warnings.warn(f"f*L/D = {term1:.2e} is very large")

    term2 = rho * V**2 / 2
    if term2 > 1e8:
        warnings.warn(f"ρV²/2 = {term2:.2e} is very large")

    dP = term1 * term2

    # Validate result
    if dP > 1e9:  # > 10,000 bar
        warnings.warn(f"Pressure drop {dP/1e6:.1f} MPa is extremely high")

    return dP
```

### Handle Square Root Arguments

Ensure arguments are non-negative:

```python
def velocity_from_pressure(dP, rho, K=1.0):
    """
    Calculate velocity from pressure drop: V = sqrt(2*ΔP/(K*ρ))

    Checks that argument to sqrt is non-negative.
    """
    if dP < 0:
        raise ValueError(f"Pressure drop {dP} Pa cannot be negative")

    if rho <= 0:
        raise ValueError(f"Density {rho} kg/m³ must be positive")

    if K <= 0:
        raise ValueError(f"Loss coefficient {K} must be positive")

    argument = 2 * dP / (K * rho)

    # Should be positive, but check for numerical issues
    if argument < 0:
        if abs(argument) < 1e-10:
            # Numerical noise - treat as zero
            argument = 0
        else:
            raise ValueError(f"sqrt argument {argument} is negative")

    V = math.sqrt(argument)
    return V
```

---

## Error Recovery

### Implement Fallback Methods

Provide alternative calculation methods:

```python
def friction_factor_with_fallbacks(Re, roughness, diameter):
    """
    Calculate friction factor with multiple fallback methods.

    Hierarchy:
    1. Colebrook-White (most accurate, iterative)
    2. Swamee-Jain (explicit approximation)
    3. Blasius (smooth pipe)
    4. Laminar formula (low Re)
    """
    # Try Colebrook-White
    try:
        return colebrook_iteration(Re, roughness, diameter)
    except ConvergenceError:
        logger.warning("Colebrook iteration failed, using Swamee-Jain")

    # Try Swamee-Jain
    try:
        return swamee_jain_formula(Re, roughness, diameter)
    except ValueError:
        logger.warning("Swamee-Jain failed, using Blasius")

    # Try Blasius (smooth pipe)
    try:
        if Re > 4000 and Re < 100000:
            return 0.316 / Re**0.25
    except:
        pass

    # Last resort: laminar formula (always works)
    logger.warning(f"All turbulent methods failed for Re={Re}, using laminar formula")
    return 64 / Re
```

### Graceful Degradation

Provide reduced functionality rather than complete failure:

```python
def calculate_system_curve(Q_range, system_params, detailed=True):
    """
    Calculate system curve with graceful degradation.

    If detailed calculation fails, fall back to simplified model.
    """
    results = []

    for Q in Q_range:
        try:
            # Try full calculation with all minor losses
            H = detailed_head_loss(Q, system_params)
            results.append({'Q': Q, 'H': H, 'method': 'detailed'})

        except Exception as e:
            if detailed:
                logger.warning(f"Detailed calculation failed at Q={Q}: {e}")
                logger.warning("Falling back to simplified model")

            try:
                # Simplified model: only major losses
                H = simplified_head_loss(Q, system_params)
                results.append({'Q': Q, 'H': H, 'method': 'simplified'})

            except Exception as e2:
                logger.error(f"Both methods failed at Q={Q}: {e2}")
                # Include NaN rather than breaking entire calculation
                results.append({'Q': Q, 'H': float('nan'), 'method': 'failed'})

    return results
```

### Iterative Solver Safeguards

Prevent infinite loops and divergence:

```python
def solve_implicit_equation(initial_guess, target, max_iter=100, tol=1e-6):
    """
    Solve implicit equation with comprehensive safeguards.
    """
    x = initial_guess
    history = [x]

    for iteration in range(max_iter):
        x_old = x

        try:
            # Newton-Raphson update
            f = implicit_function(x, target)
            df = derivative(x, target)

            # Check derivative
            if abs(df) < 1e-15:
                raise ValueError("Derivative near zero - Newton-Raphson unstable")

            x = x_old - f / df

        except Exception as e:
            logger.error(f"Iteration error at step {iteration}: {e}")
            raise

        history.append(x)

        # Check convergence
        error = abs(x - x_old)
        if error < tol:
            logger.debug(f"Converged in {iteration} iterations")
            return x, history

        # Check for oscillation
        if iteration > 5:
            if abs(x - history[-3]) < tol:
                logger.warning("Oscillating - using average")
                return (x + history[-3]) / 2, history

        # Check for divergence
        if abs(x) > 1e10:
            raise RuntimeError(f"Solution diverging: x = {x:.2e}")

        # Check for NaN
        if math.isnan(x):
            raise RuntimeError("Solution became NaN")

    raise RuntimeError(f"Failed to converge after {max_iter} iterations")
```

---

## Logging and Monitoring

### Structured Logging

Use consistent logging levels:

```python
import logging

logger = logging.getLogger(__name__)

def calculate_pump_performance(Q, N, D):
    """
    Calculate pump performance with structured logging.
    """
    # DEBUG: Detailed calculation steps
    logger.debug(f"Calculating pump performance: Q={Q}, N={N}, D={D}")

    # INFO: Important results
    H = calculate_head(Q, N, D)
    logger.info(f"Pump head: {H:.2f} m at Q={Q:.3f} m³/s")

    # WARNING: Concerning conditions
    if Q < 0.5 * Q_nominal:
        logger.warning(f"Flow rate {Q} below 50% of nominal - risk of instability")

    # ERROR: Calculation failures
    try:
        eta = calculate_efficiency(Q, N, D)
    except Exception as e:
        logger.error(f"Efficiency calculation failed: {e}", exc_info=True)
        eta = None

    # CRITICAL: System-threatening conditions
    if H > H_max_allowable:
        logger.critical(f"Head {H} exceeds maximum {H_max_allowable} - shutdown required")
        raise CriticalPressureError(f"Excessive head: {H} m")

    return {'H': H, 'eta': eta}
```

### Performance Monitoring

Track calculation performance:

```python
import time
from functools import wraps

def monitor_performance(func):
    """
    Decorator to monitor function performance.
    """
    @wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()

        try:
            result = func(*args, **kwargs)
            elapsed = time.time() - start_time

            # Log slow calculations
            if elapsed > 1.0:
                logger.warning(
                    f"{func.__name__} took {elapsed:.2f}s - consider optimization"
                )
            else:
                logger.debug(f"{func.__name__} completed in {elapsed:.4f}s")

            return result

        except Exception as e:
            elapsed = time.time() - start_time
            logger.error(
                f"{func.__name__} failed after {elapsed:.2f}s: {e}"
            )
            raise

    return wrapper


@monitor_performance
def expensive_calculation(params):
    """Calculation with automatic performance monitoring."""
    # ... complex calculation ...
    pass
```

### Calculation Auditing

Create audit trail for critical calculations:

```python
class CalculationAuditor:
    """
    Audit trail for fluid calculations.
    """
    def __init__(self, log_file='calculation_audit.log'):
        self.log_file = log_file

    def log_calculation(self, calc_type, inputs, outputs, metadata=None):
        """
        Log calculation with full context.
        """
        import json
        import datetime

        record = {
            'timestamp': datetime.datetime.now().isoformat(),
            'calculation': calc_type,
            'inputs': inputs,
            'outputs': outputs,
            'metadata': metadata or {}
        }

        with open(self.log_file, 'a') as f:
            f.write(json.dumps(record) + '\n')


# Usage
auditor = CalculationAuditor()

def audited_calculation(Q, H):
    """Calculation with audit trail."""
    result = calculate_power(Q, H)

    auditor.log_calculation(
        calc_type='pump_power',
        inputs={'Q': Q, 'H': H},
        outputs={'P': result},
        metadata={'units': 'SI', 'version': '1.0'}
    )

    return result
```

---

## Testing Edge Cases

### Comprehensive Test Suite

Test boundary conditions:

```python
import pytest

class TestReynoldsNumber:
    """Test suite for Reynolds number calculation."""

    def test_normal_values(self):
        """Test with typical values."""
        Re = reynolds_number(V=1.0, D=0.1, nu=1e-6)
        assert 90000 < Re < 110000

    def test_zero_velocity(self):
        """Zero velocity should give zero Reynolds number."""
        Re = reynolds_number(V=0, D=0.1, nu=1e-6)
        assert Re == 0

    def test_negative_velocity(self):
        """Negative velocity should raise error."""
        with pytest.raises(ValueError):
            reynolds_number(V=-1.0, D=0.1, nu=1e-6)

    def test_zero_viscosity(self):
        """Zero viscosity should raise error (division by zero)."""
        with pytest.raises(ValueError):
            reynolds_number(V=1.0, D=0.1, nu=0)

    def test_very_small_viscosity(self):
        """Very small viscosity should give warning."""
        with pytest.warns(UserWarning):
            reynolds_number(V=1.0, D=0.1, nu=1e-20)

    def test_very_large_result(self):
        """Very large Reynolds number should not cause overflow."""
        Re = reynolds_number(V=100, D=10, nu=1e-10)
        assert not math.isinf(Re)
        assert not math.isnan(Re)
```

### Property-Based Testing

Use hypothesis for property-based testing:

```python
from hypothesis import given, strategies as st

@given(
    V=st.floats(min_value=0.01, max_value=100),
    D=st.floats(min_value=0.001, max_value=10),
    nu=st.floats(min_value=1e-8, max_value=1e-3)
)
def test_reynolds_properties(V, D, nu):
    """
    Test Reynolds number properties with random valid inputs.

    Properties:
    - Result is always positive
    - Result is not NaN or Inf
    - Doubling velocity doubles Reynolds number
    """
    Re = reynolds_number(V, D, nu)

    # Should be positive
    assert Re > 0

    # Should be finite
    assert not math.isnan(Re)
    assert not math.isinf(Re)

    # Should scale linearly with velocity
    Re2 = reynolds_number(2*V, D, nu)
    assert abs(Re2 / Re - 2.0) < 1e-10
```

---

## Documentation Standards

### Docstring Requirements

Every function must document:

```python
def calculate_pressure_drop(Q, D, L, rho, mu, roughness):
    """
    Calculate pressure drop in circular pipe using Darcy-Weisbach equation.

    Uses Colebrook-White equation for friction factor with Swamee-Jain
    initial guess for iteration.

    Args:
        Q (float): Volumetric flow rate in m³/s (must be positive)
        D (float): Pipe diameter in m (must be positive)
        L (float): Pipe length in m (must be positive)
        rho (float): Fluid density in kg/m³ (must be positive)
        mu (float): Dynamic viscosity in Pa·s (must be positive)
        roughness (float): Absolute roughness in m (non-negative)

    Returns:
        float: Pressure drop in Pa (always positive for forward flow)

    Raises:
        ValueError: If any input is invalid (negative, zero, NaN, or Inf)
        RuntimeError: If friction factor iteration fails to converge

    Notes:
        - Valid for fully developed turbulent flow (Re > 4000)
        - Assumes incompressible flow (M < 0.3)
        - Uses Darcy friction factor (f_Darcy = 4 * f_Fanning)

    Examples:
        >>> # Water flow in steel pipe
        >>> dP = calculate_pressure_drop(
        ...     Q=0.01, D=0.1, L=100,
        ...     rho=1000, mu=0.001, roughness=0.045e-3
        ... )
        >>> print(f"Pressure drop: {dP/1000:.2f} kPa")
        Pressure drop: 23.45 kPa

    References:
        - Colebrook, C.F. (1939). "Turbulent flow in pipes..."
        - Moody, L.F. (1944). "Friction factors for pipe flow"
    """
    # Implementation...
```

### Error Message Quality

Provide actionable error messages:

```python
# BAD: Vague error message
if D <= 0:
    raise ValueError("Invalid diameter")

# GOOD: Specific, actionable message
if D <= 0:
    raise ValueError(
        f"Pipe diameter must be positive, got D = {D} m. "
        "Check input units (should be meters, not millimeters)."
    )

# BETTER: Include context and suggestions
if D <= 0:
    raise ValueError(
        f"Pipe diameter must be positive, got D = {D} m.\n"
        "Possible causes:\n"
        "  - Input in wrong units (mm instead of m?)\n"
        "  - Negative value from previous calculation\n"
        "  - Uninitialized variable\n"
        "Check input source and unit conversions."
    )
```

---

## Performance Considerations

### Balance Safety and Speed

```python
def reynolds_number_fast(V, D, nu):
    """
    Fast Reynolds number - minimal validation for hot paths.

    Use when:
    - Inputs already validated
    - Called in tight loops
    - Performance critical

    WARNING: No input validation! Use with caution.
    """
    return V * D / nu


def reynolds_number_safe(V, D, nu):
    """
    Safe Reynolds number - full validation.

    Use when:
    - Inputs from user or external source
    - One-off calculations
    - Safety critical
    """
    validate_positive(V, 'velocity')
    validate_positive(D, 'diameter')
    validate_positive(nu, 'viscosity')

    Re = V * D / nu
    validate_reynolds_number(Re)

    return Re
```

### Conditional Validation

Use debug mode for expensive checks:

```python
DEBUG = False  # Set True during development

def calculate_with_validation(inputs, debug=DEBUG):
    """
    Perform calculation with optional expensive validation.
    """
    # Always do cheap validation
    validate_positive(inputs['velocity'], 'velocity')

    # Only do expensive validation in debug mode
    if debug:
        validate_dimensional_consistency(inputs)
        check_numerical_stability_all(inputs)
        verify_against_correlations(inputs)

    # Perform calculation
    result = calculate(inputs)

    return result
```

---

## Summary Checklist

For every production fluid calculation function:

- [ ] Input validation at entry point
- [ ] Physical constraint checking
- [ ] Safe division (denominator checks)
- [ ] NaN/Inf checking in results
- [ ] Appropriate error messages
- [ ] Logging at appropriate levels
- [ ] Try-except blocks for risky operations
- [ ] Fallback methods where applicable
- [ ] Comprehensive docstring
- [ ] Unit tests including edge cases
- [ ] Performance profiling if in hot path
- [ ] Code review by peer

---

## Additional Resources

- See `handler.py` for complete implementation examples
- See `SKILL.md` for error catalog and strategies
- NIST guidelines for numerical accuracy
- IEEE 754 floating-point standard
- Python logging documentation
