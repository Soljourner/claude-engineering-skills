# Pump Performance Curves Skill

This skill provides comprehensive tools for generating, interpreting, and manipulating pump performance curves.

## Files

- **SKILL.md**: Main skill documentation with theory, methods, and guidelines
- **plotter.py**: Python module with verified examples and computational tools
- **reference.md**: Reference documentation with equations, standards, and nomenclature

## Quick Start

```python
from plotter import generate_typical_pump_curve, plot_pump_curves

# Generate typical pump curves
pump = generate_typical_pump_curve(
    Q_design=100,    # m³/h
    H_design=50,     # m
    specific_speed=35,
    speed_rpm=1750
)

# Plot all performance curves
plot_pump_curves(pump)
```

## Key Capabilities

1. **Performance Curves**: H-Q, η-Q, P-Q, NPSH-Q generation
2. **Affinity Laws**: Speed and diameter changes
3. **Curve Fitting**: Fit curves to test data
4. **Operating Point**: System curve intersection analysis
5. **Standards Compliance**: ISO 9906, API 610, ANSI/HI

## Dependencies

- numpy
- matplotlib  
- scipy

## Verified Examples

Run all examples:
```bash
python plotter.py
```

Individual examples:
- Example 1: Typical centrifugal pump curves
- Example 2: Affinity laws - speed change
- Example 3: Impeller trimming
- Example 4: Operating point determination
- Example 5: Curve fitting from test data
