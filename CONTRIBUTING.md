# Contributing to Claude Engineering Skills Library

We welcome contributions from mechanical engineers, aerospace engineers, fluid dynamicists, and software developers. This guide ensures consistency, accuracy, and quality across the skills library.

## Table of Contents
- [Getting Started](#getting-started)
- [Types of Contributions](#types-of-contributions)
- [Skill Development Guidelines](#skill-development-guidelines)
- [Verification Requirements](#verification-requirements)
- [Documentation Standards](#documentation-standards)
- [Code Quality](#code-quality)
- [Submission Process](#submission-process)

## Getting Started

### Prerequisites
- Engineering background (mechanical, aerospace, or related field)
- Python programming skills (for code examples)
- Familiarity with Git and GitHub
- Understanding of the domain you're contributing to

### Setting Up Development Environment

1. Fork the repository
2. Clone your fork:
```bash
git clone https://github.com/YOUR_USERNAME/claude-engineering-skills.git
cd claude-engineering-skills
```

3. Create a development branch:
```bash
git checkout -b feature/your-skill-name
```

4. Install development dependencies:
```bash
pip install numpy scipy sympy matplotlib pint fluids thermo CoolProp pytest black flake8
```

## Types of Contributions

### 1. New Skills
Add entirely new skills for databases, packages, integrations, helpers, or thinking workflows.

**Focus Areas We Need:**
- Advanced turbulence modeling techniques
- Fatigue analysis for pump components
- Noise and vibration prediction
- Multiphase flow modeling beyond cavitation
- Control system integration for pump systems
- Energy efficiency optimization algorithms
- Material selection databases for corrosive fluids

### 2. Skill Enhancements
Improve existing skills with:
- Additional examples
- Better error handling
- Extended reference materials
- Performance optimizations
- Coverage of edge cases

### 3. Verification Tests
Add test cases that validate formulas and calculations against:
- Published literature results
- Standard test problems
- Manufacturer data
- Experimental measurements

### 4. Documentation Improvements
- Fix typos and clarify explanations
- Add diagrams and illustrations
- Improve example workflows
- Update references to latest standards

### 5. Bug Fixes
Report and fix:
- Incorrect formulas
- Units inconsistencies
- Code errors
- Documentation inaccuracies

## Skill Development Guidelines

### Skill Structure

Every skill must follow this directory structure:

```
skills/[category]/[skill-name]/
├── SKILL.md           # Core instructions (REQUIRED)
├── reference.md       # Additional references (OPTIONAL)
├── examples.py        # Code examples (RECOMMENDED)
├── verification/      # Test cases (REQUIRED for calculations)
│   ├── test_cases.py
│   └── validation_data.json
└── resources/         # Additional materials (OPTIONAL)
    ├── diagrams/
    └── tables/
```

### SKILL.md Format (REQUIRED)

All `SKILL.md` files must start with YAML frontmatter and follow this structure:

```markdown
---
name: skill-name
description: Brief one-line description for Claude's activation decision (max 100 chars)
category: databases|packages|integrations|helpers|thinking
domain: mechanical|aerospace|fluids|thermodynamics|structural
complexity: basic|intermediate|advanced
dependencies: [list of required packages]
---

# Skill Name

## Overview
Brief explanation of what this skill does and when to use it (2-3 sentences).

## When to Use This Skill
- Specific scenario 1
- Specific scenario 2
- Specific scenario 3

## Prerequisites / Setup
Installation requirements, API keys, authentication, environment variables.

## Core Functions / Methods
Detailed documentation of available tools, functions, or workflows.

### Function/Concept 1
- **Purpose**: What it does
- **Parameters**: Input requirements
- **Returns**: Output description
- **Example**:
```python
# Code example with comments
```

## Example Workflows

### Workflow 1: Descriptive Title
Step-by-step example with complete code and explanation.

```python
# Complete, runnable example
import necessary_packages

# Step 1: Setup
...

# Step 2: Calculations
...

# Step 3: Results
...
```

**Expected Output:**
```
[Show what the output should look like]
```

## Engineering Considerations
Important notes about:
- Assumptions and limitations
- Validity ranges (Reynolds numbers, Mach numbers, etc.)
- Numerical stability
- Error propagation
- Safety factors

## Common Pitfalls
- Mistake 1 and how to avoid it
- Mistake 2 and how to avoid it

## Verification and Validation
How to verify results:
- Cross-check methods
- Comparison with experimental data
- Standard test problems
- Order-of-magnitude checks

## References
- [Author, Year] Standard/Paper Title, DOI or URL
- [Organization] Design Manual/Standard, Edition
- [Software] Documentation URL

## Related Skills
- `related-skill-1` - Brief description of relationship
- `related-skill-2` - Brief description of relationship
```

### Description Guidelines

The `description` field in YAML frontmatter is CRITICAL. This is what Claude uses to decide whether to activate your skill.

**Good descriptions:**
- ✅ "Calculate NPSH and assess cavitation risk in centrifugal pumps"
- ✅ "Query thermodynamic properties for 100+ fluids from CoolProp database"
- ✅ "Set up turbulent pipe flow simulations in OpenFOAM with k-omega SST"

**Bad descriptions:**
- ❌ "Pump cavitation" (too vague)
- ❌ "This skill helps with thermodynamic property calculations using various databases" (too long)
- ❌ "OpenFOAM" (not descriptive of capability)

## Verification Requirements

### For All Calculations

**CRITICAL**: Every formula or calculation must include verification against known data.

Create `verification/test_cases.py`:

```python
import pytest
import numpy as np
from your_skill_module import calculation_function

def test_formula_name_against_textbook():
    """
    Verify formula against Example X.Y from [Author, Year].

    Reference:
        Author, "Book Title", Year, p. XXX
        DOI: xxxxx (if available)
    """
    # Given conditions from reference
    input_param1 = 1000.0  # units
    input_param2 = 298.15  # K

    # Known correct result from reference
    expected_result = 42.5  # units
    expected_tolerance = 0.1  # acceptable error

    # Calculate using our implementation
    result = calculation_function(input_param1, input_param2)

    # Verify
    assert np.isclose(result, expected_result, rtol=expected_tolerance), \
        f"Expected {expected_result}, got {result}"

def test_edge_case_zero_flow():
    """Verify behavior at zero flow rate (pump shutoff head)."""
    # Edge case testing
    ...

def test_unit_consistency():
    """Verify unit conversions maintain physical consistency."""
    # Unit testing
    ...
```

### Validation Data Sources

Acceptable sources for verification (in order of preference):
1. **Primary**: Peer-reviewed journal articles with DOI
2. **Primary**: Engineering textbooks (cite edition and page)
3. **Primary**: Standards (ASME, ANSI, ISO, ASTM)
4. **Secondary**: Manufacturer technical documentation
5. **Secondary**: Government databases (NIST, NASA)
6. **Tertiary**: Reputable online calculators (document source)

**NEVER use:**
- Wikipedia as a primary source (acceptable for initial exploration only)
- Random websites without clear authorship
- Your own calculations without external verification
- "Common knowledge" without citation

### Running Verification Tests

Before submitting:
```bash
# Run all tests in your skill
pytest skills/category/your-skill/verification/

# Run with coverage
pytest --cov=skills/category/your-skill skills/category/your-skill/verification/
```

All tests must pass with 100% success rate.

## Documentation Standards

### Code Examples

- Must be complete and runnable
- Include all necessary imports
- Use realistic engineering values
- Add comments explaining non-obvious steps
- Show expected output
- Handle errors gracefully

### Units

- Default to SI units (m, kg, s, Pa, K)
- Clearly state units for all quantities
- Use `pint` library for unit-aware calculations
- Provide conversion examples for US customary units where relevant

### Mathematical Notation

Use LaTeX notation in markdown for equations:

```markdown
The Bernoulli equation for incompressible flow:

$$P_1 + \frac{1}{2}\rho v_1^2 + \rho g z_1 = P_2 + \frac{1}{2}\rho v_2^2 + \rho g z_2$$

Where:
- $P$ = pressure (Pa)
- $\rho$ = density (kg/m³)
- $v$ = velocity (m/s)
- $g$ = gravitational acceleration (9.81 m/s²)
- $z$ = elevation (m)
```

### Diagrams

- Use ASCII art for simple diagrams
- Provide SVG or PNG for complex illustrations
- Always include source files (e.g., .drawio, .svg source)
- Ensure diagrams are colorblind-friendly

## Code Quality

### Python Style

Follow PEP 8 with these specifics:
- Maximum line length: 88 characters (Black formatter)
- Use type hints for function signatures
- Docstrings in NumPy format

```python
def calculate_reynolds_number(
    velocity: float,
    diameter: float,
    kinematic_viscosity: float
) -> float:
    """
    Calculate Reynolds number for pipe flow.

    Parameters
    ----------
    velocity : float
        Flow velocity (m/s)
    diameter : float
        Pipe diameter (m)
    kinematic_viscosity : float
        Kinematic viscosity (m²/s)

    Returns
    -------
    float
        Reynolds number (dimensionless)

    Examples
    --------
    >>> calculate_reynolds_number(2.0, 0.1, 1e-6)
    200000.0

    References
    ----------
    .. [1] White, F.M. "Fluid Mechanics", 8th ed., 2016, p. 325
    """
    return (velocity * diameter) / kinematic_viscosity
```

### Error Handling

Always validate inputs:

```python
def calculate_pump_power(flow_rate: float, head: float, efficiency: float) -> float:
    """Calculate pump power requirement."""

    # Input validation
    if flow_rate < 0:
        raise ValueError("Flow rate must be non-negative")
    if head < 0:
        raise ValueError("Head must be non-negative")
    if not 0 < efficiency <= 1:
        raise ValueError("Efficiency must be between 0 and 1")

    # Constants
    rho = 1000  # kg/m³ for water
    g = 9.81    # m/s²

    # Calculate power
    power = (rho * g * flow_rate * head) / efficiency

    return power
```

## Submission Process

### 1. Pre-Submission Checklist

- [ ] SKILL.md includes complete YAML frontmatter
- [ ] All formulas have verification tests that pass
- [ ] Code examples are complete and runnable
- [ ] References include DOI/ISBN when available
- [ ] No hardcoded API keys or credentials
- [ ] Code follows PEP 8 (run `black` and `flake8`)
- [ ] All units clearly specified
- [ ] Engineering assumptions documented
- [ ] Related skills cross-referenced

### 2. Testing Your Skill

```bash
# Format code
black skills/category/your-skill/

# Check style
flake8 skills/category/your-skill/

# Run verification tests
pytest skills/category/your-skill/verification/ -v

# Test with Claude (if possible)
# Describe a task that should activate your skill and verify it works
```

### 3. Commit Guidelines

Use conventional commit format:

```
feat(pumps): add cavitation analysis skill
fix(fluids): correct Reynolds number calculation
docs(thermodynamics): improve Rankine cycle examples
test(coolprop): add verification against NIST data
```

### 4. Pull Request

1. Push to your fork
2. Create pull request with:
   - Clear title describing the contribution
   - Description of what was added/changed
   - Results of verification tests
   - Any breaking changes or dependencies
   - Screenshots/output examples if applicable

3. PR template:
```markdown
## Description
[Brief description of changes]

## Type of Change
- [ ] New skill
- [ ] Skill enhancement
- [ ] Bug fix
- [ ] Documentation improvement
- [ ] Verification tests

## Verification
- [ ] All formulas verified against published sources
- [ ] Verification tests pass (attach test output)
- [ ] Code examples tested and runnable
- [ ] Documentation reviewed for accuracy

## References
[List key references used]

## Related Issues
Closes #[issue number]

## Additional Notes
[Any special considerations]
```

### 5. Review Process

Contributors will review:
1. **Technical accuracy**: Are formulas correct?
2. **Verification**: Are tests comprehensive and passing?
3. **Documentation quality**: Is it clear and complete?
4. **Code quality**: Does it follow standards?
5. **Usefulness**: Does it fill a genuine gap?

Expect feedback and iteration. We prioritize accuracy over speed.

## Engineering Standards Compliance

When contributing skills related to regulated domains, ensure compliance with relevant standards:

### Fluid Systems
- ASME B31.1/B31.3 - Piping codes
- API 610 - Centrifugal pumps for petroleum
- ISO 5199 - Centrifugal pumps technical requirements
- HI 9.6.3 - Pump piping considerations

### Pressure Vessels
- ASME Section VIII - Pressure vessel code

### CFD and Simulation
- ASME V&V 20 - Verification and validation in CFD
- AIAA Guidelines for turbulence modeling

### Units and Measurements
- ISO 80000 - International System of Quantities
- NIST SP 811 - Guide for metric usage

## Questions?

- Open a GitHub issue with the `question` label
- Review existing skills for examples
- Consult the `docs/` directory

## Code of Conduct

- Be respectful and constructive in feedback
- Focus on technical accuracy and engineering rigor
- Credit others' work appropriately
- Assume good faith in contributions
- Prioritize safety in engineering calculations

## Recognition

Contributors will be:
- Listed in repository contributors
- Acknowledged in relevant skill documentation
- Credited in any publications using the library (if applicable)

---

Thank you for contributing to the engineering community's AI capabilities!
