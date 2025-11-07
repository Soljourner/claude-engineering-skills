# Claude Engineering Skills Library

A comprehensive collection of Agent Skills for mechanical engineering, aerospace engineering, and pump design workflows. This library equips Claude with specialized knowledge, tools, and methodologies for professional engineering analysis, design, and simulation.

## Overview

The Claude Engineering Skills Library provides 100+ skills organized across five categories:

- **Databases** (10 skills) – Access to fluid properties, material data, pump performance curves, and thermodynamic tables
- **Packages** (15 skills) – Python libraries for numerical computation, optimization, visualization, and fluid mechanics
- **Integrations** (7 skills) – Connections to CAD and simulation software (ANSYS, OpenFOAM, SolidWorks, COMSOL)
- **Helpers** (6 skills) – Utilities for unit conversion, pump selection, property calculation, and error handling
- **Thinking** (9 skills) – Structured workflows for fluid dynamics, pump design, thermodynamics, and structural analysis

## Quick Start

### Installation

1. Clone this repository:
```bash
git clone https://github.com/Soljourner/claude-engineering-skills.git
cd claude-engineering-skills
```

2. Install recommended Python packages:
```bash
pip install numpy scipy sympy matplotlib pint fluids thermo CoolProp pyvista networkx
```

3. Skills are automatically available to Claude when this repository is in your workspace.

### Using Skills

Skills follow a progressive disclosure pattern:
- **Level 1**: Claude sees skill names and descriptions automatically
- **Level 2**: Claude loads full `SKILL.md` when relevant to your task
- **Level 3**: Additional resources loaded on-demand (examples, references, code)

Simply describe your engineering task to Claude, and relevant skills will be activated automatically.

## Example Workflows

### Pump Performance Analysis
```
Analyze a centrifugal pump design for water at 20°C with a flow rate of 100 m³/h
and 50m head. Calculate efficiency, specific speed, and check for cavitation risk.
```

Claude will activate:
- `coolprop-db` or `nist-refprop` for water properties
- `pump-design/centrifugal-pumps` for design equations
- `pump-design/cavitation-analysis` for NPSH calculations
- `fluids-package` for Reynolds number and friction
- `matplotlib-visualization` for performance curves

### Fluid Dynamics Simulation Setup
```
Set up an OpenFOAM simulation for turbulent flow through a pipe bend with
water at 5 m/s inlet velocity. Use k-omega SST turbulence model.
```

Claude will activate:
- `openfoam-cfd` for case setup
- `thinking/fluid-dynamics` for boundary conditions
- `turbulence-models-db` for k-omega parameters
- `coolprop-package` for water properties

### System Design Optimization
```
Design a multi-pump system to deliver 500 m³/h against 80m total head.
Minimize energy consumption and cost.
```

Claude will activate:
- `pump-selection-helper` for pump type recommendations
- `scipy-optimization` for multi-objective optimization
- `hydraulic-components-db` for piping losses
- `pump-design/system-integration` for network modeling

## Categories Summary

### Databases (`skills/databases/`)
Access to engineering data essential for accurate analysis:
- **nist-refprop** / **coolprop-db**: Thermodynamic and transport properties of 100+ fluids
- **pump-performance-db**: Manufacturer pump curves (Grundfos, KSB, Flowserve)
- **material-properties-db**: Viscosity, density, thermal properties vs temperature
- **cavitation-risk-db**: Vapor pressure and NPSH requirements
- **turbulence-models-db**: Parameters for CFD turbulence models
- **hydraulic-components-db**: Loss coefficients for pipes, valves, fittings

### Packages (`skills/packages/`)
Python tools for engineering computation:
- **fluids-package**: Pipe flow, pump sizing, compressible flow calculations
- **thermo-package**: Thermodynamic properties, mixture calculations
- **scipy-optimization**: Design optimization, curve fitting
- **pint-units**: Consistent unit handling and conversion
- **matplotlib-visualization**: Performance curves, contour plots
- **networkx-flow-networks**: Hydraulic network analysis

### Integrations (`skills/integrations/`)
Connect with professional engineering software:
- **openfoam-cfd**: Open-source CFD setup and automation
- **ansys-simulation**: ANSYS Fluent and Workbench scripting
- **solidworks-cad**: Parametric pump design automation
- **comsol-multiphysics**: Coupled fluid-structural analysis

### Helpers (`skills/helpers/`)
Utilities for common engineering tasks:
- **unit-converter**: Flow rates, pressures, viscosities
- **pump-selection-helper**: Decision tree for pump type selection
- **fluid-property-calculator**: Quick calculations without database queries
- **engineering-context-init**: Session setup with standard constants

### Thinking (`skills/thinking/`)
Structured methodologies for engineering analysis:
- **fluid-dynamics**: Systematic CFD workflow from setup to post-processing
- **pump-design**: Complete pump design from requirements to testing
- **thermodynamics**: Cycle analysis and heat transfer calculations
- **structural-analysis**: FEA workflows for mechanical components

## Setup Instructions

### Database Access

Some database skills require authentication or installation:

**NIST REFPROP** (commercial):
```bash
# Purchase license from NIST
# Install REFPROP and Python wrapper
pip install ctREFPROP
export RPPREFIX=/path/to/refprop
```

**CoolProp** (free alternative):
```bash
pip install CoolProp
```

**NASA Earthdata**:
```bash
# Create free account at earthdata.nasa.gov
# Set credentials:
export EARTHDATA_USERNAME=your_username
export EARTHDATA_PASSWORD=your_password
```

### Simulation Software Integration

**OpenFOAM**:
```bash
# Install OpenFOAM (Ubuntu/Debian)
sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get update
sudo apt-get install openfoam10
source /opt/openfoam10/etc/bashrc
```

**ANSYS, COMSOL, SolidWorks**: Require commercial licenses. Integration skills provide Python/API scripting guidance.

## FAQs

**Q: Do I need all packages installed to use the library?**
A: No. Install only the packages relevant to your workflows. Claude will guide you if a required package is missing.

**Q: How accurate are the fluid property calculations?**
A: Database skills use validated sources (NIST, industry standards). All formulas include verification tests with known data.

**Q: Can I add custom pump data?**
A: Yes. See `skills/databases/pump-performance-db/SKILL.md` for instructions on adding manufacturer data.

**Q: What units are used?**
A: SI units by default. The `pint-units` package handles conversions. US customary units supported for legacy applications.

**Q: Are the turbulence models validated?**
A: Yes. `turbulence-models-db` includes validation cases from literature (Wilcox, Pope, Versteeg & Malalasekera).

## Engineering Domains Covered

### Mechanical Engineering
- Pump design (centrifugal, positive displacement, axial flow)
- Piping systems and hydraulic networks
- Fluid machinery efficiency and optimization
- Cavitation analysis and prevention
- Mechanical component stress analysis

### Fluid Dynamics
- Incompressible internal flows (pipes, pumps)
- Compressible flows (high-speed applications)
- Turbulence modeling (RANS, LES, DNS)
- Multiphase flows (cavitation, bubbly flows)
- Boundary layer analysis

### Aerospace Engineering
- High-speed fluid flows
- Turbopumps for rocket propulsion
- Atmospheric property modeling
- Orbital mechanics with fluid thrust

### Thermodynamics
- Rankine and Brayton cycles
- Heat transfer in fluid systems
- Isentropic efficiency calculations
- Phase equilibria

## Contributing

We welcome contributions from the engineering community! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Adding new skills
- Improving existing documentation
- Submitting verification tests
- Reporting issues

## License

This project is licensed under the MIT License - see [LICENSE.md](LICENSE.md) for details.

## Citation

If you use this skills library in your research or professional work, please cite:

```bibtex
@software{claude_engineering_skills,
  title = {Claude Engineering Skills Library},
  author = {Soljourner},
  year = {2025},
  url = {https://github.com/Soljourner/claude-engineering-skills}
}
```

## Resources

- [Anthropic Agent Skills Documentation](https://www.anthropic.com/engineering/equipping-agents-for-the-real-world-with-agent-skills)
- [Claude Code Documentation](https://docs.claude.com/claude-code)
- [OpenFOAM User Guide](https://www.openfoam.com/documentation)
- [Python Fluids Package](https://fluids.readthedocs.io/)
- [NIST REFPROP Database](https://www.nist.gov/srd/refprop)

## Support

For issues, questions, or feature requests:
- Open an issue on GitHub
- Review existing skills for examples
- Consult the `docs/` directory for detailed guides

---

Built with precision for professional mechanical and aerospace engineering workflows.
