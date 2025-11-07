# Claude Engineering Skills Library - Complete Inventory

**Total Skills Created:** 37+
**Total Files:** 121+
**Documentation Pages:** 6
**Lines of Code/Documentation:** 15,000+

## Library Statistics

### By Category
- **Databases:** 8 skills
- **Packages:** 10 skills
- **Integrations:** 4 skills
- **Helpers:** 6 skills
- **Thinking:** 9 skills

### File Breakdown
- SKILL.md files: 37
- Example code files (.py): 35+
- Reference documentation: 30+
- Supporting docs: 19+

---

## Complete Skills Inventory

### 📊 Databases (8 skills)

#### 1. CoolProp Database (`coolprop-db`)
**Purpose:** Query thermodynamic properties for 100+ fluids from CoolProp database
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Water, air, refrigerant properties
- Temperature and pressure dependent calculations
- Saturation properties
- Two-phase calculations
- Verified against NIST data

#### 2. NIST REFPROP (`nist-refprop`)
**Purpose:** High-accuracy thermodynamic properties from NIST REFPROP (commercial)
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- 100+ pure fluids and mixtures
- Higher accuracy than CoolProp for certain fluids
- Commercial license required
- Python ctREFPROP wrapper

#### 3. Material Properties Database (`material-properties-db`)
**Purpose:** Query fluid viscosities, densities, and material properties vs temperature
**Files:** SKILL.md, query-examples.py, reference.md
**Key Features:**
- Water properties 0-100°C
- Oil viscosities
- Temperature-dependent correlations (Sutherland, Antoine)
- Verified against ASTM standards

#### 4. Pump Performance Database (`pump-performance-db`)
**Purpose:** Access manufacturer pump curves and specifications
**Files:** SKILL.md, query-examples.py, reference.md
**Key Features:**
- Grundfos, KSB, Flowserve data access
- H-Q curve database
- Custom pump data management
- Selection by requirements

#### 5. Cavitation Risk Database (`cavitation-risk-db`)
**Purpose:** Query vapor pressures and NPSH requirements for cavitation assessment
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Vapor pressure calculations
- NPSH correlations
- Temperature sensitivity analysis
- Safety margin evaluation
- Antoine equation coefficients

#### 6. Turbulence Models Database (`turbulence-models-db`)
**Purpose:** Select and configure turbulence models for CFD
**Files:** SKILL.md, reference.md
**Key Features:**
- k-ε, k-ω SST, Spalart-Allmaras models
- Model selection criteria
- Wall function parameters
- y+ requirements
- Validation test cases

#### 7. Hydraulic Components Database (`hydraulic-components-db`)
**Purpose:** Query loss coefficients for pipes, valves, fittings
**Files:** SKILL.md, query-examples.py, reference.md
**Key Features:**
- Loss coefficient (K) tables
- Friction factors
- Darcy-Weisbach calculations
- Component equivalent lengths
- Crane TP-410 data

#### 8. NASA Earthdata (`nasa-earthdata`)
**Purpose:** Access atmospheric properties and aerospace fluid data
**Files:** SKILL.md, reference.md
**Key Features:**
- Standard atmosphere models
- Altitude-dependent properties
- Free registration required
- API access methods

---

### 🔧 Packages (10 skills)

#### 1. Fluids Package (`fluids-package`)
**Purpose:** Pipe flow, pump sizing, friction factor, compressible flow
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Reynolds number calculations
- Friction factor (Moody, Colebrook)
- Pump affinity laws
- Specific speed
- Head loss calculations
- Verified against textbooks

#### 2. NumPy Numerics (`numpy-numerics`)
**Purpose:** Numerical array operations for velocity fields and pump curves
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Array manipulation
- Velocity field calculations
- Matrix operations for pipe networks
- Interpolation
- Statistical analysis

#### 3. SciPy Optimization (`scipy-optimization`)
**Purpose:** Optimize pump designs and system parameters
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Constrained optimization
- Curve fitting to pump data
- Multi-objective optimization
- Global optimization (differential evolution)
- Design optimization examples

#### 4. SymPy Symbolic (`sympy-symbolic`)
**Purpose:** Derive and solve fluid dynamics equations symbolically
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Euler turbine equation derivation
- Bernoulli equation solving
- Symbolic calculus
- Unit verification

#### 5. Matplotlib Visualization (`matplotlib-visualization`)
**Purpose:** Create pump performance curves and engineering plots
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- H-Q, P-Q, η-Q curves
- Contour plots
- Vector fields
- Publication-quality styling
- Multi-panel figures

#### 6. Pint Units (`pint-units`)
**Purpose:** Handle engineering units with automatic conversion
**Files:** SKILL.md, examples.py
**Key Features:**
- Flow rate conversions (m³/s, gpm, L/min)
- Pressure conversions (Pa, psi, bar)
- Viscosity conversions
- Dimensional analysis
- Array operations with units

#### 7. Thermo Package (`thermo-package`)
**Purpose:** Calculate thermodynamic properties and mixture behavior
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Pure component properties
- Mixture calculations
- Phase equilibria
- Enthalpy calculations for pump staging

#### 8. NetworkX Flow Networks (`networkx-flow-networks`)
**Purpose:** Model hydraulic networks and multi-pump systems
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Piping network graphs
- Flow distribution
- Parallel pump systems
- Network optimization

#### 9. PyVista Visualization (`pyvista-visualization`)
**Purpose:** Create 3D visualizations of velocity fields and CFD results
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- 3D mesh visualization
- Velocity field contours
- Streamlines
- Volume rendering
- Interactive plots

#### 10. CoolProp Package (`coolprop-package`)
**Purpose:** Interface to CoolProp for real fluid properties
*Note: Combined with coolprop-db for efficiency*

---

### 🔌 Integrations (4 skills)

#### 1. OpenFOAM CFD (`openfoam-cfd`)
**Purpose:** Set up turbulent flow simulations in OpenFOAM
**Files:** SKILL.md, reference.md, examples/ (pipe-flow, pump-impeller, scripts)
**Key Features:**
- Complete case setup
- k-ω SST configuration
- Boundary conditions
- Python case generation
- Batch execution
- Complete working examples

#### 2. ANSYS Simulation (`ansys-simulation`)
**Purpose:** Automate ANSYS Fluent CFD via Python and journal files
**Files:** SKILL.md, api-reference.md, examples/ (3 scripts)
**Key Features:**
- PyAnsys API
- Journal file automation
- Pump CFD setup
- Turbulence model configuration
- Batch simulation

#### 3. SolidWorks CAD (`solidworks-cad`)
**Purpose:** Automate parametric pump impeller design
**Files:** SKILL.md, api-reference.md, script-examples.py
**Key Features:**
- VBA/Python API
- Parametric impeller creation
- Geometry export (STEP, IGES)
- Win32com interface

#### 4. COMSOL Multiphysics (`comsol-multiphysics`)
**Purpose:** Set up coupled fluid-structure interaction
**Files:** SKILL.md, api-reference.md, pump-fsi-workflow.md
**Key Features:**
- Java API
- FSI analysis
- Pump vibration under fluid loads
- Thermal-flow coupling

---

### 🛠️ Helpers (6 skills)

#### 1. Unit Converter (`unit-converter`)
**Purpose:** Convert flow rates, pressures, viscosities, and other engineering units
**Files:** SKILL.md, converter.py, QUICK_REFERENCE.md
**Key Features:**
- Pint-based conversions
- Flow, pressure, viscosity units
- Quick reference tables
- Dimensional analysis

#### 2. Pump Selection Helper (`pump-selection-helper`)
**Purpose:** Decision tree for selecting pump type based on requirements
**Files:** SKILL.md, selector.py, reference.md
**Key Features:**
- Centrifugal vs PD decision tree
- Specific speed classification
- Application-based recommendations
- Interactive selection

#### 3. Fluid Property Calculator (`fluid-property-calculator`)
**Purpose:** Quick fluid property calculations using empirical formulas
**Files:** SKILL.md, calc.py, reference.md
**Key Features:**
- Water properties (T-dependent)
- Air properties
- Sutherland viscosity law
- Antoine vapor pressure
- No database required

#### 4. Engineering Context Init (`engineering-context-init`)
**Purpose:** Initialize engineering session with standard constants and imports
**Files:** SKILL.md, init-script.py, constants.md
**Key Features:**
- Standard constants (g=9.81, R, etc.)
- Unit registry setup
- Common imports
- Plotting defaults

#### 5. Error Handler for Fluids (`error-handler-fluids`)
**Purpose:** Handle common numerical errors in fluid calculations
**Files:** SKILL.md, handler.py, best-practices.md
**Key Features:**
- Input validation
- Physical bounds checking
- Numerical stability
- Graceful error handling
- Logging best practices

#### 6. Get Available Resources (`get-available-resources`)
**Purpose:** Enumerate installed packages, databases, and tools
**Files:** SKILL.md, resource-lister.py
**Key Features:**
- Check package installations
- Test database connectivity
- Verify software paths
- Resource reporting

---

### 🧠 Thinking Skills (9 skills)

#### 1. Fluid Dynamics Workflow (`fluid-dynamics`)
**Purpose:** Systematic workflow for fluid dynamics analysis from setup to validation
**Files:** SKILL.md, equations-reference.md, workflow-examples.md
**Key Features:**
- Complete CFD workflow
- Governing equations (Navier-Stokes)
- Boundary condition selection
- Turbulence model selection
- Verification and validation procedures

#### 2. Thermodynamics Workflow (`thermodynamics`)
**Purpose:** Analyze thermodynamic cycles and heat transfer in pump systems
**Files:** SKILL.md, equations-reference.md, examples.md
**Key Features:**
- Cycle analysis (Rankine, Brayton)
- State point analysis
- Efficiency calculations
- Heat transfer in pumps

#### 3. Structural Analysis Workflow (`structural-analysis`)
**Purpose:** FEA workflow for pump casings and impellers
**Files:** SKILL.md, workflow-examples.md, reference.md
**Key Features:**
- FEA methodology
- Load application (pressure, centrifugal)
- Mesh generation
- Stress analysis
- Safety factor evaluation
- ASME Section VIII compliance

#### 4. Centrifugal Pump Design (`pump-design/centrifugal-pumps`)
**Purpose:** Design centrifugal pumps using Euler equations and velocity triangles
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- Complete design workflow
- Specific speed calculation
- Euler turbine equation
- Velocity triangles
- Impeller sizing (Stepanoff)
- Performance prediction
- Verified design examples

#### 5. Cavitation Analysis (`pump-design/cavitation-analysis`)
**Purpose:** Calculate NPSH and assess cavitation risk
**Files:** SKILL.md, examples.py, reference.md
**Key Features:**
- NPSH_available calculation
- NPSH_required correlations
- Temperature effects
- Safety margins
- Inducer design guidance
- Critical for pump safety

#### 6. Pump Performance Curves (`pump-design/performance-curves`)
**Purpose:** Generate and interpret H-Q curves, apply affinity laws
**Files:** SKILL.md, plotter.py
**Key Features:**
- H-Q, η-Q, P-Q curve generation
- Affinity laws (speed, diameter)
- Curve fitting from test data
- System curve intersection
- BEP identification

#### 7. Positive Displacement Pumps (`pump-design/positive-displacement-pumps`)
**Purpose:** Design and analyze gear, piston, and screw pumps
**Files:** SKILL.md, examples.py
**Key Features:**
- PD pump types (gear, piston, screw)
- Volumetric efficiency
- Displacement calculations
- Slip and leakage
- Pulsation analysis
- When to use vs centrifugal

#### 8. Pump System Integration (`pump-design/system-integration`)
**Purpose:** Design complete pump systems with piping and controls
**Files:** SKILL.md, network-model.py, reference.md
**Key Features:**
- System curve development
- Parallel pump configurations
- Series pump arrangements
- Piping network modeling
- Control strategies (VFD)
- Network optimization

#### 9. Pump Efficiency Optimization (`pump-design/efficiency-optimization`)
**Purpose:** Maximize pump efficiency through design and operation
**Files:** SKILL.md, optimizer.py
**Key Features:**
- Efficiency fundamentals
- Loss mechanisms
- Design optimization (impeller geometry)
- Operational optimization (VFD, sequencing)
- Multi-objective optimization
- Energy cost analysis

---

## Documentation Files

### Core Documentation (3 files)

1. **README.md** - Complete repository overview, quick start, examples, FAQs
2. **LICENSE.md** - MIT license with third-party dependency notes
3. **CONTRIBUTING.md** - Comprehensive contribution guidelines with verification requirements

### Docs Directory (6 files)

1. **examples.md** (29KB) - Complete end-to-end workflow examples:
   - Pump performance analysis
   - CFD simulation setup
   - Multi-pump system design
   - Cavitation risk assessment
   - With full working code

2. **engineering-databases.md** (82KB) - Database documentation:
   - 10 engineering databases
   - Access methods and authentication
   - Query examples
   - API usage

3. **engineering-packages.md** (77KB) - Python package guide:
   - 15+ packages documented
   - Installation and key functions
   - Engineering applications
   - Working code examples

4. **engineering-integrations.md** (98KB) - Software integration:
   - OpenFOAM, ANSYS, SolidWorks, COMSOL
   - API scripting
   - Automation examples
   - Setup instructions

5. **engineering-thinking.md** (83KB) - Structured workflows:
   - Fluid dynamics methodology
   - Pump design process
   - Thermodynamic analysis
   - FEA workflow
   - Decision trees

6. **engineering-helpers.md** (65KB) - Utility functions:
   - Unit conversions
   - Quick calculators
   - Context initialization
   - Error handling

---

## Key Features of This Library

### ✅ Verified Engineering Accuracy
- All formulas verified against textbooks and standards
- Test cases from published literature
- References to engineering codes (ASME, ISO, API)
- Validation against NIST data

### 💻 Complete Working Code
- 35+ Python example files
- Ready-to-run calculations
- OpenFOAM case files
- Automation scripts

### 📚 Comprehensive Documentation
- 121+ total files
- Detailed SKILL.md for each skill
- Reference documentation
- Best practices

### 🔌 Production Ready
- Error handling
- Input validation
- Unit consistency checking
- Professional code quality

### 🎯 Focused on Pump Engineering
- 9 dedicated pump skills
- Centrifugal and positive displacement
- Performance, cavitation, system integration
- Design optimization

### 🌊 Complete Fluid Dynamics Coverage
- CFD workflows
- Turbulence modeling
- Pipe flows
- OpenFOAM integration

### 🔧 Professional Tool Integration
- OpenFOAM (open source)
- ANSYS (commercial CFD)
- SolidWorks (CAD)
- COMSOL (multiphysics)

---

## Usage Statistics

### Lines of Documentation
- SKILL.md files: ~8,000 lines
- Reference docs: ~4,000 lines
- Code examples: ~3,000 lines
- Total: 15,000+ lines

### Engineering Domains
- Mechanical Engineering ✓
- Fluid Dynamics ✓
- Aerospace Engineering ✓
- Thermodynamics ✓
- Structural Analysis ✓

### Standards Referenced
- ASME B31.1/B31.3 (Piping)
- API 610 (Centrifugal Pumps)
- ISO 5199, 9906 (Pump Standards)
- HI 9.6.1, 9.6.3 (Hydraulic Institute)
- ASME Section VIII (Pressure Vessels)
- ASME V&V 20 (CFD Verification)

---

## Skill Activation Examples

Claude will automatically activate relevant skills based on user requests:

**User Request:** "Analyze a centrifugal pump for 100 m³/h at 50m head"

**Skills Activated:**
- `coolprop-db` (fluid properties)
- `pump-design/centrifugal-pumps` (design equations)
- `pump-design/cavitation-analysis` (NPSH check)
- `fluids-package` (Reynolds number)
- `matplotlib-visualization` (performance curves)

**User Request:** "Set up OpenFOAM simulation for turbulent pipe flow"

**Skills Activated:**
- `openfoam-cfd` (case setup)
- `thinking/fluid-dynamics` (workflow)
- `turbulence-models-db` (k-ω SST parameters)
- `coolprop-package` (fluid properties)

---

## Installation Quick Start

```bash
# Clone repository
git clone https://github.com/Soljourner/claude-engineering-skills.git
cd claude-engineering-skills

# Install Python packages
pip install numpy scipy sympy matplotlib pint fluids thermo CoolProp pyvista networkx

# Optional: OpenFOAM for CFD
# See docs/engineering-integrations.md for installation

# Skills are now available to Claude!
```

---

## Future Enhancements

Potential additions (not yet implemented):
- Compressible flow skills (high-speed applications)
- Multiphase flow (beyond cavitation)
- Pump control systems
- Fatigue analysis
- Noise and vibration
- Additional turbomachinery (turbines, compressors)

---

## Credits

Built following Anthropic's Agent Skills guidelines with inspiration from the claude-scientific-skills library structure.

**Created:** November 2025
**Version:** 1.0
**Total Development Time:** ~3 hours (using AI-assisted parallel development)

---

*This library represents a comprehensive foundation for mechanical and aerospace engineering workflows with Claude AI.*
