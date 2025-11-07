---
name: get-available-resources
description: "Enumerate installed packages, databases, and tools for engineering workflows"
category: helpers
domain: general
complexity: basic
dependencies: []
---

# Resource Availability Helper

## Purpose

Help Claude quickly determine what resources, tools, packages, and databases are available in the current environment. This skill eliminates guesswork and enables Claude to:
- Use only installed packages in code suggestions
- Recommend appropriate tools for engineering tasks
- Identify missing dependencies before starting work
- Verify database connectivity for data-driven workflows
- Discover available simulation tools (OpenFOAM, ANSYS, etc.)

This diagnostic tool should be run at the beginning of engineering sessions to establish environmental context.

## What It Checks

### Python Packages
- **Core Scientific Computing**: numpy, scipy, pandas, matplotlib
- **Engineering & Physics**: fluids, thermo, CoolProp, pint, sympy
- **Optimization**: scipy.optimize, pyomo, gekko, cvxpy
- **CFD & Simulation**: OpenFOAM bindings, PyFluent, pyMesh
- **Data Analysis**: sklearn, statsmodels, seaborn
- **Network Analysis**: networkx
- **Database Connectors**: psycopg2, pymongo, redis, sqlalchemy
- **All Installed Packages**: Complete pip list with versions

### Database Connections
- **PostgreSQL**: Connection test and version check
- **MySQL/MariaDB**: Availability and connectivity
- **SQLite**: Built-in database availability
- **MongoDB**: Server status and connection
- **Redis**: In-memory database availability
- **Database connection strings**: Environment variable detection

### Software Installations
- **OpenFOAM**: Version detection and installation path
- **ANSYS**: Workbench, Fluent, CFX availability
- **COMSOL Multiphysics**: Installation check
- **MATLAB**: Version and toolbox detection
- **SolidWorks**: CAD software availability
- **ParaView**: Visualization tool detection
- **GiD**: Pre/post-processor availability
- **Gmsh**: Mesh generator detection

### Environment Variables
- **PATH**: Executable search paths
- **LD_LIBRARY_PATH**: Dynamic library paths
- **PYTHONPATH**: Python module search paths
- **Database URLs**: Connection string variables
- **License Servers**: ANSYS_LICENSE_FILE, LM_LICENSE_FILE
- **Tool-Specific**: OPENFOAM_DIR, MATLAB_ROOT, etc.

### System Information
- **Operating System**: Linux distribution, version
- **Python Version**: Interpreter version and location
- **Architecture**: 64-bit vs 32-bit
- **Available Disk Space**: Critical for large simulations
- **Memory**: Total and available RAM
- **CPU Info**: Cores and model for parallel processing

### File Paths
- **Project Directories**: Common engineering workspace locations
- **Data Directories**: Standard data storage paths
- **Temporary Storage**: Scratch space for simulations
- **License Files**: Software license locations

## Resource Reporting Format

The tool generates a structured report with:

### Summary Section
```
=== ENVIRONMENT RESOURCE REPORT ===
Timestamp: 2025-11-07 21:30:45
Python: 3.11.14 (/usr/bin/python3)
OS: Ubuntu 22.04 LTS (Linux 5.15.0)
Architecture: x86_64
```

### Package Availability
```
[✓] numpy 1.24.3
[✓] scipy 1.10.1
[✓] pandas 2.0.2
[✗] CoolProp (not installed)
[✓] matplotlib 3.7.1
```

### Database Status
```
PostgreSQL: ✓ Connected (version 14.8)
Redis: ✓ Running (version 7.0.11)
MySQL: ✗ Not installed
MongoDB: ✗ Not running
```

### Software Tools
```
OpenFOAM: ✓ v10 (/opt/openfoam10)
ANSYS: ✗ Not found
MATLAB: ✗ Not found
ParaView: ✓ v5.11.0 (/usr/bin/paraview)
```

### Recommendations
```
MISSING PACKAGES FOR FULL ENGINEERING WORKFLOW:
- pip install CoolProp (thermodynamic properties)
- pip install fluids (fluid mechanics calculations)
- pip install thermo (chemical engineering thermodynamics)

OPTIMIZATION:
- Consider installing pyomo for optimization problems
```

## Usage

### Quick Check
```bash
python3 resource-lister.py
```

### Detailed Report
```bash
python3 resource-lister.py --detailed
```

### Check Specific Categories
```bash
# Check only Python packages
python3 resource-lister.py --packages-only

# Check only databases
python3 resource-lister.py --databases-only

# Check only software tools
python3 resource-lister.py --software-only

# Check specific package
python3 resource-lister.py --check numpy scipy pandas
```

### JSON Output
```bash
# Machine-readable format for automation
python3 resource-lister.py --json > resources.json
```

### Export Report
```bash
# Save detailed report to file
python3 resource-lister.py --detailed --output report.txt
```

## Example Output

```
=== ENGINEERING ENVIRONMENT RESOURCE REPORT ===
Generated: 2025-11-07 21:30:45
System: Linux Ubuntu 22.04 LTS (x86_64)
Python: 3.11.14 (/usr/bin/python3)
Working Directory: /home/user/projects

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PYTHON PACKAGES (38 installed)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Core Scientific Computing:
  [✓] numpy           1.24.3
  [✓] scipy           1.10.1
  [✓] matplotlib      3.7.1
  [✗] pandas          Not installed

Engineering Libraries:
  [✗] fluids          Not installed
  [✗] thermo          Not installed
  [✗] CoolProp        Not installed
  [✗] pint            Not installed
  [✓] sympy           1.12

Optimization:
  [✓] scipy.optimize  (included in scipy)
  [✗] pyomo           Not installed
  [✗] gekko           Not installed

Network Analysis:
  [✗] networkx        Not installed

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DATABASES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

[✓] PostgreSQL  14.8    /usr/bin/psql
    Connection: SUCCESS (test database)

[✓] Redis       7.0.11  /usr/bin/redis-server
    Status: Running on port 6379

[✗] MySQL       Not installed
[✗] MongoDB     Not installed
[✓] SQLite3     3.37.2  (built-in)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SOFTWARE TOOLS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Simulation & CFD:
  [✗] OpenFOAM        Not found
  [✗] ANSYS Fluent    Not found
  [✗] COMSOL          Not found

CAD & Visualization:
  [✓] ParaView        5.11.0  /usr/bin/paraview
  [✗] SolidWorks      Not found (Windows only)
  [✓] Gmsh            4.10.5  /usr/bin/gmsh

Computing:
  [✗] MATLAB          Not found

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SYSTEM RESOURCES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CPU: 8 cores (Intel Xeon E5-2680 v4)
RAM: 32.0 GB total, 28.4 GB available
Disk: 450 GB available in /home
Tmp Space: 50 GB available in /tmp

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RECOMMENDATIONS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Essential Engineering Packages:
  pip install pandas fluids thermo pint

Recommended for Fluid Mechanics:
  pip install CoolProp networkx

For Optimization Workflows:
  pip install pyomo

For Enhanced Data Analysis:
  pip install seaborn scikit-learn

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

## Integration with Claude Workflows

### Typical Workflow
1. **Start Session**: Run resource-lister.py to establish baseline
2. **Review Report**: Claude reviews available packages and tools
3. **Plan Approach**: Select methods based on available resources
4. **Code Generation**: Use only verified available packages
5. **Error Handling**: If missing dependencies detected, suggest installation

### Use Cases

#### Case 1: Pump Design Task
```
User: "Design a centrifugal pump for water transport"

Claude runs resource-lister.py and sees:
- numpy, scipy, matplotlib: Available
- fluids, thermo: Not available
- CoolProp: Not available

Claude Response:
"I can design the pump using numpy and scipy for calculations.
However, I recommend installing 'fluids' package for more
accurate hydraulic calculations:
  pip install fluids

For now, I'll use scipy for fluid properties with manual
property definitions."
```

#### Case 2: CFD Simulation
```
User: "Set up an OpenFOAM simulation for pipe flow"

Claude runs resource-lister.py and sees:
- OpenFOAM: Not installed
- ParaView: Available

Claude Response:
"OpenFOAM is not currently installed. I can:
1. Provide installation instructions for OpenFOAM
2. Create case files that will work once installed
3. Use ParaView for post-processing once simulation runs

Would you like me to proceed with case setup or help
with OpenFOAM installation?"
```

#### Case 3: Database-Backed Analysis
```
User: "Analyze pump performance data from database"

Claude runs resource-lister.py and sees:
- PostgreSQL: Connected and running
- pandas: Not installed
- numpy, scipy: Available

Claude Response:
"I can query your PostgreSQL database and perform analysis.
For better data handling, I recommend:
  pip install pandas

For now, I'll use numpy arrays to process the query results."
```

## Files Included

- **SKILL.md**: This documentation
- **resource-lister.py**: Main diagnostic script

## Notes

- The tool requires no external dependencies (uses only Python standard library)
- Database connection tests are non-destructive and read-only
- Software detection is path-based and environment variable-based
- Reports can be cached for session duration to avoid repeated checks
- Safe to run repeatedly without side effects

## Best Practices

1. **Run at Session Start**: Establish baseline before beginning work
2. **Update After Installations**: Re-run after installing new packages
3. **Share Reports**: Include in issue reports for debugging
4. **Cache Results**: Store for reference throughout session
5. **Version Tracking**: Monitor version changes between sessions

## Extensibility

Add custom checks by modifying `resource-lister.py`:
- Domain-specific packages (aerospace, chemical, structural)
- Custom software installations
- Cloud service availability
- License server status
- Network-mounted resources
