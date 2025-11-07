#!/usr/bin/env python3
"""
Resource Availability Checker for Engineering Workflows

This tool enumerates installed packages, databases, software tools, and system
resources to help Claude determine what's available for engineering tasks.

Usage:
    python3 resource-lister.py                    # Standard report
    python3 resource-lister.py --detailed         # Detailed report
    python3 resource-lister.py --json             # JSON output
    python3 resource-lister.py --packages-only    # Only Python packages
    python3 resource-lister.py --databases-only   # Only databases
    python3 resource-lister.py --software-only    # Only software tools
    python3 resource-lister.py --check pkg1 pkg2  # Check specific packages
"""

import sys
import os
import subprocess
import platform
import shutil
import json
from datetime import datetime
from pathlib import Path
import argparse


class ResourceChecker:
    """Check availability of packages, databases, software, and system resources."""

    def __init__(self):
        self.results = {
            'timestamp': datetime.now().isoformat(),
            'system': {},
            'python': {},
            'packages': {},
            'databases': {},
            'software': {},
            'environment': {},
            'resources': {},
            'recommendations': []
        }

        # Engineering packages to check
        self.engineering_packages = {
            'core_scientific': ['numpy', 'scipy', 'matplotlib', 'pandas'],
            'engineering': ['fluids', 'thermo', 'CoolProp', 'pint', 'sympy'],
            'optimization': ['pyomo', 'gekko', 'cvxpy', 'pulp'],
            'network': ['networkx'],
            'data_analysis': ['sklearn', 'scikit-learn', 'statsmodels', 'seaborn'],
            'database': ['psycopg2', 'pymongo', 'redis', 'sqlalchemy', 'pymysql'],
            'cfd': ['pyFoam', 'openfoam', 'fluiddyn']
        }

        # Database commands to check
        self.databases = {
            'PostgreSQL': 'psql',
            'MySQL': 'mysql',
            'SQLite': 'sqlite3',
            'MongoDB': 'mongod',
            'Redis': 'redis-server'
        }

        # Engineering software to check
        self.software = {
            'OpenFOAM': ['openfoam', 'simpleFoam', 'icoFoam'],
            'ANSYS Fluent': ['fluent'],
            'ANSYS Workbench': ['workbench'],
            'COMSOL': ['comsol'],
            'MATLAB': ['matlab'],
            'ParaView': ['paraview'],
            'Gmsh': ['gmsh'],
            'SolidWorks': ['SLDWORKS.exe'],  # Windows only
            'FreeCAD': ['freecad', 'FreeCAD'],
            'Salome': ['salome']
        }

    def check_system_info(self):
        """Gather system information."""
        try:
            self.results['system'] = {
                'os': platform.system(),
                'os_version': platform.version(),
                'distribution': self._get_linux_distribution(),
                'architecture': platform.machine(),
                'processor': platform.processor(),
                'hostname': platform.node()
            }
        except Exception as e:
            self.results['system']['error'] = str(e)

    def _get_linux_distribution(self):
        """Get Linux distribution info."""
        try:
            if platform.system() == 'Linux':
                # Try reading /etc/os-release
                if os.path.exists('/etc/os-release'):
                    with open('/etc/os-release', 'r') as f:
                        lines = f.readlines()
                        info = {}
                        for line in lines:
                            if '=' in line:
                                key, value = line.strip().split('=', 1)
                                info[key] = value.strip('"')
                        return f"{info.get('NAME', 'Unknown')} {info.get('VERSION', '')}"
                else:
                    return platform.platform()
            return platform.platform()
        except:
            return platform.platform()

    def check_python_info(self):
        """Gather Python interpreter information."""
        self.results['python'] = {
            'version': platform.python_version(),
            'implementation': platform.python_implementation(),
            'executable': sys.executable,
            'path': sys.path[:3]  # First 3 paths
        }

    def check_package_installed(self, package_name):
        """Check if a Python package is installed and get its version."""
        try:
            # Try importing the package
            if package_name == 'sklearn':
                # Special case for scikit-learn
                import sklearn
                return True, sklearn.__version__
            else:
                module = __import__(package_name)
                version = getattr(module, '__version__', 'unknown')
                return True, version
        except ImportError:
            # Try pip show as fallback
            try:
                result = subprocess.run(
                    ['pip3', 'show', package_name],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0:
                    for line in result.stdout.split('\n'):
                        if line.startswith('Version:'):
                            version = line.split(':', 1)[1].strip()
                            return True, version
            except:
                pass
            return False, None

    def check_all_packages(self):
        """Check all engineering packages."""
        for category, packages in self.engineering_packages.items():
            self.results['packages'][category] = {}
            for package in packages:
                installed, version = self.check_package_installed(package)
                self.results['packages'][category][package] = {
                    'installed': installed,
                    'version': version if installed else None
                }

    def get_all_installed_packages(self):
        """Get complete list of installed packages via pip."""
        try:
            result = subprocess.run(
                ['pip3', 'list', '--format=json'],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                packages = json.loads(result.stdout)
                self.results['packages']['all_installed'] = {
                    pkg['name']: pkg['version'] for pkg in packages
                }
                self.results['packages']['total_count'] = len(packages)
        except Exception as e:
            self.results['packages']['all_installed_error'] = str(e)

    def check_database(self, name, command):
        """Check if a database is installed and accessible."""
        # Check if command exists
        cmd_path = shutil.which(command)

        if not cmd_path:
            return {
                'installed': False,
                'path': None,
                'version': None,
                'status': 'Not found'
            }

        # Get version
        version = self._get_database_version(command)

        # Test connection
        connection_status = self._test_database_connection(name, command)

        return {
            'installed': True,
            'path': cmd_path,
            'version': version,
            'status': connection_status
        }

    def _get_database_version(self, command):
        """Get database version."""
        version_flags = {
            'psql': '--version',
            'mysql': '--version',
            'sqlite3': '--version',
            'mongod': '--version',
            'redis-server': '--version'
        }

        flag = version_flags.get(command, '--version')
        try:
            result = subprocess.run(
                [command, flag],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                # Extract version number from output
                output = result.stdout.split('\n')[0]
                return output
        except:
            pass
        return 'unknown'

    def _test_database_connection(self, name, command):
        """Test database connection."""
        try:
            if name == 'PostgreSQL':
                # Try to connect to default postgres database
                result = subprocess.run(
                    ['psql', '-U', 'postgres', '-c', 'SELECT version();'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                return 'Connected' if result.returncode == 0 else 'Installed, connection failed'

            elif name == 'Redis':
                # Check if Redis is running
                result = subprocess.run(
                    ['redis-cli', 'ping'],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                return 'Running' if 'PONG' in result.stdout else 'Installed, not running'

            elif name == 'SQLite':
                # SQLite is always available if command exists
                return 'Available (file-based)'

            else:
                return 'Installed'
        except:
            return 'Installed, status unknown'

    def check_all_databases(self):
        """Check all databases."""
        for name, command in self.databases.items():
            self.results['databases'][name] = self.check_database(name, command)

    def check_software(self, name, commands):
        """Check if software is installed."""
        for cmd in commands:
            cmd_path = shutil.which(cmd)
            if cmd_path:
                version = self._get_software_version(cmd)
                return {
                    'installed': True,
                    'path': cmd_path,
                    'version': version,
                    'command': cmd
                }
        return {
            'installed': False,
            'path': None,
            'version': None
        }

    def _get_software_version(self, command):
        """Get software version."""
        version_flags = ['--version', '-version', '-v', 'version']

        for flag in version_flags:
            try:
                result = subprocess.run(
                    [command, flag],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0 and result.stdout.strip():
                    # Return first line of output
                    return result.stdout.split('\n')[0][:100]
            except:
                continue

        return 'unknown'

    def check_all_software(self):
        """Check all engineering software."""
        for name, commands in self.software.items():
            self.results['software'][name] = self.check_software(name, commands)

    def check_environment_variables(self):
        """Check important environment variables."""
        important_vars = [
            'PATH', 'LD_LIBRARY_PATH', 'PYTHONPATH',
            'OPENFOAM_DIR', 'FOAM_RUN',
            'ANSYS_LICENSE_FILE', 'LM_LICENSE_FILE',
            'MATLAB_ROOT', 'DATABASE_URL'
        ]

        for var in important_vars:
            value = os.environ.get(var)
            if value:
                # Truncate long paths
                if len(value) > 200:
                    value = value[:200] + '...'
                self.results['environment'][var] = value

    def check_system_resources(self):
        """Check system resources (CPU, RAM, disk)."""
        try:
            # CPU info
            try:
                cpu_count = os.cpu_count()
                self.results['resources']['cpu_cores'] = cpu_count
            except:
                pass

            # Memory info (Linux)
            if platform.system() == 'Linux':
                try:
                    with open('/proc/meminfo', 'r') as f:
                        meminfo = f.read()
                        for line in meminfo.split('\n'):
                            if line.startswith('MemTotal:'):
                                total_kb = int(line.split()[1])
                                self.results['resources']['ram_total_gb'] = round(total_kb / 1024 / 1024, 2)
                            elif line.startswith('MemAvailable:'):
                                avail_kb = int(line.split()[1])
                                self.results['resources']['ram_available_gb'] = round(avail_kb / 1024 / 1024, 2)
                except:
                    pass

            # Disk space
            try:
                stat = shutil.disk_usage('/')
                self.results['resources']['disk_total_gb'] = round(stat.total / (1024**3), 2)
                self.results['resources']['disk_available_gb'] = round(stat.free / (1024**3), 2)
            except:
                pass

        except Exception as e:
            self.results['resources']['error'] = str(e)

    def generate_recommendations(self):
        """Generate recommendations based on missing packages."""
        recommendations = []

        # Check core scientific computing
        core_missing = []
        for pkg in self.engineering_packages['core_scientific']:
            if not self.results['packages']['core_scientific'][pkg]['installed']:
                core_missing.append(pkg)

        if core_missing:
            recommendations.append({
                'category': 'Core Scientific Computing',
                'packages': core_missing,
                'command': f"pip install {' '.join(core_missing)}",
                'priority': 'high'
            })

        # Check engineering packages
        eng_missing = []
        for pkg in self.engineering_packages['engineering']:
            if not self.results['packages']['engineering'][pkg]['installed']:
                eng_missing.append(pkg)

        if eng_missing:
            recommendations.append({
                'category': 'Engineering Libraries',
                'packages': eng_missing,
                'command': f"pip install {' '.join(eng_missing)}",
                'priority': 'medium'
            })

        # Check optimization packages
        opt_missing = []
        for pkg in self.engineering_packages['optimization']:
            if not self.results['packages']['optimization'][pkg]['installed']:
                opt_missing.append(pkg)

        if opt_missing:
            recommendations.append({
                'category': 'Optimization Tools',
                'packages': opt_missing,
                'command': f"pip install {' '.join(opt_missing)}",
                'priority': 'low'
            })

        self.results['recommendations'] = recommendations

    def run_all_checks(self):
        """Run all checks."""
        print("Checking system information...")
        self.check_system_info()

        print("Checking Python environment...")
        self.check_python_info()

        print("Checking Python packages...")
        self.check_all_packages()
        self.get_all_installed_packages()

        print("Checking databases...")
        self.check_all_databases()

        print("Checking software tools...")
        self.check_all_software()

        print("Checking environment variables...")
        self.check_environment_variables()

        print("Checking system resources...")
        self.check_system_resources()

        print("Generating recommendations...")
        self.generate_recommendations()

    def print_report(self, detailed=False):
        """Print human-readable report."""
        print("\n" + "="*60)
        print("ENGINEERING ENVIRONMENT RESOURCE REPORT")
        print("="*60)
        print(f"Generated: {self.results['timestamp']}")
        print(f"System: {self.results['system'].get('distribution', 'Unknown')}")
        print(f"Architecture: {self.results['system'].get('architecture', 'Unknown')}")
        print(f"Python: {self.results['python']['version']} ({self.results['python']['executable']})")
        print(f"Working Directory: {os.getcwd()}")

        # Python Packages (only if checked)
        if any(cat in self.results['packages'] for cat in self.engineering_packages.keys()):
            print("\n" + "━"*60)
            total = self.results['packages'].get('total_count', 0)
            print(f"PYTHON PACKAGES ({total} installed)")
            print("━"*60)

            for category, packages in self.engineering_packages.items():
                if category not in self.results['packages']:
                    continue
                category_name = category.replace('_', ' ').title()
                print(f"\n{category_name}:")

                for package in packages:
                    info = self.results['packages'][category][package]
                    if info['installed']:
                        print(f"  [✓] {package:20s} {info['version']}")
                    else:
                        print(f"  [✗] {package:20s} Not installed")

            if detailed:
                print(f"\nAll Installed Packages ({total}):")
                all_pkgs = self.results['packages'].get('all_installed', {})
                for name, version in sorted(all_pkgs.items())[:20]:  # Show first 20
                    print(f"  {name:30s} {version}")
                if len(all_pkgs) > 20:
                    print(f"  ... and {len(all_pkgs) - 20} more")

        # Databases (only if checked)
        if self.results['databases']:
            print("\n" + "━"*60)
            print("DATABASES")
            print("━"*60 + "\n")

            for name, info in self.results['databases'].items():
                if info['installed']:
                    version = info['version'].split()[0] if info['version'] else 'unknown'
                    print(f"[✓] {name:15s} {info['status']:20s} {info['path']}")
                    if detailed and info['version']:
                        print(f"    Version: {info['version']}")
                else:
                    print(f"[✗] {name:15s} Not installed")

        # Software (only if checked)
        if self.results['software']:
            print("\n" + "━"*60)
            print("SOFTWARE TOOLS")
            print("━"*60 + "\n")

            print("Simulation & CFD:")
            for name in ['OpenFOAM', 'ANSYS Fluent', 'ANSYS Workbench', 'COMSOL']:
                info = self.results['software'].get(name, {})
                if info.get('installed'):
                    print(f"  [✓] {name:20s} {info['path']}")
                    if detailed and info.get('version'):
                        print(f"      Version: {info['version']}")
                else:
                    print(f"  [✗] {name:20s} Not found")

            print("\nCAD & Visualization:")
            for name in ['ParaView', 'Gmsh', 'FreeCAD', 'SolidWorks', 'Salome']:
                info = self.results['software'].get(name, {})
                if info.get('installed'):
                    print(f"  [✓] {name:20s} {info['path']}")
                    if detailed and info.get('version'):
                        print(f"      Version: {info['version']}")
                else:
                    print(f"  [✗] {name:20s} Not found")

            print("\nComputing:")
            for name in ['MATLAB']:
                info = self.results['software'].get(name, {})
                if info.get('installed'):
                    print(f"  [✓] {name:20s} {info['path']}")
                else:
                    print(f"  [✗] {name:20s} Not found")

        # System Resources
        if self.results['resources']:
            print("\n" + "━"*60)
            print("SYSTEM RESOURCES")
            print("━"*60 + "\n")

            res = self.results['resources']
            if 'cpu_cores' in res:
                print(f"CPU Cores: {res['cpu_cores']}")
            if 'ram_total_gb' in res:
                print(f"RAM: {res['ram_total_gb']} GB total, {res.get('ram_available_gb', '?')} GB available")
            if 'disk_total_gb' in res:
                print(f"Disk: {res['disk_available_gb']} GB available (total: {res['disk_total_gb']} GB)")

        # Environment Variables
        if detailed and self.results['environment']:
            print("\n" + "━"*60)
            print("ENVIRONMENT VARIABLES")
            print("━"*60 + "\n")

            for var, value in sorted(self.results['environment'].items()):
                print(f"{var}:")
                print(f"  {value}")

        # Recommendations
        if self.results['recommendations']:
            print("\n" + "━"*60)
            print("RECOMMENDATIONS")
            print("━"*60 + "\n")

            for rec in self.results['recommendations']:
                print(f"{rec['category']} (Priority: {rec['priority']}):")
                print(f"  Missing: {', '.join(rec['packages'])}")
                print(f"  Install: {rec['command']}\n")
        else:
            print("\n" + "━"*60)
            print("All essential engineering packages are installed!")
            print("━"*60)

        print()

    def export_json(self):
        """Export results as JSON."""
        return json.dumps(self.results, indent=2)

    def check_specific_packages(self, package_names):
        """Check specific packages only."""
        print("\n" + "="*60)
        print("PACKAGE AVAILABILITY CHECK")
        print("="*60 + "\n")

        for package in package_names:
            installed, version = self.check_package_installed(package)
            if installed:
                print(f"[✓] {package:20s} {version}")
            else:
                print(f"[✗] {package:20s} Not installed")
        print()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Check available resources for engineering workflows'
    )
    parser.add_argument(
        '--detailed', '-d',
        action='store_true',
        help='Show detailed report with all packages and versions'
    )
    parser.add_argument(
        '--json', '-j',
        action='store_true',
        help='Output in JSON format'
    )
    parser.add_argument(
        '--packages-only', '-p',
        action='store_true',
        help='Check only Python packages'
    )
    parser.add_argument(
        '--databases-only', '-db',
        action='store_true',
        help='Check only databases'
    )
    parser.add_argument(
        '--software-only', '-s',
        action='store_true',
        help='Check only software tools'
    )
    parser.add_argument(
        '--check', '-c',
        nargs='+',
        metavar='PACKAGE',
        help='Check specific packages'
    )
    parser.add_argument(
        '--output', '-o',
        metavar='FILE',
        help='Save report to file'
    )

    args = parser.parse_args()

    checker = ResourceChecker()

    # Handle specific package check
    if args.check:
        checker.check_specific_packages(args.check)
        return

    # Run checks based on flags
    if args.packages_only:
        checker.check_system_info()
        checker.check_python_info()
        checker.check_all_packages()
        checker.get_all_installed_packages()
    elif args.databases_only:
        checker.check_system_info()
        checker.check_python_info()
        checker.check_all_databases()
    elif args.software_only:
        checker.check_system_info()
        checker.check_python_info()
        checker.check_all_software()
    else:
        # Run all checks
        checker.run_all_checks()

    # Output results
    if args.json:
        output = checker.export_json()
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output)
            print(f"Report saved to {args.output}")
        else:
            print(output)
    else:
        if args.output:
            # Redirect stdout to file
            import io
            old_stdout = sys.stdout
            sys.stdout = buffer = io.StringIO()

            checker.print_report(detailed=args.detailed)

            output = buffer.getvalue()
            sys.stdout = old_stdout

            with open(args.output, 'w') as f:
                f.write(output)
            print(f"Report saved to {args.output}")
        else:
            checker.print_report(detailed=args.detailed)


if __name__ == '__main__':
    main()
