---
name: solidworks-cad
description: "Automate parametric pump impeller design in SolidWorks via API"
category: integrations
domain: mechanical
complexity: advanced
dependencies: []
---

# SolidWorks CAD Integration

This skill provides guidance for automating parametric pump impeller design and other mechanical component modeling using the SolidWorks API.

## Overview of SolidWorks API

SolidWorks provides a comprehensive API that allows users to automate design tasks, create parametric models, manipulate assemblies, and generate drawings programmatically. The API exposes the full functionality of SolidWorks through COM interfaces, enabling external applications to control SolidWorks and access its features.

Key capabilities:
- Parametric part and assembly creation
- Feature manipulation and dimension control
- Sketch creation and editing
- Drawing automation and annotation
- File import/export operations
- Simulation and analysis integration
- Bill of materials (BOM) generation

## Licensing Requirements

**Important**: SolidWorks API access requires:

1. **Valid SolidWorks License**: A full SolidWorks license is required to use the API. The API cannot be used with viewer-only licenses.

2. **License Types**:
   - SolidWorks Standard, Professional, or Premium
   - Network (concurrent) or standalone licenses both support API access
   - SolidWorks PDM integration requires additional PDM licenses

3. **Windows Platform**: SolidWorks and its API are Windows-only. The API utilizes COM/ActiveX technology specific to Windows.

4. **Version Compatibility**: Ensure your API code matches your SolidWorks version. API methods may vary between versions.

## API Access Methods

### 1. VBA Macros

VBA (Visual Basic for Applications) is the built-in scripting environment in SolidWorks.

**Advantages**:
- Native integration with SolidWorks
- No external dependencies
- Easy to debug within SolidWorks
- Direct access to UI elements

**Use Cases**:
- Quick automation scripts
- User interaction through forms
- Prototyping API workflows

**Access**: Tools > Macro > New or Edit in SolidWorks

### 2. Python (via win32com)

Python can control SolidWorks through the `win32com` package (pywin32).

**Advantages**:
- Modern, readable syntax
- Extensive Python ecosystem for data processing
- Integration with scientific libraries (NumPy, pandas)
- Easier to version control and test

**Requirements**:
```bash
pip install pywin32
```

**Use Cases**:
- Batch processing of parts
- Data-driven design automation
- Integration with engineering calculations
- Export/import workflows

### 3. C# .NET

C# provides strong typing and robust development tools for SolidWorks automation.

**Advantages**:
- Best performance
- Strong typing and IntelliSense support
- Comprehensive error handling
- Professional application development

**Requirements**:
- Visual Studio
- SolidWorks Interop Assemblies (installed with SolidWorks)

**Use Cases**:
- Custom SolidWorks add-ins
- Enterprise integration applications
- Complex automation workflows
- Standalone applications controlling SolidWorks

## Common Tasks

### Parametric Part Creation

Create parts with dimensions that can be modified programmatically:
- Define sketches with dimensional constraints
- Create features (extrude, revolve, sweep, loft)
- Add relations between features
- Set global variables and equations

### Impeller Blade Geometry

Specialized tasks for pump impeller design:
- **Blade Profile Generation**: Create complex 3D curves for blade profiles
- **Swept Features**: Use guide curves and profiles to create blade surfaces
- **Circular Patterns**: Array blades around the central axis
- **Hub and Shroud Modeling**: Create the central hub and outer shroud geometries
- **Fillet Operations**: Add fillets for manufacturing and performance optimization
- **Parametric Control**: Link blade angles, thickness, and count to design variables

### Assembly Automation

Automate assembly creation and modification:
- Insert components into assemblies
- Apply mates (coincident, concentric, distance, angle)
- Create patterns of components
- Detect interferences
- Manage assembly configurations

### Drawing Generation

Automatically create and populate engineering drawings:
- Create drawing sheets from templates
- Insert model views (standard, section, detail)
- Add dimensions and annotations
- Generate BOMs and balloons
- Export to PDF or DWG

### File Export

Export SolidWorks models to neutral formats:
- **STEP (ISO 10303)**: Industry-standard for CAD data exchange
- **IGES**: Legacy format, widely supported
- **Parasolid**: High-fidelity geometry transfer
- **STL**: For 3D printing and mesh-based analysis
- **DXF/DWG**: 2D drawing export for manufacturing

## Best Practices

1. **Error Handling**: Always implement robust error handling as API operations can fail
2. **Version Control**: Keep API code separate from CAD files for better version control
3. **Documentation**: Document units, coordinate systems, and design intent
4. **Testing**: Test API scripts on simple models before applying to complex designs
5. **Performance**: Suppress screen updates and rebuild operations for better performance
6. **Cleanup**: Always release COM objects properly to avoid memory leaks

## Learning Resources

- **SolidWorks API Help**: Access from Help > API Help in SolidWorks
- **API SDK**: Includes examples and documentation (installed with SolidWorks)
- **Online Community**: SolidWorks Forums and user groups
- **Sample Code**: Available in SolidWorks installation directory

## Limitations and Considerations

- **Platform**: Windows-only, no macOS or Linux support
- **Performance**: Complex operations may be slow; consider batch processing during off-hours
- **Licensing**: Requires active SolidWorks session (consumes a license)
- **Version Dependencies**: API changes between SolidWorks versions may break scripts
- **COM Complexity**: Understanding COM object model is essential for advanced use

## Integration with Engineering Workflow

The SolidWorks API is particularly valuable for:
- **Design Optimization**: Parametric studies with multiple design variants
- **CFD Preprocessing**: Automated geometry preparation for OpenFOAM or ANSYS
- **Manufacturing**: Automated drawing generation and BOM creation
- **Product Configurators**: Customer-driven product customization
- **Data Management**: Integration with PLM/PDM systems
