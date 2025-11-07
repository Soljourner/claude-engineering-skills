# SolidWorks API Reference

Quick reference for the most commonly used SolidWorks API objects and methods for parametric pump impeller design and mechanical engineering automation.

## Core API Objects

### SldWorks Application

The main application object for accessing SolidWorks functionality.

```python
# Python
import win32com.client
sw_app = win32com.client.Dispatch("SldWorks.Application")
sw_app.Visible = True
```

```vba
' VBA
Dim swApp As SldWorks.SldWorks
Set swApp = Application.SldWorks
```

**Key Methods:**
- `NewDocument(templatePath, paperSize, width, height)` - Create new document
- `OpenDoc(filePath, docType)` - Open existing document
- `ActiveDoc` - Get currently active document
- `CloseAllDocuments(includeUnsaved)` - Close all open documents
- `GetUserPreferenceIntegerValue(preference)` - Get system preferences

### ModelDoc2

Represents a SolidWorks document (part, assembly, or drawing).

**Key Methods:**
- `FeatureManager` - Access feature creation and management
- `SelectionManager` - Access current selections
- `SketchManager` - Access sketch creation tools
- `Extension` - Access extended functionality
- `EditRebuild3()` - Rebuild the model
- `ForceRebuild3(topOnly)` - Force complete rebuild
- `SaveAs3(filePath, saveAsVersion, options)` - Save document
- `GetPathName()` - Get full file path
- `FeatureByName(name)` - Get feature by name
- `Parameter(name)` - Get dimension or parameter by name

### FeatureManager

Manages features in parts and assemblies.

**Key Methods:**
- `InsertSketch()` - Insert new sketch
- `FeatureExtrusion2()` - Create extruded boss/base or cut
- `FeatureRevolve2()` - Create revolved feature
- `FeatureSweep2()` - Create swept feature
- `FeatureLoft2()` - Create lofted feature
- `FeatureFillet3()` - Create fillet
- `FeatureChamfer()` - Create chamfer
- `FeatureCircularPattern4()` - Create circular pattern
- `FeatureLinearPattern2()` - Create linear pattern
- `InsertMirrorFeature2()` - Create mirror feature

### SketchManager

Controls sketch creation and editing.

**Key Methods:**
- `InsertSketch(markVisible)` - Start new sketch
- `CreateLine(x1, y1, z1, x2, y2, z2)` - Create line
- `CreateCircle(x, y, z, radius)` - Create circle
- `CreateCircleByRadius(x, y, z, radius)` - Create circle with radius
- `CreateArc(xc, yc, zc, x1, y1, z1, x2, y2, z2, direction)` - Create arc
- `CreateSpline2(points, addDim, closed)` - Create spline
- `CreateCornerRectangle(x1, y1, z1, x2, y2, z2)` - Create rectangle
- `AddTangentArc(x, y, z)` - Add tangent arc to line
- `CreateSketchPoint(x, y, z)` - Create sketch point
- `SketchUseEdge3()` - Convert edge to sketch entity

### SelectionMgr

Manages object selection in the model.

**Key Methods:**
- `GetSelectedObjectCount()` - Number of selected objects
- `GetSelectedObject6(index, mark)` - Get selected object
- `GetSelectedObjectType3(index, mark)` - Get type of selected object
- `DeSelect2(index, mark)` - Deselect object
- `EnableContourSelection(enable)` - Enable/disable contour selection

### ModelDocExtension

Extended functionality for ModelDoc2.

**Key Methods:**
- `SelectByID2(name, type, x, y, z, append, mark, callout, selectOption)` - Select entity by name
- `SaveAs(filePath, version, options, errors, warnings)` - Save with options
- `SetUserPreferenceInteger(preference, value)` - Set preferences
- `GetMassProperties2(accuracy, relateTo)` - Calculate mass properties
- `InsertComponent(filePath, x, y, z)` - Insert component in assembly

## Parametric Impeller Design Methods

### Creating Impeller Geometry

#### 1. Hub Creation (Central Boss)

```python
# Select plane
model.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, None, 0)

# Create sketch
sketch_mgr = model.SketchManager
sketch_mgr.InsertSketch(True)

# Draw hub profile
hub_radius = 0.025  # 25mm in meters
sketch_mgr.CreateCircleByRadius(0, 0, 0, hub_radius)

# Extrude
height = 0.040  # 40mm
feature_mgr = model.FeatureManager
feature_mgr.FeatureExtrusion2(
    True,      # Extrude both sides = False
    False,     # Flip side to extrude
    False,     # Direction type
    0,         # Direction type enum
    0,         # Direction type enum (other direction)
    height,    # Depth in direction 1
    0,         # Depth in direction 2
    False,     # Merge result
    False,     # Draft outward
    False,     # Draft type
    False,     # Auto-fillet
    0,         # Auto-fillet radius
    0,         # Angle
    False,     # Flip direction
    False,     # Translate surface
    False,     # Use offset from surface
    False,     # Offset reverse direction
    True,      # Boss or cut (True = boss)
    True,      # Join/cut (True = join)
    True,      # Optimization = None
    0,         # Start condition
    0          # End condition
)
```

#### 2. Blade Profile Using Splines

```python
# Create 3D sketch for blade curve
sketch_mgr.Create3DSketch()

# Define blade curve points (parametric blade profile)
points = []
num_points = 10

for i in range(num_points):
    t = i / (num_points - 1)
    radius = hub_radius + t * (shroud_radius - hub_radius)
    angle = inlet_angle + t * (outlet_angle - inlet_angle)
    theta = math.radians(angle)

    x = radius * math.cos(theta)
    y = radius * math.sin(theta)
    z = t * blade_height

    points.extend([x, y, z])

# Create spline through points
sketch_mgr.CreateSpline2(points, True, False)
```

#### 3. Swept Blade Surface

```python
# Select path (blade curve)
model.Extension.SelectByID2("Spline1", "SKETCHSEGMENT", 0, 0, 0, True, 0, None, 0)

# Select profile (blade cross-section)
model.Extension.SelectByID2("Sketch2", "SKETCH", 0, 0, 0, True, 1, None, 0)

# Create swept surface
feature_mgr.InsertProtrusionSwept2(
    False,  # Direction reversed
    False,  # Merge smooth faces
    0,      # Twist type
    0,      # Twist angle
    False,  # Maintain tangency
    0,      # Path alignment
    0,      # Thin feature type
    0,      # Thin feature thickness
    False   # Keep body
)
```

#### 4. Circular Pattern for Multiple Blades

```python
num_blades = 6

# Select feature to pattern
model.Extension.SelectByID2("Boss-Extrude1", "BODYFEATURE", 0, 0, 0, False, 0, None, 0)

# Select axis (Z-axis or central edge)
model.Extension.SelectByID2("Axis1", "AXIS", 0, 0, 0, True, 0, None, 0)

# Create circular pattern
feature_mgr.FeatureCircularPattern4(
    num_blades,              # Number of instances
    2 * math.pi / num_blades,  # Spacing angle (radians)
    False,                   # Equal spacing
    "",                      # Geometry pattern (unused)
    False,                   # Use feature seed
    False,                   # Vary sketch
    True,                    # Geometry pattern
    False,                   # Use instances (unused)
    False                    # Vary instances (unused)
)
```

### Modifying Parameters

#### Access and Modify Dimensions

```python
# Get feature
feature = model.FeatureByName("Boss-Extrude1")

# Get dimension by name
dimension = feature.Parameter("D1@Sketch1")

# Modify value (in meters)
new_diameter = 0.100  # 100mm
dimension.SystemValue = new_diameter

# Rebuild model to apply changes
model.EditRebuild3()
```

#### Using Global Variables

```python
# Add equation (global variable)
equation_mgr = model.GetEquationMgr()

# Add equation: "hub_diameter" = 50mm
equation_mgr.Add(0, '"hub_diameter" = 50', False)

# Link dimension to global variable
dimension.Expression = '"hub_diameter"'
```

### Export Operations

#### Export to STEP (for CFD)

```python
# Set STEP export options
step_version = 214  # AP214 or AP203

# Save as STEP
errors = 0
warnings = 0
output_path = "C:\\exports\\impeller.step"

model.Extension.SaveAs(
    output_path,
    0,          # Save as version (0 = current)
    0,          # Options
    None,       # Configuration data
    errors,
    warnings
)
```

#### Export to STL (for Meshing)

```python
# Configure STL export
stl_data = model.Extension.GetSaveAsStlOptions()

stl_data.OutputFileType = 0  # Binary (smaller file size)
stl_data.QualityType = 1     # Fine quality
stl_data.Deviation = 0.01    # Deviation in mm
stl_data.AngleTolerance = 5  # degrees

# Export
model.SaveAsSilent("C:\\exports\\impeller.stl", True)
```

## Dimension Units

**IMPORTANT**: SolidWorks API uses SI units internally:

| Quantity | SI Unit | Common Unit | Conversion |
|----------|---------|-------------|------------|
| Length | meters | millimeters | mm / 1000 |
| Angle | radians | degrees | deg × π/180 |
| Mass | kilograms | grams | g / 1000 |
| Time | seconds | seconds | 1:1 |
| Force | Newtons | Newtons | 1:1 |

**Example:**
```python
# Wrong: dimension.SystemValue = 50  # This would be 50 meters!
# Correct: dimension.SystemValue = 0.050  # 50mm = 0.050 meters
```

## Selection Types

Common selection type strings for `SelectByID2`:

- `"FACE"` - Face of a body
- `"EDGE"` - Edge
- `"VERTEX"` - Vertex
- `"PLANE"` - Reference plane
- `"AXIS"` - Reference axis
- `"SKETCH"` - Entire sketch
- `"SKETCHSEGMENT"` - Sketch entity (line, arc, etc.)
- `"BODYFEATURE"` - Solid body feature
- `"COMPONENT"` - Assembly component
- `"MATE"` - Assembly mate
- `"DIMENSION"` - Dimension

## Error Handling

Always implement error checking:

```python
# Check if object exists
if model is None:
    print("No active document")
    return False

# Check if selection succeeded
result = model.Extension.SelectByID2("Face1", "FACE", 0, 0, 0, False, 0, None, 0)
if not result:
    print("Selection failed")
    return False

# Check if feature was created
feature = model.FeatureManager.FeatureExtrusion2(...)
if feature is None:
    print("Feature creation failed")
    return False
```

## Performance Optimization

For batch operations or complex scripts:

```python
# Disable screen updates
model.FeatureManager.EnableFeatureTree = False

# Suppress rebuild during operations
model.FeatureManager.EnableFeatureTreeWindow = False

# Perform operations...

# Re-enable and rebuild
model.FeatureManager.EnableFeatureTree = True
model.FeatureManager.EnableFeatureTreeWindow = True
model.EditRebuild3()
```

## Official Documentation Links

### Primary Resources

1. **SolidWorks API Help Documentation**
   - Access from: SolidWorks > Help > API Help
   - Local path: `C:\Program Files\SOLIDWORKS Corp\SOLIDWORKS\api\`

2. **SolidWorks API SDK**
   - Includes comprehensive examples and sample code
   - Installed with SolidWorks (check installation options)

3. **Online API Help Center**
   - https://help.solidworks.com/
   - Navigate to: Developer Tools and APIs > API Documentation

4. **SolidWorks API Forum**
   - https://forum.solidworks.com/community/api
   - Community support and discussions

### Learning Resources

1. **API Fundamentals Training**
   - SolidWorks Customer Portal (requires login)
   - Search for "API training" materials

2. **YouTube Tutorials**
   - SolidWorks official channel has API examples
   - Community contributors share practical examples

3. **Sample Code Repository**
   - Local: `C:\Program Files\SOLIDWORKS Corp\SOLIDWORKS\samples\`
   - Examples for VBA, VB.NET, C#, and C++

### Python-Specific Resources

1. **pywin32 Documentation**
   - https://github.com/mhammond/pywin32
   - COM interface documentation

2. **Python for SolidWorks Guide**
   - Community-maintained guides on using Python with SW API
   - Search GitHub for "solidworks python" examples

### C# .NET Resources

1. **SolidWorks Interop Assemblies**
   - Located in SolidWorks installation directory
   - Add reference in Visual Studio: `SolidWorks.Interop.sldworks`

2. **Add-in Development Guide**
   - Comprehensive guide for creating SolidWorks add-ins
   - Available in API Help documentation

## Version Compatibility

Different SolidWorks versions may have different API methods. Check compatibility:

```python
# Get SolidWorks version
revision = sw_app.RevisionNumber()
major_version = int(revision.split('.')[0])

print(f"SolidWorks Version: {major_version}")

# Conditional API usage
if major_version >= 2023:
    # Use newer API methods
    pass
else:
    # Use legacy methods
    pass
```

## Common API Constants

### Document Types
```python
swDocPART = 1
swDocASSEMBLY = 2
swDocDRAWING = 3
```

### Feature Types
```python
swFeatureBoss = "Boss"
swFeatureCut = "Cut"
swFeatureFillet = "Fillet"
swFeatureChamfer = "Chamfer"
swFeaturePattern = "LPattern"  # or "CPattern" for circular
```

### Selection Marks
```python
# Selection marks allow organizing multiple selections
MARK_PROFILE = 1
MARK_PATH = 2
MARK_GUIDE = 4
```

## Troubleshooting

### Common Issues

1. **COM Exception: Invalid class string**
   - SolidWorks not installed or not properly registered
   - Solution: Reinstall or repair SolidWorks installation

2. **Units Confusion**
   - Remember: API uses meters, not millimeters
   - Solution: Always convert mm to m (divide by 1000)

3. **Selection Failures**
   - Objects not found by name
   - Solution: Check exact name in feature tree (case-sensitive)

4. **Memory Leaks**
   - COM objects not released
   - Solution: Set objects to `None` when done (Python)

5. **Rebuild Failures**
   - Invalid geometry or dimension constraints
   - Solution: Check sketch validity and feature parameters

## Quick Reference Card

### Most Used Operations

```python
# Connect
sw = win32com.client.Dispatch("SldWorks.Application")

# New part
model = sw.NewDocument(template_path, 0, 0, 0)

# Select plane
model.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, None, 0)

# Start sketch
model.SketchManager.InsertSketch(True)

# Draw circle (radius in meters)
model.SketchManager.CreateCircleByRadius(0, 0, 0, 0.050)

# Extrude (depth in meters)
model.FeatureManager.FeatureExtrusion2(True, False, False, 0, 0, 0.040, 0, False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False)

# Rebuild
model.EditRebuild3()

# Save
model.SaveAs3(save_path, 0, 2)
```

This reference provides the essential API information for pump impeller design automation in SolidWorks.
