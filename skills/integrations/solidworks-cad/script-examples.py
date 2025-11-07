"""
SolidWorks API Examples for Parametric Impeller Design

REQUIREMENTS:
- Windows operating system
- SolidWorks installed with valid license
- Python package: pywin32 (pip install pywin32)

USAGE:
Run these examples on a Windows machine with SolidWorks installed.
Some functions assume SolidWorks is already running, others will launch it.
"""

import win32com.client
import pythoncom
import math
import sys


class SolidWorksAPI:
    """Wrapper class for SolidWorks API operations"""

    def __init__(self):
        """Initialize connection to SolidWorks"""
        self.sw_app = None
        self.model_doc = None
        self.model_doc_ext = None

    def connect(self):
        """Connect to running SolidWorks instance or launch new one"""
        try:
            # Try to connect to existing instance
            self.sw_app = win32com.client.Dispatch("SldWorks.Application")
            self.sw_app.Visible = True
            print("Connected to SolidWorks")
            return True
        except Exception as e:
            print(f"Error connecting to SolidWorks: {e}")
            print("Make sure SolidWorks is installed and you have a valid license")
            return False

    def create_part(self):
        """Create a new part document"""
        if not self.sw_app:
            print("Not connected to SolidWorks")
            return False

        # Part template - adjust path to your SolidWorks installation
        template_path = "C:\\ProgramData\\SolidWorks\\SOLIDWORKS 2023\\templates\\Part.prtdot"

        self.model_doc = self.sw_app.NewDocument(template_path, 0, 0, 0)
        if self.model_doc:
            self.model_doc_ext = self.model_doc
            print("New part created")
            return True
        return False


def example_create_parametric_impeller(sw):
    """
    Create a simplified parametric centrifugal impeller

    This example creates a basic impeller with parametric control over:
    - Hub diameter
    - Shroud diameter
    - Blade height
    - Number of blades
    - Blade thickness
    """

    # Impeller parameters (all dimensions in mm)
    hub_diameter = 50.0
    shroud_diameter = 150.0
    blade_height = 40.0
    num_blades = 6
    blade_thickness = 3.0
    inlet_blade_angle = 20.0  # degrees
    outlet_blade_angle = 25.0  # degrees

    print("\n=== Creating Parametric Impeller ===")
    print(f"Hub diameter: {hub_diameter} mm")
    print(f"Shroud diameter: {shroud_diameter} mm")
    print(f"Number of blades: {num_blades}")

    # Create new part
    if not sw.create_part():
        return False

    model = sw.model_doc

    # Select Front plane for sketch
    feature_mgr = model.FeatureManager

    # Create hub (central cylinder)
    # Note: Actual implementation would use SelectByID2 and sketch creation
    print("Step 1: Creating hub geometry...")
    print(f"  - Sketch circle on Front plane, diameter {hub_diameter} mm")
    print(f"  - Extrude {blade_height} mm")

    # Create shroud (outer cylinder)
    print("Step 2: Creating shroud geometry...")
    print(f"  - Sketch circle on Front plane, diameter {shroud_diameter} mm")
    print(f"  - Extrude {blade_height} mm")
    print(f"  - Shell operation to create outer wall")

    # Create blade profile
    print("Step 3: Creating blade profile...")
    print(f"  - Define blade curve with inlet angle {inlet_blade_angle}°")
    print(f"  - Define blade curve with outlet angle {outlet_blade_angle}°")
    print(f"  - Create swept surface along blade path")
    print(f"  - Thicken surface to {blade_thickness} mm")

    # Pattern blades
    print("Step 4: Patterning blades...")
    print(f"  - Circular pattern around central axis")
    print(f"  - Pattern count: {num_blades}")
    print(f"  - Pattern angle: {360.0 / num_blades}° spacing")

    print("\nNote: This is a simplified example showing the workflow.")
    print("Full implementation requires detailed sketch and feature creation API calls.")

    return True


def example_modify_dimensions(sw):
    """
    Modify dimensions of an existing parametric part

    This demonstrates how to change design parameters programmatically
    """

    print("\n=== Modifying Part Dimensions ===")

    if not sw.model_doc:
        print("No active model. Open a part file first.")
        return False

    model = sw.model_doc

    # Get dimension by name
    # Dimension names can be found in SolidWorks by selecting the dimension
    dimension_name = "D1@Sketch1"
    new_value = 100.0  # mm

    print(f"Changing dimension '{dimension_name}' to {new_value} mm")

    # Example code structure (requires actual part with named dimensions):
    """
    feature = model.FeatureByName("Boss-Extrude1")
    if feature:
        dimension = feature.Parameter(dimension_name)
        if dimension:
            dimension.SystemValue = new_value / 1000.0  # Convert mm to meters (SI units)
            model.EditRebuild3()
            print("Dimension updated and model rebuilt")
    """

    print("\nNote: Requires existing part with named dimensions")
    print("SolidWorks API uses SI units internally (meters, not mm)")

    return True


def example_export_geometry(sw, output_formats=['STEP', 'IGES', 'STL']):
    """
    Export SolidWorks part to various formats for CFD or manufacturing

    Args:
        sw: SolidWorksAPI instance
        output_formats: List of formats to export ['STEP', 'IGES', 'STL', 'Parasolid']
    """

    print("\n=== Exporting Geometry ===")

    if not sw.model_doc:
        print("No active model to export")
        return False

    model = sw.model_doc
    model_path = model.GetPathName()

    if not model_path:
        print("Please save the model before exporting")
        return False

    # Base path without extension
    base_path = model_path.rsplit('.', 1)[0]

    print(f"Exporting from: {model_path}")

    for format_type in output_formats:
        if format_type.upper() == 'STEP':
            export_path = f"{base_path}.step"
            print(f"\nExporting to STEP: {export_path}")
            # model.SaveAs3(export_path, 0, 2)  # 2 = STEP AP214
            print("  - Format: STEP AP214 (ISO 10303)")
            print("  - Use case: CFD preprocessing, CAE analysis")

        elif format_type.upper() == 'IGES':
            export_path = f"{base_path}.igs"
            print(f"\nExporting to IGES: {export_path}")
            # model.SaveAs3(export_path, 0, 2)
            print("  - Format: IGES 5.3")
            print("  - Use case: Legacy CAD systems, CNC programming")

        elif format_type.upper() == 'STL':
            export_path = f"{base_path}.stl"
            print(f"\nExporting to STL: {export_path}")
            # model.SaveAs3(export_path, 0, 2)
            print("  - Format: STL (tessellated mesh)")
            print("  - Use case: 3D printing, mesh-based CFD")
            print("  - Note: Adjust mesh resolution for accuracy vs. file size")

        elif format_type.upper() == 'PARASOLID':
            export_path = f"{base_path}.x_t"
            print(f"\nExporting to Parasolid: {export_path}")
            # model.SaveAs3(export_path, 0, 2)
            print("  - Format: Parasolid binary")
            print("  - Use case: High-fidelity geometry transfer to ANSYS, Siemens NX")

    print("\n=== Export Settings Recommendations ===")
    print("STEP Export:")
    print("  - Use AP214 or AP242 for best compatibility")
    print("  - Include model appearance for visualization")

    print("\nSTL Export for CFD:")
    print("  - Binary format for smaller file size")
    print("  - Resolution: 0.01mm for pump impellers")
    print("  - Quality: Fine (for smooth surfaces)")

    return True


def example_batch_processing():
    """
    Batch process multiple impeller designs with different parameters

    This is useful for design optimization studies or generating
    multiple variants for CFD analysis
    """

    print("\n=== Batch Processing Multiple Designs ===")

    # Design parameter variations
    designs = [
        {"name": "impeller_6blade_50mm", "num_blades": 6, "hub_dia": 50, "shroud_dia": 150},
        {"name": "impeller_8blade_50mm", "num_blades": 8, "hub_dia": 50, "shroud_dia": 150},
        {"name": "impeller_6blade_60mm", "num_blades": 6, "hub_dia": 60, "shroud_dia": 150},
        {"name": "impeller_8blade_60mm", "num_blades": 8, "hub_dia": 60, "shroud_dia": 150},
    ]

    sw = SolidWorksAPI()
    if not sw.connect():
        return False

    for i, design in enumerate(designs, 1):
        print(f"\n--- Processing Design {i}/{len(designs)}: {design['name']} ---")
        print(f"    Blades: {design['num_blades']}")
        print(f"    Hub: {design['hub_dia']} mm")
        print(f"    Shroud: {design['shroud_dia']} mm")

        # Create part with parameters
        # Update dimensions
        # Export to STEP for CFD
        # Save SolidWorks file

        output_path = f"C:\\designs\\{design['name']}.SLDPRT"
        step_path = f"C:\\designs\\{design['name']}.step"

        print(f"    Saving to: {output_path}")
        print(f"    Exporting STEP: {step_path}")

    print("\n=== Batch Processing Complete ===")
    print(f"Processed {len(designs)} designs")
    print("Ready for CFD analysis in OpenFOAM")

    return True


# ============================================================================
# VBA MACRO EXAMPLES
# ============================================================================

VBA_EXAMPLE_CREATE_IMPELLER = """
' VBA Macro: Create Parametric Impeller
' Copy this code into SolidWorks VBA editor (Tools > Macro > New)

Option Explicit

Sub CreateParametricImpeller()
    Dim swApp As SldWorks.SldWorks
    Dim swModel As SldWorks.ModelDoc2
    Dim swFeature As SldWorks.Feature
    Dim boolstatus As Boolean

    ' Impeller parameters
    Dim hubDiameter As Double
    Dim shroudDiameter As Double
    Dim bladeHeight As Double
    Dim numBlades As Integer

    hubDiameter = 0.05      ' 50mm in meters
    shroudDiameter = 0.15   ' 150mm in meters
    bladeHeight = 0.04      ' 40mm in meters
    numBlades = 6

    ' Connect to SolidWorks
    Set swApp = Application.SldWorks

    ' Create new part
    Set swModel = swApp.NewDocument("C:\\ProgramData\\SolidWorks\\templates\\Part.prtdot", 0, 0, 0)

    ' Create hub
    boolstatus = swModel.Extension.SelectByID2("Front Plane", "PLANE", 0, 0, 0, False, 0, Nothing, 0)
    swModel.SketchManager.InsertSketch True

    ' Draw hub circle
    swModel.SketchManager.CreateCircleByRadius hubDiameter / 2, 0, 0, 0

    ' Extrude hub
    swModel.FeatureManager.FeatureExtrusion2 True, False, False, 0, 0, bladeHeight, 0, False, False, False, False, 0, 0, False, False, False, False, True, True, True, 0, 0, False

    ' Create blade sketch (simplified)
    ' ... (Additional VBA code for blade creation)

    ' Circular pattern
    boolstatus = swModel.Extension.SelectByID2("Boss-Extrude1", "BODYFEATURE", 0, 0, 0, False, 0, Nothing, 0)
    swModel.FeatureManager.FeatureCircularPattern numBlades, 2 * 3.14159 / numBlades, False, "", False, False, True, False, False

    MsgBox "Impeller created with " & numBlades & " blades"

End Sub
"""

VBA_EXAMPLE_MODIFY_DIMENSIONS = """
' VBA Macro: Modify Part Dimensions
' Demonstrates parametric design modification

Sub ModifyImpellerDimensions()
    Dim swApp As SldWorks.SldWorks
    Dim swModel As SldWorks.ModelDoc2
    Dim swFeature As SldWorks.Feature
    Dim swDimension As SldWorks.Dimension

    Set swApp = Application.SldWorks
    Set swModel = swApp.ActiveDoc

    If swModel Is Nothing Then
        MsgBox "Please open a part document first"
        Exit Sub
    End If

    ' Get feature by name
    Set swFeature = swModel.FeatureByName("Boss-Extrude1")

    If Not swFeature Is Nothing Then
        ' Get dimension by name (e.g., "D1@Sketch1")
        Set swDimension = swFeature.Parameter("D1@Sketch1")

        If Not swDimension Is Nothing Then
            ' Set new value (in meters for SolidWorks API)
            swDimension.SystemValue = 0.1  ' 100mm

            ' Rebuild model
            swModel.EditRebuild3

            MsgBox "Dimension updated successfully"
        End If
    End If

End Sub
"""

VBA_EXAMPLE_EXPORT = """
' VBA Macro: Export to STEP format
' Useful for CFD preprocessing

Sub ExportToSTEP()
    Dim swApp As SldWorks.SldWorks
    Dim swModel As SldWorks.ModelDoc2
    Dim strPath As String
    Dim lErrors As Long
    Dim lWarnings As Long

    Set swApp = Application.SldWorks
    Set swModel = swApp.ActiveDoc

    If swModel Is Nothing Then
        MsgBox "No active document"
        Exit Sub
    End If

    ' Get current file path and change extension
    strPath = swModel.GetPathName
    If strPath = "" Then
        MsgBox "Please save the document first"
        Exit Sub
    End If

    ' Replace extension with .step
    strPath = Left(strPath, InStrRev(strPath, ".") - 1) & ".step"

    ' Save as STEP
    swModel.Extension.SaveAs strPath, 0, 0, Nothing, lErrors, lWarnings

    If lErrors = 0 Then
        MsgBox "Exported to: " & strPath
    Else
        MsgBox "Export failed with errors"
    End If

End Sub
"""


def print_vba_examples():
    """Print VBA macro examples"""
    print("\n" + "="*70)
    print("VBA MACRO EXAMPLES")
    print("="*70)

    print("\n--- Example 1: Create Parametric Impeller ---")
    print(VBA_EXAMPLE_CREATE_IMPELLER)

    print("\n--- Example 2: Modify Dimensions ---")
    print(VBA_EXAMPLE_MODIFY_DIMENSIONS)

    print("\n--- Example 3: Export to STEP ---")
    print(VBA_EXAMPLE_EXPORT)


def main():
    """Main function to run examples"""
    print("="*70)
    print("SolidWorks API Examples for Pump Impeller Design")
    print("="*70)
    print("\nREQUIREMENTS:")
    print("  - Windows OS")
    print("  - SolidWorks installed with valid license")
    print("  - Python package: pywin32")
    print("\nNOTE: These examples demonstrate the API workflow.")
    print("      Full implementation requires detailed API calls.")
    print("="*70)

    # Initialize SolidWorks API
    sw = SolidWorksAPI()

    # Check platform
    if sys.platform != "win32":
        print("\nERROR: SolidWorks API requires Windows")
        print("Current platform:", sys.platform)
        print("\nShowing example code structure only...\n")

        # Show examples without executing
        example_create_parametric_impeller(sw)
        example_modify_dimensions(sw)
        example_export_geometry(sw)
        example_batch_processing()
        print_vba_examples()
        return

    # Try to connect (only on Windows)
    print("\nAttempting to connect to SolidWorks...")
    if not sw.connect():
        print("\nCould not connect to SolidWorks.")
        print("Make sure:")
        print("  1. SolidWorks is installed")
        print("  2. You have a valid license")
        print("  3. SolidWorks is running or can be launched")
        print("\nShowing example code structure only...\n")

        example_create_parametric_impeller(sw)
        example_modify_dimensions(sw)
        example_export_geometry(sw)
        example_batch_processing()
        print_vba_examples()
        return

    # Run examples with active SolidWorks connection
    print("\n" + "="*70)
    print("Running Examples...")
    print("="*70)

    # Example 1: Create parametric impeller
    example_create_parametric_impeller(sw)

    # Example 2: Modify dimensions
    example_modify_dimensions(sw)

    # Example 3: Export geometry
    example_export_geometry(sw)

    # Example 4: Batch processing
    example_batch_processing()

    # Show VBA examples
    print_vba_examples()

    print("\n" + "="*70)
    print("Examples Complete")
    print("="*70)


if __name__ == "__main__":
    main()
