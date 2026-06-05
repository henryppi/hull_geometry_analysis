from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.TopAbs import TopAbs_SOLID
# from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
import os

step_reader = STEPControl_Reader() # Create a STEPControl_Reader instance

#path and filename
path_str = './CAD/'
# fn = '20m_Craft.stp'
fn = 'Full_Hull_Innovator.stp'
step_file_path = path_str+fn
output_file_path = path_str+fn[:-4]+'.stl'

status = step_reader.ReadFile(step_file_path)

if status == 1:
    step_reader.TransferRoot()
    main_shape = step_reader.Shape()
    mesh = BRepMesh_IncrementalMesh(main_shape, 0.01,True,0.1,True)
    mesh.Perform()
    
    if not mesh.IsDone():
        print("Error: Meshing failed.")
    else:
        writer = StlAPI_Writer()
        writer.SetASCIIMode(True)
        writer.Write(mesh.Shape(), output_file_path)
        print(f"Shape successfully exported to {output_file_path}")