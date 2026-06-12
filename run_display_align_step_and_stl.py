import sys
import numpy as np
from my_modules import RigidBody
if sys.platform == "linux" or sys.platform == "linux2":
    print("Running on Linux")
    run_os = 'linux'
elif sys.platform == "darwin":
    print("Running on macOS")
    run_os = 'macos'

from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Display.SimpleGui import init_display
from OCC.Core.IFSelect import IFSelect_RetDone

from OCC.Core.gp import gp_Trsf
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform

from OCC.Core.StlAPI import StlAPI_Reader
from OCC.Core.TopoDS import TopoDS_Shape
import vtk

fn_stp = './data_files/rudder.stp'
fn_stl = './CAD/rudder.stl'
fn_stl_dec = fn_stl[:-4]+'_decimate.stl'

# Initialize and use the basic STEP reader
step_reader = STEPControl_Reader()
if step_reader.ReadFile(fn_stp) == IFSelect_RetDone:
    step_reader.TransferRoots()
    shape_stp = step_reader.OneShape()


if 0:
    stl_vtk_reader = vtk.vtkSTLReader()
    stl_vtk_reader.SetFileName(fn_stl)
    stl_vtk_reader.Update()

    # Initialize the quadric decimation filter
    decimate = vtk.vtkQuadricDecimation()
    decimate.SetInputConnection(stl_vtk_reader.GetOutputPort())
    decimate.SetTargetReduction(0.9) # Reduce triangle count by 90%
    decimate.Update()

    
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(fn_stl_dec)
    writer.SetInputConnection(decimate.GetOutputPort())
    writer.SetFileTypeToBinary()
    writer.Write()


stl_reader = StlAPI_Reader()
shape_stl = TopoDS_Shape()
success = stl_reader.Read(shape_stl, fn_stl_dec)

def transform(shape,axis,angle,vec):
    rb = RigidBody()
    q = rb.quaternionAxisAngle(axis,angle)
    R = rb.R(q)
    T = gp_Trsf()
    T.SetValues( R[0,0], R[0,1], R[0,2], vec[0],
                 R[1,0], R[1,1], R[1,2], vec[1],
                 R[2,0], R[2,1], R[2,2], vec[2])
    
    transformer = BRepBuilderAPI_Transform(shape, T, True)
    return transformer.Shape()

axis = np.array([0,1,0])
angle = -90*np.pi/180
vec = np.array([0,0,0],float)

shape_stl = transform(shape_stl,axis,angle,vec)

if run_os == 'macos':
    display, start_display, _, _ = init_display(backend_str="pyside6")
else:
    display, start_display, _, _ = init_display()

display.DisplayShape(
    shape_stp, 
    color='blue', 
    transparency=0.5, 
    update=True
)

display.DisplayShape(
    shape_stl, 
    color='yellow', 
    transparency=0.5, 
    update=True)

start_display()