import vtk
import numpy as np

def computeNormal(v0,v1,v2):
	# edge vector
	e01 = v1-v0
	e10 = -e01
	e02 = v2-v0
	e20 = -e02
	e12 = v2-v1
	e21 = -e12
	ab = e01
	ac = e02
	normal = np.cross(ab,ac)
	magnitude = np.sqrt(np.vdot(normal,normal))
	normal = normal/magnitude
	return normal 
	
def computeRotationAxisAndAngle(v1,v2):
    w = np.cross(v1,v2)
    ang = np.arccos(np.vdot(v1,v2))
    return w,ang

def movePolyData(polydata,vec):
    translation = vtk.vtkTransform()
    translation.Translate(vec[0], vec[1], vec[2])

    polydata0 = vtk.vtkTransformPolyDataFilter()
    polydata0.SetInputData(polydata)
    polydata0.SetTransform(translation)
    polydata0.Update()
    return polydata0
    
def createCylinder(length,radius,pos,axis='x'):
    cyl = vtk.vtkCylinderSource()
    cyl.SetResolution(360)

    trans = vtk.vtkTransform()
    trans.Translate(pos)  # translate to starting point
    trans.Scale(length, radius, radius)  # scale along the height vector
    trans.RotateZ(-90.0)  # align cylinder to x axis

    tpd = vtk.vtkTransformPolyDataFilter()
    tpd.SetInputConnection(cyl.GetOutputPort())
    tpd.SetTransform(trans) 
    tpd.Update()

    return tpd

def getVtkTransformation4x4(mat):
    matVtk = vtk.vtkMatrix4x4()
    matVtk.SetElement(0,0,mat[0,0])
    matVtk.SetElement(0,1,mat[0,1])
    matVtk.SetElement(0,2,mat[0,2])
    matVtk.SetElement(0,3,mat[0,3])
    
    matVtk.SetElement(1,0,mat[1,0])
    matVtk.SetElement(1,1,mat[1,1])
    matVtk.SetElement(1,2,mat[1,2])
    matVtk.SetElement(1,3,mat[1,3])
    
    matVtk.SetElement(2,0,mat[2,0])
    matVtk.SetElement(2,1,mat[2,1])
    matVtk.SetElement(2,2,mat[2,2])
    matVtk.SetElement(2,3,mat[2,3])
    
    matVtk.SetElement(3,0,mat[3,0])
    matVtk.SetElement(3,1,mat[3,1])
    matVtk.SetElement(3,2,mat[3,2])
    matVtk.SetElement(3,3,mat[3,3])
    
    transVtk = vtk.vtkTransform()
    transVtk.SetMatrix(matVtk)
    return transVtk

def get4x4Transformation(axisPoint,axisDir,angle):
    a = axisPoint[0]
    b = axisPoint[1]
    c = axisPoint[2]
    
    u = axisDir[0]
    v = axisDir[1]
    w = axisDir[2]
    
    ca = np.cos(angle)
    sa = np.sin(angle)
    
    mat = np.matrix(np.eye(4),float)
    mat[0,0] = u**2 + ( v**2 + w**2 ) * ca
    mat[1,1] = v**2 + ( u**2 + w**2 ) * ca
    mat[2,2] = w**2 + ( u**2 + v**2 ) * ca
    
    mat[0,1] = u * v * ( 1-ca ) - w * sa
    mat[0,2] = u * w * ( 1-ca ) + v * sa
    mat[1,2] = v * w * ( 1-ca ) - u * sa
    
    mat[1,0] = u * v * ( 1-ca ) + w * sa
    mat[2,0] = u * w * ( 1-ca ) - v * sa
    mat[2,1] = v * w * ( 1-ca ) + u * sa
    
    mat[0,3] = ( a * ( v**2 + w**2 ) - u * ( b*v +c*w ) ) * ( 1 - ca ) + ( b*w - c*v ) * sa
    mat[1,3] = ( b * ( u**2 + w**2 ) - v * ( a*u +c*w ) ) * ( 1 - ca ) + ( c*u - a*w ) * sa
    mat[2,3] = ( c * ( u**2 + v**2 ) - w * ( a*u +b*v ) ) * ( 1 - ca ) + ( a*v - b*u ) * sa
    
    mat[3,3] = 1.0
    
    return mat
    
def applyTransformation(vtkObj,transVtk):
    tpd1 = vtk.vtkTransformPolyDataFilter()
    tpd1.SetInputConnection(vtkObj.GetOutputPort())
    tpd1.SetTransform(transVtk) 
    tpd1.Update()
    
#     appendF = vtk.vtkAppendPolyData()
#     appendF.AddInput(tpd1.GetOutput())
#     return appendF
    return tpd1


def main():
    colors = vtk.vtkNamedColors()

    filename = './CAD/rudder_decimate.stl'

    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.Update()

    polydata = reader.GetOutput()
    
    p1 = np.array([-15.078249754375094,
        -21.625060092912193,
        432.98394098679887])
        
    p2 = np.array([-36.37472085843895,
        15.048355455262683,
        439.83903703936426])
        
    p3 = np.array([-42.265306795481436,
        -22.06301869611415,
        437.3713656950272])
    
    # move to center point p4
    p4 = np.array([-28.952583139973388,
    -4.8756540987,
    438.8229453669872])
    vec = -p4
    polydata0 = movePolyData(polydata,vec)
        
    map2 = vtk.vtkPolyDataMapper()
    map2.SetInputConnection(polydata0.GetOutputPort())
    act2 = vtk.vtkActor()
    act2.GetProperty().SetColor(colors.GetColor3d('Red')) 
    act2.SetMapper(map2)

    # align/rotate with x-axis
    nn = computeNormal(p1,p2,p3)
    xDir = np.array([1,0,0],float)
    axisDir, angle = computeRotationAxisAndAngle(nn,xDir)
    axisPoint = np.array([0,0,0],float)

    mat4x4 = get4x4Transformation(axisPoint,axisDir,angle)
    transVtk = getVtkTransformation4x4(mat4x4)
    prop_aligned = applyTransformation(polydata0,transVtk)
    
    map3 = vtk.vtkPolyDataMapper()
    map3.SetInputConnection(prop_aligned.GetOutputPort())
    act3 = vtk.vtkActor()
    act3.GetProperty().SetColor(colors.GetColor3d('Yellow')) 
    act3.SetMapper(map3)
    
    length = 100.0
    radius = 150.0
    pos = np.array([50,0,0],float)
    axis_dir = 'x'
    cyl = createCylinder(length,radius,pos,axis_dir)
    # create cylinder
   
    cyl_map = vtk.vtkPolyDataMapper()
    cyl_map.SetInputConnection(cyl.GetOutputPort())
    cyl_act = vtk.vtkActor()
    cyl_act.GetProperty().SetOpacity(0.2)
    cyl_act.GetProperty().SetColor(colors.GetColor3d('Green')) 
    cyl_act.SetMapper(cyl_map)


    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    
    actor = vtk.vtkActor()
    actor.GetProperty().SetLineWidth(2)
    actor.GetProperty().SetOpacity(0.7)
    actor.GetProperty().SetPointSize(10)
    actor.GetProperty().SetDiffuse(0.8)
    actor.GetProperty().SetColor(1.0,0.357,0.0) #orange
    actor.GetProperty().SetSpecular(0.3)
    actor.GetProperty().SetSpecularPower(60.0)
    actor.SetMapper(mapper)
    
    ren = vtk.vtkRenderer()
#     ren.AddActor(actor)
#     ren.AddActor(act2)
    ren.AddActor(act3)
    ren.AddActor(cyl_act)
    renWin = vtk.vtkRenderWindow()
    renWin.SetWindowName('Propeller Analyzer')
    renWin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    
    # add CSYS
    axes = vtk.vtkAxesActor();
    widget = vtk.vtkOrientationMarkerWidget()
    rgba = [0] * 4
    colors.GetColor('Carrot', rgba)
    widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
    widget.SetOrientationMarker(axes)
    widget.SetInteractor(iren)
    widget.SetViewport(0.0, 0.0, 0.4, 0.4)
    widget.SetEnabled(1)
    widget.InteractiveOn()
    # end add CSYS    
    
    ren.SetBackground(colors.GetColor3d('SlateGray'))
    ren.GetActiveCamera().Azimuth(50)
    ren.GetActiveCamera().Elevation(-30)
    ren.ResetCamera()
    renWin.Render()
    iren.Start()

if __name__ == '__main__':
    main()