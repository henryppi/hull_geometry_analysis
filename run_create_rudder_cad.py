import numpy as np
import matplotlib.pyplot as plt

from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from OCC.Core.gp import gp_Pnt
# from OCC.Core.gp import gp_Trsf, gp_Vec

# from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Display.SimpleGui import init_display

from OCC.Extend.DataExchange import write_step_file

from OCC.Core.gp import gp_Pnt, gp_XOY
from OCC.Core.TColgp import TColgp_Array1OfPnt
# from OCC.Core.Geom import Geom_BSplineCurve
# from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipeShell
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline

from my_modules import get_bbox

fn = './data_files/section_points_rudder_bezier_fit.dat'
points_raw = np.loadtxt(fn,delimiter=',',skiprows=0)
npoints = points_raw.shape[0]

bbox = get_bbox(points_raw)
print('L = ',bbox[1]-bbox[0])
print('t = ',bbox[3]-bbox[2])

z = 600

L_raw = bbox[1]-bbox[0]

points_unit = np.copy(points_raw)/L_raw*1000

bbox2 = get_bbox(points_unit)
print('L_unit = ',bbox2[1]-bbox2[0])
print('t_unit = ',bbox2[3]-bbox2[2])

L_unit = 1000

# lower profile
# blunt tail edge

L_bottom = 50+220
scale_bottom = L_bottom/L_unit
points_bottom = np.copy(points_unit)
# points_bottom *=scale_bottom
bbox3 = get_bbox(points_bottom)
print('L_bottom = ',bbox3[1]-bbox3[0])
print('t_bottom = ',bbox3[3]-bbox3[2])

p0_tail_bottom = gp_Pnt(points_bottom[-1,0],points_bottom[-1,1],0)
p1_tail_bottom = gp_Pnt(points_bottom[0,0],points_bottom[0,1],0)
edge_tail_bottom = BRepBuilderAPI_MakeEdge(p0_tail_bottom, p1_tail_bottom).Edge()

point_list_bottom = []
for i in range(npoints):
    point_list_bottom.append(gp_Pnt(points_bottom[i,0],points_bottom[i,1],0))

point_arr_bottom = TColgp_Array1OfPnt(1, npoints)
for i, p in enumerate(point_list_bottom):
    point_arr_bottom.SetValue(i + 1, p)

# 3. Interpolate the points to create a B-Spline curve geometry
bspline_geom_bottom = GeomAPI_PointsToBSpline(point_arr_bottom).Curve()

# 4. Turn the geometric curve into a topological Edge
spline_edge_bottom = BRepBuilderAPI_MakeEdge(bspline_geom_bottom).Edge()

# 5. Build the final topological Wire from the edge
profile_wire_bottom = BRepBuilderAPI_MakeWire(spline_edge_bottom,edge_tail_bottom).Wire()

# upper profile
# blunt tail edge

L_up = 125+275
scale_up = L_up/L_unit
points_up = np.copy(points_unit)
# points_up *=scale_up

bbox4 = get_bbox(points_up)
print('L_up = ',bbox4[1]-bbox4[0])
print('t_up = ',bbox4[3]-bbox4[2])

points_up = np.copy(points_unit)
p0_tail_up = gp_Pnt(points_up[-1,0],points_up[-1,1],z)
p1_tail_up = gp_Pnt(points_up[0,0],points_up[0,1],z)
edge_tail_up = BRepBuilderAPI_MakeEdge(p0_tail_up, p1_tail_up).Edge()

point_list_up = []
for i in range(npoints):
    point_list_up.append(gp_Pnt(points_up[i,0],points_up[i,1],z))

point_arr_up = TColgp_Array1OfPnt(1, npoints)
for i, p in enumerate(point_list_up):
    point_arr_up.SetValue(i + 1, p)

# spline = Geom_BSplineCurve(point_arr)
# section_wire = BRepBuilderAPI_MakeWire(\
#                BRepBuilderAPI_MakeEdge(spline).Edge() ).Wire()

# 3. Interpolate the points to create a B-Spline curve geometry
bspline_geom_up = GeomAPI_PointsToBSpline(point_arr_up).Curve()

# 4. Turn the geometric curve into a topological Edge
spline_edge_up = BRepBuilderAPI_MakeEdge(bspline_geom_up).Edge()

# 5. Build the final topological Wire from the edge
profile_wire_up = BRepBuilderAPI_MakeWire(spline_edge_up,edge_tail_up).Wire()


# 2. Perform the Loft (Variable Extrusion)
# BRepOffsetAPI_ThruSections(isSolid, ruled, pres3d)
lofter = BRepOffsetAPI_ThruSections(True, False) # True for solid, False for smooth loft

lofter.AddWire(profile_wire_bottom)
lofter.AddWire(profile_wire_up)

# Optional settings
lofter.SetSmoothing(False)
lofter.CheckCompatibility(True)

# 3. Build the shape
lofter.Build()
variable_extrusion = lofter.Shape()

if 1:
    write_step_file(variable_extrusion, "./data_files/rudder.stp")

if 1:
    display, start_display, _, _ = init_display()


    display.DisplayShape(
        variable_extrusion, 
        color='blue', 
        transparency=0.5, 
        update=True
    )
    start_display()

if 0:
    fig, ax = plt.subplots(figsize=(8,8),frameon=False)
    ax.plot(points_unit[:,0],points_unit[:,1],'.-r')
    ax.axis('equal')
    # ax.set_axis_off()
    # plt.savefig('./images/rudder_fit_bezier.png',dpi=200)
    plt.show()  