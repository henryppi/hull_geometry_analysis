import numpy as np
import matplotlib.pyplot as plt

from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_ThruSections
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_MakeEdge
from OCC.Core.gp import gp_Pnt, gp_Trsf, gp_Vec

# from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Display.SimpleGui import init_display

from OCC.Extend.DataExchange import write_step_file

from OCC.Core.gp import gp_Pnt, gp_Circ, gp_XOY
from OCC.Core.TColgp import TColgp_Array1OfPnt
from OCC.Core.Geom import Geom_BSplineCurve
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipeShell
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline

fn = './data_files/section_points_rudder_bezier_fit.dat'
points_raw = np.loadtxt(fn,delimiter=',',skiprows=0)
npoints = points_raw.shape[0]

# lower profile
# blunt tail edge
x_scale = 0.7
y_scale = 0.9
z = 60

points = np.copy(points_raw)
points[:,0] *=x_scale
points[:,1] *=y_scale
p0_tail = gp_Pnt(points[-1,0],points[-1,1],0)
p1_tail = gp_Pnt(points[0,0],points[0,1],0)
edge_tail = BRepBuilderAPI_MakeEdge(p0_tail, p1_tail).Edge()

point_list = []
for i in range(npoints):
    point_list.append(gp_Pnt(points[i,0],points[i,1],0))

point_arr = TColgp_Array1OfPnt(1, npoints)
for i, p in enumerate(point_list):
    point_arr.SetValue(i + 1, p)

# spline = Geom_BSplineCurve(point_arr)
# section_wire = BRepBuilderAPI_MakeWire(\
#                BRepBuilderAPI_MakeEdge(spline).Edge() ).Wire()

# 3. Interpolate the points to create a B-Spline curve geometry
bspline_geom = GeomAPI_PointsToBSpline(point_arr).Curve()

# 4. Turn the geometric curve into a topological Edge
spline_edge = BRepBuilderAPI_MakeEdge(bspline_geom).Edge()

# 5. Build the final topological Wire from the edge
profile_wire = BRepBuilderAPI_MakeWire(spline_edge,edge_tail).Wire()

# upper profile

# blunt tail edge
points_up = np.copy(points_raw)
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

lofter.AddWire(profile_wire)
lofter.AddWire(profile_wire_up)

# Optional settings
lofter.SetSmoothing(False)
lofter.CheckCompatibility(True)

# 3. Build the shape
lofter.Build()
variable_extrusion = lofter.Shape()

write_step_file(variable_extrusion, "./data_files/rudder.stp")

display, start_display, _, _ = init_display()


display.DisplayShape(
    variable_extrusion, 
    color='blue', 
    transparency=0.5, 
    update=True
)
start_display()


# fig, ax = plt.subplots(figsize=(8,8),frameon=False)
# ax.plot(points[:,0],points[:,1],'.-r')
# ax.axis('equal')
# # ax.set_axis_off()
# # plt.savefig('./images/rudder_fit_bezier.png',dpi=200)
# plt.show()  