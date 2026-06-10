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
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib

from my_modules import get_bbox

fn = './data_files/section_points_rudder_bezier_fit.dat'
points_raw = np.loadtxt(fn,delimiter=',',skiprows=0)
npoints = points_raw.shape[0]

bbox_raw = get_bbox(points_raw)
L_raw = bbox_raw[1]-bbox_raw[0]
t_raw = bbox_raw[3]-bbox_raw[2]
print('L_raw = ',L_raw)
print('t_raw = ',t_raw)

z = 600
L_top_fwd = 125
L_top_rev = 275
L_bot_fwd = 50
L_bot_rev = 220

L_unit = 1000
L_top = L_top_fwd + L_top_rev
L_bot = L_bot_fwd + L_bot_rev

points_unit = np.copy(points_raw)
points_unit[:,0] += -bbox_raw[0]
points_unit = np.copy(points_unit)/L_raw*L_unit

bbox_unit = get_bbox(points_unit)
print('L_unit = ',bbox_unit[1]-bbox_unit[0])
print('t_unit = ',bbox_unit[3]-bbox_unit[2])

scale_top = L_top/L_unit
scale_bot = L_bot/L_unit

points_top = np.copy(points_unit)
points_top *= scale_top
points_top[:,0] += -L_top_fwd

points_bot = np.copy(points_unit)
points_bot[:,0] *= scale_bot
points_bot[:,1] *= scale_top # keep top thickness
points_bot[:,0] += -L_bot_fwd

#==================
# upper profile
p0_tail_top = gp_Pnt(points_top[-1,0],points_top[-1,1],z)
p1_tail_top = gp_Pnt(points_top[0,0],points_top[0,1],z)
edge_tail_top = BRepBuilderAPI_MakeEdge(p0_tail_top, p1_tail_top).Edge()

point_list_top = []
for i in range(npoints):
    point_list_top.append(gp_Pnt(points_top[i,0],points_top[i,1],z))

point_arr_top = TColgp_Array1OfPnt(1, npoints)
for i, p in enumerate(point_list_top):
    point_arr_top.SetValue(i + 1, p)

bspline_geom_top = GeomAPI_PointsToBSpline(point_arr_top).Curve()
spline_edge_top = BRepBuilderAPI_MakeEdge(bspline_geom_top).Edge()
profile_wire_top = BRepBuilderAPI_MakeWire(spline_edge_top,edge_tail_top).Wire()

#==================
# lower profile
p0_tail_bot = gp_Pnt(points_bot[-1,0],points_bot[-1,1],0)
p1_tail_bot = gp_Pnt(points_bot[0,0],points_bot[0,1],0)
edge_tail_bot = BRepBuilderAPI_MakeEdge(p0_tail_bot, p1_tail_bot).Edge()

point_list_bot = []
for i in range(npoints):
    point_list_bot.append(gp_Pnt(points_bot[i,0],points_bot[i,1],0))

point_arr_bot = TColgp_Array1OfPnt(1, npoints)
for i, p in enumerate(point_list_bot):
    point_arr_bot.SetValue(i + 1, p)

bspline_geom_bot = GeomAPI_PointsToBSpline(point_arr_bot).Curve()
spline_edge_bot = BRepBuilderAPI_MakeEdge(bspline_geom_bot).Edge()
profile_wire_bot = BRepBuilderAPI_MakeWire(spline_edge_bot,edge_tail_bot).Wire()

#==================
# lofting extrusion
lofter = BRepOffsetAPI_ThruSections(True, False) # True for solid, False for smooth loft
lofter.AddWire(profile_wire_bot)
lofter.AddWire(profile_wire_top)
lofter.SetSmoothing(False)
lofter.CheckCompatibility(True)
lofter.Build()
rudder_shape = lofter.Shape()


bbox = Bnd_Box()
brepbndlib.Add(rudder_shape, bbox)
xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
print(f"L xmin xmax: [{xmax-xmin}, {xmin}, {xmax}]")
print(f"t ymin ymax: [{ymax-ymin}, {ymin}, {ymax}]")
print(f"H zmin zmax: [{zmax-zmin}, {zmin}, {zmax}]")

if 1:
    write_step_file(rudder_shape, "./data_files/rudder.stp")

if 1:
    display, start_display, _, _ = init_display()
    display.DisplayShape(
        rudder_shape, 
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
    # plt.savefig('./images/rudder_secion.png',dpi=200)
    plt.show()  