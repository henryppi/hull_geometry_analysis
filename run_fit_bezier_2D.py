import numpy as np
import matplotlib.pyplot as plt 
from my_modules import *

fn_points = './data_files/section_points_rudder.dat'

points = np.loadtxt(fn_points,delimiter='\t',skiprows=0)

bbox = get_bbox(points)
rect = np.array([[bbox[0],bbox[2]],\
                    [bbox[1],bbox[2]],\
                    [bbox[1],bbox[3]],\
                    [bbox[0],bbox[3]],\
                    [bbox[0],bbox[2]]])


x = bbox[0]
y = 0.5*(bbox[2]+bbox[3])
L = bbox[1]-bbox[0]
t = bbox[3]-bbox[2]
h = 0.05*t
xfrac = 0.33
wlead = 0.25
wfwd = 0.2
wrev = 0.2
wtail = 0.2

para_init = [x,y,L,t,h,xfrac,wlead,wfwd,wrev,wtail]

# para_coarse = coarse_fit_PCA(points)
# x_coarse = para_coarse[0]
# y_coarse = para_coarse[1]
# L_coarse = para_coarse[3]
# t_coarse = para_coarse[6]*L_coarse/100
# print(x,x_coarse)
# print(y,y_coarse)
# print(L,L_coarse)
# print(t,t_coarse)

para_fit,vert_fit,elem_fit = fit_bezier_rudder(points,para_init)
print(para_init)
print(para_fit)


fn_fit = './data_files/section_points_rudder_bezier_fit.dat'
np.savetxt(fn_fit,vert_fit,delimiter=',')

fig, ax = plt.subplots(figsize=(8,8),frameon=False)
ax.plot(points[:,0],points[:,1],'.r')
# ax.plot(nodesA[0,:],nodesA[1,:],'+--c',lw=1)
# ax.plot(nodesB[0,:],nodesB[1,:],'+--m',lw=1)
# ax.plot(nodesC[0,:],nodesC[1,:],'+--y',lw=1)
# ax.plot(nodesD[0,:],nodesD[1,:],'+--b',lw=1)
ax.plot(vert_fit[:,0],vert_fit[:,1],'-k',lw=2)
ax.axis('equal')
# ax.set_axis_off()
plt.savefig('./images/rudder_fit_bezier.png',dpi=200)
plt.show()  

