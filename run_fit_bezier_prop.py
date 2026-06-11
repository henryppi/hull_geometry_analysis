import numpy as np
import matplotlib.pyplot as plt 
from my_modules import *

nsec = 12
sections = []
for i in range(nsec):
    fn = './data_files/sections_propeller_sec{:02d}.dat'.format(i)
    fmt='%.18e'
    sec = np.loadtxt(fn,delimiter='\t',skiprows=0)
    
    sections.append(sec)

rad = np.array([72.5,87.0,101.5,116.0,145.0,174.0,203.0,232.0,261.0,275.5,282.8,287.1])
tan_lead = np.array([36,54,71,86,111,124,124,105,59,18,-13,-41])
tan_tail = np.array([35,46,57,67,87,107,128,149,165,165,145,103])

sections_fit = []

for i in range(nsec):
    sec = sections[i]
    bbox = get_bbox(sec)
    print(bbox[1]-bbox[0], tan_lead[i]+tan_tail[i])

    x = bbox[0]
    y = 0.5*(bbox[2]+bbox[3])
    L = bbox[1]-bbox[0]
    t = bbox[3]-bbox[2]
    h = 0.05*t
    d = 0.03*t
    xfrac = 0.33
    wlead = 0.2
    wmid = 0.3

    para_init = [x,y,L,t,h,d,xfrac,wlead,wmid]

    para_fit,vert_fit,elem_fit = fit_bezier_prop_section(sec,para_init)
    # vert_fit[:,0] *= -1
    sections_fit.append(vert_fit)

if 1:
    
    fig, ax = plt.subplots()
    for i in range(nsec):
        sec = sections[i]
        ax.plot(sec[:,0]+tan_lead[i],sec[:,1]+rad[i])
        fit = sections_fit[i]
        ax.plot(fit[:,0]+tan_lead[i],fit[:,1]+rad[i])


    plt.axis('equal')
    plt.show()

# bbox = get_bbox(points)
# rect = np.array([[bbox[0],bbox[2]],\
#                     [bbox[1],bbox[2]],\
#                     [bbox[1],bbox[3]],\
#                     [bbox[0],bbox[3]],\
#                     [bbox[0],bbox[2]]])


# x = bbox[0]
# y = 0.5*(bbox[2]+bbox[3])
# L = bbox[1]-bbox[0]
# t = bbox[3]-bbox[2]
# h = 0.05*t
# xfrac = 0.33
# wlead = 0.25
# wfwd = 0.2
# wrev = 0.2
# wtail = 0.2

# para_init = [x,y,L,t,h,xfrac,wlead,wfwd,wrev,wtail]

# # para_coarse = coarse_fit_PCA(points)
# # x_coarse = para_coarse[0]
# # y_coarse = para_coarse[1]
# # L_coarse = para_coarse[3]
# # t_coarse = para_coarse[6]*L_coarse/100
# # print(x,x_coarse)
# # print(y,y_coarse)
# # print(L,L_coarse)
# # print(t,t_coarse)

# para_fit,vert_fit,elem_fit = fit_bezier_rudder(points,para_init)
# print(para_init)
# print(para_fit)


# fn_fit = './data_files/section_points_rudder_bezier_fit.dat'
# np.savetxt(fn_fit,vert_fit,delimiter=',')

# fig, ax = plt.subplots(figsize=(8,8),frameon=False)
# ax.plot(points[:,0],points[:,1],'.r')
# # ax.plot(nodesA[0,:],nodesA[1,:],'+--c',lw=1)
# # ax.plot(nodesB[0,:],nodesB[1,:],'+--m',lw=1)
# # ax.plot(nodesC[0,:],nodesC[1,:],'+--y',lw=1)
# # ax.plot(nodesD[0,:],nodesD[1,:],'+--b',lw=1)
# ax.plot(vert_fit[:,0],vert_fit[:,1],'-k',lw=2)
# ax.axis('equal')
# # ax.set_axis_off()
# plt.savefig('./images/rudder_fit_bezier.png',dpi=200)
# plt.show()  

