import numpy as np 
import matplotlib.pyplot as plt
from my_modules import cubic_bezier, quadratic_bezier
nspline = 1

x = 0
y = 0

L = 1.0
xfrac = 0.5
t = 0.1
d = 0.03
h = 0.02
wlead = 0.2
wmid = 0.3

x0 = x + L
y0 = y + 0.5*h

x1 = x + xfrac*L + wmid*L
y1 = y + 0.5*t + d

x2 = x + xfrac*L
y2 = y + 0.5*t + d

x3 = x + xfrac*L - wmid*L
y3 = y + 0.5*t + d

x4 = x
y4 = y + wlead*t

x5 = x
y5 = y

x6 = x 
y6 = y - wlead*t

x7 = x + xfrac*L - wmid*L
y7 = y - 0.5*t + d

x8 = x + xfrac*L
y8 = y - 0.5*t + d

x9 = x + xfrac*L + wmid*L
y9 = y - 0.5*t + d

x10 = x + L
y10 = y - 0.5*h


cp0 = np.array([(x0,y0),(x1,y1),(x2,y2)])
cp1 = np.array([(x2,y2),(x3,y3),(x4,y4),(x5,y5)])
cp2 = np.array([(x5,y5),(x6,y6),(x7,y7),(x8,y8)])
cp3 = np.array([(x8,y8),(x9,y9),(x10,y10)])

point_list = [cp0,cp1,cp2,cp3]

nt = 21
t = np.linspace(0,1,nt)
vertices = np.empty([0,2])
nspline = len(point_list)
for i in range(nspline):
    cp = point_list[i]
    nknot = cp.shape[0]
    if nknot==3:
        vert = np.array([quadratic_bezier(ti, np.array(cp[0]), np.array(cp[1]), np.array(cp[2])) for ti in t])
    elif nknot==4:
        vert = np.array([cubic_bezier(ti, np.array(cp[0]), np.array(cp[1]), np.array(cp[2]), np.array(cp[3])) for ti in t])

    if vertices.shape[0]==0:
        vertices = np.append(vertices,vert,axis=0)
    elif np.array_equal(vertices[-1,:], vert[0,:]):
        vertices = np.append(vertices,vert[1:,:],axis=0)
    else:
        vertices = np.append(vertices,vert,axis=0)

nV = vertices.shape[0]

elements = np.zeros([nV-1,2],int)
elements[:,0] = np.arange(0,nV-1,1)
elements[:,1] = np.arange(1,nV,1)
elements = np.append(elements,np.array([[nV-1,0]]),axis=0) # close profile

fig, ax = plt.subplots()
ax.plot(vertices[:,0],vertices[:,1])
plt.axis('equal')
plt.show()