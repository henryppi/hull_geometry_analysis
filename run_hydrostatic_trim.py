import time
import matplotlib.pyplot as plt
import numpy as np
from my_modules import get_outline_xz_plane, plot3D_line, plot_line_xz
from my_modules import TriSurface
from my_modules import RigidBody


fn_hull = './CAD/20m_Craft.stl'

rho_water = 997.561 #[kg/m^3] density  (H2O starccm+ )
rho_air = 1.18415 # #[kg/m^3] density  (AIR starccm+ )
grav = np.array([0,0,-9.81],float)

tic = time.perf_counter()
outline_vert,outline_elem = get_outline_xz_plane(fn_hull)
toc = time.perf_counter()
print('runtime = {:.2f}'.format(toc-tic))

tri = TriSurface()
tri.readSTL(fn_hull)
tri.scaleSize(0.001)
tri.update()
tri.print_size()

L,B,T = tri.get_LBT()
mass = 40*1e3
cog = np.array([0.36*L,0.5*B,0.5])
tri.set_center(cog)
tri.print_size()

tri.setTrianglesToInitial()
cog0 = np.zeros(3)

rb = RigidBody()
rb.setTriangleBody(tri)
rb.setProperties(grav,rho_water)
rb.setMass(mass,cog0,L,B,T)
rb.initMatrizes()


dt = 0.01
t = 0.0
tEnd = 1000*dt

k = 0
kSave = 50

x = np.array([[0.0,0.0,-0.5]],float)#+np.array([cog])
axis = np.array([[1.0,0.0,0]],float)
# axis = vecNorm(axis)
angle = 0.0*np.pi/180.0
velo = np.array([[0,0,0]],float)
omega = np.array([[0,0,0]],float)
rb.setState(x,axis,angle,velo,omega)

rb.info()

tic = time.perf_counter()
X = []
PHI = []
T = []
damp=0.98#9#6
while t<tEnd:
    rb.computeForce()
    rb.integrate(dt,damp)
    x = rb.getX()
    phi = rb.getEulerAngles()
    T.append(t)
    X.append([x[0],x[1],x[2]])
    PHI.append([phi[0],phi[1],phi[2]])
    if np.mod(k,kSave)==0:
        print('t = %.2f' %(t),' heave (z) = %.2f [m]' %(x[2]),\
            '\t pitch = %.2f [deg]' %(phi[1]*180/np.pi))
    k+=1
    t+=dt
toc = time.perf_counter()
print('runtime = {:.2f}'.format(toc-tic))

T = np.array(T)
X = np.array(X)
PHI = np.array(PHI)*180.0/np.pi

if 1:
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot(T,X[:,2],'-b')
    ax2.plot(T,PHI[:,1],'-r')

    ax1.set_xlabel('time [s]')
    ax1.set_ylabel('heave [m]',color='b')
    ax2.set_ylabel('pitch [deg]',color='r')
    ax1.set_xlim([T[0],T[-1]])
    fig.tight_layout() 
    plt.show()

if 1:
    fn_trim = fn_hull[:-4]+'_trim'
    rb.tri.reset()
    # rb.tri.scaleSize(1000)
    rb.writeWholeSTL(fn_trim)
    outline_vert_trim,outline_elem_trim = get_outline_xz_plane(fn_trim+'.stl')
    outline_vert_trim*=1000
    fig2 = plt.figure()
    # ax2 = fig2.add_subplot(111, projection='3d') 
    ax2 = fig2.add_subplot(111) 
    color = 'r'
    linewidth=2
    linestyle='--'
    # plot3D_line(ax2,outline_vert,outline_elem,color,linewidth,linestyle)
    plot_line_xz(ax2,outline_vert,outline_elem,color,linewidth,linestyle)
    plot_line_xz(ax2,outline_vert_trim,outline_elem_trim,'k',linewidth,'-')
    ax2.axis('equal')
    plt.show()