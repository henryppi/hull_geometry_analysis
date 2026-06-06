import numpy as np
import vtk
import matplotlib.pyplot as plt

def vecMag(v):
    if len(v.shape)==1:
        return np.sqrt(np.sum(v**2))
    else:
        return np.sqrt(np.sum(v**2,axis=1))
    
def get4x4TransformationRotMat(R,t):
    M = np.matrix(np.zeros([4,4],float))
    M[0:3,0:3] = R
    M[0:3,3] = t
    M[3,3] = 1.0
    return M

def get_outline_xz_plane(fn_stl):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(fn_stl) # Replace with your STL file
    reader.Update()
    polydata = reader.GetOutput()

    points = polydata.GetPoints()
    npoints = points.GetNumberOfPoints()
    cells = polydata.GetPolys()
    ncells = cells.GetNumberOfCells()
    print('n points = ',npoints)
    print('n cells = ',ncells)

    plane_cut = vtk.vtkPlane()
    plane_cut.SetOrigin(0.0, 0.0, 0.0)  
    plane_cut.SetNormal(0.0, 1.0, 0.0)  

    cutter = vtk.vtkCutter()
    cutter.SetInputData(polydata)
    cutter.SetCutFunction(plane_cut)
    cutter.Update()

    stripper = vtk.vtkStripper()
    stripper.SetInputConnection(cutter.GetOutputPort())
    stripper.Update()
    cut_polylines = stripper.GetOutput()

    points_cut =  cut_polylines.GetPoints()
    npoints_cut = points_cut.GetNumberOfPoints()
    print(npoints_cut)

    points_vert = np.zeros([npoints_cut,3],float)
    for i in range(npoints_cut):
        point_coords = points_cut.GetPoint(i)
        points_vert[i,:] = point_coords

    points_elem = []

    vtk_lines = cut_polylines.GetLines()
    vtk_lines.InitTraversal()
    id_list = vtk.vtkIdList()
    while vtk_lines.GetNextCell(id_list):
        print("New Polyline Segment:")
        for i in range(id_list.GetNumberOfIds() - 1):
            point1_id = id_list.GetId(i)
            point2_id = id_list.GetId(i + 1)
            points_elem.append([id_list.GetId(i),id_list.GetId(i + 1)])

    points_elem = np.array(points_elem)
    return points_vert,points_elem

def plot3D_line(ax,vert,elem,color,linewidth,linestyle):
    import mpl_toolkits.mplot3d.axes3d as axes3d
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    lines = []
    for i in range(elem.shape[0]):
        e = elem[i]
        line = [vert[e[0],:],vert[e[1],:]]
        lines.append(line)
    lc = Line3DCollection(lines, colors = color,linewidth=linewidth,linestyle=linestyle)
    ax.add_collection3d(lc)

def plot_line_xz(ax,vert,elem,color,linewidth,linestyle):
    for i in range(elem.shape[0]):
        ax.plot(vert[elem[i],0],vert[elem[i],2],color=color,
                                    linewidth=linewidth,
                                    linestyle=linestyle)

class TriSurface:
    def __init__(self):
        self.vertices = np.zeros([0,3])
        self.elements = np.zeros([0,3],int)
        self.normals =  np.zeros([0,3])
        self.nElem = 0
        self.nVert = 0
        self.bbox = {'xmin':0,'xmax':0,'ymin':0,'ymax':0,'zmin':0,'zmax':0}

    def print_size(self):
        print("xmin = {:.2f},\tymin = {:.2f},\tzmin = {:.2f}"\
              .format(self.bbox['xmin'],self.bbox['ymin'],self.bbox['zmin']))
        print("xmax = {:.2f},\tymax = {:.2f},\tzmax = {:.2f}"\
              .format(self.bbox['xmax'],self.bbox['ymax'],self.bbox['zmax']))
        print("cx = {:.2f},\tcy = {:.2f},\tcz = {:.2f}"\
              .format(self.cx,self.cy,self.cz))
        print("L = {:.2f},\tB = {:.2f},\tT = {:.2f}"\
              .format(self.L,self.B,self.T))

    def get_LBT(self):
        return self.L,self.B,self.T

    def get_center(self):
        return self.cx, self.cy, self.cz

    def set_center(self,cog):
        self.update()
        vec = np.array([self.bbox['xmin'],self.bbox['ymin'],self.bbox['zmin']])
        self.translate(-vec)
        self.translate(-cog)
        self.update()

    def getBoundingBox(self):
        bbox = {}
        bbox['xmin'] = self.vertices[:,0].min()
        bbox['xmax'] = self.vertices[:,0].max()
        bbox['ymin'] = self.vertices[:,1].min()
        bbox['ymax'] = self.vertices[:,1].max()
        bbox['zmin'] = self.vertices[:,2].min()
        bbox['zmax'] = self.vertices[:,2].max()
        self.bbox = bbox
        return bbox
    
    def update(self):
        self.nVert,nDim = self.vertices.shape
        self.nElem,nEdge = self.elements.shape
        self.bbox = self.getBoundingBox()
        self.L = self.bbox['xmax']-self.bbox['xmin']
        self.B = self.bbox['ymax']-self.bbox['ymin']
        self.T = self.bbox['zmax']-self.bbox['zmin']
        self.cx = 0.5*(self.bbox['xmax']+self.bbox['xmin'])
        self.cy = 0.5*(self.bbox['ymax']+self.bbox['ymin'])
        self.cz = 0.5*(self.bbox['zmax']+self.bbox['zmin'])
        self.computeEdgeCrossProd()
        self.computeFaceDirMag()
        self.computeAreas()
        self.computeNormals()
    
    def set(self,vertices,elements):
        self.vertices = vertices
        self.elements = elements
        self.update()

    def add_triangles(self,vert,elem):
        nvert = self.vertices.shape[0]
        self.vertices = np.concatenate([self.vertices,vert])
        self.elements = np.concatenate([self.elements,elem+nvert])

    def reset(self):
        self.elements = np.copy(self.elements0)
        self.vertices = np.copy(self.vertices0)
        self.nVert,nDim = self.vertices.shape
        self.nElem,nEdge = self.elements.shape
        self.bbox = self.getBoundingBox()

    def computeEdgeCrossProd(self):
        self.e01 = self.vertices[self.elements[:,1],:]-self.vertices[self.elements[:,0],:]
        self.e02 = self.vertices[self.elements[:,2],:]-self.vertices[self.elements[:,0],:]
        self.crossProd0 = np.cross(self.e02,self.e01,axis=1)
        
    def computeFaceDirMag(self):
        self.dirMag = np.array([vecMag(self.crossProd0)]).T
        
    def computeAreas(self):
        self.areas = 0.5*self.dirMag
       
    def computeNormals(self):
        self.normals = self.crossProd0/self.dirMag
       
    def computePointPressure(self,rho,grav):
        self.pres = np.array([rho*grav*self.vertices[:,2]]).T
        
    def computeForceMoment(self,cog):
        f0 = (1.0/3.0)*self.areas*self.pres[self.elements[:,0]]*self.normals
        f1 = (1.0/3.0)*self.areas*self.pres[self.elements[:,1]]*self.normals
        f2 = (1.0/3.0)*self.areas*self.pres[self.elements[:,2]]*self.normals
        f = f0+f1+f2
        fSum = np.sum(f,axis=0)
        
        r0 = 0.25*(2.0*self.vertices[self.elements[:,0],:]+self.vertices[self.elements[:,1],:]+self.vertices[self.elements[:,2],:])-cog
        r1 = 0.25*(self.vertices[self.elements[:,0],:]+2.0*self.vertices[self.elements[:,1],:]+self.vertices[self.elements[:,2],:])-cog
        r2 = 0.25*(self.vertices[self.elements[:,0],:]+self.vertices[self.elements[:,1],:]+2.0*self.vertices[self.elements[:,2],:])-cog
        
        m0 = np.cross(r0,f0,axis=1)
        m1 = np.cross(r1,f1,axis=1)
        m2 = np.cross(r2,f2,axis=1)
        m = m0+m1+m2
        mSum = np.sum(m,axis=0)

        return fSum,mSum
    
    def getCenter(self):
        vec = np.zeros([1,3])
        vec[0,0] = np.mean([self.bbox['xmax'],self.bbox['xmin']])
        vec[0,1] = np.mean([self.bbox['ymax'],self.bbox['ymin']])
        vec[0,2] = np.mean([self.bbox['zmax'],self.bbox['zmin']])
        return vec
    
    def center(self):
        vec = self.getCenter()
        self.translate(-vec)
        self.update()
        
    def translate(self,vec):
        self.vertices = self.vertices+vec

    def move(self,mat4x4):
        n = self.vertices0.shape[0]
        vert0 = np.concatenate([self.vertices0,np.ones([n,1],float)],axis=1)
        vert = (mat4x4*vert0.T).T
        self.vertices[:,0:3] = vert[:,0:3]

    def getIndZeroAreaTriangles(self):
        n=self.elements.shape[0]
        r = np.arange(n)
        ind = r[np.where(self.areas[:,0]==0.0)[0]]

        return ind
    
    def setTrianglesToInitial(self):
        self.elements0 = np.copy(self.elements)
        self.vertices0 = np.copy(self.vertices)
        
    def removeElements(self,ind):
        ind = np.unique(ind)
        self.elements = np.delete(self.elements,ind,0)
        
    def cropPlaneZ(self):
        indBelow, indIntersect, indAbove = self.getCutIndices()
        indUnderCut, indAboveCut = self.cutTriangles(indIntersect)
        indDelete = np.concatenate([indIntersect,indAbove,indAboveCut])
        
        self.removeElements(indDelete)
        self.update()
        
        indZero = self.getIndZeroAreaTriangles()
        self.removeElements(indZero)
        self.update()

    def cutTriangles(self,ind):
        n0 = self.vertices.shape[0]
        n = ind.shape[0]
        nElem = self.elements.shape[0]
        vecNew = np.zeros([3*n,3],float)
        elemNew = np.zeros([4*n,3],int)
        indUnder = np.zeros([4*n],bool)
        ind4 = np.arange(4)
        for i in range(n):
            i0 = self.elements[ind[i],0]
            i1 = self.elements[ind[i],1]
            i2 = self.elements[ind[i],2]
            
            v0 = self.vertices[i0,:]
            v1 = self.vertices[i1,:]
            v2 = self.vertices[i2,:]
            e01 = v1-v0
            e12 = v2-v1
            e20 = v0-v2

            d0 = v0[2]
            d1 = v1[2]
            d2 = v2[2]
            
            side = np.zeros([3])
            if d0<0.0:
                side[0] = True
            else:
                side[0] = False
                
            if d1<0.0:
                side[1] = True
            else:
                side[1] = False
                            
            if d2<0.0:
                side[2] = True
            else:
                side[2] = False

            if side[0]==side[1]:
                v01 = 0.5*(v0+v1)
                elemSide = np.array([True, True, False, True])
                if side[2]:
                    elemSide = np.invert(elemSide)
            else:
                v01 = (np.abs(d0)/(np.abs(d0)+np.abs(d1)))*e01+v0

            if side[1]==side[2]:
                v12 = 0.5*(v1+v2)
                elemSide = np.array([False, True, True, True])
                if side[0]:
                    elemSide = np.invert(elemSide)
            else:
                v12 = (np.abs(d1)/(np.abs(d1)+np.abs(d2)))*e12+v1

            if side[2]==side[0]:
                v20 = 0.5*(v2+v0)
                elemSide = np.array([True, False, True, True])
                if side[1]:
                    elemSide = np.invert(elemSide)
            else:
                v20 = (np.abs(d2)/(np.abs(d2)+np.abs(d0)))*e20+v2
            
            indUnder[i*4+ind4] = elemSide
            
            i01 = n0 + i*3+0
            i12 = n0 + i*3+1
            i20 = n0 + i*3+2

            vecNew[i*3+0,:] = v01
            vecNew[i*3+1,:] = v12
            vecNew[i*3+2,:] = v20
            
            elemNew[i*4+0,:] = [i0, i01, i20]
            elemNew[i*4+1,:] = [i01, i1, i12]
            elemNew[i*4+2,:] = [i20, i12, i2]
            elemNew[i*4+3,:] = [i01, i12, i20]
            
        indAbove = np.invert(indUnder)
        elemRange = np.arange(0,n*4)
        indUnder = elemRange[indUnder]+nElem
        indAbove = elemRange[indAbove]+nElem
        
        self.vertices = np.concatenate([self.vertices,vecNew])
        self.elements = np.concatenate([self.elements,elemNew])

        return indUnder, indAbove
        

    def getCutIndices(self):
        n = self.vertices.shape[0]
        indZ = np.zeros(n,bool)
        indZ[np.where(self.vertices[:,2]<0.0)] = True
        
        ne = self.elements.shape[0]
        elemZ = np.zeros([ne,3],bool)
        elemZ[:,0] = indZ[self.elements[:,0]]
        elemZ[:,1] = indZ[self.elements[:,1]]
        elemZ[:,2] = indZ[self.elements[:,2]]
        
        indAll = np.all(elemZ,axis=1)
        indAny = np.any(elemZ,axis=1)
        
        indIntersect = np.where(indAll!=indAny)
        indBelow = np.where(indAll==True)
        indAbove = np.where(indAny==False)
        
        return indBelow[0], indIntersect[0], indAbove[0]

    def computeEdgeCrossProd(self):
        self.e01 = self.vertices[self.elements[:,1],:]-self.vertices[self.elements[:,0],:]
        self.e02 = self.vertices[self.elements[:,2],:]-self.vertices[self.elements[:,0],:]
        self.crossProd0 = np.cross(self.e02,self.e01,axis=1)
        
    def computeFaceDirMag(self):
        self.dirMag = np.array([vecMag(self.crossProd0)]).T
        
    def computeAreas(self):
        self.areas = 0.5*self.dirMag
       
    def computeNormals(self):
        self.normals = self.crossProd0/self.dirMag

    def scaleSize(self,fac):
        self.vertices[:,0] *= fac
        self.vertices[:,1] *= fac
        self.vertices[:,2] *= fac

    def readSTL(self,filename):
        fileSTL = open(filename, "r")
        normals = []
        v1 = []
        v2 = []
        v3 = []
        v = 0

        for line in fileSTL.readlines():
            line = line.lstrip()
            if line.find("facet") == 0:
                v = 0
                normal = line.split()[-3:]

                normals.append([float(normal[0]),float(normal[1]),float(normal[2])])
                continue

            if line.find("vertex") == 0:
                v = v + 1

                if v == 1:
                    vertex = line.split()[-3:]

                    v1.append([float(vertex[0]),float(vertex[1]),float(vertex[2])])

                if v == 2:
                    vertex = line.split()[-3:]

                    v2.append([float(vertex[0]),float(vertex[1]),float(vertex[2])])

                if v == 3:
                    vertex = line.split()[-3:]

                    v3.append([float(vertex[0]),float(vertex[1]),float(vertex[2])])

        fileSTL.close()

        n = len(v1)
        
        elements = np.zeros([n,3],int)
        vertices = np.concatenate([np.array(v1),np.array(v2),np.array(v3)],axis=0)
        elements[:,0] = np.arange(0,n)
        elements[:,1] = np.arange(0,n)+n
        elements[:,2] = np.arange(0,n)+2*n

        self.vertices = vertices
        self.elements = elements
        self.update() 

    def writeSTL(self,filename):
        filenamestl = filename+'.stl'
        outfile = open(filenamestl, 'w') # open file for writing
        outfile.write('solid '+filename+'\n')
        for iFace in range(self.nElem):
            vertex = np.concatenate([[np.array(self.vertices[self.elements[iFace,0],:])],\
                                     [np.array(self.vertices[self.elements[iFace,1],:])],\
                                     [np.array(self.vertices[self.elements[iFace,2],:])]],axis=0)

            # HACK to handle weird 3dim "vertex" matrix
            if ((vertex.shape[0]==3)&(vertex.shape[1]!=3)):
                vertex = vertex[:,0,:]
            iNormal = self.normals[iFace,:]
            outfile.write(' facet normal %2.5e %2.5e %2.5e \n' %(iNormal[0],iNormal[1],iNormal[2]))
            outfile.write('  outer loop\n')
            for i in range(3):
                outfile.write('   vertex %2.5e %2.5e %2.5e \n' %(vertex[i,0],vertex[i,1],vertex[i,2]))
            outfile.write('  endloop\n')
            outfile.write(' endfacet\n')
        outfile.write('endsolid '+filename+'\n')

class RigidBody():
    def __init__(self):
        self.x = np.array([[0,0,0]],float).T
        self.u = np.array([[0,0,0]],float).T
        
        self.q = np.array([[1,0,0,0]],float).T        
        self.w = np.array([[0,0,0]],float).T
        
        self.mass = 1.0
        self.inertia = np.matrix(np.eye(3,3))
        self.force = np.array([[0,0,0]],float).T
        self.tau = np.array([[0,0,0]],float).T
        self.grav = np.array([0,0,-9.81],float)
        
        self.u_ = np.zeros([6,1],float)
        self.s_ = np.zeros([7,1],float)
        self.f_ = np.zeros([6,1],float)

        self.u1_ = np.zeros([6,1],float)
        self.s1_ = np.zeros([7,1],float)

        self.tri = []

    def info(self):
        bbox = self.tri.getBoundingBox()
        print('mass = ',self.mass)
        print('bbox \t[xmin,xmax]',bbox['xmax'],bbox['xmin'],'\n \t', \
                      '[ymin,ymax]',bbox['ymax'],bbox['ymin'],'\n \t', \
                      '[zmin,ymax]',bbox['zmax'],bbox['zmin'])
                      
        print('inertia =\t',self.inertia[0,:],'\n\t\t',self.inertia[1,:],'\n\t\t',self.inertia[2,:])
        print('cog = ',self.cog )
        print('state x',self.x)

    def initMatrizes(self):
        self.S_ = np.matrix(np.zeros([7,6],float))
        self.S_[0:3,0:3] = np.eye(3,3)
        self.S_[3:7,3:6] = self.Q(self.q)
        self.M_ = np.matrix(np.eye(6,6))
        self.M_[0:3,0:3] *= self.mass
        self.M_[3:6,3:6] = self.inertia
        self.Minv_ = self.M_.I

    def setProperties(self,grav,rho):
        self.grav = grav
        self.rho = rho
    
    def setMass(self,mass,cog,L,B,T):
        self.mass = mass
        self.cog = cog
        spread = 0.5
        massOffset = spread*np.array([[L,0,0],[-L,0,0],[0,B,0],[0,-B,0],[0,0,T],[0,0,-T]],float)
        self.massList = np.ones([6,1],float)*(self.mass/6.0)
        self.massPos = np.zeros([6,3],float)+massOffset
        self.inertia = self.getInertiaTensor(self.massList,self.massPos)
        
    def computeGravitationForceMoment(self,cg,grav):
        f = grav*np.sum(self.massList[:,0],axis=0)*np.array([[0,0,-1]],float)
        r = self.massPos-cg
        m = np.sum(np.cross(r,f,axis=1),axis=0)
        return f,m
        
    def shiftCoG(self,x):
        self.cog +=x
        self.tri.translate(-np.array([self.cog]))
        
    def setState(self,x,axis,angle,velo,omega):
        self.x[:,0] = x[:]
        self.q = self.quaternionAxisAngle(axis,angle)
        
        self.u[:,0] = velo[:]
        self.w[:,0] = omega[:]
        
        self.u_[0:3,0] = self.u[:,0]
        self.u_[3:6,0] = self.w[:,0]
        
        self.s_[0:3] = self.x
        self.s_[3:7] = self.q

    def getCopy(self):
        return copy.deepcopy(self)
    
    def setTriangleBody(self,tri):
        self.tri = tri
        self.tri.setTrianglesToInitial()
        
    def getX(self):
        return self.x[:,0]
    
    def getEulerAngles(self):
        return self.eulerAngles(self.q)[:,0]
        
    def quaternionAxisAngle(self,axis,angle):
        q = np.zeros([4,1])
        q[0,0] = np.cos(angle*0.5)
        q[1:4,0] = np.sin(angle*0.5)*axis[:]
        return q
        
    def Q(self,q):
        return 0.5*np.array([[-q[1,0],-q[2,0],-q[3,0]],\
                             [q[0,0],q[3,0],-q[2,0]],\
                             [-q[3,0],q[0,0],q[1,0]],\
                             [q[2,0],-q[1,0],q[0,0]]])
        
    def eulerAngles(self,q):
        phi = np.arctan2(2.0*(q[0,0]*q[1,0]+q[2,0]*q[3,0]), 1.0-2.0*(q[1,0]**2+q[2,0]**2))
        psi = np.arcsin( 2.0*(q[0,0]*q[2,0]+q[3,0]*q[1,0]))
        chi = np.arctan2(2.0*(q[0,0]*q[3,0]+q[1,0]*q[2,0]), 1.0-2.0*(q[2,0]**2+q[3,0]**2))
        return np.array([[phi,psi,chi]]).T
        
    def R(self,q):
        w = q[0,0]
        x = q[1,0]
        y = q[2,0]
        z = q[3,0]
        
        n = w**2 + x**2 + y**2 + z**2
        if n==0.0:
            s=0.0
        else:
            s=2.0/n
        wx = s * w * x
        wy = s * w * y
        wz = s * w * z
        xx = s * x * x
        xy = s * x * y
        xz = s * x * z
        yy = s * y * y
        yz = s * y * z
        zz = s * z * z

        R = np.matrix([[1.0 - (yy + zz),       xy - wz,             xz + wy],\
                       [      xy + wz, 1.0 - (xx + zz),             yz - wx],\
                       [      xz - wy,       yz + wx,       1.0 - (xx + yy)]])
        return R
    
    def getInertiaTensor(self,M,X):
        Ixx = np.sum(M[:,0]*(X[:,1]**2+X[:,2]**2))
        Iyy = np.sum(M[:,0]*(X[:,0]**2+X[:,2]**2))
        Izz = np.sum(M[:,0]*(X[:,0]**2+X[:,1]**2))
        Ixy = Iyx = - np.sum(M[:,0]*X[:,0]*X[:,1])
        Ixz = Izx = - np.sum(M[:,0]*X[:,0]*X[:,2])
        Iyz = Izy = - np.sum(M[:,0]*X[:,1]*X[:,2])
        return np.matrix([[Ixx,Ixy,Ixz],[Iyx,Iyy,Iyz],[Izx,Izy,Izz]],float)
    
    def writeWholeSTL(self,fn2):
        # self.tri.reset()
        R = self.R(self.q)
        t = self.x
        T = get4x4TransformationRotMat(R,t)

        self.tri.move(T)
        self.tri.update()
        self.tri.writeSTL(fn2)
    
    def computeForce(self):
        self.tri.reset()
        R = self.R(self.q)
        t = self.x
        T = get4x4TransformationRotMat(R,t)

        self.tri.move(T)
        self.tri.cropPlaneZ()
        self.tri.update()
        self.tri.computePointPressure(self.rho,self.grav[2])

        cg = self.x[:,0]
        force,moment = self.tri.computeForceMoment(cg)
        fg,mg = self.computeGravitationForceMoment(cg,-self.grav[2])
        self.F = force + fg
        self.tau[:,0] = moment[:]+mg[:]

    def integrate(self,dt,damp=0.999):
        self.S_[3:7,3:6] = self.Q(self.q)
        self.f_[0:3,0] = self.F[:]
        self.f_[3:6,0] = self.tau[:,0]

        self.u1_ = self.u_ + dt*self.Minv_*self.f_
        self.s1_ = self.s_ + dt*self.S_*self.u1_

        self.u1_[0:3,0] *=damp
        self.u1_[3:6,0] *=damp

        self.s_[:] =self.s1_[:,0]
        self.u_[:] =self.u1_[:,0]

        self.u[:] = self.u1_[0:3,0]
        self.w[:] = self.u1_[3:6,0]
        self.x[:] = self.s1_[0:3,0]
        self.q[:] = self.s1_[3:7,0]

    def setCurrentState(self):
        self.tri.reset()
        R = self.R(self.q)
        t = self.x
        T = get4x4TransformationRotMat(R,t)
        self.tri.move(T)