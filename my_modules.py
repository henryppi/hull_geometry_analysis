import numpy as np
import vtk
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import leastsq
from scipy.optimize import least_squares
import bezier


def vecMag(v):
    if len(v.shape)==1:
        return np.sqrt(np.sum(v**2))
    else:
        return np.sqrt(np.sum(v**2,axis=1))
    
def get_bbox(vert):
    return np.min(vert[:,0]),np.max(vert[:,0]),np.min(vert[:,1]),np.max(vert[:,1])

    
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

def nonlinear_least_squares_scipy_bounds(ferr,v0,bounds,data):
    # print(v0)
    # test = ferr(v0,data)
    # print(test)
    tmp = []
    # fun(x, *args, **kwargs)
    sol = least_squares(ferr,v0, args=(tmp,data),bounds=bounds)
    return sol.x

def computeDotVectorized(a,b):
#     https://stackoverflow.com/questions/15616742/vectorized-way-of-calculating-row-wise-dot-product-two-matrices-with-scipy
    return np.einsum('ij,ij->i', a, b)

def computePointLineDistance(points,xa,ya,xb,yb):
    #https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    nP = points.shape[0]
    pa = np.array([[xa,ya]])
    pb = np.array([[xb,yb]])
    len2 = (xb-xa)**2 + (yb-ya)**2
    t = np.max([np.zeros(nP,float),np.min([np.ones(nP,float),computeDotVectorized(points-pa,pb - pa)/len2],axis=0)],axis=0)
    p = np.concatenate([np.array([points[:,0]-(xa + t*(xb-xa))]).T,np.array([points[:,1]-(ya + t*(yb-ya))]).T],axis=1)
    return np.linalg.norm(p,axis=1)

def fit_bezier_rudder(points,para_init):
        
    def para_to_nodes(para):
        # para = [x,y,L,t,h,xfrac,wlead,wfwd,wrev,wtail]
        #        [0.33, 0.05, 0.25, 0.2, 0.2, 0.2]
        x = para[0]
        y = para[1]
        L = para[2]
        t = para[3]
        h = para[4]
        xfrac = para[5]
        wlead = para[6]
        wfwd = para[7]
        wrev = para[8]
        wtail = para[9]

        x0 = 0.0 + L 
        y0 = 0.0 + h

        x1 = 0.0 + (1-wtail)*L
        y1 = 0.0 + h

        x2 = 0.0 + xfrac*L+wrev*L 
        y2 = 0.0 + 0.5*t

        x3 = 0.0 + xfrac*L
        y3 = 0.0 + 0.5*t

        x4 = x3-wfwd*L
        y4 = 0.0 + 0.5*t

        x5 = 0.0
        y5 = 0.0 + wlead*t

        x6 = 0.0
        y6 = 0.0
        # 
        x7 = x + x6
        y7 = y -y5

        x8 = x + x4
        y8 = y -y4

        x9 = x + x3
        y9 = y -y3

        x10 = x + x2
        y10 = y -y2

        x11 = x + x1
        y11 = y -y1

        x12 = x + x0
        y12 = y -y0

        x0+=x
        x1+=x
        x2+=x
        x3+=x
        x4+=x
        x5+=x
        x6+=x

        nodesA = np.asfortranarray([\
                [x0, x1, x2, x3],\
                [y0, y1, y2, y3],\
                ])
        nodesB = np.asfortranarray([\
                [x3, x4, x5, x6],\
                [y3, y4, y5, y6],\
                ])
        nodesC = np.asfortranarray([\
                [x6, x7, x8, x9],\
                [y6, y7, y8, y9],\
                ])
        nodesD = np.asfortranarray([\
                [x9, x10, x11, x12],\
                [y9, y10, y11, y12],\
                ])

        return [nodesA,nodesB,nodesC,nodesD]

    def gen_bezier_4point_segments(point_list,nseg=25,tail_open=True):
        nspline = len(point_list)
        
        s_vals = np.linspace(0.0, 1.0, nseg)
        vertices = np.empty([0,2])
        for i in range(nspline):
            seg_nodes = point_list[i]
            seg_curve = bezier.Curve(seg_nodes, degree=3)
            if i==nspline-1:
                seg_vertices = (seg_curve.evaluate_multi(s_vals)).T
            else:
                seg_vertices = (seg_curve.evaluate_multi(s_vals[:-1])).T
            vertices = np.append(vertices,seg_vertices,axis=0)
        
        nV = vertices.shape[0]
        
        # vertices[:,0] += -np.min(vertices[:,0])
        
        elements = np.zeros([nV-1,2],int)
        elements[:,0] = np.arange(0,nV-1,1)
        elements[:,1] = np.arange(1,nV,1)

        if not tail_open==True:
            # print(vertices.shape)
            # print(np.array([seg_nodes[:,3]]))
            vertices = np.append(vertices,np.array([seg_nodes[:,3]]),axis=0)
            elements = np.append(elements,np.array([[nV,nV+1]]),axis=0)

        return vertices,elements 
    
    def rudder_fun(para,nSeg=25):
        node_list = para_to_nodes(para)
        tail_open = True
        vert,elem = gen_bezier_4point_segments(node_list,nseg=nSeg,tail_open=tail_open)
        return vert,elem

    def error_function(para,tmp,data):
        # print('errf')
        # print(np.around(para,2).tolist())
        LARGE = 1e10

        points = data.reshape([-1,2])
        n_points = points.shape[0]

        vertices,elements = rudder_fun(para)

        # fig, ax = plt.subplots(figsize=(8,8),frameon=False)
        # ax.plot(points[:,0],points[:,1],'.r')
        # ax.plot(vertices[:,0],vertices[:,1],'+-k',lw=2)
        # ax.axis('equal')
        # plt.show()

        # vertices[:,0] += para[0]
        # vertices[:,1] += para[1]
        nE = elements.shape[0]
        # print(nE)
        # print(points.shape)
        pDist = LARGE*np.ones(n_points)
        for iE in range(nE):
            # print(iE)
            # print(elements[iE,:])
            tmpDist = computePointLineDistance(points, \
                                            vertices[elements[iE,0],0], \
                                            vertices[elements[iE,0],1], \
                                            vertices[elements[iE,1],0], \
                                            vertices[elements[iE,1],1])
            pDist = np.minimum(pDist,tmpDist)
        return pDist

    # def fit(points,rudder_fun,para_init):
    #     point_list = self.para_to_nodes(self.init_para)
    #     vert,elem = gen_bezier_4point_segments(point_list,open=True)
    #     data = points.ravel()
    #     error_function_bezier(para,tmp,data)
    #     return para

    data = points.ravel()
    # points[:,0] += -para[0]
    # points[:,1] += -para[1]

    # para = [x,y,L,t,h,xfrac,wlead,wfwd,wrev,wtail]
    #        [0,1,2,3,4,  5  ,  6  , 7  ,  8 ,  9  ] 
    b_low = len(para_init)*[-np.inf]
    b_high = len(para_init)*[np.inf]
    b_low[3] = 0
    b_low[4] = 0
    b_low[5] = 0
    b_low[6] = 0
    b_low[7] = 0
    b_low[8] = 0
    b_low[9] = 0
    b_high[2] = para_init[2]*1.0
    b_high[3] = para_init[3]*1.0
    b_high[4] = para_init[3]*1.0
    bounds = (b_low,b_high)

    para_fit = nonlinear_least_squares_scipy_bounds(error_function,para_init,bounds,data)
    # para_fine = [cx+para5[0],\
    #             cy+para5[1],\
    #             aoa]
    # para_fine.extend(para5[2:])
    # para = fit(points,rudder_fun,para_init)

    # bbox = get_bbox(points)
    
    # vert,elem = rudder_fun(para_init)
    # return vert,elem
    vert_fit,elem_fit = rudder_fun(para_fit)
    return para_fit,vert_fit,elem_fit

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