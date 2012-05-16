from __future__ import division
import component, fuse_sections
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class halfbody(component.component):

    def __init__(self, nx, ny, nz):
        Ps = []
        Ks = []

        P, K = self.createSurfaces(Ks, ny, nz, -2, 3, 0)
        Ps.extend(P)
        Ks.append(K)
         
        P, K = self.createSurfaces(Ks, nz, nx, 3, 1, 1)
        Ps.extend(P) 
        Ks.append(K)
           
        P, K = self.createSurfaces(Ks, ny, nx, -2, 1, 1)
        Ps.extend(P)  
        Ks.append(K)  
      
        P, K = self.createSurfaces(Ks, nz, nx, -3, 1, 0)
        Ps.extend(P)   
        Ks.append(K)

        P, K = self.createSurfaces(Ks, ny, nz, -2, -3, 1)
        Ps.extend(P) 
        Ks.append(K) 

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.Ps = Ps   
        self.Ks = Ks  

        self.oml0 = [] 

    def setDOFs(self):
        oml0 = self.oml0
        for f in range(5):
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,:] = True
        for f in [0]:
            for j in [0]:
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,0] = False
                    edge = oml0.surf_edge[self.Ks[f][i,j],0,0]
                    edge = abs(edge) - 1
                    oml0.edge_c1[edge,:] = True
        for f in [1]:
            for j in range(self.Ks[f].shape[1]):
                for i in [0]:
                    oml0.surf_c1[self.Ks[f][i,j],0,:] = False
                    edge = oml0.surf_edge[self.Ks[f][i,j],1,0]
                    edge = abs(edge) - 1
                    oml0.edge_c1[edge,:] = True
        for f in [3]:
            for j in range(self.Ks[f].shape[1]):
                for i in [-1]:
                    oml0.surf_c1[self.Ks[f][i,j],-1,:] = False
                    edge = oml0.surf_edge[self.Ks[f][i,j],1,1]
                    edge = abs(edge) - 1
                    oml0.edge_c1[edge,:] = True
        for f in [4]:
            for j in [-1]:
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,-1] = False
                    edge = oml0.surf_edge[self.Ks[f][i,j],0,1]
                    edge = abs(edge) - 1
                    oml0.edge_c1[edge,:] = True

    def isExteriorDOF(self, f, uType, vType):
        value = False
        if f==0:
            if uType==2 and vType==0:
                value = True
        elif f==1:
            if uType==0 and vType==2:
                value = True
        elif f==3:
            if uType==-1 and vType==2:
                value = True
        elif f==4:
            if uType==2 and vType==-1:
                value = True
        return value

    def initializeParameters(self):
        Ns = self.Ns
        self.offset = numpy.zeros(3)
        self.L = 10
        self.rz = numpy.ones(Ns[1].shape[1])
        self.ry = numpy.ones(Ns[1].shape[1])
        self.y0 = numpy.ones(Ns[1].shape[1])
        self.x0 = numpy.ones(Ns[1].shape[1])
        self.Ls = numpy.ones(self.Ks[1].shape[1]+2)
        self.sections = []
        for i in range(Ns[1].shape[1]):
            self.sections.append(fuse_sections.circular)

    def setSections(self, section, shape):
        Ns = self.Ns
        for j in range(Ns[2].shape[1]):
            if Ns[2][1,j,3]==section:
                self.sections[j] = shape

    def propagateQs(self):
        coneL = 0.05
        coneR = 0.85
        Ns = self.Ns
        Qs = self.Qs
        for f in range(1,4):
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    Qs[f][i,j,:] = self.offset
                    Qs[f][i,j,0] += sum(self.Ls[:Ns[f][i,j,3]+1])
                    Qs[f][i,j,0] += self.Ls[Ns[f][i,j,3]+1]*Ns[f][i,j,4]/(self.getni(f,1)[Ns[f][i,j,3]])
                    z,y = self.sections[j](self.rz[j],self.ry[j],self.y0[j],f-1,i/(Ns[f].shape[0]-1))
                    Qs[f][i,j,1] += y
                    Qs[f][i,j,2] += z
        for f in [0,4]:
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    Qs[f][i,j,:] = self.offset
                    if f==0:
                        x,y,z = fuse_sections.cone(self.Ls[0],coneR*self.rz[0],coneR*self.ry[0],2,2,i/(Ns[f].shape[0]-1),j/(Ns[f].shape[1]-1))
                        Qs[f][i,j,1] += self.y0[0]
                    else:
                        x,y,z = fuse_sections.cone(-self.Ls[-1],coneR*self.rz[-1],coneR*self.ry[-1],2,2,i/(Ns[f].shape[0]-1),j/(Ns[f].shape[1]-1))
                        Qs[f][i,j,0] += sum(self.Ls[:-1])
                        Qs[f][i,j,1] += self.y0[-1]
                    Qs[f][i,j,:] += [x,y,z]

    def setMain(self, y0, rz, ry):
        self.y0[:] = y0
        self.rz[:] = rz
        self.ry[:] = ry
                    
    def setNose(self, n, y0, rz, ry):
        Qs = self.Qs
        for j in range(1,n):
            r = Qs[2][1,n-1,0] - Qs[2][1,j,0]
            r /= Qs[2][1,n-1,0] - Qs[2][1,1,0]
            self.y0[j] = y0[2] - (y0[2]-y0[1])*r**2
            self.rz[j] = rz[2] - (rz[2]-rz[1])*r**2
            self.ry[j] = ry[2] - (ry[2]-ry[1])*r**2
        self.y0[0] = y0[0]
        self.rz[0] = rz[0]
        self.ry[0] = ry[0]
                    
    def setTail(self, n, y0, rz, ry):
        Qs = self.Qs
        for j in range(-n,-1):
            r = Qs[2][1,j,0] - Qs[2][1,-n,0]
            r /= Qs[2][1,-2,0] - Qs[2][1,-n,0]
            self.y0[j] = y0[2] - (y0[2]-y0[1])*r
            self.rz[j] = rz[2] - (rz[2]-rz[1])*r
            self.ry[j] = ry[2] - (ry[2]-ry[1])*r
        self.y0[-1] = y0[0]
        self.rz[-1] = rz[0]
        self.ry[-1] = ry[0]
        


if __name__ == '__main__':  

    h = halfbody([2,3],[3,4],[3,3])
    P = h.Ps
    print h.Ks
    
    ax = p3.Axes3D(pylab.figure())
    for k in range(len(P)):
        ax.plot_wireframe(P[k][:,:,0],P[k][:,:,1],P[k][:,:,2])
    pylab.show()
