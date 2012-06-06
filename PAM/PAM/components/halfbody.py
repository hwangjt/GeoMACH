from __future__ import division
from PAM.components import component, Property, fuse_sections
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class halfbody(component):

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
        self.props = {
            'posx':Property(Ns[1].shape[1]),
            'posy':Property(Ns[1].shape[1]),
            'ry':Property(Ns[1].shape[1]),
            'rz':Property(Ns[1].shape[1]),
            'noseL':0.05,
            'tailL':0.02
            }
        self.sections = []
        for i in range(Ns[1].shape[1]):
            self.sections.append(fuse_sections.circular)

    def setSections(self, section, shape):
        Ns = self.Ns
        for j in range(Ns[2].shape[1]):
            if Ns[2][1,j,3]==section:
                self.sections[j] = shape

    def propagateQs(self):
        c = 0.8
        Ns = self.Ns
        Qs = self.Qs
        for f in range(1,4):
            for j in range(Ns[f].shape[1]):
                posx = self.props['posx'].data[j]
                posy = self.props['posy'].data[j]
                rz = self.props['rz'].data[j]
                ry = self.props['ry'].data[j]
                for i in range(Ns[f].shape[0]):
                    z,y = self.sections[j](rz,ry,f-1,i/(Ns[f].shape[0]-1))
                    Qs[f][i,j,:] = self.offset + [posx,posy,0] + [0,y,z]
        for f in [0,-1]:
            posx = self.props['posx'].data[f]
            posy = self.props['posy'].data[1+3*f]
            rz = self.props['rz'].data[f]
            ry = self.props['ry'].data[f]
            if f==0:
                self.props['noseL'] = 0.5*(rz+ry)*(2**0.5-1)
                L = self.props['noseL']
            else:
                self.props['tailL'] = 0.5*(rz+ry)*(2**0.5-1)
                L = -self.props['tailL']
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    x,y,z = fuse_sections.cone(L,c*rz,c*ry,1,1,i/(Ns[f].shape[0]-1),j/(Ns[f].shape[1]-1))
                    Qs[f][i,j,:] = self.offset + [posx,posy,0] + [x,y,z] 
                    if f==-1:
                        Qs[f][i,j,0] += self.props['noseL'] + self.props['tailL']
        


if __name__ == '__main__':  

    h = halfbody([2,3],[3,4],[3,3])
    P = h.Ps
    print h.Ks
    
    ax = p3.Axes3D(pylab.figure())
    for k in range(len(P)):
        ax.plot_wireframe(P[k][:,:,0],P[k][:,:,1],P[k][:,:,2])
    pylab.show()
