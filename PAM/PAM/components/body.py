from __future__ import division
from PAM.components import Component, Property
import numpy, pylab, time
import PAM.PAMlib as PAMlib
import mpl_toolkits.mplot3d.axes3d as p3


class Body(Component):

    def __init__(self, nx, ny, nz, full=False):
        if full:
            self.faces = numpy.zeros((6,2),int)
        else:
            self.faces = numpy.zeros((5,2),int)
        self.faces[0,:] = [-2,3]
        self.faces[1,:] = [3,1]
        self.faces[2,:] = [-2,1]
        self.faces[3,:] = [-3,1]
        self.faces[4,:] = [-2,-3]
        if full:
            self.faces[5,:] = [2,1]

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
      
        P, K = self.createSurfaces(Ks, nz[::-1], nx, -3, 1, 0)
        Ps.extend(P)   
        Ks.append(K)

        P, K = self.createSurfaces(Ks, ny, nz[::-1], -2, -3, 1)
        Ps.extend(P) 
        Ks.append(K) 

        if full:
            P, K = self.createSurfaces(Ks, ny, nx, 2, 1, 0)
            Ps.extend(P) 
            Ks.append(K)             

        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.Ps = Ps   
        self.Ks = Ks  
        self.full = full

        self.oml0 = [] 

    def setDOFs(self):
        for f in range(len(self.Ks)):
            self.setC1('surf', f)
        if not self.full:
            self.setC1('surf', 0, j=0, v=0, val=False)
            self.setC1('surf', 1, i=0, u=0, val=False)
            self.setC1('surf', 3, i=-1, u=-1, val=False)
            self.setC1('surf', 4, j=-1, v=-1, val=False)
            self.setC1('edge', 0, j=0, v=0)
            self.setC1('edge', 1, i=0, u=0)
            self.setC1('edge', 3, i=-1, u=-1)
            self.setC1('edge', 4, j=-1, v=-1)

    def isExteriorDOF(self, f, uType, vType, i, j):
        check = self.check
        value = False
        if self.full:
            value = False
        elif f==0:
            value = check(uType,vType,v=0)
        elif f==1:
            value = check(uType,vType,u=0)
        elif f==3:
            value = check(uType,vType,u=-1)
        elif f==4:
            value = check(uType,vType,v=-1)
        return value

    def initializeParameters(self):
        Ns = self.Ns
        self.offset = numpy.zeros(3)
        self.props = {
            'posx':Property(Ns[1].shape[1]),
            'posy':Property(Ns[1].shape[1]),
            'ry':Property(Ns[1].shape[1]),
            'rz':Property(Ns[1].shape[1]),
            't1U':Property(Ns[1].shape[1]),
            't2U':Property(Ns[1].shape[1]),
            't1L':Property(Ns[1].shape[1]),
            't2L':Property(Ns[1].shape[1]),
            'noseL':0.1,
            'tailL':0.1
            }
        self.setSections()

    def setSections(self, sections=[], t1U=0, t2U=1, t1L=0, t2L=1):
        Ns = self.Ns
        for j in range(Ns[2].shape[1]):
            for i in range(Ns[2].shape[0]):
                val = Ns[2][i,j,3]
                if not val == -1:
                    break
            found = False
            for k in range(len(sections)):
                found = found or (val==sections[k])
            if found or sections==[]:
                self.props['t1U'].data[j] = t1U
                self.props['t2U'].data[j] = t2U
                self.props['t1L'].data[j] = t1L
                self.props['t2L'].data[j] = t2L

    def propagateQs(self):
        Ns = self.Ns
        Qs = self.Qs
        for f in range(len(Ns)):
            if f==0 or f==4:
                if f==0:
                    L = self.props['noseL']
                    i0 = 0
                    i1 = 1
                    i2 = 2
                else:
                    L = self.props['tailL']
                    i0 = -1
                    i1 = -2
                    i2 = -3
                p = self.props
                dx = abs(self.props['posx'].data[i2] - self.props['posx'].data[i1])  
                Qs[f][:,:,:] = PAMlib.computecone(self.full, Ns[f].shape[0], Ns[f].shape[1], L, p['posy'].data[i0], p['posy'].data[i1], p['posy'].data[i2], p['ry'].data[i1], p['ry'].data[i2], p['rz'].data[i1], p['rz'].data[i2], dx)
                if f==4:
                    Qs[f][:,:,0] *= -1
                    Qs[f][:,:,0] += self.props['noseL']
                    Qs[f][:,:,0] += self.props['posx'].data[-2] - self.props['posx'].data[1]
                    Qs[f][:,:,0] += self.props['tailL']
                    Qs[f][:,:,:] = Qs[f][:,::-1,:]
                Qs[f][:,:,0] += self.offset[0]
                Qs[f][:,:,1] += self.offset[1]
                Qs[f][:,:,2] += self.offset[2]
            else:
                if f==1 and self.full:
                    t1 = 3/4.0
                    t2 = 1/4.0
                elif f==1 and not self.full:
                    t1 = 1/2.0
                    t2 = 1/4.0
                elif f==2:
                    t1 = 1/4.0
                    t2 = -1/4.0
                elif f==3 and self.full:
                    t1 = 7/4.0
                    t2 = 5/4.0
                elif f==3 and not self.full:
                    t1 = 7/4.0
                    t2 = 6/4.0
                elif f==5:
                    t1 = 5/4.0
                    t2 = 3/4.0
                p = self.props
                Qs[f][:,:,:] = PAMlib.computebodysections(Ns[f].shape[0], Ns[f].shape[1], t1, t2, p['noseL'], p['posx'].data, p['posy'].data, p['rz'].data, p['ry'].data, p['t1U'].data, p['t2U'].data, p['t1L'].data, p['t2L'].data)
                Qs[f][:,:,0] += self.offset[0]
                Qs[f][:,:,1] += self.offset[1]
                Qs[f][:,:,2] += self.offset[2]

    def getFlattenedC(self, f, i, j, ni, nj):
        ii = i/(ni-1)
        jj = j/(nj-1)
        if f==1:
            return [jj,1 - 0.25*ii,0]
        elif f==2:
            return [jj,0.75 - 0.5*ii,0]
        elif f==3:
            return [jj,0.25 - 0.25*ii,0]

    def getAR(self):
        return 8

    def getSkinIndices(self):
        return [[1,2,3]]
