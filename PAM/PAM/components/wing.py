from __future__ import division
from PAM.components import Component, Property, airfoils
import numpy, pylab, time
import mpl_toolkits.mplot3d.axes3d as p3
import PAM.PAMlib as PAMlib


class Wing(Component):

    def __init__(self, nb, nc, half=False, opentip=False):
        if half:
            self.faces = numpy.zeros((1,2),int)
            self.faces[0,:] = [1,2]
        else:
            self.faces = numpy.zeros((2,2),int)
            self.faces[0,:] = [1,2]
            self.faces[1,:] = [-1,2]

        Ps = []
        Ks = []

        if not half:
            P, K = self.createSurfaces(Ks, nc[::-1], nb, -1, 3, 0)
            for k in range(len(P)):
                for v in range(P[k].shape[1]):
                    for u in range(P[k].shape[0]):
                        if (opentip or P[k][u,v,2]!=1) and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                            P[k][u,v,1] = 1
            Ps.extend(P)
            Ks.append(K)

        P, K = self.createSurfaces(Ks, nc, nb, 1, 3, 0)
        for k in range(len(P)):
            for v in range(P[k].shape[1]):
                for u in range(P[k].shape[0]):
                    if (opentip or P[k][u,v,2]!=1) and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                        P[k][u,v,1] = -1
        Ps.extend(P)
        Ks.append(K)

        self.nb = nb
        self.nc = nc
        self.Ps = Ps
        self.Ks = Ks
        self.half = half
        self.opentip = opentip

        self.oml0 = []

    def setDOFs(self):        
        half = self.half
        opentip = self.opentip
        nf = len(self.Ks)
        for f in range(nf):
            self.setC1('surf', f)
            self.setC1('surf', f, j=0, v=0, val=False)
            self.setC1('surf', f, i=-f, u=-f, val=False)
            if opentip or half:
                self.setC1('surf', f, j=-1, v=-1, val=False)
            if half:
                self.setC1('surf', f, i=f-1, u=f-1, val=False)
        for f in range(nf):
            self.setC1('edge', f, j=0, v=0)
            self.setC1('edge', f, i=-f, u=-f)
            if opentip or half:
                self.setC1('edge', f, j=-1, v=-1)
            if half:
                self.setC1('edge', f, i=f-1, u=f-1)
        for f in range(nf):
            self.setCornerC1(f, i=-f, j=0, val=False)
            if opentip or half:
                self.setCornerC1(f, i=-f, j=-1, val=False)
            if half:
                self.setCornerC1(f, i=f-1, j=0, val=False)
                self.setCornerC1(f, i=f-1, j=-1, val=False)

    def isExteriorDOF(self, f, uType, vType, i, j):
        check = self.check
        half = self.half
        opentip = self.opentip
        value = check(uType,vType,v=0) or check(uType,vType,u=-f) or check(uType,vType,u=-f,v=0)
        if opentip or half:
            value = value or check(uType,vType,v=-1) or check(uType,vType,u=-f,v=-1)
        if half:
            value = value or check(uType,vType,u=f-1) or check(uType,vType,u=f-1,v=0) or check(uType,vType,u=f-1,v=-1)
        return value

    def initializeParameters(self):
        Ns = self.Ns
        self.offset = numpy.zeros(3)
        self.SECTshape = numpy.zeros((len(self.Ks),Ns[0].shape[0],Ns[0].shape[1],3),order='F')
        self.SECTrot0 = numpy.zeros((Ns[0].shape[1],3),order='F')
        self.props = {
            'chord':Property(Ns[0].shape[1]),
            'posx':Property(Ns[0].shape[1]),
            'posy':Property(Ns[0].shape[1]),
            'posz':Property(Ns[0].shape[1]),
            'rotx':Property(Ns[0].shape[1]),
            'roty':Property(Ns[0].shape[1]),
            'rotz':Property(Ns[0].shape[1]),
            'prpx':Property(Ns[0].shape[1]),
            'prpy':Property(Ns[0].shape[1])
            }
        self.setAirfoil("naca0012.dat")

    def setAirfoil(self,filename):
        airfoil = airfoils.getAirfoil(filename)
        if self.half:
            airfoil[0][:,:] = airfoil[1][:,:]
        Ps = airfoils.fitAirfoil(self,airfoil,rev=self.half)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.SECTshape[f,:,j,:2] = Ps[f][:,:]
        if self.half:
            self.SECTshape[0,:,-1,1] = 0
        
    def propagateQs(self):
        a = 0.25
        b = 0.0
        Ns = self.Ns
        Qs = self.Qs
        for f in range(len(self.Ks)):
            p = self.props
            Qs[f][:,:,:] = PAMlib.computewingsections(Ns[f].shape[0], Ns[f].shape[1], a, b, p['posx'].data, p['posy'].data, p['posz'].data, p['rotx'].data, p['roty'].data, p['rotz'].data, p['prpx'].data, p['prpy'].data, p['chord'].data, self.SECTshape[f])
            Qs[f][:,:,0] += self.offset[0]
            Qs[f][:,:,1] += self.offset[1]
            Qs[f][:,:,2] += self.offset[2]
        if self.half:
            Qs[0][:,-1,2] = 0
            Qs[0][0,:,2] = 0
            Qs[0][-1,:,2] = 0

    def getFlattenedC(self, f, ii, jj):
        if f==0:
            return [jj,1 - ii,0]
        else:
            return [jj,ii,0]

    def getAR(self):
        return 5

    def getSkinIndices(self):
        return [[0],[1]]
