from __future__ import division
import numpy, time

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Interpolant



class Cone(Interpolant):

    def __init__(self, comp, face):
        super(Cone,self).__init__() 

        self.comp = comp
        self.face = face
        self.initializeIndices()
        self.initializeFaces()
        self.initializeSurfaces()

    def initializeIndices(self):
        ny = self.comp.ms[1].shape[0]
        nz = self.comp.ms[2].shape[0]
        if self.comp.bottom==2:
            self.ni = [int(ny/2),0,int(ny/2)]
            self.nj = [int(nz/2),0,int(nz/2)]
        else:
            self.ni = [ny,0,0]
            self.nj = [int(nz/2),0,int(nz/2)]

    def initializeSurfaces(self):
        vtx = lambda f, i, j: self.comp.Ps[self.comp.Ks[f][i,j]][i,j,:]
        verts = numpy.zeros((2,2,3),order='F')
        if self.face==0:
            verts[0,0,:] = vtx(0, -1, 0)
            verts[1,0,:] = vtx(0, 0, 0)
            verts[0,1,:] = vtx(2, 0, 0)
            verts[1,1,:] = vtx(2, -1, 0)
        else:
            verts[0,0,:] = vtx(2, 0, -1)
            verts[1,0,:] = vtx(2, -1, -1)
            verts[0,1,:] = vtx(0, -1, -1)
            verts[1,1,:] = vtx(0, 0, -1)
        ni = sum(self.ni)
        nj = sum(self.nj)

        self.Ps = []
        self.Ks = [-numpy.ones((ni,nj),int)]
        for j in range(nj):
            for i in range(ni):
                self.Ps.append(PGMlib.bilinearinterp(self.nP, ni, nj, i+1, j+1, verts))
                self.Ks[0][i,j] = j*ni + i

    def setDOFs(self):
        self.setC1('surf', 0, val=True)
        if self.comp.bottom==0:
            self.setC1('surf', 0, i=-1, u=-1, val=False)
            self.setC1('edge', 0, i=-1, u=-1, val=True)

    def computeQs(self):
        getEdge = self.getEdge
        Qs = self.comp.Qs
        zeros = numpy.zeros((1,3),order='F')
        if self.face==0:
            W = getEdge(Qs[0], j=1, d=-1)
            N = getEdge(Qs[1], j=1, d=1)
            E = getEdge(Qs[2], j=1, d=1)
            if self.comp.bottom==2:
                S = getEdge(Qs[3], j=1, d=-1)
            else:
                S = getEdge(Qs[1], j=1, d=1)
        else:
            E = getEdge(Qs[0], j=-2, d=-1)
            N = getEdge(Qs[1], j=-2, d=-1)
            W = getEdge(Qs[2], j=-2, d=1)
            if self.comp.bottom==2:
                S = getEdge(Qs[3], j=-2, d=1)
            else:
                S = getEdge(Qs[1], j=-2, d=-1)

        mu = self.getms(0,0)
        mv = self.getms(0,1)
        nu = range(3)
        nv = range(3)
        for k in range(3):
            nu[k] = sum(mu[self.si[k]:self.si[k+1]])
            nv[k] = sum(mv[self.sj[k]:self.sj[k+1]])

        v = self.variables
        self.Qs[0] = PGMlib.computecone(sum(nu)+1, sum(nv)+1, nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], v['scl']*self.comp.Qs[0].shape[1], v['fC1'], v['mC1'], W, E, N, S, v['shp'])
