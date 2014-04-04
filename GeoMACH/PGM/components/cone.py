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
        vtx = lambda f, i, j: self.comp.Ps[self.comp.faces[f].surf_indices[i,j]][i,j,:]
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
        face = self.faces[0]
        face.surf_indices[:,:] = -1
        for j in range(nj):
            for i in range(ni):
                self.Ps.append(PGMlib.bilinearinterp(self.nP, ni, nj, i+1, j+1, verts))
                face.surf_indices[i,j] = j*ni + i

    def setDOFs(self):
        faces = self.faces
        faces[0].setC1('surf', val=True)
        if self.comp.bottom==0:
            faces[0].setC1('surf', i=-1, u=-1, val=False)
            faces[0].setC1('edge', i=-1, u=-1, val=True)

    def computeQs(self):
        getEdge = self.getEdge
        comp_faces = self.comp.faces
        zeros = numpy.zeros((1,3),order='F')
        if self.face==0:
            W = getEdge(comp_faces[0].cp_array, j=1, d=-1)
            N = getEdge(comp_faces[1].cp_array, j=1, d=1)
            E = getEdge(comp_faces[2].cp_array, j=1, d=1)
            if self.comp.bottom==2:
                S = getEdge(comp_faces[3].cp_array, j=1, d=-1)
            else:
                S = getEdge(comp_faces[1].cp_array, j=1, d=1)
        else:
            E = getEdge(comp_faces[0].cp_array, j=-2, d=-1)
            N = getEdge(comp_faces[1].cp_array, j=-2, d=-1)
            W = getEdge(comp_faces[2].cp_array, j=-2, d=1)
            if self.comp.bottom==2:
                S = getEdge(comp_faces[3].cp_array, j=-2, d=1)
            else:
                S = getEdge(comp_faces[1].cp_array, j=-2, d=-1)

        mu, mv = self.faces[0].num_cp_list[:]
        nu = range(3)
        nv = range(3)
        for k in range(3):
            nu[k] = sum(mu[self.si[k]:self.si[k+1]])
            nv[k] = sum(mv[self.sj[k]:self.sj[k+1]])

        v = self.variables
        self.faces[0].cp_array[:,:,:] = PGMlib.computecone(sum(nu)+1, sum(nv)+1, nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], v['scl']*self.comp.faces[0].num_cp[1], v['fC1'], v['mC1'], W, E, N, S, v['shp'])
