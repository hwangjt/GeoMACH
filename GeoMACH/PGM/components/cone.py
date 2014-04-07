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
            verts[0,0,:] = vtx('rgt', -1, 0)
            verts[1,0,:] = vtx('rgt', 0, 0)
            verts[0,1,:] = vtx('lft', 0, 0)
            verts[1,1,:] = vtx('lft', -1, 0)
        else:
            verts[0,0,:] = vtx('lft', 0, -1)
            verts[1,0,:] = vtx('lft', -1, -1)
            verts[0,1,:] = vtx('rgt', -1, -1)
            verts[1,1,:] = vtx('rgt', 0, -1)
        ni = sum(self.ni)
        nj = sum(self.nj)

        self.Ps = []
        face = self.faces['def']
        face.surf_indices[:,:] = -1
        for j in range(nj):
            for i in range(ni):
                self.Ps.append(PGMlib.bilinearinterp(self.nP, ni, nj, i+1, j+1, verts))
                face.surf_indices[i,j] = j*ni + i

    def setDOFs(self):
        face = self.faces['def']
        face.setC1('surf', val=True)
        if self.comp.bottom==0:
            face.setC1('surf', i=-1, u=-1, val=False)
            face.setC1('edge', i=-1, u=-1, val=True)

    def computeQs(self):
        getEdge = self.getEdge
        comp_faces = self.comp.faces
        face = self.faces['def']
        zeros = numpy.zeros((1,3),order='F')
        if self.face==0:
            W = getEdge(comp_faces['rgt'].cp_array, j=1, d=-1)
            N = getEdge(comp_faces['top'].cp_array, j=1, d=1)
            E = getEdge(comp_faces['lft'].cp_array, j=1, d=1)
            if self.comp.bottom==2:
                S = getEdge(comp_faces['bot'].cp_array, j=1, d=-1)
            else:
                S = getEdge(comp_faces['top'].cp_array, j=1, d=1)
        else:
            E = getEdge(comp_faces['rgt'].cp_array, j=-2, d=-1)
            N = getEdge(comp_faces['top'].cp_array, j=-2, d=-1)
            W = getEdge(comp_faces['lft'].cp_array, j=-2, d=1)
            if self.comp.bottom==2:
                S = getEdge(comp_faces['bot'].cp_array, j=-2, d=1)
            else:
                S = getEdge(comp_faces['top'].cp_array, j=-2, d=-1)

        mu, mv = face.num_cp_list[:]
        nu = range(3)
        nv = range(3)
        for k in range(3):
            nu[k] = sum(mu[self.si[k]:self.si[k+1]])
            nv[k] = sum(mv[self.sj[k]:self.sj[k+1]])

        p = self.properties
        face.cp_array[:,:,:] = PGMlib.computecone(sum(nu)+1, sum(nv)+1, nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], p['scl']*self.comp.faces['rgt'].num_cp[1], p['fC1'], p['mC1'], W, E, N, S, p['shp','def'])
