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

    def setDOFs(self):
        face = self.faces['def']
        face.setC1('surf', val=True)
        if self.comp.bottom==0:
            face.setC1('surf', i=-1, u=-1, val=False)
            face.setC1('edge', i=-1, u=-1, val=True)

    def compute_cp_wireframe(self):
        comp_faces = self.comp.faces
        face = self.faces['def']
        
        if self.face == 0:
            W = comp_faces['rgt'].cp_indices[::-1,:2]
            N = comp_faces['top'].cp_indices[:,:2]
            E = comp_faces['lft'].cp_indices[:,:2]
            S = comp_faces['bot'].cp_indices[::-1,:2]
        else:
            E = comp_faces['rgt'].cp_indices[::-1,-1:-3:-1]
            N = comp_faces['top'].cp_indices[::-1,-1:-3:-1]
            W = comp_faces['lft'].cp_indices[:,-1:-3:-1]
            S = comp_faces['bot'].cp_indices[:,-1:-3:-1]

        nu, nv = face.num_cp[:]

        nD = 0
        nD += 2 * nv + 2 * (nu - 2)
        nD += 8 * (nv - 3 + nu - 3)
        nD += 8

        p = self.properties
        Da0, Di0, Dj0 = PGMlib.computeconewireframe(nD, nu, nv, 1.0, p['fC1'], p['mC1'], W, E, N, S, face.cp_indices)
        Da, Di, Dj = numpy.zeros(3*nD), numpy.zeros(3*nD, int), numpy.zeros(3*nD, int)
        for coord in xrange(3):
            Da[coord::3] = Da0
            Di[coord::3] = 3*Di0 + coord
            Dj[coord::3] = 3*Dj0 + coord
        return Da, Di, Dj

    def compute_cp_surfs(self):
        nu, nv = self.faces['def'].num_cp[:]

        nD = 8 * 4 * ((nu-1)/2 -1) * ((nv-1)/2 -1)

        Da0, Di0, Dj0 = PGMlib.computeconecoons(nD, nu, nv, self.faces['def'].cp_indices)
        Da, Di, Dj = numpy.zeros(3*nD), numpy.zeros(3*nD, int), numpy.zeros(3*nD, int)
        for coord in xrange(3):
            Da[coord::3] = Da0
            Di[coord::3] = 3*Di0 + coord
            Dj[coord::3] = 3*Dj0 + coord
        return Da, Di, Dj
