from __future__ import division
import numpy, time

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Interpolant, Property



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
        self.ni = [int(ny/2),0,int(ny/2)]
        self.nj = [int(nz/2),0,int(nz/2)]

    def setDOFs(self):
        face = self.faces['def']
        face.setC1('surf', val=True)

    def compute_cp_wireframe(self):
        comp_faces = self.comp.faces
        face = self.faces['def']
        
        if self.face == 0:
            W = comp_faces['rgt'].cp_indices[::-1,:2,:]
            N = comp_faces['top'].cp_indices[:,:2,:]
            E = comp_faces['lft'].cp_indices[:,:2,:]
            S = comp_faces['bot'].cp_indices[::-1,:2,:]
        else:
            E = comp_faces['rgt'].cp_indices[::-1,-1:-3:-1,:]
            N = comp_faces['top'].cp_indices[::-1,-1:-3:-1,:]
            W = comp_faces['lft'].cp_indices[:,-1:-3:-1,:]
            S = comp_faces['bot'].cp_indices[:,-1:-3:-1,:]

        nu, nv = face.num_cp[:]

        nD = 0
        nD += 3 * 2 * nv + 3 * 2 * (nu - 2)
        nD += 3 * 8 * (nv - 3 + nu - 3)
        nD += 3 * 8

        fC1 = self.props['fC1'].prop_vec
        mC1 = self.props['mC1'].prop_vec
        Da, Di, Dj = PGMlib.computeconewireframe(nD, nu, nv, 1.0, fC1, mC1, W, E, N, S, face.cp_indices)
        return Da, Di, Dj

    def compute_cp_surfs(self):
        nu, nv = self.faces['def'].num_cp[:]

        nD = 3 * 8 * 4 * ((nu-1)/2 -1) * ((nv-1)/2 -1)

        Da, Di, Dj = PGMlib.computeconecoons(nD, nu, nv, self.faces['def'].cp_indices)
        return Da, Di, Dj
