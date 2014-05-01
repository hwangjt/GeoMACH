from __future__ import division
import numpy, time

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Interpolant, Property



class Tip(Interpolant):

    def __init__(self, comp, side):
        super(Tip, self).__init__()

        self.comp = comp
        self.side = side
        self.initializeIndices()
        self.initializeFaces()
        self.initializeSurfaces()

    def initializeIndices(self):
        ny = 2
        nx = self.comp.ms[0].shape[0]
        self.ni = [1, 0, 1]
        self.nj = [0, nx, 0]

    def setDOFs(self):
        face = self.faces['def']
        face.setC1('surf', val=True)
        if self.side == 1:
            face.setC1('surf', j=-1, v=-1, val=False)
            face.setC1('edge', j=-1, v=-1, val=True)
        elif self.side == 0:
            face.setC1('surf', j=0, v=0, val=False)
            face.setC1('edge', j=0, v=0, val=True)

    def compute_cp_wireframe(self):
        comp_faces = self.comp.faces
        face = self.faces['def']

        nu, nv = face.num_cp[:]

        if self.side == 0:
            N = comp_faces['upp'].cp_indices[:,:2,:]
            S = comp_faces['low'].cp_indices[::-1,:2,:]
        elif self.side == 1:
            N = comp_faces['upp'].cp_indices[::-1,-1:-3:-1,:]
            S = comp_faces['low'].cp_indices[:,-1:-3:-1,:]

        nD = 3 * 2 * nv + 3 * 4 * (nu-2) * nv

        fC1 = self.props['fC1'].prop_vec
        mC1 = self.props['mC1'].prop_vec
        Da, Di, Dj = PGMlib.computetip(nD, nu, nv, fC1, mC1, N, S, face.cp_indices)
        return Da, Di, Dj
