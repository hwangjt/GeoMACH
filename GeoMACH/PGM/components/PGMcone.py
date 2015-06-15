"""
GeoMACH cone class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components.PGMinterpolant import PGMinterpolant
from GeoMACH.PGM.core.PGMface import PGMface


class PGMcone(PGMinterpolant):
    """ Cone component """

    def __init__(self, config, comp, side, weight=1):
        """
        Parameters
        ----------
        config : ``PGMconfiguration``
           Pointer to the configuration containing this component
        comp : ``str``
           Name of the body component
        side : ``str``
           The side of the body being closed off:
           'front' for the v=0 side,
           'rear' for the v=1 side
        weight : ``float``
           The weight given to the tangent vectors
           of the tip of the body when interpolating
        """
        super(PGMcone, self).__init__()

        self._comp = config.comps[comp]
        self._side = side
        self._weight = weight

        ny = self._comp.faces['rgt']._num_surf['u']
        nz = self._comp.faces['top']._num_surf['u']

        self._num_surf['y'] = ny
        self._num_surf['z'] = nz

        self.faces[''] = PGMface(ny, nz)

    def set_diff(self):
        face = self.faces['']
        face.set_diff_surf(True)

    def compute(self, name):
        face = self.faces['']
        num_u = face._num_cp_total['u']
        num_v = face._num_cp_total['v']

        if name == 'cp_prim': #If we are at the cp_prim step...
            return super(PGMcone, self).compute(name) #Call the function that sets up the normal properties
        elif name == 'cp_bez': #If we are at the cp_bez step...
            rgt = self._comp.faces['rgt']
            top = self._comp.faces['top']
            lft = self._comp.faces['lft']
            bot = self._comp.faces['bot']
            if self._side == 'front':
                W = rgt.vec_inds['cp_prim'][::-1,:2,:]
                N = top.vec_inds['cp_prim'][:,:2,:]
                E = lft.vec_inds['cp_prim'][:,:2,:]
                S = bot.vec_inds['cp_prim'][::-1,:2,:]
            elif self._side == 'rear':
                E = rgt.vec_inds['cp_prim'][::-1,-1:-3:-1,:]
                N = top.vec_inds['cp_prim'][::-1,-1:-3:-1,:]
                W = lft.vec_inds['cp_prim'][:,-1:-3:-1,:]
                S = bot.vec_inds['cp_prim'][:,-1:-3:-1,:]
            nD = 0
            nD += 3 * 2 * num_v + 3 * 2 * (num_u - 2)
            nD += 3 * 8 * (num_v - 3 + num_u - 3)
            nD += 3 * 8

            Da, Di, Dj = PGMlib.computeconewireframe(nD, num_u, num_v, 1.0, self._weight, W, E, N, S, face.vec_inds['cp_bez'])
            Das, Dis, Djs = super(PGMcone, self).compute(name) #We will recover identity matrices just to carry over the normal parameters (Check PGMinterpolant.py)
            return Das + [Da], Dis + [Di], Djs + [Dj]
        elif name == 'cp_coons': #If we are at the cp_coons step...
            nD = 3 * 8 * 4 * ((num_u-1)/2 -1) * ((num_v-1)/2 -1)
            Da, Di, Dj = PGMlib.computeconecoons(nD, num_u, num_v, face.vec_inds['cp_coons'])
            Das, Dis, Djs = super(PGMcone, self).compute(name) #We will recover identity matrices just to carry over the normal parameters (Check PGMinterpolant.py)
            return Das + [Da], Dis + [Di], Djs + [Dj]
            
