"""
GeoMACH tip class
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


class PGMtip(PGMinterpolant):
    """ Tip component """

    def __init__(self, config, comp, side, weight=0.1):
        """
        Parameters
        ----------
        config : ``PGMconfiguration``
           Pointer to the configuration containing this component
        comp : ``str``
           Name of the wing component
        side : ``str``
           The side of the wing being closed off by this tip:
           'right' for the v=0 side,
           'left' for the v=1 side
        weight : ``float``
           The weight given to the tangent vectors
           of the tip of the wing when interpolating
        """
        super(PGMtip, self).__init__()

        self._comp = config.comps[comp]
        self._side = side
        self._weight = weight

        self._num_surf_wing = self._comp.faces['upp']._num_surf['u']

        self.faces[''] = PGMface(2, self._num_surf_wing)

    def set_diff(self):
        face = self.faces['']
        face.set_diff_surf(True)

        if self._side == 'right':
            face.set_diff_surf(False, ind_j=0, ind_v=0)
            face.set_diff_edge(True, 'v0', ind_j=0)
        elif self._side == 'left':
            face.set_diff_surf(False, ind_j=-1, ind_v=2)
            face.set_diff_edge(True, 'v1', ind_j=-1)

    def compute(self, name):
        upp = self._comp.faces['upp']
        low = self._comp.faces['low']
        face = self.faces['']

        num_u = face._num_cp_total['u']
        num_v = face._num_cp_total['v']

        if self._side == 'right':
            N = upp.vec_inds['cp_prim'][:,:2,:]
            S = low.vec_inds['cp_prim'][::-1,:2,:]
        elif self._side == 'left':
            N = upp.vec_inds['cp_prim'][::-1,-1:-3:-1,:]
            S = low.vec_inds['cp_prim'][:,-1:-3:-1,:]

        nD = 3 * 2 * num_v + 3 * 4 * (num_u-2) * num_v
        Da, Di, Dj = PGMlib.computetip(nD, num_u, num_v, 
                                       self._weight, N, S, 
                                       face.vec_inds['cp_bez'])

        Das, Dis, Djs = [Da], [Di], [Dj]
        if name == 'cp_bez':
            return Das, Dis, Djs
        elif name == 'cp_coons':
            return Das, Dis, Djs
        elif name == 'cp_prim':
            return [], [], []
