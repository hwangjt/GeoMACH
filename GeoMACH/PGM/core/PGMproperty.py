"""
GeoMACH property class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.PGM.core.PGMobject import PGMobject


class PGMproperty(PGMobject):

    def __init__(self):
        super(PGMproperty, self).__init__()

        self._num_pt = {'u': 0, 'v': 0}
        self.params = OrderedDict()

    def assemble_sizes(self, num_u, num_v):
        self._num_pt['u'] = num_u
        self._num_pt['v'] = num_v

        for param in self.params.values():
            param.assemble_sizes(self, num_u, num_v)

        self._shape = (num_u, num_v)

    def add_param(self, name, num_cp, val=None, 
                  rng_u=None, rng_v=None, 
                  order_u=None, order_v=None):
        limits = {'u': rng_u if rng_u is not None else [0, 1],
                  'v': rng_v if rng_v is not None else [0, 1]}
        order = {'u': order_u if order_u is not None 
                 else min(4, shape[0]),
                 'v': order_v if order_v is not None
                 else min(4, shape[1])}
        num_cp = self._num_cp
        self.params[name] = PGMparameter(num_pt, num_cp, 
                                         limits, order, val)
