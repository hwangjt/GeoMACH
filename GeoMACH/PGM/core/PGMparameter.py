"""
GeoMACH parameter class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.core.PGMobject import PGMobject


class PGMparameter(PGMobject):

    def __init__(self, num_u, num_v,
                 pos_u=None, pos_v=None, 
                 order_u=None, order_v=None):
        """
        Parameters
        ----------
        num_u, num_v : ``int``
           Number of parameters
           (B-spline control points)
        pos_u, pos_v : ``list`` of num_u, num_v ``float``s
           Parametric positions of the control points
           in the u and v directions
           normalized to the interval [0,1]
        order_u, order_v : ``int``
           Order of the B-spline in each direction
        """
        super(PGMparameter, self).__init__()

        self._pos = {'u': numpy.array(pos_u)
                     if pos_u is not None
                     else numpy.linspace(0, 1, num_u),
                     'v': numpy.array(pos_v)
                     if pos_v is not None
                     else numpy.linspace(0, 1, num_v)}
        self._order = {'u': order_u if order_u is not None 
                       else min(2, num_u),
                       'v': order_v if order_v is not None
                       else min(2, num_v)}
        self._prop = None
        self._num_cp = {'u': num_u,
                        'v': num_v}
        self._num_pt = {'u': 0, 'v': 0}
        self._shape = (num_u, num_v)

    def assemble_sizes(self, prop, num_u, num_v):
        self._prop = prop
        self._num_pt['u'] = num_u
        self._num_pt['v'] = num_v

    def val(self, val):
        val = numpy.array(val).reshape(self._shape, 
                                       order='F')
        self.vec_data['param'][:] = val

    def compute(self, name):
        order = self._order
        num_cp = self._num_cp
        num_pt = self._num_pt
        pos = self._pos
        cp_indices = self.vec_inds['param']
        pt_indices = self._prop.vec_inds['prop']

        nnz = num_pt['u'] * num_pt['v'] * order['u'] * order['v']
        vals, rows, cols \
            = PGMlib.computebspline(nnz, 
                                    order['u'], order['v'],
                                    num_cp['u'], num_cp['v'],
                                    num_pt['u'], num_pt['v'],
                                    cp_indices, pt_indices,
                                    pos['u'], pos['v'])
        return [vals], [rows], [cols]
