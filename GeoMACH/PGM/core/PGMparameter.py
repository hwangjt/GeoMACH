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
                 limits_u=None, limits_v=None, 
                 order_u=None, order_v=None):
        """
        Parameters
        ----------
        num_u, num_v : ``int``
           Number of parameters
           (B-spline control points)
        limits_u, limits_v : ``list`` of 2 ``int``s
           The intervals for u and v of the properties
           controlled by this parameter.
           Normalized so limits_u, limits_v
           :math:`\subseteq [0, 1]`
        order_u, order_v : ``int``
           Order of the B-spline in each direction
        """
        super(PGMparameter, self).__init__()

        self._limits = {'u': limits_u 
                        if limits_u is not None else [0, 1],
                        'v': limits_v
                        if limits_v is not None else [0, 1]}
        self._order = {'u': order_u if order_u is not None 
                       else min(2, num_u),
                       'v': order_v if order_v is not None
                       else min(2, num_v)}
        self._num_cp = {'u': num_u,
                        'v': num_v}
        self._num_prop = {'u': 0, 'v': 0}
        self._prop = None
        self._shape = (num_u, num_v)

    def assemble_sizes(self, prop, num_u, num_v):
        self._prop = prop
        self._num_prop['u'] = num_u
        self._num_prop['v'] = num_v

    def val(self, val):
        val = numpy.array(val).reshape(self._shape, 
                                       order='F')
        self.vec_data['param'][:] = val

    def compute(self, name):
        tu1, tu2 = self._limits['u']
        tv1, tv2 = self._limits['v']
        num_u, num_v = self._num_prop['u'], self._num_prop['v']

        if num_u > 1:
            range_u = [ind_i for ind_i in xrange(num_u)
                       if tu1 <= ind_i / (num_u-1) <= tu2]
        else:
            range_u = [0]
        if num_v > 1:
            range_v = [ind_j for ind_j in xrange(num_v)
                       if tv1 <= ind_j / (num_v-1) <= tv2]
        else:
            range_v = [0]
        order_u = self._order['u']
        order_v = self._order['v']
        num_cp_u = self._num_cp['u']
        num_cp_v = self._num_cp['v']
        num_pt_u = len(range_u)
        num_pt_v = len(range_v)
        ind_u1, ind_u2 = range_u[0], range_u[-1] + 1
        ind_v1, ind_v2 = range_v[0], range_v[-1] + 1

        cp_indices = self.vec_inds['param']
        pt_indices = self._prop.vec_inds['prop']
        pt_indices = pt_indices[ind_u1: ind_u2, ind_v1: ind_v2]

        nnz = num_pt_u * num_pt_v * order_u * order_v
        vals, rows, cols \
            = PGMlib.computebspline(nnz, order_u, order_v,
                                    num_cp_u, num_cp_v,
                                    num_pt_u, num_pt_v,
                                    cp_indices, pt_indices)
        return [vals], [rows], [cols]
