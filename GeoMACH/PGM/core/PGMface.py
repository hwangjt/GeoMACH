"""
GeoMACH face class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
from collections import OrderedDict

from GeoMACH.PGM.core.PGMobject import PGMobject
from GeoMACH.PGM.core.PGMsurf import PGMsurf


class PGMface(PGMobject):

    def __init__(self, num_i, num_j):
        """
        Attributes
        ----------
        **surfs** : ``OrderedDict``
           Dictionary of ``PGMsurf`` objects
        **_name** : ``str``
           Face name
        **_num** : ``int``
           Unique face index within the containing component
        **_bse** : ``BSEmodel``
           BSE model of the configuration
        **_num_surf** : dictionary of ``int``s
           Number of surfaces in the 'u' and 'v' directions
        **_num_cp_list** : dictionary of ``numpy.ndarray``, ``int``[num]
           Vector of number of control points minus one
           for each surface in the 'u' and 'v' directions
        **_num_pt_list** : dictionary of ``numpy.ndarray``, ``int``[num]
           Vector of number of points minus one
           for each surface in the 'u' and 'v' directions
        **_num_cp_list** : dictionary of ``int``s
           Total number of control points
           in the 'u' and 'v' directions
        **_surf_indices** : ``numpy.ndarray``, ``int``[ni, nj]
           Array of global surface indices at the configuration level
        """
        super(PGMface, self).__init__()

        self.surfs = OrderedDict()

        self._name = None
        self._num = None
        self._bse = None

        ones = numpy.ones
        zeros = numpy.zeros
        self._num_surf = {'u': num_i, 
                          'v': num_j}
        self._num_cp_list = {'u': 3*ones(num_i, int),
                             'v': 3*ones(num_j,int)}
        self._num_pt_list = {'u': zeros(num_i, int),
                             'v': zeros(num_j,int)}
        self._num_cp_total = {'u': 0,
                              'v': 0}
        self._surf_indices = numpy.zeros((num_i, num_j), int, order='F')
        for ind_j in xrange(num_j):
            for ind_i in xrange(num_i):
                self.surfs[ind_i, ind_j] = PGMsurf()

    def _set_bspline_option(self, option, axis, vals):
        """
        Arguments
        ---------
        **option** : ``str``
           'order', 'num_cp', or 'num_pt'
        **axis** : ``str``
           'u' or 'v'
        **vals** : ``list`` of ``int``s
           Length must be number of surfaces
           or 1, to copy to all
        """
        num_surf = self._num_surf
        vals = vals * num_surf[axis] if len(vals) == 1 else vals

        for ind_j in xrange(num_surf['v']):
            for ind_i in xrange(num_surf['u']):
                isurf = self._surf_indices[ind_i, ind_j]
                if isurf != -1:
                    if axis == 'u':
                        val = vals[ind_i]
                    elif axis == 'v':
                        val = vals[ind_j]
                    self._bse.set_bspline_option(option, isurf, 
                                                 axis, val)

    def _get_surf_indices(self, ind_i, ind_j):
        if ind_i is None:
            ilist = range(self._num_surf['u'])
        else:
            ilist = [ind_i]

        if ind_j is None:
            jlist = range(self._num_surf['v'])
        else:
            jlist = [ind_j]

        surf_list = [self._surf_indices[ind_i, ind_j]
                     for ind_i in ilist
                     for ind_j in jlist]
        return [isurf for isurf in surf_list if isurf != -1]

    def set_info(self, name, num, bse):
        self._name = name
        self._num = num
        self._bse = bse

    def assemble_sizes(self, bse):
        num_cp_list = self._num_cp_list
        num_pt_list = self._num_pt_list
        num_cp_total = self._num_cp_total
        num_surf = self._num_surf
        surf_inds = self._surf_indices

        for ind_j in xrange(num_surf['v']):
            for ind_i in xrange(num_surf['u']):
                isurf = surf_inds[ind_i, ind_j]
                if isurf != -1:
                    for dim in ['u', 'v']:
                        if dim == 'u':
                            ind = ind_i
                        elif dim == 'v':
                            ind = ind_j
                        if bse is not None:
                            num_cp_list[dim][ind] \
                                = bse.get_bspline_option('num_cp', isurf, dim) - 1
                            num_pt_list[dim][ind] \
                                = bse.get_bspline_option('num_pt', isurf, dim) - 1
                        else:
                            num_cp_list[dim][ind] = 3
                            num_pt_list[dim][ind] = 9

                    num_u = num_cp_list['u'][ind_i]
                    num_v = num_cp_list['v'][ind_j]
                    surf = self.surfs[ind_i, ind_j] 
                    surf.set_shape(num_u, num_v)

        for dim in ['u', 'v']:
            num_cp_total[dim] = sum(num_cp_list[dim]) + 1

        num_u = self._num_cp_total['u']
        num_v = self._num_cp_total['v']
        self._shape = (num_u, num_v, 3)

    def initialize_bse(self, surfs_list, surf_index):
        num_surf = self._num_surf
        surf_inds = self._surf_indices

        for ind_j in xrange(num_surf['v']):
            for ind_i in xrange(num_surf['u']):
                if surf_inds[ind_i, ind_j] != -1:
                    surf = self.surfs[ind_i, ind_j]
                    data = surf.vec_data['df_surf']
                    surfs_list.append(data)

                    surf_inds[ind_i, ind_j] = surf_index
                    surf_index += 1
        return surf_index

    def set_option(self, option, dim, vals, both=True):
        """
        Parameters
        ----------
        option : ``str``
           'order', 'num_cp', or 'num_pt'
        dim : ``str``
           'u' or 'v'
        vals : ``list`` of ``int``(s)
           Length 1 to copy to all or
           length number of surfaces in the given direction
        both : ``bool``
           If option is 'num_cp' and both is ``True``,
           3 * vals is assigned to 'num_pt'.
           If option is 'num_pt' and both is ``True``,
           ``max``(4, ``int``(val/3)) is assigned to 'num_pt'
        """
        set_option = self._set_bspline_option

        set_option(option, dim, vals)
        if option == 'num_cp' and both:
            set_option('num_pt', dim, [3*val for val in vals])
        if option == 'num_pt' and both:
            set_option('num_cp', dim, [max(4, int(val/3))
                                       for val in vals])

    def set_diff_surf(self, val, ind_i=None, ind_j=None, 
                      ind_u=None, ind_v=None):
        for isurf in self._get_surf_indices(ind_i, ind_j):
            if isurf != -1:
                self._bse.set_diff_surf(val, isurf, ind_u, ind_v)

    def set_diff_edge(self, val, edge, 
                      ind_i=None, ind_j=None, side=None):
        for isurf in self._get_surf_indices(ind_i, ind_j):
            if isurf != -1:
                self._bse.set_diff_edge(val, isurf, edge, side)

    def set_diff_corner(self, val, ind_i, ind_j):
        if ind_i == 0 and ind_j == 0:
            self.set_diff_edge(val, 'u0', ind_i, ind_j, side=0)
            self.set_diff_edge(val, 'v0', ind_i, ind_j, side=0)
        elif ind_i == 0 and ind_j == -1:
            self.set_diff_edge(val, 'u0', ind_i, ind_j, side=1)
            self.set_diff_edge(val, 'v1', ind_i, ind_j, side=0)
        elif ind_i == -1 and ind_j == 0:
            self.set_diff_edge(val, 'u1', ind_i, ind_j, side=0)
            self.set_diff_edge(val, 'v0', ind_i, ind_j, side=1)
        elif ind_i == -1 and ind_j == -1:
            self.set_diff_edge(val, 'u1', ind_i, ind_j, side=1)
            self.set_diff_edge(val, 'v1', ind_i, ind_j, side=1)

    def compute(self, name):
        vals_list, rows_list, cols_list = [], [], []

        v_start, v_end = 0, 0
        for ind_j in xrange(self._num_surf['v']):
            num_v = self._num_cp_list['v'][ind_j]
            v_end += num_v

            u_start, u_end = 0, 0
            for ind_i in xrange(self._num_surf['u']):
                num_u = self._num_cp_list['u'][ind_i]
                u_end += num_u

                if self._surf_indices[ind_i, ind_j] != -1:
                    surf = self.surfs[ind_i, ind_j]
                    inds = self.vec_inds['cp_coons']
                    vals = numpy.ones(3 * (num_u+1) * (num_v+1))
                    rows = surf.vec_inds[name][:, :, :]
                    cols = inds[u_start : u_end + 1, 
                                v_start : v_end + 1, :]

                    vals_list.append(vals)
                    rows_list.append(rows.flatten())
                    cols_list.append(cols.flatten())

                u_start += num_u
            v_start += num_v

        return vals_list, rows_list, cols_list   
