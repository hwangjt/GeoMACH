"""
B-spline surface-modeling engine (BSE)
"""

# pylint: disable=E1101

from __future__ import division
import numpy
import time
import scipy.sparse
from GeoMACH.BSE import BSElib
from BSEvec import BSEvecUns, BSEvecStr


class BSEmodel(object):
    """ A continuous union of untrimmed B-spline surfaces """

    def __init__(self, initial_surfaces):
        """ Creates a new surface model

        Arguments
        ---------
        initial_surfaces : [list(nsurf) of double(nu, nv, 3)]

        Attributes
        ----------
        _num['surf'] : [int] number of surfaces
        _num['edge'] : [int] number of unique edges
        _num['vert'] : [int] number of unique vertices
        _num['group'] : [int] number of groups of edges that share properties

        """

        nsurf = len(initial_surfaces)
        nvert, nedge, ngroup, \
            surf_ptrs, edge_ptrs, \
            surf_group, edge_group \
            = self._compute_topology(initial_surfaces)

        self._num = {
            'surf': len(initial_surfaces),
            'vert': nvert,
            'edge': nedge,
            'group': ngroup,
        }

        self._topo = {
            'surf_ptrs': surf_ptrs,
            'edge_ptrs': edge_ptrs,
            'surf_group': surf_group,
            'edge_group': edge_group,
        }

        self._mult = {
            'vert': numpy.zeros(nvert, int),
            'edge': numpy.zeros(nedge, int),
            'diff_vert': numpy.zeros(nvert, int),
            'diff_edge': numpy.zeros(nedge, int),
        }

        self._bspline = {
            'order': 4 * numpy.ones(ngroup, int),
            'num_cp': 4 * numpy.ones(ngroup, int),
            'num_pt': 10 * numpy.ones(ngroup, int),
        }

        self._surf_indices = {
            'df': numpy.zeros((nsurf, 2), int, 'F'),
            'cp': numpy.zeros((nsurf, 2), int, 'F'),
            'pt': numpy.zeros((nsurf, 2), int, 'F'),
        }

        self._edge_indices = {
            'df': numpy.zeros((nedge, 2), int, 'F'),
            'cp': numpy.zeros((nedge, 2), int, 'F'),
            'pt': numpy.zeros((nedge, 2), int, 'F'),
        }

        self._str_indices = {
            'cp': numpy.zeros((nsurf, 2), int, 'F'),
            'pt': numpy.zeros((nsurf, 2), int, 'F'),
        }

        self._vert_indices = numpy.zeros(nvert, int)

        self._size = {
            'df_str': 0,
            'df': 0,
            'cp': 0,
            'cp_str': 0,
            'pt_str': 0,
            'pt': 0,
        }

        self.diff = {
            'surf': numpy.zeros((nsurf, 3, 3), bool, 'F'),
            'edge': numpy.zeros((nedge, 2), bool, 'F'),
        }

        self.hidden = numpy.zeros(nsurf, bool)

        self.jac = {
            'd(df)/d(df_str)': None,
            'd(cp)/d(df)': None,
            'd(cp_str)/d(cp)': None,
            'd(pt_str)/d(cp_str)': None,
            'd(pt)/d(pt_str)': None,
        }

        self.vec = {
            'df_str': None,
            'df': None,
            'cp': None,
            'cp_str': None,
            'pt_str': None,
            'pt': None,
        }
            

    def _compute_topology(self, initial_surfaces):
        """ Load an initial set of surfaces and compute the topology

        Arguments
        ---------
        initial_surfaces : list(nsurf) of double(nu, nv, 3)

        Returns
        -------
        nvert, nedge, ngroup : see __init__
        surf_vert : [int(nsurf, 2, 2)] surface to vertex mapping, 1-based
           (isurf, i, j) : u=i, v=j corner
        surf_edge : [int(nsurf, 2, 2)] surface to edge mapping, 1-based
           (isurf, 1, 1) : v=0 edge
           (isurf, 1, 2) : v=1 edge
           (isurf, 2, 1) : u=0 edge
           (isurf, 2, 2) : u=1 edge
        surf_group : [int(nsurf, 2)] surface to group mapping, 1-based
           (isurf, d) : d=0 for v=0, v=1 edges; d=1 for u=0, u=1 edges
        """

        nsurf = len(initial_surfaces)
        surfaces = numpy.zeros((nsurf, 3, 3, 3), float, 'F')

        for isurf in xrange(nsurf):
            surface = initial_surfaces[isurf]
            num_u, num_v = surface.shape[:2]
            mid_u1 = int(numpy.floor((num_u - 1) / 2.0))
            mid_u2 = int(numpy.ceil((num_u - 1) / 2.0))
            mid_v1 = int(numpy.floor((num_v - 1) / 2.0))
            mid_v2 = int(numpy.ceil((num_v - 1) / 2.0))

            for ind_u in xrange(2):
                for ind_v in xrange(2):
                    surfaces[isurf, -ind_u, -ind_v] = surface[-ind_u, -ind_v]

            for ind_u in xrange(2):
                surfaces[isurf, -ind_u, 1] += 0.5 * surface[-ind_u, mid_v1] + \
                                              0.5 * surface[-ind_u, mid_v2]

            for ind_v in xrange(2):
                surfaces[isurf, 1, -ind_v] += 0.5 * surface[mid_u1, -ind_v] + \
                                              0.5 * surface[mid_u2, -ind_v]

        nvert, nedge, surf_ptrs \
            = BSElib.computesurfconnectivities(nsurf, 1e-16, 1e-10, surfaces)

        edge_ptrs \
            = BSElib.computeedgeconnectivities(nsurf, nedge, surf_ptrs)

        ngroup, surf_group, edge_group \
            = BSElib.computegroups(nsurf, nedge, surf_ptrs)

        topology = [nvert, nedge, ngroup, \
                    surf_ptrs, edge_ptrs, \
                    surf_group, edge_group, \
                ]

        return topology

    def _compute_indices(self):
        num = self._num
        topo = self._topo
        diff = self.diff
        bspline = self._bspline
        vert_indices = self._vert_indices
        edge_indices = self._edge_indices
        surf_indices = self._surf_indices
        str_indices = self._str_indices
        mult = self._mult
        size = self._size

        vert_indices[:] \
            = BSElib.computevertindices(num['surf'], num['edge'], num['vert'],
                                        topo['surf_ptrs'], topo['edge_ptrs'],
                                        diff['surf'], diff['edge'])

        edge_indices['df'], edge_indices['cp'], edge_indices['pt'], \
            = BSElib.computeedgeindices(num['surf'], num['edge'],
                                        num['group'], topo['surf_ptrs'],
                                        diff['surf'], topo['edge_group'],
                                        bspline['num_cp'], bspline['num_pt'])

        surf_indices['df'], surf_indices['cp'], surf_indices['pt'], \
            = BSElib.computesurfindices(num['surf'], num['group'],
                                        topo['surf_group'],
                                        bspline['num_cp'], bspline['num_pt'])

        str_indices['df'], str_indices['cp'], str_indices['pt'], \
            = BSElib.computestrindices(num['surf'], num['group'],
                                       topo['surf_group'],
                                       bspline['num_cp'], bspline['num_pt'])

        mult['vert'], mult['edge'], mult['diff_vert'], mult['diff_edge'] \
            = BSElib.computemults(num['surf'], num['edge'], num['vert'],
                                  topo['surf_ptrs'], topo['edge_ptrs'],
                                  diff['surf'], diff['edge'])

        nvert_df = numpy.max(vert_indices)

        edge_indices['df'][:] += nvert_df
        edge_indices['cp'][:] += num['vert']
        edge_indices['pt'][:] += num['vert']

        surf_indices['df'][:] += numpy.max(edge_indices['df'])
        surf_indices['cp'][:] += numpy.max(edge_indices['cp'])
        surf_indices['pt'][:] += numpy.max(edge_indices['pt'])

        size['df_str'] = str_indices['df'][-1, -1]
        size['df'] = surf_indices['df'][-1, -1]
        size['cp'] = surf_indices['cp'][-1, -1]
        size['cp_str'] = str_indices['cp'][-1, -1]
        size['pt_str'] = str_indices['pt'][-1, -1]
        size['pt'] = surf_indices['pt'][-1, -1]

    def _compute_jacobians(self, nbs_jac=1):
        num = self._num
        topo = self._topo
        diff = self.diff
        bspline = self._bspline
        mult = self._mult
        vert_indices = self._vert_indices
        edge_indices = self._edge_indices
        surf_indices = self._surf_indices
        str_indices = self._str_indices
        size = self._size
        jac = self.jac

        nnz = BSElib.computeennz(num['surf'], num['edge'],
                                 num['vert'], num['group'],
                                 topo['surf_ptrs'], topo['surf_group'],
                                 mult['diff_vert'], mult['diff_edge'],
                                 bspline['num_cp'])
        data, rows, cols \
            = BSElib.computeemtx(nnz, num['surf'], num['edge'],
                                 num['vert'], num['group'],
                                 topo['surf_ptrs'], topo['surf_group'],
                                 surf_indices['df'], edge_indices['df'],
                                 str_indices['df'], vert_indices,
                                 mult['vert'], mult['edge'],
                                 mult['diff_vert'], mult['diff_edge'],
                                 bspline['num_cp'])
        jac['d(df)/d(df_str)'] \
            = scipy.sparse.csr_matrix((data, (rows-1, cols-1)),
                                      shape=(size['df'],
                                             size['df_str']))
            
                                 
        nnz = BSElib.computednnz(num['surf'], num['edge'],
                                num['vert'], num['group'],
                                topo['surf_group'], topo['edge_group'],
                                diff['surf'], diff['edge'],
                                mult['diff_vert'], mult['diff_edge'],
                                bspline['num_cp'])
        data, rows, cols \
            = BSElib.computedmtx(nnz, num['surf'], num['edge'],
                                 num['vert'], num['group'],
                                 topo['surf_group'], topo['edge_group'],
                                 topo['surf_ptrs'], topo['edge_ptrs'],
                                 diff['surf'], diff['edge'],
                                 mult['diff_vert'], mult['diff_edge'],
                                 surf_indices['df'], surf_indices['cp'],
                                 edge_indices['df'], edge_indices['cp'],
                                 vert_indices, bspline['num_cp'])
        jac['d(cp)/d(df)'] \
            = scipy.sparse.csr_matrix((data, (rows-1, cols-1)),
                                      shape=(size['cp'],
                                             size['df']))

        nnz = size['cp_str']
        data, rows, cols \
            = BSElib.computecmtx(nnz, num['surf'], num['edge'], num['group'],
                                 topo['surf_ptrs'], topo['surf_group'],
                                 surf_indices['cp'], edge_indices['cp'],
                                 str_indices['cp'], bspline['num_cp'])
        jac['d(cp_str)/d(cp)'] \
            = scipy.sparse.csr_matrix((data, (rows-1, cols-1)),
                                      shape=(size['cp_str'],
                                             size['cp']))

        nnz = BSElib.computebnnz(num['surf'], num['group'],
                                topo['surf_group'],
                                bspline['order'], bspline['num_pt'])
        uder = [0, 1, 0, 2, 1, 0]
        vder = [0, 0, 1, 0, 1, 2]
        name = ['', '_du', '_dv', '_duu', '_duv', '_dvv']
        for k in xrange(nbs_jac):
            data, rows, cols \
                = BSElib.computebmtx(nnz, num['surf'], num['group'], 
                                     uder[k], vder[k],
                                     str_indices['cp'], str_indices['pt'],
                                     topo['surf_group'], bspline['order'],
                                     bspline['num_cp'], bspline['num_pt'])
            jac['d(pt_str)/d(cp_str)' + name[k]] \
                = scipy.sparse.csr_matrix((data, (rows-1, cols-1)),
                                          shape=(size['pt_str'],
                                                 size['cp_str']))

        nnz = size['pt_str']
        data, rows, cols \
            = BSElib.computepmtx(nnz, num['surf'], num['edge'],
                                 num['vert'], num['group'],
                                 topo['surf_ptrs'], topo['surf_group'],
                                 surf_indices['pt'], edge_indices['pt'],
                                 str_indices['pt'], mult['vert'],
                                 mult['edge'], bspline['num_pt'])
        jac['d(pt)/d(pt_str)'] \
            = scipy.sparse.csr_matrix((data, (rows-1, cols-1)),
                                      shape=(size['pt'],
                                             size['pt_str']))

    def set_bspline_option(self, option, isurf, dim, val):
        if dim == 'u':
            idim = 0
        elif dim == 'v':
            idim = 1

        igroup = self._topo['surf_group'][isurf, idim] - 1
        self._bspline[option][igroup] = val

    def get_bspline_option(self, option, isurf, dim):
        if dim == 'u':
            idim = 0
        elif dim == 'v':
            idim = 1

        igroup = self._topo['surf_group'][isurf, idim] - 1
        return self._bspline[option][igroup]

    def set_diff_surf(self, val, isurf, ind_u=None, ind_v=None):
        if ind_u is None:
            ind_u1, ind_u2 = 0, 3
        else:
            ind_u1, ind_u2 = ind_u, ind_u + 1

        if ind_v is None:
            ind_v1, ind_v2 = 0, 3
        else:
            ind_v1, ind_v2 = ind_v, ind_v + 1

        self.diff['surf'][isurf, ind_u1:ind_u2, ind_v1:ind_v2] = val

    def set_diff_edge(self, val, isurf, edge, side=None):
        diff = self.diff
        surf_ptrs = self._topo['surf_ptrs']

        if edge == 'u0':
            iedge = surf_ptrs[isurf, 0, 1]
        elif edge == 'u1':
            iedge = surf_ptrs[isurf, 2, 1]
        elif edge == 'v0':
            iedge = surf_ptrs[isurf, 1, 0]
        elif edge == 'v1':
            iedge = surf_ptrs[isurf, 1, 2]
        
        if side is None:
            self.diff['edge'][abs(iedge) - 1, :] = val
        else:
            iside = side if iedge > 0 else 1 - side
            self.diff['edge'][abs(iedge) - 1, iside] = val

    def assemble(self):
        self._compute_indices()
        self._compute_jacobians()

        for vec_type in ['df_str', 'df', 'cp', 
                         'cp_str', 'pt_str', 'pt']:
            self.initialize_vec(vec_type, vec_type, 3)

    def print_info(self):
        num = self._num
        size = self._size

        print
        print 'Topology information'
        print '--------------------'
        print '# surfaces: ', num['surf']
        print '# edges:    ', num['edge']
        print '# vertices: ', num['vert']
        print '# groups:   ', num['group']
        print
        print 'Vector sizes (unique, structured)'
        print '---------------------------------'
        print '# free control points:', size['df'], size['df_str']
        print '# control points:     ', size['cp'], size['cp_str']
        print '# discretized points: ', size['pt'], size['pt_str']
        print

    def initialize_vec(self, name, vec_type, ndim=1):
        vec = self.vec
        size = self._size
        surf_group = self._topo['surf_group']
        num_cp = self._bspline['num_cp'][surf_group-1]
        num_pt = self._bspline['num_pt'][surf_group-1]

        hidden = self.hidden
        if vec_type == 'df_str':
            vec[name] = BSEvecStr(name, size['df_str'],
                                  ndim, num_cp, hidden)
        elif vec_type == 'df':
            vec[name] = BSEvecUns(name, size['df'],
                                  ndim, hidden)
        elif vec_type == 'cp':
            vec[name] = BSEvecUns(name, size['cp'],
                                  ndim, hidden)
        elif vec_type == 'cp_str':
            vec[name] = BSEvecStr(name, size['cp_str'],
                                  ndim, num_cp, hidden, self)
        elif vec_type == 'pt_str':
            vec[name] = BSEvecStr(name, size['pt_str'],
                                  ndim, num_pt, hidden)
        elif vec_type == 'pt':
            vec[name] = BSEvecUns(name, size['pt'],
                                  ndim, hidden)
        else:
            raise Exception('vec_type not recognized')

    def apply_jacobian(self, vec2_name, jac_name, vec1_name):        
        vec1 = self.vec[vec1_name].array
        vec2 = self.vec[vec2_name].array

        jac = self.jac[jac_name]
        for ind in xrange(vec1.shape[1]):
            vec2[:, ind] = jac.dot(vec1[:, ind])

    def compute_projection(self, name, pts, surf_pts=None, ndim=1):
        num = self._num
        size = self._size
        str_indices = self._str_indices
        topo = self._topo
        bspline = self._bspline
        cp_str = self.vec['cp_str'].array
        pt_str = self.vec['pt_str'].array

        pts = numpy.array(pts, order='F')
        if surf_pts is None:
            surf_pts = numpy.linspace(1, num['surf'], num['surf'])
        else:
            surf_pts = numpy.array(surf_pts, int) + 1

        npts = pts.shape[0]
        nsurf_pts = surf_pts.shape[0]
        nrefine = 10
        surfs, ind_u, ind_v \
            = BSElib.computeproj(npts, nsurf_pts, 
                                 size['cp_str'], size['pt_str'], 
                                 num['surf'], num['group'], nrefine,
                                 surf_pts,
                                 str_indices['cp'], str_indices['pt'],
                                 topo['surf_group'], bspline['order'],
                                 bspline['num_cp'], bspline['num_pt'],
                                 cp_str, pt_str, pts)
        self.add_jacobian(name, surfs-1, ind_u, ind_v, ndim)

    def add_jacobian(self, name, surfs, ind_u, ind_v, ndim=1):
        num = self._num
        size = self._size
        str_indices = self._str_indices
        topo = self._topo
        bspline = self._bspline
        vec = self.vec
        jac = self.jac

        surfs = numpy.array(surfs, int) + 1

        npt = surfs.shape[0]
        nnz = BSElib.computesnnz(npt, num['surf'], num['group'],
                                 topo['surf_group'], bspline['order'],
                                 surfs)
        
        quant = ['', '_du', '_dv']
        der_u = [0, 1, 0]
        der_v = [0, 0, 1]
        for ind in xrange(3):
            data, rows, cols \
                = BSElib.computesmtx(der_u[ind], der_v[ind], 
                                     nnz, npt,
                                     num['surf'], num['group'],
                                     str_indices['cp'], topo['surf_group'],
                                     bspline['order'], bspline['num_cp'],
                                     surfs, ind_u, ind_v)
            qname = name + quant[ind]
            jac['d(' + qname + ')/d(cp_str)'] \
                = scipy.sparse.csr_matrix((data, (rows-1, cols-1)),
                                          shape=(npt,
                                                 size['cp_str']))
            vec[qname] = BSEvecUns(qname, npt, ndim, self.hidden)
