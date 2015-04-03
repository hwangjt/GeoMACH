"""
GeoMACH PGM configuration class
John Hwang, July 2014
"""
# pylint: disable=E1101

from __future__ import division
import numpy
import time
import scipy.sparse
from collections import OrderedDict

from GeoMACH.BSE.BSEmodel import BSEmodel
from GeoMACH.PGM.core.PGMvec import PGMvec
from GeoMACH.PGM.core.PGMdv import PGMdv


class PGMconfiguration(object):
    """
    Base class for an aircraft configuration

    Attributes
    ----------
    **comps** : ``OrderedDict``
       Dictionary of ``PGMcomp`` objects
    **dvs** : ``OrderedDict``
       Dictionary of ``PGMdv`` objects
    **_vecs** : ``OrderedDict``
       Dictionary of ``PGMvec`` objects

       == ======== ===================== ========= ================
       #  Vec key   Description          Stored in Computed in
       == ======== ===================== ========= ================
       1  dv       design variables      PGMdv     ..
       2  param    parameters            PGMparam  PGMconfiguration
       3  prop     properties            PGMprop   PGMparam
       4  cp_prim  control points        PGMface   PGMcomponent
       .. ..       (primitives)          ..        ..
       5  cp_bez   control points        PGMface   PGMcomponent
       .. ..       (prims. + wireframes) ..        ..
       6  cp_coons control points        PGMface   PGMcomponent
       .. ..       (all)                 ..        ..
       7  df_surf  ordered by surface    PGMsurf   PGMface
       == ======== ===================== ========= ================

    **_jacs** : ``OrderedDict``
       Dictionary of sparse Jacobian matrices

       == ======================
       #  Jac key
       == ======================
       1  d(param)/d(dv)
       2  d(prop)/d(param)
       3  d(cp_prim)/d(prop)
       4  d(cp_bez)/d(cp_prim)
       5  d(cp_coons)/d(cp_bez)
       6  d(df_surf)/d(cp_coons)
       7  d(df_str)/d(df_surf)
       == ======================

    **_bse** : ``BSEmodel``
       GeoMACH B-spline Surface-modeling Engine (BSE) object
    """

    def __init__(self):
        """ Initialize attributes """
        self.comps = OrderedDict()
        self.dvs = OrderedDict()

        self._vecs = OrderedDict()
        self._jacs = OrderedDict()
        self._bse = None

    def initialize(self):
        """
        Builds the PGM model and computes it.

        Returns
        -------
        **bse** : ``BSEmodel``
           A pointer to the BSE object 
        """
        self._define_comps()

        for comp in self.comps.values():
            comp.initialize_props()

        self._define_params()
        self._define_dvs()
        for dv_name in self.dvs:
            self.dvs[dv_name].name = dv_name

        self._initialize_pgm()
        self._initialize_bse()
        self._initialize_pgm()
        self._set_airfoils()
        self._compute_bse()

        return self._bse

    def _define_comps(self):
        """ 
        **Must be implemented in derived class.**

        Adds components to 
        ``self.comps[comp]``.
        """
        pass

    def _define_params(self):
        """ 
        **Must be implemented in derived class.**

        Adds parameters to 
        ``self.comps[comp].props[prop].params``.
        """
        pass

    def _define_dvs(self):
        """
        **Must be implemented in derived class.**

        Adds design variables to 
        ``self.comps[comp].props[prop].params[param].dvs``.
        """
        self.dvs[''] = PGMdv((1))

    def _compute_params(self):
        """
        **Must be implemented in derived class.**

        Computes parameters with given design variable 
        values (if applicable) and returns Jacobian of
        derivatives of parameters w.r.t. DVs
  
        Returns:
        --------
        vals : ``list`` of ``numpy.ndarray``s, float[nnz]
           Values of non-zeros in the Jacobian
        rows : ``list`` of ``numpy.ndarray``s, int[nnz]
           Row indices of non-zeros (0-based)
        cols : ``list`` of ``numpy.ndarray``s, int[nnz]
           Column indices of non-zeros (0-based)
        """
        return [], [], []

    def _set_bspline_options(self):
        """
        **Must be implemented in derived class.**

        Set B-spline order, number of control points,
        and number of points for parallel edges using
        ``self.comps[comp].faces[face].set_num_cp`` and
        ``self.comps[comp].faces[face].set_num_pt``
        """
        pass

    def _set_airfoils(self):
        """
        Optional method to specify airfoils for Wing objects
        """
        pass

    def _initialize_pgm(self):
        """
        Initializes data/computation part of PGM. 3 steps:

        1. Initiates recursive call to compute
           number of control points
           for each component, face, property, etc.
        2. With sizes known, creates vecs and
           computes constant Jacobians
        3. Runs PGM
        """
        for comp in self.comps.values():
            comp.assemble_sizes(self._bse)

        faces = [face
                 for comp in self.comps.values()
                 for face in comp.faces.values()]
        surfs = [surf
                 for comp in self.comps.values()
                 for face in comp.faces.values()
                 for surf in face.surfs.values()
                 if surf.get_vec_shape('') is not None]
        props = [prop
                 for comp in self.comps.values()
                 for prop in comp.props.values()]
        params = [param
                  for comp in self.comps.values()
                  for prop in comp.props.values()
                  for param in prop.params.values()]
        comps = self.comps.values()
        dvs = self.dvs.values()

        vecs = self._vecs
        vecs['dv'] = PGMvec('dv', dvs, [])
        vecs['param'] = PGMvec('param', params, [self])
        vecs['prop'] = PGMvec('prop', props, params)
        vecs['cp_prim'] = PGMvec('cp_prim', faces, comps)
        vecs['cp_bez'] = PGMvec('cp_bez', faces, comps)
        vecs['cp_coons'] = PGMvec('cp_coons', faces, comps)
        vecs['df_surf'] = PGMvec('df_surf', surfs, faces)

        for ind in xrange(6):
            vec1 = vecs.values()[ind]
            vec2 = vecs.values()[ind + 1]
            vec2.set_col_size(vec1.get_size())

        jacs = self._jacs
        jacs['d(param)/d(dv)'] = None
        jacs['d(prop)/d(param)'] = vecs['prop'].compute()
        jacs['d(cp_prim)/d(prop)'] = None
        jacs['d(cp_bez)/d(cp_prim)'] = vecs['cp_bez'].compute()
        jacs['d(cp_coons)/d(cp_bez)'] = vecs['cp_coons'].compute()
        jacs['d(df_surf)/d(cp_coons)'] = vecs['df_surf'].compute()

        self._compute_pgm()

        linspace = numpy.linspace
        size = vecs['df_surf'].get_size()
        vals = numpy.ones(size)

        rows = numpy.zeros(size, int)
        cols = numpy.zeros(size, int)
        irow1, irow2 = 0, 0
        icol1, icol2 = 0, 0
        for surf in surfs:
            num = numpy.prod(surf.get_vec_shape('')[:2])
            lins = numpy.linspace(0, num-1, num).astype(int)

            irow2 += num
            for dim in xrange(3):
                icol2 += num

                rows[icol1:icol2] = dim*size/3 + irow1 + lins
                cols[icol1:icol2] = icol1 + lins

                icol1 += num
            irow1 += num

        jacs['d(df_str)/d(df_surf)'] \
            = scipy.sparse.csr_matrix((vals, (rows, cols)),
                                      shape=(size, size))

    def _initialize_bse(self):
        """ 
        Instantiates and sets up BSE model. 5 steps:

        1. Recursively loops over comps' faces' surfs' to
           collect initial list of surfaces from each surf
           and give each surface a unique global index
        2. Instantiates BSE B-spline surface model
        3. Components set differentiability options
        4. Set user-specified B-spline options
        5. Assemble BSE Jacobians
        """
        surf_index = 0
        surfs_list = []
        for comp in self.comps.values():
            surf_index = comp.initialize_bse(surfs_list, surf_index)

        self._bse = BSEmodel(surfs_list)

        bse = self._bse

        comp_num = 0
        for comp_name in self.comps:
            comp = self.comps[comp_name]
            comp.set_info(comp_name, comp_num, bse)
            comp.set_diff()
            comp.set_hidden_surfaces()
            comp_num += 1
        self._set_bspline_options()

        for isurf in xrange(len(surfs_list)):
            surf = surfs_list[isurf]
            if numpy.average(surf[:, :, 2]) < 0:
                bse.hidden[isurf] = True

        bse.assemble()

    def _apply_jacobian(self, vec2, jac, vec1):
        """
        Performs matrix-vector multiplication : 
        :math:`\textbf{y} = \textbf{A} \textbf{x}`

        Arguments
        ---------
        **vec2** : ``str``
           Key for vector :math:`\textbf{y}`
        **jac** : ``str``
           Key for Jacobian matrix :math:`\textbf{A}`
        **vec1** : ``str``
           Key for vector :math:`\textbf{x}`
        """
        vecs = self._vecs
        jacs = self._jacs
        vecs[vec2]()[:] = jacs[jac].dot(vecs[vec1]())

    def _compute_pgm(self):
        """
        Computes the vectors in sequence.

        - *Nonlinear maps*: the ``PGMvec``s compute method is called
        - *Linear maps*: simply a matrix-vector product
        """
        jacs = self._jacs
        vecs = self._vecs
        mult = self._apply_jacobian

        jacs['d(param)/d(dv)'] = vecs['param'].compute()
        mult('prop', 'd(prop)/d(param)', 'param')
        jacs['d(cp_prim)/d(prop)'] = vecs['cp_prim'].compute()
        mult('cp_bez', 'd(cp_bez)/d(cp_prim)', 'cp_prim')
        mult('cp_coons', 'd(cp_coons)/d(cp_bez)', 'cp_bez')
        mult('df_surf', 'd(df_surf)/d(cp_coons)', 'cp_coons')

    def _compute_bse(self):
        """
        Copies the output of PGM (``vec['df_str']``) to BSE
        and runs the BSE of the computation
        """
        bse = self._bse
        size = bse.vec['df_str'].size

        jac = self._jacs['d(df_str)/d(df_surf)']
        vec = self._vecs['df_surf']()
        bse.vec['df_str'].array[:, :] \
            = jac.dot(vec).reshape((size, 3), order='F')

        bse.apply_jacobian('df', 'd(df)/d(df_str)', 'df_str')
        bse.apply_jacobian('cp', 'd(cp)/d(df)', 'df')
        bse.apply_jacobian('cp_str', 'd(cp_str)/d(cp)', 'cp')
        bse.apply_jacobian('pt_str', 'd(pt_str)/d(cp_str)', 'cp_str')
        bse.apply_jacobian('pt', 'd(pt)/d(pt_str)', 'pt_str')

    def compute_all(self):
        """
        Runs PGM and BSE given the current set of shape variables
        (parameter and design variable values in PGM)
        """
        self._compute_pgm()
        self._compute_bse()

    def compute(self, name):
        """ Compute the *dv* to *param* mapping """
        Das, Dis, Djs = self._compute_params()
        for dv_name in self.dvs:
            dv = self.dvs[dv_name]
            if dv.identity_param is not None:
                comp_name, prop_name, param_name, inds = dv.identity_param
                param = self.comps[comp_name].props[prop_name].params[param_name]
                if inds is None:
                    param.data[:, :] = dv.data[:, :]
                    Das.append([1.0]*param._num_cp['u']*param._num_cp['v'])
                    Dis.append(param.inds.flatten())
                    Djs.append(dv.inds.flatten())
                else:
                    param.data[inds] = dv.data[0]
                    Das.append([1.0])
                    Dis.append([param.inds[inds]])
                    Djs.append([dv.inds[0]])
        return Das, Dis, Djs

    def compute_normals(self):
        for comp in self.comps.values():
            comp.compute_normals()
