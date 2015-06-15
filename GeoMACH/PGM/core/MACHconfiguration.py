from __future__ import division
import numpy
import scipy.sparse
from collections import OrderedDict
from mpi4py import MPI

from GeoMACH.PGM.core.PGMconfiguration import PGMconfiguration


class MACHconfiguration(PGMconfiguration):

    def __init__(self):
        super(MACHconfiguration, self).__init__()
        self.points = OrderedDict()
        self.diff = OrderedDict()
        self.jacobians = OrderedDict()
        self.updated = {}

    def addPointSet(self, points, pt_name, origConfig=True, **kwargs):
        bse = self._bse
        points = numpy.array(points).real.astype('d')
        self.points[pt_name] = points

        if self.points[pt_name].shape[0] > 0:
	    bse.compute_projection(pt_name, points, ndim=3)
	    self.jacobians[pt_name] = bse.jac['d(' + pt_name + ')/d(cp_str)']
            self.diff[pt_name] = points - self.jacobians[pt_name].dot(bse.vec['cp_str'].array)

        ### FOR DEBUGGING
        bse.apply_jacobian(pt_name, 'd(' + pt_name + ')/d(cp_str)', 'cp_str')
        bse.vec[pt_name].export_tec_scatter('projected_points.dat')
        bse.vec[pt_name].array[:, :] = points
        bse.vec[pt_name].export_tec_scatter('CFD_surf_points.dat')

        self.updated[pt_name] = False

    def setDesignVars(self, dv_dict):
        for dv in self.dvs.values():
            dv.data[:] = numpy.array(numpy.atleast_1d(dv_dict[dv.name]).reshape(dv._shape, order='F')).real.astype('D')

        for pt_name in self.updated:
            self.updated[pt_name] = False

    def getValues(self):
        dv_dict = {}
        for dv in self.dvs.values():
            dv_dict[dv.name] = numpy.array(dv.data.reshape(numpy.prod(dv._shape), order='F'))
        return dv_dict

    def update(self, pt_name, childDelta=True, config=None):
        bse = self._bse
        self.compute_all()
        if self.points[pt_name].shape[0] > 0:
            self.points[pt_name] = self.jacobians[pt_name].dot(bse.vec['cp_str'].array) + self.diff[pt_name]
        self.updated[pt_name] = True
        return numpy.array(self.points[pt_name])

    def pointSetUpToDate(self, pt_name):
        return self.updated[pt_name]

    def getVarNames(self):
        return list(self.dvs.keys())

    def getNDV(self):
        num_dv = 0
        for dv in self.dvs.values():
            num_dv += numpy.prod(dv._shape)
        return num_dv

    def totalSensitivity(self, dfunc_dpt, ptSetName, comm=None, child=False, nDVStore=0, config=None):
        bse = self._bse
        jacs = self._jacs
        pt_name = ptSetName

        if len(dfunc_dpt.shape) == 2:
            dfunc_dpt = numpy.array([dfunc_dpt])
        num_func = dfunc_dpt.shape[0]
        num_dv = self.getNDV()
        num_pt = self.points[pt_name].shape[0]
        ndf = bse._size['df_str']

        dfunc_ddv = numpy.zeros((num_func, num_dv))
        if self.points[pt_name].shape[0] > 0:
	    ddf_ddv = jacs['d(df_str)/d(df_surf)'] * \
                      jacs['d(df_surf)/d(cp_coons)'] * \
                      jacs['d(cp_coons)/d(cp_bez)'] * \
                      jacs['d(cp_bez)/d(cp_prim)'] * \
                      jacs['d(cp_prim)/d(prop)'] * \
                      jacs['d(prop)/d(param)'] * \
                      jacs['d(param)/d(dv)']
            dpt_ddf_s = self.jacobians[pt_name] * \
                bse.jac['d(cp_str)/d(cp)'] * \
                bse.jac['d(cp)/d(df)'] * \
                bse.jac['d(df)/d(df_str)']
            dpt_ddf = scipy.sparse.bmat(
                [
                    [dpt_ddf_s, None, None],
                    [None, dpt_ddf_s, None],
                    [None, None, dpt_ddf_s]
                ],
                format='csc')
            dpt_ddv_T = (dpt_ddf * ddf_ddv).transpose()

            for i in xrange(num_func):
                for k in xrange(3):
                    dfunc_ddv[i, :] += dpt_ddv_T[:, k*num_pt:(k+1)*num_pt].dot(dfunc_dpt[i, :, k])

        if comm:
            dfunc_ddv = comm.allreduce(dfunc_ddv, op=MPI.SUM)

        return self.convertSensitivityToDict(dfunc_ddv)

    def convertSensitivityToDict(self, dfunc_ddv, out1D=False):
        dfunc_ddv_dict = {}

        start, end = 0, 0
        for dv in self.dvs.values():
            end += numpy.prod(dv._shape)
            dfunc_ddv_dict[dv.name] = dfunc_ddv[:, start:end]
            start += numpy.prod(dv._shape)

        return dfunc_ddv_dict

    def addVariablesPyOpt(self, optProb):
        for dv in self.dvs.values():
            optProb.addVarGroup(dv.name, numpy.prod(dv._shape), 'c', 
                                value=dv.value, lower=dv.lower, upper=dv.upper,
                                scale=dv.scale)

    def addConstraintsPyOpt(self, optProb):
        for comp in self.comps.values():
	    for func in comp.funcs.values():
		func.initialize()

	funcsSens = {}
	self.evalFunctionsSens(funcsSens)
	for comp in self.comps.values():
	    for func in comp.funcs.values():
		optProb.addConGroup(func.name, func.size, upper=0,
				    wrt=[dv_name for dv_name in self.dvs],
				    jac={dv_name:funcsSens[func.name][dv_name] for dv_name in self.dvs})

    def evalFunctions(self, funcs):
        self.compute_all()
	for comp in self.comps.values():
	    for func in comp.funcs.values():
		funcs[func.name] = func.get_func()

    def evalFunctionsSens(self, funcsSens):
        bse = self._bse
        jacs = self._jacs

	num_cp = bse.vec['cp_str'].size
	num_dv = self.getNDV()
        ddf_ddv = jacs['d(df_str)/d(df_surf)'] * \
                  jacs['d(df_surf)/d(cp_coons)'] * \
                  jacs['d(cp_coons)/d(cp_bez)'] * \
                  jacs['d(cp_bez)/d(cp_prim)'] * \
                  jacs['d(cp_prim)/d(prop)'] * \
                  jacs['d(prop)/d(param)'] * \
                  jacs['d(param)/d(dv)']
        dcp_ddf_s = bse.jac['d(cp_str)/d(cp)'] * \
                  bse.jac['d(cp)/d(df)'] * \
                  bse.jac['d(df)/d(df_str)']
	dcp_ddf = scipy.sparse.bmat(
            [
                [dcp_ddf_s, None, None],
                [None, dcp_ddf_s, None],
                [None, None, dcp_ddf_s]
            ],
            format='csc')
	dcp_ddv = dcp_ddf * ddf_ddv

	for comp in self.comps.values():
	    for func in comp.funcs.values():
		if func.name not in funcsSens:
		    funcsSens[func.name] = {}

        start = 0
	for dv in self.dvs.values():
            size = numpy.prod(dv._shape)
            ones = numpy.ones(size)
            lins = numpy.array(numpy.linspace(0, size-1, size), int)
            P = scipy.sparse.csr_matrix((ones, (start+lins, lins)),
                                        shape=(num_dv, size))
            # P extracts the columns of the current DV
            start += size
	    for comp in self.comps.values():
		for func in comp.funcs.values():
                    df_dcp = func.get_jacobian()
                    df_ddv = df_dcp * dcp_ddv
		    funcsSens[func.name][dv.name] = (df_ddv * P).todense()
