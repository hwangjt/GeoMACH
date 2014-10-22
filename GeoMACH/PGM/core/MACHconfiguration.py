from __future__ import division
import numpy
import scipy.sparse
from collections import OrderedDict

from GeoMACH.PGM.core.PGMconfiguration import PGMconfiguration


class MACHconfiguration(PGMconfiguration):

    def __init__(self):
        super(MACHconfiguration, self).__init__()
        self.points = OrderedDict()
        self.jacobians = OrderedDict()
        self.updated = {}

    def addPointSet(self, points, pt_name, origConfig=True, **kwargs):
        bse = self._bse
        points = numpy.array(points).real.astype('d')
        self.points[pt_name] = points

        if self.points[pt_name].shape[0] > 0:
	    bse.compute_projection(pt_name, points, ndim=3)
	    self.jacobians[pt_name] = bse.jac['d(' + pt_name + ')/d(cp_str)']

        self.updated[pt_name] = False

    def setDesignVars(self, dv_dict): 
        for dv in self.dvs.values():
            dv.data[:] = numpy.atleast_1d(dv_dict[dv.name]).reshape(dv._shape, order='F')

        for pt_name in self.updated:
            self.updated[pt_name] = False

    def getValues(self):
        dv_dict = {}
        for dv in self.dvs.values():
            dv_dict[dv.name] = dv.data.reshape(numpy.prod(dv._shape), order='F')
        return dv_dict

    def update(self, pt_name, childDelta=True):
        bse = self._bse
        self.compute_all()
        if self.points[pt_name].shape[0] > 0:
            self.points[pt_name] = self.jacobians[pt_name].dot(bse.vec['cp_str'])
        self.updated[pt_name] = True
        return self.points[pt_name]

    def pointSetUpToDate(self, pt_name):
        return self.updated[pt_name]

    def getVarNames(self):
        return list(self.dvs.keys())

    def getNDV(self):
        num_dv = 0
        for dv in self.dvs.values():
            num_dv += numpy.prod(dv._shape)
        return num_dv

    def totalSensitivity(self, dfunc_dpt_T, pt_name, comm=None, child=False, nDVStore=0):
        bse = self._bse

        dfunc_dpt_T = numpy.atleast_3d(dfunc_dpt_T)
        num_func = dfunc_dpt_T.shape[2]
        num_dv = self.getNDV()
        num_pt = self.points[pt_name].shape[0]
        ndf = bse._size['df_str']

        dfunc_ddv_T = numpy.zeros((num_dv, num_func))
        if self.points[pt_name].shape[0] > 0:
	    ddf_ddv = pgm._jacs['d(df_str)/d(df_surf)'] * \
                      pgm._jacs['d(df_surf)/d(cp_coons)'] * \
                      pgm._jacs['d(cp_coons)/d(cp_bez)'] * \
                      pgm._jacs['d(cp_bez)/d(cp_prim)'] * \
                      pgm._jacs['d(cp_prim)/d(prop)'] * \
                      pgm._jacs['d(prop)/d(param)'] * \
                      pgm._jacs['d(param)/d(dv)']
            dpt_ddf = self.jacobians[pt_name] * \
                      bse.jac['d(cp_str)/d(cp)'] * \
                      bse.jac['d(cp)/d(df)'] * \
                      bse.jac['d(df)/d(df_str)']
            dpt_ddv_T = (dpt_ddf * ddf_ddv).transpose()
            dfunc_ddv_T[:, :] = dpt_ddv_T * dfunc_dpt_T

        if comm:
            dfunc_ddv_T = comm.allreduce(dfunc_ddv_T)

        return self.convertSensitivityToDict(dfunc_ddv_T)

    def convertSensitivityToDict(self, dfunc_ddv_T):
        dfunc_ddv = {}

        start, end = 0, 0
        for dv in self.dvs.values():
            end += numpy.prod(dv._shape)
            dfunc_ddv[dv.name] = dfunc_ddv_T[start:end].squeeze().T
            start += numpy.prod(dv._shape)

        return dfunc_ddv

    def addVariablesPyOpt(self, optProb):
        for dv in self.dvs.values():
            optProb.addVarGroup(dv.name, numpy.prod(dv._shape), 'c', 
                                value=dv.value, lower=dv.lower, upper=dv.upper,
                                scale=dv.scale)

    def addConstraintsPyOpt(self, optProb):
        pass

    def evalFunctions(self, funcs):
        pass

    def evalFunctionsSens(self, funcsSens):
        pass
