from __future__ import division
import numpy
import scipy.sparse
from collections import OrderedDict

from GeoMACH.PGM.configurations import Configuration


class ConfigurationMACH(Configuration):

    def __init__(self):
        super(ConfigurationMACH, self).__init__()
        self.points = OrderedDict()
        self.jacobians = OrderedDict()
        self.updated = {}

    def addPointSet(self, points, pt_name, origConfig=True, **kwargs):
        points = numpy.array(points).real.astype('d')
        self.points[pt_name] = points

        surf, indu, indv = self.oml0.evaluateProjection(points)
        self.jacobians[pt_name] = self.oml0.evaluateBases(surf, indu, indv)

        self.updated[pt_name] = False

    def setDesignVars(self, dv_dict):
        for dv in self.dvs.values():
            dv.vec[:] = numpy.atleast_1d(dv_dict[dv.name]).reshape(dv.shape, order='F')

        for pt_name in self.updated:
            self.updated[pt_name] = False

    def getValues(self):
        dv_dict = {}
        for dv in self.dvs.values():
            dv_dict[dv.name] = dv.vec.reshape(dv.size, order='F')
        return dv_dict

    def update(self, pt_name, childDelta=True):
        self.compute()
        self.points[pt_name] = self.jacobians[pt_name].dot(self.oml0.C[:,:3])
        self.updated[pt_name] = True
        return self.points[pt_name]

    def pointSetUpToDate(self, pt_name):
        return self.updated[pt_name]

    def getVarNames(self):
        return list(self.dvs.keys())

    def getNDV(self):
        num_dv = 0
        for dv in self.dvs.values():
            num_dv += dv.size
        return num_dv

    def totalSensitivity(self, dfunc_dpt_T, ptSetName, comm=None, child=False, nDVStore=0):
        pt_name = ptSetName

        dfunc_dpt_T = numpy.atleast_3d(dfunc_dpt_T)
        num_func = dfunc_dpt_T.shape[2]
        num_dv = self.getNDV()
        num_pt = self.points[pt_name].shape[0]
        nQ = self.oml0.nQ

        dfunc_ddv_T = numpy.zeros((num_dv, num_func))
        ones = numpy.ones(nQ)
        lins = numpy.array(numpy.linspace(0, nQ-1, nQ), int)
        for k in range(3):
            P = scipy.sparse.csr_matrix((ones, (lins, 3*lins + k)), 
                                        shape=(nQ, 3*nQ))
            dfunc_ddv_T[:,:] += self.jacobians[pt_name].dot(self.oml0.M.dot(P.dot(self.jac))).transpose().dot(dfunc_dpt_T[:,k,:])

        if comm:
            dfunc_ddv_T = comm.allreduce(dfunc_ddv_T)           

        return self.convertSensitivityToDict(dfunc_ddv_T)

    def convertSensitivityToDict(self, dfunc_ddv_T):
        dfunc_ddv = {}

        start, end = 0, 0
        for dv in self.dvs.values():
            end += dv.size
            dfunc_ddv[dv.name] = dfunc_ddv_T[start:end].squeeze().T
            start += dv.size

        return dfunc_ddv

    def addVariablesPyOpt(self, optProb):
        for dv in self.dvs.values():
            optProb.addVarGroup(dv.name, dv.size, 'c', 
                                value=dv.val, lower=dv.lower, upper=dv.upper,
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
        for comp in self.comps.values():
            for func in comp.funcs.values():
                funcs[func.name] = func.get_func()

    def evalFunctionsSens(self, funcsSens):
        M = scipy.sparse.bmat([[self.oml0.M, None, None],
                               [None, self.oml0.M, None],
                               [None, None, self.oml0.M]], format='csc')
        for comp in self.comps.values():
            for func in comp.funcs.values():
                if func.name not in funcsSens:
                    funcsSens[func.name] = {}

        start = 0
        for dv in self.dvs.values():
            ones = numpy.ones(dv.size)
            lins = numpy.array(numpy.linspace(0, dv.size-1, dv.size), int)
            P = scipy.sparse.csr_matrix((ones, (start+lins, lins)), 
                                        shape=(self.num_dv, dv.size))
            for comp in self.comps.values():
                for func in comp.funcs.values():
                    funcsSens[func.name][dv.name] = func.get_jacobian().dot(M.dot(self.jac.dot(P)))
            start += dv.size
