"""
GeoMACH configuration class
John Hwang, April 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.PUBS import PUBS


class Configuration(object):
    """ Base class for an aircraft configuration """
    
    def __init__(self):
        """ Initializes the outer mold line (OML) and parametrization """
        self.comps = OrderedDict()
        self.dvs = OrderedDict()

        # Adds primitive components and separate them
        self.primitive_comps = self.define_primitive_comps()
        self.comps.update(self.primitive_comps)

        # Adds interpolant components
        self.interpolant_comps = self.define_interpolant_comps()
        self.comps.update(self.interpolant_comps)

        for name in self.comps:
            self.comps[name].name = name

        self.initialize()
        self.initialize_dvs()
        self.compute_cp()

        self.initialize_oml()
        for comp in self.comps.values():
            comp.set_oml(self.oml0)
            comp.setDOFs()
        self.set_oml_resolution()
        self.oml0.update()
        for comp in self.comps.values():
            for face in comp.faces.values():
                face.compute_num_cp()

        self.initialize()
        self.initialize_cp_jacobian()
        self.compute()

    def initialize(self):
        self.initialize_cp_vecs()
        self.initialize_properties()
        self.define_oml_parameters()
        self.initialize_parameters()
        self.prop_vec[:] = self.dprop_dparam.dot(self.param_vec) 
        self.compute_cp_prim()
        self.initialize_cp_bezier()
        self.initialize_cp_coons()  

    def add_dv(self, name, shape, val=0.0, lower=None, upper=None, scale=1.0):
        self.dvs[name] = DV(name, shape, val, lower, upper, scale)

    def initialize_dvs(self):
        self.define_dvs()

        num_dv = 0
        for dv in self.dvs.values():
            num_dv += dv.size

        dv_vec = numpy.zeros(num_dv)
        dv_ind = numpy.array(numpy.linspace(0, num_dv-1, num_dv), int)
            
        start, end = 0, 0
        for dv in self.dvs.values():
            end += dv.size
            dv.vec = dv_vec[start:end].reshape(dv.shape, order='F')
            dv.ind = dv_ind[start:end].reshape(dv.shape, order='F')
            start += dv.size
            dv.vec[:] = dv.val

        self.num_dv = num_dv
        self.dv_vec = dv_vec

    def compute_dvs(self):
        Das, Dis, Djs = self.apply_dvs()
        for i in xrange(len(Das)):
            Das[i] = numpy.atleast_1d(Das[i])
        for i in xrange(len(Dis)):
            Dis[i] = numpy.atleast_1d(Dis[i])
        for i in xrange(len(Djs)):
            Djs[i] = numpy.atleast_1d(Djs[i])

        Da = numpy.concatenate(Das)
        Di = numpy.concatenate(Dis)
        Dj = numpy.concatenate(Djs)

        self.dparam_ddv = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                                  shape=(self.num_param,
                                                         self.num_dv))

    def initialize_oml(self):
        index_offset = 0
        surfs_list = []
        for comp in self.comps.values():
            for face in comp.faces.values():
                v_start, v_end = 0, 0
                for ind_j in range(face.num_surf[1]):
                    v_end += face.num_cp_list[1][ind_j]
                    u_start, u_end = 0, 0
                    for ind_i in range(face.num_surf[0]):
                        u_end += face.num_cp_list[0][ind_i]
                        if face.surf_indices[ind_i, ind_j] != -1:
                            face.surf_indices[ind_i, ind_j] += index_offset
                            surfs_list.append(face.cp_array[u_start:u_end+1,v_start:v_end+1])
                        u_start += face.num_cp_list[0][ind_i]
                    v_start += face.num_cp_list[1][ind_j]
            index_offset = numpy.max(comp.faces.values()[-1].surf_indices) + 1
        self.oml0 = PUBS.PUBS(surfs_list)

    def initialize_cp_vecs(self):
        # Creates global face-wise cp and index vectors
        num_cp = 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                num_cp += 3 * face.num_cp[0] * face.num_cp[1]

        cp_vec = numpy.zeros(num_cp)
        cp_indices = numpy.array(numpy.linspace(0, num_cp-1, num_cp), int)
        index_vec = -numpy.ones(num_cp, int)

        # Passes views to each face
        start, end = 0, 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                end += 3 * face.num_cp[0] * face.num_cp[1]
                face.initialize_cp_data(cp_vec[start:end],
                                        cp_indices[start:end],
                                        index_vec[start:end])
                start += 3 * face.num_cp[0] * face.num_cp[1]

        self.cp_vec = cp_vec
        self.cp_indices = cp_indices
        self.index_vec = index_vec
        self.num_cp = num_cp

    def initialize_cp_jacobian(self):
        for comp in self.comps.values():
            for face in comp.faces.values():
                face.initializeDOFmappings()

        self.num_cp_dof = 3 * self.oml0.nQ

        num_cp = self.num_cp
        num_cp_dof = self.num_cp_dof

        # Sets up face-wise cp to oml's free cp vec mapping
        data = numpy.minimum(1, self.index_vec + 1)
        rows = numpy.maximum(0, self.index_vec)
        cols = numpy.array(numpy.linspace(0, num_cp-1, num_cp), int)
        cp_jacobian = scipy.sparse.csr_matrix((data, (rows, cols)), 
                                              shape=(num_cp_dof,
                                                     num_cp))
        row_sums = cp_jacobian.dot(numpy.ones(num_cp, int))
        lins = numpy.array(numpy.linspace(0, num_cp_dof-1, num_cp_dof), int)
        inv_row_sum = scipy.sparse.csr_matrix((1.0/row_sums, (lins, lins)),
                                              shape=(num_cp_dof, num_cp_dof))
        self.ddof_dcp = inv_row_sum.dot(cp_jacobian)

        self.dof_vec0 = numpy.zeros((num_cp_dof))
        self.dof_vec = self.dof_vec0.reshape((self.oml0.nQ, 3), order='C')

    def initialize_properties(self):
        num_prop = 0
        for comp in self.comps.values():
            comp.declare_properties()
            comp.count_properties()
            num_prop += comp.size_prop

        prop_vec = numpy.zeros(num_prop)
        prop_ind = numpy.array(numpy.linspace(0, num_prop-1, num_prop), int)
            
        start, end = 0, 0
        for comp in self.comps.values():
            end += comp.size_prop
            comp.initialize_properties(prop_vec[start:end],
                                       prop_ind[start:end])
            start += comp.size_prop

        self.prop_vec = prop_vec
        self.prop_ind = prop_ind
        self.num_prop = num_prop

    def initialize_parameters(self):
        num_param = 0
        for comp in self.comps.values():
            for prop in comp.props.values():
                for param in prop.params.values():
                    num_param += param.count_parameters()

        param_vec = numpy.zeros(num_param)
        param_ind = numpy.array(numpy.linspace(0, num_param-1, num_param), int)
            
        start, end = 0, 0
        for comp in self.comps.values():
            for prop in comp.props.values():
                for param in prop.params.values():
                    end += param.count_parameters()
                    param.initialize_parameters(param_vec[start:end],
                                                param_ind[start:end])
                    start += param.count_parameters()

        self.param_vec = param_vec
        self.param_ind = param_ind
        self.num_param = num_param

        Das, Dis, Djs = [], [], []
        for comp in self.comps.values():
            for prop in comp.props.values():
                for param in prop.params.values():
                    Da, Di, Dj = param.compute_jacobian(prop.prop_ind)
                    Das.append(Da)
                    Dis.append(Di)
                    Djs.append(Dj)
        Da = numpy.concatenate(Das)
        Di = numpy.concatenate(Dis)
        Dj = numpy.concatenate(Djs)

        self.dprop_dparam = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                                shape=(self.num_prop,
                                                       self.num_param))

    def initialize_cp_bezier(self):
        face_cps = self.face_cps

        Das, Dis, Djs = [], [], []
        for comp in self.interpolant_comps.values():
            Da, Di, Dj = comp.compute_cp_wireframe()
            Das.append(Da)
            Dis.append(Di)
            Djs.append(Dj)
        Da = numpy.concatenate(Das + [numpy.ones(face_cps.shape[0])])
        Di = numpy.concatenate(Dis + [face_cps])
        Dj = numpy.concatenate(Djs + [face_cps])
        self.dbezier_dprim = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                                     shape=(self.num_cp,
                                                            self.num_cp))
        self.face_grid_cps = numpy.unique(Di)

    def initialize_cp_coons(self):
        face_grid_cps = self.face_grid_cps

        Das, Dis, Djs = [], [], []
        for comp in self.interpolant_comps.values():
            Da, Di, Dj = comp.compute_cp_surfs()
            Das.append(Da)
            Dis.append(Di)
            Djs.append(Dj)
        Da = numpy.concatenate(Das + [numpy.ones(face_grid_cps.shape[0])])
        Di = numpy.concatenate(Dis + [face_grid_cps])
        Dj = numpy.concatenate(Djs + [face_grid_cps])
        self.dcoons_dbezier = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                                      shape=(self.num_cp,
                                                             self.num_cp))

    def define_primitive_comps(self):
        """ Virtual method; must be implemented in derived class """
        pass

    def define_interpolant_comps(self):
        """ Virtual method; must be implemented in derived class """
        pass

    def set_oml_resolution(self):
        """ Virtual method; must be implemented in derived class """
        pass

    def define_oml_parameters(self):
        """ Virtual method; must be implemented in derived class """
        pass

    def compute(self):
        """ Computes OML points based on current parameter values """
        self.compute_cp()
        self.compute_oml()
        self.jac = self.ddof_dcp.dot(self.dcoons_dbezier.dot(self.dbezier_dprim.dot(self.dprim_dprop.dot(self.dprop_dparam.dot(self.dparam_ddv)))))

    def compute_cp(self):
        time0 = time.time()
        self.compute_dvs()
        self.prop_vec[:] = self.dprop_dparam.dot(self.param_vec) 
        self.cp_vec[:] = 0.0
        self.compute_cp_prim()
        self.cp_vec[:] = self.dbezier_dprim.dot(self.cp_vec)
        self.cp_vec[:] = self.dcoons_dbezier.dot(self.cp_vec)
        #print 'Compute CP:', time.time() - time0

    def compute_oml(self):
        time0 = time.time()
        self.dof_vec0[:] = self.ddof_dcp.dot(self.cp_vec)
        self.oml0.Q[:, :3] = self.dof_vec
        self.oml0.computePoints()
        #print 'Compute OML:', time.time() - time0

    def compute_cp_prim(self):
        """ Computes face control points from section properties """ 
        Das, Dis, Djs = [], [], []
        for comp in self.primitive_comps.values():
            Da, Di, Dj = comp.computeQs()
            Das.extend(Da)
            Dis.extend(Di)
            Djs.extend(Dj)
        Da = numpy.concatenate(Das)
        Di = numpy.concatenate(Dis)
        Dj = numpy.concatenate(Djs)
        self.dprim_dprop = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                                   shape=(self.num_cp,
                                                          self.num_prop))
        self.face_cps = numpy.unique(Di)

    def test_derivatives_dv(self):
        self.compute()

        dof = numpy.array(self.dof_vec0)
        deriv_AN = numpy.array(self.dof_vec0)
        deriv_FD = numpy.array(self.dof_vec0)

        h = 1e-3
        in_vec = numpy.zeros(self.num_dv)
        ind = 0
        for dv in self.dvs.values():
            for k in xrange(dv.size):
                in_vec[ind] = 1.0
                deriv_AN[:] = self.jac.dot(in_vec)
                in_vec[ind] = 0.0

                self.dv_vec[ind] += h
                self.compute()
                deriv_FD[:] = (self.dof_vec0 - dof) / h
                self.dv_vec[ind] -= h

                ind += 1

                print '%6s %3i %17.10e %17.10e %17.10e' % \
                    (dv.name, k, 
                     numpy.linalg.norm(deriv_FD),
                     numpy.linalg.norm(deriv_AN),
                     numpy.linalg.norm(deriv_FD - deriv_AN) / numpy.linalg.norm(deriv_FD))

    def test_derivatives(self, comp_names=None, prop_names=None):
        if comp_names is None:
            comp_names = self.comps.keys()
        if prop_names is None:
            prop_names = ['pos', 'rot', 'scl']

        self.compute()
        dof = numpy.array(self.dof_vec0)
        deriv_AN = numpy.array(self.dof_vec0)
        deriv_FD = numpy.array(self.dof_vec0)
        jac = self.ddof_dcp.dot(self.dcoons_dbezier.dot(self.dbezier_dprim.dot(self.dprim_dprop.dot(self.dprop_dparam))))

        h = 1e-3
        in_vec = numpy.zeros(self.num_param)
        for comp_name in comp_names:
            comp = self.comps[comp_name]
            for prop_name in prop_names + \
                    [('shX', name) for name in comp.faces] + \
                    [('shY', name) for name in comp.faces] + \
                    [('shZ', name) for name in comp.faces]:
                if prop_name in comp.props:
                    prop = comp.props[prop_name]
                    for param_name in prop.params:
                        param = prop.params[param_name]
                        mu, mv = param.mu, param.mv
                        for i in range(mu):
                            for j in range(mv):
                                in_vec[:] = 0.0
                                in_vec[param.param_ind[i,j,0]] = 1.0
                                deriv_AN[:] = jac.dot(in_vec)
                                param.param_vec[i,j,0] += h     
                                self.compute()
                                deriv_FD[:] = (self.dof_vec0 - dof) / h
                                param.param_vec[i,j,0] -= h            
                                print '%6s %6s %6s %3i %3i %17.10e %17.10e %17.10e' % \
                                    (comp_name, prop_name, param_name, i, j, 
                                     numpy.linalg.norm(deriv_FD),
                                     numpy.linalg.norm(deriv_AN),
                                     numpy.linalg.norm(deriv_FD - deriv_AN) / numpy.linalg.norm(deriv_FD))


class DV(object):

    def __init__(self, name, shape, val, lower, upper, scale):
        self.name = name
        self.shape = shape
        self.val = val
        self.lower = lower
        self.upper = upper
        self.scale = scale
        self.size = numpy.prod(shape)
