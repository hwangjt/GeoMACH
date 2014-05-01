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

        # Adds primitive components and separate them
        self.primitive_comps = self.define_primitive_comps()
        self.comps.update(self.primitive_comps)

        # Adds interpolant components
        self.interpolant_comps = self.define_interpolant_comps()
        self.comps.update(self.interpolant_comps)

        for name in self.comps:
            self.comps[name].name = name

        self.initialize_cp_vecs()
        self.initialize_properties()
        self.define_oml_parameters()
        self.initialize_parameters()
        self.prop_vec[:] = self.dprop_dparam.dot(self.param_vec) 
        self.compute_cp_prim()
        self.initialize_cp_bezier()
        self.initialize_cp_coons()
        self.compute_cp()

        # Initializes surfaces and instantiates OML object
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

        # Sets comp data and OML properties
        for comp in self.comps.values():
            comp.set_oml(self.oml0)
            comp.setDOFs()
        self.set_oml_resolution()
        self.oml0.update()
        for comp in self.comps.values():
            for face in comp.faces.values():
                face.compute_num_cp()

        self.initialize_cp_vecs()
        for comp in self.comps.values():
            for face in comp.faces.values():
                face.initializeDOFmappings()
        self.initialize_cp_jacobian()
        self.initialize_properties()
        self.define_oml_parameters()
        self.initialize_parameters()
        self.prop_vec[:] = self.dprop_dparam.dot(self.param_vec) 
        self.compute_cp_prim()
        self.initialize_cp_bezier()
        self.initialize_cp_coons()
        self.compute()

    def initialize_cp_vecs(self):
        # Creates global face-wise cp and index vectors
        num_cp = 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                num_cp += 3 * face.num_cp[0] * face.num_cp[1]
        num_cp_prim = 0
        for comp in self.primitive_comps.values():
            for face in comp.faces.values():
                num_cp_prim += 3 * face.num_cp[0] * face.num_cp[1]

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
        self.num_cp_prim = num_cp_prim

    def initialize_cp_jacobian(self):
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
        inv_row_sum = scipy.sparse.diags(1.0/row_sums, 0, format='csr')
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
                    num_param += 3 * param.mu * param.mv

        param_vec = numpy.zeros(num_param)
        param_ind = numpy.array(numpy.linspace(0, num_param-1, num_param), int)
            
        start, end = 0, 0
        for comp in self.comps.values():
            for prop in comp.props.values():
                for param in prop.params.values():
                    end += 3 * param.mu * param.mv
                    param.initialize_parameters(param_vec[start:end],
                                                param_ind[start:end])
                    start += 3 * param.mu * param.mv

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

    def compute_cp(self):
        time0 = time.time()
        self.prop_vec[:] = self.dprop_dparam.dot(self.param_vec) 
        self.cp_vec[:] = 0.0
        self.compute_cp_prim()
        self.cp_vec[:] = self.dbezier_dprim.dot(self.cp_vec)
        self.cp_vec[:] = self.dcoons_dbezier.dot(self.cp_vec)
        print 'Compute CP:', time.time() - time0

    def compute_oml(self):
        time0 = time.time()
        self.dof_vec0[:] = self.ddof_dcp.dot(self.cp_vec)
        self.oml0.Q[:, :3] = self.dof_vec
        self.oml0.computePoints()
        print 'Compute OML:', time.time() - time0

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
        num_cp_prim = self.num_cp_prim
        self.face_cps = numpy.array(numpy.linspace(0,num_cp_prim-1,num_cp_prim), int)

    def get_derivatives(self, comp_name, par_name, ind,
                        clean=True, useFD=False, step=1e-5):
        comp = self.comps[comp_name]
        par = comp.params[par_name]
        var = par.var
        self.compute_properties()
        self.compute_face_ctrlpts()
        self.compute_free_ctrlpt_vector()
        V0 = numpy.array(comp.properties[var])
        Q0 = numpy.array(self.oml0.Q[:, :3])
        if useFD:
            par.P[ind[0], ind[1], 0] += step
            self.compute_properties()
            self.compute_face_ctrlpts()
            par.P[ind[0], ind[1], 0] -= step
        else:
            step = 1.0
            par.P[ind[0], ind[1], 0] += step
            self.compute_properties()
            par.P[ind[0], ind[1], 0] -= step
            dV = comp.properties[var] - V0
            self.compute_properties()
            comp.setDerivatives(var, dV)
            self.compute_face_ctrlpts(False, comp_name)
        self.compute_free_ctrlpt_vector()
        res = (self.oml0.Q[:, :3] - Q0)/step
        if clean:
            self.compute()
        return res

    def test_derivatives(self, c, ps=[]):
        self.compute()

        comp = self.comps[c]
        if ps == []:
            ps = comp.params.keys()

        self.compute()
        step = 1e-5
        for p in ps:
            par = comp.params[p]
            var = par.var
            if not (var in ['nor', 'ogn', 'flt']):
                ni, nj = par.P.shape[:2]
                for i in range(ni):
                    for j in range(nj):
                        ind = (i, j)
                        drv1 = self.get_derivatives(c, p, ind, clean=False)
                        drv2 = self.get_derivatives(c, p, ind, clean=False,
                                                    useFD=True, step=step)
                        norm0 = numpy.linalg.norm(drv1)
                        norm0 = 1.0 if norm0 == 0 else norm0
                        error = numpy.linalg.norm(drv2-drv1)/norm0
                        good = 'O' if error < 1e-4 else 'X'
                        print good, ' ', c, ' ', p, ' ', ind, ' ', error
        self.compute()

"""
    def get_derivatives0(self, comp, var, ind, clean=True, useFD=False, step=1e-5):
        self.compute_face_ctrlpts()
        self.compute_free_ctrlpt_vector()
        Q0 = numpy.array(self.oml0.Q[:, :3])
        if useFD:
            self.comps[comp].properties[var][ind] += step
            self.compute_face_ctrlpts()
            self.comps[comp].properties[var][ind] -= step
        else:
            self.comps[comp].setDerivatives(var, ind)
            self.compute_face_ctrlpts(False, comp)
            step = 1.0
        self.compute_free_ctrlpt_vector()
        res = (self.oml0.Q[:, :3] - Q0)/step
        if clean:
            self.computePoints()
        return res

    def test_derivatives0(self, comp, properties=[]):
        self.compute()
        if properties == []:
            properties = self.comps[comp].properties.keys()
        step = 1e-5
        for var in properties:
            if not (var in ['nor', 'origin', 'fillet']):
                dat = self.comps[comp].properties[var]
                for ind, val in numpy.ndenumerate(dat):
                    ind = ind[0] if len(ind) == 1 else ind
                    drv1 = self.get_derivatives(comp, var, ind, clean=False)
                    drv2 = self.get_derivatives(comp, var, ind, clean=False,
                                               useFD=True, step=step)
                    norm0 = numpy.linalg.norm(drv2)
                    norm0 = 1.0 if norm0 == 0 else norm0
                    error = numpy.linalg.norm(drv2-drv1)/norm0
                    good = 'O' if error < 1e-4 else 'X'
                    print good, ' ', comp, ' ', var, ' ', ind, ' ', error
        self.compute()
"""
