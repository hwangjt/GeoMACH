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
        self.compute_properties()

        # Assembles initial surfaces and instantiates OML object
        index_offset = 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                for ind_j in xrange(face.num_surf[1]):
                    for ind_i in xrange(face.num_surf[0]):
                        if face.surf_indices[ind_i, ind_j] != -1:
                            face.surf_indices[ind_i, ind_j] += index_offset
            index_offset = numpy.max(comp.faces.values()[-1].surf_indices) + 1

        self.compute_face_ctrlpts()

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
                            surfs_list.append(face.cp_array[u_start:u_end+1,v_start:v_end+1])
                        u_start += face.num_cp_list[0][ind_i]
                    v_start += face.num_cp_list[1][ind_j]
        self.oml0 = PUBS.PUBS(surfs_list)

        for comp in self.comps.values():
            comp.set_oml(self.oml0)
            #comp.removeHiddenSurfaces()

        self.oml0.write2Tec('test2.dat')
        #exit()
        #self.oml0.plot()

        # Sets comp data and OML properties
        for name in self.comps:
            comp = self.comps[name]
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
        self.compute()

    def initialize_cp_vecs(self):
        # Creates global face-wise cp and index vectors
        num_cp_total = 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                num_cp_total += face.num_cp[0] * face.num_cp[1]
        num_cp_prim = 0
        for comp in self.primitive_comps.values():
            for face in comp.faces.values():
                num_cp_prim += face.num_cp[0] * face.num_cp[1]
        self.num_cp_prim = num_cp_prim
        cp_vec = numpy.zeros(3*num_cp_total)
        index_vec = -numpy.ones(num_cp_total, int)
        cp_indices = numpy.array(
            numpy.linspace(0, num_cp_total-1, num_cp_total), int)

        # Passes views to each face
        start, end = 0, 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                end += face.num_cp[0] * face.num_cp[1]
                face.initialize_cp_data(cp_vec[3*start:3*end],
                                        index_vec[start:end],
                                        cp_indices[start:end])
                start += face.num_cp[0] * face.num_cp[1]

        self.cp_vec = cp_vec
        self.index_vec = index_vec
        self.cp_indices = cp_indices
        self.num_cp_total = num_cp_total
        self.num_cp_prim = num_cp_prim

    def initialize_cp_jacobian(self):
        num_cp_total = self.num_cp_total

        # Sets up face-wise cp to oml's free cp vec mapping
        data0 = numpy.minimum(1, self.index_vec + 1)
        rows0 = numpy.maximum(0, self.index_vec)
        cols0 = numpy.linspace(0, num_cp_total-1, num_cp_total)
        data = numpy.zeros(3*num_cp_total, int)
        rows = numpy.zeros(3*num_cp_total, int)
        cols = numpy.zeros(3*num_cp_total, int)
        for coord in xrange(3):
            data[coord::3] = data0
            rows[coord::3] = 3*rows0 + coord
            cols[coord::3] = 3*cols0 + coord
        cp_jacobian = scipy.sparse.csr_matrix((data, (rows, cols)), 
                                              shape=(3*self.oml0.nQ, 
                                                     3*num_cp_total))
        row_sums = cp_jacobian.dot(numpy.ones(3*num_cp_total, int))
        inv_row_sum = scipy.sparse.diags(1.0/row_sums, 0, format='csr')
        self.cp_jacobian = inv_row_sum.dot(cp_jacobian)

        self.q_vec0 = numpy.zeros((3*self.oml0.nQ))
        self.q_vec = self.q_vec0.reshape((self.oml0.nQ, 3), order='C')

    def initialize_properties(self):
        num_prop_total = 0
        for comp in self.comps.values():
            comp.declare_properties()
            comp.count_properties()
            num_prop_total += comp.size_prop

        prop_vec = numpy.zeros(num_prop_total)
        prop_ind = numpy.array(
            numpy.linspace(0, num_prop_total-1, num_prop_total), int)
            
        start, end = 0, 0
        for comp in self.comps.values():
            end += comp.size_prop
            comp.initialize_properties(prop_vec[start:end],
                                       prop_ind[start:end])
            start += comp.size_prop

        self.prop_vec = prop_vec
        self.prop_ind = prop_ind

    def initialize_parameters(self):
        num_param_total = 0
        for comp in self.comps.values():
            for prop in comp.props.values():
                for param in prop.params.values():
                    num_param_total += 3 * param.mu * param.mv

        param_vec = numpy.zeros(num_param_total)
        param_ind = numpy.array(
            numpy.linspace(0, num_param_total-1, num_param_total), int)
            
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

        Das = []
        Dis = []
        Djs = []
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

        self.prop_jac = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                                shape=(self.prop_vec.shape[0],
                                                       self.param_vec.shape[0]))


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
        time0 = time.time()
        self.compute_properties()
        self.compute_face_ctrlpts()
        self.compute_free_ctrlpt_vector()
        time1 = time.time()
        self.oml0.computePoints()
        time2 = time.time()
        print time2-time1, time1-time0

    def compute_properties(self):
        """ Computes section properties from parameters """
        self.prop_vec[:] = self.prop_jac.dot(self.param_vec)

    def compute_face_ctrlpts(self, full=True, name0=None):
        """ Computes face control points from section properties """
        def linspace(n):
            return numpy.array(numpy.linspace(0,n-1,n), int)
        
        self.cp_vec[:] = 0.0

        # Step 1: compute primitives' CPs
        for comp in self.primitive_comps.values():
            comp.computeQs()
        face_cps = linspace(3*self.num_cp_prim)

        # Step 2: compute interpolants' wireframe CPs
        Das = []
        Dis = []
        Djs = []
        for comp in self.interpolant_comps.values():
            Da, Di, Dj = comp.compute_cp_wireframe()
            Das.append(Da)
            Dis.append(Di)
            Djs.append(Dj)
        Da = numpy.concatenate(Das + [numpy.ones(face_cps.shape[0])])
        Di = numpy.concatenate(Dis + [face_cps])
        Dj = numpy.concatenate(Djs + [face_cps])
        D1 = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                     shape=(self.cp_vec.shape[0],
                                            self.cp_vec.shape[0]))
        self.cp_vec[:] = D1.dot(self.cp_vec)
        face_grid_cps = numpy.unique(Di)

        # Step 3: compute interpolants' surface CPs
        Das = []
        Dis = []
        Djs = []
        for comp in self.interpolant_comps.values():
            Da, Di, Dj = comp.compute_cp_surfs()
            Das.append(Da)
            Dis.append(Di)
            Djs.append(Dj)
        Da = numpy.concatenate(Das + [numpy.ones(face_grid_cps.shape[0])])
        Di = numpy.concatenate(Dis + [face_grid_cps])
        Dj = numpy.concatenate(Djs + [face_grid_cps])
        D2 = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                     shape=(self.cp_vec.shape[0],
                                            self.cp_vec.shape[0]))
        self.cp_vec[:] = D2.dot(self.cp_vec)

    def compute_free_ctrlpt_vector(self):
        """ Computes vector of free control points from face control points """
        self.q_vec0[:] = self.cp_jacobian.dot(self.cp_vec)
        self.oml0.Q[:, :3] = self.q_vec

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
