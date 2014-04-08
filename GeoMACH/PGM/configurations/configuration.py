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
        z_coord = 0
        for comp in self.primitive_comps.values():
            for surfs in comp.Ps:
                surfs[:, :, :] += z_coord
            z_coord += 4

        # Adds interpolant components
        self.interpolant_comps = self.define_interpolant_comps()
        self.comps.update(self.interpolant_comps)

        # Assembles initial surfaces and instantiates OML object
        surfs_list = []
        index_offset = 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                for ind_j in xrange(face.num_surf[1]):
                    for ind_i in xrange(face.num_surf[0]):
                        if face.surf_indices[ind_i, ind_j] != -1:
                            face.surf_indices[ind_i, ind_j] += index_offset
            index_offset = numpy.max(comp.faces.values()[-1].surf_indices) + 1
            surfs_list.extend(comp.Ps)
            comp.Ps = []
        self.oml0 = PUBS.PUBS(surfs_list)

        # Sets comp data and OML properties
        for name in self.comps:
            comp = self.comps[name]
            comp.name = name
            comp.set_oml(self.oml0)
            comp.setDOFs()
        self.set_oml_resolution()
        self.oml0.update()
        for comp in self.comps.values():
            for face in comp.faces.values():
                face.compute_num_cp()

        # Creates global face-wise cp and index vectors
        num_cp_total = 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                num_cp_total += face.num_cp[0] * face.num_cp[1]
        cp_vec = numpy.zeros(3*num_cp_total)
        index_vec = -numpy.ones(num_cp_total, int)

        # Passes views to each face
        start, end = 0, 0
        for comp in self.comps.values():
            for face in comp.faces.values():
                end += face.num_cp[0] * face.num_cp[1]
                face.initializeDOFmappings(cp_vec[3*start:3*end],
                                           index_vec[start:end])
                start += face.num_cp[0] * face.num_cp[1]

        # Sets up face-wise cp to oml's free cp vec mapping
        data0 = numpy.minimum(1, index_vec + 1)
        rows0 = numpy.maximum(0, index_vec)
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

        self.cp_vec = cp_vec
        self.q_vec0 = numpy.zeros((3*self.oml0.nQ))
        self.q_vec = self.q_vec0.reshape((self.oml0.nQ, 3), order='C')

        for comp in self.comps.values():
            comp.removeHiddenDOFs()

        num_prop_total = 0
        for comp in self.comps.values():
            comp.declare_properties()
            comp.count_properties()
            num_prop_total += comp.size_prop

        prop_vec = numpy.zeros(num_prop_total)
        prop_index_vec = numpy.array(
            numpy.linspace(0, num_prop_total-1, num_prop_total), int)
            
        start, end = 0, 0
        for comp in self.comps.values():
            end += comp.size_prop
            comp.initialize_properties(prop_vec[start:end],
                                       prop_index_vec[start:end])
            start += comp.size_prop

        self.define_oml_parameters()
        self.compute()

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
        for comp in self.comps.values():
            comp.computeVs()

    def compute_face_ctrlpts(self, full=True, name0=None):
        """ Computes face control points from section properties """
        for comp in self.primitive_comps.values():
            comp.computeQs()
        for comp in self.interpolant_comps.values():
            comp.computeQs()
        #if full:
        #    for comp in self.comps.values():
        #        comp.computeQs()
        #else:
        #    for name in self.interpolant_comps.keys():
        #        if name != name0:
        #            self.comps[name].computeQs()

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
