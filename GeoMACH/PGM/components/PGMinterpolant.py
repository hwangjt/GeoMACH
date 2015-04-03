"""
GeoMACH interpolant class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.core.PGMcomponent import PGMcomponent
from GeoMACH.PGM.core.PGMproperty import PGMproperty


class PGMinterpolant(PGMcomponent):
    """ Base class for interpolant components """

    def initialize_props(self):
        """ Adds the *X*, *Y*, and *Z* shape properties """
        super(PGMinterpolant, self).initialize_props()

        props = self.props
        for name in self.faces:
            props['shN', name] = PGMproperty()

        self.normals = {}

    def assemble_sizes(self, bse):
        super(PGMinterpolant, self).assemble_sizes(bse)

        props = self.props
        for name in self.faces:
            num_u = self.faces[name]._num_cp_total['u']
            num_v = self.faces[name]._num_cp_total['v']
            props['shN', name].assemble_sizes(num_u, num_v)
            #self.normals[name] = numpy.zeros((num_u, num_v,3)) Check this line

    def compute(self, name):
        props = self.props

        vals_list, rows_list, cols_list = [], [], []
        if name == 'cp_prim':
	    #At this step we write the Jacobian that consists of the normals of each control point
            if len(self.normals) is not 0:
                for fname in self.faces:
                    face = self.faces[fname]
                    num_u = face._num_cp_total['u']
                    num_v = face._num_cp_total['v']

                    vals = self.normals[fname][1:-1, 1:-1, :].flatten() #We cannot apply the normal disturbance to the wireframe points. Otherwise, cp_coons would be disturbed twice as we will add a disturbance over an already disturbed Coons patch. The wireframe points will be regenerated with the averaging.
                    rows = face.vec_inds['cp_prim'][1:-1, 1:-1, :].flatten()
                    cols = numpy.zeros((num_u-2, num_v-2, 3), int)
                    for ind in xrange(3):
                        cols[:, :, ind] = props['shN', fname].vec_inds['prop'][1:-1, 1:-1] #Use the same normal weight for three directions (x, y, and z)
                    cols = cols.flatten()

                    vals_list.append(vals)
                    rows_list.append(rows)
                    cols_list.append(cols)

                    for ind in xrange(3):
                        face.vec_data['cp_prim'][1:-1, 1:-1, ind] = self.normals[fname][1:-1, 1:-1, ind] * \
                                                                    props['shN', fname].vec_data['prop'][1:-1, 1:-1]
            return vals_list, rows_list, cols_list
        elif name == 'cp_bez' or name == 'cp_coons':
	    #At these steps we just need to return an identity matrix as Jacobian, just to carry over the normal parameters to the next step
            for face in self.faces.values():
                num_u = face._num_cp_total['u']
                num_v = face._num_cp_total['v']

                vals = numpy.ones(3 * num_u * num_v)
                inds = face.vec_inds['cp_prim'].flatten()

                vals_list.append(vals)
                rows_list.append(inds)
                cols_list.append(inds)
            return vals_list, rows_list, cols_list

    def compute_normals(self):
        for face in self.faces.values():
            num_u = face._num_cp_total['u']
            num_v = face._num_cp_total['v']

            face.vec_data['cp_coons']
            pts = face.vec_data['cp_coons'].reshape((num_u*num_v, 3), order='C')

            name = 'proj_cp'
            surf_pts = face._surf_indices.flatten()
            surf_pts = surf_pts[surf_pts>=0]
            self._bse.compute_projection(name, pts, surf_pts, ndim=3)

            for qname in [name + '_du', name + '_dv']:
                self._bse.apply_jacobian(qname, 'd(' + qname + ')/d(cp_str)', 'cp_str')

            cross = numpy.cross(self._bse.vec[name + '_du'].array, self._bse.vec[name + '_dv'].array)
            norms = numpy.linalg.norm(cross, axis=1)

            for ind in xrange(3):
                cross[:, ind] /= norms

            self.normals[face._name] = cross.reshape((num_u, num_v, 3), order='C')
