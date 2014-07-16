"""
GeoMACH primitive class
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


class PGMprimitive(PGMcomponent):
    """ Base class for non-interpolating components """

    def __init__(self):
        super(PGMprimitive, self).__init__()

    def initialize_props(self):
        super(PGMprimitive, self).initialize_props()

        props = self.props
        props['scl'] = PGMproperty()
        props['pos'] = PGMproperty()
        props['rot'] = PGMproperty()
        props['ogn'] = PGMproperty()
        props['nor'] = PGMproperty()
        props['flt'] = PGMproperty()

    def assemble_sizes(self, bse):
        super(PGMprimitive, self).assemble_sizes(bse)

        num = self.faces.values()[0]._num_cp_total['v']
        props = self.props
        props['scl'].assemble_sizes(num, 3)
        props['pos'].assemble_sizes(num, 3)
        props['rot'].assemble_sizes(num, 3)
        props['ogn'].assemble_sizes(num, 3)
        props['nor'].assemble_sizes(num, 3)
        props['flt'].assemble_sizes(num, 4)

    def compute(self, name):
        props = self.props

        if name == 'cp_prim':
            data = {}
            inds = {}
            for pname in ['pos', 'rot', 'nor',
                          'ogn', 'scl']:
                data[pname] = props[pname].vec_data['prop']
                inds[pname] = props[pname].vec_inds['prop']

            vals_list, rows_list, cols_list = [], [], []
            for fname in self.faces:
                face = self.faces[fname]

                sh_data = {}
                sh_inds = {}
                
                pnames = ['shX', 'shY', 'shZ']
                for ind in xrange(3):
                    pname = pnames[ind]
                    sh_data[pname] \
                        = props[pname, fname].vec_data['prop'] \
                        + self._shapes[fname][:, :, ind]
                    sh_inds[pname] \
                        = props[pname, fname].vec_inds['prop']

                num_u = face._num_cp_total['u']
                num_v = face._num_cp_total['v']

                cp_data = face.vec_data['cp_prim']
                cp_inds = face.vec_inds['cp_prim']

                nnz = num_u*num_v*3*3*6
                cp_data[:, :, :], vals, rows, cols \
                    = PGMlib.computesections(self._ax1, self._ax2, num_u, num_v, nnz, data['ogn'], data['nor'], data['pos'], data['rot'], data['scl'], sh_data['shX'], sh_data['shY'], sh_data['shZ'], inds['ogn'], inds['nor'], inds['pos'], inds['rot'], inds['scl'], sh_inds['shX'], sh_inds['shY'], sh_inds['shZ'], cp_inds)

                vals_list.append(vals)
                rows_list.append(rows)
                cols_list.append(cols)
            return vals_list, rows_list, cols_list

        elif name == 'cp_bez' or name == 'cp_coons':
            vals_list, rows_list, cols_list = [], [], []
            for face in self.faces.values():
                num_u = face._num_cp_total['u']
                num_v = face._num_cp_total['v']

                vals = numpy.ones(3 * num_u * num_v)
                inds = face.vec_inds['cp_prim'].flatten()

                vals_list.append(vals)
                rows_list.append(inds)
                cols_list.append(inds)
            return vals_list, rows_list, cols_list
                
