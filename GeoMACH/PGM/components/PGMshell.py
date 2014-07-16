"""
GeoMACH shell class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components.PGMprimitive import PGMprimitive
from GeoMACH.PGM.core.PGMface import PGMface
from GeoMACH.PGM.core.PGMproperty import PGMproperty


class PGMshell(PGMprimitive):
    """ Shell component """

    def __init__(self, num_x=1, num_y=1, num_z=1):
        """
        Parameters
        ----------
        num_x : ``float``
           Number of surfaces in the stream-wise direction
        num_z : ``float``
           Number of surfaces over the height
        num_z : ``float``
           Number of surfaces across the width
        """
        super(PGMshell, self).__init__()

        self._num_surf['x'] = num_x
        self._num_surf['y'] = num_y
        self._num_surf['z'] = num_z

        self._ax1 = 3
        self._ax2 = 1

        self.faces['rt0'] = PGMface(num_y, num_x)
        self.faces['tp0'] = PGMface(num_z, num_x)
        self.faces['lt0'] = PGMface(num_y, num_x)
        self.faces['bt0'] = PGMface(num_z, num_x)

        self.faces['rt1'] = PGMface(num_y, num_x)
        self.faces['tp1'] = PGMface(num_z, num_x)
        self.faces['lt1'] = PGMface(num_y, num_x)
        self.faces['bt1'] = PGMface(num_z, num_x)

    def set_diff(self):
        for face in self.faces.values():
            face.set_diff_surf(True)

    def initialize_props(self):
        super(PGMshell, self).initialize_props()

        props = self.props
        props['thk'] = PGMproperty()

    def assemble_sizes(self, bse):
        super(PGMshell, self).assemble_sizes(bse)

        num = self.faces.values()[0]._num_cp_total['v']
        props = self.props
        props['thk'].assemble_sizes(num, 3)

    def compute(self, name):
        if name == 'cp_prim':
            theta1 = {'rt0': -1/4.0,
                      'tp0': 1/4.0,
                      'lt0': 3/4.0,
                      'bt0': 5/4.0,
                      'rt1': 1/4.0,
                      'tp1': 3/4.0,
                      'lt1': 5/4.0,
                      'bt1': 7/4.0,
            }
            theta2 = {'rt0': 1/4.0,
                      'tp0': 3/4.0,
                      'lt0': 5/4.0,
                      'bt0': 7/4.0,
                      'rt1': -1/4.0,
                      'tp1': 1/4.0,
                      'lt1': 3/4.0,
                      'bt1': 5/4.0,
            }

            flt = self.props['flt'].vec_data['prop']
            thk = self.props['thk'].vec_data['prop']
            for fname in self.faces:
                face = self.faces[fname]
                num_u = face._num_cp_total['u']
                num_v = face._num_cp_total['v']
                if fname[2] == '0':
                    sgn = 1.0
                elif fname[2] == '1':
                    sgn = -1.0
                self._shapes[fname][:,:,:] \
                    = PGMlib.computeshape(num_u, num_v, 
                                          theta1[fname], theta2[fname],
                                          numpy.ones((num_v, 3),
                                                     order='F') + sgn*thk/2.0, 
                                          flt, numpy.zeros((num_u, num_v),
                                                           order='F'))
            output = super(PGMshell, self).compute(name)

            for fname in ['rt', 'tp', 'lt', 'bt']:
                for ind in range(2):
                    outer = self.faces[fname+'0'].vec_data['cp_prim'][:, -ind, :]
                    inner = self.faces[fname+'1'].vec_data['cp_prim'][::-1, -ind, :]
                    outer[:, :] = 0.5 * (outer + inner)
                    inner[:, :] = outer[:, :]

            return output

        else:
            return super(PGMshell, self).compute(name)
