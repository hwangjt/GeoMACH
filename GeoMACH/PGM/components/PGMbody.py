"""
GeoMACH body class
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


class PGMbody(PGMprimitive):
    """ Body component """

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
        super(PGMbody, self).__init__()

        self._num_surf['x'] = num_x
        self._num_surf['y'] = num_y
        self._num_surf['z'] = num_z

        self._ax1 = 3
        self._ax2 = 1

        self.faces['rgt'] = PGMface(num_y, num_x)
        self.faces['top'] = PGMface(num_z, num_x)
        self.faces['lft'] = PGMface(num_y, num_x)
        self.faces['bot'] = PGMface(num_z, num_x)

    def set_diff(self):
        for face in self.faces.values():
            face.set_diff_surf(True)

    def compute(self, name):
        if name == 'cp_prim':
            theta1 = {'rgt': -1/4.0, 
                      'top': 1/4.0,
                      'lft': 3/4.0,
                      'bot': 5/4.0,
            }
            theta2 = {'rgt': 1/4.0, 
                      'top': 3/4.0,
                      'lft': 5/4.0,
                      'bot': 7/4.0,
            }

            flt = self.props['flt'].vec_data['prop']
            for fname in self.faces:
                face = self.faces[fname]
                num_u = face._num_cp_total['u']
                num_v = face._num_cp_total['v']
                self._shapes[fname][:,:,:] \
                    = PGMlib.computeshape(num_u, num_v, 
                                          theta1[fname], theta2[fname],
                                          numpy.ones((num_v, 3),
                                                     order='F'), 
                                          flt, numpy.zeros((num_u, num_v),
                                                           order='F'))

        return super(PGMbody, self).compute(name)
