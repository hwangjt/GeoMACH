"""
GeoMACH PGM object class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse



class PGMobject(object):

    def __init__(self):
        self._shape = None
        self.vec_data = {}
        self.vec_inds = {}
        self.value = 0.0

    def get_vec_shape(self, name):
        return self._shape

    def initialize_vec(self, name, data, inds):
        self.vec_data[name] = data
        self.vec_inds[name] = inds
        self.data = data
        self.inds = inds
