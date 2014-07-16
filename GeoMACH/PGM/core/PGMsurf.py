"""
GeoMACH surface class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy

from GeoMACH.PGM.core.PGMobject import PGMobject


class PGMsurf(PGMobject):

    def set_shape(self, num_u, num_v):
        self._shape = (num_u+1, num_v+1, 3)
