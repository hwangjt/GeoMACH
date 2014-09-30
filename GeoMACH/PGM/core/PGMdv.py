"""
GeoMACH design variable class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy

from GeoMACH.PGM.core.PGMobject import PGMobject


class PGMdv(PGMobject):

    def __init__(self, shape, value=0.0, lower=None, upper=None, scale=None):
        super(PGMdv, self).__init__()

        self._shape = shape
        self.value = value
        self.lower = lower
        self.upper = upper
        self.scale = scale
