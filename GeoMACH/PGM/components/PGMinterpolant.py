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

    pass
