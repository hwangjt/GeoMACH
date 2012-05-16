import os.path
#import setuptools
import sys

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

CGNS_INCLUDE='/usr/local/include'
LIB_CGNS='cgns'

include_dirs = [CGNS_INCLUDE]
library_dirs = []
libraries = [LIB_CGNS]
sources = [
    'bspline2d/src/bspline/basis.f90',
    'bspline2d/src/bspline/knotopen.f90',
    'bspline2d/src/bspline/paramuni.f90',
    'bspline2d/src/tensor/curve.f90',
    'bspline2d/src/tensor/surface.f90',
    'bspline2d/src/fileIO/readcgnsvol.f90',
    'bspline2d/src/fileIO/readcgnssurf.f90',
    'bspline2d/src/patchwork/initializeTopology.f90',
    'bspline2d/src/patchwork/initializeBsplines.f90',
    'bspline2d/src/patchwork/initializePoints.f90',
    'bspline2d/src/patchwork/computeIndices.f90',
    'bspline2d/src/patchwork/computeKnots.f90',
    'bspline2d/src/patchwork/computeParameters.f90',
    'bspline2d/src/patchwork/computeJacobian.f90',
    'bspline2d/src/patchwork/computeDOFmapping.f90', 
    'bspline2d/src/patchwork/computeDerivative.f90', 
    'bspline2d/src/patchwork/computeProjection.f90', 
    'bspline2d/src/patchwork/plotSurfaces.f90',
    'bspline2d/src/bspline2d.py'
    ]

config = Configuration(name='bspline2d')
config.add_extension('bspline2d_lib', sources=sources, include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries)

kwds = {'install_requires':['numpy','scipy','pylab'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'public domain',
        }
kwds.update(config.todict())

setup(**kwds)
