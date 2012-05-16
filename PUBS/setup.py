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
    'PUBS/src/bspline/basis.f90',
    'PUBS/src/bspline/knotopen.f90',
    'PUBS/src/bspline/paramuni.f90',
    'PUBS/src/tensor/curve.f90',
    'PUBS/src/tensor/surface.f90',
    'PUBS/src/fileIO/readcgnsvol.f90',
    'PUBS/src/fileIO/readcgnssurf.f90',
    'PUBS/src/patchwork/initializeTopology.f90',
    'PUBS/src/patchwork/initializeBsplines.f90',
    'PUBS/src/patchwork/initializePoints.f90',
    'PUBS/src/patchwork/computeIndices.f90',
    'PUBS/src/patchwork/computeKnots.f90',
    'PUBS/src/patchwork/computeParameters.f90',
    'PUBS/src/patchwork/computeJacobian.f90',
    'PUBS/src/patchwork/computeDOFmapping.f90', 
    'PUBS/src/patchwork/computeDerivative.f90', 
    'PUBS/src/patchwork/computeProjection.f90', 
    'PUBS/src/patchwork/plotSurfaces.f90',
    'PUBS/src/PUBS.py'
    ]

config = Configuration(name='PUBS')
config.add_extension('PUBSlib', sources=sources, include_dirs=include_dirs, library_dirs=library_dirs, libraries=libraries)

kwds = {'install_requires':['numpy','scipy','pylab'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'public domain',
        }
kwds.update(config.todict())

setup(**kwds)
