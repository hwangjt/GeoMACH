import setuptools
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

sources = [
    'PUBS/src/bspline/basis.f90',
    'PUBS/src/bspline/knotopen.f90',
    'PUBS/src/bspline/paramuni.f90',
    'PUBS/src/tensor/curve.f90',
    'PUBS/src/tensor/surface.f90',
    'PUBS/src/patchwork/initializeTopology.f90',
    'PUBS/src/patchwork/initializeBsplines.f90',
    'PUBS/src/patchwork/initializePoints.f90',
    'PUBS/src/patchwork/computeIndices.f90',
    'PUBS/src/patchwork/computeKnots.f90',
    'PUBS/src/patchwork/computeParameters.f90',
    'PUBS/src/patchwork/computeJacobian.f90',
    'PUBS/src/patchwork/computeDOFmapping.f90', 
    'PUBS/src/patchwork/evaluatePoints.f90',
    'PUBS/src/patchwork/evaluateProjections.f90',
    'PUBS/src/patchwork/getSurface.f90',
    ]

config = Configuration(name='PUBS')
config.add_extension('PUBSlib', sources=sources)

kwds = {'install_requires':['numpy','scipy'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'LGPL',
        }
kwds.update(config.todict())

setup(**kwds)
