import setuptools
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

sources = [
    'PAM/src/layout/importEdges.f90',
    'PAM/src/layout/computeIntersections.f90',
    'PAM/src/layout/deleteDuplicates.f90',
    'PAM/src/layout/addConnectors.f90',
    'PAM/src/layout/computePolygons.f90',
    'PAM/src/layout/splitPolygons.f90',
    'PAM/src/layout/extractSurface.f90',
    'PAM/src/component/initialize.f90',
    ]

config = Configuration(name='PAM')
config.add_extension('PAMlib', sources=sources)

kwds = {'install_requires':['numpy','scipy'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'LGPL',
        }
kwds.update(config.todict())

setup(**kwds)
