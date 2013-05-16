import setuptools
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

PUBSsources = [
    'src/PUBS/bspline/basis.f90',
    'src/PUBS/bspline/knotopen.f90',
    'src/PUBS/bspline/paramuni.f90',
    'src/PUBS/tensor/curve.f90',
    'src/PUBS/tensor/surface.f90',
    'src/PUBS/patchwork/initializeTopology.f90',
    'src/PUBS/patchwork/initializeBsplines.f90',
    'src/PUBS/patchwork/initializePoints.f90',
    'src/PUBS/patchwork/computeIndices.f90',
    'src/PUBS/patchwork/computeKnots.f90',
    'src/PUBS/patchwork/computeParameters.f90',
    'src/PUBS/patchwork/computeJacobian.f90',
    'src/PUBS/patchwork/computeDOFmapping.f90', 
    'src/PUBS/patchwork/evaluatePoints.f90',
    'src/PUBS/patchwork/evaluateProjections.f90',
    'src/PUBS/patchwork/getSurface.f90',
    ]

PAMsources = [
    'src/PAM/layout/importEdges.f90',
    'src/PAM/layout/computeIntersections.f90',
    'src/PAM/layout/deleteDuplicates.f90',
    'src/PAM/layout/addConnectors.f90',
    'src/PAM/layout/computePolygons.f90',
    'src/PAM/layout/splitPolygons.f90',
    'src/PAM/layout/extractSurface.f90',
    'src/PAM/component/component.f90',
    'src/PAM/component/parameter.f90',
    'src/PAM/component/primitive/computeAngles.f90',
    'src/PAM/component/primitive/computeRotations.f90',
    'src/PAM/component/primitive/computeSections.f90',
    'src/PAM/component/primitive/computeShape.f90',
    'src/PAM/component/interpolant/computeCone.f90',
    'src/PAM/component/interpolant/computeJunction.f90',
    'src/PAM/component/interpolant/coonsPatch.f90',
    'src/PAM/component/interpolant/interpolateGrid.f90',
    'src/PAM/component/structures/computeInternalStructure.f90',
    'src/PAM/component/structures/shapeTypes.f90',
    ]

entry_points = """
[openmdao.parametric_geometry]
GeoMACH.PAM.configurations.conventional.Conventional = GeoMACH.PAM.configurations.conventional:Conventional

[openmdao.binpub]
GeoMACH.PAM.configurations.configuration.GeoMACHSender = GeoMACH.PAM.configurations.configuration:GeoMACHSender
"""

config = Configuration(name='GeoMACH')
config.add_extension('PUBS.PUBSlib', sources=PUBSsources)
config.add_extension('PAM.PAMlib', sources=PAMsources)

kwds = {'install_requires':['numpy','scipy'],
        'version': '0.1',
        'zip_safe': False,
        'license': 'LGPL',
        'include_package_data': True,
        'package_dir': {'': '.'},
        'packages': setuptools.find_packages('.'),
        'package_data': {
            'GeoMACH': [
                'sphinx_build/html/*.html',
                'sphinx_build/html/*.js',
                'sphinx_build/html/*.inv',
                'sphinx_build/html/_static/*',
                'sphinx_build/html/_sources/*.txt',
                'sphinx_build/html/_modules/index.html',
                'sphinx_build/html/_modules/GeoMACH/PAM/components/*.html',
                'sphinx_build/html/_modules/GeoMACH/PAM/configurations/*.html',
                'sphinx_build/html/_modules/GeoMACH/PUBS/*.html',
            ]
        },
        'entry_points': entry_points
        }
kwds.update(config.todict())

setup(**kwds)
