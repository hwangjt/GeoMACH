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

PGMsources = [
    'src/PGM/component.f90',
    'src/PGM/parameter.f90',
    'src/PGM/primitive/computeAngles.f90',
    'src/PGM/primitive/computeRotations.f90',
    'src/PGM/primitive/computeSections.f90',
    'src/PGM/primitive/computeShape.f90',
    'src/PGM/interpolant/computeCone.f90',
    'src/PGM/interpolant/computeJunction.f90',
    'src/PGM/interpolant/coonsPatch.f90',
    'src/PGM/interpolant/interpolateGrid.f90',
    ]

PSMsources = [
    'src/PSM/computeMemberEdges.f90',
    'src/PSM/computeProjtnInputs.f90',
    'src/PSM/computePreviewSurfaces.f90',
    'src/PSM/computeEdgeLengths.f90',
    'src/PSM/computeFaceDimensions.f90',
    'src/PSM/importMembers.f90',
    'src/PSM/computePreviewMembers.f90',
    'src/PSM/computeMemberTopology.f90',
    'src/PSM/computeAdjoiningEdges.f90',
    'src/PSM/computeFaceEdges.f90',
    'src/PSM/computeGroupIntersections.f90',
    'src/PSM/computeIntersectionVerts.f90',
    'src/PSM/computeSurfaces.f90',
    'src/PSM/computeSurfaceProjections.f90',
    ]

QUADsources = [
    'src/QUAD/computeAdjMap.f90',
    'src/QUAD/computeDivisions.f90',
    'src/QUAD/computeGrid.f90',
    'src/QUAD/computeIntersections.f90',
    'src/QUAD/computeQuads.f90',
    'src/QUAD/computeTri2Quad.f90',
    'src/QUAD/computeTriangles.f90',
    'src/QUAD/deleteDuplicateEdges.f90',
    'src/QUAD/deleteDuplicateQuads.f90',
    'src/QUAD/deleteDuplicateTriangles.f90',
    'src/QUAD/deleteDuplicateVerts.f90',
    'src/QUAD/importEdges.f90',
    'src/QUAD/reorderCollinear.f90',
    'src/QUAD/splitEdges.f90',
    'src/QUAD/trianglesToEdges.f90',
    ]

CDTsources = [
    'src/CDT/addNode.f90',
    'src/CDT/computeCDT.f90',
    'src/CDT/constraints.f90',
    'src/CDT/delaunay.f90',
    'src/CDT/delete.f90',
    'src/CDT/misc.f90',
    'src/CDT/nearest.f90',
    'src/CDT/output.f90',
    'src/CDT/postProcess.f90',
    ]

entry_points = """
[openmdao.parametric_geometry]
GeoMACH.PGM.configurations.conventional.Conventional = GeoMACH.PGM.configurations.conventional:Conventional

[openmdao.binpub]
GeoMACH.PGM.configurations.configuration.GeoMACHSender = GeoMACH.PGM.configurations.configuration:GeoMACHSender
"""

config = Configuration(name='GeoMACH')
config.add_extension('PUBS.PUBSlib', sources=PUBSsources)
config.add_extension('PGM.PGMlib', sources=PGMsources)
config.add_extension('PSM.PSMlib', sources=PSMsources)
config.add_extension('PSM.QUADlib', sources=QUADsources)
config.add_extension('PSM.CDTlib', sources=CDTsources)

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
                'sphinx_build/html/_modules/GeoMACH/PGM/components/*.html',
                'sphinx_build/html/_modules/GeoMACH/PGM/configurations/*.html',
                'sphinx_build/html/_modules/GeoMACH/PUBS/*.html',
            ]
        },
        'entry_points': entry_points
        }
kwds.update(config.todict())

setup(**kwds)
