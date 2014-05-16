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
    'src/PGM/interpolant/computeTip.f90',
    'src/PGM/interpolant/coonsPatch.f90',
    'src/PGM/interpolant/interpolant.f90',
    'src/PGM/interpolant/misc.f90',
    ]

PSMsources = [
    'src/PSM/GFEM/computeProjtnInputs.f90',
    'src/PSM/GFEM/computePreviewSurfaces.f90',
    'src/PSM/GFEM/computeEdgeLengths.f90',
    'src/PSM/GFEM/computeFaceDimensions.f90',
    'src/PSM/GFEM/importMembers.f90',
    'src/PSM/GFEM/computePreviewMembers.f90',
    'src/PSM/GFEM/computeMemberTopology.f90',
    'src/PSM/GFEM/computeAdjoiningEdges.f90',
    'src/PSM/GFEM/computeFaceEdges.f90',
    'src/PSM/GFEM/computeGroupIntersections.f90',
    'src/PSM/GFEM/computeGroupSplits.f90',
    'src/PSM/GFEM/computeIntersectionVerts.f90',
    'src/PSM/GFEM/computeSurfaces.f90',
    'src/PSM/GFEM/computeSurfaceProjections.f90',
    'src/PSM/GFEM/computeMemberEdges.f90',
    'src/PSM/GFEM/computeMemberNodes.f90',
    'src/PSM/GFEM/computeMembers.f90',
    'src/PSM/GFEM/removeDuplicateNodes.f90',
    'src/PSM/GFEM/removeRightQuads.f90',
    'src/PSM/GFEM/identifySymmNodes.f90',
    ]

QUADsources = [
    'src/PSM/QUAD/importEdges.f90',
    'src/PSM/QUAD/reorderCollinear.f90',
    'src/PSM/QUAD/addIntersectionPts.f90',
    'src/PSM/QUAD/addEdgePts.f90',
    'src/PSM/QUAD/addInteriorPts.f90',
    'src/PSM/QUAD/splitEdges.f90',
    'src/PSM/QUAD/removeDegenerateEdges.f90',
    'src/PSM/QUAD/removeDuplicateEdges.f90',
    'src/PSM/QUAD/removeDuplicateQuads.f90',
    'src/PSM/QUAD/removeDuplicateTriangles.f90',
    'src/PSM/QUAD/removeDuplicateVerts.f90',
    'src/PSM/QUAD/computeAdjMap.f90',
    'src/PSM/QUAD/computeTriangles.f90',
    'src/PSM/QUAD/computeQuads.f90',
    'src/PSM/QUAD/computeConstraints.f90',
    'src/PSM/QUAD/computeQuadDominant.f90',
    'src/PSM/QUAD/splitTrisNQuads.f90',
    'src/PSM/QUAD/computeQuad2Edge.f90',
    'src/PSM/QUAD/removeInvalidQuads.f90',
    ]

CDTsources = [
    'src/PSM/CDT/addNode.f90',
    'src/PSM/CDT/computeCDT.f90',
    'src/PSM/CDT/constraints.f90',
    'src/PSM/CDT/delaunay.f90',
    'src/PSM/CDT/delete.f90',
    'src/PSM/CDT/misc.f90',
    'src/PSM/CDT/nearest.f90',
    'src/PSM/CDT/output.f90',
    'src/PSM/CDT/postProcess.f90',
    ]

BLSsources = [
    'src/PSM/BLS/assembleMtx.f90',
    ]

entry_points = """
[openmdao.parametric_geometry]
GeoMACH.PGM.configurations.conventional.Conventional = GeoMACH.PGM.configurations.conventional:Conventional

[openmdao.binpub]
GeoMACH.PGM.configurations.configuration.GeoMACHSender = GeoMACH.PGM.configurations.configuration:GeoMACHSender
"""

addLib = lambda name, srcs: config.add_extension(name, sources=srcs, extra_compile_args=['-fbounds-check'])#, f2py_options=['--debug-capi'])

config = Configuration(name='GeoMACH')
addLib('PUBS.PUBSlib', PUBSsources)
addLib('PGM.PGMlib', PGMsources)
addLib('PSM.PSMlib', PSMsources)
addLib('PSM.QUADlib', QUADsources)
addLib('PSM.CDTlib', CDTsources)
addLib('PSM.BLSlib', BLSsources)

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
