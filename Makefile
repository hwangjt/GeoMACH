#/bin/sh
LIB_CGNS=-lcgns 
CGNS_INCLUDE=/usr/local/include

SRC_FILES = src/basis/basis.f90 \
	src/basis/knotopen.f90 \
	src/basis/paramuni.f90 \
	src/tensor/curve.f90 \
	src/tensor/surface.f90 \
	src/fileIO/readcgnsvol.f90 \
	src/fileIO/readcgnssurf.f90 \
	src/oml/initializeTopology.f90 \
	src/oml/initializeBsplines.f90 \
	src/oml/initializePoints.f90 \
	src/oml/computeIndices.f90 \
	src/oml/computeKnots.f90 \
	src/oml/computeParameters.f90 \
        src/oml/computeJacobian.f90 \
        src/oml/computeDOFmapping.f90 \
	src/oml/computeDerivative.f90 \
	src/oml/computeProjection.f90 \
	src/oml/plotSurfaces.f90

default:	
	f2py --fcompiler=intelem -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m GeoMACH ${LIB_CGNS}
	-mv GeoMACH.so ./python/

g95:	
	f2py --fcompiler=g95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m GeoMACH ${LIB_CGNS}
	-mv GeoMACH.so ./python/

gfortran:
	f2py --fcompiler=gnu95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m GeoMACH ${LIB_CGNS}
	-mv GeoMACH.so ./python/

intel:
	f2py --fcompiler=intel -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m GeoMACH ${LIB_CGNS}
	-mv GeoMACH.so ./python/

clean:
	-rm python/GeoMACH.so
