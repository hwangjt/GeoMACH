#/bin/sh
LIB_CGNS=#-lcgns 
CGNS_INCLUDE=#/usr/local/include

SRC_FILES = PUBS/src/bspline/basis.f90 \
	PUBS/src/bspline/knotopen.f90 \
	PUBS/src/bspline/paramuni.f90 \
	PUBS/src/tensor/curve.f90 \
	PUBS/src/tensor/surface.f90 \
	PUBS/src/patchwork/initializeTopology.f90 \
	PUBS/src/patchwork/initializeBsplines.f90 \
	PUBS/src/patchwork/initializePoints.f90 \
	PUBS/src/patchwork/computeIndices.f90 \
	PUBS/src/patchwork/computeKnots.f90 \
	PUBS/src/patchwork/computeParameters.f90 \
        PUBS/src/patchwork/computeJacobian.f90 \
        PUBS/src/patchwork/computeDOFmapping.f90 \
	PUBS/src/patchwork/evaluatePoints.f90 \
	PUBS/src/patchwork/evaluateProjections.f90 \
	PUBS/src/patchwork/getSurface.f90

default:	
	f2py --fcompiler=gnu95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PUBSlib ${LIB_CGNS}
	-mv PUBSlib.so ./PUBS/

g95:	
	f2py --fcompiler=g95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PUBSlib ${LIB_CGNS}
	-mv PUBSlib.so ./PUBS/

gfortran:
	f2py --fcompiler=gnu95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PUBSlib ${LIB_CGNS}
	-mv PUBSlib.so ./PUBS/

intel:
	f2py --fcompiler=intel -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PUBSlib ${LIB_CGNS}
	-mv PUBSlib.so ./PUBS/

clean:
	-rm PUBS/PUBSlib.so
