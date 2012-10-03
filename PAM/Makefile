#/bin/sh
LIB_CGNS=#-lcgns 
CGNS_INCLUDE=#/usr/local/include

SRC_FILES = PAM/src/layout/addConnectors.f90 \
	PAM/src/layout/computeIntersections.f90 \
	PAM/src/layout/computePolygons.f90 \
	PAM/src/layout/deleteDuplicates.f90 \
	PAM/src/layout/extractSurface.f90 \
	PAM/src/layout/importEdges.f90 \
	PAM/src/layout/splitPolygons.f90 \
	PAM/src/component/initialize.f90 \
	PAM/src/component/updateQs.f90 \
	PAM/src/component/computeInternalStructure.f90 \
	PAM/src/component/shapeTypes.f90 \
	PAM/src/component/bodysections.f90 \
	PAM/src/component/wingsections.f90 \
	PAM/src/component/interpolation.f90

default:	
	f2py --fcompiler=intelem -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PAMlib ${LIB_CGNS}
	-mv PAMlib.so ./PAM/

g95:	
	f2py --fcompiler=g95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PAMlib ${LIB_CGNS}
	-mv PAMlib.so ./PAM/

gfortran:
	f2py --fcompiler=gnu95 -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PAMlib ${LIB_CGNS}
#	f2py --fcompiler=gnu95 --f90flags="-fbounds-check" -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PAMlib ${LIB_CGNS}
	-mv PAMlib.so ./PAM/

intel:
	f2py --fcompiler=intel -I$(CGNS_INCLUDE) -c ${SRC_FILES} -m PAMlib ${LIB_CGNS}
	-mv PAMlib.so ./PAM/

clean:
	-rm PAM/PAMlib.so
