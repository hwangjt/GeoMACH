from __future__ import division
import sys
import PUBS
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3
from mayavi import mlab


n = [40,40]
P0 = []
P0.append(numpy.zeros((n[0],n[1],3)))
P0.append(numpy.zeros((n[0],n[1],3)))
P0.append(numpy.zeros((n[0],n[1],3)))
P0.append(numpy.zeros((n[0],n[1],3)))
P0.append(numpy.zeros((n[0],n[1],3)))

for i in range(n[0]):
    for j in range(n[1]):
        P0[0][i,j,0] = -1
        P0[0][i,j,1] = -i/(n[0]-1)
        P0[0][i,j,2] = -1+2*j/(n[1]-1)
        P0[1][i,j,0] = 1
        P0[1][i,j,1] = -i/(n[0]-1)
        P0[1][i,j,2] = -1+2*j/(n[1]-1)
        P0[2][i,j,0] = -1+2*i/(n[0]-1)
        P0[2][i,j,1] = -1
        P0[2][i,j,2] = -1+2*j/(n[1]-1)
        P0[3][i,j,0] = 1-2*i/(n[0]-1)
        P0[3][i,j,1] = -j/(n[1]-1)
        P0[3][i,j,2] = 1
        P0[4][i,j,0] = -1+2*i/(n[0]-1)
        P0[4][i,j,1] = -j/(n[1]-1)
        P0[4][i,j,2] = -1

oml1 = PUBS.PUBS()
oml1.importSurfaces(P0)

oml1.C[10,0] -= 0.3
oml1.computePoints()
oml1.edge_c1[0,0] = True
oml1.surf_c1[-1,:] = True
oml1.edge_c1[3,1] = True
oml1.updateBsplines()

#oml1.plotm(mlab.figure(),False)
#mlab.show()

oml1.plot(pylab.figure(),False)
pylab.show()
