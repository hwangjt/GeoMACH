from __future__ import division
import numpy

from GeoMACH.PUBS import PUBS


n = [6,4]
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

oml1 = PUBS.PUBS(P0)
oml1.write2Tec('test')
oml1.plot()
