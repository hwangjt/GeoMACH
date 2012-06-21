from __future__ import division
import PUBS
import numpy, pylab
import copy, time
import mpl_toolkits.mplot3d.axes3d as p3

n = [20,40]

P0 = []
P0.append(numpy.zeros((n[1],n[1],3),order='F'))
P0.append(numpy.zeros((n[0],n[1],3),order='F'))
P0.append(numpy.zeros((n[0],n[1],3),order='F'))
P0.append(numpy.zeros((n[1],n[0],3),order='F'))
P0.append(numpy.zeros((n[1],n[0],3),order='F'))

r = 1
dx = 1
dy = 1
dz = 1
            
a1 = [ 1]
a2 = [-1]
b1 = [ 1]
b2 = [-1]
s  = [-1]
for k in range(1):
    for i in range(n[1]):
        for j in range(n[1]):
            a = a1[k] + (a2[k]-a1[k])*i/(n[1]-1)
            b = b1[k] + (b2[k]-b1[k])*j/(n[1]-1)
            y = s[k]*r/(1+a**2+b**2)**0.5
            P0[0][i,j,0] = a*y*dx
            P0[0][i,j,1] = y*dy
            P0[0][i,j,2] = b*y*dz
a1 = [ 0, -1]
a2 = [ 1,  0]
b1 = [ 1, -1]
b2 = [-1,  1]
s  = [-1,  1]
for k in range(2):
    for i in range(n[0]):
        for j in range(n[1]):
            a = a1[k] + (a2[k]-a1[k])*i/(n[0]-1)
            b = b1[k] + (b2[k]-b1[k])*j/(n[1]-1)
            x = s[k]*r/(1+a**2+b**2)**0.5
            P0[k+1][i,j,0] = x*dx
            P0[k+1][i,j,1] = a*x*dy
            P0[k+1][i,j,2] = b*x*dz
a1 = [-1,  1]
a2 = [ 1, -1]
b1 = [-1,  0]
b2 = [ 0,  1]
s  = [ 1, -1]
for k in range(2):
    for i in range(n[1]):
        for j in range(n[0]):
            a = a1[k] + (a2[k]-a1[k])*i/(n[1]-1)
            b = b1[k] + (b2[k]-b1[k])*j/(n[0]-1)
            z = s[k]*r/(1+a**2+b**2)**0.5
            P0[k+3][i,j,0] = a*z*dx
            P0[k+3][i,j,1] = b*z*dy
            P0[k+3][i,j,2] = z*dz

oml1 = PUBS.PUBS()
oml1.importSurfaces(P0)

oml1.computePoints()

t0 = time.time()
for k in range(100):
    P,s,u,v = oml1.computeProjection(numpy.array([0.5,-1,0.1],order='F'))
print time.time()-t0
print 'Projection test:'
print P[0,0],s[0,0],u[0,0],v[0,0]

h=1e-5
print '1st parametric derivative test:'
print (oml1.computePt(0,0.1,0.4+h) - oml1.computePt(0,0.1,0.4-h))/2/h
print oml1.computePt(0,0.1,0.4,0,1)
print '2nd parametric derivative test:'
print (oml1.computePt(0,0.1,0.4+2*h) - 2*oml1.computePt(0,0.1,0.4+h) + oml1.computePt(0,0.1,0.4))/h**2
print oml1.computePt(0,0.1,0.4,0,2)

for i in range(oml1.ngroup):
    oml1.group_n[i] -= 4
oml1.updateEvaluation()
oml1.plot(pylab.figure(),False)
pylab.show()
