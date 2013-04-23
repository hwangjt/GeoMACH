from __future__ import division
import numpy, pylab, time, copy
import scipy.sparse

from GeoMACH.PAM import PAMlib


class Parameter(object):

    def __init__(self, var, shp, shp0, P, T, Tdim, D, Ddim, B, Bdim):

        mu,mv = shp
        zeros = numpy.zeros
        self.var = var
        self.shp0 = shp0
        self.P = zeros((shp[0],shp[1],5), order='F')
        self.T = [numpy.linspace(0,1,mu), numpy.linspace(0,1,mv)]

        if P != None:
            self.setP(P)
        if T != None:
            self.setT(T, Tdim)
        if D != None:
            self.setD(D, Ddim)
        if B != None:
            self.setB(B, Bdim)

    def set(self, A, B, dim):
        B = numpy.array(B,order='F')
        if len(B.shape)==1:
            if dim==0:
                for j in range(A.shape[1]):
                    A[:,j] = B
            else:
                for i in range(A.shape[0]):
                    A[i,:] = B
        else:
            A[:,:] = B.reshape(A.shape)

    def setP(self, P):
        self.P[:,:,0] = numpy.array(P,order='F').reshape(self.P.shape[:2])

    def setT(self, T, dim=0):
        self.T[dim][:] = numpy.array(T,order='F').reshape(self.T[dim].shape)

    def setD(self, D, dim=0):
        self.set(self.P[:,:,1+dim], D, dim)

    def setB(self, B, dim=0):
        self.set(self.P[:,:,3+dim], B, dim)

    def compute(self):
        mu,mv = self.P.shape[:2]
        nu,nv = self.shp0
        P = self.P
        Tu,Tv = self.T
        return PAMlib.computeparameter(mu, mv, nu, nv, P, Tu, Tv)
