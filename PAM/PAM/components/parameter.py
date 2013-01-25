from __future__ import division
import numpy, pylab, time, copy
import scipy.sparse
import PAM.PAMlib as PAMlib


class Parameter(object):

    def __init__(self, var, shp, shp0, P, Tdim, T, Ddim, D, Bdim, B):

        mu,mv = shp
        zeros = numpy.zeros
        self.var = var
        self.shp0 = shp0
        self.P = zeros(shp, order='F')
        self.Ts = [numpy.linspace(0,1,mu), numpy.linspace(0,1,mv)]
        self.Ds = [zeros(shp,order='F'), zeros(shp,order='F')]
        self.Bs = [zeros(shp,dtype=bool,order='F'), zeros(shp,dtype=bool,order='F')]

        if P != None:
            self.setP(P)
        if T != None:
            self.setT(Tdim, T)
        if D != None:
            self.setD(Ddim, D)
        if B != None:
            self.setB(Bdim, B)

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
        self.P[:,:] = numpy.array(P,order='F').reshape(self.P.shape)

    def setT(self, dim, T):
        self.Ts[dim][:] = numpy.array(T,order='F').reshape(self.Ts[dim].shape)

    def setD(self, dim, D):
        self.set(self.Ds[dim], D, dim)

    def setB(self, dim, B):
        self.set(self.Bs[dim], B, dim)

    def compute(self):
        mu,mv = self.P.shape
        nu,nv = self.shp0
        P = self.P
        Tu,Tv = self.Ts
        Du,Dv = self.Ds
        Bu,Bv = self.Bs
        return PAMlib.computeparameter(mu, mv, nu, nv, P, Tu, Tv, Du, Dv, Bu, Bv)


