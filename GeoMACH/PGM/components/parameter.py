from __future__ import division
import numpy

from GeoMACH.PGM import PGMlib


class Parameter(object):

    def __init__(self, shp, shp0, P, T, Tdim, D, Ddim, B, Bdim):
        self.mu, self.mv = shp
        self.nu, self.nv = shp0
        self.P = numpy.zeros((shp[0],shp[1],3), order='C')
        self.B = numpy.zeros((shp[0],shp[1],2), order='C', dtype=bool)
        self.T = [numpy.linspace(0,1,self.mu), numpy.linspace(0,1,self.mv)]

        if P is not None:
            self.setP(P)
        if T is not None:
            self.setT(T, Tdim)
        if D is not None:
            self.setD(D, Ddim)
        if B is not None:
            self.setB(B, Bdim)

        self.param_vec = None
        self.param_ind = None

    def initialize_parameters(self, param_vec, param_ind):
        self.param_vec = param_vec.reshape((self.mu, self.mv, 3), order='C')
        self.param_ind = param_ind.reshape((self.mu, self.mv, 3), order='C')

        self.param_vec[:,:,:] = self.P[:,:,:3]

    def set(self, A, B, dim):
        B = numpy.array(B,order='C')
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
        self.P[:,:,0] = numpy.array(P,order='C').reshape(self.P.shape[:2])

    def setT(self, T, dim=0):
        self.T[dim][:] = numpy.array(T,order='C').reshape(self.T[dim].shape)

    def setD(self, D, dim=0):
        self.set(self.P[:,:,1+dim], D, dim)

    def setB(self, B, dim=0):
        self.set(self.B[:,:,dim], B, dim)

    def compute_jacobian(self, prop_ind):
        Tu,Tv = self.T
        nD = self.nu * self.nv * 12
        Da, Di, Dj = PGMlib.computeprops(nD, self.mu, self.mv, self.nu, self.nv, 
                                         self.B, Tu, Tv, self.param_ind, prop_ind)
        return Da, Di, Dj
