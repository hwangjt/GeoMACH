from __future__ import division
import numpy
import scipy.sparse

from GeoMACH.PGM import PGMlib
from GeoMACH.PUBS import PUBS


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

    def count_parameters(self):
        return 3 * self.mu * self.mv

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


class Parameter2(object):

    def __init__(self, shp_param, shp_prop, P, Tu, Tv, ku, kv):
        self.mu, self.mv = shp_param
        self.ku, self.kv = ku, kv

        nu, nv = shp_prop
        self.urange = [i for i in range(nu) if Tu[0] <= i/(nu-1) <= Tu[1]]
        self.vrange = [j for j in range(nv) if Tv[0] <= j/(nv-1) <= Tv[1]]
        self.nu, self.nv = len(self.urange), len(self.vrange)

        self.param0 = numpy.zeros((self.mu, self.mv), order='F')
        if P is not None:
            self.set(P)

        self.param_vec = None
        self.param_ind = None

    def count_parameters(self):
        return self.mu * self.mv

    def set(self, P):
        self.param0[:,:] = numpy.array(P,order='F').reshape((self.mu,self.mv), order='F')

    def initialize_parameters(self, param_vec, param_ind):
        self.param_vec = param_vec.reshape((self.mu, self.mv), order='F')
        self.param_ind = param_ind.reshape((self.mu, self.mv), order='F')

        self.param_vec[:,:] = self.param0[:,:]

    def compute_jacobian(self, prop_ind):
        ku0, kv0 = max(2,self.ku), max(2,self.kv)
        mu0, mv0 = max(2,self.mu), max(2,self.mv)

        oml = PUBS.PUBS([numpy.zeros((self.nu,self.nv,3), order='F')])
        oml.edgeProperty(0, 0, 0, ku0)
        oml.edgeProperty(0, 0, 1, kv0)
        oml.edgeProperty(0, 1, 0, mu0)
        oml.edgeProperty(0, 1, 1, mv0)
        oml.updateBsplines(True)

        nC = mu0 * mv0
        Ca = numpy.ones(nC)
        Ci = numpy.zeros(nC, int)
        Cj = numpy.array(numpy.linspace(0,nC-1,nC), int)
        indices = Cj.reshape((mu0, mv0), order='F')
        for indv in xrange(mv0):
            for indu in xrange(mu0):
                Ci[indices[indu,indv]] = oml.getIndex(0, indu, indv, 1)
        Cmtx = scipy.sparse.csr_matrix((Ca, (Ci, Cj)), shape=(nC, nC))

        nP = self.nu * self.nv
        Pa = numpy.ones(nP)
        Pi = numpy.array(numpy.linspace(0,nP-1,nP), int)
        Pj = numpy.zeros(nP, int)
        indices = Pi.reshape((self.nu, self.nv), order='F')
        for indv in xrange(self.nv):
            for indu in xrange(self.nu):
                Pj[indices[indu,indv]] = oml.getIndex(0, indu, indv, 0)
        Pmtx = scipy.sparse.csr_matrix((Pa, (Pi, Pj)), shape=(nP, nP))

        indices = Cj.reshape((mu0, mv0), order='F')
        Cmap = numpy.zeros((nC,2), int)
        for indv in xrange(mv0):
            for indu in xrange(mu0):
                Cmap[indices[indu,indv],:] = [indu if self.ku!=1 else 0, 
                                              indv if self.kv!=1 else 0]

        indices = Pi.reshape((self.nu, self.nv), order='F')
        Pmap = numpy.zeros((nP,2), int)
        for indv in xrange(self.nv):
            for indu in xrange(self.nu):
                Pmap[indices[indu,indv],:] = [self.urange[indu], 
                                              self.vrange[indv]]

        jac = Pmtx.dot(oml.J.dot(Cmtx)).tocoo(True)
        Da = numpy.array(jac.data)
        Di = numpy.array(prop_ind[Pmap[jac.row,0], Pmap[jac.row,1]])
        Dj = numpy.array(self.param_ind[Cmap[jac.col,0], Cmap[jac.col,1]])
        return Da, Di, Dj
        
