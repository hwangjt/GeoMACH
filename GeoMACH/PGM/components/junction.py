from __future__ import division
import numpy, time

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Interpolant, Property
        


class Junction(Interpolant):

    def __init__(self, fComp, fFace, fRot, fNW, mComp, mComp0=None, mComp1=None, mSide=-1, ni=None, nj=None):
        super(Junction,self).__init__() 

        self.fComp = fComp
        self.fFace = fFace
        self.fRot = fRot
        self.fNW = fNW
        self.mComp = mComp
        self.mSide = mSide
        self.mComp0 = mComp0
        self.mComp1 = mComp1
        self.initializeIndices(ni, nj)
        self.initializeFaces()
        self.initializeSurfaces()

    def setDOFs(self):
        face = self.faces['def']

        face.setC1('surf', val=True)

        si = self.si
        sj = self.sj
        for j in range(sj[1],sj[2]):
            face.setC1('surf', i=si[1]-1, j=j, u=-1, val=False)
            face.setC1('surf', i=si[2], j=j, u=0, val=False)
        if self.mSide==-1:
            for i in range(si[1],si[2]):
                face.setC1('surf', i=i, j=sj[1]-1, v=-1, val=False)
                face.setC1('surf', i=i, j=sj[2], v=0, val=False)

        face.setC1('surf', i=si[1]-1, j=sj[1]-1, u=-1, v=-1, val=False)
        face.setC1('surf', i=si[2], j=sj[1]-1, u=0, v=-1, val=False)
        face.setC1('surf', i=si[1]-1, j=sj[2], u=-1, v=0, val=False)
        face.setC1('surf', i=si[2], j=sj[2], u=0, v=0, val=False)

        self.removeHiddenSurfaces()

    def initializeIndices(self, ni, nj):
        if not ni==None:
            self.ni = numpy.array(ni,int)
        if not nj==None:
            self.nj = numpy.array(nj,int)
        if self.mSide==-1:
            if ni==None:
                self.ni = [1,self.mComp.ms[0].shape[0],1]
            if nj==None:
                self.nj = [1,self.mComp.ms[2].shape[0],1]
        else:
            if ni==None:
                self.ni = [1,1,1]
            if nj==None:
                self.nj = [1,self.mComp.ms[0].shape[0],1]

        if self.fRot==0:
            self.rotate = lambda P: P
            self.flip = lambda nu, nv: [nu,nv]
        elif self.fRot==1:
            self.rotate = lambda P: numpy.swapaxes(P,0,1)[::-1,:]
            self.flip = lambda nu, nv: [nv[::-1],nu]
        elif self.fRot==2:
            self.rotate = lambda P: P[::-1,::-1]
            self.flip = lambda nu, nv: [nu[::-1],nv[::-1]]
        elif self.fRot==3:
            self.rotate = lambda P: numpy.swapaxes(P,0,1)[:,::-1]
            self.flip = lambda nu, nv: [nv,nu[::-1]]

    def initializeSurfaces(self):        
        nP = self.nP
        si = self.si
        sj = self.sj

        face = self.faces['def']
        face.surf_indices[:,:] = -1
        counter = 0
        for j in range(sj[3]):
            for i in range(si[3]):
                if i<si[1] or i>=si[2] or ((j<sj[1] or j>=sj[2]) and self.mSide==-1):
                    face.surf_indices[i,j] = counter
                    counter += 1

        return

    def removeHiddenSurfaces(self):
        fK = self.rotate(self.fComp.faces[self.fFace].surf_indices)[self.fNW[0]:self.fNW[0]+sum(self.ni),self.fNW[1]:self.fNW[1]+sum(self.nj)]
        if self.mSide != -1:
            fK = numpy.vstack((fK[0,:], -numpy.ones(fK.shape[1], int), fK[1,:]))

        for j in range(fK.shape[1]):
            for i in range(fK.shape[0]):
                surf = fK[i,j]
                if surf >= 0:
                    self.oml0.visible[surf] = False

    def compute_cp_wireframe(self):
        fu, fv = self.fComp.faces[self.fFace].num_cp_list[:]
        fu,fv = self.flip(fu,fv)
        fu1 = sum(fu[:self.fNW[0]])
        fu2 = sum(fu[:self.fNW[0]+self.si[3]])
        if self.mSide!=-1:
            fu2 = sum(fu[:self.fNW[0]+2])
        fv1 = sum(fv[:self.fNW[1]])
        fv2 = sum(fv[:self.fNW[1]+self.sj[3]])
        fFace_inds = self.rotate(self.fComp.faces[self.fFace].cp_indices)[fu1:fu2+1,fv1:fv2+1]
        
        #getEdge = self.getEdge
        if self.mSide==-1:
            #W = getEdge(self.mComp.faces['lft'].cp_indices, i=-1, d=1)
            #E = getEdge(self.mComp.faces['rgt'].cp_indices, i=0, d=1)
            #N = getEdge(self.mComp0.faces['def'].cp_indices, i=-1, d=-1)
            #S = getEdge(self.mComp1.faces['def'].cp_indices, i=-1, d=1)
            pass
        elif self.mSide==0:
            W = numpy.zeros((4,2,3),order='F')
            E = numpy.zeros((4,2,3),order='F')
            N = self.mComp.faces['upp'].cp_indices[::-1,:2,:]
            S = self.mComp.faces['low'].cp_indices[:,:2,:]
        elif self.mSide==1:
            W = numpy.zeros((4,2,3),order='F')
            E = numpy.zeros((4,2,3),order='F')
            N = self.mComp.faces['upp'].cp_indices[:,-1:-3:-1,:]
            S = self.mComp.faces['low'].cp_indices[::-1,-1:-3:-1,:]

        mu, mv = self.faces['def'].num_cp_list[:]
        nu = range(3)
        nv = range(3)
        for k in range(3):
            nu[k] = sum(mu[self.si[k]:self.si[k+1]]) + 1
            nv[k] = sum(mv[self.sj[k]:self.sj[k+1]]) + 1
        nu0 = sum(nu)-2
        nv0 = sum(nv)-2
        fInds = -numpy.ones((nu0, nv0, 3), dtype=int, order='F')
        fInds[:nu[0],:] = fFace_inds[:nu[0],:]
        fInds[-nu[2]:,:] = fFace_inds[-nu[2]:,:]

        nD = 3 * 2 * nv0 + 3 * 2 * (nu0 - 2)
        nD += 3 * 2
        if nu[1] != 1:
            nD += 3 * 2
        nD += 3 * 2 * (nv[1] - 2)
        if nu[1] != 1:
            nD += 3 * 2 * (nu[1] - 2)
        nD += 3 * 4 * 2 * (nu[0] - 2)
        nD += 3 * 4 * 2 * (nu[2] - 2)
        nD += 3 * 4 * (nv[0] - 2)
        nD += 3 * 4 * (nv[2] - 2)
        if nu[1] != 1:
            nD += 3 * 4 * (nv[0] - 2)
            nD += 3 * 4 * (nv[2] - 2)

        fC1 = self.props['fC1'].prop_vec
        mC1 = self.props['mC1'].prop_vec
        Da, Di, Dj = PGMlib.computejunctionwireframe(nD, nu0, nv0, nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], fC1, mC1, W, E, N, S, fInds, self.faces['def'].cp_indices)
        Da = Da * (-1 != Dj)
        Dj = numpy.maximum(0, Dj)
        return Da, Di, Dj

    def compute_cp_surfs(self):
        mu, mv = self.faces['def'].num_cp_list[:]
        nu = range(3)
        nv = range(3)
        for k in range(3):
            nu[k] = sum(mu[self.si[k]:self.si[k+1]]) + 1
            nv[k] = sum(mv[self.sj[k]:self.sj[k+1]]) + 1
        nu0 = sum(nu)-2
        nv0 = sum(nv)-2

        nD = 0
        for i in range(3):
            for j in range(3):
                if (nu[1] != 1) or (i !=1):
                    nD += 3 * 8 * (nu[i]-2) * (nv[j]-2)
        
        Da, Di, Dj = PGMlib.computejunctioncoons(nD, nu0, nv0, nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], self.faces['def'].cp_indices)
        return Da, Di, Dj
