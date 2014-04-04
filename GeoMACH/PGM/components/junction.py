from __future__ import division
import numpy, time

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Interpolant
        


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
        faces = self.faces

        faces[0].setC1('surf', val=True)

        si = self.si
        sj = self.sj
        for j in range(sj[1],sj[2]):
            faces[0].setC1('surf', i=si[1]-1, j=j, u=-1, val=False)
            faces[0].setC1('surf', i=si[2], j=j, u=0, val=False)
        for i in range(si[1],si[2]):
            faces[0].setC1('surf', i=i, j=sj[1]-1, v=-1, val=False)
            faces[0].setC1('surf', i=i, j=sj[2], v=0, val=False)

        faces[0].setC1('surf', i=si[1]-1, j=sj[1]-1, u=-1, v=-1, val=False)
        faces[0].setC1('surf', i=si[2], j=sj[1]-1, u=0, v=-1, val=False)
        faces[0].setC1('surf', i=si[1]-1, j=sj[2], u=-1, v=0, val=False)
        faces[0].setC1('surf', i=si[2], j=sj[2], u=0, v=0, val=False)

        self.removeSurfaces()

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
                self.ni = [1,0,1]
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

        self.fK = self.rotate(self.fComp.faces[self.fFace].surf_indices)[self.fNW[0]:self.fNW[0]+sum(self.ni),self.fNW[1]:self.fNW[1]+sum(self.nj)]

    def initializeVerts(self):
        vtx = lambda i, j, u, v: self.rotate(self.fComp.Ps[self.fK[i,j]])[u,v,:]
        vtxM = lambda f, i, j, u, v: self.mComp.Ps[self.mComp.faces[f].surf_indices[i,j]][u,v,:]

        verts = numpy.zeros((4,4,3),order='F')
        for i in [0,-1]:
            for j in range(3):
                verts[i,j,:] = vtx(i, self.sj[j], i, 0)
                verts[i,j+1,:] = vtx(i, self.sj[j+1]-1, i, -1)
        for j in [0,-1]:
            for i in range(3):
                verts[i,j,:] = vtx(self.si[i], j, 0, j)
                verts[i+1,j,:] = vtx(self.si[i+1]-1, j, -1, j)
        if self.mSide==-1:
            verts[1,1,:] = vtxM(2, -1, 0, -1, 0)
            verts[2,1,:] = vtxM(2, -1, -1, -1, -1)
            verts[1,2,:] = vtxM(0, 0, 0, 0, 0)
            verts[2,2,:] = vtxM(0, 0, -1, 0, -1)
        else:
            if self.mSide==0:
                L = vtxM(0, -1, 0, -1, 0)
                R = vtxM(0, 0, 0, 0, 0)
            else:
                L = vtxM(0, 0, -1, 0, -1)
                R = vtxM(0, -1, -1, -1, -1)
            verts[1,1,:] = L
            verts[2,1,:] = L
            verts[1,2,:] = R
            verts[2,2,:] = R

        return verts

    def initializeSurfaces(self):
        def get(P, u, v, d):
            edge = P[u,:,:] if not u==None else P[:,v,:]
            return edge if d==1 else edge[::-1,:]

        getM = lambda f, i, j: self.mComp.Ps[self.mComp.faces[f].surf_indices[i,j]]
        getI = lambda i, j: self.Ps[self.faces[0].surf_indices[i,j]]
        getF = lambda i, j: self.rotate(self.fComp.Ps[self.fK[i,j]])

        def copy(iI, jI, fM, iM, jM, uI=None, vI=None, uM=None, vM=None, d=1):
            edgeI = get(getI(iI, jI), uI, vI, 1)
            edgeM = get(getM(fM, iM, jM), uM, vM, d)
            edgeI[:,:] = edgeM[:,:]

        def copy2(i, j, u=None, v=None):
            edgeI = get(getI(i, j), u, v, 1)
            edgeF = get(getF(i, j), u, v, 1)
            edgeI[:,:] = edgeF[:,:]

        verts = self.initializeVerts()
        
        nP = self.nP
        ni = self.ni
        nj = self.nj
        si = self.si
        sj = self.sj

        self.Ps = []
        face = self.faces[0]
        face.surf_indices[:,:] = -1
        counter = 0
        for j in range(sj[3]):
            for i in range(si[3]):
                if i<si[1] or j<sj[1] or i>=si[2] or j>=sj[2]:
                    self.Ps.append(numpy.zeros((nP,nP,3),order='F'))
                    face.surf_indices[i,j] = counter
                    counter += 1

        for b in range(3):
            for a in range(3):
                for j in range(nj[b]):
                    jj = sj[b] + j
                    for i in range(ni[a]):
                        ii = si[a] + i
                        if not face.surf_indices[ii,jj]==-1:
                            self.Ps[face.surf_indices[ii,jj]][:,:,:] = PGMlib.bilinearinterp(nP, ni[a], nj[b], i+1, j+1, verts[a:a+2,b:b+2,:])

        for i in [0,-1]:
            for j in range(face.num_surf[1]):
                copy2(i, j, u=i)
        for j in [0,-1]:
            for i in range(face.num_surf[0]):
                copy2(i, j, v=j)

        if not self.mSide==-1:
            mComp = self.mComp
            ii = si[1] - 1
            for j in range(nj[1]):
                jj = sj[1] + j
                if self.mSide==0:
                    copy(ii, jj, 0, -1-j, 0, uI=-1, vM=0, d=-1)
                else:
                    copy(ii, jj, 0, j, -1, uI=-1, vM=-1, d=1)
            ii = si[2]
            for j in range(nj[1]):
                jj = sj[1] + j
                if self.mSide==0:
                    copy(ii, jj, 1, j, 0, uI=0, vM=0, d=1)
                else:
                    copy(ii, jj, 1, -1-j, -1, uI=0, vM=-1, d=-1)

    def removeSurfaces(self):
        fK0 = self.fK

        for j in range(fK0.shape[1]):
            for i in range(fK0.shape[0]):
                self.oml0.visible[fK0[i,j]] = False

    def initializeDOFmappings(self):
        super(Junction,self).initializeDOFmappings()
        mu, mv = self.faces[0].num_cp_list[:]
        nu = range(4)
        nv = range(4)
        for k in range(4):
            nu[k] = sum(mu[:self.si[k]])
            nv[k] = sum(mv[:self.sj[k]])
        N = self.faces[0].index_array
        N[nu[1],nv[1]:nv[2]+1] = -1
        N[nu[2],nv[1]:nv[2]+1] = -1
        N[nu[1]:nu[2]+1,nv[1]] = -1
        N[nu[1]:nu[2]+1,nv[2]] = -1
        
    def computeQs(self):
        fu, fv = self.fComp.faces[self.fFace].num_cp_list[:]
        fu,fv = self.flip(fu,fv)
        fu1 = sum(fu[:self.fNW[0]])
        fu2 = sum(fu[:self.fNW[0]+self.si[3]])
        fv1 = sum(fv[:self.fNW[1]])
        fv2 = sum(fv[:self.fNW[1]+self.sj[3]])
        fQ = self.rotate(self.fComp.faces[self.fFace].cp_array)[fu1:fu2+1,fv1:fv2+1,:]

        getEdge = self.getEdge
        if self.mSide==-1:
            W = getEdge(self.mComp.faces[2].cp_array, i=-1, d=1)
            E = getEdge(self.mComp.faces[0].cp_array, i=0, d=1)
            N = getEdge(self.mComp0.faces[0].cp_array, i=-1, d=-1)
            S = getEdge(self.mComp1.faces[0].cp_array, i=-1, d=1)
        elif self.mSide==0:
            W = numpy.zeros((1,2,3),order='F')
            E = numpy.zeros((1,2,3),order='F')
            N = getEdge(self.mComp.faces[0].cp_array, j=0, d=-1)
            S = getEdge(self.mComp.faces[1].cp_array, j=0, d=1)
        elif self.mSide==1:
            W = numpy.zeros((1,2,3),order='F')
            E = numpy.zeros((1,2,3),order='F')
            N = getEdge(self.mComp.faces[0].cp_array, j=-1, d=1)
            S = getEdge(self.mComp.faces[1].cp_array, j=-1, d=-1)

        mu, mv = self.faces[0].num_cp_list[:]
        nu = range(3)
        nv = range(3)
        for k in range(3):
            nu[k] = sum(mu[self.si[k]:self.si[k+1]]) + 1
            nv[k] = sum(mv[self.sj[k]:self.sj[k+1]]) + 1

        v = self.variables
        self.faces[0].cp_array = PGMlib.computejunction(sum(nu)-2, sum(nv)-2, nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], v['fC1'], v['mC1'], W, E, N, S, fQ, v['shp'])
