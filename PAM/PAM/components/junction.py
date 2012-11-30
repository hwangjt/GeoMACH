from __future__ import division
from PAM.components import Component
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib



class Junction(Component):

    def __init__(self, mComp, mSide, fComp, fFace, NW, SE):
        super(Junction,self).__init__() 

        self.n = 10
        self.mComp = mComp
        self.mSide = mSide
        self.fComp = fComp
        self.fFace = fFace
        self.NW = NW
        self.SE = SE
        self.initializeIndices()
        self.initializePoints()
        self.removeSurfaces()

        self.ms = []
        self.ms.append(numpy.zeros(self.nu,int))
        self.ms.append(numpy.zeros(self.nv,int))
        self.ms.append(None)
        self.faces.append([1,2])

    def initializeIndices(self):
        NW = self.NW
        SE = self.SE
        uInc = SE[0] > NW[0]
        vInc = SE[1] > NW[1]
        if uInc and vInc:
            self.rotate1 = lambda P: P
            self.rotate2 = lambda P: P
        elif (not uInc) and (not vInc):
            self.rotate1 = lambda P: P[::-1,::-1,:]
            self.rotate2 = lambda P: P[::-1,::-1,:]
        elif (not uInc) and vInc:
            self.rotate1 = lambda P: numpy.swapaxes(P,0,1)[::-1,:,:]
            self.rotate2 = lambda P: numpy.swapaxes(P,0,1)[:,::-1,:]
        elif uInc and (not vInc):
            self.rotate1 = lambda P: numpy.swapaxes(P,0,1)[:,::-1,:]
            self.rotate2 = lambda P: numpy.swapaxes(P,0,1)[::-1,:,:]
        u1 = min(SE[0],NW[0])
        u2 = max(SE[0],NW[0])
        v1 = min(SE[1],NW[1])
        v2 = max(SE[1],NW[1])
        self.fK = self.rotate2(self.fComp.Ks[self.fFace][u1:u2+1,v1:v2+1])
        self.nu = self.fK.shape[0]
        self.nv = self.fK.shape[1]

    def initializePoints(self):
        def getEdge(f, i, j, u=0, v=0):
            if u==-1:
                i = -1-i
            if v==-1:
                j = -1-j
            surf = mKs[f][i,j]
            if u==1:
                return mPs[surf][:,j,:]
            elif u==-1:
                return mPs[surf][::-1,j,:]
            elif v==1:
                return mPs[surf][i,:,:]
            elif v==-1:
                return mPs[surf][i,::-1,:]

        rot = self.rotate1
        n = self.n
        nu = self.nu
        nv = self.nv
        fK = self.fK
        fPs = self.fComp.Ps
        mPs = self.mComp.Ps
        mKs = self.mComp.Ks

        self.Ps = []
        self.Ks = [-numpy.ones((nu,nv),int)]
        counter = 0
        for j in range(nv):
            for i in range(nu):
                if i==0 or j==0 or i==nu-1 or j==nv-1:
                    self.Ps.append(numpy.zeros((n,n,3),order='F'))
                    self.Ks[0][i,j] = counter
                    counter += 1

        iK = self.Ks[0]
        for j in range(1,nv-1):
            i = 0
            e1 = rot(fPs[fK[i,j]])[0,:,:]
            if nu > 2:
                e2 = getEdge(0, -1, j-1, v=-1)
            elif self.mSide==0:
                e2 = getEdge(0, j-1, 0, u=-1)
            elif self.mSide==1:
                e2 = getEdge(0, j-1, -1, u=1)
            self.Ps[iK[i,j]][:,:,:] = PAMlib.createinterface(n, n, e1, e2)
            i = -1
            if nu > 2:
                e1 = getEdge(1, -1, j-1, v=1)
            elif self.mSide==0:
                e1 = getEdge(1, j-1, 0, u=1)
            elif self.mSide==1:
                e1 = getEdge(1, j-1, -1, u=-1)
            e2 = rot(fPs[fK[i,j]])[-1,:,:]
            self.Ps[iK[i,j]][:,:,:] = PAMlib.createinterface(n, n, e1, e2)
        for i in range(1,nu-1):
            j = 0
            e1 = rot(fPs[fK[i,j]])[:,0,:]
            e2 = getEdge(4, -1, i-1, v=1)
            self.Ps[iK[i,j]][:,:,:] = PAMlib.createinterface(n, n, e1, e2)
            j = -1
            e1 = getEdge(2, 0, i-1, v=1)
            e2 = rot(fPs[fK[i,j]])[:,-1,:]
            self.Ps[iK[i,j]][:,:,:] = PAMlib.createinterface(n, n, e1, e2)

        for i in [0,-1]:
            e1 = rot(fPs[fK[i,0]])[:,0,:]
            e2 = self.Ps[iK[i,1]][:,0,:]
            self.Ps[iK[i,0]][:,:,:] = PAMlib.createinterface(n, n, e1, e2)

            e1 = self.Ps[iK[i,-2]][:,-1,:]
            e2 = rot(fPs[fK[i,-1]])[:,-1,:]
            self.Ps[iK[i,-1]][:,:,:] = PAMlib.createinterface(n, n, e1, e2)

    def initializeIndices2(self, NW, SE):
        uInc = SE[0] > NW[0]
        vInc = SE[1] > NW[1]
        if uInc and vInc:
            E = [ 0, 1]
            N = [-1, 0]
            self.rotate = lambda P: P
        elif (not uInc) and (not vInc):
            E = [ 0,-1]
            N = [ 1, 0]
            self.rotate = lambda P: P[::-1,::-1,:]
        elif (not uInc) and vInc:
            E = [-1, 0]
            N = [ 0,-1]
            self.rotate = lambda P: numpy.swapaxes(P,0,1)[::-1,:,:]
        elif uInc and (not vInc):
            E = [ 1, 0]
            N = [ 0, 1]
            self.rotate = lambda P: numpy.swapaxes(P,0,1)[:,::-1,:]
        self.N = numpy.array(N,int)
        self.E = numpy.array(E,int)
        self.NW = numpy.array(NW,int)
        self.SE = numpy.array(SE,int)
        self.SW = self.NW - self.N*(self.SE-self.NW)
        self.NE = self.SE + self.N*(self.SE-self.NW)
        self.nu = int(sum(abs(self.NW-self.SW))) + 1
        self.nv = int(sum(abs(self.NW-self.NE))) + 1   

    def removeSurfaces(self):
        fPs = self.fComp.Ps
        fKs = self.fComp.Ks
        f1 = self.fFace
        NW = self.NW
        SE = self.SE
        umin = min(NW[0],SE[0])
        umax = max(NW[0],SE[0])
        vmin = min(NW[1],SE[1])
        vmax = max(NW[1],SE[1])
        for j1 in range(vmax,vmin-1,-1):
            for i1 in range(umax,umin-1,-1):
                fPs.pop(fKs[f1][i1,j1])
                for f2 in range(len(fKs)):
                    for j2 in range(fKs[f2].shape[1]):
                        for i2 in range(fKs[f2].shape[0]):
                            if fKs[f2][i2,j2] > fKs[f1][i1,j1]:
                                fKs[f2][i2,j2] -= 1
                fKs[f1][i1,j1] = -1

    def setDOFs(self):
        self.setC1('surf', 0, val=True)

    def initializeVariables(self):
        Ns = self.Ns

    def propagateQs(self):
        shortenEdge = lambda Q, ni, i: Q[sum(ni[:i]):sum(ni[:i+1])+1,:]
        def writeEdge(i, j, d):
            if d==0:
                vCurves[i][j] = shortenEdge(fQ[:,j,:], ni, i)
            else:
                hCurves[i][j] = shortenEdge(fQ[i,:,:], nj, j)
        def getCorner(i, j, d):
            if d==0:
                return fQ[sum(ni[:i]),j,:]
            else:
                return fQ[i,sum(nj[:j]),:]
        def getEdge(f, i, j, d, u=0, v=0):
            if u==1:
                P = mQs[f][:,j,:]
            elif u==-1:
                P = mQs[f][::-1,j,:]
            elif v==1:
                P = mQs[f][i,:,:]
            elif v==-1:
                P = mQs[f][i,::-1,:]
            if d==0:
                return shortenEdge(P, ni, i)
            else:
                return shortenEdge(P, nj, j)

        n = self.n
        nu = self.nu
        nv = self.nv

        hCurves = [[0 for v in range(nv)] for u in range(nu+1)]
        vCurves = [[0 for v in range(nv+1)] for u in range(nu)]
        corners = numpy.zeros((nu+1,nv+1,2,3))

        fQ = self.rotate(self.fComp.Qs[self.fFace])
        mQs = self.mComp.Qs
        ni = self.getms(0,0)
        nj = self.getms(0,1)
        
        for j in range(nv):
            writeEdge( 0, j, 1)
            writeEdge(-1, j, 1)
        for i in range(nu):
            writeEdge( i, 0, 0)
            writeEdge( i,-1, 0)
            
        for j in range(nv+1):
            corners[0,j,0,:] = getCorner(0,j,1)
            corners[0,j,1,:] = getCorner(1,j,1)
            corners[-1,j,0,:] = getCorner(-1,j,1)
            corners[-1,j,1,:] = getCorner(-2,j,1)
        for i in range(nu+1):
            corners[i,0,0,:] = getCorner(i,0,0)
            corners[i,0,1,:] = getCorner(i,1,0)
            corners[i,-1,0,:] = getCorner(i,-1,0)
            corners[i,-1,1,:] = getCorner(i,-2,0)

        if nu > 2:
            for j in range(nv-2):
                hCurves[1][1+j] = getEdge(0, -1, j, 1, v=-1)
                hCurves[-2][1+j] = getEdge(1, -1, j, 1, v=1)
            for i in range(nu-2):
                vCurves[1+i][1] = getEdge(4, -1, i, 0, v=1)
                vCurves[1+i][-2] = getEdge(2, 0, i, 0, v=1)       
        elif self.mSide==0:
            for j in range(nv-2):
                hCurves[1][1+j] = getEdge(0, j, 0, 1, u=-1)
                hCurves[-2][1+j] = getEdge(1, j, 0, 1, u=1)
        elif self.mSide==1:
            for j in range(nv-2):
                hCurves[1][1+j] = getEdge(0, j, -1, 1, u=1)
                hCurves[-2][1+j] = getEdge(1, j, -1, 1, u=-1)  

        Z = numpy.zeros(3)
        for j in range(nv-1):
            if j==0:
                P1 = hCurves[0][j+1][0,:]
                P2 = hCurves[1][j+1][0,:]
            else:
                P1 = hCurves[0][j][-1,:]
                P2 = hCurves[1][j][-1,:]
            vCurves[0][1+j] = PAMlib.beziercurve(vCurves[0][0].shape[0], P1, P2, Z, Z)
            if j==0:
                P1 = hCurves[-2][j+1][0,:]
                P2 = hCurves[-1][j+1][0,:]
            else:
                P1 = hCurves[-2][j][-1,:]
                P2 = hCurves[-1][j][-1,:]
            vCurves[-1][1+j] = PAMlib.beziercurve(vCurves[-1][0].shape[0], P1, P2, Z, Z)
        for i in range(nu-1):
            if i==0:
                P1 = vCurves[i+1][0][0,:]
                P2 = vCurves[i+1][1][0,:]
            else:
                P1 = vCurves[i][0][-1,:]
                P2 = vCurves[i][1][-1,:]
            hCurves[1+i][0] = PAMlib.beziercurve(hCurves[0][0].shape[0], P1, P2, Z, Z)
            if i==0:
                P1 = vCurves[i+1][-2][0,:]
                P2 = vCurves[i+1][-1][0,:]
            else:
                P1 = vCurves[i][-2][-1,:]
                P2 = vCurves[i][-1][-1,:]
            hCurves[1+i][-1] = PAMlib.beziercurve(hCurves[0][-1].shape[0], P1, P2, Z, Z) 

        for j in range(nv):
            for i in range(nu):
                if i==0 or j==0 or i==nu-1 or j==nv-1:
                    P = PAMlib.coonspatch(vCurves[i][j].shape[0], hCurves[i][j].shape[0], \
                                              hCurves[i][j][:,:], hCurves[i+1][j][:,:], \
                                              vCurves[i][j][:,:], vCurves[i][j+1][:,:])
                    self.Qs[0][sum(ni[:i]):sum(ni[:i+1])+1,sum(nj[:j]):sum(nj[:j+1])+1,:] = P
