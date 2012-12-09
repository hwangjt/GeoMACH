from __future__ import division
from PAM.components import Component
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib



class Junction(Component):

    def __init__(self, mComp, mSide, fComp, fFace, fNW, fSE):
        super(Junction,self).__init__() 

        self.nP = 10
        self.mComp = mComp
        self.mSide = mSide
        self.fComp = fComp
        self.fFace = fFace
        self.fNW = fNW
        self.fSE = fSE
        self.initializeIndices()
        self.initializePoints()
        self.removeSurfaces()

        self.ms = []
        self.ms.append(numpy.zeros(self.ni,int))
        self.ms.append(numpy.zeros(self.nj,int))
        self.ms.append(None)
        self.faces.append([1,2])

    def initializeIndices(self):
        fNW = self.fNW
        fSE = self.fSE
        uInc = fSE[0] > fNW[0]
        vInc = fSE[1] > fNW[1]
        if uInc and vInc:
            self.rotate = lambda P: P
        elif (not uInc) and (not vInc):
            self.rotate = lambda P: P[::-1,::-1]
        elif (not uInc) and vInc:
            self.rotate = lambda P: numpy.swapaxes(P,0,1)[:,::-1]
        elif uInc and (not vInc):
            self.rotate = lambda P: numpy.swapaxes(P,0,1)[::-1,:]
        self.fi1 = min(fSE[0],fNW[0])
        self.fi2 = max(fSE[0],fNW[0])
        self.fj1 = min(fSE[1],fNW[1])
        self.fj2 = max(fSE[1],fNW[1])
        self.fK = self.rotate(self.fComp.Ks[self.fFace][self.fi1:self.fi2+1,self.fj1:self.fj2+1])
        self.ni = self.fK.shape[0]
        self.nj = self.fK.shape[1]

    def initializePoints(self):
        def getmEdge(f, i, j, du=0, dv=0):
            if du==-1:
                i = -1-i
            if dv==-1:
                j = -1-j
            surf = mKs[f][i,j]
            if du==1:
                return mPs[surf][:,j,:]
            elif du==-1:
                return mPs[surf][::-1,j,:]
            elif dv==1:
                return mPs[surf][i,:,:]
            elif dv==-1:
                return mPs[surf][i,::-1,:]

        rot = self.rotate
        nP = self.nP
        ni = self.ni
        nj = self.nj
        fK = self.fK
        fPs = self.fComp.Ps
        mPs = self.mComp.Ps
        mKs = self.mComp.Ks

        self.Ps = []
        self.Ks = [-numpy.ones((ni,nj),int)]
        counter = 0
        for j in range(nj):
            for i in range(ni):
                if i==0 or j==0 or i==ni-1 or j==nj-1:
                    self.Ps.append(numpy.zeros((nP,nP,3),order='F'))
                    self.Ks[0][i,j] = counter
                    counter += 1

        iK = self.Ks[0]
        for j in range(1,nj-1):
            i = 0
            e1 = rot(fPs[fK[i,j]])[0,:,:]
            if ni > 2:
                e2 = getmEdge(0, -1, j-1, dv=-1)
            elif self.mSide==0:
                e2 = getmEdge(0, j-1, 0, du=-1)
            elif self.mSide==1:
                e2 = getmEdge(0, j-1, -1, du=1)
            self.Ps[iK[i,j]][:,:,:] = PAMlib.createinterface(nP, nP, e1, e2)
            i = -1
            if ni > 2:
                e1 = getmEdge(1, -1, j-1, dv=1)
            elif self.mSide==0:
                e1 = getmEdge(1, j-1, 0, du=1)
            elif self.mSide==1:
                e1 = getmEdge(1, j-1, -1, du=-1)
            e2 = rot(fPs[fK[i,j]])[-1,:,:]
            self.Ps[iK[i,j]][:,:,:] = PAMlib.createinterface(nP, nP, e1, e2)
        for i in range(1,ni-1):
            j = 0
            e1 = rot(fPs[fK[i,j]])[:,0,:]
            e2 = getmEdge(4, -1, i-1, dv=1)
            self.Ps[iK[i,j]][:,:,:] = numpy.swapaxes(PAMlib.createinterface(nP, nP, e1, e2),0,1)
            j = -1
            e1 = getmEdge(2, 0, i-1, dv=1)
            e2 = rot(fPs[fK[i,j]])[:,-1,:]
            self.Ps[iK[i,j]][:,:,:] = numpy.swapaxes(PAMlib.createinterface(nP, nP, e1, e2),0,1)

        for i in [0,-1]:
            e1 = rot(fPs[fK[i,0]])[:,0,:]
            e2 = self.Ps[iK[i,1]][:,0,:]
            self.Ps[iK[i,0]][:,:,:] = numpy.swapaxes(PAMlib.createinterface(nP, nP, e1, e2),0,1)

            e1 = self.Ps[iK[i,-2]][:,-1,:]
            e2 = rot(fPs[fK[i,-1]])[:,-1,:]
            self.Ps[iK[i,-1]][:,:,:] = numpy.swapaxes(PAMlib.createinterface(nP, nP, e1, e2),0,1)

    def removeSurfaces(self):
        fPs = self.fComp.Ps
        fKs = self.fComp.Ks
        f0 = self.fFace
        for j0 in range(self.fj2,self.fj1-1,-1):
            for i0 in range(self.fi2,self.fi1-1,-1):
                fPs.pop(fKs[f0][i0,j0])
                for f in range(len(fKs)):
                    for j in range(fKs[f].shape[1]):
                        for i in range(fKs[f].shape[0]):
                            if fKs[f][i,j] > fKs[f0][i0,j0]:
                                fKs[f][i,j] -= 1
                fKs[f0][i0,j0] = -1

    def setDOFs(self):
        self.setC1('surf', 0, val=True)

    def initializeVariables(self):
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        self.variables = {
            'f0':1.0,
            'm0':1.0,
            'shape':numpy.zeros((ni,nj),order='F')
            }

    def computeQs(self):
        getmPt = lambda f,i,j: mQs[f][i,j,:]
        def getmEdge(f, i=None, j=None, d=1):
            if j==None:
                P = mQs[f][i,:,:]
            else:
                P = mQs[f][:,j,:]
            if d==1:
                return P
            elif d==-1:
                return P[::-1,:]

        nfu = self.fComp.getms(self.fFace,0)
        nfv = self.fComp.getms(self.fFace,1)
        nfu1 = sum(nfu[:self.fi1])
        nfv1 = sum(nfv[:self.fj1])
        nfu2 = sum(nfu[:self.fi2+1])
        nfv2 = sum(nfv[:self.fj2+1])
        fQ = self.rotate(self.fComp.Qs[self.fFace][nfu1:nfu2+1,nfv1:nfv2+1,:])

        nu = self.getms(0,0)
        nv = self.getms(0,1)
        nu = numpy.array([nu[0]+1, sum(nu[1:-1])+1, nu[-1]+1])
        nv = numpy.array([nv[0]+1, sum(nv[1:-1])+1, nv[-1]+1])
        mQs = self.mComp.Qs
        mA = numpy.zeros((2,2,3),order='F')
        mB = numpy.zeros((2,2,3),order='F')

        if self.ni > 2:
            mQT = getmEdge(0, i=-1, d=-1)
            mQB = getmEdge(1, i=-1, d= 1)
            mQL = getmEdge(4, i=-1, d= 1)
            mQR = getmEdge(2, i= 0, d= 1)
            mA[0,0,:] = getmPt(0,-1,-1)
            mA[0,1,:] = getmPt(0,-1, 0)
            mA[1,0,:] = getmPt(1,-1, 0)
            mA[1,1,:] = getmPt(1,-1,-1)
            mB[0,0,:] = getmPt(0,-2,-1)
            mB[0,1,:] = getmPt(0,-2, 0)
            mB[1,0,:] = getmPt(1,-2, 0)
            mB[1,1,:] = getmPt(1,-2,-1)
        elif self.mSide==0:
            mQT = getmEdge(0, j= 0, d=-1)
            mQB = getmEdge(1, j= 0, d= 1)
            mQL = numpy.zeros((1,3),order='F')
            mQR = numpy.zeros((1,3),order='F')
            mA[0,0,:] = getmPt(1, 0, 0)
            mA[0,1,:] = getmPt(1,-1, 0)
            mA[1,0,:] = getmPt(1, 0, 0)
            mA[1,1,:] = getmPt(1,-1, 0)
            mB[0,0,:] = getmPt(1, 0, 1)
            mB[0,1,:] = getmPt(1,-1, 1)
            mB[1,0,:] = getmPt(1, 0, 1)
            mB[1,1,:] = getmPt(1,-1, 1)
        elif self.mSide==1:
            mQT = getmEdge(0, j=-1, d= 1)
            mQB = getmEdge(1, j=-1, d=-1)
            mQL = numpy.zeros((1,3),order='F')
            mQR = numpy.zeros((1,3),order='F')
            mA[0,0,:] = getmPt(1,-1,-1)
            mA[0,1,:] = getmPt(1, 0,-1)
            mA[1,0,:] = getmPt(1,-1,-1)
            mA[1,1,:] = getmPt(1, 0,-1)
            mB[0,0,:] = getmPt(1,-1,-2)
            mB[0,1,:] = getmPt(1, 0,-2)
            mB[1,0,:] = getmPt(1,-1,-2)
            mB[1,1,:] = getmPt(1, 0,-2)

        v = self.variables
        self.Qs[0] = PAMlib.computejunction(fQ.shape[0], fQ.shape[1], nu[0], nu[1], nu[2], nv[0], nv[1], nv[2], v['f0'], v['m0'], mQT, mQB, mQL, mQR, mA, mB, fQ, v['shape'])
