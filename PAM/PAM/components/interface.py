from __future__ import division
from PAM.components import Component
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3



class Interface(Component):

    def tup(self, v):
        return v[0], v[1]

    def mSpliceP(self, face, ind, rev=False):
        mPs = self.mComp.Ps
        mKs = self.mComp.Ks[face]
        if rev:
            face = 1-face
        mSide = self.mSide
        if mSide==face:
            return mPs[mKs[-1-ind,-mSide]][::-1,-mSide,:]
        else:
            return mPs[mKs[ind,-mSide]][:,-mSide,:]

    def fSpliceP(self, ind, u, v):
        E = self.E
        fPs = self.fComp.Ps
        fKs = self.fComp.Ks[self.fFace]
        P = fPs[fKs[self.tup(ind)]]
        nu = P.shape[0]
        nv = P.shape[1]
        if E[1]==1:
            u,v = u,v
            du = 1
            dv = 1
        elif E[1]==-1:
            u,v = u,v
            du = -1
            dv = -1
        elif E[0]==-1: 
            u,v = v,u
            du = -1
            dv = 1
        elif E[0]==1: 
            u,v = v,u
            du = 1
            dv = -1
            
        if v==None and dv==1:
            P = P
        elif v==None and dv==-1:
            P = P[:,::-1]
        elif dv==1:
            P = P[:,v]
        elif dv==-1:
            if v<0:
                v += nv
            P = P[:,nv-1 - v]

        if u==None and du==1:
            P = P
        elif u==None and du==-1:
            P = P[::-1]
        elif du==1:
            P = P[u]
        elif du==-1:
            if u<0:
                u += nu
            P = P[nu-1 - u]

        return P

    def mSpliceQ(self, face, ind, rev=False):
        mQs = self.mComp.Qs[face]
        mSide = self.mSide
        if rev:
            face = 1-face
        if mSide==face:
            return mQs[-1-ind,-mSide,:]
        else:
            return mQs[ind,-mSide,:]        

    def fSpliceQ(self, ref, i, j, u, v):
        fi = self.fComp.getni(self.fFace,0)
        fj = self.fComp.getni(self.fFace,1)    
        fQs = self.fComp.Qs[self.fFace]
        E = self.E
        if E[1]==1:
            u,v = u,v
            i,j = i,j
            i += 0
            j += 0
        elif E[1]==-1:
            u,v = -u,-v
            i,j = -i,-j
            i += 1
            j += 1
        elif E[0]==-1: 
            u,v = -v,u
            i,j = -j,i
            i += 1
            j += 0
        elif E[0]==1: 
            u,v = v,-u
            i,j = j,-i
            i += 0
            j += 1
        uu = sum(fi[:ref[0]+i])+u
        vv = sum(fj[:ref[1]+j])+v
        if uu == -1:
            uu = 0
        elif uu == fQs.shape[0]:
            uu = -1
        if vv == -1:
            vv = 0
        elif vv == fQs.shape[1]:
            vv = -1
        return fQs[uu,vv]

    def initializeParameters(self):
        Ns = self.Ns



class FullInterface(Interface):

    def __init__(self, comps, mComp, mSide, fComp, fFace, NW, SE):
        mComp = comps[mComp]
        fComp = comps[fComp]
        self.faces = numpy.zeros((1,2),int)
        self.faces[0,:] = [1,2]

        NW = numpy.array(NW)
        SE = numpy.array(SE)
        uInc = SE[0] > NW[0]
        vInc = SE[1] > NW[1]
        if uInc and vInc:
            E = numpy.array([ 0, 1],int)
            N = numpy.array([-1, 0],int)
        elif (not uInc) and (not vInc):
            E = numpy.array([ 0,-1],int)
            N = numpy.array([ 1, 0],int)
        elif (not uInc) and vInc:
            E = numpy.array([-1, 0],int)
            N = numpy.array([ 0,-1],int)
        elif uInc and (not vInc):
            E = numpy.array([ 1, 0],int)
            N = numpy.array([ 0, 1],int)
        SW = NW - N
        NE = SE + N
        
        self.mComp = mComp
        self.mSide = mSide
        self.fComp = fComp
        self.fFace = fFace
        self.NW = NW
        self.SE = SE
        self.SW = SW
        self.NE = NE
        self.E = E
        self.N = N

        Ps = []
        Ks = []

        fPs = fComp.Ps
        fKs = fComp.Ks[fFace]
        mPs = mComp.Ps
        mKs = mComp.Ks
        for i in range(mComp.Ks[0].shape[0]):
            n = fPs[fKs[self.tup(NW)]].shape[abs(N[1])]
            edge1 = self.fSpliceP(NW+(1+i)*E, 0, None)
            edge2 = self.mSpliceP(0,i)
            Ps.append(self.createInterface(n, edge1, edge2))

            n = fPs[fKs[self.tup(SW)]].shape[abs(N[1])]
            edge2 = self.fSpliceP(SW+(1+i)*E, -1, None)
            edge1 = self.mSpliceP(1,i)
            Ps.append(self.createInterface(n, edge1, edge2))

        n = fPs[fKs[self.tup(NW)]].shape[abs(N[0])]
        edge1 = self.fSpliceP(NW, None, 0)
        edge2 = Ps[0][:,0,:]
        Ps.append(self.createInterface(n, edge1, edge2, True))

        n = fPs[fKs[self.tup(SW)]].shape[abs(N[0])]
        edge1 = self.fSpliceP(SW, None, 0)
        edge2 = Ps[1][:,0,:]
        Ps.append(self.createInterface(n, edge1, edge2, True))

        n = fPs[fKs[self.tup(NE)]].shape[abs(N[0])]
        edge2 = self.fSpliceP(NE, None, -1)
        edge1 = Ps[2*mComp.Ks[0].shape[0]-2][:,-1,:]
        Ps.append(self.createInterface(n, edge1, edge2, True))

        n = fPs[fKs[self.tup(SE)]].shape[abs(N[0])]
        edge2 = self.fSpliceP(SE, None, -1)
        edge1 = Ps[2*mComp.Ks[0].shape[0]-1][:,-1,:]
        Ps.append(self.createInterface(n, edge1, edge2, True))

        Ps[-4][0,:,:] = self.fSpliceP(NW, 0, None)
        Ps[-3][-1,:,:] = self.fSpliceP(SW, -1, None)
        Ps[-2][0,:,:] = self.fSpliceP(NE, 0, None)
        Ps[-1][-1,:,:] = self.fSpliceP(SE, -1, None)

        K = numpy.zeros((2,2+mComp.Ks[0].shape[0]),int)
        counter = 0
        for j in range(mComp.Ks[0].shape[0]):
            for i in range(2):
                K[i,j+1] = counter
                counter += 1
        for j in range(2):
            for i in range(2):
                K[i,-j] = counter
                counter += 1
        Ks.append(K)

        u1 = min(NW[0],SE[0])
        u2 = max(NW[0],SE[0])
        v1 = min(NW[1],SE[1])
        v2 = max(NW[1],SE[1])
        for j in range(v2,v1-1,-1):
            for i in range(u2,u1-1,-1):
                fPs.pop(fKs[i,j])
                for f in range(len(fComp.Ks)):
                    for v in range(fComp.Ks[f].shape[1]):
                        for u in range(fComp.Ks[f].shape[0]):
                            if fComp.Ks[f][u,v] > fKs[i,j]:
                                fComp.Ks[f][u,v] -= 1
                fKs[i,j] = -1
            
        self.Ps = Ps
        self.Ks = Ks

        self.oml0 = []

    def setDOFs(self):
        oml0 = self.oml0
        Ks = self.Ks
        self.setSurfC1(0, val=True)        
        for j in range(1,Ks[0].shape[1]-1):
            for i in range(2):
                oml0.surf_c1[Ks[0][i,j],i-1,:] = False
        oml0.surf_c1[Ks[0][0,0],-1,-1] = False
        oml0.surf_c1[Ks[0][1,0],0,-1] = False
        oml0.surf_c1[Ks[0][0,-1],-1,0] = False
        oml0.surf_c1[Ks[0][1,-1],0,0] = False

    def isExteriorDOF(self, f, uType, vType, i, j):
        return False

    def propagateQs(self):
        si = self.getni(0,0)
        sj = self.getni(0,1)

        Ns = self.Ns
        Qs = self.Qs
        Ks = self.Ks
        fQs = self.fComp.Qs[self.fFace]    
        Qs[0][:,:,:] = 0
        for j in range(int(Ns[0].shape[1]-sj[0]-sj[-1])):
            Q1 = self.fSpliceQ(self.NW, 0, 1, -1, j)
            Q2 = self.mSpliceQ(0, j)
            Qs[0][:si[0]+1,j+sj[0],:] += numpy.outer(numpy.linspace(1,0,si[0]+2)[1:],Q1)
            Qs[0][:si[0]+1,j+sj[0],:] += numpy.outer(numpy.linspace(0,1,si[0]+2)[1:],Q2)

            Q2 = self.fSpliceQ(self.NW, 2, 1, 1, j)
            Q1 = self.mSpliceQ(1, j)
            Qs[0][si[0]:,j+sj[0],:] += numpy.outer(numpy.linspace(1,0,si[1]+2)[:-1],Q1)
            Qs[0][si[0]:,j+sj[0],:] += numpy.outer(numpy.linspace(0,1,si[1]+2)[:-1],Q2) 
        for i in range(si[0]):
            Q1 = self.fSpliceQ(self.NW, 0, 0, i, -1)
            Q2 = Qs[0][i,sj[0]+1]
            Qs[0][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(1,0,sj[0]+3)[1:-1],Q1)
            Qs[0][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(0,1,sj[0]+3)[1:-1],Q2)

            Q2 = self.fSpliceQ(self.NW, 0, Ks[0].shape[1], i, 1)
            Q1 = Qs[0][i,sum(sj[:-1])-1]
            Qs[0][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(1,0,sj[-1]+3)[1:-1],Q1)
            Qs[0][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(0,1,sj[-1]+3)[1:-1],Q2)
        for i in range(si[0],sum(si)):
            Q1 = self.fSpliceQ(self.NW, 1, 0, i-si[0], -1)
            Q2 = Qs[0][i,sj[0]+1]
            Qs[0][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(1,0,sj[0]+3)[1:-1],Q1)
            Qs[0][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(0,1,sj[0]+3)[1:-1],Q2)

            Q2 = self.fSpliceQ(self.NW, 1, Ks[0].shape[1], i-si[0], -1)
            Q1 = Qs[0][i,sum(sj[:-1])-1]
            Qs[0][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(1,0,sj[-1]+3)[1:-1],Q1)
            Qs[0][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(0,1,sj[-1]+3)[1:-1],Q2)



class HalfInterface(Interface):

    def __init__(self, comps, mComp, mSide, fComp, fFace, SW, SE):
        mComp = comps[mComp]
        fComp = comps[fComp]
        self.faces = numpy.zeros((1,2),int)
        self.faces[0,:] = [1,2]
        
        SW = numpy.array(SW,int)
        SE = numpy.array(SE,int)
        E = SE-SW
        E /= abs(sum(E))
        
        self.mComp = mComp
        self.mSide = mSide
        self.fComp = fComp
        self.fFace = fFace
        self.SE = SE
        self.SW = SW
        self.E = E

        Ps = []
        Ks = []

        fPs = fComp.Ps
        fKs = fComp.Ks[fFace]
        mPs = mComp.Ps
        mKs = mComp.Ks

        for i in range(mComp.Ks[0].shape[0]):
            n = fPs[fKs[self.tup(SW)]].shape[abs(E[0])]
            edge2 = self.fSpliceP(SW+(1+i)*E, -1, None)
            edge1 = self.mSpliceP(0,i,True)
            Ps.append(self.createInterface(n, edge1, edge2))

        n = fPs[fKs[self.tup(SW)]].shape[abs(E[1])]
        edge1 = self.fSpliceP(SW, None, 0)
        edge2 = Ps[0][:,0,:]
        Ps.append(self.createInterface(n, edge1, edge2, True))

        n = fPs[fKs[self.tup(SE)]].shape[abs(E[1])]
        edge2 = self.fSpliceP(SE, None, -1)
        edge1 = Ps[-2][:,-1,:]
        Ps.append(self.createInterface(n, edge1, edge2, True))

        Ps[-2][-1,:,:] = self.fSpliceP(SW, -1, None)
        Ps[-1][-1,:,:] = self.fSpliceP(SE, -1, None)

        K = numpy.zeros((1,2+mComp.Ks[0].shape[0]),int)
        counter = 0
        for j in range(mComp.Ks[0].shape[0]):
            K[0,j+1] = counter
            counter += 1
        for i in range(2):
            K[0,-i] = counter
            counter += 1
        Ks.append(K)

        u1 = min(SW[0],SE[0])
        u2 = max(SW[0],SE[0])
        v1 = min(SW[1],SE[1])
        v2 = max(SW[1],SE[1])
        for j in range(v2,v1-1,-1):
            for i in range(u2,u1-1,-1):
                fPs.pop(fKs[i,j])
                for f in range(len(fComp.Ks)):
                    for v in range(fComp.Ks[f].shape[1]):
                        for u in range(fComp.Ks[f].shape[0]):
                            if fComp.Ks[f][u,v] > fKs[i,j]:
                                fComp.Ks[f][u,v] -= 1
                fKs[i,j] = -1
            
        self.Ps = Ps
        self.Ks = Ks

        self.oml0 = []

    def setDOFs(self):
        oml0 = self.oml0
        Ks = self.Ks
        self.setSurfC1(0, val=True)
        self.setSurfC1(0, i=0)
        edge = oml0.surf_edge[self.Ks[0][0,0],1,0]
        if edge > 0:
            oml0.edge_c1[abs(edge)-1,0] = True
        else:
            oml0.edge_c1[abs(edge)-1,1] = True
        edge = oml0.surf_edge[self.Ks[0][0,-1],1,0]
        if edge > 0:
            oml0.edge_c1[abs(edge)-1,1] = True
        else:
            oml0.edge_c1[abs(edge)-1,0] = True

    def isExteriorDOF(self, f, uType, vType, i, j):
        if j==0 or j==self.Ks[0].shape[1]-1:
            return self.check(uType, vType, u=0)
        else:
            return False

    def propagateQs(self):
        si = self.getni(0,0)
        sj = self.getni(0,1)

        Ns = self.Ns
        Qs = self.Qs
        Ks = self.Ks
        fQs = self.fComp.Qs[self.fFace]    
        Qs[0][:,:,:] = 0
        for j in range(int(Ns[0].shape[1]-sj[0]-sj[-1])):
            Q2 = self.fSpliceQ(self.SW, 1, 1, 1, j)
            Q1 = self.mSpliceQ(0, j, True)
            Qs[0][:,j+sj[0],:] += numpy.outer(numpy.linspace(1,0,si[0]+2)[:-1],Q1)
            Qs[0][:,j+sj[0],:] += numpy.outer(numpy.linspace(0,1,si[0]+2)[:-1],Q2)
        for i in range(si[0]):
            Q1 = self.fSpliceQ(self.SW, 0, 0, i, -1)
            Q2 = Qs[0][i,sj[0]]
            Qs[0][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(1,0,sj[0]+2)[1:],Q1)
            Qs[0][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(0,1,sj[0]+2)[1:],Q2)

            Q2 = self.fSpliceQ(self.SW, 0, Ks[0].shape[1], i, -1)
            Q1 = Qs[0][i,sum(sj[:-1])]
            Qs[0][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(1,0,sj[-1]+2)[:-1],Q1)
            Qs[0][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(0,1,sj[-1]+2)[:-1],Q2)
