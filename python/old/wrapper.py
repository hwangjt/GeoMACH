from __future__ import division
import GeoMACH
import numpy, copy


class curve:

    def initialize(self,k,m,n):
        self.k = k
        self.m = m
        self.n = n
        self.d = GeoMACH.knotopen(k,k+m)
        self.t = GeoMACH.paramuni(k+m,m,n,self.d) 
        self.evaluateBases()      

    def create(self,C,n,k=4):
        self.initialize(k,C.shape[0],n)
        self.C = copy.copy(C)
        return self.evaluatePts()

    def fit(self,P0,m,k=4):
        self.initialize(k,m,P0.shape[0])
        dPdC,P = GeoMACH.curvefit(self.k,self.m,self.n,self.B,self.i0,P0)
        self.P1 = P
        self.C = numpy.zeros((m,3))
        self.C[0,:] = P0[0,:]
        self.C[-1,:] = P0[-1,:]
        for i in range(3):
            self.C[1:-1,i] = numpy.linalg.lstsq(dPdC,P[:,i])[0]
        return self.evaluatePts()

    def fit0(self,P0,m,k=4):
        self.initialize(k,m,P0.shape[0])
        self.evaluateJacob()
        self.C = numpy.zeros((m,3))
        for i in range(3):
            self.C[:,i] = numpy.linalg.lstsq(self.dPdC,P0[:,i])[0]
        return self.evaluatePts()

    def evaluateBases(self):
        self.B = numpy.zeros((self.n,self.k))
        self.i0 = numpy.zeros(self.n)
        for u in range(self.n):
            self.B[u,:], self.i0[u] = GeoMACH.basis(self.k,self.k+self.m,self.t[u],self.d)
        return copy.copy(self.B)

    def evaluatePts(self):
        self.P = numpy.zeros((self.n,3))
        for i in range(3):
            self.P[:,i] = GeoMACH.evalcurve(self.k, self.m, self.n, self.B, self.i0, self.C[:,i])
        return copy.copy(self.P)

    def evaluateJacob(self):
        self.dPdC = GeoMACH.curvejacob(self.k,self.m,self.n,self.B,self.i0)
        return copy.copy(self.dPdC)

    def evaluatedPdu(self,u):
        Bu, i0 = GeoMACH.basis1(self.k,self.k+self.m,self.t[u],self.d)
        dPdu = numpy.zeros(3)
        for i in range(3):
            dPdu[i] = GeoMACH.evalcurve(self.k,self.m,1,Bu,i0,self.C[:,i])
        return dPdu
        

class surface:  

    def initialize(self,k,m,n):
        self.k = k
        self.m = m
        self.n = n
        self.d1 = GeoMACH.knotopen(k[0],k[0]+m[0])
        self.d2 = GeoMACH.knotopen(k[1],k[1]+m[1])
        self.t = numpy.zeros((n[0],n[1],2))
        self.t[:,0,0] = GeoMACH.paramuni(k[0]+m[0],m[0],n[0],self.d1)
        self.t[0,:,1] = GeoMACH.paramuni(k[1]+m[1],m[1],n[1],self.d2)
        for i in range(1,n[1]):
            self.t[:,i,0] = self.t[:,0,0]
        for i in range(1,n[0]):
            self.t[i,:,1] = self.t[0,:,1]
        self.evaluateBases()    

    def create(self,C,n,k=[4,4]):
        self.initialize(k,C.shape[0:2],n)
        self.C = copy.copy(C)
        return self.evaluatePts()

    def fit(self,P0,C0,m,k=[4,4]):
        self.initialize(k,m,P0.shape[0:2])
        dPdC,P = GeoMACH.surfacefit(self.k[0],self.k[1],self.m[0],self.m[1],self.n[0],self.n[1],self.B1,self.B2,self.i0,self.j0,P0,C0)
        self.C = numpy.zeros((m[0],m[1],3))
        for i in range(3):
            X = numpy.linalg.lstsq(dPdC,P[:,i])[0]
            self.C[:,:,i] = GeoMACH.expand2d(m[0],m[1],X,C0[:,:,i])
        return self.evaluatePts()

    def fit0(self,P0,m,k=[4,4]):
        self.initialize(k,m,P0.shape[0:2])
        self.evaluateJacob()
        self.C = numpy.zeros((m[0],m[1],3))
        for i in range(3):
            Q = GeoMACH.flatten2d(self.n[0],self.n[1],P0[:,:,i])
            D = numpy.linalg.lstsq(self.dPdC,Q)[0]
            self.C[:,:,i] = GeoMACH.expand2d0(m[0],m[1],D)
        return self.evaluatePts()

    def evaluateBases(self):
        self.B1 = numpy.zeros((self.n[0],self.n[1],self.k[0]))
        self.B2 = numpy.zeros((self.n[0],self.n[1],self.k[1]))
        self.i0 = numpy.zeros((self.n[0],self.n[1]))
        self.j0 = numpy.zeros((self.n[0],self.n[1]))
        for u in range(self.n[0]):
            for v in range(self.n[1]):
                self.B1[u,v,:], self.i0[u,v] = GeoMACH.basis(self.k[0],self.k[0]+self.m[0],self.t[u,v,0],self.d1)
                self.B2[u,v,:], self.j0[u,v] = GeoMACH.basis(self.k[1],self.k[1]+self.m[1],self.t[u,v,1],self.d2)

    def evaluatePts(self):
        self.P = numpy.zeros((self.n[0],self.n[1],3))
        for i in range(3):
            self.P[:,:,i] = GeoMACH.evalsurface(self.k[0], self.k[1], self.m[0], self.m[1], self.n[0], self.n[1], self.B1, self.B2, self.i0, self.j0, self.C[:,:,i])
        return copy.copy(self.P)

    def evaluateJacob(self):
        self.dPdC = GeoMACH.surfacejacob(self.k[0], self.k[1], self.m[0], self.m[1], self.n[0], self.n[1], self.B1, self.B2, self.i0, self.j0)
        return copy.copy(self.dPdC)

    def evaluatePder(self,u,v,du,dv):
        k = self.k
        m = self.m
        n = self.n
        t = self.t
        d1 = self.d1
        d2 = self.d2
        Bu = numpy.zeros((1,1,k[0]))
        Bv = numpy.zeros((1,1,k[1]))
        i0 = numpy.zeros((1,1))
        j0 = numpy.zeros((1,1))
        if du==0:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis(k[0],k[0]+m[0],t[u,v,0],d1)
        elif du==1:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis1(k[0],k[0]+m[0],t[u,v,0],d1)
        elif du==2:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis2(k[0],k[0]+m[0],t[u,v,0],d1)
        if dv==0:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis(k[1],k[1]+m[1],t[u,v,1],d2)
        elif dv==1:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis1(k[1],k[1]+m[1],t[u,v,1],d2)
        elif dv==2:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis2(k[1],k[1]+m[1],t[u,v,1],d2)
        Pder = numpy.zeros(3)
        for i in range(3):
            Pder[i] = GeoMACH.evalsurface(k[0], k[1],m[0],m[1],1,1,Bu,Bv,i0,j0,self.C[:,:,i])
        return Pder

    def evaluatePderAt(self,t0,t1,du,dv):
        k = self.k
        m = self.m
        n = self.n
        d1 = self.d1
        d2 = self.d2
        t0 = GeoMACH.param(k[0]+m[0],m[0],n[0],d1,t0)
        t1 = GeoMACH.param(k[1]+m[1],m[1],n[1],d2,t1)
        Bu = numpy.zeros((1,1,k[0]))
        Bv = numpy.zeros((1,1,k[1]))
        i0 = numpy.zeros((1,1))
        j0 = numpy.zeros((1,1))
        if du==0:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis(k[0],k[0]+m[0],t0,d1)
        elif du==1:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis1(k[0],k[0]+m[0],t0,d1)
        elif du==2:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis2(k[0],k[0]+m[0],t0,d1)
        if dv==0:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis(k[1],k[1]+m[1],t1,d2)
        elif dv==1:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis1(k[1],k[1]+m[1],t1,d2)
        elif dv==2:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis2(k[1],k[1]+m[1],t1,d2)
        Pder = numpy.zeros(3)
        for i in range(3):
            Pder[i] = GeoMACH.evalsurface(k[0],k[1],m[0],m[1],1,1,Bu,Bv,i0,j0,self.C[:,:,i])
        return Pder

    def evaluatePderAt2(self,t0,t1,du,dv):
        k = self.k
        m = self.m
        n = self.n
        d1 = self.d1
        d2 = self.d2
        Bu = numpy.zeros((1,1,k[0]))
        Bv = numpy.zeros((1,1,k[1]))
        i0 = numpy.zeros((1,1))
        j0 = numpy.zeros((1,1))
        if du==0:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis(k[0],k[0]+m[0],t0,d1)
        elif du==1:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis1(k[0],k[0]+m[0],t0,d1)
        elif du==2:
            Bu[0,0,:], i0[0,0] = GeoMACH.basis2(k[0],k[0]+m[0],t0,d1)
        if dv==0:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis(k[1],k[1]+m[1],t1,d2)
        elif dv==1:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis1(k[1],k[1]+m[1],t1,d2)
        elif dv==2:
            Bv[0,0,:], j0[0,0] = GeoMACH.basis2(k[1],k[1]+m[1],t1,d2)
        Pder = numpy.zeros(3)
        for i in range(3):
            Pder[i] = GeoMACH.evalsurface(k[0],k[1],m[0],m[1],1,1,Bu,Bv,i0,j0,self.C[:,:,i])
        return Pder

    def getCentroid(self):
        cen = numpy.zeros(3)
        for i in range(3):
            cen[i] = numpy.sum(self.P[:,:,i])/self.P.shape[0]/self.P.shape[1]
        return cen

    def getMid(self):
        u = int(self.P.shape[0]/2.0)
        v = int(self.P.shape[1]/2.0)
        return self.P[u,v,:]

    def projectPt(self, P0, t0, t1):
        k = self.k
        m = self.m
        n = self.n
        d1 = self.d1
        d2 = self.d2
        u = GeoMACH.param(k[0]+m[0],m[0],n[0],d1,t0)
        v = GeoMACH.param(k[1]+m[1],m[1],n[1],d2,t1)

        P = self.evaluatePderAt2(u,v,0,0)
        dPdu = self.evaluatePderAt2(u,v,1,0)
        dPdv = self.evaluatePderAt2(u,v,0,1)
        d2Pdu2 = self.evaluatePderAt2(u,v,2,0)
        d2Pdv2 = self.evaluatePderAt2(u,v,0,2)
        d2Pdudv = self.evaluatePderAt2(u,v,1,1)
        dP = P-P0
        Ru = numpy.dot(dP,dPdu)
        Rv = numpy.dot(dP,dPdv)
        R = Ru**2 + Rv**2
        h = 1e-5
        counter = 0
        while R > 1e-10 and counter < 100:
            dRdu = 0
            dRdv = 0
            for i in range(3):
                dRdu += 2*Ru*(dPdu[i]**2 + dP[i]*d2Pdu2[i])
                dRdv += 2*Rv*(dPdv[i]**2 + dP[i]*d2Pdv2[i])
                dRdu += 2*Rv*(dPdu[i]*dPdv[i] + dP[i]*d2Pdudv[i])
                dRdv += 2*Ru*(dPdu[i]*dPdv[i] + dP[i]*d2Pdudv[i])

            dRdu = (R1-R)/h
            dRdv = (R2-R)/h
            
            u -= dRdu
            v -= dRdv
            if u<0 or u>1:
                u = 1 - u%1
            if v<0 or v>1:
                v = 1 - v%1
            P = self.evaluatePderAt2(u,v,0,0)
            dPdu = self.evaluatePderAt2(u,v,1,0)
            dPdv = self.evaluatePderAt2(u,v,0,1)
            d2Pdu2 = self.evaluatePderAt2(u,v,2,0)
            d2Pdv2 = self.evaluatePderAt2(u,v,0,2)
            d2Pdudv = self.evaluatePderAt2(u,v,1,1)
            dP = P-P0
            Ru = numpy.dot(dP,dPdu)
            Rv = numpy.dot(dP,dPdv)
            R = Ru**2 + Rv**2
            counter += 1
            print counter,R,u,v

        return P,u,v

        
        
