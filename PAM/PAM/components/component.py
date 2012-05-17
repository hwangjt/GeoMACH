from __future__ import division
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class component(object):
        
    def createSurfaces(self, Ks, nu, nv, du, dv, d):
        r = [0,0,0]
        if du<0:
            r[abs(du)-1] = 1
        if dv<0:
            r[abs(dv)-1] = 1
        du = abs(du)-1
        dv = abs(dv)-1

        Ps = []
        for j in range(len(nv)):
            for i in range(len(nu)):
                P = d*numpy.ones((nu[i],nv[j],3),order='F')
                u1 = (sum(nu[:i])-i)/(sum(nu)-len(nu))
                u2 = (sum(nu[:i+1])-i-1)/(sum(nu)-len(nu))
                v1 = (sum(nv[:j])-j)/(sum(nv)-len(nv))
                v2 = (sum(nv[:j+1])-j-1)/(sum(nv)-len(nv))
                for v in range(nv[j]):
                    P[:,v,du] = numpy.linspace(u1,u2,nu[i])
                for u in range(nu[i]):
                    P[u,:,dv] = numpy.linspace(v1,v2,nv[j])
                Ps.append(P)
        for k in range(len(Ps)):
            for i in range(3):
                if r[i]:
                    Ps[k][:,:,i] *= -1
                    Ps[k][:,:,i] += 1    

        K = numpy.zeros((len(nu),len(nv)),int)
        counter = 0
        if len(Ks) > 0:
            counter = numpy.max(Ks[-1]) + 1
        for j in range(len(nv)):
            for i in range(len(nu)):
                K[i,j] = counter
                counter += 1            
        return Ps, K

    def createInterface(self, n, edge1, edge2, swap=False):
        nu = n
        nv = edge1.shape[0]
        P = numpy.zeros((nu,nv,3),order='F')
        for j in range(nv):
            P[:,j,:] += numpy.outer(numpy.linspace(0,1,nu),edge2[j,:])
            P[:,j,:] += numpy.outer(numpy.linspace(1,0,nu),edge1[j,:])
        if swap:
            P = numpy.swapaxes(P,0,1)            
        return P

    def translatePoints(self, dx, dy, dz):
        for k in range(len(self.Ps)):
            self.Ps[k][:,:,0] += dx
            self.Ps[k][:,:,1] += dy
            self.Ps[k][:,:,2] += dz

    def updateQs(self):
        Qs = self.Qs
        Ns = self.Ns
        oml0 = self.oml0

        for f in range(len(Ns)):
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    if Ns[f][i,j,0] != -1:
                        oml0.Q[Ns[f][i,j,0],:] = Qs[f][i,j,:].real

    def initializeDOFs(self):
        oml0 = self.oml0

        Qs = []
        Ns = []

        for f in range(len(self.Ks)):
            ni = self.getni(f,0)
            nj = self.getni(f,1)
            Qs.append(numpy.zeros((sum(ni)+1,sum(nj)+1,3),complex))
            Ns.append(numpy.zeros((sum(ni)+1,sum(nj)+1,5),int))
            Ns[f][:,:,:] = -1
            for j in range(Ns[f].shape[1]):
                v,jj = self.divide(j,nj)
                for i in range(Ns[f].shape[0]):
                    u,ii = self.divide(i,ni)
                    surf = self.Ks[f][ii,jj]
                    uType = self.classifyC(u,i,Ns[f].shape[0])
                    vType = self.classifyC(v,j,Ns[f].shape[1])
                    isInteriorDOF = (uType==2 and vType==2)
                    if surf != -1 and (isInteriorDOF or self.isExteriorDOF(f,uType,vType)):
                        Ns[f][i,j,0] = oml0.computeIndex(surf,u,v,2)
                        Ns[f][i,j,1] = ii
                        Ns[f][i,j,2] = u
                        Ns[f][i,j,3] = jj
                        Ns[f][i,j,4] = v

        self.Qs = Qs
        self.Ns = Ns

    def classifyC(self, u, i, leni):
        if i==0:
            return 0
        elif i==leni-1:
            return -1
        elif u==-1:
            return 1
        else:
            return 2

    def getni(self, f, a):  
        oml0 = self.oml0
        Ks = self.Ks
        ni = numpy.zeros(Ks[f].shape[a],int)
        for i in range(Ks[f].shape[a]):
            surf = Ks[f][i*(1-a),i*a]
            edge = oml0.surf_edge[surf,a,0]
            group = oml0.edge_group[abs(edge)-1] - 1
            m = oml0.group_m[group] - 1
            ni[i] = int(m)
        return ni

    def divide(self, i, ni):
        u = i
        for ii in range(ni.shape[0]):
            if u > ni[ii]:
                u -= ni[ii]
            elif u == ni[ii]:
                u = -1
                break
            else:
                break
        return u, ii
