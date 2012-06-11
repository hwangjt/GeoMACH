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
                    if surf != -1 and (isInteriorDOF or self.isExteriorDOF(f,uType,vType,ii,jj)):
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

    def computeDims(self):
        ndims = int(numpy.max(abs(self.faces)))
        oml0 = self.oml0
        Ks = self.Ks
        dims = []
        for d in range(ndims):
            dims.append([])
        for f in range(self.faces.shape[0]):
            for k in range(2):
                d = abs(self.faces[f,k]) - 1
                dims[d] = numpy.zeros(Ks[f].shape[k],int)
        for f in range(self.faces.shape[0]):
            for i in range(Ks[f].shape[0]):
                for j in range(Ks[f].shape[1]):
                    surf = Ks[f][i,j]
                    if not surf==-1:
                        for k in range(2):
                            edge = oml0.surf_edge[surf,k,0]
                            group = oml0.edge_group[abs(edge)-1] - 1
                            m = oml0.group_m[group] - 1
                            d = abs(self.faces[f,k]) - 1
                            if k==0:
                                index = i
                            else:
                                index = j
                            if self.faces[f,k] > 0:
                                dims[d][index] = int(m)
                            else:
                                dims[d][-index-1] = int(m)
        self.dims = dims
    
    def getni(self, f, a):
        d = abs(self.faces[f,a]) - 1
        if self.faces[f,a] > 0:
            return self.dims[d]
        else:
            return self.dims[d][::-1]

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

    def setSurfC1(self, f, i=None, j=None, val=False):
        oml0 = self.oml0
        if not i==None:
            for j in range(self.Ks[f].shape[1]):
                oml0.surf_c1[self.Ks[f][i,j],i,:] = val
        elif not j==None:
            for i in range(self.Ks[f].shape[0]):
                oml0.surf_c1[self.Ks[f][i,j],:,j] = val
        else:
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,:] = val 

    def setEdgeC1(self, f, i=None, j=None, val=True):
        oml0 = self.oml0
        if not i==None:
            for j in range(self.Ks[f].shape[1]):
                edge = oml0.surf_edge[self.Ks[f][i,j],1,i]
                edge = abs(edge) - 1
                oml0.edge_c1[edge,:] = val
        elif not j==None:
            for i in range(self.Ks[f].shape[0]):
                edge = oml0.surf_edge[self.Ks[f][i,j],0,j]
                edge = abs(edge) - 1
                oml0.edge_c1[edge,:] = val

    def setCornerC1(self, f, i, j, val=False):
        oml0 = self.oml0
        edge = oml0.surf_edge[self.Ks[f][i,j],0,j]
        if edge > 0:
            oml0.edge_c1[abs(edge)-1,i] = val
        else:
            oml0.edge_c1[abs(edge)-1,i+1] = val 
        edge = oml0.surf_edge[self.Ks[f][i,j],1,i]
        if edge > 0:
            oml0.edge_c1[abs(edge)-1,j] = val
        else:
            oml0.edge_c1[abs(edge)-1,j+1] = val 

    def check(self, uType, vType, u=None, v=None):
        if u==None:
            uVal = uType==2
        else:
            uVal = uType==u
        if v==None:
            vVal = vType==2
        else:
            vVal = vType==v
        return uVal and vVal



class Property(object):

    def __init__(self, n0):
        self.data = numpy.zeros(n0)
        self.set([0.0,1.0],[0,1])

    def set(self, val, ind, p=None, w=None, d=None):
        self.G = numpy.zeros((numpy.array(val).shape[0],5))
        self.G[:,0] = val
        self.G[:,1] = ind
        if p==None:
            self.G[:,2] = 2
        else:
            self.G[:,2] = p
        if w==None:
            self.G[:,3] = 0
        else:
            self.G[:,3] = w
        if d==None:
            self.G[:,4] = 0
        else:
            self.G[:,4] = d
        self.evaluate()

    def evaluate(self):
        indices = numpy.round((self.data.shape[0]-1)*self.G[:,1])
        for i in range(self.G.shape[0]-1):
            v1 = self.G[i,0]
            v2 = self.G[i+1,0]
            i1 = int(indices[i])
            i2 = int(indices[i+1])
            p = self.G[i,2]
            w = self.G[i,3]
            x = numpy.linspace(0,1,i2-i1+1)
            if self.G[i,4]==0:
                self.data[i1:i2+1] = v1 + (v2-v1)*(1-w)*x + (v2-v1)*w*x**p
            else:
                self.data[i1:i2+1] = v2 + (v1-v2)*(1-w)*x[::-1] + (v1-v2)*w*x[::-1]**p
