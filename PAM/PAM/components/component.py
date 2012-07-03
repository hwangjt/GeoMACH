from __future__ import division
import numpy, pylab, time
import PAM.PAMlib as PAMlib
import mpl_toolkits.mplot3d.axes3d as p3


class Component(object):
        
    def createSurfaces(self, Ks, nu, nv, du, dv, d):
        Ps = []
        for j in range(len(nv)):
            for i in range(len(nu)):
                u1 = (sum(nu[:i])-i)/(sum(nu)-len(nu))
                u2 = (sum(nu[:i+1])-i-1)/(sum(nu)-len(nu))
                v1 = (sum(nv[:j])-j)/(sum(nv)-len(nv))
                v2 = (sum(nv[:j+1])-j-1)/(sum(nv)-len(nv))
                P = PAMlib.createsurfaces(nu[i],nv[j],du,dv,d,u1,u2,v1,v2)
                Ps.append(P)  

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
        P = PAMlib.createinterface(n, edge1.shape[0], edge1, edge2)
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
            PAMlib.updateqs(oml0.nQ, Ns[f].shape[0], Ns[f].shape[1], Ns[f], Qs[f], oml0.Q)

    def initializeDOFs(self):
        oml0 = self.oml0
        Ks = self.Ks

        Qs = []
        Ns = []

        for f in range(len(Ks)):
            ni = self.getni(f,0)
            nj = self.getni(f,1)
            Qs.append(numpy.zeros((sum(ni)+1,sum(nj)+1,3),complex))
            Ns.append(numpy.zeros((sum(ni)+1,sum(nj)+1,5),int))
            Ns[f][:,:,:] = -1
            for j in range(Ks[f].shape[1]):
                for i in range(Ks[f].shape[0]):
                    surf = Ks[f][i,j]
                    if surf != -1:
                        for v in range(nj[j]+1):
                            jj = sum(nj[:j]) + v
                            for u in range(ni[i]+1):
                                ii = sum(ni[:i]) + u
                                uType = self.classifyC(u,ii,ni[i]+1,Ns[f].shape[0])
                                vType = self.classifyC(v,jj,nj[j]+1,Ns[f].shape[1])
                                isInteriorDOF = (uType==2 and vType==2)
                                if isInteriorDOF or self.isExteriorDOF(f,uType,vType,i,j):
                                    Ns[f][ii,jj,0] = oml0.computeIndex(surf,u,v,2)
                                    Ns[f][ii,jj,1] = i
                                    Ns[f][ii,jj,2] = u
                                    Ns[f][ii,jj,3] = j
                                    Ns[f][ii,jj,4] = v
        self.Qs = Qs
        self.Ns = Ns

    def classifyC(self, u, i, lenu, leni):
        if i==0:
            return 0
        elif i==leni-1:
            return -1
        elif u==0:
            return 1
        elif u==lenu-1:
            return 1
        else:
            return 2

    def computeDims(self, aircraft):
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

    def setC1(self, t, f, i=None, j=None, u=None, v=None, d=None, val=True):
        if t=='surf':
            func = self.setSurfC1
        elif t=='edge':
            func = self.setEdgeC1
        if (not i==None) and (not j==None):
            func(f, i, j, u, v, d, val)
        elif not i==None:
            for j in range(self.Ks[f].shape[1]):
                func(f, i, j, u, v, d, val)
        elif not j==None:
            for i in range(self.Ks[f].shape[0]):
                func(f, i, j, u, v, d, val)
        else:
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    func(f, i, j, u, v, d, val)

    def setSurfC1(self, f, i, j, u, v, d, val):
        oml0 = self.oml0
        surf = self.Ks[f][i,j]
        if not surf==-1:
            if u==None and v==None:
                oml0.surf_c1[surf,:,:] = val                    
            elif u==None:
                oml0.surf_c1[surf,:,v] = val
            elif v==None:
                oml0.surf_c1[surf,u,:] = val
            else:
                oml0.surf_c1[surf,u,v] = val


    def setEdgeC1(self, f, i, j, u, v, d, val):
        oml0 = self.oml0
        surf = self.Ks[f][i,j]
        if not surf==-1:
            if u==None:
                edge = oml0.surf_edge[surf,0,v]
            else:
                edge = oml0.surf_edge[surf,1,u]
            if d==None:
                oml0.edge_c1[abs(edge)-1,:] = val
            elif edge>0:
                oml0.edge_c1[abs(edge)-1,d] = val
            else:
                oml0.edge_c1[abs(edge)-1,1-abs(d)] = val

    def setCornerC1(self, f, i=0, j=0, val=True):
        self.setC1('edge', f, i=i, j=j, u=i, d=j, val=val)
        self.setC1('edge', f, i=i, j=j, v=j, d=i, val=val)

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

    def computeMs(self):
        oml0 = self.oml0
        Ks = self.Ks

        Ms = []
        for f in range(len(Ks)):
            ni = self.getni(f,0)
            nj = self.getni(f,1)
            Ms.append(numpy.zeros((sum(ni)+1,sum(nj)+1),int))
            Ms[f][:,:] = -1
            for j in range(Ks[f].shape[1]):
                for i in range(Ks[f].shape[0]):
                    surf = Ks[f][i,j]
                    if surf != -1:
                        for v in range(nj[j]+1):
                            jj = sum(nj[:j]) + v
                            for u in range(ni[i]+1):
                                ii = sum(ni[:i]) + u
                                Ms[f][ii,jj] = oml0.computeIndex(surf,u,v,1)
        self.Ms = Ms

    def findJunctions(self):
        oml0 = self.oml0
        Ks = self.Ks
        Ms = self.Ms

        Js = []
        for f in range(len(Ks)):
            ni = self.getni(f,0)
            nj = self.getni(f,1)
            k0 = Ks[f].shape[0]
            k1 = Ks[f].shape[1]
            SPs = []
            EPs = []
            for j in range(k1):
                for i in range(k0):
                    if Ks[f][i,j]==-1:
                        SPs.append([i,j])
                        EPs.append([i,j])
                        if (i > 0 and Ks[f][i-1,j] == -1) or (j > 0 and Ks[f][i,j-1] == -1):
                            SPs.pop()
                        if (i < k0-1 and Ks[f][i+1,j] == -1) or (j < k1-1 and Ks[f][i,j+1] == -1):
                            EPs.pop()
            J = numpy.zeros((len(SPs),4),int)
            for s in range(len(SPs)):
                SP = SPs[s]
                for e in range(len(EPs)):
                    EP = EPs[e]
                    if SP[0] <= EP[0] and SP[1] <= EP[1]:
                        if numpy.linalg.norm(1+Ks[f][SP[0]:EP[0]+1,SP[1]:EP[1]+1]) < 1e-14:
                            J[s,:] = [sum(ni[:SP[0]]),sum(nj[:SP[1]]),sum(ni[:EP[0]+1]),sum(nj[:EP[1]+1])]
            Js.append(J)
        self.Js = Js

    def addJunctionEdges(self):
        oml0 = self.oml0
        Ms = self.Ms
        Js = self.Js
        ctr = 0
        for f in range(len(Js)):
            for k in range(Js[f].shape[0]):
                C11 = oml0.C[Ms[f][Js[f][k,0],Js[f][k,1]],:2]
                C12 = oml0.C[Ms[f][Js[f][k,0],Js[f][k,3]],:2]
                C21 = oml0.C[Ms[f][Js[f][k,2],Js[f][k,1]],:2]
                C22 = oml0.C[Ms[f][Js[f][k,2],Js[f][k,3]],:2]      
                self.structure.addMembers('Junction'+str(ctr+0), 1, 1, SP1=C11, EP1=C12)
                self.structure.addMembers('Junction'+str(ctr+1), 1, 1, SP1=C21, EP1=C22)
                self.structure.addMembers('Junction'+str(ctr+2), 1, 1, SP1=C11, EP1=C21)
                self.structure.addMembers('Junction'+str(ctr+3), 1, 1, SP1=C12, EP1=C22) 
                ctr += 4

    def findJunctionQuadsAndEdges(self):
        nquad = self.structure.nquad
        nedge = self.structure.nedge
        edges = self.structure.edges
        verts = self.structure.verts
        poly_vert = self.structure.poly_vert

        oml0 = self.oml0
        Ms = self.Ms
        Js = self.Js

        C = numpy.zeros((4,2))
        ctd = numpy.zeros(2)

        JQs = []
        JEs = []
        for f in range(len(Js)):
            JQ = []
            for q in range(nquad):
                ctd[:] = 0
                for k in range(4):
                    ctd[:] += 0.25*verts[poly_vert[q,k]-1,:]
                for k in range(Js[f].shape[0]):
                    C[0,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,1]],:2]
                    C[1,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,3]],:2]
                    C[2,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,1]],:2]
                    C[3,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,3]],:2]
                    if min(C[:,0]) < ctd[0] and ctd[0] < max(C[:,0]) and min(C[:,1]) < ctd[1] and ctd[1] < max(C[:,1]):
                        JQ.append(q+1)
            JQs.append(JQ)
            JE = []
            for e in range(nedge):
                ctd[:] = 0
                for k in range(2):
                    ctd[:] += 0.5*verts[edges[e,k]-1,:]
                for k in range(Js[f].shape[0]):
                    C[0,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,1]],:2]
                    C[1,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,3]],:2]
                    C[2,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,1]],:2]
                    C[3,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,3]],:2]
                    if min(C[:,0]) < ctd[0] and ctd[0] < max(C[:,0]) and min(C[:,1]) < ctd[1] and ctd[1] < max(C[:,1]):
                        JE.append(e+1)
            JEs.append(JE)
        self.JQs = JQs
        self.JEs = JEs


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
