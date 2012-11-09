from __future__ import division
from layout import Layout
import numpy, pylab, time, copy
import scipy.sparse
import PAM.PAMlib as PAMlib
import mpl_toolkits.mplot3d.axes3d as p3


class Component(object):
    """ Base class for Wing, Body, and Junction components. """

    def __init__(self):
        self.Ps = []
        self.Ks = []   
        self.oml0 = []
        self.faces = []
        
    def addFace(self, du, dv, d, ru=0.4, rv=0.4):
        """ Creates a set of rectangular surfaces, their IDs, and face dims.
        nu,nv: number of surfaces in the u and v directions
        du,dv: {1,2,3} maps to {x,y,z}; negative sign means reverse order
        d: position of the surfaces in the remaining coordinate axis
        ru,rv: surfaces span -ru to +ru in u dir. and -rv to +rv in v dir.

        Adds to self.Ps and self.Ks
        """
        self.faces.append([du,dv])

        n = 10
        nu = self.ms[abs(du)-1].shape[0]
        nv = self.ms[abs(dv)-1].shape[0]
        Ps = self.Ps
        Ks = self.Ks
        for j in range(nv):
            for i in range(nu):
                u1 = ru*(-1+2*i/nu)
                u2 = ru*(-1+2*(i+1)/nu)
                v1 = rv*(-1+2*j/nv)
                v2 = rv*(-1+2*(j+1)/nv)
                P = PAMlib.createsurfaces(n,du,dv,d,u1,u2,v1,v2)
                Ps.append(P)  

        K = numpy.zeros((nu,nv),int)
        counter = 0
        if len(Ks) > 0:
            counter = numpy.max(Ks[-1]) + 1
        for j in range(nv):
            for i in range(nu):
                K[i,j] = counter
                counter += 1
        Ks.append(K)

    def connectEdges(self, f1=0, u1=None, v1=None, f2=0, u2=None, v2=None):
        def edge(f, u, v, kk, P0=None):
            Ks = self.Ks
            Ps = self.Ps
            d = 0 if u==None else 1
            r = self.faces[f][d]
            k = kk if r > 0 else -1-kk
            surf = Ks[f][k,v] if d==0 else Ks[f][u,k]
            if not (P0==None):
                if d == 0:
                    if r > 0:
                        Ps[surf][:,v] = P0
                    else:
                        Ps[surf][::-1,v] = P0
                else:
                    if r > 0:
                        Ps[surf][u,:] = P0
                    else:
                        Ps[surf][u,::-1] = P0
            P = Ps[surf][:,v] if d==0 else Ps[surf][u,:]
            return P[::-1] if r==-1 else P                
            
        for k in range(self.Ks[f1].shape[v1==None]):
            avg = 0.5*edge(f1,u1,v1,k) + 0.5*edge(f2,u2,v2,k)
            edge(f1,u1,v1,k,P0=avg)
            edge(f2,u2,v2,k,P0=avg)  

    def setC1(self, t, f, i=None, j=None, u=None, v=None, d=None, val=True):
        """ Set C1 continuity 
        t: {string} surface or edge C1
        f: face index
        i,j: surface index
            both given: only consider [i,j] surface
            one given: loop through and apply to all of the other index
            none given: apply to all surfaces
        u,v: edge/vert index (for surfC1)
            both given: only consider [u,v] corner/side
            one given: loop through and apply to all of the other index
            none given: apply to all corners/sides
        u,v,d: side index (for edgeC1)
        """

        def setSurfC1(f, i, j, u, v, d, val):
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


        def setEdgeC1(f, i, j, u, v, d, val):
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

        if t=='surf':
            func = setSurfC1
        elif t=='edge':
            func = setEdgeC1
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

    def setCornerC1(self, f, i=0, j=0, val=True):
        self.setC1('edge', f, i=i, j=j, u=i, d=j, val=val)
        self.setC1('edge', f, i=i, j=j, v=j, d=i, val=val)

    def computems(self):
        oml0 = self.oml0
        Ks = self.Ks
        for f in range(len(Ks)):
            for i in range(Ks[f].shape[0]):
                for j in range(Ks[f].shape[1]):
                    surf = Ks[f][i,j]
                    if not surf==-1:
                        for k in range(2):
                            edge = oml0.surf_edge[surf,k,0]
                            group = oml0.edge_group[abs(edge)-1] - 1
                            m = oml0.group_m[group] - 1
                            if k==0:
                                self.getms(f,k)[i] = int(m)
                            else:
                                self.getms(f,k)[j] = int(m)

    def getms(self, f, d):
        dim = self.faces[f][d]
        if dim > 0:
            return self.ms[abs(dim)-1]
        else:
            return self.ms[abs(dim)-1][::-1]

    def initializeDOFmappings(self):
        def classify(i, n):
            if i==0:
                return 0
            elif i==n-1:
                return 2
            else:
                return 1

        def getC1(surf, u=None, v=None, d=0):
            if u==None:
                edge = self.oml0.surf_edge[surf,0,v]
            else:
                edge = self.oml0.surf_edge[surf,1,u]
            if edge > 0:
                return self.oml0.edge_c1[abs(edge)-1,d]
            else:
                return self.oml0.edge_c1[abs(edge)-1,1-abs(d)]

        oml0 = self.oml0
        Ks = self.Ks
        surf_c1 = oml0.surf_c1
        edge_c1 = oml0.edge_c1

        Qs = []
        Ns = []

        for f in range(len(Ks)):
            ni = self.getms(f,0)
            nj = self.getms(f,1)
            Qs.append(numpy.zeros((sum(ni)+1,sum(nj)+1,3)))
            Ns.append(numpy.zeros((sum(ni)+1,sum(nj)+1,5),int))
            Ns[f][:,:,:] = -1
            for j in range(Ks[f].shape[1]):
                for i in range(Ks[f].shape[0]):
                    surf = Ks[f][i,j]
                    if surf != -1:
                        mu,mv = oml0.getEdgeProperty(surf,1)
                        for v in range(mv):
                            jj = sum(nj[:j]) + v
                            vType = classify(v,mv)
                            for u in range(mu):
                                ii = sum(ni[:i]) + u
                                uType = classify(u,mu)
                                DOF = True
                                if uType==0 or uType==2 or vType==0 or vType==2:
                                    DOF = DOF and not surf_c1[surf,uType,vType]
                                    if (not uType==1) and (not vType==1):
                                        DOF = DOF and not getC1(surf,u=int(uType/2),d=int(vType/2))
                                        DOF = DOF and not getC1(surf,v=int(vType/2),d=int(uType/2))
                                if DOF:
                                    Ns[f][ii,jj,0] = oml0.computeIndex(surf,u,v,2)
                                    Ns[f][ii,jj,1] = i
                                    Ns[f][ii,jj,2] = u
                                    Ns[f][ii,jj,3] = j
                                    Ns[f][ii,jj,4] = v
        self.Qs = Qs
        self.Ns = Ns

    def updateQs(self):
        Qs = self.Qs
        Ns = self.Ns
        oml0 = self.oml0
        for f in range(len(Ns)):
            PAMlib.updateqs(oml0.nQ, Ns[f].shape[0], Ns[f].shape[1], Ns[f], Qs[f], oml0.Q)

    def translatePoints(self, dx, dy, dz):
        for k in range(len(self.Ps)):
            self.Ps[k][:,:,0] += dx
            self.Ps[k][:,:,1] += dy
            self.Ps[k][:,:,2] += dz


class Component2(object):
    """ Base class for Wing, Body, and Junction components. """

    def __init__(self):
        self.Ps = []
        self.Ks = []   
        self.oml0 = []      
        
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
        """ Set C1 continuity 
        t: {string} surface or edge C1
        f: face index
        i,j: surface index
            both given: only consider [i,j] surface
            one given: loop through and apply to all of the other index
            none given: apply to all surfaces
        u,v: edge/vert index (for surfC1)
            both given: only consider [u,v] corner/side
            one given: loop through and apply to all of the other index
            none given: apply to all corners/sides
        u,v,d: side index (for edgeC1)
        """

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

    def flattenGeometry(self):
        oml0 = self.oml0
        Ms = self.Ms
        for f in range(len(Ms)):      
            ni = Ms[f].shape[0]
            nj = Ms[f].shape[1]
            for i in range(ni):
                for j in range(nj):
                    oml0.C[Ms[f][i,j],:] = self.getFlattenedC(f, i, j, ni, nj)
        oml0.computePointsC()

    def initializeLayout(self):
        def getv(r, v1, v2):
            return numpy.array(v1)*(1-r) + numpy.array(v2)*r

        oml0 = self.oml0
        Ms = self.Ms
        Js = self.Js

        edges = []
        edges.append([0,0,0,1])
        edges.append([0,0,1,0])
        edges.append([0,1,1,1])
        edges.append([1,0,1,1])

        skinIndices = self.getSkinIndices()
        for s in range(len(skinIndices)):
            for f in skinIndices[s]:
                for k in range(Js[f].shape[0]):
                    C11 = oml0.C[Ms[f][Js[f][k,0],Js[f][k,1]],:2]
                    C12 = oml0.C[Ms[f][Js[f][k,0],Js[f][k,3]],:2]
                    C21 = oml0.C[Ms[f][Js[f][k,2],Js[f][k,1]],:2]
                    C22 = oml0.C[Ms[f][Js[f][k,2],Js[f][k,3]],:2]
                    edges.append([C11[0],C11[1],C12[0],C12[1]])
                    edges.append([C11[0],C11[1],C21[0],C21[1]])
                    edges.append([C12[0],C12[1],C22[0],C22[1]])
                    edges.append([C21[0],C21[1],C22[0],C22[1]])
        for m in self.keys:
            member = self.members[m]
            for i in range(member.nmem):
                if member.nmem==1:
                    ii = 0
                else:
                    ii = i/(member.nmem-1)
                A = getv(ii, member.A1, member.A2)
                B = getv(ii, member.B1, member.B2)
                C = getv(ii, member.C1, member.C2)
                D = getv(ii, member.D1, member.D2)
                for j in range(member.ndiv):
                    if (B[2]==0 and C[2]==0) or (B[2]==1 and C[2]==1):
                        V0 = getv(j/member.ndiv, B, C)
                        V1 = getv((j+1)/member.ndiv, B, C)
                        edges.append([V0[0],V0[1],V1[0],V1[1]])
                    if (A[2]==0 and D[2]==0) or (A[2]==1 and D[2]==1):
                        V0 = getv(j/member.ndiv, A, D)
                        V1 = getv((j+1)/member.ndiv, A, D)
                        edges.append([V0[0],V0[1],V1[0],V1[1]])
        self.layout = Layout(6, self.getAR(),1,numpy.array(edges))

    def findJunctionQuadsAndEdges(self):
        oml0 = self.oml0
        Ms = self.Ms
        Js = self.Js

        C = numpy.zeros((4,2))
        ctd = numpy.zeros(2)

        nquad = self.layout.nquad
        nedge = self.layout.nedge
        edges = self.layout.edges
        verts = self.layout.verts
        poly_vert = self.layout.poly_vert

        skinIndices = self.getSkinIndices()

        JQs = []
        JEs = []
        quad_indices = []
        for s in range(len(skinIndices)):
            JQ = []
            for f in skinIndices[s]:
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
            if JQ==[]:
                JQ.append(-1)
            JQs.append(JQ)

            JE = []
            for f in skinIndices[s]:
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
                            JEs.append(e+1)
            JEs.append(JE)

            quad_indices.append(self.layout.getQuadIndices(JQ))
        
        if len(quad_indices)==1:
            quad_indices.append(quad_indices[0])

        self.JQs = JQs
        self.JEs = JEs
        self.quad_indices = numpy.array(quad_indices, order='F').T

    def projectSkins(self):
        oml0 = self.oml0
        Ms = self.Ms
        Js = self.Js
        skinIndices = self.getSkinIndices()

        C = numpy.zeros((4,2),order='F')
        Bs = []
        Ss = []
        quad_indices = []
        for s in range(len(skinIndices)):
            Ps = []
            P, S = self.layout.extractFlattened(self.JQs[s], max(self.quad_indices[:,s]))
            Ps.append(P)
            Ss.append(S)
            for f in skinIndices[s]:
                for k in range(Js[f].shape[0]):
                    C[0,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,1]],:2]
                    C[1,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,3]],:2]
                    C[2,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,1]],:2]
                    C[3,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,3]],:2]
                    P = self.layout.extractEdges(self.JQs[s], min(C[:,0]), max(C[:,0]), min(C[:,1]), max(C[:,1]))
                    Ps.append(P)
            P = numpy.vstack(Ps)

            surfindices = []
            for f in skinIndices[s]:
                surfindices.append(self.Ks[f].flatten())
            surfs = numpy.unique(numpy.hstack(surfindices))
            if surfs[0] == -1:            
                surfs = numpy.delete(surfs,0)
            Q = numpy.zeros((P.shape[0],3))
            Q[:,2] = 1.0
            ss, u, v = oml0.computeProjection(P, surfs=surfs, Q=Q)
            Bs.append(oml0.computeBases(ss,u,v))

        B = oml0.vstackSparse(Bs)

        BM = B.dot(oml0.M)
        P = BM.dot(oml0.Q)

        As,S = self.computeStructure(P)
        Ss.append(S)
        A = oml0.vstackSparse(As)

        ABM = A.dot(BM)
        S = numpy.vstack(Ss)

        self.strABM = ABM
        self.strS = S

    def computeStructure(self, P):
        oml0 = self.oml0
        Ms = self.Ms
        Js = self.Js
        layout = self.layout
        skinIndices = self.getSkinIndices()

        C = numpy.zeros((4,2),order='F')
        As = []
        Jns = []
        Jus = []
        col = 0
        for s in range(len(skinIndices)):
            n = max(self.quad_indices[:,s])*self.layout.n**2
            Aa = numpy.ones(n)
            Ai = numpy.linspace(0, n-1, n)
            Aj = numpy.linspace(0, n-1, n) + col
            As.append(scipy.sparse.csr_matrix((Aa,(Ai,Aj)), shape=(n,P.shape[0])))
            Jn = []
            Ju = []
            col += n
            for f in skinIndices[s]:
                for k in range(Js[f].shape[0]):
                    C[0,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,1]],:2]
                    C[1,:] = oml0.C[Ms[f][Js[f][k,0],Js[f][k,3]],:2]
                    C[2,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,1]],:2]
                    C[3,:] = oml0.C[Ms[f][Js[f][k,2],Js[f][k,3]],:2]
                    nu1, nu2, nv1, nv2 = layout.countJunctionEdges(self.JQs[s], min(C[:,0]), max(C[:,0]), min(C[:,1]), max(C[:,1]))
                    Jn.append([nu1, nu2, nv1, nv2])
                    Ju.append([min(C[:,0]), max(C[:,0]), min(C[:,1]), max(C[:,1])])
                    nP = layout.n*(nu1 + nu2 + nv1 + nv2)
                    col += nP
            if Jn==[]:
                Jn.append([-1,-1,-1,-1])
                Ju.append([-1,-1,-1,-1])
            Jns.append(numpy.array(Jn, int, order='F'))
            Jus.append(numpy.array(Ju, order='F'))
            
        nM = len(self.keys)
        members = numpy.zeros((nM,34),order='F')
        for i in range(nM):
            member = self.members[self.keys[i]]
            members[i,:4] = [member.domain, member.shape, member.nmem, member.ndiv]
            members[i, 4:16] = member.A1 + member.B1 + member.C1 + member.D1
            members[i,16:28] = member.A2 + member.B2 + member.C2 + member.D2
            members[i,28:] = [member.tx, member.ty, member.rx, member.ry, member.n1, member.n2]

        nP, nS = PAMlib.countinternalnodes(self.layout.n, nM, members)
        P2, M, S = PAMlib.computeinternalnodes(nP, nS, self.layout.n, nM, members)
        nA = PAMlib.countannz(nP, layout.nvert, layout.nquad, layout.verts, layout.poly_vert, self.quad_indices, P2, M)
        if len(skinIndices)==1:
            Aa, Ai, Aj = PAMlib.assembleamtx(nA, self.layout.n, nP, Jns[0].shape[0], Jns[0].shape[0], self.layout.nvert, self.layout.nquad, Jns[0], Jns[0], Jus[0], Jus[0], self.quad_indices, self.layout.verts, self.layout.poly_vert, P2, M)
        else:
            Aa, Ai, Aj = PAMlib.assembleamtx(nA, self.layout.n, nP, Jns[0].shape[0], Jns[1].shape[0], self.layout.nvert, self.layout.nquad, Jns[0], Jns[1], Jus[0], Jus[1], self.quad_indices, self.layout.verts, self.layout.poly_vert, P2, M)
        As.append(scipy.sparse.csr_matrix((Aa,(Ai,Aj)), shape=(max(Ai)+1,P.shape[0])))
            
        return As, S
        
        
    def buildStructure(self):
        self.computeMs()
        self.findJunctions()
        self.flattenGeometry()
        self.initializeLayout()
        self.findJunctionQuadsAndEdges()
        self.projectSkins()

    def addMembers(self, name, domain, shape, nmem, ndiv, 
                   A1=None, B1=None, C1=None, D1=None, 
                   A2=None, B2=None, C2=None, D2=None, 
                   tx=0.5, ty=0.5, rx=1.0, ry=1.0, n1=5, n2=5):

        if A2==None:
            A2 = [-1,-1,-1]
            B2 = [-1,-1,-1]
            C2 = [-1,-1,-1]
            D2 = [-1,-1,-1]
        if B1==None:
            B1 = [-1,-1,-1]
            B1[:2] = A1[:2]
            B1[2] = C1[2]
        if D1==None:
            D1 = [-1,-1,-1]
            D1[:2] = C1[:2]
            D1[2] = A1[2]
        if B2==None:
            B2 = [-1,-1,-1]
            B2[:2] = A2[:2]
            B2[2] = C2[2]
        if D2==None:
            D2 = [-1,-1,-1]
            D2[:2] = C2[:2]
            D2[2] = A2[2]
        self.members[name] = Member(domain, shape, nmem, ndiv, A1, B1, C1, D1, 
                                    A2, B2, C2, D2, tx, ty, rx, ry, n1, n2)
        self.keys.append(name)


class Member(object):

    def __init__(self, domain, shape, nmem, ndiv, 
                 A1, B1, C1, D1, A2, B2, C2, D2, tx, ty, rx, ry, n1, n2):
        self.domain = domain
        self.shape = shape
        self.nmem = nmem
        self.ndiv = ndiv
        self.A1 = A1
        self.B1 = B1
        self.C1 = C1
        self.D1 = D1
        self.A2 = A2
        self.B2 = B2
        self.C2 = C2
        self.D2 = D2
        self.tx = tx
        self.ty = ty
        self.rx = rx
        self.ry = ry
        self.n1 = n1
        self.n2 = n2


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
