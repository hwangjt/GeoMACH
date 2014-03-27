from __future__ import division
import numpy
import scipy.sparse

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Parameter


class Component(object):
    """ Base class for Wing, Body, and Junction components. """

    def __init__(self):
        self.Ps = []
        self.Ks = []   
        self.outers = []
        self.oml0 = None
        self.faces = []
        self.variables = {}
        self.params = {}
        self.Cv0 = 2
        self.Cv1 = 2

    def computeLeft(self):
        for f in range(len(self.Ks)):
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    #print self.oml0.visible[self.Ks[f][i,j]]
                    try:
                        self.oml0.visible[self.Ks[f][i,j]] = self.outers[f][i,j]
                    except:
                        self.oml0.visible[self.Ks[f][i,j]] = True
                    #try:
                    #    self.oml0.visible[self.Ks[f][i,j]] = avg > 0 and self.outers[f][i,j]
                    #except:
                    #    self.oml0.visible[self.Ks[f][i,j]] = avg > 0

    def setm(self, f, d, val, ind=[], both=True):
        self.set(1, f, d, val, ind)
        if both:
            self.set(2, f, d, [x*3 for x in val], ind)

    def setn(self, f, d, val, ind=[], both=True):
        self.set(2, f, d, val, ind)
        if both:
            self.set(1, f, d, [max(4,int(x/3.0)) for x in val], ind)

    def set(self, p, f, d, val, ind):
        ind = range(self.getms(f,d).shape[0]) if len(ind)==0 else ind
        val = [val[0] for x in range(len(ind))] if len(val) < len(ind) else val

        Ks = self.Ks
        ilist = range(Ks[f].shape[0]) if d==1 else ind
        jlist = range(Ks[f].shape[1]) if d==0 else ind
        for i in range(len(ilist)):
            for j in range(len(jlist)):
                surf = Ks[f][ilist[i],jlist[j]]
                pos = i if d==0 else j
                if not surf==-1:
                    self.oml0.edgeProperty(surf,p,d,val[pos])

    def computeVs(self):
        vs = self.variables
        ps = self.params
        for v in vs:
            vs[v][:,:] = 0.0
        for p in ps:
            vs[ps[p].var][:,:] += ps[p].compute()

    def addParam(self, name, var, shp, P=None, T=None, Tdim=0, D=None, Ddim=0, B=None, Bdim=0):
        self.params[name] = Parameter(var, shp, self.variables[var].shape, P, T, Tdim, D, Ddim, B, Bdim)
        
    def addFace(self, du, dv, d, ru=0.5, rv=0.5):
        """ Creates a set of rectangular surfaces, their IDs, and face dims.
        nu,nv: number of surfaces in the u and v directions
        du,dv: {1,2,3} maps to {x,y,z}; negative sign means reverse order
        d: position of the surfaces in the remaining coordinate axis
        ru,rv: surfaces span -ru to +ru in u dir. and -rv to +rv in v dir.

        Adds to self.Ps and self.Ks
        """
        self.faces.append([du,dv])
        nP = 10
        ni = self.ms[abs(du)-1].shape[0]
        nj = self.ms[abs(dv)-1].shape[0]
        verts = numpy.zeros((2,2,3),order='F')
        verts[:,:,:] = d
        verts[0,:,abs(du)-1] = -ru*numpy.sign(du)
        verts[1,:,abs(du)-1] = ru*numpy.sign(du)
        verts[:,0,abs(dv)-1] = -rv*numpy.sign(dv)
        verts[:,1,abs(dv)-1] = rv*numpy.sign(dv)
        for j in range(nj):
            for i in range(ni):
                self.Ps.append(PGMlib.bilinearinterp(nP, ni, nj, i+1, j+1, verts))

        if len(self.Ks) > 0:
            counter = numpy.max(self.Ks[-1]) + 1
        else:
            counter = 0
        K = numpy.zeros((ni,nj),int)
        for j in range(nj):
            for i in range(ni):
                K[i,j] = counter
                counter += 1
        self.Ks.append(K)
            
        self.outers.append(numpy.ones((ni,nj),bool))

    def averageEdges(self, edge1, edge2):
        avg = 0.5*edge1 + 0.5*edge2
        edge1[:,:] = avg
        edge2[:,:] = avg

    def connectEdges(self, f1=0, u1=None, v1=None, f2=0, u2=None, v2=None):
        def edge(f, u, v, kk):
            Ks = self.Ks
            Ps = self.Ps
            d = 0 if u==None else 1
            r = self.faces[f][d]
            k = kk if r > 0 else -1-kk
            surf = Ks[f][k,v] if d==0 else Ks[f][u,k]
            P = Ps[surf][:,v] if d == 0 else Ps[surf][u,:]
            return P[::-1] if r < 0 else P                
            
        for k in range(self.Ks[f1].shape[v1==None]):
            self.averageEdges(edge(f1,u1,v1,k), edge(f2,u2,v2,k))

    def setC1(self, t, f, i=None, j=None, u=None, v=None, d=None, val=True):
        """Set C1 continuity 

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

    def computeEdgeInfo(self):
        edgeProperty = self.oml0.edgeProperty
        Ks = self.Ks
        for f in range(len(Ks)):
            for i in range(Ks[f].shape[0]):
                for j in range(Ks[f].shape[1]):
                    surf = Ks[f][i,j]
                    if surf != -1:
                        self.getms(f,0)[i] = edgeProperty(surf,1)[0] - 1
                        self.getms(f,1)[j] = edgeProperty(surf,1)[1] - 1
                        self.getns(f,0)[i] = edgeProperty(surf,2)[0] - 1
                        self.getns(f,1)[j] = edgeProperty(surf,2)[1] - 1

    def getms(self, f, d):
        dim = self.faces[f][d]
        if dim > 0:
            return self.ms[abs(dim)-1]
        else:
            return self.ms[abs(dim)-1][::-1]

    def getns(self, f, d):
        dim = self.faces[f][d]
        if dim > 0:
            return self.ns[abs(dim)-1]
        else:
            return self.ns[abs(dim)-1][::-1]

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
                        mu,mv = oml0.edgeProperty(surf,1)
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
                                    Ns[f][ii,jj,0] = oml0.getIndex(surf,u,v,2)
                                    Ns[f][ii,jj,1] = i
                                    Ns[f][ii,jj,2] = u
                                    Ns[f][ii,jj,3] = j
                                    Ns[f][ii,jj,4] = v
        self.Qs = Qs
        self.Ns = Ns

    def propagateQs(self):
        Qs = self.Qs
        Ns = self.Ns
        oml0 = self.oml0
        for f in range(len(Ns)):
            PGMlib.updateqs(oml0.nQ, Ns[f].shape[0], Ns[f].shape[1], oml0.nvar, Ns[f], Qs[f], oml0.Q)

    def translatePoints(self, dx, dy, dz):
        for k in range(len(self.Ps)):
            self.Ps[k][:,:,0] += dx
            self.Ps[k][:,:,1] += dy
            self.Ps[k][:,:,2] += dz
