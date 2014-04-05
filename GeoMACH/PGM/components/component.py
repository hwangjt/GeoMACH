from __future__ import division
import numpy
import scipy.sparse

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Parameter


class Component(object):
    """ Base class for Wing, Body, and Junction components. """

    def __init__(self):
        self.name = None
        self.Ps = []
        self.outers = []
        self.oml0 = None
        self.faces = []
        self.variables = {}
        self.params = {}
        self.Cv0 = 2
        self.Cv1 = 2

    def set_oml(self, oml):
        self.oml0 = oml
        for face in self.faces:
            face.oml = oml

    def computeVs(self):
        vs = self.variables
        ps = self.params
        for v in vs:
            vs[v][:,:] = 0.0
        for p in ps:
            vs[ps[p].var][:,:] += ps[p].compute()

    def addParam(self, name, var, shp, P=None, T=None, Tdim=0, D=None, Ddim=0, B=None, Bdim=0):
        self.params[name] = Parameter(var, shp, self.variables[var].shape, P, T, Tdim, D, Ddim, B, Bdim)
        
    def addFace(self, axis_u, axis_v, d, ru=0.5, rv=0.5):
        """ Creates a set of rectangular surfaces, their IDs, and face dims.
        nu,nv: number of surfaces in the u and v directions
        axis_u,axis_v: {1,2,3} maps to {x,y,z}; negative sign means reverse order
        d: position of the surfaces in the remaining coordinate axis
        ru,rv: surfaces span -ru to +ru in u dir. and -rv to +rv in v dir.
        """
        nP = 10
        ni = self.ms[abs(axis_u)-1].shape[0]
        nj = self.ms[abs(axis_v)-1].shape[0]
        verts = numpy.zeros((2,2,3),order='F')
        verts[:,:,:] = d
        verts[0,:,abs(axis_u)-1] = -ru*numpy.sign(axis_u)
        verts[1,:,abs(axis_u)-1] = ru*numpy.sign(axis_u)
        verts[:,0,abs(axis_v)-1] = -rv*numpy.sign(axis_v)
        verts[:,1,abs(axis_v)-1] = rv*numpy.sign(axis_v)
        for j in range(nj):
            for i in range(ni):
                self.Ps.append(PGMlib.bilinearinterp(nP, ni, nj, i+1, j+1, verts))

        if len(self.faces) > 0:
            counter = numpy.max(self.faces[-1].surf_indices) + 1
        else:
            counter = 0

        face = Face(len(self.faces), axis_u, axis_v, ni, nj)
        for j in range(nj):
            for i in range(ni):
                face.surf_indices[i,j] = counter
                counter += 1
        self.faces.append(face)
            
        self.outers.append(numpy.ones((ni,nj),bool))

    def averageEdges(self, edge1, edge2):
        avg = 0.5*edge1 + 0.5*edge2
        edge1[:,:] = avg
        edge2[:,:] = avg

    def connectEdges(self, f1=0, u1=None, v1=None, f2=0, u2=None, v2=None):
        def edge(f, u, v, kk):
            surf_indices = self.faces[f].surf_indices
            Ps = self.Ps
            d = 0 if u==None else 1
            r = self.faces[f].axes[d]
            k = kk if r > 0 else -1-kk
            surf = surf_indices[k,v] if d==0 else surf_indices[u,k]
            P = Ps[surf][:,v] if d == 0 else Ps[surf][u,:]
            return P[::-1] if r < 0 else P                
            
        for k in range(self.faces[f1].num_surf[v1==None]):
            self.averageEdges(edge(f1,u1,v1,k), edge(f2,u2,v2,k))

    def removeHiddenDOFs(self):
        pass


class Face(object):

    def __init__(self, num, axis_u, axis_v, ni, nj):
        self.num = num
        self.axes = [axis_u, axis_v]
        self.surf_indices = numpy.zeros((ni,nj),int)
        self.num_cp_list = [numpy.zeros(ni, int), numpy.zeros(nj, int)]
        self.num_pt_list = [numpy.zeros(ni, int), numpy.zeros(nj, int)]
        self.num_surf = [ni, nj]

        self.oml = None
        self.cp_array = None
        self.index_array = None
        self.num_cp = [None, None]

    def initializeDOFmappings(self, cp_vec, index_vec):
        def classify(i, n):
            if i==0:
                return 0
            elif i==n-1:
                return 2
            else:
                return 1

        def getC1(self, surf, u=None, v=None, d=0):
            if u==None:
                edge = self.oml.surf_edge[surf,0,v]
            else:
                edge = self.oml.surf_edge[surf,1,u]
            if edge > 0:
                return self.oml.edge_c1[abs(edge)-1,d]
            else:
                return self.oml.edge_c1[abs(edge)-1,1-abs(d)]

        oml = self.oml
        ni, nj = self.num_cp_list[:]
        self.cp_array = cp_vec.reshape((sum(ni)+1,sum(nj)+1,3), order='C')
        self.index_array = index_vec.reshape((sum(ni)+1,sum(nj)+1), order='C')
        for j in range(self.num_surf[1]):
            for i in range(self.num_surf[0]):
                surf = self.surf_indices[i,j]
                if surf != -1:
                    mu,mv = oml.edgeProperty(surf,1)
                    for v in range(mv):
                        jj = sum(nj[:j]) + v
                        vType = classify(v,mv)
                        for u in range(mu):
                            ii = sum(ni[:i]) + u
                            uType = classify(u,mu)
                            DOF = True
                            if uType==0 or uType==2 or vType==0 or vType==2:
                                DOF = DOF and not oml.surf_c1[surf,uType,vType]
                                if (not uType==1) and (not vType==1):
                                    DOF = DOF and not getC1(self,surf,u=int(uType/2),d=int(vType/2))
                                    DOF = DOF and not getC1(self,surf,v=int(vType/2),d=int(uType/2))
                            if DOF:
                                self.index_array[ii,jj] = oml.getIndex(surf,u,v,2)

    def compute_num_cp(self):
        edgeProperty = self.oml.edgeProperty
        for j in range(self.num_surf[1]):
            for i in range(self.num_surf[0]):
                surf = self.surf_indices[i,j]
                if surf != -1:
                    self.num_cp_list[0][i] = edgeProperty(surf,1)[0] - 1
                    self.num_cp_list[1][j] = edgeProperty(surf,1)[1] - 1
                    self.num_pt_list[0][i] = edgeProperty(surf,2)[0] - 1
                    self.num_pt_list[1][j] = edgeProperty(surf,2)[1] - 1
        for d in xrange(2):
            self.num_cp[d] = sum(self.num_cp_list[d]) + 1

    def propagateQs(self):
        PGMlib.updateqs(self.oml.nQ, self.num_cp[0], self.num_cp[1], self.oml.nvar, 
                        self.index_array, self.cp_array, self.oml.Q)

    def setm(self, d, val, ind=[], both=True):
        self.set(1, d, val, ind)
        if both:
            self.set(2, d, [x*3 for x in val], ind)

    def setn(self, d, val, ind=[], both=True):
        self.set(2, d, val, ind)
        if both:
            self.set(1, d, [max(4,int(x/3.0)) for x in val], ind)

    def set(self, p, d, val, ind):
        ind = range(self.num_surf[d]) if len(ind)==0 else ind
        val = [val[0] for x in range(len(ind))] if len(val) < len(ind) else val

        ilist = range(self.num_surf[0]) if d==1 else ind
        jlist = range(self.num_surf[1]) if d==0 else ind
        for i in range(len(ilist)):
            for j in range(len(jlist)):
                surf = self.surf_indices[ilist[i],jlist[j]]
                pos = i if d==0 else j
                if not surf==-1:
                    self.oml.edgeProperty(surf,p,d,val[pos])

    def setC1(self, t, i=None, j=None, u=None, v=None, d=None, val=True):
        """Set C1 continuity 

        t: {string} surface or edge C1

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

        def setSurfC1(i, j, u, v, d, val):
            oml = self.oml
            surf = self.surf_indices[i,j]
            if not surf==-1:
                if u==None and v==None:
                    oml.surf_c1[surf,:,:] = val                    
                elif u==None:
                    oml.surf_c1[surf,:,v] = val
                elif v==None:
                    oml.surf_c1[surf,u,:] = val
                else:
                    oml.surf_c1[surf,u,v] = val


        def setEdgeC1(i, j, u, v, d, val):
            oml = self.oml
            surf = self.surf_indices[i,j]
            if not surf==-1:
                if u==None:
                    edge = oml.surf_edge[surf,0,v]
                else:
                    edge = oml.surf_edge[surf,1,u]
                if d==None:
                    oml.edge_c1[abs(edge)-1,:] = val
                elif edge>0:
                    oml.edge_c1[abs(edge)-1,d] = val
                else:
                    oml.edge_c1[abs(edge)-1,1-abs(d)] = val

        if t=='surf':
            func = setSurfC1
        elif t=='edge':
            func = setEdgeC1
        if (not i==None) and (not j==None):
            func(i, j, u, v, d, val)
        elif not i==None:
            for j in range(self.num_surf[1]):
                func(i, j, u, v, d, val)
        elif not j==None:
            for i in range(self.num_surf[0]):
                func(i, j, u, v, d, val)
        else:
            for j in range(self.num_surf[1]):
                for i in range(self.num_surf[0]):
                    func(i, j, u, v, d, val)

    def setCornerC1(self, i=0, j=0, val=True):
        self.setC1('edge', i=i, j=j, u=i, d=j, val=val)
        self.setC1('edge', i=i, j=j, v=j, d=i, val=val)
