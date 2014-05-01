from __future__ import division
import numpy
import scipy.sparse
from collections import OrderedDict

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Parameter


class Component(object):
    """ Base class for Wing, Body, and Junction components. """

    def __init__(self):
        self.name = None
        self.Ps = []
        self.outers = []
        self.oml0 = None
        self.faces = OrderedDict()
        self.props = {}
        self.params = {}
        self.Cv0 = 2
        self.Cv1 = 2

    def count_properties(self):
        self.size_prop = 0
        for prop in self.props.values():
            self.size_prop += prop.ni * prop.nj

    def declare_properties(self):
        self.shapes = {}
        for name in self.faces:
            ni, nj = self.faces[name].num_cp
            self.props['shN', name] = Property(ni, nj)
            self.props['shX', name] = Property(ni, nj)
            self.props['shY', name] = Property(ni, nj)
            self.props['shZ', name] = Property(ni, nj)
            self.shapes[name] = numpy.zeros((ni, nj, 3), order='F')

    def initialize_properties(self, prop_vec, prop_ind):
        start, end = 0, 0
        for prop in self.props.values():
            end += prop.ni * prop.nj
            prop.initialize_properties(prop_vec[start:end], 
                                       prop_ind[start:end])
            start += prop.ni * prop.nj

    def set_oml(self, oml):
        self.oml0 = oml
        for face in self.faces.values():
            face.oml = oml

    def computeVs(self):
        for prop in self.props.values():
            prop.compute()
        
    def addFace(self, name, axis_u, axis_v, d):
        """ Creates a set of rectangular surfaces, their IDs, and face dims.
        nu,nv: number of surfaces in the u and v directions
        axis_u,axis_v: {1,2,3} maps to {x,y,z}; negative sign means reverse order
        d: position of the surfaces in the remaining coordinate axis
        ru,rv: surfaces span -ru to +ru in u dir. and -rv to +rv in v dir.
        """
        ni = self.ms[abs(axis_u)-1].shape[0]
        nj = self.ms[abs(axis_v)-1].shape[0]

        if len(self.faces) > 0:
            counter = numpy.max(self.faces.values()[-1].surf_indices) + 1
        else:
            counter = 0

        self.faces[name] = Face(len(self.faces), axis_u, axis_v, ni, nj)
        for j in range(nj):
            for i in range(ni):
                self.faces[name].surf_indices[i,j] = counter
                counter += 1
            
        self.outers.append(numpy.ones((ni,nj),bool))

    def removeHiddenSurfaces(self):
        pass

    def compute_cp_wireframe(self):
        return numpy.zeros(0), numpy.zeros(0), numpy.zeros(0)

    def compute_cp_surfs(self):
        return numpy.zeros(0), numpy.zeros(0), numpy.zeros(0)


class Face(object):

    def __init__(self, num, axis_u, axis_v, ni, nj):
        self.num = num
        self.axes = [axis_u, axis_v]
        self.surf_indices = numpy.zeros((ni,nj),int)
        self.num_cp_list = [numpy.zeros(ni, int), numpy.zeros(nj, int)]
        self.num_pt_list = [numpy.zeros(ni, int), numpy.zeros(nj, int)]
        self.num_cp_list = [3*numpy.ones(ni, int), 3*numpy.ones(nj, int)]
        self.num_surf = [ni, nj]

        self.oml = None
        self.cp_array = None
        self.index_array = None
        self.num_cp = [None, None]
        for d in xrange(2):
            self.num_cp[d] = sum(self.num_cp_list[d]) + 1

    def initialize_cp_data(self, cp_vec, cp_indices, index_vec):
        num_cp = self.num_cp
        self.cp_array = cp_vec.reshape((num_cp[0],num_cp[1],3), order='C')
        self.index_array = index_vec.reshape((num_cp[0],num_cp[1],3), order='C')
        self.cp_indices = cp_indices.reshape((num_cp[0],num_cp[1],3), order='C')

    def initializeDOFmappings(self):
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
                                if oml.getIndex(surf,u,v,2) >= 0:
                                    self.index_array[ii,jj,0] = 3*oml.getIndex(surf,u,v,2) + 0
                                    self.index_array[ii,jj,1] = 3*oml.getIndex(surf,u,v,2) + 1
                                    self.index_array[ii,jj,2] = 3*oml.getIndex(surf,u,v,2) + 2

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



class Property(object):

    def __init__(self, ni, nj):
        self.ni = ni
        self.nj = nj
        self.prop_vec = None
        self.prop_ind = None
        self.params = {}

    def initialize_properties(self, prop_vec, prop_ind):
        ni, nj = self.ni, self.nj
        self.prop_vec = prop_vec.reshape((ni, nj), order='C')
        self.prop_ind = prop_ind.reshape((ni, nj), order='C')

    def compute(self):
        self.prop_vec[:,:] = 0.0
        for param in self.params.values():
            self.prop_vec[:,:] += param.compute()

    def addParam(self, name, shp, P=None, T=None, Tdim=0, D=None, Ddim=0, B=None, Bdim=0):
        self.params[name] = Parameter(shp, [self.ni, self.nj], P, T, Tdim, D, Ddim, B, Bdim)
