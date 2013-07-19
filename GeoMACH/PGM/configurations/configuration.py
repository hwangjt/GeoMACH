from __future__ import division
import numpy
import time

from GeoMACH.PUBS import PUBS, PUBSlib
from GeoMACH.PSM import PSMlib

try: 
    from pyV3D.sender import WV_Sender
    from openmdao.main.interfaces import IParametricGeometry, implements, IStaticGeometry
    USE_OPENDMAO = True

except ImportError: 
    USE_OPENDMAO = False

class Configuration(object):

    if USE_OPENDMAO: 
        implements(IParametricGeometry)
    
    def __init__(self):
        self.comps = {}
        self.keys = []
        self.inds = {}

        self._callbacks = []

    def regen_model(self):
        self.computePoints()

    def list_parameters(self):
        self.tris = self.oml0.exportPtri()

        params = []
        for k in range(len(self.comps)):
            c = self.keys[k]
            comp = self.comps[c]
            for p in comp.params.keys():
                par = comp.params[p]
                meta = {}
                meta['value'] = par.P[:,:,0]
                meta['iotype'] = 'in'
                meta['shape'] = par.P.shape[:2]
                params.append((c+'.'+p, meta))
        return params

    def set_parameter(self, name, value):
        c, p = name.split('.', 1)
        self.comps[c].params[p].setP(value[:])

    def get_parameters(self, names):
        vals = []
        for name in names:
            c, p = name.split('.', 1)
            val = self.comps[c].params[p].P[:,:,0]
            vals.append(val)
        return vals

    def get_static_geometry(self):
        g = GeoMACHGeometry()
        g.tris = self.tris
        g._model = self
        return g
    
    def register_param_list_changedCB(self, callback):
        """Register a callback that will be called when self.invoke_callbacks() is called.
        self.invoke_callbacks() should be called from the inheriting class whenever any
        parameters are added, removed, or change their type.
        """
        self._callbacks.append(callback)

    def invoke_callbacks(self):
        """Invokes any callbacks that have been registered via register_param_list_changedCB."""

        for cb in self._callbacks:
            cb()

    def get_attributes(self, io_only=True):
        """Return an attribute dict for use by the openmdao GUI.
        """
        
        return {
        }

    def addComp(self, name, comp):
        self.comps[name] = comp
        self.keys.append(name)
        self.inds[name] = len(self.inds)

    def separateComps(self):
        self.nprim = len(self.comps)
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].translatePoints(0,0,k*4)   

    def assembleComponents(self):
        Ps = []
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            if k > 0:
                comp0 = self.comps[self.keys[k-1]]
                maxk = numpy.max(comp0.Ks[-1]) + 1
                for s in range(len(comp.Ks)):
                    for j in range(comp.Ks[s].shape[1]):
                        for i in range(comp.Ks[s].shape[0]):
                            if comp.Ks[s][i,j] != -1:
                                comp.Ks[s][i,j] += maxk
            Ps.extend(comp.Ps)
            comp.Ps = []

        self.oml0 = PUBS.PUBS(Ps)

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.oml0 = self.oml0
            comp.setDOFs()
        self.oml0.updateBsplines()
        self.updateParametrization()

    def updateParametrization(self):
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.computeEdgeInfo()
            comp.initializeDOFmappings()
            comp.initializeVariables()
        self.computePoints()

    def update(self):
        self.oml0.update()
        self.updateParametrization()

    def computePoints(self):
        t0 = time.time()
        self.computeVs()
        self.computeQs()
        self.propagateQs()
        self.oml0.computePoints()
        print time.time()-t0

    def computeVs(self):
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].computeVs()

    def computeQs(self, full=True, comp=None):
        if full:
            for k in range(len(self.comps)):
                self.comps[self.keys[k]].computeQs()
        else:
            for k in range(self.nprim,len(self.comps)):
                if self.keys[k] != comp:
                    self.comps[self.keys[k]].computeQs()

    def propagateQs(self):
        self.oml0.Q[:,:3] = 0.0
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].propagateQs()

    def getDerivatives(self, c, p, ind, clean=True, FD=False, h=1e-5):
        comp = self.comps[c]
        par = comp.params[p]
        var = par.var
        self.computeVs()
        self.computeQs()
        self.propagateQs()
        V0 = numpy.array(comp.variables[var])
        Q0 = numpy.array(self.oml0.Q[:,:3])
        if FD:
            par.P[ind[0],ind[1],0] += h
            self.computeVs()
            self.computeQs()
            par.P[ind[0],ind[1],0] -= h
        else:
            h = 1.0
            par.P[ind[0],ind[1],0] += h
            self.computeVs()
            par.P[ind[0],ind[1],0] -= h
            dV = comp.variables[var] - V0
            self.computeVs()
            comp.setDerivatives(var,dV)
            self.computeQs(False, c)
        self.propagateQs()
        res = (self.oml0.Q[:,:3] - Q0)/h
        if clean:
            self.computePoints()
        return res

    def runDerivativeTest(self, c, ps=[]):
        self.computePoints()

        comp = self.comps[c]
        if ps==[]:
            ps = comp.params.keys()

        self.computePoints()
        h = 1e-5
        for p in ps:
            par = comp.params[p]
            var = par.var
            if not (var in ['nor','ogn','flt']):
                ni,nj = par.P.shape[:2]
                for i in range(ni):
                    for j in range(nj):
                        ind = (i,j)
                        t0 = time.time()
                        d1 = self.getDerivatives(c,p,ind,clean=False)
                        t1 = time.time()
                        d2 = self.getDerivatives(c,p,ind,clean=False,FD=True,h=h)
                        t2 = time.time()
                        norm0 = numpy.linalg.norm(d2)
                        norm0 = 1.0 if norm0==0 else norm0
                        error = numpy.linalg.norm(d2-d1)/norm0
                        good = 'O' if error < 1e-4 else 'X'
                        print good, ' ', c, ' ', p, ' ', ind, ' ', error #t1-t0, t2-t1
        self.computePoints()

    def getDerivatives0(self, comp, var, ind, clean=True, FD=False, h=1e-5):
        self.computeQs()
        self.propagateQs()
        Q0 = numpy.array(self.oml0.Q[:,:3])
        if FD:
            self.comps[comp].variables[var][ind] += h
            self.computeQs()
            self.comps[comp].variables[var][ind] -= h
        else:
            self.comps[comp].setDerivatives(var,ind)
            self.computeQs(False, comp)
            h = 1.0
        self.propagateQs()
        res = (self.oml0.Q[:,:3] - Q0)/h
        if clean:
            self.computePoints()
        return res

    def runDerivativeTest0(self, comp, variables=[]):
        self.computePoints()
        if variables==[]:
            variables = self.comps[comp].variables.keys()
        h = 1e-5
        for var in variables:
            if not (var in ['nor','origin','fillet']):
                dat = self.comps[comp].variables[var]
                for ind,x in numpy.ndenumerate(dat):
                    ind = ind[0] if len(ind)==1 else ind
                    t0 = time.time()
                    d1 = self.getDerivatives(comp,var,ind,clean=False)
                    t1 = time.time()
                    d2 = self.getDerivatives(comp,var,ind,clean=False,FD=True,h=h)
                    t2 = time.time()
                    norm0 = numpy.linalg.norm(d2)
                    norm0 = 1.0 if norm0==0 else norm0
                    error = numpy.linalg.norm(d2-d1)/norm0
                    good = 'O' if error < 1e-4 else 'X'
                    print good, ' ', comp, ' ', var, ' ', ind, ' ', error #t1-t0, t2-t1
        self.computePoints()

    def meshStructure(self, members, lengths):
        oml0 = self.oml0

        nmem = len(members)
        faces = -numpy.ones((nmem,4,2),int,order='F')
        coords = numpy.zeros((nmem,4,2,2,3),order='F')
        for imem in range(nmem):
            key = members.keys()[imem]
            for icontr in range(len(members[key])):
                faces[imem,icontr,0] = self.inds[members[key][icontr][0]]
                faces[imem,icontr,1] = members[key][icontr][1]
                coords[imem,icontr,0,0,:] = members[key][icontr][2]
                coords[imem,icontr,1,0,:] = members[key][icontr][3]
                coords[imem,icontr,0,1,:] = members[key][icontr][4]
                coords[imem,icontr,1,1,:] = members[key][icontr][5]

        surfs = numpy.linspace(0,oml0.nsurf-1,oml0.nsurf)
        s = numpy.zeros(4*oml0.nsurf)
        s[0::4] = surfs
        s[1::4] = surfs
        s[2::4] = surfs
        s[3::4] = surfs
        u = numpy.zeros(4*oml0.nsurf)
        u[1::4] = 1.
        u[3::4] = 1.
        v = numpy.zeros(4*oml0.nsurf)
        v[2::4] = 1.
        v[3::4] = 1.
        quadsS = numpy.zeros((oml0.nsurf,4),order='F')
        quadsS[:,0] = 4*surfs + 0
        quadsS[:,1] = 4*surfs + 1
        quadsS[:,2] = 4*surfs + 3
        quadsS[:,3] = 4*surfs + 2
        B = oml0.evaluateBases(s,u,v)
        nodesS = B.dot(oml0.C[:,:3])
        oml0.export.write2TecQuads('john.dat',nodesS,quadsS)

        #oml0.C[:,:] = -1.0
        #for k in range(len(self.comps)):
        #    c = self.keys[k]
        #    comp = self.comps[c]
        #    for f in range(len(comp.Ks)):
        #        ni, nj = comp.Ks[f].shape
        #        for i in range(ni):
        #            for j in range(nj):
        #                surf = comp.Ks[f][i,j]
        #                mu, mv = oml0.edgeProperty(surf,1)
        #                ugroup = oml0.edge_group[abs(oml0.surf_edge[surf,0,0])-1]
        #                vgroup = oml0.edge_group[abs(oml0.surf_edge[surf,1,0])-1]
        #                mu = oml0.group_m[ugroup-1]
        #                mv = oml0.group_m[vgroup-1]
        #                for u in range(mu):
        #                    for v in range(mv):
        #                        oml0.C[oml0.getIndex(surf,u,v,1),:3] = [u/(mu-1), v/(mv-1), 0]

        oml0.C[:,:] = -1.0
        for surf in range(oml0.nsurf):
            mu, mv = oml0.edgeProperty(surf,1)
            ugroup = oml0.edge_group[abs(oml0.surf_edge[surf,0,0])-1]
            vgroup = oml0.edge_group[abs(oml0.surf_edge[surf,1,0])-1]
            mu = oml0.group_m[ugroup-1]
            mv = oml0.group_m[vgroup-1]
            for u in range(mu):
                for v in range(mv):
                    oml0.C[oml0.getIndex(surf,u,v,1),:3] = [u/(mu-1), v/(mv-1), 0]
                    print u/(mu-1), v/(mv-1)

        oml0.computePointsC()
        #oml0.export.write2TecQuads('john2.dat',nodesM,quadsM)
        oml0.write2Tec('test2')
        oml0.write2TecC('test2')

        self.computePoints()
        print oml0.C[oml0.getIndex(147,-1, 0,1),:3]
        print oml0.C[oml0.getIndex(144,-1,-1,1),:3]
        print oml0.C[oml0.getIndex( 49, 0, 0,1),:3]


        s = numpy.zeros(4*nmem)
        P = numpy.zeros((4*nmem,3),order='F')
        Q = numpy.zeros((4*nmem,3),order='F')
        w = numpy.zeros((4*nmem,4),order='F')
        mems = numpy.linspace(0,nmem-1,nmem)
        Q[:,2] = 1.
        Bs = []
        ws = []
        for icontr in range(4):
            for imem in range(nmem):
                k = faces[imem,icontr,0]
                f = faces[imem,icontr,1]
                c = self.keys[k]
                comp = self.comps[c]
                ni, nj = comp.Ks[f].shape
                for i in range(2):
                    for j in range(2):
                        u, v = coords[imem,icontr,i,j,:2]
                        ii = int(numpy.floor(u*ni))
                        jj = int(numpy.floor(v*nj))
                        s[4*imem+2*j+i] = comp.Ks[f][ii,jj]
                        P[4*imem+2*j+i,0] = u*ni - ii
                        P[4*imem+2*j+i,1] = v*nj - jj
                        w[4*imem+2*j+i,icontr] = coords[imem,icontr,i,j,2]
            surf,u,v = oml0.evaluateProjection(P, Q=Q)
            Bs.append(oml0.evaluateBases(surf, u, v))


        quadsM = numpy.zeros((nmem,4),order='F')
        quadsM[:,0] = 4*mems + 0
        quadsM[:,1] = 4*mems + 1
        quadsM[:,2] = 4*mems + 3
        quadsM[:,3] = 4*mems + 2
        nodesM = numpy.zeros((4*nmem,3),order='F')
        for icontr in range(4):
            for k in range(3):
                nodesM[:,k] += w[:,icontr] * Bs[icontr].dot(oml0.C[:,k])
            
        oml0.export.write2TecQuads('john2.dat',nodesM,quadsM)

        exit()



        Ps = numpy.zeros((oml0.nsurf,3,3,3),order='F')
        for s in range(oml0.nsurf):
            for i in range(3):
                for j in range(3):
                    Ps[s,i,j] = oml0.evaluatePoint(s,i/2.0,j/2.0)[:3]
        nvertS,ngroupS,surf_vert,surf_group = PUBSlib.initializeconnectivities(oml0.nsurf,1e-13,1e-5,Ps)

        nvertM,ngroupM,mem_vert,mem_group = PSMlib.computemembertopology(nmem, faces, coords)
        mem_group[:,:,:] += ngroupS
        ngroup = ngroupS + ngroupM

        nint = PSMlib.countfaceintersections(nmem, coords)
        intFaces, intCoords = PSMlib.computefaceintersections(nmem, nint, faces, coords)

        groupIntCount = numpy.zeros(ngroup,int)
        meshesS = []
        for k in range(len(self.comps)):
            c = self.keys[k]
            comp = self.comps[c]
            meshes = []
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                nedge = PSMlib.countfaceedges(k, f, ni, nj, nint, intFaces)
                edges, edge_group = PSMlib.computefaceedges(k, f, ni, nj, oml0.nsurf, nint, nmem, nedge, comp.Ks[f], surf_group, mem_group, intFaces, intCoords)
                mesh = QuadMesh(0.2,edges)
                mesh.computeIntersections()
                mesh.computeDivisions()
                mesh.deleteDuplicateVerts()
                groupIntCount = PSMlib.countgroupintersections(mesh.verts.shape[0], mesh.edges.shape[0], ngroup, mesh.verts, mesh.edges, edge_group, groupIntCount)
                meshes.append([mesh,edge_group])
            meshesS.append(meshes)

        meshesM = []
        for imem in range(nmem):
            edges, edge_group = PSMlib.computememberedges(imem+1, nmem, mem_group)
            mesh = QuadMesh(1e6,edges)
            meshesM.append([mesh,edge_group])

        groupIntPtr = PSMlib.computegroupintptr(ngroup, groupIntCount)
        nint = groupIntPtr[-1,-1]
        groupInts = numpy.zeros(nint)
        for k in range(len(self.comps)):
            c = self.keys[k]
            comp = self.comps[c]
            for f in range(len(comp.Ks)):
                mesh, edge_group = meshesS[k][f]
                groupInts = PSMlib.computegroupintersections(mesh.verts.shape[0], mesh.edges.shape[0], ngroup, nint, mesh.verts, mesh.edges, edge_group, groupIntPtr, groupInts)
                #print groupInts

        for k in range(len(self.comps)):
            c = self.keys[k]
            comp = self.comps[c]
            for f in range(len(comp.Ks)):
                mesh, edge_group = meshesS[k][f]
                nvert = PSMlib.countintersectionverts(mesh.edges.shape[0], ngroup, edge_group, groupIntPtr)
                mesh.verts = PSMlib.computeintersectionverts(mesh.verts.shape[0], mesh.edges.shape[0], ngroup, nint, nvert + mesh.verts.shape[0], mesh.verts, mesh.edges, edge_group, groupIntPtr, groupInts)
                print 'QM1', k, f
                mesh.mesh()
                #edges[:,:,0] *= lengths[k,0]
                #edges[:,:,1] *= lengths[k,1]
                if k==1 and f==0 and 0:
                    import pylab
                    mesh.plot(111,pt=False,pq=False)
                    pylab.show()
                    exit()

        for imem in range(nmem):
            mesh, edge_group = meshesM[imem]
            nvert = PSMlib.countintersectionverts(mesh.edges.shape[0], ngroup, edge_group, groupIntPtr)
            mesh.verts = PSMlib.computeintersectionverts(mesh.verts.shape[0], mesh.edges.shape[0], ngroup, nint, nvert + mesh.verts.shape[0], mesh.verts, mesh.edges, edge_group, groupIntPtr, groupInts)
            print 'QM2', imem
            mesh.mesh()
            if 0:
                import pylab
                mesh.plot(111,pt=False,pq=False)
                pylab.show()
                exit()
            meshesM.append([mesh,edge_group])

        oml0.C[:,:] = -1.0
        for k in range(len(self.comps)):
            c = self.keys[k]
            comp = self.comps[c]
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                for i in range(ni):
                    for j in range(nj):
                        surf = comp.Ks[f][i,j]
                        mu, mv = oml0.edgeProperty(surf,1)
                        for u in range(mu):
                            uu = u/(mu-1)
                            for v in range(mv):
                                vv = v/(mv-1)
                                oml0.C[oml0.getIndex(surf,u,v,1),:3] = [(uu+i)/ni, (vv+j)/nj, 0]
        oml0.computePointsC()
        #oml0.write2Tec('test2')
        #oml0.write2TecC('test2')
        #exit()

        Bs = []
        quads = []
        nquad0 = 0
        for k in range(len(self.comps)):
            c = self.keys[k]
            comp = self.comps[c]
            for f in range(len(comp.Ks)):
                mesh, edge_group = meshesS[k][f]
                ni, nj = comp.Ks[f].shape
                print mesh.verts.shape[0]
                P0, surfs, Q = PSMlib.computeprojtninputs(mesh.verts.shape[0], ni, nj, mesh.verts, comp.Ks[f])
                surf,u,v = oml0.evaluateProjection(P0, comp.Ks[f].flatten(), Q)
                Bs.append(oml0.evaluateBases(surf,u,v))
                quads.append(nquad0 + mesh.quads - 1)
                #quads.append(mesh.quads - 1)
                nquad0 += mesh.verts.shape[0]

        self.computePoints()

        #i = 0
        #for k in range(len(self.comps)):
        #    c = self.keys[k]
        #    comp = self.comps[c]
        #    for f in range(len(comp.Ks)):
        #        P = Bs[i].dot(oml0.C[:,:3])
        #        oml0.export.write2TecQuads('data'+str(k)+'-'+str(f)+'.dat',P,quads[i]-1)
        #        i += 1

        import scipy.sparse
        B = scipy.sparse.vstack(Bs)
        P = B.dot(oml0.C[:,:3])
        quads = numpy.vstack(quads)
        oml0.export.write2TecQuads('john.dat',P,quads)#mesh.quads-1)

                

class GeoMACHGeometry(object):
    '''A wrapper for a GeoMACH object that respresents a specific instance of a
    geometry at a designated set of parameters. Parameters are not modifiable.
    This object is able to provide data for visualization.
    '''  

    if USE_OPENDMAO: 
        implements(IStaticGeometry)
        
    def get_visualization_data(self, wv):
        '''Fills the given WV_Wrapper object with data for faces,
        edges, colors, etc.
        
        wv: WV_Wrapper object
 
        '''
        if self._model is None:
            return []

        t0 = time.time()
        xyzs = self._model.oml0.P0[:,:3]
        tris = self.tris

        mins = numpy.min(xyzs, axis=0)
        maxs = numpy.max(xyzs, axis=0)

        box = [mins[0], mins[1], mins[2], maxs[0], maxs[1], maxs[2]]

        print 'xyz shape = %s' % list(xyzs.shape)
        #xyzs = xyzs.astype(numpy.float32).flatten(order='C')

        print 'len(tris) = ',len(tris)

#        wv.set_face_data(xyzs.astype(numpy.float32).flatten(order='C'), 
#                             tris.astype(numpy.int32).flatten(order='C'), bbox=box, name="oml_surf")

        for i,tri in enumerate(tris):
            min_idx = int(numpy.min(tri))
            max_idx = int(numpy.max(tri))
            #print "type = %s" % type(tri[0,0])
            #print 'i: %d    len(tri): %d' % (i,len(tri))
            new_a = xyzs[min_idx:max_idx+1]
            new_tri = tri - min_idx

            wv.set_face_data(new_a.astype(numpy.float32).flatten(order='C'), 
                             new_tri.astype(numpy.int32).flatten(order='C'), bbox=box, name="oml_surf%d" % i)
#        print 'T2', time.time()-t0

if USE_OPENDMAO: 
    class GeoMACHSender(WV_Sender):

        def initialize(self, **kwargs):
            eye    = numpy.array([0.0, 0.0, 7.0], dtype=numpy.float32)
            center = numpy.array([0.0, 0.0, 0.0], dtype=numpy.float32)
            up     = numpy.array([0.0, 1.0, 0.0], dtype=numpy.float32)
            fov   = 30.0
            zNear = 1.0
            zFar  = 10.0

            bias  = 0
            self.wv.createContext(bias, fov, zNear, zFar, eye, center, up)

        @staticmethod
        def supports(obj):
            return isinstance(obj, GeoMACHGeometry) or isinstance(obj, Configuration)

        def geom_from_obj(self, obj):
            if isinstance(obj, Configuration):
                obj = obj.get_geometry()
                if obj is None:
                    raise RuntimeError("can't get Geometry object from GeoMACHParametricGeometry")
            elif not isinstance(obj, GeoMACHGeometry):
                raise TypeError("object must be a GeoMACHParametricGeometry or GeoMACHGeometry but is a '%s' instead" %
                    str(type(obj)))
            obj.get_visualization_data(self.wv)
