from __future__ import division
import numpy
import time

from GeoMACH.PUBS import PUBS

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


class Configuration2(object):
    
    def __init__(self):
        self.comps = {}
        self.keys = []

    def addComp(self, name, comp):
        self.comps[name] = comp
        self.keys.append(name)

    def separateComps(self):
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].translatePoints(k*4,0,0)   

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
        self.export = PUBS.PUBSexport(self.oml0)

        #self.oml0.plotm(mlab.figure())
        #self.oml0.write2Tec('test')
        #mlab.show()
        #exit()

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.oml0 = self.oml0
            comp.members = {}
            comp.keys = []
            comp.computeDims(self)
            comp.setDOFs()
        self.oml0.updateBsplines(True)

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.initializeDOFs()
            comp.initializeParameters()
            comp.propagateQs()
            comp.updateQs()
        self.oml0.computePoints()

    def updateComponents(self):
        self.oml0.updateBsplines()
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.computeDims(self)
            comp.initializeDOFs()
        self.computePoints()

    def computePoints(self):
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.propagateQs()
            comp.updateQs()
        self.oml0.computePoints()

    def buildStructure(self):
        ABMs = []        
        Ss = []
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            if not comp.keys==[]:
                print 'Building structure for',self.keys[k]
                comp.buildStructure()     
                ABMs.append(comp.strABM)
                Ss.append(comp.strS)
        self.ABM = self.oml0.vstackSparse(ABMs)
        self.S = numpy.vstack(Ss)
        self.computePoints()

    def writeStructure(self, name):
        S = self.S
        P = self.ABM.dot(self.oml0.Q)
        f = open(name+'_str.dat','w')
        f.write('title = "PUBSlib output"\n')
        f.write('variables = "x", "y", "z"\n')
        iP = 0
        for surf in range(S.shape[0]):    
            nu = int(S[surf,0])
            nv = int(S[surf,1])
            f.write('zone i='+str(nu)+', j='+str(nv)+', DATAPACKING=POINT\n')
            for v in range(nv):
                for u in range(nu):
                    f.write(str(P[iP,0]) + ' ' + str(P[iP,1]) + ' ' + str(P[iP,2]) + '\n')
                    iP += 1
        f.close()

    def plot(self):
        #self.oml0.plot(pylab.figure(),False)
        #pylab.show()
        self.oml0.plotm(mlab.figure(),False)
        mlab.show()



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
