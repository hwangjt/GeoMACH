from __future__ import division
import numpy
import time

from GeoMACH.PUBS import PUBS, PUBSlib
from GeoMACH.PSM import PSMlib
from GeoMACH.PGM.configurations import Configuration

try: 
    from pyV3D.sender import WV_Sender
    from openmdao.main.interfaces import IParametricGeometry, implements, IStaticGeometry
    USE_OPENDMAO = True

except ImportError: 
    USE_OPENDMAO = False


class ConfigurationOpenMDAO(Configuration):
    ''' Implements OpenMDAO geometry API for Configuration class '''

    if USE_OPENDMAO: 
        implements(IParametricGeometry)

    def __init__(self):
        super(ConfigurationOpenMDAO, self).__init__()
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
