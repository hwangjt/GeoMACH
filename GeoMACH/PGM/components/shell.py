from __future__ import division
import numpy, time

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Primitive


class Shell(Primitive):
    """ A component used to model hollow bodies. """

    def __init__(self, nx=1, ny=1, nz=1, bottom=2):
        """ Initialization method
        nx: integer
            Number of surfaces in x direction
        ny: integer
            Number of surfaces in y direction
        nz: integer
            Number of surfaces in z direction
        bottom: integer
            Bottom of the shell
            0: open, C0
            1: open, C1
            2: closed
        """ 

        super(Shell,self).__init__(nx,ny,nz)

        self.addFace('rt0', 2, 1, -0.5)
        self.addFace('tp0', 3, 1, 0.5)
        self.addFace('lt0', -2, 1, 0.5)
        self.addFace('lt1', 2, 1, 0.4, 0.4, 0.4)
        self.addFace('tp1', -3, 1, 0.4, 0.4, 0.4)
        self.addFace('rt1', -2, 1, -0.4, 0.4, 0.4)
        if bottom==2:
            self.addFace('bt0', -3, 1, -0.5)
            self.addFace('bt1', 3, 1, -0.4, 0.4, 0.4)
        self.connectEdges(f1=0,v1= 0,f2=5,v2= 0)
        self.connectEdges(f1=0,v1=-1,f2=5,v2=-1)
        self.connectEdges(f1=1,v1= 0,f2=4,v2= 0)
        self.connectEdges(f1=1,v1=-1,f2=4,v2=-1)
        self.connectEdges(f1=2,v1= 0,f2=3,v2= 0)
        self.connectEdges(f1=2,v1=-1,f2=3,v2=-1)
        if bottom==2:
            self.connectEdges(f1=6,v1= 0,f2=7,v2= 0)
            self.connectEdges(f1=6,v1=-1,f2=7,v2=-1)

        self.bottom = bottom
        self.ax1 = 3
        self.ax2 = 1

    def setDOFs(self):
        faces = self.faces
        for face in faces.values():
            face.setC1('surf', val=True)
        if self.bottom==0:
            faces['rt0'].setC1('surf', i= 0, u= 0, val=False)
            faces['rt0'].setC1('edge', i= 0, u= 0, val=True)
            faces['lt0'].setC1('surf', i=-1, u=-1, val=False)
            faces['lt0'].setC1('edge', i=-1, u=-1, val=True)

    def declare_properties(self):
        super(Shell, self).declare_properties()
        n = self.faces['rt0'].num_cp[1]
        self.properties['thk'] = [n,3]

    def computeQs(self):
        faces = self.faces
        nx = faces['rt0'].num_cp[1]
        ny = faces['rt0'].num_cp[0]
        nz = faces['tp0'].num_cp[0]
        p = self.properties
        b = self.bottom==2

        r0 = p['scl'] + p['thk']/2.0
        r1 = p['scl'] - p['thk']/2.0

        shapes = range(8)
        shapes[0] = PGMlib.computeshape(ny, nx,-b/4.0, 1/4.0, p['flt'], p['shp','rt0'])
        shapes[1] = PGMlib.computeshape(nz, nx, 1/4.0, 3/4.0, p['flt'], p['shp','tp0'])
        shapes[2] = PGMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, p['flt'], p['shp','lt0'])
        shapes[6] = PGMlib.computeshape(nz, nx, 5/4.0, 7/4.0, p['flt'], p['shp','bt0'])
        shapes[5] = PGMlib.computeshape(ny, nx, 1/4.0,-b/4.0, p['flt'], p['shp','rt1'])
        shapes[4] = PGMlib.computeshape(nz, nx, 3/4.0, 1/4.0, p['flt'], p['shp','tp1'])
        shapes[3] = PGMlib.computeshape(ny, nx, (4+b)/4.0, 3/4.0, p['flt'], p['shp','lt1'])
        shapes[7] = PGMlib.computeshape(nz, nx, 7/4.0, 5/4.0, p['flt'], p['shp','bt1'])
        
        nQ = nx*(9+12*ny+12*nz) if self.bottom==2 else nx*(9+12*ny+6*nz)
        radii = [r0,r0,r0,r1,r1,r1,r0,r1]
        self.computeSections(nQ, shapes, radii=radii)
