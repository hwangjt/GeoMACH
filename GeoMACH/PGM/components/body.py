from __future__ import division
import numpy

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Primitive


class Body(Primitive):
    """ A component used to model blunt bodies. """

    def __init__(self, nx=1, ny=1, nz=1, bottom=2):
        """ Initialization method
        nx: integer
            Number of surfaces in x direction
        ny: integer
            Number of surfaces in y direction
        nz: integer
            Number of surfaces in z direction
        bottom: integer
            Bottom face of the body
            0: open, C0
            1: open, C1
            2: closed
        """ 

        super(Body,self).__init__(nx,ny,nz)

        self.addFace('rgt', 2, 1, -0.5)
        self.addFace('top', 3, 1, 0.5)
        self.addFace('lft', -2, 1, 0.5)
        if bottom==2:
            self.addFace('bot', -3, 1, -0.5)

        self.bottom = bottom
        self.ax1 = 3
        self.ax2 = 1

    def setDOFs(self):
        faces = self.faces
        for face in self.faces.values():
            face.setC1('surf', val=True)
        if self.bottom==0:
            faces['rgt'].setC1('surf', i= 0, u= 0, val=False)
            faces['rgt'].setC1('edge', i= 0, u= 0, val=True)
            faces['lft'].setC1('surf', i=-1, u=-1, val=False)
            faces['lft'].setC1('edge', i=-1, u=-1, val=True)

    def declare_properties(self):
        super(Body, self).declare_properties()

    def computeQs(self):
        faces = self.faces
        nx = faces['rgt'].num_cp[1]
        ny = faces['rgt'].num_cp[0]
        nz = faces['top'].num_cp[0]
        p = self.properties
        b = self.bottom==2

        #p['pos'][0] = 2*p['pos'][1] - p['pos'][2]
        #p['pos'][-1] = 2*p['pos'][-2] - p['pos'][-3]

        shapes = range(4)
        shapes[0] = PGMlib.computeshape(ny, nx,-b/4.0, 1/4.0, p['flt'], p['shp','rgt'])
        shapes[1] = PGMlib.computeshape(nz, nx, 1/4.0, 3/4.0, p['flt'], p['shp','top'])
        shapes[2] = PGMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, p['flt'], p['shp','lft'])
        shapes[3] = PGMlib.computeshape(nz, nx, 5/4.0, 7/4.0, p['flt'], p['shp','bot'])

        nQ = nx*(9+6*ny+6*nz) if self.bottom==2 else nx*(9+6*ny+3*nz)
        self.computeSections(nQ, shapes)
