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

        self.addFace( 2, 1,-0.5)
        self.addFace( 3, 1, 0.5)
        self.addFace(-2, 1, 0.5)
        if bottom==2:
            self.addFace(-3, 1,-0.5)

        self.bottom = bottom
        self.ax1 = 3
        self.ax2 = 1

    def setDOFs(self):
        faces = self.faces
        for f in range(len(faces)):
            faces[f].setC1('surf', val=True)
        if self.bottom==0:
            faces[0].setC1('surf', i= 0, u= 0, val=False)
            faces[0].setC1('edge', i= 0, u= 0, val=True)
            faces[2].setC1('surf', i=-1, u=-1, val=False)
            faces[2].setC1('edge', i=-1, u=-1, val=True)

    def initializeVariables(self):
        super(Body,self).initializeVariables()
        faces = self.faces
        nx = faces[0].num_cp[1]
        ny = faces[0].num_cp[0]
        nz = faces[1].num_cp[0]
        zeros = numpy.zeros
        v = self.variables
        a = self.addParam

        v['shR'] = zeros((ny,nx),order='F')
        v['shT'] = zeros((nz,nx),order='F')
        v['shL'] = zeros((ny,nx),order='F')
        v['shB'] = zeros((nz,nx),order='F')

        a('shR','shR',(1,1),P=[0.0])
        a('shT','shT',(1,1),P=[0.0])
        a('shL','shL',(1,1),P=[0.0])
        a('shB','shB',(1,1),P=[0.0])
        self.params['pos'].setP(P=[[0.,0.,0.],[1.,0.,0.]])

    def computeQs(self):
        faces = self.faces
        nx = faces[0].num_cp[1]
        ny = faces[0].num_cp[0]
        nz = faces[1].num_cp[0]
        v = self.variables
        b = self.bottom==2

        #v['pos'][0] = 2*v['pos'][1] - v['pos'][2]
        #v['pos'][-1] = 2*v['pos'][-2] - v['pos'][-3]

        shapes = range(4)
        shapes[0] = PGMlib.computeshape(ny, nx,-b/4.0, 1/4.0, v['flt'], v['shR'])
        shapes[1] = PGMlib.computeshape(nz, nx, 1/4.0, 3/4.0, v['flt'], v['shT'])
        shapes[2] = PGMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, v['flt'], v['shL'])
        shapes[3] = PGMlib.computeshape(nz, nx, 5/4.0, 7/4.0, v['flt'], v['shB'])

        nQ = nx*(9+6*ny+6*nz) if self.bottom==2 else nx*(9+6*ny+3*nz)
        self.computeSections(nQ, shapes)
