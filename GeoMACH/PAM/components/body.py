from __future__ import division
import numpy

from GeoMACH.PAM import PAMlib
from GeoMACH.PAM.components import Primitive


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
        setC1 = self.setC1
        for f in range(len(self.Ks)):
            self.setC1('surf', f, val=True)
        if self.bottom==0:
            setC1('surf', 0, i= 0, u= 0, val=False)
            setC1('edge', 0, i= 0, u= 0, val=True)
            setC1('surf', 2, i=-1, u=-1, val=False)
            setC1('edge', 2, i=-1, u=-1, val=True)

    def initializeVariables(self):
        super(Body,self).initializeVariables()
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
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

        self.setSections()

    def computeQs(self):
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        v = self.variables
        b = self.bottom==2

        #v['pos'][0] = 2*v['pos'][1] - v['pos'][2]
        #v['pos'][-1] = 2*v['pos'][-2] - v['pos'][-3]

        shapes = range(4)
        shapes[0] = PAMlib.computeshape(ny, nx,-b/4.0, 1/4.0, v['flt'], v['shR'])
        shapes[1] = PAMlib.computeshape(nz, nx, 1/4.0, 3/4.0, v['flt'], v['shT'])
        shapes[2] = PAMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, v['flt'], v['shL'])
        shapes[3] = PAMlib.computeshape(nz, nx, 5/4.0, 7/4.0, v['flt'], v['shB'])

        nQ = nx*(9+6*ny+6*nz) if self.bottom==2 else nx*(9+6*ny+3*nz)
        self.computeSections(nQ, shapes)


if __name__ == '__main__':
    b = Body(nx=8,ny=4,nz=4,bottom=0)
    import PUBS
    from mayavi import mlab
    b.oml0 = PUBS.PUBS(b.Ps)
    b.setDOFs()
    b.oml0.updateBsplines()
    b.computems()
    b.initializeDOFmappings()
    b.initializeVariables()
    b.variables['pos'][:,0] = numpy.linspace(0,4,b.Qs[2].shape[1])
    #b.variables['pos'][:,1] = numpy.linspace(0,4,b.Qs[2].shape[1])
    b.variables['pos'][:,2] = numpy.linspace(0,4,b.Qs[2].shape[1])
    b.variables['pos'][4,1] = -0.5
    #b.variables['shapeL'][:,:] = 0.5
    #b.variables['shapeR'][:,:] = 0.5
    b.parameters['fillet'][:,0] = 0.5
    b.parameters['fillet'][:,1] = 0.5
    b.variables['noseL'] = 0.5
    b.variables['tailL'] = 0.5
    #b.variables['shapeR'][:10,3:-3] = -0.5
    #b.variables['shapeL'][-10:,3:-3] = -0.5
    b.propagateQs()
    b.updateQs()
    b.oml0.computePoints()
    b.oml0.plot()
    name='body'
    b.oml0.write2Tec(name)
