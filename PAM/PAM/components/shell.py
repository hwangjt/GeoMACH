from __future__ import division
from PAM.components import Primitive, Variable
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib


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

        self.addFace( 2, 1,-0.5)
        self.addFace( 3, 1, 0.5)
        self.addFace(-2, 1, 0.5)
        self.addFace( 2, 1, 0.4,0.4,0.4)
        self.addFace(-3, 1, 0.4,0.4,0.4)
        self.addFace(-2, 1,-0.4,0.4,0.4)
        if bottom==2:
            self.addFace(-3, 1,-0.5)
            self.addFace( 3, 1,-0.4,0.4,0.4)
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
        setC1 = self.setC1
        setCornerC1 = self.setCornerC1
        for f in range(len(self.Ks)):
            self.setC1('surf', f, val=True)
        if self.bottom==0:
            setC1('surf', 0, i= 0, u= 0, val=False)
            setC1('edge', 0, i= 0, u= 0, val=True)
            setC1('surf', 2, i=-1, u=-1, val=False)
            setC1('edge', 2, i=-1, u=-1, val=True)

    def initializeVariables(self):
        super(Shell,self).initializeVariables()
        nx = self.Qs[2].shape[1]
        ny = self.Qs[2].shape[0]
        nz = self.Qs[3].shape[0]
        v = self.variables
        v['shapeR0'] = Variable((ny,nx))
        v['shapeT0'] = Variable((nz,nx))
        v['shapeL0'] = Variable((ny,nx))
        v['shapeB0'] = Variable((nz,nx))
        v['shapeR1'] = Variable((ny,nx))
        v['shapeT1'] = Variable((nz,nx))
        v['shapeL1'] = Variable((ny,nx))
        v['shapeB1'] = Variable((nz,nx))
        v['thickness'] = Variable((nx,3))
        v['fillet'] = Variable((nx,4),False)
        self.setSections()

    def computeQs(self):
        val = lambda var: self.variables[var]()
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        b = self.bottom==2

        r0 = val('scl') + val('thickness')/2.0
        r1 = val('scl') - val('thickness')/2.0

        shapes = range(8)
        shapes[0] = PAMlib.computeshape(ny, nx,-b/4.0, 1/4.0, val('fillet'), val('shapeR0'))
        shapes[1] = PAMlib.computeshape(nz, nx, 1/4.0, 3/4.0, val('fillet'), val('shapeT0'))
        shapes[2] = PAMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, val('fillet'), val('shapeL0'))
        shapes[6] = PAMlib.computeshape(nz, nx, 5/4.0, 7/4.0, val('fillet'), val('shapeB0'))
        shapes[5] = PAMlib.computeshape(ny, nx, 1/4.0,-b/4.0, val('fillet'), val('shapeR1'))
        shapes[4] = PAMlib.computeshape(nz, nx, 3/4.0, 1/4.0, val('fillet'), val('shapeT1'))
        shapes[3] = PAMlib.computeshape(ny, nx, (4+b)/4.0, 3/4.0, val('fillet'), val('shapeL1'))
        shapes[7] = PAMlib.computeshape(nz, nx, 7/4.0, 5/4.0, val('fillet'), val('shapeB1'))
        
        nQ = nx*(9+12*ny+12*nz) if self.bottom==2 else nx*(9+12*ny+6*nz)
        radii = [r0,r0,r0,r1,r1,r1,r0,r1]
        self.computeSections(nQ, shapes, radii=radii)


if __name__ == '__main__':
    s = Shell(nx=4,ny=4,nz=4,bottom=2)
    import PUBS
    from mayavi import mlab
    s.oml0 = PUBS.PUBS(s.Ps)
    s.setDOFs()
    s.oml0.updateBsplines()
    s.computems()
    s.initializeDOFmappings()
    s.initializeVariables()
    s.variables['pos'][:,0] = numpy.linspace(0,4,s.Qs[0].shape[1])
    #s.variables['pos'][:,1] = numpy.linspace(0,4,s.Qs[0].shape[1])
    s.variables['pos'][:,2] = numpy.linspace(0,4,s.Qs[0].shape[1])
    #s.variables['pos'][4,1] = -0.5
    #s.variables['fillet'][:,0] = 0.4
    #s.variables['fillet'][:,1] = 0.6
    #s.variables['shapeR'][:10,3:-3] = -0.5
    #s.variables['shapeL'][-10:,3:-3] = -0.5
    s.propagateQs()
    s.updateQs()
    s.oml0.computePoints()
    s.oml0.plot()
    name='shell'
    s.oml0.write2Tec(name)
