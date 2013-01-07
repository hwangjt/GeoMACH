from __future__ import division
from PAM.components import Primitive, Property
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
        self.variables['shapeR0'] = numpy.zeros((ny,nx),order='F')
        self.variables['shapeT0'] = numpy.zeros((nz,nx),order='F')
        self.variables['shapeL0'] = numpy.zeros((ny,nx),order='F')
        self.variables['shapeB0'] = numpy.zeros((nz,nx),order='F')
        self.variables['shapeR1'] = numpy.zeros((ny,nx),order='F')
        self.variables['shapeT1'] = numpy.zeros((nz,nx),order='F')
        self.variables['shapeL1'] = numpy.zeros((ny,nx),order='F')
        self.variables['shapeB1'] = numpy.zeros((nz,nx),order='F')
        self.variables['thickness'] = 0.1*numpy.ones((nx,3),order='F')
        self.parameters['fillet'] = numpy.zeros((nx,4),order='F')
        self.setSections()

    def computeQs(self):
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        v = self.variables
        p = self.parameters
        b = self.bottom==2

        rot, self.drot0_dpos = self.computeRotations()

        r0 = v['scale'] + v['thickness']
        r1 = v['scale'] - v['thickness']

        shapes = range(8)
        shapes[0] = PAMlib.computeshape(ny, nx,-b/4.0, 1/4.0, p['fillet'], v['shapeR0'])
        shapes[1] = PAMlib.computeshape(nz, nx, 1/4.0, 3/4.0, p['fillet'], v['shapeT0'])
        shapes[2] = PAMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, p['fillet'], v['shapeL0'])
        shapes[6] = PAMlib.computeshape(nz, nx, 5/4.0, 7/4.0, p['fillet'], v['shapeB0'])
        shapes[5] = PAMlib.computeshape(ny, nx, 1/4.0,-b/4.0, p['fillet'], v['shapeR1'])
        shapes[4] = PAMlib.computeshape(nz, nx, 3/4.0, 1/4.0, p['fillet'], v['shapeT1'])
        shapes[3] = PAMlib.computeshape(ny, nx, (4+b)/4.0, 3/4.0, p['fillet'], v['shapeL1'])
        shapes[7] = PAMlib.computeshape(nz, nx, 7/4.0, 5/4.0, p['fillet'], v['shapeB1'])

        if self.bottom==2:
            nQ = nx*(6+12*ny+12*nz)
        else:
            nQ = nx*(6+12*ny+6*nz)
        radii = [r0,r0,r0,r1,r1,r1,r0,r1]
        self.computeSections(nQ, rot, shapes,radii=radii)

    def setDerivatives(self, var, ind):
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        if var=='offset':
            for f in range(len(self.Qs)):
                self.Qs[f][:,:,ind] += 1.0
        elif var=='thickness':
            p = 0
        elif var=='radii':
            p = 0
        elif var=='pos':
            p = 0
        elif var=='rot':
            j = ind[0]
            k = ind[1]
            for f in range(len(self.Qs)):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(3*nj+nj*k+j).todense()*numpy.pi/180.0)
        elif var=='shapeR0':
            p = 0
        elif var=='shapeT0':
            p = 0
        elif var=='shapeL0':
            p = 0
        elif var=='shapeB0':
            p = 0
        elif var=='shapeR1':
            p = 0
        elif var=='shapeT1':
            p = 0
        elif var=='shapeL1':
            p = 0
        elif var=='shapeB1':
            p = 0


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
