from __future__ import division
from PAM.components import Component, Property
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib


class Shell(Component):
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

        super(Shell,self).__init__() 

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(numpy.zeros(ny,int))
        self.ms.append(numpy.zeros(nz,int))

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
        nx = self.Qs[2].shape[1]
        ny = self.Qs[2].shape[0]
        nz = self.Qs[3].shape[0]
        zeros = numpy.zeros
        ones = numpy.ones
        self.variables = {
            'offset':zeros(3),
            'thickness':0.1*ones((nx,3),order='F'),
            'radii':ones((nx,3),order='F'),
            'pos':zeros((nx,3),order='F'),
            'rot':zeros((nx,3),order='F'),
            'shapeR0':zeros((ny,nx),order='F'),
            'shapeT0':zeros((nz,nx),order='F'),
            'shapeL0':zeros((ny,nx),order='F'),
            'shapeB0':zeros((nz,nx),order='F'),
            'shapeR1':zeros((ny,nx),order='F'),
            'shapeT1':zeros((nz,nx),order='F'),
            'shapeL1':zeros((ny,nx),order='F'),
            'shapeB1':zeros((nz,nx),order='F')
            }
        self.parameters = {
            'nor':ones((nx,3),order='F'),
            'fillet':zeros((nx,4))
            }
        self.setSections()

    def setSections(self, sections=[], t1U=0, t2U=0, t1L=0, t2L=0):
        Ns = self.Ns
        v = self.variables
        p = self.parameters
        for j in range(Ns[0].shape[1]):
            for i in range(Ns[0].shape[0]):
                val = Ns[0][i,j,3]
                if not val == -1:
                    break
            found = False
            for k in range(len(sections)):
                found = found or (val==sections[k])
            if found or sections==[]:
                p['fillet'][j,0] = t1U
                p['fillet'][j,1] = t2U
                p['fillet'][j,2] = t1L
                p['fillet'][j,3] = t2L

    def computeQs(self):
        r = numpy.zeros(3)
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        v = self.variables
        p = self.parameters
        b = self.bottom==2
        ax1 = self.ax1
        ax2 = self.ax2

        rot0, Da, Di, Dj = PAMlib.computerotations(ax1, ax2, nx, 9*(nx*3-2), v['pos'], p['nor'])
        drot0_dpos = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(nx*3,nx*3))
        rot = v['rot']*numpy.pi/180.0 + rot0
        r0 = v['radii'] + v['thickness']
        r1 = v['radii'] - v['thickness']

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
        self.dQs_dv = range(len(self.Qs))

        counter = 0
        for f in range(len(self.Qs)):
            if f in [0,1,2,6]:
                radii = r0
            else:
                radii = r1
            ni, nj = self.Qs[f].shape[:2]
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computesections(ax1, ax2, ni, nj, ni*nj*27, counter, r, v['offset'], radii, v['pos'], rot, shapes[f])
            self.dQs_dv[f] = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nQ))
            counter += 3*ni*nj

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