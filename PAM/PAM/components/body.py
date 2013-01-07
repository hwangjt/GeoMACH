from __future__ import division
from PAM.components import Primitive, Property
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib


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
        setCornerC1 = self.setCornerC1
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
        self.variables['shapeR'] = numpy.zeros((ny,nx),order='F')
        self.variables['shapeT'] = numpy.zeros((nz,nx),order='F')
        self.variables['shapeL'] = numpy.zeros((ny,nx),order='F')
        self.variables['shapeB'] = numpy.zeros((nz,nx),order='F')
        self.parameters['fillet'] = numpy.zeros((nx,4),order='F')
        self.setSections()

    def computeQs(self):
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        v = self.variables
        p = self.parameters
        b = self.bottom==2

        #v['pos'][0] = 2*v['pos'][1] - v['pos'][2]
        #v['pos'][-1] = 2*v['pos'][-2] - v['pos'][-3]

        rot, self.drot0_dpos = self.computeRotations()

        shapes = range(4)
        shapes[0] = PAMlib.computeshape(ny, nx,-b/4.0, 1/4.0, p['fillet'], v['shapeR'])
        shapes[1] = PAMlib.computeshape(nz, nx, 1/4.0, 3/4.0, p['fillet'], v['shapeT'])
        shapes[2] = PAMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, p['fillet'], v['shapeL'])
        shapes[3] = PAMlib.computeshape(nz, nx, 5/4.0, 7/4.0, p['fillet'], v['shapeB'])

        nQ = nx*(6+6*ny+6*nz) if self.bottom==2 else nx*(6+6*ny+3*nz)
        self.computeSections(nQ, rot, shapes)

    def setDerivatives(self, var, ind):
        nx = self.Qs[0].shape[1]
        ny = self.Qs[0].shape[0]
        nz = self.Qs[1].shape[0]
        v = self.variables
        if var=='offset':
            for f in range(len(self.Qs)):
                self.Qs[f][:,:,ind] += 1.0
        elif var=='radii':
            j = ind[0]
            k = ind[1]
        elif var=='pos':
            j = ind[0]
            k = ind[1]
            A = scipy.sparse.csc_matrix((3*nx,3*nx))
            B = self.drot0_dpos
            C = scipy.sparse.csc_matrix((self.dQs_dv[2].shape[1]-6*nx,3*nx))
            D = scipy.sparse.vstack([A,B,C],format='csc')
            for f in range(2,len(self.Qs)):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,j,k] += 1.0
                Q = self.dQs_dv[f].dot(D).getcol(nj*k+j).todense()
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, Q)
            if j==1:
                ni, nj = self.Qs[0].shape[:2]
                self.Qs[0][:,:,:] += 1.0
                for l in range(3):
                    self.Qs[0][:,:,:] += self.dQ_drot[0][:,:,:,l]*self.drot0_dpos[nj*l+1,nj*k+j]
            elif j==nx-2:
                ni, nj = self.Qs[0].shape[:2]
                self.Qs[1][:,:,:] += 1.0
                for l in range(3):
                    self.Qs[1][:,:,:] += self.dQ_drot[1][:,:,:,l]*self.drot0_dpos[nj*l+1,nj*k+j]
        elif var=='rot':
            j = ind[0]
            k = ind[1]
            for f in range(2,len(self.Qs)):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(3*nj+nj*k+j).todense()*numpy.pi/180.0)
            #if j==1:
            #    self.Qs[0][:,:,:] += self.dQ_drot[0][:,:,:,k]*numpy.pi/180.0
            #elif j==nx-2:
            #    self.Qs[1][:,:,:] += self.dQ_drot[1][:,:,:,k]*numpy.pi/180.0
        elif var=='shapeR':
            p = 0
        elif var=='shapeT':
            p = 0
        elif var=='shapeL':
            p = 0
        elif var=='shapeB':
            p = 0


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
