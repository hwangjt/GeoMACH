from __future__ import division
from PAM.components import Component, Property
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib
import mpl_toolkits.mplot3d.axes3d as p3


class Body(Component):
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

        super(Body,self).__init__() 

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(numpy.zeros(ny,int))
        self.ms.append(numpy.zeros(nz,int))

        self.addFace(-2, 3,-0.5)
        self.addFace(-2,-3, 0.5)
        self.addFace( 2, 1,-0.5)
        self.addFace( 3, 1, 0.5)
        self.addFace(-2, 1, 0.5)
        if bottom==2:
            self.addFace(-3, 1,-0.5)

        self.bottom = bottom

    def setDOFs(self):
        setC1 = self.setC1
        setCornerC1 = self.setCornerC1
        for f in range(len(self.Ks)):
            self.setC1('surf', f, val=True)
        if self.bottom==0:
            setC1('surf', 0, i=-1, u=-1, val=False)
            setC1('edge', 0, i=-1, u=-1, val=True)
            setC1('surf', 1, i=-1, u=-1, val=False)
            setC1('edge', 1, i=-1, u=-1, val=True)
            setC1('surf', 2, i= 0, u= 0, val=False)
            setC1('edge', 2, i= 0, u= 0, val=True)
            setC1('surf', 4, i=-1, u=-1, val=False)
            setC1('edge', 4, i=-1, u=-1, val=True)

    def initializeVariables(self):
        nx = self.Qs[2].shape[1]
        ny = self.Qs[2].shape[0]
        nz = self.Qs[3].shape[0]
        zeros = numpy.zeros
        ones = numpy.ones
        self.variables = {
            'noseL':0.1,
            'tailL':0.1,
            'offset':zeros(3),
            'radii':ones((nx,3),order='F'),
            'pos':zeros((nx,3),order='F'),
            'rot':zeros((nx,3),order='F'),
            'nor':ones((nx,3),order='F'),
            'fillet':zeros((nx,4)),
            'shapeR':zeros((ny,nx),order='F'),
            'shapeT':zeros((nz,nx),order='F'),
            'shapeL':zeros((ny,nx),order='F'),
            'shapeB':zeros((nz,nx),order='F')
            }
        self.setSections()

    def setSections(self, sections=[], t1U=0, t2U=1, t1L=0, t2L=1):
        Ns = self.Ns
        v = self.variables
        for j in range(Ns[2].shape[1]):
            for i in range(Ns[2].shape[0]):
                val = Ns[2][i,j,3]
                if not val == -1:
                    break
            found = False
            for k in range(len(sections)):
                found = found or (val==sections[k])
            if found or sections==[]:
                v['fillet'][j,0] = t1U
                v['fillet'][j,1] = t2U
                v['fillet'][j,2] = t1L
                v['fillet'][j,3] = t2L

    def propagateQs(self):
        r = numpy.zeros(3)
        nx = self.Qs[2].shape[1]
        ny = self.Qs[2].shape[0]
        nz = self.Qs[3].shape[0]
        v = self.variables
        b = self.bottom==2

        Ns = self.Ns
        Qs = self.Qs

        rot0, Da, Di, Dj = PAMlib.computerotations(nx, 9*(nx*3-2), v['pos'])
        drot0_dpos = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(nx*3,nx*3))
        rot = v['rot']*numpy.pi/180.0 + rot0*v['nor']
        shapeR = PAMlib.computeshape(ny, nx, (-b)/4.0, 1/4.0, v['radii'], v['fillet'], v['shapeR'])
        shapeT = PAMlib.computeshape(nz, nx, 1/4.0, 3/4.0, v['radii'], v['fillet'], v['shapeT'])
        shapeL = PAMlib.computeshape(ny, nx, 3/4.0, (4+b)/4.0, v['radii'], v['fillet'], v['shapeL'])
        shapeB = PAMlib.computeshape(nz, nx, 5/4.0, 7/4.0, v['radii'], v['fillet'], v['shapeB'])
        chord = numpy.ones(nx)

        self.dQs_dv = range(len(self.Qs))

        self.Qs[2][:,:,:], Da, Di, Dj = PAMlib.computesections(-1, ny, nx, nx*ny*24, 0, r, v['offset'], chord, v['pos'], rot, shapeR)
        self.dQs_dv[2] = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(3*nx*ny,nx*(9+6*ny+6*nz)))

        self.Qs[3][:,:,:], Da, Di, Dj = PAMlib.computesections(-1, nz, nx, nx*nz*24, 3*nx*ny, r, v['offset'], chord, v['pos'], rot, shapeT)
        self.dQs_dv[3] = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(3*nx*nz,nx*(9+6*ny+6*nz)))

        self.Qs[4][:,:,:], Da, Di, Dj = PAMlib.computesections(-1, ny, nx, nx*ny*24, 3*nx*(ny+nz), r, v['offset'], chord, v['pos'], rot, shapeL)
        self.dQs_dv[4] = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(3*nx*ny,nx*(9+6*ny+6*nz)))

        if b:
            self.Qs[5][:,:,:], Da, Di, Dj = PAMlib.computesections(-1, nz, nx, nx*nz*24, 3*nx*(2*ny+nz), r, v['offset'], chord, v['pos'], rot, shapeB)
            self.dQs_dv[5] = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(3*nx*nz,nx*(9+6*ny+6*nz)))

        #Qs[f][:,:,:] = PAMlib.computecone(self.bottom, ny, nz, L, p['posy'].data[i0], p['posy'].data[i1], p['posy'].data[i2], p['ry'].data[i1], p['ry'].data[i2], p['rz'].data[i1], p['rz'].data[i2], dx)
        

        for f in []:#range(len(Ns)):
            if f==0 or f==4:
                if f==0:
                    L = self.props['noseL']
                    i0 = 0
                    i1 = 1
                    i2 = 2
                else:
                    L = self.props['tailL']
                    i0 = -1
                    i1 = -2
                    i2 = -3
                p = self.props
                dx = abs(self.props['posx'].data[i2] - self.props['posx'].data[i1])  
                Qs[f][:,:,:] = PAMlib.computecone(self.full, Ns[f].shape[0], Ns[f].shape[1], L, p['posy'].data[i0], p['posy'].data[i1], p['posy'].data[i2], p['ry'].data[i1], p['ry'].data[i2], p['rz'].data[i1], p['rz'].data[i2], dx)
                if f==4:
                    Qs[f][:,:,0] *= -1
                    Qs[f][:,:,0] += self.props['noseL']
                    Qs[f][:,:,0] += self.props['posx'].data[-2] - self.props['posx'].data[1]
                    Qs[f][:,:,0] += self.props['tailL']
                    Qs[f][:,:,:] = Qs[f][:,::-1,:]
                Qs[f][:,:,0] += self.offset[0]
                Qs[f][:,:,1] += self.offset[1]
                Qs[f][:,:,2] += self.offset[2]

if __name__ == '__main__':
    #b = Body(nx=2,ny=2,nz=2,bottom=0)
    b = Body(nx=4,ny=2,nz=2,bottom=2)
    import PUBS
    from mayavi import mlab
    b.oml0 = PUBS.PUBS(b.Ps)
    b.setDOFs()
    b.oml0.updateBsplines()
    b.computems()
    b.initializeDOFmappings()
    b.initializeVariables()
    b.variables['pos'][:,0] = numpy.linspace(0,1,b.Qs[2].shape[1])
    #b.variables['shapeL'][:,:] = 0.5
    #b.variables['shapeR'][:,:] = 0.5
    b.propagateQs()
    b.updateQs()
    b.oml0.computePoints()
    b.oml0.plot(pylab.figure(),False)
    export = PUBS.PUBSexport(b.oml0)
    name='body'
    export.write2Tec(name)
    export.write2TecC(name+'_C')
    pylab.show()
