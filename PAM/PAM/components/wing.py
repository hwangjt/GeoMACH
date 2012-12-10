from __future__ import division
from PAM.components import Component, Property, airfoils
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib


class Wing(Component):
    """ A component used to model lifting surfaces. """

    def __init__(self, nx=1, nz=1, left=2, right=2):
        """ Initialization method
        nx: integer
            Number of surfaces in x (chord-wise) direction
        nz: integer
            Number of surfaces in z (span-wise) direction
        left, right: integer
            The v[0] and v[-1] sections of the wing
            0: open tip, C0
            1: open tip, C1
            2: closed tip
        """ 

        super(Wing,self).__init__() 

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(None)
        self.ms.append(numpy.zeros(nz,int))

        self.addFace(-1, 3, 0.5)
        self.addFace( 1, 3,-0.5)
        self.connectEdges(f1=0,u1=0,f2=1,u2=-1)
        self.connectEdges(f1=0,u1=-1,f2=1,u2=0)
        if left==2:
            self.connectEdges(f1=0,v1=-1,f2=1,v2=-1)
        if right==2:
            self.connectEdges(f1=0,v1=0,f2=1,v2=0)

        self.left = left
        self.right = right
        self.ax1 = 3
        self.ax2 = 2

    def setDOFs(self):
        setC1 = self.setC1
        setCornerC1 = self.setCornerC1
        for f in range(2):
            setC1('surf', f, val=True) #C1 Everywhere
            setC1('surf', f, i=-f, u=-f, val=False) #C0 trailing edge
            setC1('edge', f, i=-f, u=-f, val=True) #C0 trailing edge
            if self.left==0:  
                setC1('surf', f, j=-1, v=-1, val=False) #C0 left edge
                setC1('edge', f, j=-1, v=-1, val=True) #C0 left edge
                setCornerC1(f, i=-f, j=-1, val=False) #C0 left TE corner     
            if self.right==0:         
                setC1('surf', f, j=0, v=0, val=False) #C0 right edge
                setC1('edge', f, j=0, v=0, val=True) #C0 right edge
                setCornerC1(f, i=-f, j=0, val=False) #C0 right TE corner

    def initializeVariables(self):
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        zeros = numpy.zeros
        ones = numpy.ones
        self.variables = {
            'offset':zeros(3),
            'chord':ones(nj),
            'pos':zeros((nj,3),order='F'),
            'rot':zeros((nj,3),order='F'),
            'shape':zeros((2,ni,nj,3),order='F')
            }
        self.parameters = {
            'nor':ones((nj,3),order='F')
            }
        self.setAirfoil()

    def setAirfoil(self,filename="naca0012"):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.variables['shape'][f,:,j,:2] = Ps[f][:,:]
        
    def computeQs(self):
        r = numpy.array([0.25,0,0])
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        v = self.variables
        p = self.parameters
        ax1 = self.ax1
        ax2 = self.ax2

        #if self.left==2:
        #    v['pos'][-1] = 2*v['pos'][-2] - v['pos'][-3]
        #if self.right==2:
        #    v['pos'][0] = 2*v['pos'][1] - v['pos'][2]
        rot0, Da, Di, Dj = PAMlib.computerotations(ax1, ax2, nj, 9*(nj*3-2), v['pos'], p['nor'])
        self.drot0_dpos = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(nj*3,nj*3))
        rot = v['rot']*numpy.pi/180.0 + rot0

        self.dQs_dv = range(2)
        for f in range(2):
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computesections(ax1, ax2, f, ni, nj, ni*nj*21, f*3*ni*nj, r, v['offset'], v['chord'], v['pos'], rot, v['shape'][f,:,:,:])
            self.dQs_dv[f] = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nj*(4+6*ni)))

    def setDerivatives(self, var, ind):
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        if var=='offset':
            for f in range(2):
                self.Qs[f][:,:,ind] += 1.0
        elif var=='chord':
            for f in range(2):
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(ind).todense())
        elif var=='pos':
            j = ind[0]
            k = ind[1]
            A = scipy.sparse.csc_matrix((nj,3*nj))
            B = self.drot0_dpos
            C = scipy.sparse.csc_matrix((self.dQs_dv[0].shape[1]-4*nj,3*nj))
            D = scipy.sparse.vstack([A,B,C],format='csc')
            for f in range(2):
                self.Qs[f][:,j,k] += 1.0
                Q = self.dQs_dv[f].dot(D).todense()[:,nj*k+j]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, Q)
        elif var=='rot':
            j = ind[0]
            k = ind[1]
            for f in range(2):
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(nj+nj*k+j).todense()*numpy.pi/180.0)
        elif var=='shape':
            f = ind[0]
            i = ind[1]
            j = ind[2]
            k = ind[3]
            self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(4*nj+f*3*ni*nj+ni*nj*k+ni*j+i).todense())


if __name__ == '__main__':
    w = Wing(nx=2,nz=2,left=0)
    import PUBS
    from mayavi import mlab
    w.oml0 = PUBS.PUBS(w.Ps)
    w.setDOFs()
    w.oml0.updateBsplines()
    w.computems()
    w.initializeDOFmappings()
    w.initializeVariables()
    w.variables['pos'][:,2] = numpy.linspace(0,1,w.Qs[0].shape[1])
    #for j in range(w.Qs[0].shape[1]):
    #    w.variables['shape'][0,:,j,0] = 1 - numpy.linspace(0,1,w.Qs[0].shape[0])
    #    w.variables['shape'][1,:,j,0] = numpy.linspace(0,1,w.Qs[0].shape[0])
    w.variables['pos'][:,0] = numpy.linspace(0,1,w.Qs[0].shape[1])
    w.variables['pos'][:,1] = numpy.linspace(0,1,w.Qs[0].shape[1])
    #w.variables['rot'][:,2] = 20
    w.parameters['nor'][:,:] = 1.0
    w.setAirfoil("naca0012.dat")
    w.computeQs()
    w.propagateQs()
    w.oml0.computePoints()
    w.oml0.plot()
    name='wing'
    w.oml0.write2Tec(name)
