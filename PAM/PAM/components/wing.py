from __future__ import division
from PAM.components import Component, Property, airfoils
import numpy, pylab, time, scipy.sparse
import mpl_toolkits.mplot3d.axes3d as p3
import PAM.PAMlib as PAMlib

import PUBS


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
            self.connectEdges(f1=0,v1=0,f2=1,v2=0)
        if right==2:
            self.connectEdges(f1=0,v1=-1,f2=1,v2=-1)

        self.left = left
        self.right = right

    def setDOFs(self):
        setC1 = self.setC1
        setCornerC1 = self.setCornerC1
        for f in range(2):
            setC1('surf', f, val=True) #C1 Everywhere
            setC1('surf', f, i=-f, u=-f, val=False) #C0 trailing edge
            setC1('edge', f, i=-f, u=-f, val=True) #C0 trailing edge
            if self.left==0:                
                setC1('surf', f, j=0, v=0, val=False) #C0 left edge
                setC1('edge', f, j=0, v=0, val=True) #C0 left edge
                setCornerC1(f, i=-f, j=0, val=False) #C0 left TE corner
            if self.right==0:
                setC1('surf', f, j=-1, v=-1, val=False) #C0 right edge
                setC1('edge', f, j=-1, v=-1, val=True) #C0 right edge
                setCornerC1(f, i=-f, j=-1, val=False) #C0 right TE corner

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
            'nor':zeros((nj,3),order='F'),
            'shape':zeros((2,ni,nj,3),order='F')
            }
        self.setAirfoil()

    def setAirfoil(self,filename="naca0012"):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.variables['shape'][f,:,j,:2] = Ps[f][:,:]
        
    def propagateQs(self):
        r = numpy.array([0.25,0,0])
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        v = self.variables

        rot0, Da, Di, Dj = PAMlib.computerotations(nj, 9*(nj*3-2), v['pos'])
        drot0_dpos = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(nj*3,nj*3))
        rot = v['rot']*numpy.pi/180.0 + rot0*v['nor']

        self.dQs_dv = range(2)
        for f in range(2):
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computesections(f, ni, nj, ni*nj*24, f*3*ni*nj, r, v['offset'], v['chord'], v['pos'], rot, v['shape'][f,:,:,:])
            self.dQs_dv[f] = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nj*(7+6*ni)))


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
    for j in range(w.Qs[0].shape[1]):
        w.variables['shape'][0,:,j,0] = 1 - numpy.linspace(0,1,w.Qs[0].shape[0])
        w.variables['shape'][1,:,j,0] = numpy.linspace(0,1,w.Qs[0].shape[0])
    w.variables['pos'][:,0] = numpy.linspace(0,1,w.Qs[0].shape[1])
    w.variables['pos'][:,1] = numpy.linspace(0,1,w.Qs[0].shape[1])
    #w.variables['rot'][:,2] = 20
    w.variables['nor'][:,:] = 1.0
    w.setAirfoil("naca0012.dat")
    w.propagateQs()
    w.updateQs()
    w.oml0.computePoints()
    w.oml0.plot(pylab.figure(),False)
    pylab.show()
