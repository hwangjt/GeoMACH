from __future__ import division
from PAM.components import Component, Property, airfoils
import numpy, pylab, time
import mpl_toolkits.mplot3d.axes3d as p3
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

        self.addFace(-1, 3, 1, 1, 1, left==2, right==2)
        self.addFace( 1, 3,-1, 1, 1, left==2, right==2)

        self.left = left
        self.right = right

    def setDOFs(self):
        left = self.left
        right = self.right

        for f in range(2):
            self.setC1('surf', f, val=True) #C1 Everywhere
            self.setC1('surf', f, i=-f, u=-f, val=False) #C0 trailing edge
            self.setC1('edge', f, i=-f, u=-f, val=True) #C0 trailing edge
            if left==0:                
                self.setC1('surf', f, j=0, v=0, val=False) #C0 left edge
                self.setC1('edge', f, j=0, v=0, val=True) #C0 left edge
                self.setCornerC1(f, i=-f, j=0, val=False) #C0 left TE corner
            if right==0:
                self.setC1('surf', f, j=-1, v=-1, val=False) #C0 right edge
                self.setC1('edge', f, j=-1, v=-1, val=True) #C0 right edge
                self.setCornerC1(f, i=-f, j=-1, val=False) #C0 right TE corner

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
        
    def propagateQs(self):
        r = numpy.array([0.25,0,0])
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        v = self.variables
        for f in range(2):
            self.Qs[f][:,:,:] = PAMlib.computewingsections(ni, nj, r, v['offset'], v['chord'], v['pos'], v['rot'], v['nor'], v['shape'][f,:,:,:])

    def setAirfoil(self,filename):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.variables['shape'][f,:,j,:2] = Ps[f][:,:]


if __name__ == '__main__':
    w = Wing(nx=2,nz=2)#,left=0)
    import PUBS
    from mayavi import mlab
    w.oml0 = PUBS.PUBS(w.Ps)
    w.setDOFs()
    w.oml0.updateBsplines()
    w.computems()
    w.initializeDOFmappings()
    w.initializeVariables()
    #for i in range(w.Qs[0].shape[0]):
    #    w.Qs[0][i,:,2] = numpy.linspace(0,1,w.Qs[0].shape[1])
    #    w.Qs[1][i,:,2] = numpy.linspace(0,1,w.Qs[0].shape[1])
    #for i in range(w.Qs[0].shape[1]):
    #    w.Qs[0][::-1,i,0] = numpy.linspace(0,1,w.Qs[0].shape[0])
    #    w.Qs[1][:,i,0] = numpy.linspace(0,1,w.Qs[0].shape[0])
    #w.Qs[0][:,:,1] = 0.1
    #w.Qs[1][:,:,1] = -0.1
    #w.Qs[0][0,:,1] = 0.0
    #w.Qs[1][-1,:,1] = 0.0
    w.variables['pos'][:,2] = numpy.linspace(0,1,w.Qs[0].shape[1])
    for j in range(w.Qs[0].shape[1]):
        w.variables['shape'][0,:,j,0] = 1 - numpy.linspace(0,1,w.Qs[0].shape[0])
        w.variables['shape'][1,:,j,0] = numpy.linspace(0,1,w.Qs[0].shape[0])
    #w.variables['pos'][:,0] = numpy.linspace(0,0.3,w.Qs[0].shape[1])
    w.setAirfoil("naca0012.dat")
    w.propagateQs()
    w.updateQs()
    w.oml0.computePoints()
    w.oml0.plot(pylab.figure(),False)
    pylab.show()
