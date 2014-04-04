from __future__ import division
import numpy

from GeoMACH.PGM.components import Primitive, airfoils


class Wing(Primitive):
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

        super(Wing,self).__init__(nx,0,nz)

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
        faces = self.faces
        for f in range(2):
            faces[f].setC1('surf', val=True) #C1 Everywhere
            faces[f].setC1('surf', i=-f, u=-f, val=False) #C0 trailing edge
            faces[f].setC1('edge', i=-f, u=-f, val=True) #C0 trailing edge
            if self.left==0:  
                faces[f].setC1('surf', j=-1, v=-1, val=False) #C0 left edge
                faces[f].setC1('edge', j=-1, v=-1, val=True) #C0 left edge
                faces[f].setCornerC1(i=-f, j=-1, val=False) #C0 left TE corner
            if self.right==0:
                faces[f].setC1('surf', j=0, v=0, val=False) #C0 right edge
                faces[f].setC1('edge', j=0, v=0, val=True) #C0 right edge
                faces[f].setCornerC1(i=-f, j=0, val=False) #C0 right TE corner

    def initializeVariables(self):
        super(Wing,self).initializeVariables()
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        zeros = numpy.zeros
        v = self.variables
        a = self.addParam

        v['shU'] = zeros((ni,nj),order='F')
        v['shL'] = zeros((ni,nj),order='F')

        a('shU','shU',(1,1),P=[0.0])
        a('shL','shL',(1,1),P=[0.0])
        self.params['pos'].setP([[0.,0.,0.],[0.,0.,1.]])
        self.params['ogn'].setP([0.25,0.,0.])

        self.shapeU = zeros((ni,nj,3),order='F')
        self.shapeL = zeros((ni,nj,3),order='F')
        self.setAirfoil()

    def setAirfoil(self,filename="naca0012"):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.faces)):
            for j in range(self.Ns[f].shape[1]):
                shape = self.shapeU if f==0 else self.shapeL
                shape[:,j,:2] = Ps[f][:,:]
        
    def computeQs(self):
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        v = self.variables

        #if self.left==2:
        #    v['pos'][-1] = 2*v['pos'][-2] - v['pos'][-3]
        #if self.right==2:
        #    v['pos'][0] = 2*v['pos'][1] - v['pos'][2]

        shapes = [self.shapeU, self.shapeL]
        shapes[0][:,:,1] += v['shU']
        shapes[1][:,:,1] -= v['shL']
        nQ = nj*(9+6*ni)
        self.computeSections(nQ, shapes)
        shapes[0][:,:,1] -= v['shU']
        shapes[1][:,:,1] += v['shL']


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
    w.oml0.write2TecC(name)
