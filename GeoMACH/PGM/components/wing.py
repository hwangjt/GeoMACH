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

        self.addFace('upp', -1, 3, 0.5)
        self.addFace('low', 1, 3, -0.5)

        self.left = left
        self.right = right
        self.ax1 = 3
        self.ax2 = 2

    def setDOFs(self):
        for f in xrange(len(self.faces)):
            face = self.faces.values()[f]
            face.setC1('surf', val=True) #C1 Everywhere
            face.setC1('surf', i=-f, u=-f, val=False) #C0 trailing edge
            face.setC1('edge', i=-f, u=-f, val=True) #C0 trailing edge
            if self.left==0:  
                face.setC1('surf', j=-1, v=-1, val=False) #C0 left edge
                face.setC1('edge', j=-1, v=-1, val=True) #C0 left edge
                face.setCornerC1(i=-f, j=-1, val=False) #C0 left TE corner
            if self.right==0:
                face.setC1('surf', j=0, v=0, val=False) #C0 right edge
                face.setC1('edge', j=0, v=0, val=True) #C0 right edge
                face.setCornerC1(i=-f, j=0, val=False) #C0 right TE corner

    def declare_properties(self):
        super(Wing, self).declare_properties()

        ni = self.faces['upp'].num_cp[0]
        nj = self.faces['upp'].num_cp[1]
        self.airfoils = {'upp': numpy.zeros((ni,nj,3),order='F'),
                         'low': numpy.zeros((ni,nj,3),order='F'),
                         }
        if self.oml0 is not None:
            self.setAirfoil()
        else:
            self.airfoils['upp'][:,:,:] = 0.0
            self.airfoils['low'][:,:,:] = 0.0
            self.airfoils['upp'][1:-1,:,1] = 0.05
            self.airfoils['low'][1:-1,:,1] = -0.05
            n = self.airfoils['upp'].shape[0]
            for i in range(n):
                self.airfoils['upp'][i,:,0] = 1 - i/(n-1)
                self.airfoils['low'][i,:,0] = i/(n-1)

    def setAirfoil(self,filename="naca0012"):
        Ps = airfoils.fitAirfoil(self, filename)
        for name in self.airfoils:
            for j in range(self.faces[name].num_cp[1]):
                self.airfoils[name][:,j,:2] = Ps[name][:,:]
        
    def computeQs(self):
        faces = self.faces
        ni = faces['upp'].num_cp[0]
        nj = faces['upp'].num_cp[1]
        p = self.properties

        #if self.left==2:
        #    v['pos'][-1] = 2*v['pos'][-2] - v['pos'][-3]
        #if self.right==2:
        #    v['pos'][0] = 2*v['pos'][1] - v['pos'][2]

        shapes = self.shapes
        for name in shapes:
            shapes[name][:,:,:] = self.airfoils[name][:,:,:]
            shapes[name][:,:,1] += p['shp', name]
        self.computeSections()
