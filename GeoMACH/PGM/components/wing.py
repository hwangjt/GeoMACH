from __future__ import division
import numpy

from GeoMACH.PGM.components import Primitive, airfoils, Property


class Wing(Primitive):
    """ A component used to model lifting surfaces. """

    def __init__(self, nx=1, nz=1, left=1, right=1):
        """ Initialization method
        nx: integer
            Number of surfaces in x (chord-wise) direction
        nz: integer
            Number of surfaces in z (span-wise) direction
        left, right: integer
            The v[0] and v[-1] sections of the wing
            0: open tip, C0
            1: open tip, C1
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

        if self.oml0 is not None:
            self.setAirfoil()
        else:
            self.shapes['upp'][:,:,:] = 0.0
            self.shapes['low'][:,:,:] = 0.0
            self.shapes['upp'][1:-1,:,1] = 0.05
            self.shapes['low'][1:-1,:,1] = -0.05
            n = self.shapes['upp'].shape[0]
            for i in range(n):
                self.shapes['upp'][i,:,0] = 1 - i/(n-1)
                self.shapes['low'][i,:,0] = i/(n-1)

    def setAirfoil(self,filename="naca0012"):
        Ps = airfoils.fitAirfoil(self, filename)
        for name in self.shapes:
            for j in range(self.faces[name].num_cp[1]):
                self.shapes[name][:,j,:2] = Ps[name][:,:]
        
    def computeQs(self):
        return self.computeSections()

    def add_thickness_constraints_grid(self, nu, nv, u1, u2, v1, v2):
        locations = numpy.zeros((nu,nv,2))
        for i in range(nu):
            for j in range(nv):
                locations[i,j,0] = u1 + (u2-u1) * i/(nu-1)
                locations[i,j,1] = v1 + (v2-v1) * j/(nv-1)
        self.add_thickness_constraints(locations)

    def add_thickness_constraints(self, locations):
        ni, nj = locations.shape[:2]

        face = self.faces['upp']
        increments = [[None, None], [None, None]]
        for f in xrange(2):
            for d in xrange(2):
                n = face.num_surf[d]
                increments[f][d] = numpy.zeros(n + 1)
                for i in xrange(n+1):
                    increments[f][d][i] = sum(face.num_cp_list[d][:i]) / sum(face.num_cp_list[d])
            if f==0:
                increments[f][0][:] = 1 - increments[f][0][::-1]

        surf = numpy.zeros((2,ni,nj), dtype=int, order='F')
        locs = [numpy.zeros((2,ni,nj), order='F'),
                numpy.zeros((2,ni,nj), order='F')]
        for f in xrange(2):
            for i in range(ni):
                for j in range(nj):
                    if f==0:
                        loc_face = [1-locations[i,j,0], locations[i,j,1]]
                    else:
                        loc_face = [locations[i,j,0], locations[i,j,1]]
                    loc_surf = [-1,-1]
                    for d in xrange(2):
                        n = face.num_surf[d]
                        for k in range(n):
                            if increments[f][d][k] <= loc_face[d] and loc_face[d] <= increments[f][d][k+1]:
                                locs[d][f,i,j] = (loc_face[d] - increments[f][d][k]) / \
                                    (increments[f][d][k+1] - increments[f][d][k])
                                loc_surf[d] = k
                    if loc_surf[0] == -1 or loc_surf[1] == -1:
                        raise Exception('Invalid thickness constraint locations')
                    surf[f,i,j] = self.faces.values()[f].surf_indices[loc_surf[0], loc_surf[1]]

        J = self.oml0.evaluateBases(surf.flatten(order='F'), locs[0].flatten(order='F'), locs[1].flatten(order='F'))
        self.oml0.export.write2TecScatter('thk.dat', J.dot(self.oml0.C[:,:3]), ['x','y','z'])
