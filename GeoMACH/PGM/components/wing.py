from __future__ import division
import numpy

from GeoMACH.PGM.components import Primitive, airfoils, Property
import scipy.sparse


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

    def add_thk_con(self, name, nu, nv, urange, vrange, factor):
        self.funcs[name] = WingThicknessFunction(self, nu, nv, urange, vrange, factor)



class WingFunction(object):

    def __init__(self, comp):
        self.oml = comp.oml0
        self.comp = comp

    def get_grid(self, nu, nv, urange, vrange):
        locations = numpy.zeros((nu,nv,2))
        for i in range(nu):
            for j in range(nv):
                locations[i,j,0] = urange[0] + (urange[1]-urange[0]) * i/(nu-1)
                locations[i,j,1] = vrange[0] + (vrange[1]-vrange[1]) * j/(nv-1)
        ni, nj = nu, nv

        face = self.comp.faces['upp']
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
                    surf[f,i,j] = self.comp.faces.values()[f].surf_indices[loc_surf[0], loc_surf[1]]

        J = self.oml.evaluateBases(surf.flatten(order='F'), locs[0].flatten(order='F'), locs[1].flatten(order='F'))
        J = scipy.sparse.block_diag((J, J, J), format='csc')
        return J


class WingThicknessFunction(WingFunction):

    def __init__(self, comp, nu, nv, urange, vrange, factor):
        super(WingThicknessFunction, self).__init__(comp)
        self.J = self.get_grid(nu, nv, urange, vrange)
        self.M = scipy.sparse.block_diag((self.oml.M, self.oml.M, self.oml.M), format='csc')
        self.factor = factor
        self.nu, self.nv = nu, nv

        self.pts = numpy.zeros(2*nu*nv*3,)
        self.pts_array = self.pts.reshape((2,nu,nv,3), order='F')

        self.func = numpy.zeros(nu*nv)
        self.func_array = self.func.reshape((nu,nv), order='F')

        self.func0 = numpy.array(self.func)
        
#        self.oml.export.write2TecScatter('thk.dat', self.pts.reshape((2*nu*nv,3),order='F'), ['x','y','z'])

    def initialize(self):
        self.compute()
        self.func0[:] = self.func[:]

    def compute(self):
        self.pts[:] = self.J.dot(self.oml.C[:,:3].reshape(3*self.oml.nC, order='F'))

        self.func_array[:,:] = 0.0
        for k in xrange(3):
            self.func_array[:,:] += (self.pts_array[0,:,:,k] - self.pts_array[1,:,:,k])**2
        self.func_array[:,:] = self.func_array**0.5

    def get_func(self):
        self.compute()
        return self.factor * self.func0[:] - self.func[:]

    def get_jacobian(self):
        self.compute()
        nu, nv = self.nu, self.nv

        nD = 6 * nu * nv
        Da = numpy.zeros((2,nu,nv,3), order='F')
        Di = numpy.zeros((2,nu,nv,3), dtype=int, order='F')
        Dj = numpy.zeros((2,nu,nv,3), dtype=int, order='F')
        
        pts_indices = numpy.array(numpy.linspace(0, 2*nu*nv*3-1, 2*nu*nv*3), int).reshape((2,nu,nv,3), order='F')

        for f in xrange(2):
            for k in xrange(3):
                Da[f,:,:,k] = -(self.pts_array[f,:,:,k] - self.pts_array[1-f,:,:,k]) / self.func_array[:,:]
                Di[f,:,:,k] = numpy.linspace(0, nu*nv-1, nu*nv).reshape((nu,nv), order='F')
                Dj[f,:,:,k] = pts_indices[f,:,:,k]

        Da = Da.reshape(2*nu*nv*3, order='F')
        Di = Di.reshape(2*nu*nv*3, order='F')
        Dj = Dj.reshape(2*nu*nv*3, order='F')
        D = scipy.sparse.csr_matrix((Da, (Di, Dj)), 
                                    shape=(nu*nv, 2*nu*nv*3))
        return D.dot(self.J)
