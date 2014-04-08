from __future__ import division
import numpy, time, scipy.sparse

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Component



class Primitive(Component):

    def __init__(self, nx, ny, nz):
        super(Primitive,self).__init__()

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(numpy.zeros(ny,int))
        self.ms.append(numpy.zeros(nz,int))

        self.ns = []
        self.ns.append(numpy.zeros(nx,int))
        self.ns.append(numpy.zeros(ny,int))
        self.ns.append(numpy.zeros(nz,int))

    def declare_properties(self):
        super(Primitive, self).declare_properties()

        n = self.faces.values()[0].num_cp[1]
        props = self.properties
        props['scl'] = [n,3]
        props['pos'] = [n,3]
        props['rot'] = [n,3]
        props['ogn'] = [n,3]
        props['nor'] = [n,3]
        props['flt'] = [n,4]

    def computeSections(self, radii=None):
        nf = len(self.faces)
        n = self.faces.values()[0].num_cp[1]
        p = self.properties
        nQ = 9*n
        for face in self.faces.values():
            nQ += 3 * face.num_cp[0] * face.num_cp[1]

        rot0, Da, Di, Dj = PGMlib.computerotations(self.ax1, self.ax2, n, 9*(n*3-2), p['pos'], p['nor'])
        rot = p['rot']*numpy.pi/180.0 + rot0
        dv_dpos0 = scipy.sparse.csc_matrix((Da,(Di+6*n,Dj+3*n)),shape=(nQ,nQ))

        self.dQs_dv = {}
        counter = 0
        for name in self.faces:
            if radii==None:
                scale = p['scl']
            else:
                scale = radii[name]
            ni, nj = self.faces[name].num_cp
            self.faces[name].cp_array[:,:,:], Da, Di, Dj = PGMlib.computesections(self.ax1, self.ax2, ni, nj, ni*nj*27, counter, p['ogn'], scale, p['pos'], rot, self.shapes[name])
            self.dQs_dv[name] = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nQ))
            self.dQs_dv[name] = self.dQs_dv[name] + self.dQs_dv[name].dot(dv_dpos0)
            counter += 3*ni*nj

    def setDerivatives(self, var, dV0):
        nf = len(self.faces)
        n = self.faces.values()[0].num_cp[1]
        nv = self.dQs_dv.values()[0].shape[1]
        dV = numpy.zeros(nv)
        if var=='scl':
            dV[:3*n] = dV0.T.flatten()
            for name in self.faces:
                ni, nj = self.faces[name].num_cp
                self.faces[name].cp_array[:,:,:] += PGMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[name].dot(dV))
        elif var=='pos':
            dV[3*n:6*n] = dV0.T.flatten()
            for name in self.faces:
                ni, nj = self.faces[name].num_cp
                for i in range(ni):
                    self.faces[name].cp_array[i,:,:] += dV0
                self.faces[name].cp_array[:,:,:] += PGMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[name].dot(dV))
        elif var=='rot':
            dV[6*n:9*n] = dV0.T.flatten()
            for name in self.faces:
                ni, nj = self.faces[name].num_cp
                self.faces[name].cp_array[:,:,:] += PGMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[name].dot(dV)*numpy.pi/180.0)
        else:
            self.properties[var] += dV0
            self.computeQs()
            self.properties[var] -= dV0
