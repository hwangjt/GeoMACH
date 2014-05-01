from __future__ import division
import numpy, time, scipy.sparse

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Component, Property



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
        props = self.props
        props['scl'] = Property(n,3)
        props['pos'] = Property(n,3)
        props['rot'] = Property(n,3)
        props['ogn'] = Property(n,3)
        props['nor'] = Property(n,3)
        props['flt'] = Property(n,4)

    def computeSections(self):
        nf = len(self.faces)
        n = self.faces.values()[0].num_cp[1]

        pos = self.props['pos'].prop_vec
        nor = self.props['nor'].prop_vec
        rot = self.props['rot'].prop_vec
        ogn = self.props['ogn'].prop_vec
        scl = self.props['scl'].prop_vec

        pos_ind = self.props['pos'].prop_ind
        nor_ind = self.props['nor'].prop_ind
        rot_ind = self.props['rot'].prop_ind
        ogn_ind = self.props['ogn'].prop_ind
        scl_ind = self.props['scl'].prop_ind

        for name in self.faces:
            shX = self.props['shX',name].prop_vec + self.shapes[name][:,:,0]
            shY = self.props['shY',name].prop_vec + self.shapes[name][:,:,1]
            shZ = self.props['shZ',name].prop_vec + self.shapes[name][:,:,2]

            shX_ind = self.props['shX',name].prop_ind
            shY_ind = self.props['shY',name].prop_ind
            shZ_ind = self.props['shZ',name].prop_ind

            ni, nj = self.faces[name].num_cp
            
            cp_ind = self.faces[name].cp_indices

            self.faces[name].cp_array[:,:,:], Da, Di, Dj = PGMlib.computesections2(self.ax1, self.ax2, ni, nj, ni*nj*3*3*4, ogn, nor, pos, rot, scl, shX, shY, shZ, ogn_ind, nor_ind, pos_ind, rot_ind, scl_ind, shX_ind, shY_ind, shZ_ind, cp_ind)

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
