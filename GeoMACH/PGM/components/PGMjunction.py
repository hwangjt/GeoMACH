"""
GeoMACH junction class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.BSE.BSEmodel import BSEmodel
from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components.PGMinterpolant import PGMinterpolant
from GeoMACH.PGM.core.PGMface import PGMface


class PGMjunction(PGMinterpolant):
    """ Junction component """

    def __init__(self, config, fcomp, face, 
                 vdir, loc, mcomp, side,
                 fweight=0, mweight=0):
        """
        Parameters
        ----------
        config : ``PGMconfiguration``
           Pointer to the configuration containing this component
        fcomp : ``str``
           Name of the 'female' component
           (the one to which the wing is attached)
        face : ``str``
           Name of the face on the female component
           to which the wing is attached
        vdir : ``str``
           Direction the v axis (of the face of the female component)
           is pointing: 'E', 'N', 'W', or 'S'
        loc : ``list`` of 2 ``int``s
           The 'i' and 'j' indices of the northwest-most face
           of the junction when viewed
           with the wing's upper surface facing up
        mcomp : ``str``
           Name of the 'male' component
           (the one being attached to the female component)
        side : ``str``
           The side of wing being attached:
           'right' for the v=0 side,
           'left' for the v=1 side
        fweight : ``float``
           The weight applied to the tangent vectors
           of the female component when interpolating
        mweight : ``float``
           The weight given to the tangent vectors
           of the male component when interpolating
        """
        super(PGMjunction, self).__init__()

        self._fcomp = config.comps[fcomp].faces[face]
        self._loc = {'u': loc[0], 'v': loc[1]}
        self._mcomp = config.comps[mcomp]
        self._side = side

        #Check if the user gave a single weight for all the edges or if he specified one weight for each edge
	self._fweight = numpy.zeros(6)
	self._mweight = numpy.zeros(6)
	self._fweight[:] = fweight
	self._mweight[:] = mweight

        self._num_surf_wing = self._mcomp.faces['upp']._num_surf['u']

        self.faces[''] = PGMface(3, 2 + self._num_surf_wing)
        self.faces['']._surf_indices[1, :] = -1

        if vdir == 'E':
            self._rotate = lambda P: P
            self._flip = lambda nu, nv: [nu,nv]
        elif vdir == 'N':
            self._rotate = lambda P: numpy.swapaxes(P,0,1)[::-1,:]
            self._flip = lambda nu, nv: [nv[::-1],nu]
        elif vdir == 'W':
            self._rotate = lambda P: P[::-1,::-1]
            self._flip = lambda nu, nv: [nu[::-1],nv[::-1]]
        elif vdir == 'S':
            self._rotate = lambda P: numpy.swapaxes(P,0,1)[:,::-1]
            self._flip = lambda nu, nv: [nv,nu[::-1]]

    def set_diff(self):
        face = self.faces['']
        face.set_diff_surf(True)
        
        for ind_j in xrange(1, 1 + self._num_surf_wing):
            face.set_diff_surf(False, ind_i=0, ind_j=ind_j, ind_u=2)
            face.set_diff_surf(False, ind_i=-1, ind_j=ind_j, ind_u=0)

        face.set_diff_surf(False, ind_i=0, ind_j=0, ind_u=2, ind_v=2)
        face.set_diff_surf(False, ind_i=-1, ind_j=0, ind_u=0, ind_v=2)
        face.set_diff_surf(False, ind_i=0, ind_j=-1, ind_u=2, ind_v=0)
        face.set_diff_surf(False, ind_i=-1, ind_j=-1, ind_u=0, ind_v=0)

    def set_hidden_surfaces(self):
        loc = self._loc
        fInds = self._rotate(self._fcomp._surf_indices)
        for ind_v in xrange(2 + self._num_surf_wing):
            for ind_u in xrange(2):
                isurf = fInds[ind_u + loc['u'], ind_v + loc['v']]
                self._bse.hidden[isurf] = True

    def compute(self, name):
        loc = self._loc

        fu = self._fcomp._num_cp_list['u']
        fv = self._fcomp._num_cp_list['v']
        fu,fv = self._flip(fu,fv)
        fu1 = sum(fu[:loc['u']])
        fu2 = sum(fu[:loc['u']+2])
        fv1 = sum(fv[:loc['v']])
        fv2 = sum(fv[:loc['v']+2+self._num_surf_wing])
        fFace_inds = self._rotate(self._fcomp.vec_inds['cp_bez'])
        fFace_inds = fFace_inds[fu1:fu2+1,fv1:fv2+1]
        
        mcomp = self._mcomp
        if self._side == 'right':
            W = numpy.zeros((4,2,3),order='F')
            E = numpy.zeros((4,2,3),order='F')
            N = mcomp.faces['upp'].vec_inds['cp_prim'][::-1,:2,:]
            S = mcomp.faces['low'].vec_inds['cp_prim'][:,:2,:]
        elif self._side == 'left':
            W = numpy.zeros((4,2,3),order='F')
            E = numpy.zeros((4,2,3),order='F')
            N = mcomp.faces['upp'].vec_inds['cp_prim'][:,-1:-3:-1,:]
            S = mcomp.faces['low'].vec_inds['cp_prim'][::-1,-1:-3:-1,:]

        num_u = [fu[loc['u']] + 1, 4, fu[loc['u']+1] + 1]
        num_v = [fv[loc['v']] + 1, 
                 sum(fv[loc['v']+1:loc['v'] + 1 + self._num_surf_wing]) + 1, 
                 fv[loc['v'] + 1 + self._num_surf_wing] + 1]
        nu0 = sum(num_u) - 2
        nv0 = sum(num_v) - 2
        fInds = -numpy.ones((nu0, nv0, 3), dtype=int, order='F')
        fInds[:num_u[0],:] = fFace_inds[:num_u[0],:]
        fInds[-num_u[2]:,:] = fFace_inds[-num_u[2]:,:]

        nD = 3 * 2 * nv0 + 3 * 2 * (nu0 - 2)
        nD += 3 * 2
        if num_u[1] != 1:
            nD += 3 * 2
        nD += 3 * 2 * (num_v[1] - 2)
        if num_u[1] != 1:
            nD += 3 * 2 * (num_u[1] - 2)
        nD += 3 * 4 * 2 * (num_u[0] - 2)
        nD += 3 * 4 * 2 * (num_u[2] - 2)
        nD += 3 * 4 * (num_v[0] - 2)
        nD += 3 * 4 * (num_v[2] - 2)
        if num_u[1] != 1:
            nD += 3 * 4 * (num_v[0] - 2)
            nD += 3 * 4 * (num_v[2] - 2)

        Da, Di, Dj = PGMlib.computejunctionwireframe(nD, nu0, nv0, num_u[0], num_u[1], num_u[2], num_v[0], num_v[1], num_v[2], self._fweight, self._mweight, W, E, N, S, fInds, self.faces[''].vec_inds['cp_bez'])
        Da = Da * (-1 != Dj)
        Dj = numpy.maximum(0, Dj)

        Das, Dis, Djs = [Da], [Di], [Dj]
        if name == 'cp_bez':
            return Das, Dis, Djs
        elif name == 'cp_coons':
            nD = 0
            for i in range(3):
                for j in range(3):
                    if (num_u[1] != 1) or (i !=1):
                        nD += 3 * 8 * (num_u[i]-2) * (num_v[j]-2)
        
            Da, Di, Dj = PGMlib.computejunctioncoons(nD, nu0, nv0, num_u[0], num_u[1], num_u[2], num_v[0], num_v[1], num_v[2], self.faces[''].vec_inds['cp_coons'])
            Das.append(Da)
            Dis.append(Di)
            Djs.append(Dj)
            return Das, Dis, Djs
        elif name == 'cp_prim':
            return [], [], []

            
