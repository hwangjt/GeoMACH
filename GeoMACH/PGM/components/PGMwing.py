"""
GeoMACH wing class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse
import time
from collections import OrderedDict

from GeoMACH.BSE.BSEmodel import BSEmodel
from GeoMACH.PGM.components.PGMprimitive import PGMprimitive
from GeoMACH.PGM.core.PGMface import PGMface


class PGMwing(PGMprimitive):
    """ Wing component """

    def __init__(self, num_x=1, num_z=1,
                 left_closed=False, right_closed=False):
        """
        Parameters
        ----------
        num_x : ``float``
           Number of surfaces in the chord-wise direction
        num_z : ``float``
           Number of surfaces in the span-wise direction
        left_closed, right_closed : ``bool``
           ``True`` if closed and 
           ``False`` if attached to another component
        """
        super(PGMwing, self).__init__()

        self._num_surf['x'] = num_x
        self._num_surf['z'] = num_z

        self._left_closed = left_closed
        self._right_closed = right_closed

        self._ax1 = 3
        self._ax2 = 2

        self.faces['upp'] = PGMface(num_x, num_z)
        self.faces['low'] = PGMface(num_x, num_z)

    def assemble_sizes(self, bse):
        super(PGMwing, self).assemble_sizes(bse)

        if self._bse is None or True:
            upp = self._shapes['upp']
            low = self._shapes['low']
            upp[:, :, :] = 0.0
            low[:, :, :] = 0.0
            upp[1:-1, :, 1] = 0.05
            low[1:-1, :, 1] = -0.05
            num = upp.shape[0]
            for ind in range(num):
                upp[ind, :, 0] = 1 - ind / (num-1)
                low[ind, :, 0] = ind / (num-1)
        else:
            self.set_airfoil()

    def set_diff(self):
        faces = self.faces.values()
        for ind in xrange(2):
            face = self.faces.values()[ind]
             #C1 Everywhere
            face.set_diff_surf(True)
            #C0 trailing edge
            face.set_diff_surf(False, ind_i=-ind, ind_u=2*ind)
            face.set_diff_edge(True, 'u' + str(ind), ind_i=-1)
            #C0 left edge
            if not self._left_closed:
                face.set_diff_surf(False, ind_j=-1, ind_v=2)
                face.set_diff_edge(True, 'v1', ind_j=-1)
                face.set_diff_corner(False, ind_i=-ind, ind_j=-1)
                #face.set_diff_corner(False, ind_i=ind, ind_j=-1)
            #C0 right edge
            if not self._right_closed:
                face.set_diff_surf(False, ind_j=0, ind_v=0)
                face.set_diff_edge(True, 'v0', ind_j=0)
                face.set_diff_corner(False, ind_i=-ind, ind_j=0)
                #face.set_diff_corner(False, ind_i=ind, ind_j=0)

    def set_airfoil(self, filename='naca0012'):
        if filename[:4]=='naca' or filename[:4]=='NACA':
            airfoil = self._get_airfoil_naca(filename[:4])
        else:
            airfoil = self._get_airfoil_file(filename)

        Qs = {}
        for name in ['upp', 'low']:
            nsurf = self.faces[name].num_surf[0]
            ms = self.faces[name].num_cp_list[0]
            ns = self.faces[name].num_pt_list[0]
            nP = sum(ns) + 1

            P = _get_P(nP, airfoil[name])
            Q = _get_Q(ms, ns, P)
            Qs[name] = Q

        for name in self.shapes:
            for j in range(self.faces[name].num_cp[1]):
                self.shapes[name][:,j,:2] = Qs[name][:,:]

    def _get_P(self, nP, airfoil):
        P = numpy.zeros((airfoil.shape[0],4,3),order='F')
        for j in range(4):
            P[:,j,:2] = airfoil[:,:]
            P[:,j,2] = j

        bse = BSEmodel([P])
        bse.set_bspline_option('num_pt', 0, 'u', 
                               airfoil.shape[0])
        bse.vec['pt_str'].surfs[0][:, :, :] = P
        jac1 = bse.jac['d(df)/d(df_str)']
        jac2 = bse.jac['d(cp)/d(df)']
        jac3 = bse.jac['d(cp_str)/d(cp)']
        jac4 = bse.jac['d(pt_str)/d(cp_str)']
        jac = jac4.dot(jac3.dot(jac2.dot(jac1)))
        jacT = jac.transpose()
        jacTjac = jacT.dot(jac)
        bse.assemble()

        oml0 = PUBS.PUBS([P])
        oml0.edgeProperty(0,2,0,nP)
        oml0.updateEvaluation()

        P = numpy.zeros((nP,2),order='F')
        for i in range(nP):
            P[i,:] = oml0.P[oml0.getIndex(0,i,0,0),:2]
        return P

    def _get_airfoil_naca(self, naca):
        num = 50
        max_cmb = int(naca[0]) / 100.0
        pos_cmb = int(naca[1]) / 10.0
        thickness = int(naca[2:4]) / 100.0
        x = 0.5 * (1 - numpy.cos(numpy.pi * 
                                 numpy.linspace(0, 1, num)))
        y_sym = thickness/0.2 * (0.2969*x**0.5 - 0.1260*x - 
                                 0.3516*x**2 + 0.2843*x**3 - 
                                 0.1036*x**4)
        y_cmb = numpy.zeros(num)
        if p != 0:
            y1 = max_cmb * x / pos_cmb**2 * \
                 (2 * pos_cmb - x)
            y2 = max_cmb * (1-x) / (1-pos_cmb)**2 * \
                 (1 + x - 2 * pos_cmb)
            for i in range(num):
                if x[i] < p:
                    y_cmb[i] = y1[i]
                else:
                    y_cmb[i] = y2[i]
        upper = numpy.zeros((num, 2), order='F')
        lower = numpy.zeros((num, 2), order='F')
        upper[:,0] = x[::-1]
        lower[:,0] = x
        upper[:,1] = y_sym[::-1] + y_cmb[::-1]
        lower[:,1] = -y_sym + y_cmb
        return {'upp': upper, 'low': lower}

    def _get_airfoil_file(self, filename):
        path = __import__(__name__).__file__
        index_slash = path[::-1].index('/')
        path = path[:-index_slash]
        data = numpy.genfromtxt(path+'PGM/airfoils/'+filename)

        if data[0,0] > data[1,0]:
            mark = numpy.argmin(data,0)[0]
            upper = data[:mark+1,:]
            lower = data[mark:,:]
        else:
            for i in range(data.shape[0]-1):
                if abs(data[i+1,0]-data[i,0]) > 0.8:
                    mark = i
            upper = data[mark::-1,:]
            lower = data[mark+1:,:]
        return {'upp': upper, 'low': lower}
        
            
