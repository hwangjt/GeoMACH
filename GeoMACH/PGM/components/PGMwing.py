"""
GeoMACH wing class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse, scipy.sparse.linalg
from collections import OrderedDict

from GeoMACH.BSE.BSEmodel import BSEmodel
from GeoMACH.PGM.components.PGMprimitive import PGMprimitive
from GeoMACH.PGM.core.PGMface import PGMface


class PGMwing(PGMprimitive):
    """ Wing component """

    def __init__(self, num_x=1, num_z=1,
                 left_closed=False, right_closed=False, blunt_te=False):
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
        self._blunt_te = blunt_te

        self._ax1 = 3
        self._ax2 = 2

        self.faces['upp'] = PGMface(num_x, num_z)
        self.faces['low'] = PGMface(num_x, num_z)

        self.funcs = {}

    def assemble_sizes(self, bse):
        super(PGMwing, self).assemble_sizes(bse)

        if self._bse is None:
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
            if not self._blunt_te:
                face.set_diff_surf(False, ind_i=-ind, ind_u=2*ind)
                face.set_diff_edge(True, 'u' + str(ind), ind_i=-ind)
            #C0 left edge
            if not self._left_closed:
                face.set_diff_surf(False, ind_j=-1, ind_v=2)
                face.set_diff_edge(True, 'v1', ind_j=-1)
                if not self._blunt_te:
                    face.set_diff_corner(False, ind_i=-ind, ind_j=-1)
                #face.set_diff_corner(False, ind_i=ind, ind_j=-1)
            #C0 right edge
            if not self._right_closed:
                face.set_diff_surf(False, ind_j=0, ind_v=0)
                face.set_diff_edge(True, 'v0', ind_j=0)
                if not self._blunt_te:
                    face.set_diff_corner(False, ind_i=-ind, ind_j=0)
                #face.set_diff_corner(False, ind_i=ind, ind_j=0)

    def set_airfoil(self, filename='naca0012', blunt_thk=0.0, blunt_pos=0.95, bunch_LE=1.0, bunch_TE=1.0):
        if filename[:4]=='naca' or filename[:4]=='NACA':
            airfoils = self._get_airfoil_naca(filename[4:])
        else:
            airfoils = self._get_airfoil_file(filename)

        for name in ['upp', 'low']:
            ms = self.faces[name]._num_cp_list['u']
            ns = self.faces[name]._num_pt_list['u']
            nP = sum(ns) + 1

            P = self._get_P(nP, airfoils, name, bunch_LE, bunch_TE)
            P[:, 0] /= numpy.max(P[:, 0])
            P[:, 1] -= numpy.linspace(P[0,1], P[-1,1], P.shape[0])

            if name == 'upp':
                sign = 1.0
            elif name == 'low':
                sign = -1.0
            t, p, x = blunt_thk, blunt_pos, P[:, 0]
            P[:, 1] += sign * t * (-2*(x/p)**3 + 3*(x/p)**2) * (x < p)
            P[:, 1] += sign * t * numpy.sqrt(1 - (numpy.maximum(x,2*p-1) - p)**2/(1-p)**2) * (x >= p)

            Q = self._get_Q(ms, ns, P)
            for j in range(self.faces[name]._num_cp_total['v']):
                self._shapes[name][:,j,:] = Q[:,:]
                self._shapes[name][:,j,2] = 0.0

    def _get_Q(self, ms, ns, P0):
        nsurf = ns.shape[0]

        Ps = []
        for i in xrange(nsurf):
            P = numpy.zeros((ns[i]+1,4,3), order='F')
            for j in xrange(4):
                P[:,j,:] = P0[sum(ns[:i]):sum(ns[:i+1])+1,:]
                P[:,j,2] = j
            Ps.append(P)

        bse = BSEmodel(Ps)
        for i in xrange(nsurf):
            bse.set_bspline_option('num_cp', i, 'u', ms[i]+1)
            bse.set_bspline_option('num_pt', i, 'u', ns[i]+1)
            bse.set_bspline_option('num_pt', i, 'v', 4)
        bse.assemble()
        cp = bse.vec['cp_str'].array
        pt = bse.vec['pt_str'].array
        jac = bse.jac['d(pt_str)/d(cp_str)']

        for i in xrange(nsurf):
            bse.vec['pt_str'].surfs[i][:,:,:] = Ps[i]
        mtx = jac.T.dot(jac)
        rhs = jac.T.dot(pt)
        for dim in xrange(3):
            cp[:, dim] = scipy.sparse.linalg.cg(mtx, rhs[:, dim])[0]

        Q = numpy.zeros((sum(ms) + 1,3),order='F')
        for i in xrange(nsurf):
            Q[sum(ms[:i]):sum(ms[:i+1])+1,:] = \
                bse.vec['cp_str'].surfs[i][:,0,:]
        return Q

    def _get_P(self, nP, airfoils, name, bunch_LE, bunch_TE):
        def bunch_start(ind, a):
            ind[:] = ind**a
        def bunch_end(ind, a):
            ind[:] = 1 - (1-ind)**a

        airfoil = airfoils[name]

        P = numpy.zeros((airfoil.shape[0],4,3),order='F')
        for j in range(4):
            P[:,j,:2] = airfoil[:,:]
            P[:,j,2] = j
        bse = BSEmodel([P])

        bse.set_bspline_option('num_pt', 0, 'u', 
                               airfoil.shape[0])
        bse.set_bspline_option('num_cp', 0, 'u', 
                               int(airfoil.shape[0]/3))
        bse.set_bspline_option('num_pt', 0, 'v', 4)
        bse.set_bspline_option('num_cp', 0, 'v', 4)
        bse.assemble()
        cp = bse.vec['cp_str'].array
        pt = bse.vec['pt_str'].array
        jac = bse.jac['d(pt_str)/d(cp_str)']

        bse.vec['pt_str'].surfs[0][:, :, :] = P
        mtx = jac.T.dot(jac)
        rhs = jac.T.dot(pt)
        for dim in xrange(3):
            cp[:, dim] = scipy.sparse.linalg.cg(mtx, rhs[:, dim])[0]
        fit = numpy.array(cp)

        bse.set_bspline_option('num_pt', 0, 'u', nP)
        bse.assemble()
        cp = bse.vec['cp_str'].array
        pt = bse.vec['pt_str'].array

        cp[:, :] = fit[:, :]
        surfs = numpy.zeros(nP).astype(int)
        ind_u = numpy.linspace(0, 1, nP)
        if name == 'upp':
            bunch_start(ind_u, bunch_TE)
            bunch_end(ind_u, bunch_LE)
        elif name == 'low':
            bunch_start(ind_u, bunch_LE)
            bunch_end(ind_u, bunch_TE)
        ind_v = numpy.zeros(nP)
        bse.add_jacobian('new', surfs, ind_u, ind_v, 3)
        bse.apply_jacobian('new', 'd(new)/d(cp_str)', 'cp_str')

        return bse.vec['new'].array[:, :]

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
        if pos_cmb != 0:
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

    def add_thk_con(self, name, urange, vrange, factor):
	self.funcs[name] = WingThicknessFunction(self, name, urange, vrange, factor)



class WingFunction(object):

    def __init__(self, comp, name, urange, vrange, factor):
	self.bse = comp._bse
	self.comp = comp
	self.name = name
	self.num_u, self.num_v = len(urange), len(vrange)
	self.factor = factor

	num_u, num_v = self.num_u, self.num_v
	self.pt = numpy.zeros(2*num_u*num_v*3)
	self.pt_array = self.pt.reshape((2,num_u,num_v,3), order='F')

	ni, nj = self.num_u, self.num_v
	locations = {'u':numpy.zeros((ni,nj)), 'v':numpy.zeros((ni,nj))}
	for i in range(ni):
	    for j in range(nj):
		locations['u'][i,j] = urange[i]
		locations['v'][i,j] = vrange[j]

	face = self.comp.faces['upp']

        increments = {'upp': {'u': None, 'v': None}, 'low': {'u': None, 'v': None}}
	for f in ['upp', 'low']:
	    for d in ['u', 'v']:
		n = face._num_surf[d]
		increments[f][d] = numpy.zeros(n+1)
		for i in xrange(n+1):
		    increments[f][d][i] = sum(face._num_cp_list[d][:i]) / sum(face._num_cp_list[d])
	    if f=='upp':
		increments[f]['u'][:] = 1 - increments[f]['u'][::-1]

	surf = {'upp':numpy.zeros((ni,nj), dtype=int, order='F'), 'low':numpy.zeros((ni,nj), dtype=int, order='F')}
	locs = {'upp':{'u':numpy.zeros((ni,nj), order='F'), 'v':numpy.zeros((ni,nj), order='F')}, 'low':{'u':numpy.zeros((ni,nj), order='F'), 'v':numpy.zeros((ni,nj), order='F')}}


	loc_face = {'u': None, 'v': None}
	for f in ['upp', 'low']:
	    for i in range(ni):
		for j in range(nj):
		    if f=='upp':
			loc_face['u'] = 1-locations['u'][i,j]
			loc_face['v'] = locations['v'][i,j]
		    else:
			loc_face['u'] = locations['u'][i,j]
			loc_face['v'] = locations['v'][i,j]
		    loc_surf = {'u':-1, 'v':-1}
		    for d in ['u', 'v']:
			n = face._num_surf[d]
			for k in range(n):
			    if increments[f][d][k] <= loc_face[d] <= increments[f][d][k+1]:
				locs[f][d][i,j] = (loc_face[d] - increments[f][d][k]) / \
				    (increments[f][d][k+1] - increments[f][d][k])
				loc_surf[d] = k
		    if loc_surf['u'] == -1 or loc_surf['v'] == -1:
			raise Exception('Invalid thickness constraint locations')
		    ind = 0 if f=='upp' else 1
		    surf[f][i,j] = self.comp.faces.values()[ind]._surf_indices[loc_surf['u'], loc_surf['v']]

	surf_flat = numpy.zeros(2*ni*nj)
	locs_u_flat = numpy.zeros(2*ni*nj)
	locs_v_flat = numpy.zeros(2*ni*nj)

	a = surf['upp'].flatten(order='F')
	b = surf['low'].flatten(order='F')
	for k in range(ni*nj):
	    surf_flat[2*k] = a[k]
	    surf_flat[2*k + 1] = b[k]

	c = locs['upp']['u'].flatten(order='F')
	d = locs['low']['u'].flatten(order='F')
	for k in range(ni*nj):
	    locs_u_flat[2*k] = c[k]
	    locs_u_flat[2*k + 1] = d[k]

	e = locs['upp']['v'].flatten(order='F')
	f = locs['low']['v'].flatten(order='F')
	for k in range(ni*nj):
	    locs_v_flat[2*k] = e[k]
	    locs_v_flat[2*k + 1] = f[k]

	self.bse.add_jacobian('constr_pts', surf_flat, locs_u_flat, locs_v_flat, ndim=3)
	J = self.bse.jac['d(constr_pts)/d(cp_str)']

#	self.bse.apply_jacobian('constr_pts', 'd(constr_pts)/d(cp_str)', 'cp_str')
#	self.bse.vec['thickness'].export_tec_scatter()
#	exit()

	self.dpt_dcp = scipy.sparse.bmat(
            [
                [J, None, None],
                [None, J, None],
                [None, None, J]
            ],
            format = 'csc')
	

    def initialize(self):
	self.compute_all()
	self.func0[:] = self.func[:]

    def get_func(self):
	self.compute_all()
	return self.factor**2 * self.func0[:] - self.func[:]



class WingThicknessFunction(WingFunction):

    def __init__(self, comp, name, urange, vrange, factor):
	super(WingThicknessFunction, self).__init__(comp, name, urange, vrange, factor)
	self.size = self.num_u * self.num_v
	num_u, num_v = self.num_u, self.num_v

	self.func = numpy.zeros(num_u*num_v)
	self.func_array = self.func.reshape((num_u, num_v), order='F')

	self.func0 = numpy.array(self.func)

    def compute_all(self):
	self.pt[:] = self.dpt_dcp.dot(self.bse.vec['cp_str'].array.reshape(3*self.bse.vec['cp_str'].size, order='F'))
	num_u, num_v = self.num_u, self.num_v

	self.func_array[:,:] = 0.0
	for k in xrange(3):
	    self.func_array[:,:] += (self.pt_array[0,:,:,k] - self.pt_array[1,:,:,k])**2
#	self.func_array[:,:] += self.pt_array[0,:,:,1]


    def get_jacobian(self):
	self.compute_all()
	num_u, num_v = self.num_u, self.num_v

	nD = 6 * num_u * num_v
	Da = numpy.zeros((2,num_u,num_v,3), order = 'F')
	Di = numpy.zeros((2,num_u,num_v,3), dtype=int, order='F')
	Dj = numpy.zeros((2,num_u,num_v,3), dtype=int, order='F')

	pt_indices = numpy.array(numpy.linspace(0, 2*num_u*num_v*3-1, 2*num_u*num_v*3), int).reshape((2,num_u,num_v,3), order='F')

#	for f in [0]:
#	    for k in [1]:
#		Da[f,:,:,k] = -1
#		Di[f,:,:,k] = numpy.linspace(0, num_u*num_v-1, num_u*num_v).reshape((num_u,num_v), order='F')
#		Dj[f,:,:,k] = pt_indices[f,:,:,k]
	for f in xrange(2):
	    for k in xrange(3):
		Da[f,:,:,k] = -2 * (self.pt_array[f,:,:,k] - self.pt_array[1-f,:,:,k])
		Di[f,:,:,k] = numpy.linspace(0, num_u*num_v-1, num_u*num_v).reshape((num_u,num_v), order='F')
		Dj[f,:,:,k] = pt_indices[f,:,:,k]

	Da = Da.reshape(2*num_u*num_v*3, order='F')
	Di = Di.reshape(2*num_u*num_v*3, order='F')
	Dj = Dj.reshape(2*num_u*num_v*3, order='F')
	df_dpt = scipy.sparse.csr_matrix((Da, (Di, Dj)), shape=(num_u*num_v, 2*num_u*num_v*3))
        df_dcp = df_dpt * self.dpt_dcp
	return df_dcp
