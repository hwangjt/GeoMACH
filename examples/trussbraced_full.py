# Davide Ivaldi, John Hwang, 2014

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone


class Trussbraced(PGMconfiguration):

    def _define_comps(self):
        self.comps['fuse'] = PGMbody(num_x=17, num_y=6, num_z=4)
        self.comps['lwing'] = PGMwing(num_x=7, num_z=7, left_closed=True)
	self.comps['lstrut'] = PGMwing(num_x=4, num_z=4, left_closed=True, right_closed=True)
 	self.comps['lv'] = PGMwing(num_z=4)
        self.comps['ltail'] = PGMwing(left_closed=True)
        self.comps['vtail'] = PGMwing(num_x=5, num_z=4, left_closed=True)

        self.comps['fuse_f'] = PGMcone(self, 'fuse', 'front', 18)
        self.comps['fuse_r'] = PGMcone(self, 'fuse', 'rear', 2)
        self.comps['lwing_t'] = PGMtip(self, 'lwing', 'left', 0.1)
        self.comps['ltail_t'] = PGMtip(self, 'ltail', 'left', 0.1)
        self.comps['vtail_t'] = PGMtip(self, 'vtail', 'left', 0.1)
        self.comps['lwing_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [0,1], 'lwing', 'right', fweight=4, mweight=2)
        self.comps['lstrut_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [4,2], 'lstrut', 'right')
        self.comps['lstrut_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [4,1], 'lstrut', 'left', fweight=3, mweight=3)
	self.comps['lv_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [1,3], 'lv', 'left')
	self.comps['lv_lstrut'] = PGMjunction(self, 'lstrut', 'upp', 'S', [1,0], 'lv', 'right')
        self.comps['vtail_fuse'] = PGMjunction(self, 'fuse', 'top', 'E', [1,10], 'vtail', 'right')
        self.comps['ltail_vtail'] = PGMjunction(self, 'vtail', 'low', 'N', [0,1], 'ltail', 'right')

    def _define_params(self):
        fuse = self.comps['fuse'].props
	fuse['nor'].params[''] = PGMparameter(1, 3)
	fuse['pos'].params[''] = PGMparameter(2, 3)
	fuse['pos'].params['nose'] = PGMparameter(2, 3, pos_u=[0,0.12])
	fuse['pos'].params['tail'] = PGMparameter(2, 3, pos_u=[0.76,1.0])
#	fuse['scl'].params['rad'] = PGMparameter(4, 1, pos_u=[0,0.14,0.75,1.0])
	fuse['scl'].params['rad1'] = PGMparameter(4, 1, pos_u=[0,0.01,0.05,0.12], order_u=4)
	fuse['scl'].params['rad2'] = PGMparameter(2, 1, pos_u=[0.12,0.76])
	fuse['scl'].params['rad3'] = PGMparameter(4, 1, pos_u=[0.76,0.83,0.99,1], order_u=4)
	fuse['scl'].params['tail'] = PGMparameter(2, 3, pos_u=[0.76,1.0])
	fuse['flt'].params['flt1a'] = PGMparameter(4, 2, pos_u=[0.24,0.27,0.33,0.36], pos_v=[0.5,1], order_u=4)
	fuse['flt'].params['flt1b'] = PGMparameter(2, 2, pos_u=[0.36,0.41], pos_v=[0.5,1])
	fuse['flt'].params['flt1c'] = PGMparameter(4, 2, pos_u=[0.41,0.44,0.49,0.52], pos_v=[0.5,1], order_u=4)
	fuse['flt'].params['flt2a'] = PGMparameter(4, 2, pos_u=[0.24,0.27,0.33,0.36], pos_v=[0,0.5], order_u=4)
	fuse['flt'].params['flt2b'] = PGMparameter(2, 2, pos_u=[0.36,0.41], pos_v=[0,0.5])
	fuse['flt'].params['flt2c'] = PGMparameter(4, 2, pos_u=[0.41,0.44,0.49,0.52], pos_v=[0,0.5], order_u=4)

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['scl'].params[''] = PGMparameter(2, 1)
        lwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=3)

        lstrut = self.comps['lstrut'].props
	lstrut['pos'].params[''] = PGMparameter(1, 3)
	lstrut['pos'].params['lin'] = PGMparameter(2, 3)
#	lstrut['pos'].params['curv'] = PGMparameter(3, 3, pos_u=[0.9,0.95,1.0], order_u=3)
	lstrut['scl'].params[''] = PGMparameter(2, 1)
	lstrut['nor'].params[''] = PGMparameter(1, 1)
#	lstrut['rot'].params[''] = PGMparameter(2, 3, pos_u=[0.85,1.0])
#	lstrut['rot'].params[''] = PGMparameter(2, 3, pos_u=[0,1.0])

        lv = self.comps['lv'].props
	lv['pos'].params[''] = PGMparameter(1, 3)
	lv['pos'].params['lin'] = PGMparameter(2, 3)
	lv['scl'].params[''] = PGMparameter(2, 1)
	lv['nor'].params[''] = PGMparameter(1, 1)
	lv['rot'].params[''] = PGMparameter(2, 3, pos_u=[0,1])

        ltail = self.comps['ltail'].props
        ltail['pos'].params[''] = PGMparameter(1, 3)
        ltail['pos'].params['lin'] = PGMparameter(2, 3)
        ltail['scl'].params[''] = PGMparameter(2, 1)

        vtail = self.comps['vtail'].props
        vtail['pos'].params[''] = PGMparameter(1, 3)
        vtail['pos'].params['lin'] = PGMparameter(2, 3)
        vtail['scl'].params[''] = PGMparameter(2, 1)
        vtail['nor'].params[''] = PGMparameter(1, 3)

    def _compute_params(self):
        fuse = self.comps['fuse'].props
	fuse['nor'].params[''].val([1.0,0.0,1.0])
	fuse['pos'].params[''].val([[0,0,0],[36,0,0]])
	fuse['pos'].params['nose'].val([[0,-0.4,0],[0,0,0]])
	fuse['pos'].params['tail'].val([[0,0,0],[0,1.6,0]])
#	fuse['scl'].params['rad'].val([1.4,2,2,0.4])
	fuse['scl'].params['rad1'].val([1,1.2,1.9,2])
	fuse['scl'].params['rad2'].val([2,2])
	fuse['scl'].params['rad3'].val([2,1.7,0.6,0.4])
	fuse['scl'].params['tail'].val([[0,0,0],[-0.3,0,0]])
	fuse['flt'].params['flt1a'].val([[0,0],[0,0],[0.6,0.6],[0.6,0.6]])
	fuse['flt'].params['flt1b'].val([[0.6,0.6],[0.6,0.6]])
	fuse['flt'].params['flt1c'].val([[0.6,0.6],[0.6,0.6],[0,0],[0,0]])
	fuse['flt'].params['flt2a'].val([[0,0],[0,0],[1,1],[1,1]])
	fuse['flt'].params['flt2b'].val([[1,1],[1,1]])
	fuse['flt'].params['flt2c'].val([[1,1],[1,1],[0,0],[0,0]])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([12,1.7,2.8])
        lwing['scl'].params[''].val([3.6,0.8])
        lwing['pos'].params['lin'].val([[0,0,0],[2.5,-0.1,11],[5,-0.8,22]])

        lstrut = self.comps['lstrut'].props
	lstrut['pos'].params[''].val([13.4,-1.6,2.6])
	lstrut['pos'].params['lin'].val([[0,0,0],[1.6,2.6,11.8]])
#	lstrut['pos'].params['curv'].val([[0,0,0],[0,0.25,0],[0,0.5,0]])
	lstrut['scl'].params[''].val([1.8,1.6])
	lstrut['nor'].params[''].val([1.0])
#	lstrut['rot'].params[''].val([[0,0,0],[-30,0,0]])
#	lstrut['rot'].params[''].val([[0,0,0],[-90,0,0]])

        lv = self.comps['lv'].props
	lv['pos'].params[''].val([14.3,-0.12,8.8])
	lv['pos'].params['lin'].val([[0,0,0],[0,1.58,0]])
	lv['scl'].params[''].val([1.1,1.1])
	lv['nor'].params[''].val([1.0])
	lv['rot'].params[''].val([[0,2,0],[0,-2,0]])

        ltail = self.comps['ltail'].props
        ltail['pos'].params[''].val([35.3,6.6,0.25])
        ltail['pos'].params['lin'].val([[0,0,0],[2.6,0,5]])
        ltail['scl'].params[''].val([3.3,1])

        vtail = self.comps['vtail'].props
        vtail['pos'].params[''].val([30.7,2.1,0])
        vtail['pos'].params['lin'].val([[0,0,0],[4.6,5,0]])
        vtail['scl'].params[''].val([5,4.5])
        vtail['nor'].params[''].val([1.0,0.0,0.0])

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

        comps['fuse'].faces['lft'].set_option('num_cp', 'u', [4,4,14,14,4,4])
        comps['fuse'].faces['rgt'].set_option('num_cp', 'v', [85,4,4,4,4,4,4,4,4,4,102,4,4,16,8,4,6])
	comps['vtail'].faces['low'].set_option('num_cp', 'u', [6,4,30,4,4])
	comps['vtail'].faces['low'].set_option('num_cp', 'v', [10,10,10,4])
        comps['lwing'].faces['upp'].set_option('num_cp', 'v', [20,4,4,20,5,4,31])
        comps['lwing'].faces['low'].set_option('num_cp', 'u', [12,12,20,4,4,4,4])
	comps['lstrut'].faces['upp'].set_option('num_cp', 'u', [4,8,12,4])
	comps['lstrut'].faces['upp'].set_option('num_cp', 'v', [4,5,4,4])
        

if __name__ == '__main__':

    pgm = Trussbraced()
    bse = pgm.initialize()

    pgm.comps['lwing'].set_airfoil('rae2822.dat')
    pgm.comps['ltail'].set_airfoil('naca0012')
    pgm.comps['vtail'].set_airfoil('naca0010')
    pgm.compute_all()

    bse.vec['pt_str'].export_tec_str()
    bse.vec['df'].export_tec_scatter()
    bse.vec['cp'].export_tec_scatter()
    bse.vec['pt'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()
