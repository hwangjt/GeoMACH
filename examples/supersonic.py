# Davide Ivaldi, 2014

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone


class Supersonic(PGMconfiguration):

    def _define_comps(self):
        self.comps['fuse'] = PGMbody(num_x=8, num_y=6, num_z=4)
        self.comps['lwing'] = PGMwing(num_x=5, num_z=7, left_closed=True)
        self.comps['pylon'] = PGMwing()
        self.comps['nac'] = PGMshell(num_x=3, num_y=1, num_z=4)

        self.comps['fuse_f'] = PGMcone(self, 'fuse', 'front', 4)
        self.comps['fuse_r'] = PGMcone(self, 'fuse', 'rear', 2)
        self.comps['lwing_t'] = PGMtip(self, 'lwing', 'left', 0.1)
        self.comps['lwing_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [0,1], 'lwing', 'right')
        self.comps['pylon_fuse'] = PGMjunction(self, 'fuse', 'top', 'E', [1,5], 'pylon', 'right')
        self.comps['pylon_nac'] = PGMjunction(self, 'nac', 'bt0', 'W', [1,0], 'pylon', 'left')


    def _define_params(self):
        fuse = self.comps['fuse'].props
	fuse['pos'].params[''] = PGMparameter(2, 3)
	fuse['pos'].params['nose'] = PGMparameter(2, 3, pos_u=[0,0.45])
	fuse['pos'].params['tail'] = PGMparameter(2, 3, pos_u=[0.85,1.0])
	fuse['scl'].params['rad_nose'] = PGMparameter(4, 1, pos_u=[0,0.05,0.41,0.45], order_u=4)
	fuse['scl'].params['rad_middle'] = PGMparameter(2, 1, pos_u=[0.45,0.85])
	fuse['scl'].params['rad_tail'] = PGMparameter(4, 1, pos_u=[0.85,0.87,0.98,1.0], order_u=4)
	fuse['nor'].params[''] = PGMparameter(1, 3)


        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=2)
        lwing['scl'].params[''] = PGMparameter(2, 1)
        lwing['rot'].params[''] = PGMparameter(2, 3)

	pylon = self.comps['pylon'].props
	pylon['pos'].params[''] = PGMparameter(1, 3)
	pylon['pos'].params['lin'] = PGMparameter(2, 3)
	pylon['scl'].params[''] = PGMparameter(2, 1)
	pylon['nor'].params[''] = PGMparameter(1, 3)

	nac = self.comps['nac'].props
	nac['pos'].params[''] = PGMparameter(1, 3)
	nac['pos'].params['lin'] = PGMparameter(2, 3)
	nac['nor'].params[''] = PGMparameter(1, 1)
	nac['scl'].params[''] = PGMparameter(3, 1, pos_u=[0, 0.8, 1], order_u=3)
	nac['thk'].params[''] = PGMparameter(3, 1)


    def _compute_params(self):
        fuse = self.comps['fuse'].props
	fuse['pos'].params[''].val([[0,0,0],[24,0,0]])
	fuse['pos'].params['nose'].val([[0,-0.9,0],[0,0,0]])
	fuse['pos'].params['tail'].val([[0,0,0],[0,0.8,0]])
	fuse['scl'].params['rad_nose'].val([0.2,0.25,1.2,1.3])
	fuse['scl'].params['rad_middle'].val([1.3,1.28])
	fuse['scl'].params['rad_tail'].val([1.28,1.2,0.4,0.3])
	fuse['nor'].params[''].val([1.0,0.0,1.0])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([12,0.54,1.29])
        lwing['pos'].params['lin'].val([[0,0,0],[4,0,1],[12,0,8]])
        lwing['scl'].params[''].val([10,1.2])
	lwing['rot'].params[''].val([[-40,0,0],[0,0,0]])

        pylon = self.comps['pylon'].props
        pylon['pos'].params[''].val([19.8,1.35,0])
        pylon['pos'].params['lin'].val([[0,0,0],[0.1,0.1,0]])
        pylon['scl'].params[''].val([2,1.9])
	pylon['nor'].params[''].val([1,0,0])

        nac = self.comps['nac'].props
        nac['pos'].params[''].val([19,2,0])
        nac['pos'].params['lin'].val([[0,0,0],[5,0,0]])
        nac['nor'].params[''].val([1])
        nac['scl'].params[''].val([0.4,0.6,0.5])
        nac['thk'].params[''].val([0.08,0.2,0.08])

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

        comps['fuse'].faces['lft'].set_option('num_cp', 'u', [4,5,4,4,4,4])
        comps['fuse'].faces['lft'].set_option('num_cp', 'v', [36,4,10,6,6,6,10,8])
        comps['lwing'].faces['upp'].set_option('num_cp', 'u', [10,6,6,6,10])
        comps['nac'].faces['bt0'].set_option('num_cp', 'v', [6,10,6])
        comps['nac'].faces['bt1'].set_option('num_cp', 'v', [6,10,6])
        

if __name__ == '__main__':

    pgm = Supersonic()
    bse = pgm.initialize()

    pgm.comps['lwing'].set_airfoil('n64206.dat')
    pgm.comps['pylon'].set_airfoil('naca0010')
    pgm.compute_all()

    bse.vec['pt_str'].export_tec_str()
    bse.vec['df'].export_tec_scatter()
    bse.vec['cp'].export_tec_scatter()
    bse.vec['pt'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()
