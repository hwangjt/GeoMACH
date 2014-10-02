# Davide Ivaldi, John Hwang, 2014

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone


class Trussbraced(PGMconfiguration):

    def _define_comps(self):
        self.comps['lwing'] = PGMwing(num_x=9, num_z=7, left_closed=True, blunt_te=True)
	self.comps['lstrut'] = PGMwing(num_x=4, num_z=4, left_closed=True, blunt_te=True)
 	self.comps['lv'] = PGMwing(num_z=4, blunt_te=True)

        self.comps['lwing_t'] = PGMtip(self, 'lwing', 'left', 0.01)
        self.comps['lstrut_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [4,2], 'lstrut', 'left', fweight=3, mweight=0.2)
	self.comps['lv_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [1,4], 'lv', 'left', fweight=0.00000001)
	self.comps['lv_lstrut'] = PGMjunction(self, 'lstrut', 'upp', 'S', [1,0], 'lv', 'right', fweight=3)

    def _define_params(self):

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=2)
        lwing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.52,1])


        lstrut = self.comps['lstrut'].props
	lstrut['pos'].params[''] = PGMparameter(1, 3)
	lstrut['pos'].params['lin'] = PGMparameter(2, 3)
	lstrut['scl'].params[''] = PGMparameter(1, 1)
	lstrut['scl'].params['2'] = PGMparameter(2, 3)
#	lstrut['nor'].params[''] = PGMparameter(2, 1)
	lstrut['rot'].params[''] = PGMparameter(2, 3, pos_u=[0,1])

        lv = self.comps['lv'].props
	lv['pos'].params[''] = PGMparameter(1, 3)
	lv['pos'].params['lin'] = PGMparameter(2, 3)
	lv['scl'].params[''] = PGMparameter(2, 1)
	lv['nor'].params[''] = PGMparameter(3, 1)
	lv['rot'].params[''] = PGMparameter(2, 3, pos_u=[0,1])

    def _compute_params(self):

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([0,0,0])
        lwing['pos'].params['lin'].val([[0,0,0],[3.106,-0.3932,15.015],[6.3416,-0.65567,25.8989]])
        lwing['scl'].params[''].val([3.4,2.3,1.115])

        lstrut = self.comps['lstrut'].props
	lstrut['pos'].params[''].val([2,-4.524,0])
	#lstrut['pos'].params['lin'].val([[0,0,0],[1.5,3.9,14.753]])
	lstrut['pos'].params['lin'].val([[0,0,0],[1.5,4.00005,15.153]])
	lstrut['scl'].params[''].val([1.6])
	lstrut['scl'].params['2'].val([[0,0,0],[0,2,0]])
#	lstrut['nor'].params[''].val([0.0,1.0])
	lstrut['rot'].params[''].val([[0,0,0],[-90,0,0]])

        lv = self.comps['lv'].props
	lv['pos'].params[''].val([3.1,-1.92,9.376])
	#lv['pos'].params['lin'].val([[0,0,0],[0,1.5,0]])
	lv['pos'].params['lin'].val([[0,0,0],[-0.2,1.45,0]])
	lv['scl'].params[''].val([0.9,0.9])
	lv['nor'].params[''].val([0,0.,0.])
	lv['rot'].params[''].val([[-15-90,2,0],[-90,-8,0]])

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

        comps['lwing'].faces['upp'].set_option('num_cp', 'u', [15,15, 4, 4, 4, 4, 4, 4, 12], both=False)
        comps['lwing'].faces['upp'].set_option('num_pt', 'u', [24,24,12,12,12,12,12,12,30], both=False)
        comps['lwing'].faces['upp'].set_option('num_cp', 'v', [20,4,4,8,5,4,32])
	comps['lstrut'].faces['upp'].set_option('num_cp', 'u', [6,12,12,4])
	comps['lstrut'].faces['upp'].set_option('num_cp', 'u', [6,8,12,8])
	comps['lstrut'].faces['upp'].set_option('num_cp', 'v', [8,4,4,4])
        

if __name__ == '__main__':

    pgm = Trussbraced()
    bse = pgm.initialize()

    pgm.comps['lwing'].set_airfoil('rae2822.dat', blunt_thk=0.002, blunt_pos=0.9995, bunch_LE=1.0, bunch_TE=1.5)
    pgm.comps['lstrut'].set_airfoil('naca0012', blunt_thk=0.002, blunt_pos=0.9995, bunch_LE=1.0, bunch_TE=2.0)
    pgm.comps['lv'].set_airfoil('naca0010', blunt_thk=0.002, blunt_pos=0.9995, bunch_LE=1.0, bunch_TE=2.0)
    pgm.compute_all()

    bse.vec['pt_str'].export_tec_str()
    bse.vec['df'].export_tec_scatter()
    bse.vec['cp'].export_tec_scatter()
    bse.vec['pt'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()
