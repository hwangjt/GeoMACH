# Davide Ivaldi, John Hwang, Ney Secco, 2015
#Using junction design variables

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone


class Trussbraced(PGMconfiguration):

    def _define_comps(self):
        self.comps['lwing'] = PGMwing(num_x=9, num_z=7, left_closed=True, blunt_te=True)
	self.comps['lstrut'] = PGMwing(num_x=5, num_z=4, left_closed=True, blunt_te=True)
 	self.comps['lv'] = PGMwing(blunt_te=True)

        self.comps['lwing_t'] = PGMtip(self, 'lwing', 'left', 0.01)
        #self.comps['lstrut_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [4,1], 'lstrut', 'left', fweight=[15,10,9,0,0,0], mweight=[0, 0, 0, 0.2, 0.2, 13.5])
	self.comps['lstrut_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [4,1], 'lstrut', 'left', fweight=[10,10,10,10,10,10], mweight=[5, 5, 5, 5, 5, 10])
	#self.comps['lstrut_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [4,1], 'lstrut', 'left', fweight=10, mweight=5)
	self.comps['lv_lwing'] = PGMjunction(self, 'lwing', 'low', 'S', [1,2], 'lv', 'left')#, fweight=0.00000001)
	self.comps['lv_lstrut'] = PGMjunction(self, 'lstrut', 'upp', 'S', [1,1], 'lv', 'right', fweight=3)

    def _define_params(self):

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['pos'].params['lin'] = PGMparameter(3, 3, order_u=2)
        lwing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.52,1])
	lwing['shY','upp'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)
	lwing['shY','low'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)

        lstrut = self.comps['lstrut'].props
	lstrut['pos'].params[''] = PGMparameter(1, 3)
	lstrut['pos'].params['lin'] = PGMparameter(2, 3)
	lstrut['scl'].params[''] = PGMparameter(1, 1)
	lstrut['scl'].params['2'] = PGMparameter(2, 3)
	lstrut['nor'].params[''] = PGMparameter(2, 1)
	lstrut['rot'].params[''] = PGMparameter(3, 3, pos_u=[0,0.15,1])
	lstrut['shY','upp'].params[''] = PGMparameter(10, 3, order_u=4, order_v=3)
	lstrut['shY','low'].params[''] = PGMparameter(10, 3, order_u=4, order_v=3)

        lv = self.comps['lv'].props
	lv['pos'].params[''] = PGMparameter(1, 3)
	lv['pos'].params['lin'] = PGMparameter(2, 3)
	lv['scl'].params[''] = PGMparameter(2, 1)
	lv['nor'].params[''] = PGMparameter(3, 1)
	lv['rot'].params[''] = PGMparameter(2, 3, pos_u=[0,1])
	lv['shY','upp'].params[''] = PGMparameter(10, 2, order_u=4, order_v=2)
	lv['shY','low'].params[''] = PGMparameter(10, 2, order_u=4, order_v=2)

	lstrut_lwing = self.comps['lstrut_lwing'].props
	lstrut_lwing['shN',''].params[''] = PGMparameter(7,7)

    def _compute_params(self):

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([0,0,0])
        lwing['pos'].params['lin'].val([[0,0,0],[3.106,-0.3932,15.015],[6.3416,-0.65567,25.8989]])
        lwing['scl'].params[''].val([3.4,2.3,1.115])

        lstrut = self.comps['lstrut'].props
	lstrut['pos'].params[''].val([2,-4.524,0])
	strut_factor = 0.95
	#lstrut['pos'].params['lin'].val([[0,0,0],[1.5,3.9,14.753]])
	lstrut['pos'].params['lin'].val([[0,0,0],[1.5*strut_factor,4.00005*strut_factor,15.153*strut_factor]])
	lstrut['scl'].params[''].val([1.6])
	lstrut['scl'].params['2'].val([[0,0,0],[0,0,0]])
	lstrut['nor'].params[''].val([0.0,0.0])
	lstrut['rot'].params[''].val([[0,0,0],[-30,0,0],[-30,0,0]])

        lv = self.comps['lv'].props
	lv['pos'].params[''].val([3.1,-1.92,9.376])
	lv['pos'].params['lin'].val([[0,0,0],[-0.2,1.5,0]])
	#lv['pos'].params['lin'].val([[0,0,0],[-0.2,1.45,0]])
	lv['scl'].params[''].val([0.9,0.9])
	lv['nor'].params[''].val([0,0.,0.])
	lv['rot'].params[''].val([[-15-90,2,0],[-90,-5,0]])

	lstrut_lwing = self.comps['lstrut_lwing'].props
	lstrut_lwing['shN',''].params[''].data[:,:] = 0

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

	import numpy
	ary = numpy.array

	factor = 2
#        comps['lwing'].faces['upp'].set_option('num_pt', 'u', factor*ary([24,24,12,12,12,12,12,12,30], both=False)
        comps['lwing'].faces['upp'].set_option('num_cp', 'u', factor*ary([4,10, 2, 4, 4, 4, 4, 3, 3]))
        comps['lwing'].faces['upp'].set_option('num_cp', 'v', factor*ary([21,3,3,9,5,4,32]))
	#comps['lstrut'].faces['upp'].set_option('num_cp', 'u', factor*ary([6,16,12,2,4]))
	#comps['lstrut'].faces['upp'].set_option('num_cp', 'u', factor*ary([6,8,20,2,4]))
	#comps['lstrut'].faces['upp'].set_option('num_cp', 'u', factor*ary([6,4,12,4,6))
	comps['lstrut'].faces['upp'].set_option('num_cp', 'u', factor*ary([2,6,6,2,2]))
	comps['lstrut'].faces['upp'].set_option('num_cp', 'v', factor*ary([24,3,3,12]))
	comps['lv'].faces['upp'].set_option('num_cp', 'v', factor*ary([8]))

#        comps['lwing'].faces['upp'].set_option('num_pt', 'v', 3*ary([20,4,4,8,5,4,32]), both=False)
#	comps['lstrut'].faces['upp'].set_option('num_pt', 'u', 3*ary([6,8,12,8,6]), both=False)
#	comps['lstrut'].faces['upp'].set_option('num_pt', 'v', 3*ary([27,4,4,12]), both=False)
#	comps['lv'].faces['upp'].set_option('num_pt', 'v', 3*ary([8]), both=False)
 #       comps['lwing'].faces['upp'].set_option('num_pt', 'u', 3*factor*ary([4,4, 2, 4, 4, 4, 4, 4, 4]), both=False)
        

    def _define_dvs(self):
        dvs = self.dvs
	dvs['shape_wing_upp'] = PGMdv((10,6)).set_identity_param('lwing', ('shY', 'upp'), '')
	dvs['shape_wing_low'] = PGMdv((10,6)).set_identity_param('lwing', ('shY', 'low'), '')
	dvs['shape_strut_upp'] = PGMdv((10,3)).set_identity_param('lstrut', ('shY', 'upp'), '')
	dvs['shape_strut_low'] = PGMdv((10,3)).set_identity_param('lstrut', ('shY', 'low'), '')
	dvs['shape_vstrut_upp'] = PGMdv((10,2)).set_identity_param('lv', ('shY', 'upp'), '')
	dvs['shape_vstrut_low'] = PGMdv((10,2)).set_identity_param('lv', ('shY', 'low'), '')
	dvs['lstrut_lwing_normal'] = PGMdv((7,7)).set_identity_param('lstrut_lwing', ('shN', ''), '')

if __name__ == '__main__':

    pgm = Trussbraced()
    bse = pgm.initialize()

    pgm.comps['lwing'].set_airfoil('rae2822.dat', blunt_thk=0.002, blunt_pos=0.999, bunch_LE=1.2, bunch_TE=1.5)
    pgm.comps['lstrut'].set_airfoil('naca0012', blunt_thk=0.002, blunt_pos=0.9995, bunch_LE=1.2, bunch_TE=2.0)
    pgm.comps['lv'].set_airfoil('naca0010', blunt_thk=0.002, blunt_pos=0.9995, bunch_LE=1.2, bunch_TE=2.0)
    pgm.dvs['lstrut_lwing_normal'].data[:,:] = 0.0 #Select normal perturbation for wing-fuselage junction

    pgm.compute_all()
    pgm.compute_normals() #We have to compute normals with the baseline geometry
    pgm.compute_all()

    bse.vec['pt_str'].export_tec_str()
    #bse.vec['df'].export_tec_scatter()
    bse.vec['cp'].export_tec_scatter()
    #bse.vec['pt'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()

    import os
    #Rename the geometry output file
    os.system('mv pt_str_surf.dat modified.dat')
    os.system('tec360 layout_le.lay')
    os.system('mv cp_str.igs modified.igs')
