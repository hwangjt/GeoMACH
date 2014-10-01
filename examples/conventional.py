# John Hwang

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone


class Conventional(PGMconfiguration):

    def _define_comps(self):
        self.comps['fuse'] = PGMbody(num_x=12, num_y=4, num_z=2)
        self.comps['lwing'] = PGMwing(num_x=4, num_z=4, left_closed=True)
        self.comps['lpylon'] = PGMwing()
        self.comps['lnac'] = PGMshell(num_x=4, num_y=1, num_z=4)
        self.comps['ltail'] = PGMwing(left_closed=True)
        self.comps['vtail'] = PGMwing(num_x=2, left_closed=True)

        self.comps['fuse_f'] = PGMcone(self, 'fuse', 'front', 2)
        self.comps['fuse_r'] = PGMcone(self, 'fuse', 'rear', 2)
        self.comps['lwing_t'] = PGMtip(self, 'lwing', 'left', 0.1)
        self.comps['ltail_t'] = PGMtip(self, 'ltail', 'left', 0.1)
        self.comps['vtail_t'] = PGMtip(self, 'vtail', 'left', 0.1)
        self.comps['lwing_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [2,1], 'lwing', 'right')
        self.comps['lpylon_lwing'] = PGMjunction(self, 'lwing', 'low', 'N', [1,0], 'lpylon', 'right')
        self.comps['lpylon_lnac'] = PGMjunction(self, 'lnac', 'tp0', 'W', [1,0], 'lpylon', 'left')
        self.comps['ltail_fuse'] = PGMjunction(self, 'fuse', 'lft', 'E', [1,9], 'ltail', 'right')
        self.comps['vtail_fuse'] = PGMjunction(self, 'fuse', 'top', 'E', [0,8], 'vtail', 'right')

    def _define_params(self):
        fuse = self.comps['fuse'].props
        fuse['pos'].params[''] = PGMparameter(2, 3)
        fuse['nor'].params[''] = PGMparameter(1, 1)
        fuse['scl'].params[''] = PGMparameter(1, 1)
        fuse['flt'].params[''] = PGMparameter(2, 4, pos_u=[0.28,0.53])
        fuse['pos'].params['nose'] = PGMparameter(3, 3, pos_u=[0,0.065,0.13], order_u=3)
        fuse['scl'].params['nose'] = PGMparameter(3, 1, pos_u=[0,0.07,0.14], order_u=3)
        fuse['scl'].params['tail'] = PGMparameter(2, 1, pos_u=[0.7,1.0])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''] = PGMparameter(1, 3)
        lwing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.35,1.0])
        lwing['pos'].params['lin'] = PGMparameter(2, 3)
	lwing['shY','upp'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)
	lwing['shY','low'].params[''] = PGMparameter(10, 6, order_u=4, order_v=4)

        lpylon = self.comps['lpylon'].props
        lpylon['pos'].params[''] = PGMparameter(1, 3)
        lpylon['pos'].params['lin'] = PGMparameter(2, 3)
        lpylon['scl'].params[''] = PGMparameter(2, 1)
        lpylon['nor'].params[''] = PGMparameter(1, 3)

        lnac = self.comps['lnac'].props
        lnac['pos'].params[''] = PGMparameter(1, 3)
        lnac['pos'].params['lin'] = PGMparameter(2, 3)
        lnac['nor'].params[''] = PGMparameter(1, 1)
        lnac['scl'].params[''] = PGMparameter(1, 1)
        lnac['thk'].params[''] = PGMparameter(3, 1)

        ltail = self.comps['ltail'].props
        ltail['pos'].params[''] = PGMparameter(1, 3)
        ltail['pos'].params['lin'] = PGMparameter(2, 3)
        ltail['scl'].params[''] = PGMparameter(2, 1)
        ltail['rot'].params[''] = PGMparameter(2, 3)
        ltail['ogn'].params[''] = PGMparameter(1, 3)

        vtail = self.comps['vtail'].props
        vtail['pos'].params[''] = PGMparameter(1, 3)
        vtail['pos'].params['lin'] = PGMparameter(2, 3)
        vtail['nor'].params[''] = PGMparameter(1, 3)
        vtail['scl'].params[''] = PGMparameter(2, 1)
        vtail['rot'].params[''] = PGMparameter(2, 3)
        vtail['ogn'].params[''] = PGMparameter(1, 3)

    def _define_dvs(self):
        dvs = self.dvs
        dvs['span'] = PGMdv((1), 23.3).set_identity_param('lwing', 'pos', 'lin', (1,2))
        dvs['mid_chord'] = PGMdv((1), 4.5).set_identity_param('lwing', 'scl', '', (1,0))
        dvs['tip_chord'] = PGMdv((1), 1.2).set_identity_param('lwing', 'scl', '', (2,0))
	dvs['shape_wing_upp'] = PGMdv((10,6)).set_identity_param('lwing', ('shY', 'upp'), '')
	dvs['shape_wing_low'] = PGMdv((10,6)).set_identity_param('lwing', ('shY', 'low'), '')

    def _compute_params(self):
        fuse = self.comps['fuse'].props
        fuse['pos'].params[''].val([[0,0,0],[50,0,0]])
        fuse['nor'].params[''].val([1.0])
        fuse['scl'].params[''].val([2.6])
        fuse['flt'].params[''].val([[0,0,0.5,0.5],[0,0,0.5,0.5]])
        fuse['pos'].params['nose'].val([[0,-1.1,0],[0,0,0],[0,0,0]])
        fuse['scl'].params['nose'].val([-2, 0, 0])
        fuse['scl'].params['tail'].val([0, -2.3])

        lwing = self.comps['lwing'].props
        lwing['pos'].params[''].val([16,-1,2.6])
        lwing['scl'].params[''].val([10,4.5,1.2])
        lwing['pos'].params['lin'].val([[0,0,0],[16.5,4.4,23.3]])

        lpylon = self.comps['lpylon'].props
        lpylon['pos'].params[''].val([21.2,-0.5,9])
        lpylon['pos'].params['lin'].val([[0,0,0],[-2,-0.5,0]])
        lpylon['scl'].params[''].val([2.1,2.5])
        lpylon['nor'].params[''].val([1,0,0])

        lnac = self.comps['lnac'].props
        lnac['pos'].params[''].val([16.4,-2.4,9])
        lnac['pos'].params['lin'].val([[0,0,0],[4.5,0,0]])
        lnac['nor'].params[''].val([1])
        lnac['scl'].params[''].val([1.25])
        lnac['thk'].params[''].val([0.08,0.2,0.08])

        ltail = self.comps['ltail'].props
        ltail['pos'].params[''].val([44,0,1.3])
        ltail['pos'].params['lin'].val([[0,0,0],[6,1.4,8]])
        ltail['scl'].params[''].val([4,1])
        ltail['rot'].params[''].val([[0,10,0],[0,0,0]])
        ltail['ogn'].params[''].val([0.25,0,0])

        vtail = self.comps['vtail'].props
        vtail['pos'].params[''].val([42,1.7,0])
        vtail['pos'].params['lin'].val([[0,0,0],[6,8,0]])
        vtail['nor'].params[''].val([1,0,0])
        vtail['scl'].params[''].val([5.8,2])
        vtail['rot'].params[''].val([[0,10,0],[0,0,0]])
        vtail['ogn'].params[''].val([0.25,0,0])

        return [], [], []

    def _set_bspline_options(self):
        comps = self.comps

        comps['fuse'].faces['rgt'].set_option('num_cp', 'u', [4,4,4,4])
        comps['fuse'].faces['rgt'].set_option('num_cp', 'v', [18,4,4,4,4,8,4,15,4,4,10,4])
        comps['fuse'].faces['rgt'].set_option('num_pt', 'v', [40,16,16,16,16,60,16,60,16,16,70,16], both=False)
        comps['fuse'].faces['top'].set_option('num_cp', 'u', [8,8])
        comps['lwing'].faces['upp'].set_option('num_cp', 'v', [6,4,4,20])
        
        

if __name__ == '__main__':

    pgm = Conventional()
    bse = pgm.initialize()

    pgm.comps['lwing'].set_airfoil('rae2822.dat')
    pgm.comps['ltail'].set_airfoil()
    pgm.dvs['shape_wing_upp'].data[2,2] = 0.0
    pgm.dvs['span'].data[0] = 23.3
    pgm.dvs['tip_chord'].data[0] = 1.2
    pgm.compute_all()

    bse.vec['pt_str'].export_tec_str()
    bse.vec['df'].export_tec_scatter()
    bse.vec['cp'].export_tec_scatter()
    bse.vec['pt'].export_tec_scatter()
    bse.vec['cp_str'].export_IGES()
