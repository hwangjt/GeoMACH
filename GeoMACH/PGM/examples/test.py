from __future__ import division

from GeoMACH.PGM.core.PGMconfiguration import PGMconfiguration
from GeoMACH.PGM.core.PGMparameter import PGMparameter
from GeoMACH.PGM.core.PGMdv import PGMdv
from GeoMACH.PGM.components.PGMwing import PGMwing
from GeoMACH.PGM.components.PGMbody import PGMbody
from GeoMACH.PGM.components.PGMshell import PGMshell
from GeoMACH.PGM.components.PGMjunction import PGMjunction
from GeoMACH.PGM.components.PGMtip import PGMtip
from GeoMACH.PGM.components.PGMcone import PGMcone



class Test(PGMconfiguration):

    def _define_comps(self):
        self.comps['wing'] = PGMwing(left_closed=True)
        self.comps['body'] = PGMbody(5,4,2)
        self.comps['shell'] = PGMshell(2,2,2)
        self.comps['junction'] = PGMjunction(self, 'body', 'lft', 'E',
                                             [1, 1], 'wing', 'right')
        self.comps['tip'] = PGMtip(self, 'wing', 'left')
        self.comps['cone1'] = PGMcone(self, 'body', 'front', 10)
        self.comps['cone2'] = PGMcone(self, 'body', 'rear', 10)

    def _define_params(self):
        wing = self.comps['wing'].props
        wing['pos'].params[''] = PGMparameter(2, 3)
        wing['scl'].params[''] = PGMparameter(1, 1)

        body = self.comps['body'].props
        body['pos'].params[''] = PGMparameter(2, 3)
        body['scl'].params[''] = PGMparameter(1, 3)
        body['nor'].params[''] = PGMparameter(1, 1)    

        shell = self.comps['shell'].props
        shell['pos'].params[''] = PGMparameter(2, 3)
        shell['scl'].params[''] = PGMparameter(1, 1)
        shell['nor'].params[''] = PGMparameter(1, 1)    
        shell['thk'].params[''] = PGMparameter(1, 1)

    def _compute_params(self):
        wing = self.comps['wing'].props
        wing['pos'].params[''].val([[-0.5,0,1],[1,0,4]])
        wing['scl'].params[''].val([1])

        body = self.comps['body'].props
        body['pos'].params[''].val([[-3,0,0],[3,0,0]])
        body['scl'].params[''].val([1, 0.4, 0])
        body['nor'].params[''].val([1])

        shell = self.comps['shell'].props
        shell['pos'].params[''].val([[6,0,0],[7,0,0]])
        shell['scl'].params[''].val([1])
        shell['nor'].params[''].val([1])
        shell['thk'].params[''].val([0.1])
        return [], [], []

    def _set_bspline_options(self):
        wing = self.comps['wing'].faces
        wing['low'].set_option('num_cp', 'u', [10])
        wing['upp'].set_option('num_cp', 'u', [10])
        wing['upp'].set_option('num_cp', 'v', [10])

        body = self.comps['body'].faces
        body['rgt'].set_option('num_cp', 'v', [10])
        body['rgt'].set_option('num_cp', 'u', [10])
        body['top'].set_option('num_cp', 'u', [10])

        shell = self.comps['shell'].faces
        shell['rt0'].set_option('num_cp', 'v', [10])
        shell['rt0'].set_option('num_cp', 'u', [10])
        shell['tp0'].set_option('num_cp', 'u', [10])
        shell['rt1'].set_option('num_cp', 'v', [10])
        shell['rt1'].set_option('num_cp', 'u', [10])
        shell['tp1'].set_option('num_cp', 'u', [10])
        

if __name__ == '__main__':

    t = Test()
    bse = t.initialize()
    bse.vec['pt_str'].export_tec_str()
    bse.vec['cp'].export_tec_scatter()
    bse.print_info()
