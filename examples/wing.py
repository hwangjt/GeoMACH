# John Hwang, 2014

from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMtip
from GeoMACH.PSM import Airframe
import numpy


class Wing(PGMconfiguration):

    def _define_comps(self):
        self.comps['wing'] = PGMwing(num_x=1, num_z=1, left_closed=True)
        self.comps['tip'] = PGMtip(self, 'wing', 'left', 0.1)

    def _define_params(self):
        wing = self.comps['wing'].props
        wing['pos'].params[''] = PGMparameter(3, 3, pos_u=[0,0.37,1.0])
        wing['scl'].params[''] = PGMparameter(3, 1, pos_u=[0,0.37,1.0])

    def _compute_params(self):
        wing = self.comps['wing'].props
        wing['pos'].params[''].data[0, :] = [904.294, 174.126, 0.0]
        wing['pos'].params[''].data[1, :] = [1225.82, 181.071, 427.999]
        wing['pos'].params[''].data[2, :] = [1780.737, 263.827, 1156.753]
        wing['scl'].params[''].data[:, 0] = [536.181, 285.782, 107.4]
        return [], [], []

    def _set_bspline_options(self):
        wing = self.comps['wing'].faces
        wing['upp'].set_option('num_cp', 'u', [40])
        wing['upp'].set_option('num_cp', 'v', [40])

    def meshStructure(self, res, filename):
        afm = Airframe(self, res)
        idims = numpy.linspace(0.45,0.9,7)
        jdims = numpy.linspace(0,1,16)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlw_1:'+str(i)+':'+str(j),'wing',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                if i is 0 or i is idims.shape[0]-1:
                    afm.addVertFlip('Mlw_2:'+str(i)+':'+str(j),'wing',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                else:
                    afm.addVertFlip('Mlw_2a:'+str(i)+':'+str(j),'wing',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mlw_2b:'+str(i)+':'+str(j),'wing',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[0.15,0])
        idims = numpy.linspace(0.18,0.45,6)
        for j in range(idims.shape[0]-1):
            afm.addVertFlip('Mlw_sec1:'+str(j),'wing',[idims[j],jdims[j]],[idims[j+1],jdims[j+1]])
            afm.addVertFlip('Mlw_sec2:'+str(j),'wing',[idims[j],jdims[j]],[0.45,jdims[j]])
        afm.preview('wing_pvw.dat')
        afm.mesh()
        afm.computeMesh(filename + '_str.dat')

        

if __name__ == '__main__':

    pgm = Wing()
    bse = pgm.initialize()
    pgm.comps['wing'].set_airfoil('rae2822.dat')
    bse.vec['pt_str'].export_tec_str()
    bse.vec['cp_str'].export_IGES()
    exit()
    for num in [20]:#[5, 10, 20, 50, 100]:
        pgm.meshStructure(num, 'wing_'+str(num))
