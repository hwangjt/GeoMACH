from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys

from GeoMACH.PGM.components import Wing, Body, Shell, Junction, Cone, Tip
from GeoMACH.PGM.configurations import Configuration
from GeoMACH.PSM import Airframe


class Trussbraced(Configuration):

    def define_primitive_comps(self):
        comps = {
            'fu': Body(nx=17, ny=6, nz=4),
            'lw': Wing(nx=7, nz=7, right=0),
            'rw': Wing(nx=7, nz=7, left=0),
            'ls': Wing(nx=4, nz=4, left=0, right=0),
            'lv': Wing(left=0, right=0),
            'vt': Wing(nx=4,nz=4,right=0),
            'lt': Wing(right=0),
        }
        return comps

    def define_interpolant_comps(self):
        c = self.comps
        comps = {
            'fu_n': Cone(c['fu'], 0),
            'fu_t': Cone(c['fu'], 1),
            'lw_T': Tip(c['lw'], 1),
            'rw_T': Tip(c['rw'], 0),
            'lt_T': Tip(c['lt'], 1),
            'vt_T': Tip(c['vt'], 1),
            'lw_fu': Junction(c['fu'], 'lft', 0, [0,1], c['lw'], mSide=0),
            'rw_fu': Junction(c['fu'], 'rgt', 2, [0,7], c['rw'], mSide=1),
            'ls_fu': Junction(c['fu'], 'lft', 0, [4,2], c['ls'], mSide=0),
            'ls_lw': Junction(c['lw'], 'low', 3, [4,1], c['ls'], mSide=1),
            'lv_lw': Junction(c['lw'], 'low', 3, [1,3], c['lv'], mSide=1),
            'lv_ls': Junction(c['ls'], 'upp', 3, [1,0], c['lv'], mSide=0),
            'vt_fu': Junction(c['fu'], 'top', 0, [1,11], c['vt'], mSide=0),
            'lt_vt': Junction(c['vt'], 'low', 1, [0,0], c['lt'], mSide=0),
            #'lp_ln': Junction(c['ln'], 'tp0', 2, [1,0], c['lp'], mSide=1),
            #'rp_rw': Junction(c['rw'], 'low', 1, [1,0], c['rp'], mSide=0),
            #'rp_rn': Junction(c['rn'], 'tp0', 2, [1,0], c['rp'], mSide=1),
            #'lw_fu': Junction(c['fu'], 'lft', 0, [2,1], c['lw'], mSide=0),
            #'rw_fu': Junction(c['fu'], 'rgt', 2, [2,5], c['rw'], mSide=1),
            #'lt_fu': Junction(c['fu'], 'lft', 0, [1,9], c['lt'], mSide=0),
            #'rt_fu': Junction(c['fu'], 'rgt', 2, [1,0], c['rt'], mSide=1),
                 }
        return comps

    def define_dvs(self):
        self.add_dv('dv1', [1], 19)

    def apply_dvs(self):
        Das = [numpy.zeros(1)]
        Dis = [numpy.zeros(1)]
        Djs = [numpy.zeros(1)]

        return Das, Dis, Djs

    def set_oml_resolution(self):
        comps = self.comps
        comps['lv'].faces['upp'].setm(0,[4])
        comps['lv'].faces['upp'].setm(1,[4])
        comps['lv'].faces['low'].setm(0,[4])
        comps['lv'].faces['low'].setm(1,[4])
        comps['lt'].faces['upp'].setm(0,[4])
        comps['lt'].faces['upp'].setm(1,[4])
        comps['lt'].faces['low'].setm(0,[4])
        comps['lt'].faces['low'].setm(1,[4])

        comps['fu'].faces['rgt'].setm(1,[42,4,4,4,7,4,4,4,4,4,33,4,4,12,4,4,4])
        comps['fu'].faces['lft'].setm(0,[4,4,8,8,4,4])
        comps['lw'].faces['upp'].setm(1,[8,4,4,5,4,4,22])
        comps['rw'].faces['upp'].setm(1,[8,4,4,5,4,4,22][::-1])
#        comps['rw'].faces['low'].setm(0,[4,4,4,7,4,4,4])
        comps['ls'].faces['upp'].setm(1,[10,4,4,10])

        comps['lw_T'].faces['def'].setm(0,[10])
        comps['lt_T'].faces['def'].setm(0,[10])
        comps['vt_T'].faces['def'].setm(0,[10])

        #comps['fu'].faces['rgt'].setm(1,[18,4,4,4,4,8,4,15,4,4,10,4])
        #comps['fu'].faces['rgt'].setm(0,[4,4,4,4])
        #comps['fu'].faces['top'].setm(0,[8,8])
        #comps['lw'].faces['upp'].setm(1,[6,4,4,20])

    def define_oml_parameters(self):
        c = self.comps

        c['lw'].setAirfoil('rae2822.dat')
        c['rw'].setAirfoil('rae2822.dat')

        c['fu'].props['nor'].addParam('nor1',[1,3],P=[1.0,0.0,1.0])
        c['fu'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[38,0,0]])
        c['fu'].props['pos'].addParam('nose',[2,3],P=[[0,-1.1,0],[0,0,0]],T=[0,0.13],B=[False,True])
        c['fu'].props['pos'].addParam('tail',[2,3],P=[[0,0,0],[0,1.6,0]],T=[0.75,1.0],B=[False,False])
        c['fu'].props['scl'].addParam('rad',[4,1],P=[0.65,2,2,0.4],T=[0,0.14,0.75,1.0],B=[False,True,False,False])
        c['fu'].props['scl'].addParam('tail',[2,3],P=[[0,0,0],[-0.3,0,0]],T=[0.75,1.0])
        c['fu'].props['flt'].addParam('flt1',[2,4],P=[[0,0,0.9,0.9],[0,0,0.9,0.9]],T=[0.33,0.50])
        c['fu'].props['flt'].addParam('flt2',[2,4],P=[[0.9,0.9,0,0],[0.9,0.9,0,0]],T=[0.32,0.55])

        c['lw'].props['pos'].addParam('offset',[1,3],P=[14,1.5,2.1])
        c['lw'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[3.8,1.5,30]])
        c['lw'].props['scl'].addParam('scl1',[2,1],P=[4.9,0.9])

        c['rw'].props['pos'].addParam('offset',[1,3],P=[14,1.5,-2.1])
        c['rw'].props['pos'].addParam('pos1',[2,3],P=[[3.8,1.5,-30],[0,0,0]])
        c['rw'].props['scl'].addParam('scl1',[2,1],P=[0.9,4.9])

        c['ls'].props['pos'].addParam('offset',[1,3],P=[15,-1.6,2.3])
        c['ls'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[1,3.0,30*0.45]])
        c['ls'].props['pos'].addParam('curv',[2,3],P=[[0,0,0],[0,0.6,0]],B=[True,False],T=[0.8,1.0])
        c['ls'].props['scl'].addParam('scl1',[2,1],P=[2.5,1.9])
        c['ls'].props['nor'].addParam('nor',[1,1],P=[1.0])
        c['ls'].props['rot'].addParam('rot',[2,3],P=[[0,0,0],[-60,0,0]],T=[0.7,1])

        c['lv'].props['pos'].addParam('offset',[1,3],P=[15.5,0,8.8])
        c['lv'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[0,1.58,0]])
        c['lv'].props['scl'].addParam('scl1',[2,1],P=[1,1.4])
        c['lv'].props['nor'].addParam('nor',[1,1],P=[1.0])

        c['vt'].props['pos'].addParam('offset',[1,3],P=[31.3,2.1,0])
        c['vt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[3,6,0]])
        c['vt'].props['scl'].addParam('scl1',[2,1],P=[5.8,3.3])
        c['vt'].props['nor'].addParam('nor',[1,3],P=[1.0,0.0,0.0])

        c['lt'].props['pos'].addParam('offset',[1,3],P=[33.9,6.7,0.25])
        c['lt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[3.3,0,5]])
        c['lt'].props['scl'].addParam('scl1',[2,1],P=[2.9,1])

        c['fu_n'].props['fC1'].params['fC1'].setP([4])
        c['fu_t'].props['fC1'].params['fC1'].setP([4])
        c['lw_T'].props['fC1'].params['fC1'].setP([0.5])
        c['vt_T'].props['fC1'].params['fC1'].setP([1.0])
        c['lt_T'].props['fC1'].params['fC1'].setP([0.5])
        c['lw_fu'].props['mC1'].params['mC1'].setP([0.1])
        c['lw_fu'].props['fC1'].params['fC1'].setP([1.0])
        c['ls_fu'].props['mC1'].params['mC1'].setP([0.02])
        c['ls_fu'].props['fC1'].params['fC1'].setP([0.02])
        c['ls_lw'].props['mC1'].params['mC1'].setP([0.1])
        c['ls_lw'].props['fC1'].params['fC1'].setP([0.1])
        c['lv_ls'].props['mC1'].params['mC1'].setP([0.1])
        c['lv_ls'].props['fC1'].params['fC1'].setP([0.1])
        c['lv_lw'].props['mC1'].params['mC1'].setP([0.1])
        c['lv_lw'].props['fC1'].params['fC1'].setP([0.1])
        c['vt_fu'].props['mC1'].params['mC1'].setP([0.1])
        c['vt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['lt_vt'].props['mC1'].params['mC1'].setP([0.1])
        c['lt_vt'].props['fC1'].params['fC1'].setP([0.1])

    def meshStructure(self):
        afm = Airframe(self, 1) #0.2)

        idims = numpy.linspace(0.3,0.85,7)
        jdims = numpy.linspace(0,1,17)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlw_1:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                if i is 0 or i is idims.shape[0]-1:
                    afm.addVertFlip('Mlw_2:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                else:
                    afm.addVertFlip('Mlw_2a:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mlw_2b:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[0.15,0])

        idims = numpy.linspace(0.3,0.85,7)
        jdims = numpy.linspace(0,0.9,17)
        for i in range(idims.shape[0]):
            if i is 0 or i is idims.shape[0]-1:
                afm.addCtrVert('Mcw_2:'+str(i)+':'+str(j),'lw','rw',idims[i])
            else:
                afm.addCtrVert('Mcw_2a:'+str(i)+':'+str(j),'lw','rw',idims[i],w=[1,0.85])
                afm.addCtrVert('Mcw_2b:'+str(i)+':'+str(j),'lw','rw',idims[i],w=[0.15,0])
        for i in range(idims.shape[0]-1):
            afm.addCtr('Mcw_u:','lw','rw',0,[idims[i],idims[i+1]])
        for i in range(idims.shape[0]-1):
            afm.addCtr('Mcw_l:','lw','rw',1,[1-idims[i],1-idims[i+1]])
        afm.addCtrVert('Mcw_sec:'+str(i)+':'+str(j),'lw','rw',0.18)

        idims = numpy.linspace(0.25,0.65,2)
        jdims = numpy.linspace(0,1,11)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlt_1:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mvt_1:'+str(i)+':'+str(j),'vt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mls_1:'+str(i)+':'+str(j),'ls',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mlv_1:'+str(i)+':'+str(j),'lv',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVertFlip('Mlt_2:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mvt_2:'+str(i)+':'+str(j),'vt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mls_2:'+str(i)+':'+str(j),'ls',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mlv_2:'+str(i)+':'+str(j),'lv',[idims[i],jdims[j]],[idims[i],jdims[j+1]])

        idims = numpy.linspace(0,1,5)
        jdims = numpy.linspace(0,1,21)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVert('Mfu_F1:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[0,2])
                afm.addVert('Mfu_F2:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[1,3])
                afm.addVert('Mfu_F3:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[2,0])
                afm.addVert('Mfu_F4:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[3,1])
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]-1):
                afm.addVert('Mfu_L1:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[0,2])
                afm.addVert('Mfu_L2:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[1,3])
                afm.addVert('Mfu_L3:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[2,0])
                afm.addVert('Mfu_L4:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[3,1])
        for j in range(jdims.shape[0]-1):
            afm.addVertFlip('Mfu_0a:'+str(j),'fu',[0.4,jdims[j]],[0.4,jdims[j+1]],w=[1.0,0.5],i=[0,2])
            afm.addVertFlip('Mfu_0b:'+str(j),'fu',[0.4,jdims[j]],[0.4,jdims[j+1]],w=[0.5,0.0],i=[0,2])

        afm.preview('trussbraced_pvw.dat')
        afm.mesh()
        afm.computeMesh('trussbraced_str.dat')


if __name__ == '__main__':

    import cProfile

    name = 'trussbraced'
    aircraft = Trussbraced()

    #der = aircraft.get_derivatives('lw', 'scl1', (0,0), useFD=False)
    #aircraft.oml0.addVars(['der'])
    ##aircraft.oml0.addVars(['dx','dy','dz'])
    ##aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)
    #aircraft.oml0.Q[:,6] = numpy.sum(der*der,1)**0.5
    aircraft.oml0.computePoints()
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)
    aircraft.oml0.write2IGES(name)

    aircraft.comps['lw'].add_thk_con('name', 
                                     numpy.linspace(0.1,0.9,10),
                                     numpy.linspace(0.1,0.9,10),
                                     0.25)

    aircraft.meshStructure()

    #aircraft.test_derivatives()
    #aircraft.oml0.plot()
    #aircraft.meshStructure()

#    cProfile.run('aircraft.meshStructure()')
