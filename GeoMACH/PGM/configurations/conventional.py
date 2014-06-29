from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys

from GeoMACH.PGM.components import Wing, Body, Shell, Junction, Cone, Tip
from GeoMACH.PGM.configurations import Configuration
from GeoMACH.PSM import Airframe


class Conventional(Configuration):

    def define_primitive_comps(self):
        comps = {'fu': Body(nx=12, ny=4, nz=2),
                 'lw': Wing(nx=4, nz=4, right=0),
                 'rw': Wing(nx=4, nz=4, left=0),
                 'lp': Wing(left=0, right=0),
                 'rp': Wing(left=0, right=0),
                 'ln': Shell(nx=4, ny=1, nz=4),
                 'rn': Shell(nx=4, ny=1, nz=4),
                 'lt': Wing(right=0),
                 'rt': Wing(left=0),
                 'vt': Wing(nx=2,right=0),
                 }
        return comps

    def define_interpolant_comps(self):
        c = self.comps
        comps = {'fu_n': Cone(c['fu'], 0),
                 'fu_t': Cone(c['fu'], 1),
                 'lw_T': Tip(c['lw'], 1),
                 'rw_T': Tip(c['rw'], 0),
                 'lt_T': Tip(c['lt'], 1),
                 'rt_T': Tip(c['rt'], 0),
                 'vt_T': Tip(c['vt'], 1),
                 'lp_lw': Junction(c['lw'], 'low', 1, [1,0], c['lp'], mSide=0),
                 'lp_ln': Junction(c['ln'], 'tp0', 2, [1,0], c['lp'], mSide=1),
                 'rp_rw': Junction(c['rw'], 'low', 1, [1,0], c['rp'], mSide=0),
                 'rp_rn': Junction(c['rn'], 'tp0', 2, [1,0], c['rp'], mSide=1),
                 'lw_fu': Junction(c['fu'], 'lft', 0, [2,1], c['lw'], mSide=0),
                 'rw_fu': Junction(c['fu'], 'rgt', 2, [2,5], c['rw'], mSide=1),
                 'lt_fu': Junction(c['fu'], 'lft', 0, [1,9], c['lt'], mSide=0),
                 'rt_fu': Junction(c['fu'], 'rgt', 2, [1,0], c['rt'], mSide=1),
                 'vt_fu': Junction(c['fu'], 'top', 0, [0,8], c['vt'], mSide=0),
                 }
        return comps

    def define_dvs(self):
        self.add_dv('sweep', [1], val=16.5, lower=15, upper=23)
        self.add_dv('rot', [1], val=0.1, lower=0.0, upper=1.0)

    def apply_dvs(self):
        pos1 = self.comps['lw'].props['pos'].params['pos1']
        rot1 = self.comps['lt'].props['rot'].params['rot2']
        sweep = self.dvs['sweep']
        rot = self.dvs['rot']

        Das, Dis, Djs = [], [], []

        pos1.param_vec[1,0,0] = sweep.vec[0]
        Das.append(1.0)
        Dis.append(pos1.param_ind[1,0,0])
        Djs.append(sweep.ind[0])

        rot1.param_vec[1,2] = rot.vec[0]
        Das.append(1.0)
        Dis.append(rot1.param_ind[1,2])
        Djs.append(rot.ind[0])

        return Das, Dis, Djs

    def set_oml_resolution(self):
        comps = self.comps
        comps['fu'].faces['rgt'].setm(1,[18,4,4,4,4,8,4,15,4,4,10,4])
        comps['fu'].faces['rgt'].setm(0,[4,4,4,4])
        comps['fu'].faces['top'].setm(0,[8,8])
        comps['lw'].faces['upp'].setm(1,[6,4,4,20])
        comps['rw'].faces['upp'].setm(1,[20,4,4,6])
        comps['lt'].faces['upp'].setm(1,[15])
        comps['rt'].faces['upp'].setm(1,[15])
        comps['vt'].faces['upp'].setm(1,[15])
        comps['ln'].faces['rt0'].setm(1,[4])
        comps['ln'].faces['rt1'].setm(1,[4])
        comps['rn'].faces['rt0'].setm(1,[4])
        comps['rn'].faces['rt1'].setm(1,[4])
        comps['ln'].faces['rt0'].setm(0,[4])
        comps['ln'].faces['lt0'].setm(0,[4])
        comps['rn'].faces['rt0'].setm(0,[4])
        comps['rn'].faces['lt0'].setm(0,[4])
        comps['lw_T'].faces['def'].setm(0,[10])
        comps['rw_T'].faces['def'].setm(0,[10])
        comps['lt_T'].faces['def'].setm(0,[10])
        comps['rt_T'].faces['def'].setm(0,[10])
        comps['vt_T'].faces['def'].setm(0,[10])

    def define_oml_parameters(self):
        c = self.comps

        c['lw'].setAirfoil('rae2822.dat')
        c['rw'].setAirfoil('rae2822.dat')

        c['fu'].props['nor'].addParam('nor1',[1,1],P=[1.0])
        c['fu'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[50,0,0]])
        c['fu'].props['pos'].addParam('nose',[2,3],P=[[0,-1.1,0],[0,0,0]],T=[0,0.13],B=[False,True])
        c['fu'].props['scl'].addParam('rad',[4,1],P=[0.65,2.6,2.6,0.35],T=[0,0.14,0.7,1.0],B=[False,True,True,False])
        c['fu'].props['flt'].addParam('flt1',[2,4],P=[[0,0,0.5,0.5],[0,0,0.5,0.5]],T=[0.28,0.53])

        #c['lw'].props['pos'].addParam('offset',[1,3],P=[16,-1,3])
        c['lw'].props['pos'].addParam('offset',[1,3],P=[16,-1,2.6])
        c['lw'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[16.5,4.4,23.3]],T=[0,1.0])
        c['lw'].props['scl'].addParam('scl1',[3,1],P=[10,4.5,1.2],T=[0,0.35,1.0])
#        c['lw'].props['pos'].addParam('rkd_p',[2,3],P=[[0,0,0],[3,2,0]],T=[0.85,1.0],B=[[True,True,False],[False,False,False]])
#        c['lw'].props['scl'].addParam('rkd_s',[2,1],P=[0,-1.5],T=[0.85,1.0],B=[True,False])
        #c['lw'].addParam('aileron1','shU',[2,4],P=[[0,0.16,0.16,0],[0,0,0,0]])
        #c['lw'].params['aileron1'].setT([0.0,0.5],0)
        #c['lw'].params['aileron1'].setT([0.6,0.63,0.8,0.83],1)
        #c['lw'].params['aileron1'].setB([[0,0,0,0],[1,1,1,1]],0)
        #c['lw'].params['aileron1'].setB([[1,1,1,1],[1,1,1,1]],1)
        #c['lw'].addParam('aileron2','shL',[2,4],P=[[0,0,0,0],[0,-0.16,-0.16,0]])
        #c['lw'].params['aileron2'].setT([0.5,1.0],0)
        #c['lw'].params['aileron2'].setT([0.6,0.63,0.8,0.83],1)
        #c['lw'].params['aileron2'].setB([[1,1,1,1],[0,0,0,0]],0)
        #c['lw'].params['aileron2'].setB([[1,1,1,1],[1,1,1,1]],1)

        c['rw'].props['pos'].addParam('offset',[1,3],P=[16,-1,-3])
        c['rw'].props['pos'].addParam('pos1',[2,3],P=[[19,3,-24],[0,0,0]],T=[0,1.0])
        c['rw'].props['scl'].addParam('scl1',[3,1],P=[1.8,4.5,10],T=[0,0.65,1.0])
        c['rw'].props['pos'].addParam('rkd_p',[2,3],P=[[3,2,0],[0,0,0]],T=[0,0.15],B=[[False,False,False],[True,True,False]])
        c['rw'].props['scl'].addParam('rkd_s',[2,1],P=[-1.5,0],T=[0,0.15],B=[False,True])
        #c['rw'].addParam('aileron1','shU',[2,4],P=[[0,-0.16,-0.16,0],[0,0,0,0]])
        #c['rw'].params['aileron1'].setT([0.0,0.5],0)
        #c['rw'].params['aileron1'].setT([0.17,0.2,0.37,0.4],1)
        #c['rw'].params['aileron1'].setB([[0,0,0,0],[1,1,1,1]],0)
        #c['rw'].params['aileron1'].setB([[1,1,1,1],[1,1,1,1]],1)
        #c['rw'].addParam('aileron2','shL',[2,4],P=[[0,0,0,0],[0,0.16,0.16,0]])
        #c['rw'].params['aileron2'].setT([0.5,1.0],0)
        #c['rw'].params['aileron2'].setT([0.17,0.2,0.37,0.4],1)
        #c['rw'].params['aileron2'].setB([[1,1,1,1],[0,0,0,0]],0)
        #c['rw'].params['aileron2'].setB([[1,1,1,1],[1,1,1,1]],1)

        c['lt'].props['ogn'].addParam('ogn',[1,3],P=[0.25,0,0])
        c['lt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[6,1.4,8]])
        c['lt'].props['pos'].addParam('offset',[1,3],P=[44,0,1.7])
        c['lt'].props['scl'].addParam('scl1',[2,1],P=[4,1])
        c['lt'].props['rot'].addParam('rot1',[2,3],P=[[0,10,0],[0,0,0]])

        c['lt'].props['rot'].addParam2('rot2', [6,3])
        c['lt'].props['rot'].params['rot2'].param0[1,2] = 15
#        c['lt'].props['shY','upp'].addParam2('shp', [10,10])
#        c['lt'].props['shY','upp'].params['shp'].param0[5,5] = 1.0

        c['rt'].props['ogn'].addParam('ogn',[1,3],P=[0.25,0,0])
        c['rt'].props['pos'].addParam('pos1',[2,3],P=[[6,0,-8],[0,0,0]])
        c['rt'].props['pos'].addParam('offset',[1,3],P=[44,0,-1.7])
        c['rt'].props['scl'].addParam('scl1',[2,1],P=[1,4])
        c['rt'].props['rot'].addParam('rot1',[2,3],P=[[0,0,0],[0,-10,0]])

        c['vt'].props['ogn'].addParam('ogn',[1,3],P=[0.25,0,0])
        c['vt'].props['nor'].addParam('nor1',[1,3],P=[1,0,0])
        c['vt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[6,8,0]])
        c['vt'].props['pos'].addParam('offset',[1,3],P=[42,2.0,0])
        c['vt'].props['scl'].addParam('scl1',[2,1],P=[5.8,2])
        c['vt'].props['rot'].addParam('rot1',[2,3],P=[[0,10,0],[0,0,0]])
        #c['vt'].addParam('aileron1','shU',[2,1],P=[-0.12,0.0])
        #c['vt'].params['aileron1'].setT([0.0,0.25],0)
        #c['vt'].addParam('aileron2','shL',[2,1],P=[0.0,0.12])
        #c['vt'].params['aileron2'].setT([0.75,1.0],0)

        c['lp'].props['ogn'].addParam('ogn',[1,3],P=[0.25,0,0])
        c['lp'].props['nor'].addParam('nor1',[1,3],P=[1,0,0])
        c['lp'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[-2,-0.5,0]])
        c['lp'].props['pos'].addParam('offset',[1,3],P=[21.7,-0.5,9])
        c['lp'].props['scl'].addParam('scl1',[2,1],P=[2.1,2.5])

        c['ln'].props['nor'].addParam('nor1',[1,1],P=[1])
        c['ln'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[4.5,0,0]])
        c['ln'].props['pos'].addParam('offset',[1,3],P=[16.4,-2.4,9])
        c['ln'].props['scl'].addParam('scl1',[1,1],P=[1.25])
        c['ln'].props['thk'].addParam('thk1',[3,1],P=[0.08,0.2,0.08],B=[False,True,False])

        c['rp'].props['ogn'].addParam('ogn',[1,3],P=[0.25,0,0])
        c['rp'].props['nor'].addParam('nor1',[1,3],P=[1,0,0])
        c['rp'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[-2,-0.3,0]])
        c['rp'].props['pos'].addParam('offset',[1,3],P=[21.2,-0.7,-9])
        c['rp'].props['scl'].addParam('scl1',[2,1],P=[2.1,3])

        c['rn'].props['nor'].addParam('nor1',[1,1],P=[1])
        c['rn'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[4.5,0,0]])
        c['rn'].props['scl'].addParam('scl1',[1,1],P=[1.25])
        c['rn'].props['pos'].addParam('offset',[1,3],P=[16,-2.4,-9])
        c['rn'].props['thk'].addParam('thk1',[3,1],P=[0.08,0.2,0.08],B=[False,True,False])

        #c['lw_fu'].props['mC1'].params['mC1'].setP([1.5])
        c['lw_fu'].props['mC1'].params['mC1'].setP([0.01])
        c['lw_fu'].props['fC1'].params['fC1'].setP([0.5])
        c['rw_fu'].props['mC1'].params['mC1'].setP([1.5])
        c['rw_fu'].props['fC1'].params['fC1'].setP([0.5])
        c['lt_fu'].props['mC1'].params['mC1'].setP([0.1])
        c['lt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['rt_fu'].props['mC1'].params['mC1'].setP([0.1])
        c['rt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['vt_fu'].props['mC1'].params['mC1'].setP([0.01])
        c['vt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['lp_lw'].props['mC1'].params['mC1'].setP([0])
        c['lp_ln'].props['mC1'].params['mC1'].setP([0])
        c['rp_rn'].props['mC1'].params['mC1'].setP([0])
        c['fu_n'].props['fC1'].params['fC1'].setP([2])
        c['fu_t'].props['fC1'].params['fC1'].setP([2])
        c['lw_T'].props['fC1'].params['fC1'].setP([0.5])
        c['rw_T'].props['fC1'].params['fC1'].setP([0.5])
        c['lt_T'].props['fC1'].params['fC1'].setP([0.5])
        c['rt_T'].props['fC1'].params['fC1'].setP([0.5])
        c['vt_T'].props['fC1'].params['fC1'].setP([0.5])

        #c['rw'].params['pos'].setP([[0,0,0],[0,3,30]])
        #c['rw'].params['ogn'].setP([0,0,0])
        #c['rw'].addParam('offset','pos',[1,3],P=[18,-1,3])
        #c['rw'].addParam('pos1','pos',[3,3],P=[[0,0,0],[18,0,25],[22,0,29]],T=[0,0.9,1.0])
        #c['rw'].params['scl'].setP([1])
        #c['rw'].addParam('scl1','scl',[3,1],P=[10,5,0.8],T=[0,0.35,1.0])

    def meshStructure(self):
        afm = Airframe(self, 1) #0.2)

        idims = numpy.linspace(0.45,0.9,7)
        jdims = numpy.linspace(0,1,16)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlw_1:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mrw_1:'+str(i)+':'+str(j),'rw',[idims[i],1-jdims[j]],[idims[i+1],1-jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                if i is 0 or i is idims.shape[0]-1:
                    afm.addVertFlip('Mlw_2:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                    afm.addVertFlip('Mrw_2:'+str(i)+':'+str(j),'rw',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]])
                else:
                    afm.addVertFlip('Mlw_2a:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mlw_2b:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[0.15,0])
                    afm.addVertFlip('Mrw_2a:'+str(i)+':'+str(j),'rw',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mrw_2b:'+str(i)+':'+str(j),'rw',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]],w=[0.15,0])
        idims = numpy.linspace(0.18,0.45,6)
        for j in range(idims.shape[0]-1):
            afm.addVertFlip('Mlw_sec1:'+str(j),'lw',[idims[j],jdims[j]],[idims[j+1],jdims[j+1]])
            afm.addVertFlip('Mrw_sec1:'+str(j),'rw',[idims[j],1-jdims[j]],[idims[j+1],1-jdims[j+1]])
            afm.addVertFlip('Mlw_sec2:'+str(j),'lw',[idims[j],jdims[j]],[0.45,jdims[j]])
            afm.addVertFlip('Mrw_sec2:'+str(j),'rw',[idims[j],1-jdims[j]],[0.45,1-jdims[j]])

        idims = numpy.linspace(0.45,0.85,7)
        jdims = numpy.linspace(0,1,16)
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
        jdims = numpy.linspace(0,0.9,10)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlt_1:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mrt_1:'+str(i)+':'+str(j),'rt',[idims[i],1-jdims[j]],[idims[i+1],1-jdims[j]])
                afm.addVertFlip('Mvt_1:'+str(i)+':'+str(j),'vt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVertFlip('Mlt_2:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mrt_2:'+str(i)+':'+str(j),'rt',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]])
                afm.addVertFlip('Mvt_2:'+str(i)+':'+str(j),'vt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
        for i in range(idims.shape[0]):
                afm.addCtrVert('Mct_2:'+str(i)+':'+str(j),'lt','rt',idims[i])
        for i in range(idims.shape[0]-1):
            afm.addCtr('Mct_u:','lt','rt',0,[idims[i],idims[i+1]])
        for i in range(idims.shape[0]-1):
            afm.addCtr('Mct_l:','lt','rt',1,[1-idims[i],1-idims[i+1]])

        idims = numpy.linspace(0,1,4)
        jdims = numpy.linspace(0,1,20)
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

        afm.preview('conventional_pvw.dat')
        afm.mesh()
        afm.computeMesh('conventional_str.dat')

if __name__ == '__main__':

    import cProfile

    name = 'conventional'
    aircraft = Conventional()

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

    #aircraft.test_derivatives_dv()
    #aircraft.oml0.plot()
#    aircraft.meshStructure()

#    cProfile.run('aircraft.meshStructure()')
