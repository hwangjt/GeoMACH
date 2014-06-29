from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys

from GeoMACH.PGM.components import Wing, Body, Shell, Junction, Cone, Tip
from GeoMACH.PGM.configurations import Configuration
from GeoMACH.PSM import Airframe


class D8(Configuration):

    def define_primitive_comps(self):
        comps = {'fu': Body(nx=9, ny=4, nz=12),
                 'lw': Wing(right=0),
                 'rw': Wing(left=0),
                 'lt': Wing(nx=2,right=0,left=0),
                 'rt': Wing(nx=2,right=0,left=0),
                 'ht': Wing(nx=6, nz=6),
                 }
        return comps

    def define_interpolant_comps(self):
        c = self.comps
        comps = {'fu_n': Cone(c['fu'], 0),
                 'fu_t': Cone(c['fu'], 1),
                 'lw_T': Tip(c['lw'], 1),
                 'rw_T': Tip(c['rw'], 0),
                 'ht_R': Tip(c['ht'], 0),
                 'ht_L': Tip(c['ht'], 1),
                 'lw_fu': Junction(c['fu'], 'lft', 0, [2,1], c['lw'], mSide=0),
                 'rw_fu': Junction(c['fu'], 'rgt', 2, [2,5], c['rw'], mSide=1),
                 'lt_fu': Junction(c['fu'], 'top', 0, [10,5], c['lt'], mSide=0),
                 'rt_fu': Junction(c['fu'], 'top', 0, [0,5], c['rt'], mSide=0),
                 'lt_ht': Junction(c['ht'], 'low', 3, [3,1], c['lt'], mSide=1),
                 'rt_ht': Junction(c['ht'], 'low', 3, [1,1], c['rt'], mSide=1),
                 }
        return comps

    def define_dvs(self):
        self.add_dv('sweep', [1], val=19, lower=15, upper=23)
        self.add_dv('shape', [1], val=0.1, lower=0.0, upper=1.0)
        self.add_dv('rot', [1], val=0.1, lower=0.0, upper=1.0)

    def apply_dvs(self):
        return [0], [0], [0]

    def set_oml_resolution(self):
        c = self.comps

        c['fu'].faces['rgt'].setm(1,[25,4,7,4,13,4,4,10,4])
        c['fu'].faces['rgt'].setm(0,[4,4,4,4])
        c['fu'].faces['top'].setm(0,[4,4,4,4])
        c['lw'].faces['upp'].setm(1,[30])
        c['rw'].faces['upp'].setm(1,[30])
#        c['lt'].setm(0,1,[15])
#        c['rt'].setm(0,1,[15])
        c['lt'].faces['upp'].setm(1,[15])
        c['rt'].faces['upp'].setm(1,[15])
        c['ht'].faces['low'].setm(0,[5])
        c['ht'].faces['low'].setm(1,[5])
        c['lw_T'].faces['def'].setm(0,[10])
        c['rw_T'].faces['def'].setm(0,[10])
        c['ht_L'].faces['def'].setm(0,[10])
        c['ht_R'].faces['def'].setm(0,[10])
#        c['le'].setm(1,0,[5])
#        c['le'].setm(1,1,[5])
#        c['ln'].setm(0,1,[7])
#        c['ln'].setm(5,1,[7])
#        c['rn'].setm(0,1,[7])
#        c['rn'].setm(5,1,[7])
#        c['ln'].setm(0,0,[7])
#        c['ln'].setm(2,0,[7])
#        c['rn'].setm(0,0,[7])
#        c['rn'].setm(2,0,[7])


    def define_oml_parameters(self):
        c = self.comps

        c['fu'].props['nor'].addParam('nor1',[1,3],P=[1,0,1])
        c['fu'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[50,0,0]])
        c['fu'].props['pos'].addParam('nose',[2,3],P=[[0,0.75,0],[0,0,0]],T=[0,0.13],B=[False,True])
        c['fu'].props['scl'].addParam('rad',[4,3],P=[[0.7,0.3,0],[4.2,1.73,0],[4.2,1.73,0],[4.2,0.3,0]],T=[0,0.14,0.85,1.0],B=[False,True,False,False])
        c['fu'].props['flt'].addParam('flt1',[2,4],P=[[0,0.2,0,0.2],[0,0.2,0,0.2]],T=[0.0,1.0])
        c['fu'].props['flt'].addParam('flt2',[2,4],P=[[0,0,0.6,0],[0,0,0.6,0]],T=[0.4,0.6])

        c['lw'].props['pos'].addParam('offset',[1,3],P=[21,-0.7,4.2])
        c['lw'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[1.0,3,24]])
        c['lw'].props['scl'].addParam('scl1',[2,1],P=[5,1.8],T=[0,1.0])

        c['rw'].props['pos'].addParam('offset',[1,3],P=[21,-0.7,-4.4])
        c['rw'].props['pos'].addParam('pos1',[2,3],P=[[1.0,3,-24],[0,0,0]])
        c['rw'].props['scl'].addParam('scl1',[2,1],P=[1.8,5],T=[0,1.0])

        c['lt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[4.5,5.9,0]])
        c['lt'].props['nor'].addParam('nor1',[1,3],P=[1,0,0])
        c['lt'].props['pos'].addParam('offset',[1,3],P=[42,1.5,2.6])
        c['lt'].props['scl'].addParam('scl1',[2,1],P=[5.8,2])
        c['lt'].props['rot'].addParam('rot1',[2,3],P=[[15,8,0],[0,-5,0]])

        c['rt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[4.5,6,0]])
        c['rt'].props['nor'].addParam('nor1',[1,3],P=[1,0,0])
        c['rt'].props['pos'].addParam('offset',[1,3],P=[42,1.7,-2.5])
        c['rt'].props['scl'].addParam('scl1',[2,1],P=[5.8,2])
        c['rt'].props['rot'].addParam('rot1',[2,3],P=[[0,5,0],[0,0,0]])

        c['ht'].props['pos'].addParam('offset',[1,3],P=[44.3,7.7,0])
        c['ht'].props['pos'].addParam('pos1',[3,3],P=[[3,0,-8],[0,0,0],[3,0,8]],T=[0,0.5,1.0])
        c['ht'].props['scl'].addParam('scl1',[3,1],P=[2,5,2])

        c['lw_fu'].props['mC1'].params['mC1'].setP([0.01])
        c['lw_fu'].props['fC1'].params['fC1'].setP([1])
        c['rw_fu'].props['mC1'].params['mC1'].setP([1.5])
        c['rw_fu'].props['fC1'].params['fC1'].setP([0.5])
        c['lt_fu'].props['mC1'].params['mC1'].setP([0.01])
        c['lt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['rt_fu'].props['mC1'].params['mC1'].setP([0.01])
        c['rt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['lt_ht'].props['mC1'].params['mC1'].setP([0.01])
        c['lt_ht'].props['fC1'].params['fC1'].setP([0.1])
        c['rt_ht'].props['mC1'].params['mC1'].setP([0.01])
        c['rt_ht'].props['fC1'].params['fC1'].setP([0.1])
        c['fu_n'].props['fC1'].params['fC1'].setP([1])
        c['fu_t'].props['fC1'].params['fC1'].setP([2])
        c['lw_T'].props['fC1'].params['fC1'].setP([0.9])
        c['rw_T'].props['fC1'].params['fC1'].setP([0.9])
        c['ht_L'].props['fC1'].params['fC1'].setP([0.9])
        c['ht_R'].props['fC1'].params['fC1'].setP([0.9])

    def meshStructure(self):
        afm = Airframe(self, 0.4)

        idims = numpy.linspace(0.25,0.65,7)
        jdims = numpy.linspace(0,0.99,16)
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
#        idims = numpy.linspace(0.18,0.45,6)
#        for j in range(idims.shape[0]-1):
#            afm.addVertFlip('Mlw_sec1:'+str(j),'lw',[idims[j],jdims[j]],[idims[j+1],jdims[j+1]])
#            afm.addVertFlip('Mrw_sec1:'+str(j),'rw',[idims[j],1-jdims[j]],[idims[j+1],1-jdims[j+1]])
#            afm.addVertFlip('Mlw_sec2:'+str(j),'lw',[idims[j],jdims[j]],[0.45,jdims[j]])
#            afm.addVertFlip('Mrw_sec2:'+str(j),'rw',[idims[j],1-jdims[j]],[0.45,1-jdims[j]])

        idims = numpy.linspace(0.25,0.65,7)
        jdims = numpy.linspace(0,0.9,16)
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
#        afm.addCtrVert('Mcw_sec:'+str(i)+':'+str(j),'lw','rw',0.18)

        idims = numpy.linspace(0.25,0.65,2)
        jdims = numpy.linspace(0,0.9,10)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlt_1:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mrt_1:'+str(i)+':'+str(j),'rt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVertFlip('Mlt_2:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mrt_2:'+str(i)+':'+str(j),'rt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
        jdims = numpy.linspace(0.1,0.9,10)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mvt_1:'+str(i)+':'+str(j),'ht',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVertFlip('Mvt_2:'+str(i)+':'+str(j),'ht',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
#        for i in range(idims.shape[0]):
#                afm.addCtrVert('Mct_2:'+str(i)+':'+str(j),'lt','rt',idims[i])
#        for i in range(idims.shape[0]-1):
#            afm.addCtr('Mct_u:','lt','rt',0,[idims[i],idims[i+1]])
#        for i in range(idims.shape[0]-1):
#            afm.addCtr('Mct_l:','lt','rt',1,[1-idims[i],1-idims[i+1]])

        idims = numpy.linspace(0,1,4)
        jdims = numpy.linspace(0,1,20)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVert('Mfu_1:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[0,2])
                afm.addVert('Mfu_2:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[1,3])
                afm.addVert('Mfu_3:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[2,0])
                afm.addVert('Mfu_4:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[3,1])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVert('Mfu_5:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[0,2])
                afm.addVert('Mfu_6:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[1,3])
                afm.addVert('Mfu_7:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[2,0])
                afm.addVert('Mfu_8:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[3,1])
        for j in range(jdims.shape[0]-1):
            afm.addVertFlip('Mfu_0:'+str(j),'fu',[0.4,jdims[j]],[0.4,jdims[j+1]],w=[1.0,0.0],i=[0,2])

        afm.preview('D8_pvw.dat')
        afm.mesh()
        afm.computeMesh('D8_str.dat')


if __name__ == '__main__':

    name = 'D8'
    aircraft = D8()

    #der = aircraft.getDerivatives('lw', 'scl1', (1,0), FD=False)
    #aircraft.oml0.addVars(['der'])
    #aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)
    #aircraft.oml0.Q[:,6] = numpy.sum(der*der,1)**0.5

    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)
    aircraft.oml0.write2IGES(name)

    aircraft.meshStructure()

    #aircraft.oml0.plot()
