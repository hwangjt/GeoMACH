from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys

from GeoMACH.PGM.components import Wing, Body, Shell, Junction, Cone
from GeoMACH.PGM.configurations import Configuration
from GeoMACH.PSM import Airframe


class Conventional(Configuration):

    def __init__(self):
        super(Conventional,self).__init__()
        c = self.comps

        self.addComp('fu', Body(nx=9, ny=4, nz=12))
        self.addComp('lw', Wing(right=0))
        self.addComp('rw', Wing(left=0))
#        self.addComp('lp', Wing(left=0, right=0))
#        self.addComp('rp', Wing(left=0, right=0))
#        self.addComp('ln', Shell(nx=4, ny=1, nz=4))
#        self.addComp('rn', Shell(nx=4, ny=1, nz=4))
#        self.addComp('lt', Wing(right=0))
#        self.addComp('rt', Wing(left=0))
        self.addComp('lt', Wing(nx=2,right=0))
        self.addComp('rt', Wing(nx=2,right=0))
        self.addComp('ht', Wing(nx=6, nz=6))
#        self.addComp('le', Wing(nx=6, nz=8))

        self.separateComps()

        self.addComp('fu_n', Cone(c['fu'], 0))
        self.addComp('fu_t', Cone(c['fu'], 1))
#        self.addComp('lp_ht', Junction(c['ht'], 1, 1, [0,0], c['lp'], mSide=0))
#        self.addComp('lp_ln', Junction(c['ln'], 1, 2, [1,0], c['lp'], mSide=1))
#        self.addComp('rp_rw', Junction(c['rw'], 1, 1, [1,0], c['rp'], mSide=0))
#        self.addComp('rp_rn', Junction(c['rn'], 1, 2, [1,0], c['rp'], mSide=1))
        self.addComp('lw_fu', Junction(c['fu'], 2, 0, [2,1], c['lw'], mSide=0))
        self.addComp('rw_fu', Junction(c['fu'], 0, 2, [2,5], c['rw'], mSide=1))
#        self.addComp('lt_fu', Junction(c['fu'], 2, 0, [1,9], c['lt'], mSide=0))
#        self.addComp('rt_fu', Junction(c['fu'], 0, 2, [1,0], c['rt'], mSide=1))
        self.addComp('lt_fu', Junction(c['fu'], 1, 0, [10,5], c['lt'], mSide=0))
        self.addComp('rt_fu', Junction(c['fu'], 1, 0, [0,5], c['rt'], mSide=0))
        self.addComp('lt_ht', Junction(c['ht'], 1, 3, [3,1], c['lt'], mSide=1))
        self.addComp('rt_ht', Junction(c['ht'], 1, 3, [1,1], c['rt'], mSide=1))
#        self.addComp('le_lf', Junction(c['fu'], 1, 0, [9,5], c['le'], mSide=0))

        self.assembleComponents()

        c['fu'].setm(0,1,[25,4,7,4,13,4,4,10,4])
        c['fu'].setm(0,0,[4,4,4,4])
        c['fu'].setm(1,0,[4,4,4,4])
        c['lw'].setm(0,1,[30])
        c['rw'].setm(0,1,[30])
#        c['lt'].setm(0,1,[15])
#        c['rt'].setm(0,1,[15])
        c['lt'].setm(0,1,[15])
        c['rt'].setm(0,1,[15])
        c['ht'].setm(1,0,[5])
        c['ht'].setm(1,1,[5])
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

        self.update()

        c['fu'].params['pos'].setP([[0,0,0],[50,0,0]])
        c['fu'].params['scl'].setP([0])
        c['fu'].addParam('nose','pos',[2,3],P=[[0,-1.1,0],[0,0,0]],T=[0,0.13],B=[False,True])
        c['fu'].addParam('rad','scl',[4,3],P=[[1.0,0.7,0],[4.2,1.73,0],[4.2,1.73,0],[4.2,0.4,0]],T=[0,0.14,0.85,1.0],B=[False,True,False,False])
        c['fu'].addParam('flt1','flt',[2,4],P=[[0,0.2,0,0.2],[0,0.2,0,0.2]],T=[0.0,1.0])
        c['fu'].addParam('flt2','flt',[2,4],P=[[0,0,0.6,0],[0,0,0.6,0]],T=[0.35,0.55])

        c['lw'].params['nor'].setP([0])
        c['lw'].params['pos'].setP([[0,0,0],[0,3,0]])
        c['lw'].params['ogn'].setP([0,0,0])
        c['lw'].params['scl'].setP([0])
        c['lw'].addParam('offset','pos',[1,3],P=[21,-0.7,4.4])
        c['lw'].addParam('pos1','pos',[2,3],P=[[0,0,0],[1.0,0,24]],T=[0,1.0])
        c['lw'].addParam('scl1','scl',[2,1],P=[5,1.8],T=[0,1.0])
#        c['lw'].addParam('rakedp1','pos',[2,3],P=[[0,0,0],[3,2,0]],T=[0.85,1.0],B=[[True,True,False],[False,False,False]])
#        c['lw'].addParam('rakeds','scl',[2,1],P=[0,-1.5],T=[0.85,1.0],B=[True,False])
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

        c['rw'].params['nor'].setP([0])
        c['rw'].params['pos'].setP([[0,3,0],[0,0,0]])
        c['rw'].params['ogn'].setP([0,0,0])
        c['rw'].params['scl'].setP([0])
        c['rw'].addParam('offset','pos',[1,3],P=[21,-0.7,-4.4])
        c['rw'].addParam('pos1','pos',[2,3],P=[[1.0,0,-24],[0,0,0]],T=[0,1.0])
        c['rw'].addParam('scl1','scl',[2,1],P=[1.8,5],T=[0,1.0])
#        c['rw'].addParam('rakedp1','pos',[2,3],P=[[3,2,0],[0,0,0]],T=[0,0.15],B=[[False,False,False],[True,True,False]])
#        c['rw'].addParam('rakeds','scl',[2,1],P=[-1.5,0],T=[0,0.15],B=[False,True])
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

#        c['lt'].params['nor'].setP([0])
#        c['lt'].params['pos'].setP([[0,0,0],[6,0,8]])
#        c['lt'].params['scl'].setP([0])
#        c['lt'].addParam('offset','pos',[1,3],P=[44,0,1.7])
#        c['lt'].addParam('scl1','scl',[2,1],P=[4,1])
#        c['lt'].addParam('rot1','rot',[2,3],P=[[0,10,0],[0,0,0]])

#        c['rt'].params['nor'].setP([0])
#        c['rt'].params['pos'].setP([[6,0,-8],[0,0,0]])
#        c['rt'].params['scl'].setP([0])
#        c['rt'].addParam('offset','pos',[1,3],P=[44,0,-1.7])
#        c['rt'].addParam('scl1','scl',[2,1],P=[1,4])
#        c['rt'].addParam('rot1','rot',[2,3],P=[[0,0,0],[0,-10,0]])

        c['lt'].params['nor'].setP([0])
        c['lt'].params['pos'].setP([[0,0,0],[4.5,6,0]])
        c['lt'].params['scl'].setP([0])
        c['lt'].addParam('nor1','nor',[1,3],P=[1,0,0])
        c['lt'].addParam('offset','pos',[1,3],P=[42,1.7,2.5])
        c['lt'].addParam('scl1','scl',[2,1],P=[5.8,2])
        c['lt'].addParam('rot1','rot',[2,3],P=[[0,5,0],[0,0,0]])

        c['rt'].params['nor'].setP([0])
        c['rt'].params['pos'].setP([[0,0,0],[4.5,6,0]])
        c['rt'].params['scl'].setP([0])
        c['rt'].addParam('nor1','nor',[1,3],P=[1,0,0])
        c['rt'].addParam('offset','pos',[1,3],P=[42,1.7,-2.5])
        c['rt'].addParam('scl1','scl',[2,1],P=[5.8,2])
        c['rt'].addParam('rot1','rot',[2,3],P=[[0,5,0],[0,0,0]])

        c['ht'].params['nor'].setP([0])
        c['ht'].params['pos'].setP([[0,0,0],[0,0,0]])
        c['ht'].params['ogn'].setP([0,0,0])
        c['ht'].params['scl'].setP([0])
        c['ht'].addParam('offset','pos',[1,3],P=[44,7.7,0])
        c['ht'].addParam('pos1','pos',[3,3],P=[[3,0,-8],[0,0,0],[3,0,8]],T=[0,0.5,1.0])
        c['ht'].addParam('scl1','scl',[3,1],P=[2,4,2])

#        c['le'].params['nor'].setP([0])
#        c['le'].params['pos'].setP([[0,0,0],[0,0,0]])
#        c['le'].params['scl'].setP([0])
#        c['le'].addParam('nor1','nor',[1,3],P=[1,0,0])
#        c['le'].addParam('offset','pos',[1,3],P=[42,2,2.2])
#        c['le'].addParam('pos','pos',[3,3],P=[[0,0,0],[0,2,-1],[0,0,-2]],D=[[0,1,0],[0,0,1],[0,1,0]],T=[0,0.5,1])
#        c['le'].addParam('scl1','scl',[1,1],P=[2])
#        c['le'].addParam('rot1','rot',[2,3],P=[[0,10,0],[0,0,0]])

        #c['vt'].addParam('aileron1','shU',[2,1],P=[-0.12,0.0])
        #c['vt'].params['aileron1'].setT([0.0,0.25],0)
        #c['vt'].addParam('aileron2','shL',[2,1],P=[0.0,0.12])
        #c['vt'].params['aileron2'].setT([0.75,1.0],0)

#        c['lp'].params['nor'].setP([0])
#        c['lp'].params['pos'].setP([[0,0,0],[-2,-0.3,0]])
#        c['lp'].params['scl'].setP([0])
#        c['lp'].addParam('nor1','nor',[1,3],P=[1,0,0])
#        c['lp'].addParam('offset','pos',[1,3],P=[20,-0.5,9])
#        c['lp'].addParam('scl1','scl',[2,1],P=[4,3])

#        c['ln'].params['pos'].setP([[0,0,0],[4.5,0,0]])
#        c['ln'].params['scl'].setP([1.25])
#        c['ln'].params['thk'].setP([0])
#        c['ln'].addParam('offset','pos',[1,3],P=[15,-2.2,9])
#        c['ln'].addParam('thk1','thk',[3,1],P=[0.15,0.4,0.15],B=[False,True,False])

#        c['rp'].params['nor'].setP([0])
#        c['rp'].params['pos'].setP([[0,0,0],[-2,-0.3,0]])
#        c['rp'].params['scl'].setP([0])
#        c['rp'].addParam('nor1','nor',[1,3],P=[1,0,0])
#        c['rp'].addParam('offset','pos',[1,3],P=[20,-0.5,-9])
#        c['rp'].addParam('scl1','scl',[2,1],P=[4,3])

#        c['rn'].params['pos'].setP([[0,0,0],[4.5,0,0]])
#        c['rn'].params['scl'].setP([1.25])
#        c['rn'].params['thk'].setP([0])
#        c['rn'].addParam('offset','pos',[1,3],P=[15,-2.2,-9])
#        c['rn'].addParam('thk1','thk',[3,1],P=[0.15,0.4,0.15],B=[False,True,False])

        c['lw_fu'].params['mC1'].setP([1.5])
        c['lw_fu'].params['fC1'].setP([0.5])
        c['rw_fu'].params['mC1'].setP([1.5])
        c['rw_fu'].params['fC1'].setP([0.5])
#        c['lt_fu'].params['mC1'].setP([0.1])
#        c['lt_fu'].params['fC1'].setP([0.1])
#        c['rt_fu'].params['mC1'].setP([0.1])
#        c['rt_fu'].params['fC1'].setP([0.1])
        c['lt_fu'].params['mC1'].setP([0.01])
        c['lt_fu'].params['fC1'].setP([0.1])
        c['rt_fu'].params['mC1'].setP([0.01])
        c['rt_fu'].params['fC1'].setP([0.1])
#        c['lt_ht'].params['mC1'].setP([0.01])
#        c['lt_ht'].params['fC1'].setP([0.1])
#        c['rt_ht'].params['mC1'].setP([0.01])
#        c['rt_ht'].params['fC1'].setP([0.1])
#        c['lp_ln'].params['mC1'].setP([0])
#        c['rp_rn'].params['mC1'].setP([0])
        c['fu_n'].params['scl'].setP([0.02])
        c['fu_t'].params['scl'].setP([0.02])


        #c['rw'].params['pos'].setP([[0,0,0],[0,3,30]])
        #c['rw'].params['ogn'].setP([0,0,0])
        #c['rw'].addParam('offset','pos',[1,3],P=[18,-1,3])
        #c['rw'].addParam('pos1','pos',[3,3],P=[[0,0,0],[18,0,25],[22,0,29]],T=[0,0.9,1.0])
        #c['rw'].params['scl'].setP([1])
        #c['rw'].addParam('scl1','scl',[3,1],P=[10,5,0.8],T=[0,0.35,1.0])

        self.computePoints()

    def meshStructure(self):
        afm = Airframe(self, 1)

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

#    import cProfile
#    cProfile.run('Conventional()')

    name = 'D8'
    aircraft = Conventional()

    der = aircraft.getDerivatives('lw', 'scl1', (1,0), FD=False)
    aircraft.oml0.addVars(['der'])
    #aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)
    aircraft.oml0.Q[:,6] = numpy.sum(der*der,1)**0.5

    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)

    aircraft.meshStructure()

    #aircraft.oml0.plot()
