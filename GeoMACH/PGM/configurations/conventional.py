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

        self.addComp('fu', Body(nx=12, ny=4, nz=2))
        self.addComp('lw', Wing(nx=4, nz=4, right=0))
        self.addComp('rw', Wing(nx=4, nz=4, left=0))
        self.addComp('lp', Wing(left=0, right=0))
        self.addComp('rp', Wing(left=0, right=0))
        self.addComp('ln', Shell(nx=4, ny=1, nz=4))
        self.addComp('rn', Shell(nx=4, ny=1, nz=4))
        self.addComp('lt', Wing(right=0))
        self.addComp('rt', Wing(left=0))
        self.addComp('vt', Wing(nx=2,right=0))

        self.separateComps()

        self.addComp('fu_n', Cone(c['fu'], 0))
        self.addComp('fu_t', Cone(c['fu'], 1))
        self.addComp('lp_lw', Junction(c['lw'], 1, 1, [1,0], c['lp'], mSide=0))
        self.addComp('lp_ln', Junction(c['ln'], 1, 2, [1,0], c['lp'], mSide=1))
        self.addComp('rp_rw', Junction(c['rw'], 1, 1, [1,0], c['rp'], mSide=0))
        self.addComp('rp_rn', Junction(c['rn'], 1, 2, [1,0], c['rp'], mSide=1))
        self.addComp('lw_fu', Junction(c['fu'], 2, 0, [2,1], c['lw'], mSide=0))
        self.addComp('rw_fu', Junction(c['fu'], 0, 2, [2,5], c['rw'], mSide=1))
        self.addComp('lt_fu', Junction(c['fu'], 2, 0, [1,9], c['lt'], mSide=0))
        self.addComp('rt_fu', Junction(c['fu'], 0, 2, [1,0], c['rt'], mSide=1))
        self.addComp('vt_fu', Junction(c['fu'], 1, 0, [0,8], c['vt'], mSide=0))

        self.assembleComponents()

        c['fu'].setm(0,1,[17,4,4,4,4,7,4,13,4,4,10,4])
        c['fu'].setm(0,0,[4,4,4,4])
        c['fu'].setm(1,0,[8,8])
        c['lw'].setm(0,1,[6,4,4,20])
        c['rw'].setm(0,1,[20,4,4,6])
        c['lt'].setm(0,1,[15])
        c['rt'].setm(0,1,[15])
        c['vt'].setm(0,1,[15])
        c['ln'].setm(0,1,[7])
        c['ln'].setm(5,1,[7])
        c['rn'].setm(0,1,[7])
        c['rn'].setm(5,1,[7])
        c['ln'].setm(0,0,[7])
        c['ln'].setm(2,0,[7])
        c['rn'].setm(0,0,[7])
        c['rn'].setm(2,0,[7])

        self.update()

        c['fu'].params['pos'].setP([[0,0,0],[50,0,0]])
        c['fu'].params['scl'].setP([0])
        c['fu'].addParam('nose','pos',[2,3],P=[[0,-1.1,0],[0,0,0]],T=[0,0.13],B=[False,True])
        c['fu'].addParam('rad','scl',[4,1],P=[0.65,2.6,2.6,0.35],T=[0,0.14,0.7,1.0],B=[False,True,True,False])
        c['fu'].addParam('flt1','flt',[2,4],P=[[0,0,0.5,0.5],[0,0,0.5,0.5]],T=[0.28,0.53])

        c['lw'].params['nor'].setP([0])
        c['lw'].params['pos'].setP([[0,0,0],[0,3,0]])
        c['lw'].params['ogn'].setP([0,0,0])
        c['lw'].params['scl'].setP([0])
        c['lw'].addParam('offset','pos',[1,3],P=[15,-1,3])
        c['lw'].addParam('pos1','pos',[2,3],P=[[0,0,0],[19,0,24]],T=[0,1.0])
        c['lw'].addParam('scl1','scl',[3,1],P=[10,4.5,1.8],T=[0,0.35,1.0])
        c['lw'].addParam('rakedp1','pos',[2,3],P=[[0,0,0],[3,2,0]],T=[0.85,1.0],B=[[True,True,False],[False,False,False]])
        c['lw'].addParam('rakeds','scl',[2,1],P=[0,-1.5],T=[0.85,1.0],B=[True,False])
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
        c['rw'].addParam('offset','pos',[1,3],P=[15,-1,-3])
        c['rw'].addParam('pos1','pos',[2,3],P=[[19,0,-24],[0,0,0]],T=[0,1.0])
        c['rw'].addParam('scl1','scl',[3,1],P=[1.8,4.5,10],T=[0,0.65,1.0])
        c['rw'].addParam('rakedp1','pos',[2,3],P=[[3,2,0],[0,0,0]],T=[0,0.15],B=[[False,False,False],[True,True,False]])
        c['rw'].addParam('rakeds','scl',[2,1],P=[-1.5,0],T=[0,0.15],B=[False,True])
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

        c['lt'].params['nor'].setP([0])
        c['lt'].params['pos'].setP([[0,0,0],[6,0,8]])
        c['lt'].params['scl'].setP([0])
        c['lt'].addParam('offset','pos',[1,3],P=[44,0,1.7])
        c['lt'].addParam('scl1','scl',[2,1],P=[4,1])
        c['lt'].addParam('rot1','rot',[2,3],P=[[0,10,0],[0,0,0]])

        c['rt'].params['nor'].setP([0])
        c['rt'].params['pos'].setP([[6,0,-8],[0,0,0]])
        c['rt'].params['scl'].setP([0])
        c['rt'].addParam('offset','pos',[1,3],P=[44,0,-1.7])
        c['rt'].addParam('scl1','scl',[2,1],P=[1,4])
        c['rt'].addParam('rot1','rot',[2,3],P=[[0,0,0],[0,-10,0]])

        c['vt'].params['nor'].setP([0])
        c['vt'].params['pos'].setP([[0,0,0],[6,8,0]])
        c['vt'].params['scl'].setP([0])
        c['vt'].addParam('nor1','nor',[1,3],P=[1,0,0])
        c['vt'].addParam('offset','pos',[1,3],P=[42,2.0,0])
        c['vt'].addParam('scl1','scl',[2,1],P=[5.8,2])
        c['vt'].addParam('rot1','rot',[2,3],P=[[0,10,0],[0,0,0]])
        #c['vt'].addParam('aileron1','shU',[2,1],P=[-0.12,0.0])
        #c['vt'].params['aileron1'].setT([0.0,0.25],0)
        #c['vt'].addParam('aileron2','shL',[2,1],P=[0.0,0.12])
        #c['vt'].params['aileron2'].setT([0.75,1.0],0)

        c['lp'].params['nor'].setP([0])
        c['lp'].params['pos'].setP([[0,0,0],[-2,-0.3,0]])
        c['lp'].params['scl'].setP([0])
        c['lp'].addParam('nor1','nor',[1,3],P=[1,0,0])
        c['lp'].addParam('offset','pos',[1,3],P=[20,-0.5,9])
        c['lp'].addParam('scl1','scl',[2,1],P=[4,3])

        c['ln'].params['pos'].setP([[0,0,0],[4.5,0,0]])
        c['ln'].params['scl'].setP([1.25])
        c['ln'].params['thk'].setP([0])
        c['ln'].addParam('offset','pos',[1,3],P=[15,-2.2,9])
        c['ln'].addParam('thk1','thk',[3,1],P=[0.15,0.4,0.15],B=[False,True,False])

        c['rp'].params['nor'].setP([0])
        c['rp'].params['pos'].setP([[0,0,0],[-2,-0.3,0]])
        c['rp'].params['scl'].setP([0])
        c['rp'].addParam('nor1','nor',[1,3],P=[1,0,0])
        c['rp'].addParam('offset','pos',[1,3],P=[20,-0.5,-9])
        c['rp'].addParam('scl1','scl',[2,1],P=[4,3])

        c['rn'].params['pos'].setP([[0,0,0],[4.5,0,0]])
        c['rn'].params['scl'].setP([1.25])
        c['rn'].params['thk'].setP([0])
        c['rn'].addParam('offset','pos',[1,3],P=[15,-2.2,-9])
        c['rn'].addParam('thk1','thk',[3,1],P=[0.15,0.4,0.15],B=[False,True,False])

        c['lw_fu'].params['mC1'].setP([1.5])
        c['lw_fu'].params['fC1'].setP([0.5])
        c['rw_fu'].params['mC1'].setP([1.5])
        c['rw_fu'].params['fC1'].setP([0.5])
        c['lt_fu'].params['mC1'].setP([0.1])
        c['lt_fu'].params['fC1'].setP([0.1])
        c['rt_fu'].params['mC1'].setP([0.1])
        c['rt_fu'].params['fC1'].setP([0.1])
        c['vt_fu'].params['mC1'].setP([0.01])
        c['vt_fu'].params['fC1'].setP([0.1])
        c['lp_ln'].params['mC1'].setP([0])
        c['rp_rn'].params['mC1'].setP([0])
        c['fu_n'].params['scl'].setP([0.02])
        c['fu_t'].params['scl'].setP([0.02])

        #c['rw'].params['pos'].setP([[0,0,0],[0,3,30]])
        #c['rw'].params['ogn'].setP([0,0,0])
        #c['rw'].addParam('offset','pos',[1,3],P=[18,-1,3])
        #c['rw'].addParam('pos1','pos',[3,3],P=[[0,0,0],[18,0,25],[22,0,29]],T=[0,0.9,1.0])
        #c['rw'].params['scl'].setP([1])
        #c['rw'].addParam('scl1','scl',[3,1],P=[10,5,0.8],T=[0,0.35,1.0])

        self.computePoints()
        ind = self.inds


        members = [
            [
                [[ind['lw'], 0,0], [1,0.5,0], [1,0.5,0.5], [0,0.5,0], [0,0.5,0.5]],
                [[ind['lw'], 1,0], [0,0.5,0], [0,0.5,0.5], [1,0.5,0], [1,0.5,0.5]],
                ],
            [
                [[ind['lw'], 0,0], [1,0.2,0.5], [1,0.5,0.5], [0,0.2,0.5], [0,0.5,0.5]],
                [[ind['lw'], 1,0], [0,0.8,0.5], [0,0.5,0.5], [1,0.8,0.5], [1,0.5,0.5]],
                ],
            ]

        afm = Airframe(self, members, 1.0)
        #afm.mesh()

        #self.meshStructure(members, lengths)

if __name__ == '__main__':

#    import cProfile
#    cProfile.run('Conventional()')

    name = 'conventional'
    aircraft = Conventional()

    der = aircraft.getDerivatives('lw', 'scl1', (1,0), FD=False)
    aircraft.oml0.addVars(['der'])
    aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)

    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)

    #aircraft.oml0.plot()
