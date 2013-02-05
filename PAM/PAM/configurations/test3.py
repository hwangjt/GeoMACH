from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys
from PAM.components import Wing, Body, Shell, Junction, Cone
from PAM.configurations import Configuration


class Test(Configuration):

    def __init__(self):
        super(Test,self).__init__()
        c = self.comps

        self.addComp('body', Body(nx=5,ny=4,nz=4,bottom=2))
        self.addComp('wing', Wing())

        self.separateComps()

        self.addComp('nose', Cone(c['body'], 0))
        self.addComp('tail', Cone(c['body'], 1))

        self.assembleComponents()

        c['wing'].setm(0,0,[20])
        c['wing'].setn(0,0,[40])
        c['wing'].setm(0,1,[30])
        c['wing'].setn(0,1,[60])
        self.update()

        c['wing'].params['pos'].setP([[0,0,0],[0,0,0]])
        c['wing'].addParam('offset','pos',[1,1],P=[2])
        c['wing'].addParam('winglet','pos',[3,3],P=[[0,0,0],[0,0,3.5],[0,1,4]],T=[0,0.6,1],D=[[0,0,1],[0,0,1],[0,1,0]],B=numpy.ones((3,3),bool))
        #c['wing'].addParam('aileron1','shU',[1,1],P=[0.03])
        #c['wing'].addParam('aileron2','shL',[1,1],P=[-0.03])
        if 1:
            c['wing'].addParam('aileron1','shU',[2,4],P=[[0,0.04,0.04,0],[0,0,0,0]])
            c['wing'].params['aileron1'].setT([0.0,0.5],0)
            c['wing'].params['aileron1'].setT([0.1,0.13,0.3,0.33],1)
            c['wing'].params['aileron1'].setB([[0,0,0,0],[1,1,1,1]],0)
            c['wing'].params['aileron1'].setB([[1,1,1,1],[1,1,1,1]],1)
        if 1:
            c['wing'].addParam('aileron2','shL',[2,4],P=[[0,0,0,0],[0,-0.04,-0.04,0]])
            c['wing'].params['aileron2'].setT([0.5,1.0],0)
            c['wing'].params['aileron2'].setT([0.1,0.13,0.3,0.33],1)
            c['wing'].params['aileron2'].setB([[1,1,1,1],[0,0,0,0]],0)
            c['wing'].params['aileron2'].setB([[1,1,1,1],[1,1,1,1]],1)

        self.computePoints()



if __name__ == '__main__':

    name = 'test'
    aircraft = Test()
    aircraft.oml0.addVars(['der'])

    d1 = aircraft.getDerivatives('wing', 'aileron1', (0,1), FD=True)
    d2 = aircraft.getDerivatives('wing', 'aileron1', (0,2), FD=True)
    d3 = aircraft.getDerivatives('wing', 'aileron2', (1,1), FD=True)
    d4 = aircraft.getDerivatives('wing', 'aileron2', (1,2), FD=True)
    aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(d1+d2+d3+d4)

    #aircraft.runDerivativeTest('body')
    #aircraft.runDerivativeTest('nose')
    #aircraft.runDerivativeTest('wing')
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)
    aircraft.oml0.write2STL(name)
    aircraft.oml0.write2IGES(name)

    if 0:
        t = Tecplot()
        t.importDataSet(name)
        t.setTransparency(False)
        t.setRelCameraPosition(-20, 10, 20, -140)
        t.writeImage('screen')
        t.runTecplot()
        exit()

    aircraft.oml0.plot()
