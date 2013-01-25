from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys
from PAM.components import Wing, Body, Shell, Junction, Cone
from PAM.configurations import Configuration


class Test(Configuration):

    def __init__(self):
        super(Test,self).__init__()

        self.addComp('body', Body(nx=5,ny=4,nz=4,bottom=2))
        self.addComp('wing', Wing())

        self.separateComps()

        c = self.comps
        self.addComp('nose', Cone(c['body'], 0))
        self.addComp('tail', Cone(c['body'], 1))

        self.assembleComponents()

        c = self.comps

        c['wing'].params['pos'].setP([[0,0,0],[0,0,4]])
        c['wing'].addParam('offset','pos',[1,1],P=[2])
        c['wing'].addParam('winglet','pos',[3,3],P=[[0,0,0],[0,0,0],[0,1,0]],Tdim=0,T=[0,0.8,1],Ddim=0,D=[[0,0,1],[0,0,1],[0,1,0]],Bdim=0,B=numpy.ones((3,3),bool))

        self.computePoints()
        print c['wing'].variables['pos']


if __name__ == '__main__':

    name = 'test'
    aircraft = Test()
    aircraft.oml0.addVars(['der'])

    d1 = aircraft.getDerivatives('body', 'shT', [0,0], FD=False)
    d2 = aircraft.getDerivatives('body', 'shT', [0,0], FD=True) 
    aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(d1)
    #aircraft.computePoints()

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
