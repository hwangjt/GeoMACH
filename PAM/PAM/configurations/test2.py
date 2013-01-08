from __future__ import division
import numpy, time
import sys
from PAM.components import Wing, Body, Shell, Junction, Cone
from PAM.configurations import Configuration


class Test(Configuration):

    def __init__(self):
        super(Test,self).__init__()

        self.addComp('body', Body(nx=5,ny=4,nz=4,bottom=2))
        self.addComp('body2', Body(nx=1,ny=1,nz=2,bottom=1))
        self.addComp('wingL', Wing(nx=1,nz=1,right=1))
        self.addComp('wingR', Wing(nx=1,nz=1,left=1))
        self.addComp('shell', Shell())

        self.separateComps()

        c = self.comps
        self.addComp('nose', Cone(c['body'], 0))
        self.addComp('tail', Cone(c['body'], 1))
        self.addComp('nose2', Cone(c['body2'], 0))
        self.addComp('tail2', Cone(c['body2'], 1))
        self.addComp('wbL', Junction(c['body'], 2, 0, [1,1], c['wingL'], mSide=0))
        self.addComp('wbR', Junction(c['body'], 0, 2, [1,1], c['wingR'], mSide=1))
        self.addComp('bb', Junction(c['body'], 1, 3, [1,0], c['body2'], c['nose2'], c['tail2']))

        self.assembleComponents()

        c = self.comps

        c['wingL'].variables['offset'] = [0.5,0,1]
        c['wingL'].variables['scale'][:,:] = 0.2
        c['wingL'].variables['pos'][:,2] = numpy.linspace(0,1,c['wingL'].Qs[0].shape[1])
        c['wingR'].variables['offset'] = [0.5,0,-1]
        c['wingR'].variables['scale'][:,:] = 0.2
        c['wingR'].variables['pos'][:,2] = numpy.linspace(-1,0,c['wingR'].Qs[0].shape[1])
        c['body'].variables['pos'][:,0] = numpy.linspace(0,1,c['body'].Qs[0].shape[1])
        c['nose'].variables['scale'][0] = 0.25
        c['nose'].variables['f0'][0] = 3.0
        c['tail'].variables['scale'][0] = 0.25
        c['tail'].variables['f0'][0] = 3.0
        c['body2'].variables['offset'] = [0.5,1.1,0]
        c['body2'].variables['pos'][:,0] = numpy.linspace(-0.2,0.2,c['body2'].Qs[0].shape[1])
        c['body2'].variables['scale'][:,:] = 0.2
        c['nose2'].variables['scale'][0] = 0.25
        c['nose2'].variables['f0'][0] = 3.0
        c['tail2'].variables['scale'][0] = 0.25
        c['tail2'].variables['f0'][0] = 3.0
        c['shell'].variables['offset'] = [0.5,-1.5,0]
        c['shell'].variables['pos'][:,0] = numpy.linspace(-0.5,0.5,c['shell'].Qs[0].shape[1])
        c['shell'].variables['pos'][:,2] = numpy.linspace(-0.5,0.5,c['shell'].Qs[0].shape[1])
        c['shell'].variables['scale'][:,:] = 0.2

        self.computePoints()


if __name__ == '__main__':

    name = 'test'
    aircraft = Test()

    if 0:
        aircraft.oml0.addVars(['dQxdc','dQydc','dQzdc'])
        aircraft.oml0.addVars(['dQxdc2','dQydc2','dQzdc2'])
        d1 = aircraft.getDerivatives('bb','shape',(5,10),FD=True)
        d2 = aircraft.getDerivatives('bb','shape',(5,10),FD=True)
        print numpy.linalg.norm(d1)
        print numpy.linalg.norm(d2)
        print numpy.linalg.norm(d2-d1)
        aircraft.oml0.Q[:,6:9] = d1
        aircraft.oml0.Q[:,9:12] = d2
        aircraft.computePoints()
        aircraft.oml0.write2Tec(name)
        exit()

    aircraft.runDerivativeTest('wingL')
    aircraft.runDerivativeTest('body')
    aircraft.runDerivativeTest('shell')
    aircraft.runDerivativeTest('wbL')
    aircraft.runDerivativeTest('bb')
    aircraft.runDerivativeTest('nose')
    aircraft.runDerivativeTest('nose2')
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)
    aircraft.oml0.write2STL(name)
    aircraft.oml0.write2IGES(name)
    aircraft.oml0.plot()
