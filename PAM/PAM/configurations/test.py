from __future__ import division
import numpy, time
import sys
from PAM.components import Wing, Body, Junction
from PAM.configurations import Configuration


class Test(Configuration):

    def __init__(self):
        super(Test,self).__init__() 

        self.addComp('body', Body(nx=5,ny=4,nz=5))
        self.addComp('wingL', Wing(nz=3,right=1))
        self.addComp('wingR', Wing(nz=3,left=1))
        self.addComp('body2', Body(bottom=1))

        self.separateComps()

        c = self.comps
        self.addComp('juncL', Junction(c['wingL'],0,c['body'],4,[1,1],[2,3]))
        self.addComp('juncR', Junction(c['wingR'],1,c['body'],2,[2,3],[1,1]))
        self.addComp('juncB', Junction(c['body2'],2,c['body'],3,[3,1],[1,3]))

        self.assembleComponents()

        c = self.comps
        c['wingR'].variables['offset'] = [1,0,-1.2]
        c['wingR'].variables['pos'][:,2] = numpy.linspace(-1,0,c['wingR'].Qs[0].shape[1])
        c['wingR'].variables['chord'][:] = 0.2
        c['wingL'].variables['offset'] = [1,0,1.2]
        c['wingL'].variables['pos'][:,2] = numpy.linspace(0,1,c['wingL'].Qs[0].shape[1])
        c['wingL'].variables['chord'][:] = 0.2
        c['body'].variables['pos'][:,0] = numpy.linspace(0,2,c['body'].Qs[2].shape[1])
        c['body'].variables['fillet'][:,2] = 1.0
        c['body'].variables['fillet'][:,3] = 1.0
        c['body'].variables['noseL'] = 0.5
        c['body'].variables['tailL'] = 0.5
        c['body2'].variables['offset'] = [0.5,1.2,0]
        c['body2'].variables['radii'][:] = 0.2
        c['body2'].variables['pos'][:,0] = numpy.linspace(0,1,c['body2'].Qs[2].shape[1])
        #c['juncB'].variables['shape'][3:8,3:8] -= 0.2

        self.computePoints()   


if __name__ == '__main__':

    name = 'test'
    aircraft = Test()
    #aircraft.buildStructure()
    #aircraft.writeStructure(name)
    aircraft.export.write2Tec(name)
    aircraft.export.write2STL(name)
    aircraft.export.write2TecC(name+'_C')
    #aircraft.export.write2IGES(name)
    #aircraft.export.write2EGADS(name+'_EGADS')
    aircraft.plot()
