from __future__ import division
import numpy, time
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration
from tecplot import Tecplot


class Conventional(Configuration):

    def __init__(self):
        super(Conventional,self).__init__() 

        self.addComp('fuse', Body([100,10,10,35,10,10,10,50,10,25,10,10],[25,25,25,25],[15]))
        self.addComp('wing', Wing([10,10,10,50],[10,35,10,10]))
        self.addComp('tail', Wing([30],[25]))
        self.addComp('fin', Wing([30],[25],half=True))

        self.separateComps()

        self.addComp('wingfuse', FullInterface(self.comps, 'wing', 0, 'fuse', 2, [2,1], [3,6]))
        self.addComp('tailfuse', FullInterface(self.comps, 'tail', 0, 'fuse', 2, [1,8], [2,10]))
        self.addComp('finfuse', HalfInterface(self.comps, 'fin', 0, 'fuse', 1, [0,8], [0,10]))

        self.assembleComponents()

        self.computePoints()

        c = self.comps
        c['fuse'].props['posx'].set([0,10],[0,1])
        c['fuse'].props['posy'].set([0,0],[0,1])
        c['fuse'].props['ry'].set([1,1],[0,1])
        c['fuse'].props['rz'].set([1,1],[0,1])
        
        c['wing'].offset = [4,-.4,1]
        c['wing'].props['posx'].set([0,0],[0,1])
        c['wing'].props['posy'].set([0,0],[0,1])
        c['wing'].props['posz'].set([0,8],[0,1])
        c['wing'].props['chord'].set([2,2],[0,1])



        t0 = time.time()
        self.computePoints()        
        print time.time() - t0


if __name__ == '__main__':

    name = 'conventional'
    aircraft = Conventional()
    #aircraft.buildStructure()
    #aircraft.writeStructure(name)
    aircraft.oml0.write2Tec(name)
    #aircraft.oml0.write2TecC(name+'_C')
    #aircraft.oml0.write2IGES(name)
    #aircraft.oml0.write2EGADS(name+'_EGADS')
    aircraft.plot()
