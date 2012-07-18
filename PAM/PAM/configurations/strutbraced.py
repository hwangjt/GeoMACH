from __future__ import division
import numpy
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration


class Strutbraced(Configuration):

    def __init__(self):
        self.comps = {}
        self.keys = []

        self.addComp('fuse', Body([90,10,10,10,10,10,10,45,10,10,10,40,10,10],[10,25,25,25,25],[15]))
        self.addComp('wing', Wing([30,10,10,50],[10,10,10,10]))
        self.addComp('strut', Wing([30],[10],opentip=True))
        self.addComp('vtail', Wing([40,10,10,10],[10,10,40,10],half=True))
        self.addComp('htail', Wing([20],[40]))

        self.separateComps()

        self.addComp('wingfuse', FullInterface(self.comps, 'wing', 0, 'fuse', 2, [0,1], [1,6]))
        self.addComp('strutfuse', FullInterface(self.comps, 'strut', 0, 'fuse', 2, [3,2], [4,4]))
        self.addComp('strutwing', FullInterface(self.comps, 'strut', 1, 'wing', 1, [2,1], [0,2]))
        self.addComp('vtailfuse', HalfInterface(self.comps, 'vtail', 0, 'fuse', 1, [0,8], [0,13]))
        self.addComp('htailvtail', FullInterface(self.comps, 'htail', 0, 'vtail', 0, [1,2], [3,1]))

        self.assembleComponents()

        self.comps['htailvtail'].setC1('surf', 0, j=-1, v=-1, val=False)
        self.comps['htailvtail'].setC1('edge', 0, j=-1, v=-1)
        self.comps['htailvtail'].extraDOFs.append([None,-1])
        self.updateComponents()

        c = self.comps
        c['fuse'].setSections(sections=[2,3,4], t1U=0.5, t2U=0.5)
        c['fuse'].props['posx'].set([0,10],[0,1])
        c['fuse'].props['posy'].set([0.3,0.5,0.5,0.9],[0,0.15,0.75,1.0],w=[1.0,0,0,0],d=[1,0,0,0])
        c['fuse'].props['ry'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])
        c['fuse'].props['rz'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])

        c['wing'].offset[:] = [3.7, 0.94, 0.5]
        c['wing'].setAirfoil("rae2822.dat")
        c['wing'].props['posx'].set([0,0],[0,1])
        c['wing'].props['posy'].set([0,0],[0,1])
        c['wing'].props['posz'].set([0,8],[0,1])
        c['wing'].props['prpx'].set([0,0],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([1,0.3],[0,1],w=[1,0])

        c['strut'].offset[:] = [3.95, 0.25, 0.5]
        c['strut'].props['posx'].set([0,0],[0,1])
        c['strut'].props['posy'].set([0,0.6],[0,1])
        c['strut'].props['posz'].set([0,3.2],[0,1])
        c['strut'].props['prpx'].set([1,1],[0,1])
        c['strut'].props['prpy'].set([0,0],[0,1])
        c['strut'].props['chord'].set([0.5,0.5],[0,1])

        c['vtail'].offset[:] = [8.0, 1, 0] #0.85, 0]
        c['vtail'].props['posx'].set([0,1.5],[0,1])
        c['vtail'].props['posy'].set([0,1.4],[0,1])
        c['vtail'].props['posz'].set([0,0],[0,1])
        c['vtail'].props['prpx'].set([1,1],[0,1])
        c['vtail'].props['prpy'].set([0,0],[0,1])
        c['vtail'].props['chord'].set([1.3,0.9],[0,1])

        c['htail'].offset[:] = [9.25, 2.1, 0.1]
        c['htail'].props['posx'].set([0,0.6],[0,1],w=[0.2,0])
        c['htail'].props['posy'].set([0,0],[0,1])
        c['htail'].props['posz'].set([0,1.7],[0,1])
        c['htail'].props['prpx'].set([0,0],[0,1])
        c['htail'].props['prpy'].set([0,0],[0,1])
        c['htail'].props['chord'].set([0.65,0.3],[0,1])

        self.computePoints()

if __name__ == '__main__':

    name = 'strutbraced'
    aircraft = Strutbraced()
    #aircraft.buildStructure(name)
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name+'_C')
    aircraft.oml0.write2IGES(name)
    #aircraft.oml0.write2EGADS(name+'_EGADS')
    #aircraft.plot()
