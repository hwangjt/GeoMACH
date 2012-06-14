from __future__ import division
import numpy
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration


class DoubleBubble(Configuration):

    def __init__(self):
        self.comps = {}
        self.keys = []

        self.addComp('fuse', Body([10,60,10,20,10,40,10,15,10],[10,10,30,10,10],[20]))
        self.addComp('wing', Wing([40],[20]))
        self.addComp('vtail', Wing([40],[15], opentip=True))
        self.addComp('htail', Wing([10,10,10,20],[10,10,15,10]))

        self.separateComps()

        self.addComp('wingfuse', FullInterface(self.comps, 'wing', 0, 'fuse', 2, [3,2], [4,4]))
        self.addComp('vtailfuse', FullInterface(self.comps, 'vtail', 0, 'fuse', 2, [0,6], [1,8]))
        self.addComp('vtailhtail', FullInterface(self.comps, 'vtail', 1, 'htail', 1, [3,1], [1,2]))

        self.assembleComponents()

        self.comps['vtailhtail'].setC1('surf', 0, j=0, v=0, val=False)
        self.comps['vtailhtail'].setC1('edge', 0, j=0, v=0)
        self.comps['vtailhtail'].extraDOFs.append([None,0])
        self.updateComponents()

        c = self.comps
        c['fuse'].setSections(sections=range(1,10), t1U=0, t2U=0.6, t1L=0, t2L=0.6)
        c['fuse'].props['posx'].set([0,10],[0,1])
        c['fuse'].props['posy'].set([0.7,0.4,0.4],[0,0.15,1])
        c['fuse'].props['ry'].set([0.1,0.4,0.4,0.1],[0,0.15,0.75,1.0],w=[0,0,0,0],d=[1,0,0,0])
        c['fuse'].props['rz'].set([0.3,0.8,0.8,0.8],[0,0.15,0.75,1.0])

        c['wing'].offset[:] = [4.7, 0.05, 0.8]
        c['wing'].setAirfoil("rae2822.dat")
        c['wing'].props['posx'].set([0,0.25],[0,1])
        c['wing'].props['posy'].set([0,0.4],[0,1])
        c['wing'].props['posz'].set([0,5],[0,1])
        c['wing'].props['prpx'].set([0,1],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([0.7,0.2],[0,1],w=[1,0])

        c['vtail'].offset[:] = [9.3, 0.53, 0.65]
        c['vtail'].props['posx'].set([0,0.45],[0,1])
        c['vtail'].props['posy'].set([0,1.2],[0,1])
        c['vtail'].props['posz'].set([0,0],[0,1])
        c['vtail'].props['rotz'].set([5,0],[0,1])
        c['vtail'].props['prpx'].set([1,1],[0,1])
        c['vtail'].props['prpy'].set([0,0],[0,1])
        c['vtail'].props['chord'].set([0.5,0.25],[0,1])

        c['htail'].offset[:] = [9.55, 1.8, 0]
        c['htail'].props['posx'].set([0,0.5],[0,1])
        c['htail'].props['posy'].set([0,0],[0,1])
        c['htail'].props['posz'].set([0,2],[0,1])
        c['htail'].props['prpx'].set([0,0],[0,1])
        c['htail'].props['prpy'].set([0,0],[0,1])
        c['htail'].props['chord'].set([0.7,0.2],[0,1])

        self.computePoints()

if __name__ == '__main__':

    name = 'd8'
    aircraft = DoubleBubble()
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name+'_C')
    aircraft.oml0.write2IGES(name)
    aircraft.plot()
