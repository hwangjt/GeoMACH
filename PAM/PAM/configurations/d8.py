from __future__ import division
import numpy
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration


class DoubleBubble(Configuration):

    def __init__(self):
        self.comps = {}
        self.keys = []

        self.addComp('fuse', Body([10,60,10,20,10,40,10,15,10],[20,10,10,10],[20]))
        self.addComp('wing', Wing([40],[20]))
        self.addComp('vtail', Wing([40],[15], opentip=True))
        self.addComp('htail', Wing([10,10,10,35],[10,10,15,15]))

        self.separateComps()

        self.addComp('wingfuse', FullInterface(self.comps, 'wing', 0, 'fuse', 2, [1,2], [2,4]))
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
        c['fuse'].props['rz'].set([0.3,0.8,0.8,0.8,0.7],[0,0.15,0.75,0.96,1.0],w=[0,0,0,0,1],d=[0,0,0,0,1])
        c['fuse'].props['tailL'] = 0.2

        c['wing'].offset[:] = [4.7, 0.3, 0.8]
        c['wing'].setAirfoil("rae2822.dat")
        c['wing'].props['posx'].set([0,0.25],[0,1])
        c['wing'].props['posy'].set([0,0.4],[0,1])
        c['wing'].props['posz'].set([0,5],[0,1])
        c['wing'].props['prpx'].set([0,1],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([0.7,0.2],[0,1],w=[1,0])

        c['vtail'].offset[:] = [8.8, 0.54, 0.61]
        c['vtail'].props['posx'].set([0,0.45],[0,1])
        c['vtail'].props['posy'].set([0,1.2],[0,1])
        c['vtail'].props['posz'].set([0,0],[0,1])
        c['vtail'].props['rotz'].set([5,0],[0,1])
        c['vtail'].props['prpx'].set([1,1],[0,1])
        c['vtail'].props['prpy'].set([0,0],[0,1])
        c['vtail'].props['chord'].set([0.5,0.25],[0,1])

        c['htail'].offset[:] = [9.0, 1.78, 0]
        c['htail'].props['posx'].set([0,0.5],[0,1])
        c['htail'].props['posy'].set([0,0],[0,1])
        c['htail'].props['posz'].set([0,2],[0,1])
        c['htail'].props['prpx'].set([0,0],[0,1])
        c['htail'].props['prpy'].set([0,0],[0,1])
        c['htail'].props['chord'].set([0.7,0.2],[0,1])

        self.computePoints()
  
        c['fuse'].addMembers('Longerons', 2, 1, 12, 15, A1=[0,0,0.95], C1=[1,0,1], A2=[0,1,0.95], C2=[1,1,1])
        c['fuse'].addMembers('Frames', 2, 2, 16, 11, A1=[0,0,0.85], C1=[0,1,1], A2=[1,0,0.85], C2=[1,1,1])

        c['wing'].addMembers('RibsLE', 1, 1, 13, 1, A1=[0,0,0], C1=[0,0.125,1], A2=[1,0,0], C2=[1,0.125,1])
        c['wing'].addMembers('Ribs', 1, 2, 13, 5, A1=[0,0.125,0], C1=[0,0.75,1], A2=[1,0.125,0], C2=[1,0.75,1])
        c['wing'].addMembers('RibsTE', 1, 1, 13, 1, A1=[0,0.75,0], C1=[0,1,1], A2=[1,0.75,0], C2=[1,1,1])
        c['wing'].addMembers('Spars', 1, 2, 2, 12, A1=[0,0.125,0], C1=[1,0.125,1], A2=[0,0.75,0], C2=[1,0.75,1])
        c['wing'].addMembers('Ustiff', 1, 1, 4, 12, A1=[0,0.25,0.9], C1=[1,0.25,1], A2=[0,0.625,0.9], C2=[1,0.625,1])
        c['wing'].addMembers('UstiffL', 1, 1, 4, 12, A1=[0,0.25,0.9], B1=[0,0.255,0.9], C1=[1,0.255,0.9], D1=[1,0.25,0.9], A2=[0,0.625,0.9], B2=[0,0.63,0.9], C2=[1,0.63,0.9], D2=[1,0.625,0.9])
        c['wing'].addMembers('Lstiff', 1, 1, 4, 12, A1=[0,0.25,0], C1=[1,0.25,0.1], A2=[0,0.625,0], C2=[1,0.625,0.1])
        c['wing'].addMembers('LstiffL', 1, 1, 4, 12, A1=[0,0.25,0.1], B1=[0,0.255,0.1], C1=[1,0.255,0.1], D1=[1,0.25,0.1], A2=[0,0.625,0.1], B2=[0,0.63,0.1], C2=[1,0.63,0.1], D2=[1,0.625,0.1])

        c['htail'].addMembers('RibsLE', 1, 1, 7, 1, A1=[0,0,0], C1=[0,0.125,1], A2=[1,0,0], C2=[1,0.125,1])
        c['htail'].addMembers('Ribs', 1, 2, 7, 5, A1=[0,0.125,0], C1=[0,0.75,1], A2=[1,0.125,0], C2=[1,0.75,1])
        c['htail'].addMembers('RibsTE', 1, 1, 7, 1, A1=[0,0.75,0], C1=[0,1,1], A2=[1,0.75,0], C2=[1,1,1])
        c['htail'].addMembers('Spars', 1, 2, 2, 6, A1=[0,0.125,0], C1=[1,0.125,1], A2=[0,0.75,0], C2=[1,0.75,1])
        c['htail'].addMembers('Ustiff', 1, 1, 4, 6, A1=[0,0.25,0.9], C1=[1,0.25,1], A2=[0,0.625,0.9], C2=[1,0.625,1])
        c['htail'].addMembers('UstiffL', 1, 1, 4, 6, A1=[0,0.25,0.9], B1=[0,0.255,0.9], C1=[1,0.255,0.9], D1=[1,0.25,0.9], A2=[0,0.625,0.9], B2=[0,0.63,0.9], C2=[1,0.63,0.9], D2=[1,0.625,0.9])
        c['htail'].addMembers('Lstiff', 1, 1, 4, 6, A1=[0,0.25,0], C1=[1,0.25,0.1], A2=[0,0.625,0], C2=[1,0.625,0.1])
        c['htail'].addMembers('LstiffL', 1, 1, 4, 6, A1=[0,0.25,0.1], B1=[0,0.255,0.1], C1=[1,0.255,0.1], D1=[1,0.25,0.1], A2=[0,0.625,0.1], B2=[0,0.63,0.1], C2=[1,0.63,0.1], D2=[1,0.625,0.1])

        c['vtail'].addMembers('Ribs', 1, 2, 13, 5, A1=[0,0.125,0], C1=[0,0.75,1], A2=[1,0.125,0], C2=[1,0.75,1])
        c['vtail'].addMembers('Spars', 1, 2, 2, 12, A1=[0,0.125,0], C1=[1,0.125,1], A2=[0,0.75,0], C2=[1,0.75,1])

if __name__ == '__main__':

    name = 'd8'
    aircraft = DoubleBubble()
    #aircraft.buildStructure(name)
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name+'_C')
    aircraft.oml0.write2IGES(name)
    #aircraft.oml0.write2EGADS(name+'_EGADS')
    #aircraft.plot()
