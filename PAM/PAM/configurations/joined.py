from __future__ import division
import numpy
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration


class Joined(Configuration):

    def __init__(self):
        self.comps = {}
        self.keys = []

        self.addComp('fuse', Body([70,10,25,10,40,10,10,10,25,10,20],[25,25,25,25],[15]))
        self.addComp('wing', Wing([80],[25],opentip=True))
        self.addComp('nacelle', Body([50,10,25,30],[20,10,10,20],[30],full=True))
        self.addComp('pylon', Wing([30],[25],opentip=True))
        self.addComp('fin', Wing([30,10,10,30],[10,10,25,10],half=True))

        self.separateComps()

        self.addComp('wingfuse', FullInterface(self.comps, 'wing', 0, 'fuse', 2, [2,1], [3,3]))
        self.addComp('pylonfuse', FullInterface(self.comps, 'pylon', 0, 'fuse', 2, [1,7], [2,9]))
        self.addComp('pylonnacelle', FullInterface(self.comps, 'pylon', 1, 'nacelle', 5, [2,3], [1,1]))
        self.addComp('finfuse', HalfInterface(self.comps, 'fin', 0, 'fuse', 1, [0,5], [0,10]))
        self.addComp('wingfin', FullInterface(self.comps, 'wing', 1, 'fin', 0, [3,1], [1,2]))

        self.assembleComponents()

        self.comps['wingfin'].setC1('surf', 0, j=0, v=0, val=False)
        self.comps['wingfin'].setC1('edge', 0, j=0, v=0)
        self.comps['wingfin'].extraDOFs.append([None,0])
        self.updateComponents()

        c = self.comps
        c['fuse'].setSections(sections=[1,2,3], t1L=0.35, t2L=0.65)
        c['fuse'].props['posx'].set([0,10],[0,1])
        c['fuse'].props['posy'].set([0.3,0.5,0.5],[0,0.15,1],w=[1.0,0,0],d=[1,0,0])
        c['fuse'].props['ry'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])
        c['fuse'].props['rz'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])

        c['wing'].offset[:] = [3.75, 0.3, 0.5]
        c['wing'].props['posx'].set([0,2.4,2.9,5.35],[0,0.45,0.55,1])
        c['wing'].props['posy'].set([0,0.3,1.0,1.19],[0,0.45,0.55,1])
        c['wing'].props['posz'].set([0,5.5,5.5,-0.45],[0,0.45,0.55,1])
        c['wing'].props['prpx'].set([1,1],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([1.4,0.5,0.5,0.8],[0,0.45,0.55,1])

        e = numpy.zeros((11,4))
        e[0,:] = [0.20, 0.03, 1, 1]
        e[1,:] = [0.30, 0.10, 0, 0]
        e[2,:] = [0.30, 0.29, 0.5, 0]
        e[3,:] = [0.00, 0.29, 1, 1]
        e[4,:] = [0.00, 0.31, 1, 1]
        e[5,:] = [0.50, 0.35, 1, 0]
        e[6,:] = [1.00, 0.30, 0, 0]
        e[7,:] = [1.00, 0.27, 0, 0]
        e[8,:] = [0.65, 0.27, 0, 0]
        e[9,:] = [0.65, 0.20, 1, 0]
        e[10,:] = [1.2, 0.05, 0, 0]
        e[:,:2] *= 1
        l = numpy.linspace(0,1,e.shape[0])
        c['nacelle'].offset[:] = [7.6, 0.65, 0.87]
        c['nacelle'].props['posx'].set(e[:,0],l)
        c['nacelle'].props['posy'].set([0,0],[0,1])
        c['nacelle'].props['ry'].set(e[:,1],l,e[:,2],e[:,3])
        c['nacelle'].props['rz'].set(e[:,1],l,e[:,2],e[:,3])

        c['pylon'].offset[:] = [7.9, 0.6, 0.43]
        c['pylon'].setAirfoil("naca0010")
        c['pylon'].props['posx'].set([0,0],[0,1])
        c['pylon'].props['posy'].set([0,0.05],[0,1])
        c['pylon'].props['posz'].set([0,0.1],[0,1])
        c['pylon'].props['roty'].set([10,0],[0,1])
        c['pylon'].props['prpx'].set([1,1],[0,1])
        c['pylon'].props['prpy'].set([0,0],[0,1])
        c['pylon'].props['chord'].set([0.75,0.75],[0,1])

        c['fin'].offset[:] = [8.4, 0.88, 0]
        c['fin'].props['posx'].set([0,1.5],[0,1])
        c['fin'].props['posy'].set([0,1.4],[0,1])
        c['fin'].props['posz'].set([0,0],[0,1])
        c['fin'].props['rotz'].set([10,0],[0,1])
        c['fin'].props['prpx'].set([1,1],[0,1])
        c['fin'].props['prpy'].set([0,0],[0,1])
        c['fin'].props['chord'].set([1.2,0.65],[0,1])

        self.computePoints()

if __name__ == '__main__':

    name = 'joined'
    aircraft = Joined()
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name+'_C')
    aircraft.oml0.write2IGES(name)
    aircraft.oml0.write2EGADS(name+'_EGADS')
    aircraft.plot()
