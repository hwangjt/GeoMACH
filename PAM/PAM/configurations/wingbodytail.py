from __future__ import division
import numpy
import sys
from PAM.components import fullplate, halfbody, fulljunction, fuse_sections
from PAM.configurations import configuration


class wingbodytail(configuration):

    def __init__(self):
        fuse = halfbody([40,20,40,20,20,15,15,15,10],[30,15,15,20,20,10],[15])
        fuse.translatePoints(0,0,0)
        wing = fullplate([30,30],[40])
        wing.translatePoints(4,0,0)
        tail = fullplate([30],[15])
        tail.translatePoints(8,0,0)
        wingfuse = fulljunction(fuse, 2, 0, [3,1], [4,3], wing, 0, 0, 0)
        tailfuse = fulljunction(fuse, 2, 0, [1,5], [2,7], tail, 0, 0, 0)

        self.components = []
        self.components.append(fuse)
        self.components.append(wing)
        self.components.append(tail)
        self.components.append(wingfuse)
        self.components.append(tailfuse)

        self.assembleComponents()

        self.components[0].setSections(1,fuse_sections.rounded2)
        self.components[0].setSections(2,fuse_sections.rounded2)
        self.components[0].setSections(3,fuse_sections.rounded2)
        self.components[0].Ls = [0.1, 2.75, 0.5, 2.5, 0.5, 2, 0.3, 1.35, 0.3, 0.5, 0.05]

        self.components[0].propagateQs()
        self.components[0].setMain(0.7, 0.7, 0.7)

        n = 8
        y0 = [0.42, 0.49, 0.7]
        rz = [0.25, 0.42, 0.7]
        ry = [0.25, 0.42, 0.7]
        self.components[0].setNose(n,y0,rz,ry)

        n = 15
        y0 = [0.7, 0.7, 0.7]
        rz = [0.1, 0.15, 0.7]
        ry = [0.1, 0.15, 0.7]
        self.components[0].setTail(n,y0,rz,ry)

        self.components[1].offset[:] = [3.75, 0.45, 0.7]
        self.components[1].setAirfoil("rae2822.dat")
        self.components[1].props['posx'].set([0,3.2,4],[0,0.8,1],w=[0.4,1,0])
        self.components[1].props['posy'].set([0,0.5,1.7],[0,0.8,1],w=[1,1,0])
        self.components[1].props['posz'].set([0,5.5,6],[0,0.8,1],w=[0,1,0])
        self.components[1].props['prpx'].set([1,1],[0,1])
        self.components[1].props['prpy'].set([0,0,0,0],[0,0.4,0.8,1])
        self.components[1].props['chord'].set([2,0.25],[0,1])

        self.components[2].offset[:] = [9.1, 0.75, 0.55]
        self.components[2].props['posx'].set([0,1.6],[0,1],w=[0.2,0])
        self.components[2].props['posy'].set([0,0.3],[0,1],w=[0,0])
        self.components[2].props['posz'].set([0,1.7],[0,1])
        self.components[2].props['roty'].set([15,0],[0,1])
        self.components[2].props['prpx'].set([0,0],[0,1])
        self.components[2].props['prpy'].set([0,0],[0,1])
        self.components[2].props['chord'].set([1,0.15],[0,1])

        self.computePoints()

        self.components[1].computeCproj()
        self.components[2].computeCproj()

if __name__ == '__main__':

    aircraft = wingbodytail()
    aircraft.oml0.write2Tec('wingbodytail')
    aircraft.oml0.write2TecC('wingbodytailC')
    aircraft.plot()
