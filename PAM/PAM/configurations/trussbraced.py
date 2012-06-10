from __future__ import division
import numpy
import sys
from PAM.components import fullplate, halfbody, fulljunction, fuse_sections
from PAM.configurations import configuration


class trussbraced(configuration):

    def __init__(self):
        fuse = halfbody([70,10,10,20,10,10,10,50,10,25,10,10],[25,25,25,25,25,25],[15])
        fuse.translatePoints(0,0,0)
        wing = fullplate([30,10,10,50],[10,20,10,10])
        wing.translatePoints(4,0,0)
        tail = fullplate([30],[25])
        tail.translatePoints(8,0,0)
        strut = fullplate([30],[20],opentip=True)
        strut.translatePoints(12,0,0)
        wingfuse = fulljunction(wing, 0, fuse, 2, [0,1], [1,6])
        strutfuse = fulljunction(strut, 0, fuse, 2, [4,2], [5,4])
        strutwing = fulljunction(strut, 1, wing, 1, [2,1], [0,2])
        tailfuse = fulljunction(tail, 0, fuse, 2, [2,8], [3,10])

        self.components = []
        self.components.append(fuse)
        self.components.append(wing)
        self.components.append(strut)
        self.components.append(tail)
        self.components.append(wingfuse)
        self.components.append(strutfuse)
        self.components.append(strutwing)
        self.components.append(tailfuse)

        self.assembleComponents()

        self.components[0].setSections(2,fuse_sections.rounded4)
        self.components[0].setSections(3,fuse_sections.rounded4)
        self.components[0].setSections(4,fuse_sections.rounded4)
        self.components[0].props['posx'].set([0,10],[0,1])
        self.components[0].props['posy'].set([0.3,0.5,0.5],[0,0.15,1],w=[1.0,0,0],d=[1,0,0])
        self.components[0].props['ry'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])
        self.components[0].props['rz'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])

        self.components[1].offset[:] = [3.75, 0.85, 0.5]
        self.components[1].setAirfoil("rae2822.dat")
        self.components[1].props['posx'].set([0,0],[0,1])
        self.components[1].props['posy'].set([0,0],[0,1])
        self.components[1].props['posz'].set([0,8],[0,1])
        self.components[1].props['prpx'].set([0,0],[0,1])
        self.components[1].props['prpy'].set([0,0],[0,1])
        self.components[1].props['chord'].set([1,0.3],[0,1],w=[1,0])

        self.components[2].offset[:] = [3.75, 0.15, 0.5]
        self.components[2].props['posx'].set([0,0],[0,1])
        self.components[2].props['posy'].set([0,0.65],[0,1])
        self.components[2].props['posz'].set([0,3.2],[0,1])
        self.components[2].props['prpx'].set([1,1],[0,1])
        self.components[2].props['prpy'].set([0,0],[0,1])
        self.components[2].props['chord'].set([0.5,0.5],[0,1])

        self.components[3].offset[:] = [8.5, 0.5, 0.35]
        self.components[3].props['posx'].set([0,0.5],[0,1],w=[0.2,0])
        self.components[3].props['posy'].set([0,0.1],[0,1],w=[0,0])
        self.components[3].props['posz'].set([0,2],[0,1])
        self.components[3].props['roty'].set([10,0],[0,1])
        self.components[3].props['prpx'].set([0,0],[0,1])
        self.components[3].props['prpy'].set([0,0],[0,1])
        self.components[3].props['chord'].set([0.85,0.15],[0,1])

        self.computePoints()

if __name__ == '__main__':

    aircraft = trussbraced()
    aircraft.oml0.write2Tec('trussbraced')
    aircraft.oml0.write2TecC('trussbracedC')
    aircraft.plot()
