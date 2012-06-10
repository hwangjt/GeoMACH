from __future__ import division
import numpy
import sys
from PAM.components import fullplate, halfbody, fulljunction, fuse_sections
from PAM.configurations import configuration


class wingbody(configuration):

    def __init__(self):
#        fuse = halfbody([50,15,20,30,15,40,10,25,10,10],[30,15,15,20,20,10],[15])
        fuse = halfbody([50,15,20,30,15,40,10,25,10,10],[10,40,10],[15])
        fuse.translatePoints(0,0,0)
#        wing = fullplate([30,30],[20,30])
        wing = fullplate([30,30],[40])
        wing.translatePoints(4,0,0)
        wingfuse = fulljunction(wing, 0, fuse, 2, [2,2], [0,3])

        self.components = []
        self.components.append(fuse)
        self.components.append(wing)
        self.components.append(wingfuse)

        self.assembleComponents()

        self.components[0].setSections(1,fuse_sections.rounded2)
        self.components[0].setSections(2,fuse_sections.rounded2)
        self.components[0].setSections(3,fuse_sections.rounded2)
        self.components[0].props['posx'].set([0,10],[0,1])
        self.components[0].props['posy'].set([0.3,0.5,0.5],[0,0.15,1],w=[1.0,0,0],d=[1,0,0])
        self.components[0].props['ry'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])
        self.components[0].props['rz'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])

        self.components[1].offset[:] = [3.75, 0.3, 0.5]
        self.components[1].setAirfoil("rae2822.dat")
        self.components[1].props['posx'].set([0,3.2,4],[0,0.8,1],w=[0.4,1,0])
        self.components[1].props['posy'].set([0,0.5,1.7],[0,0.8,1],w=[1,1,0])
        self.components[1].props['posz'].set([0,4.5,5],[0,0.8,1],w=[0,1,0])
        self.components[1].props['prpx'].set([1,1],[0,1])
        self.components[1].props['prpy'].set([0,0],[0,1])
        self.components[1].props['rotz'].set([-90,-90],[0,1])
        self.components[1].props['chord'].set([2,0.25],[0,1])

        self.computePoints()

if __name__ == '__main__':

    aircraft = wingbody()
    aircraft.oml0.write2Tec('wingbody')
    aircraft.oml0.write2TecC('wingbodyC')
    aircraft.plot()
