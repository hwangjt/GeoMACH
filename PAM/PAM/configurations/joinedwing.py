from __future__ import division
import numpy
import sys
from PAM.components import fullplate, halfbody, fulljunction, fuse_sections
from PAM.configurations import configuration


class joined(configuration):

    def __init__(self):
        fuse = halfbody([50,15,40,15,40,10,40,10,10],[30,15,15,20,20,10],[15])
        fuse.translatePoints(0,0,0)
        wing = fullplate([30,30],[40])
        wing.translatePoints(4,0,0)
        wingfuse = fulljunction(wing, 0, fuse, 2, [3,1], [4,3])
        tailfuse = fulljunction(wing, 1, fuse, 2, [2,7], [1,5])
        #wingfuse = fulljunction(fuse, 2, 0, [3,1], [4,3], wing, 0, 0, 0)
        #tailfuse = fulljunction(fuse, 2, 0, [1,5], [2,7], tail, 0, 0, 0)

        self.components = []
        self.components.append(fuse)
        self.components.append(wing)
        self.components.append(wingfuse)
        self.components.append(tailfuse)

        self.assembleComponents()

        self.components[0].setSections(1,fuse_sections.rounded2)
        self.components[0].setSections(2,fuse_sections.rounded2)
        self.components[0].setSections(3,fuse_sections.rounded2)
        self.components[0].props['posx'].set([0,10],[0,1])
        self.components[0].props['posy'].set([0.3,0.5,0.5],[0,0.15,1],w=[1.0,0,0],d=[1,0,0])
        self.components[0].props['ry'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])
        self.components[0].props['rz'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])

        self.components[1].offset[:] = [3.75, 0.3, 0.5]
#        self.components[1].setAirfoil("rae2822.dat")
        self.components[1].props['posx'].set([0,4.2],[0,1])
        self.components[1].props['posy'].set([0,0,0.5,0.5],[0,0.45,0.55,1])
        self.components[1].props['posz'].set([0,5.5,5.5,0],[0,0.45,0.55,1])
        self.components[1].props['prpx'].set([1,1],[0,1])
        self.components[1].props['prpy'].set([0,0],[0,1])
        self.components[1].props['chord'].set([1.4,0.5,0.5,1],[0,0.45,0.55,1])

        self.computePoints()

if __name__ == '__main__':

    aircraft = joined()
    aircraft.oml0.write2Tec('joined')
    aircraft.oml0.write2TecC('joinedC')
    aircraft.plot()
