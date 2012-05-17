from __future__ import division
import numpy
import sys
sys.path.append(sys.path[0]+'/components/')
import fullplate, halfbody, fulljunction, fuse_sections
import configuration


class wingbodytail(configuration.configuration):

    def __init__(self):
        fuse = halfbody.halfbody([40,20,20,20,20,15,15,15,10],[30,15,15,20,20,10],[15])
        fuse.translatePoints(0,0,0)
        wing = fullplate.fullplate([30,30],[20])
        wing.translatePoints(4,0,0)
        tail = fullplate.fullplate([30],[15])
        tail.translatePoints(8,0,0)
        wingfuse = fulljunction.fulljunction(fuse, 2, 0, [3,1], [4,3], wing, 0, 0, 0)
        tailfuse = fulljunction.fulljunction(fuse, 2, 0, [1,5], [2,7], tail, 0, 0, 0)

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

        self.components[1].offset[0] = 3.75
        self.components[1].offset[1] += 0.45
        self.components[1].offset[2] += 0.7
        self.components[1].chord[:] = numpy.linspace(2,0.3,self.components[1].chord.shape[0])
        self.components[1].sweep[:] = numpy.linspace(0,3.75,self.components[1].sweep.shape[0])
        self.components[1].span = 6
        self.components[1].thickness = 0.3

        self.components[2].offset[0] = 8.8
        self.components[2].offset[1] += 0.75
        self.components[2].offset[2] += 0.65
        self.components[2].chord[:] = numpy.linspace(1,0.15,self.components[2].chord.shape[0])
        self.components[2].sweep[:] = numpy.linspace(0,1.6,self.components[2].sweep.shape[0])
        self.components[2].span = 1.7
        self.components[2].thickness = 0.1

        self.computePoints()


if __name__ == '__main__':

    aircraft = wingbodytail()
    aircraft.oml0.write2Tec('wingbodytail')
    aircraft.oml0.write2TecC('wingbodytailC')
    aircraft.plot()
