from __future__ import division
import numpy
import sys
sys.path.append(sys.path[0]+'/components/')
import fullplate, halfbody, fulljunction, fuse_sections
import configuration



class wingbody(configuration.configuration):

    def __init__(self):
        #fuse = halfbody.halfbody([12,6,10,6,12],[8,6,6,8],[10])
        #fuse.translatePoints(0,0,0)
        #wing = fullplate.fullplate([8,8],[10])
        fuse = halfbody.halfbody([40,20,20,20,40],[60,20,20,10],[15])
        fuse.translatePoints(0,0,0)
        wing = fullplate.fullplate([30,30],[20])
        wing.translatePoints(4,0.5,2)
        wingfuse = fulljunction.fulljunction(fuse, 2, 0, [1,1], [2,3], wing, 0, 0, 0)

        self.components = []
        self.components.append(fuse)
        self.components.append(wing)
        self.components.append(wingfuse)

        self.assembleComponents()

        self.components[0].setSections(1,fuse_sections.rounded2)
        self.components[0].setSections(2,fuse_sections.rounded2)
        self.components[0].setSections(3,fuse_sections.rounded2)
        self.components[0].Ls = [0.5, 2.75, 0.5, 2.5, 0.5, 3.25, 0.5]
        self.components[0].y0 *= 1

        self.components[1].offset[0] = 3.75
        self.components[1].offset[1] += 0.6
        self.components[1].offset[2] += 1
        self.components[1].chord[:] = numpy.linspace(2,0.5,self.components[1].chord.shape[0])
        self.components[1].sweep[:] = numpy.linspace(0,2.5,self.components[1].sweep.shape[0])

        self.computePoints()


if __name__ == '__main__':

    aircraft = wingbody()
    aircraft.oml0.write2Tec('wingbody')
    aircraft.oml0.write2TecC('wingbodyC')
    aircraft.plot()
