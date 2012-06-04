from __future__ import division
import numpy
import sys
from PAM.components import halfbody
import configuration


class fuse(configuration.configuration):

    def __init__(self):
        #fuse = halfbody.halfbody([12,6,10,6,12],[8,6,6,8],[10])
        fuse = halfbody([40,20,20,20,40],[30,20,20,30],[30])


        self.components = []
        self.components.append(fuse)

        self.assembleComponents()

        self.components[0].setSections(1,fuse_sections.rounded2)
        self.components[0].setSections(2,fuse_sections.rounded2)
        self.components[0].setSections(3,fuse_sections.rounded2)
        self.components[0].Ls = [0.5, 3.0, 0.5, 2.0, 0.5, 4, 0.5]
        self.components[0].y0 *= 0.5

        self.computePoints()


if __name__ == '__main__':

    aircraft = fuse()
    aircraft.oml0.write2Tec('fuse')
    aircraft.plot()
