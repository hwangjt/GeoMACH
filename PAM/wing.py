from __future__ import division
import numpy
import sys
sys.path.append(sys.path[0]+'/components/')
import fullplate, halfbody, fulljunction
import configuration


class wing(configuration.configuration):

    def __init__(self):
        wing = fullplate.fullplate([8,8],[10])

        self.components = []
        self.components.append(wing)

        self.assembleComponents()


if __name__ == '__main__':

    aircraft = wing()
