from __future__ import division
import numpy
from PAM.components import fullplate
from PAM.configurations import configuration


class wing(configuration):

    def __init__(self):
        wing = fullplate([8,8],[10])

        self.components = []
        self.components.append(wing)

        self.assembleComponents()


if __name__ == '__main__':

    aircraft = wing()
