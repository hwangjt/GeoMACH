from __future__ import division
import numpy

from GeoMACH.PGM import PGMlib
from GeoMACH.PGM.components import Primitive, Property


class Body(Primitive):
    """ A component used to model blunt bodies. """

    def __init__(self, nx=1, ny=1, nz=1):
        """ Initialization method
        nx: integer
            Number of surfaces in x direction
        ny: integer
            Number of surfaces in y direction
        nz: integer
            Number of surfaces in z direction
        """ 

        super(Body,self).__init__(nx,ny,nz)

        self.addFace('rgt', 2, 1, -0.5)
        self.addFace('top', 3, 1, 0.5)
        self.addFace('lft', -2, 1, 0.5)
        self.addFace('bot', -3, 1, -0.5)

        self.ax1 = 3
        self.ax2 = 1

    def setDOFs(self):
        faces = self.faces
        for face in self.faces.values():
            face.setC1('surf', val=True)

    def declare_properties(self):
        super(Body, self).declare_properties()

    def computeQs(self):
        faces = self.faces
        nx = faces['rgt'].num_cp[1]
        ny = faces['rgt'].num_cp[0]
        nz = faces['top'].num_cp[0]

        theta1 = {'rgt': -1/4.0, 
                  'top': 1/4.0,
                  'lft': 3/4.0,
                  'bot': 5/4.0,
                  }
        theta2 = {'rgt': 1/4.0, 
                  'top': 3/4.0,
                  'lft': 5/4.0,
                  'bot': 7/4.0,
                  }

        flt = self.props['flt'].prop_vec
        for name in self.faces:
            ni, nj = self.faces[name].num_cp
            self.shapes[name][:,:,:] = \
                PGMlib.computeshape(ni, nj, theta1[name], theta2[name],
                                    numpy.ones((nj,3),order='F'), 
                                    flt, numpy.zeros((ni,nj),order='F'))

        return self.computeSections()
