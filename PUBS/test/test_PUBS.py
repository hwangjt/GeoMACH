from __future__ import division
import unittest
import numpy
import PUBS

class TestPUBS(unittest.TestCase):

    def __init__(self, *args):
        unittest.TestCase.__init__(self, *args)
        P0 = []
        P0.append(numpy.zeros((10,10,3)))
        P0.append(numpy.zeros((10,10,3)))
        
        for i in range(10):
            for j in range(10):
                P0[0][i,j,0] = i/9
                P0[0][i,j,1] = j/9
                P0[1][i,j,0] = 1 + i/9
                P0[1][i,j,1] = j/9

        self.P0 = P0   
        
    def test_topology(self):
        """ Check that the numbers of vertices, edges, 
        and surfaces are correct """

        oml1 = PUBS.PUBS()
        oml1.initializeTopology(self.P0)

        self.assertEqual(oml1.nsurf,2)
        self.assertEqual(oml1.nedge,7)
        self.assertEqual(oml1.nvert,6)
        
if __name__ == "__main__":
    unittest.main()
