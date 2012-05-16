from __future__ import division
import sys
sys.path.append(sys.path[0]+'/../')
import PUBS
import numpy, pylab
import copy, time
import mpl_toolkits.mplot3d.axes3d as p3

oml1 = oml.oml()
oml1.importCGNSsurf('WBSolKW.cgns')
oml1.write2Tec('dlrf6')
oml1.plot(pylab.figure())
pylab.show()
