from __future__ import division
import numpy, pylab
import PUBS
import mpl_toolkits.mplot3d.axes3d as p3
from mayavi import mlab


class configuration(object):
    
    def __init__(self): 
        self.components = {}
        
    def __gettattr__(self,name):
        return self.components[name]
    
    def add(self,name,component): 
        self.components[name] = component    

    def assembleComponents(self):
        Ps = []
        for k in range(len(self.components)):
            if k > 0:
                maxk = numpy.max(self.components[k-1].Ks[-1]) + 1
                for s in range(len(self.components[k].Ks)):
                    for j in range(self.components[k].Ks[s].shape[1]):
                        for i in range(self.components[k].Ks[s].shape[0]):
                            if self.components[k].Ks[s][i,j] != -1:
                                self.components[k].Ks[s][i,j] += maxk
            Ps.extend(self.components[k].Ps)
            self.components[k].Ps = []

        self.oml0 = PUBS.PUBS()
        self.oml0.importSurfaces(Ps)
        #self.oml0.plotm(mlab.figure())
        #self.oml0.write2Tec('test')
        #mlab.show()
        #exit()

        for k in range(len(self.components)):
            self.components[k].oml0 = self.oml0
            self.components[k].setDOFs()
        self.oml0.updateBsplines()

        for k in range(len(self.components)):
            self.components[k].initializeDOFs()
            self.components[k].initializeParameters()
            self.components[k].propagateQs()
            self.components[k].updateQs()
        self.oml0.computePoints()

    def computePoints(self):
        for k in range(len(self.components)):
            self.components[k].propagateQs()
            self.components[k].updateQs()
        self.oml0.computePoints()

    def plot(self):
        #self.oml0.plot(pylab.figure(),False)
        #pylab.show()
        self.oml0.plotm(mlab.figure(),True)
        mlab.show()



