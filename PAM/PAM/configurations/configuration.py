from __future__ import division
import numpy, pylab, time
import PUBS
import mpl_toolkits.mplot3d.axes3d as p3
from mayavi import mlab


class Configuration(object):
    
    def __init__(self): 
        self.comps = {}
        self.keys = []

    def addComp(self, name, comp):
        self.comps[name] = comp
        self.keys.append(name)

    def separateComps(self):
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].translatePoints(k*4,0,0)   

    def assembleComponents(self):
        Ps = []        
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            if k > 0:
                comp0 = self.comps[self.keys[k-1]]
                maxk = numpy.max(comp0.Ks[-1]) + 1
                for s in range(len(comp.Ks)):
                    for j in range(comp.Ks[s].shape[1]):
                        for i in range(comp.Ks[s].shape[0]):
                            if comp.Ks[s][i,j] != -1:
                                comp.Ks[s][i,j] += maxk
            Ps.extend(comp.Ps)
            comp.Ps = []

        self.oml0 = PUBS.PUBS()
        self.oml0.importSurfaces(Ps)

        #self.oml0.plotm(mlab.figure())
        #self.oml0.write2Tec('test')
        #mlab.show()
        #exit()

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.oml0 = self.oml0
            comp.members = {}
            comp.keys = []
            comp.computeDims(self)
            comp.setDOFs()
        self.oml0.updateBsplines()

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.initializeDOFs()
            comp.initializeParameters()
            comp.propagateQs()
            comp.updateQs()
        self.oml0.computePoints()

    def updateComponents(self):
        self.oml0.updateBsplines()
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.computeDims(self)
            comp.initializeDOFs()
        self.computePoints()

    def computePoints(self):
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.propagateQs()
            comp.updateQs()
        t0 = time.time()
        self.oml0.computePoints()
        print '0',time.time()-t0

    def plot(self):
        #self.oml0.plot(pylab.figure(),False)
        #pylab.show()
        self.oml0.plotm(mlab.figure(),True)
        mlab.show()



