from __future__ import division
import numpy, pylab, time
import scipy.sparse
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
        self.nprim = len(self.comps)
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].translatePoints(0,0,k*4)   

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

#        a = PUBS.PUBSexport()
#        a.write2TecStruct('test.dat',Ps,['x','y','z'])
#        a.plot(Ps)
#        exit()
        self.oml0 = PUBS.PUBS(Ps)

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.oml0 = self.oml0
            comp.computems()
            comp.setDOFs()
        self.oml0.updateBsplines()

        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.initializeDOFmappings()
            comp.initializeVariables()
        self.computePoints()

    def updateComponents(self):
        self.oml0.updateBsplines()
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            comp.computeDims(self)
            comp.initializeDOFs()
        self.computePoints()

    def computePoints(self):
        self.computeQs()
        self.propagateQs()
        self.oml0.computePoints()

    def computeQs(self, full=True):
        if full:
            for k in range(len(self.comps)):
                self.comps[self.keys[k]].computeQs()
        else:
            for k in range(self.nprim,len(self.comps)):
                self.comps[self.keys[k]].computeQs()

    def propagateQs(self):
        self.oml0.Q[:,:3] = 0.0
        for k in range(len(self.comps)):
            self.comps[self.keys[k]].propagateQs()

    def getDerivatives(self, comp, var, ind, clean=True, FD=False, h=1e-5):
        self.computeQs()
        self.propagateQs()
        Q0 = numpy.array(self.oml0.Q[:,:3])
        if FD:
            self.comps[comp].variables[var][ind] += h
            self.computeQs()
            self.comps[comp].variables[var][ind] -= h
        else:
            self.comps[comp].setDerivatives(var,ind)
            self.computeQs(False)
            h = 1.0
        self.propagateQs()
        res = (self.oml0.Q[:,:3] - Q0)/h
        if clean:
            self.computePoints()
        return res

    def runDerivativeTest(self, comp, variables=[]):
        self.computePoints()
        if variables==[]:
            variables = self.comps[comp].variables.keys()
        h = 1e-5
        for var in variables:
            dat = self.comps[comp].variables[var]
            for ind,x in numpy.ndenumerate(dat):
                ind = ind[0] if len(ind)==1 else ind
                d1 = self.getDerivatives(comp,var,ind,clean=False)
                d2 = self.getDerivatives(comp,var,ind,clean=False,FD=True,h=h)
                norm0 = numpy.linalg.norm(d2)
                norm0 = 1.0 if norm0==0 else norm0
                error = numpy.linalg.norm(d2-d1)/norm0
                good = 'O' if error < 1e-4 else 'X'
                print good, ' ', comp, ' ', var, ' ', ind, ' ', error
        self.computePoints()


class Configuration2(object):
    
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

        self.oml0 = PUBS.PUBS(Ps)
        self.export = PUBS.PUBSexport(self.oml0)

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
        self.oml0.computePoints()

    def buildStructure(self):
        ABMs = []        
        Ss = []
        for k in range(len(self.comps)):
            comp = self.comps[self.keys[k]]
            if not comp.keys==[]:
                print 'Building structure for',self.keys[k]
                comp.buildStructure()     
                ABMs.append(comp.strABM)
                Ss.append(comp.strS)
        self.ABM = self.oml0.vstackSparse(ABMs)
        self.S = numpy.vstack(Ss)
        self.computePoints()

    def writeStructure(self, name):
        S = self.S
        P = self.ABM.dot(self.oml0.Q)
        f = open(name+'_str.dat','w')
        f.write('title = "PUBSlib output"\n')
        f.write('variables = "x", "y", "z"\n')
        iP = 0
        for surf in range(S.shape[0]):    
            nu = int(S[surf,0])
            nv = int(S[surf,1])
            f.write('zone i='+str(nu)+', j='+str(nv)+', DATAPACKING=POINT\n')
            for v in range(nv):
                for u in range(nu):
                    f.write(str(P[iP,0]) + ' ' + str(P[iP,1]) + ' ' + str(P[iP,2]) + '\n')
                    iP += 1
        f.close()

    def plot(self):
        #self.oml0.plot(pylab.figure(),False)
        #pylab.show()
        self.oml0.plotm(mlab.figure(),False)
        mlab.show()



