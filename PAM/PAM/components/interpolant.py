from __future__ import division
from PAM.components import Component
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib



class Interpolant(Component):

    def __init__(self):
        super(Interpolant,self).__init__() 
        self.nP = 10

    def initializeFaces(self):
        self.ms = []
        self.ms.append(numpy.zeros(sum(self.ni),int))
        self.ms.append(numpy.zeros(sum(self.nj),int))
        self.ms.append(None)
        self.faces.append([1,2])

        self.si = numpy.zeros(4,int)
        self.sj = numpy.zeros(4,int)
        for k in range(4):
            self.si[k] = sum(self.ni[:k])
            self.sj[k] = sum(self.nj[:k])

    def setDOFs(self):
        self.setC1('surf', 0, val=True)

    def initializeVariables(self):
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        self.variables = {
            'scale':numpy.array([0.15]),
            'f0':numpy.array([1.0]),
            'm0':numpy.array([1.0]),
            'shape':numpy.zeros((ni,nj),order='F')
            }

    def getEdge(self, Q, i=None, j=None, d=1):
        if j==None:
            if i==0:
                P = Q[:2,:,:]
            elif i==1:
                P = Q[1:3,:,:]
            elif i==-1:
                P = Q[-1:-3:-1,:,:]
            elif i==-2:
                P = Q[-2:-4:-1,:,:]
            P = numpy.swapaxes(P,0,1)
        else:
            if j==0:
                P = Q[:,:2,:]
            elif j==1:
                P = Q[:,1:3,:]
            elif j==-1:
                P = Q[:,-1:-3:-1,:]
            elif j==-2:
                P = Q[:,-2:-4:-1,:]
        if d==1:
            return P
        elif d==-1:
            return P[::-1,:,:]

    def setDerivatives(self, var, ind):
        self.variables[var][ind] += 1
        self.computeQs()
        self.variables[var][ind] -= 1
