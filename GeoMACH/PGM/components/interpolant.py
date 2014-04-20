from __future__ import division
import numpy

from GeoMACH.PGM.components import Component, Face



class Interpolant(Component):

    def __init__(self):
        super(Interpolant,self).__init__() 
        self.nP = 10

    def initializeFaces(self):
        self.ms = []
        self.ms.append(numpy.zeros(sum(self.ni),int))
        self.ms.append(numpy.zeros(sum(self.nj),int))
        self.ms.append(None)
        self.ns = []
        self.ns.append(numpy.zeros(sum(self.ni),int))
        self.ns.append(numpy.zeros(sum(self.nj),int))
        self.ns.append(None)
        self.faces['def'] = Face(0, 1, 2, sum(self.ni), sum(self.nj))
            
        self.si = numpy.zeros(4,int)
        self.sj = numpy.zeros(4,int)
        for k in range(4):
            self.si[k] = sum(self.ni[:k])
            self.sj[k] = sum(self.nj[:k])

    def declare_properties(self):
        super(Interpolant, self).declare_properties()

        props = self.properties
        props['scl'] = [1,1]
        props['fC1'] = [1,1]
        props['mC1'] = [1,1]

    def initialize_properties(self, prop_vec, prop_index_vec):
        super(Interpolant, self).initialize_properties(prop_vec, prop_index_vec)

        add = self.addParam
        add('scl','scl',(1,1),P=[0.15])
        add('fC1','fC1',(1,1),P=[1.0])
        add('mC1','mC1',(1,1),P=[1.0])

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

    def getEdge2(self, Q, i=None, j=None, d=1):
        if j==None:
            if i==0:
                P = Q[:2,:]
            elif i==1:
                P = Q[1:3,:]
            elif i==-1:
                P = Q[-1:-3:-1,:]
            elif i==-2:
                P = Q[-2:-4:-1,:]
            P = numpy.swapaxes(P,0,1)
        else:
            if j==0:
                P = Q[:,:2]
            elif j==1:
                P = Q[:,1:3]
            elif j==-1:
                P = Q[:,-1:-3:-1]
            elif j==-2:
                P = Q[:,-2:-4:-1]
        if d==1:
            return P
        elif d==-1:
            return P[::-1,:]

    def setDerivatives(self, var, dV0):
        self.properties[var] += dV0
        self.computeQs()
        self.properties[var] -= dV0
