from __future__ import division
from PAM.components import Component
import numpy, pylab, time, scipy.sparse
import PAM.PAMlib as PAMlib



class Primitive(Component):

    def __init__(self, nx, ny, nz):
        super(Primitive,self).__init__() 

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(numpy.zeros(ny,int))
        self.ms.append(numpy.zeros(nz,int))

    def initializeVariables(self):
        n = self.Qs[0].shape[1]
        zeros = numpy.zeros
        v = self.variables
        a = self.addParam

        v['scl'] = zeros((n,3),order='F')
        v['pos'] = zeros((n,3),order='F')
        v['rot'] = zeros((n,3),order='F')
        v['ogn'] = zeros((n,3),order='F')
        v['nor'] = zeros((n,3),order='F')
        v['flt'] = zeros((n,4),order='F')

        a('scl','scl',(1,1),P=[1.0])
        a('pos','pos',(2,3),P=[[0.,0.,0.],[1.,1.,1.]])
        a('rot','rot',(1,1),P=[0.0])
        a('ogn','ogn',(1,3),P=[0.,0.,0.])
        a('nor','nor',(1,1),P=[1.0])
        a('flt','flt',(1,1),P=[0.0])

    def setSections(self, sections=[], t1U=0, t2U=0, t1L=0, t2L=0):
        Ns = self.Ns
        v = self.variables
        for f in range(len(Ns)):
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    val = Ns[f][i,j,3]
                    if not val == -1:
                        break
                found = False
                for k in range(len(sections)):
                    found = found or (val==sections[k])
                if found or sections==[]:
                    v['flt'][j,0] = t1U
                    v['flt'][j,1] = t2U
                    v['flt'][j,2] = t1L
                    v['flt'][j,3] = t2L

    def computeSections(self, nQ, shapes, radii=None):
        nf = len(self.Qs)
        n = self.Qs[0].shape[1]
        v = self.variables

        rot0, Da, Di, Dj = PAMlib.computerotations(self.ax1, self.ax2, n, 9*(n*3-2), v['pos'], v['nor'])
        rot = v['rot']*numpy.pi/180.0 + rot0
        dv_dpos0 = scipy.sparse.csc_matrix((Da,(Di+6*n,Dj+3*n)),shape=(nQ,nQ))

        self.dQs_dv = range(nf)
        counter = 0
        for f in range(nf):
            if radii==None:
                scale = v['scl']
            else:
                scale = radii[f]
            ni, nj = self.Qs[f].shape[:2]
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computesections(self.ax1, self.ax2, ni, nj, ni*nj*27, counter, v['ogn'], scale, v['pos'], rot, shapes[f])
            self.dQs_dv[f] = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nQ))
            self.dQs_dv[f] = self.dQs_dv[f] + self.dQs_dv[f].dot(dv_dpos0)
            counter += 3*ni*nj

    def setDerivatives(self, var, dV0):
        nf = len(self.Qs)
        nv = self.dQs_dv[0].shape[1]
        n = self.Qs[0].shape[1]
        dV = numpy.zeros(nv)
        if var=='scl':
            dV[:3*n] = dV0.T.flatten()
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].dot(dV))
        elif var=='pos':
            dV[3*n:6*n] = dV0.T.flatten()
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                for i in range(ni):
                    self.Qs[f][i,:,:] += dV0
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].dot(dV))
        elif var=='rot':
            dV[6*n:9*n] = dV0.T.flatten()
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].dot(dV)*numpy.pi/180.0)
        else:
            self.variables[var] += dV0
            self.computeQs()
