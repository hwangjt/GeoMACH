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
        ones = numpy.ones
        self.variables = {
            'offset':zeros(3),
            'scale':ones((n,3)),
            'pos':zeros((n,3),order='F'),
            'rot':zeros((n,3),order='F'),
            }
        self.parameters = {
            'origin': zeros(3),
            'nor':ones((n,3),order='F'),
            }

    def setSections(self, sections=[], t1U=0, t2U=0, t1L=0, t2L=0):
        Ns = self.Ns
        v = self.variables
        p = self.parameters
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
                    p['fillet'][j,0] = t1U
                    p['fillet'][j,1] = t2U
                    p['fillet'][j,2] = t1L
                    p['fillet'][j,3] = t2L

    def computeRotations(self):
        n = self.Qs[0].shape[1]
        v = self.variables
        p = self.parameters

        rot0, Da, Di, Dj = PAMlib.computerotations(self.ax1, self.ax2, n, 9*(n*3-2), v['pos'], p['nor'])
        drot0_dpos = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*n,3*n))
        rot = v['rot']*numpy.pi/180.0 + rot0
        return rot, drot0_dpos

    def computeSections(self, nQ, rot, shapes, radii=None):
        nf = len(self.Qs)
        n = self.Qs[0].shape[1]
        v = self.variables
        p = self.parameters
        self.dQs_dv = range(nf)
        counter = 0
        for f in range(nf):
            if radii==None:
                scale = v['scale']
            else:
                scale = radii[f]
            ni, nj = self.Qs[f].shape[:2]
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computesections(self.ax1, self.ax2, ni, nj, ni*nj*27, counter, p['origin'], v['offset'], scale, v['pos'], rot, shapes[f])
            self.dQs_dv[f] = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nQ))
            counter += 3*ni*nj

    def setDerivatives(self, var, ind):
        nf = len(self.Qs)
        n = self.Qs[0].shape[1]
        if var=='offset':
            for f in range(nf):
                self.Qs[f][:,:,ind] += 1.0
        elif var=='scale':
            j,k = ind[:2]
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(nj*k+j).todense())
        elif var=='pos':
            j,k = ind[:2]
            A = scipy.sparse.csc_matrix((3*n,3*n))
            B = self.drot0_dpos
            C = scipy.sparse.csc_matrix((self.dQs_dv[0].shape[1]-6*n,3*n))
            D = scipy.sparse.vstack([A,B,C],format='csc')
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,j,k] += 1.0
                Q = self.dQs_dv[f].dot(D).getcol(nj*k+j).todense()
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, Q)
        elif var=='rot':
            j,k = ind[:2]
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(3*nj+nj*k+j).todense()*numpy.pi/180.0)
        else:
            self.variables[var][ind] += 1
            self.computeQs()
            self.variables[var][ind] -= 1
            
        
        
