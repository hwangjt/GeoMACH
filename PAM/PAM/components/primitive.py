from __future__ import division
from PAM.components import Component, Variable
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
        v = self.variables
        v['scl'] = Variable((n,3))
        v['pos'] = Variable((n,3))
        v['rot'] = Variable((n,3))
        v['ogn'] = Variable((n,3),False)
        v['nor'] = Variable((n,3),False)

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
                    v['fillet'].V[j,0] = t1U
                    v['fillet'].V[j,1] = t2U
                    v['fillet'].V[j,2] = t1L
                    v['fillet'].V[j,3] = t2L

    def computeRotations(self):
        val = lambda var: self.variables[var]()
        n = self.Qs[0].shape[1]

        rot0, Da, Di, Dj = PAMlib.computerotations(self.ax1, self.ax2, n, 9*(n*3-2), val('pos'), val('nor'))
        drot0_dpos = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*n,3*n))
        rot = val('rot')*numpy.pi/180.0 + rot0
        return rot, drot0_dpos

    def computeSections(self, nQ, shapes, radii=None):
        val = lambda var: self.variables[var]()
        nf = len(self.Qs)
        n = self.Qs[0].shape[1]

        rot0, Da, Di, Dj = PAMlib.computerotations(self.ax1, self.ax2, n, 9*(n*3-2), val('pos'), val('nor'))
        rot = val('rot')*numpy.pi/180.0 + rot0
        dv_dpos0 = scipy.sparse.csc_matrix((Da,(Di+6*n,Dj+3*n)),shape=(nQ,nQ))

        self.dQs_dv = range(nf)
        counter = 0
        for f in range(nf):
            if radii==None:
                scale = val('scl')
            else:
                scale = radii[f]
            ni, nj = self.Qs[f].shape[:2]
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computesections(self.ax1, self.ax2, ni, nj, ni*nj*27, counter, val('ogn'), scale, val('pos'), rot, shapes[f])
            self.dQs_dv[f] = scipy.sparse.csc_matrix((Da,(Di,Dj)),shape=(3*ni*nj,nQ))
            self.dQs_dv[f] = self.dQs_dv[f] + self.dQs_dv[f].dot(dv_dpos0)
            counter += 3*ni*nj

    def setDerivatives(self, var, ind):
        nf = len(self.Qs)
        n = self.Qs[0].shape[1]
        if var=='scl':
            j,k = ind[:2]
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(nj*k+j).todense())
        elif var=='pos':
            j,k = ind[:2]
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,j,k] += 1.0
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(3*nj+nj*k+j).todense())
        elif var=='rot':
            j,k = ind[:2]
            for f in range(nf):
                ni, nj = self.Qs[f].shape[:2]
                self.Qs[f][:,:,:] += PAMlib.inflatevector(ni, nj, 3*ni*nj, self.dQs_dv[f].getcol(6*nj+nj*k+j).todense()*numpy.pi/180.0)
        else:
            self.variables[var][ind] += 1
            self.computeQs()
            self.variables[var][ind] -= 1
            
        
        
