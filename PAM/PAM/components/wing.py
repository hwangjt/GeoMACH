from __future__ import division
from PAM.components import Component, Property, airfoils
import numpy, pylab, time, scipy.sparse
import mpl_toolkits.mplot3d.axes3d as p3
import PAM.PAMlib as PAMlib

import PUBS


class Wing(Component):
    """ A component used to model lifting surfaces. """

    def __init__(self, nx=1, nz=1, left=2, right=2):
        """ Initialization method
        nx: integer
            Number of surfaces in x (chord-wise) direction
        nz: integer
            Number of surfaces in z (span-wise) direction
        left, right: integer
            The v[0] and v[-1] sections of the wing
            0: open tip, C0
            1: open tip, C1
            2: closed tip
        """ 

        super(Wing,self).__init__() 

        self.ms = []
        self.ms.append(numpy.zeros(nx,int))
        self.ms.append(None)
        self.ms.append(numpy.zeros(nz,int))

        self.addFace(-1, 3, 1, 1, 1, left==2, right==2)
        self.addFace( 1, 3,-1, 1, 1, left==2, right==2)

        self.left = left
        self.right = right

    def setDOFs(self):
        left = self.left
        right = self.right

        for f in range(2):
            self.setC1('surf', f, val=True) #C1 Everywhere
            self.setC1('surf', f, i=-f, u=-f, val=False) #C0 trailing edge
            self.setC1('edge', f, i=-f, u=-f, val=True) #C0 trailing edge
            if left==0:                
                self.setC1('surf', f, j=0, v=0, val=False) #C0 left edge
                self.setC1('edge', f, j=0, v=0, val=True) #C0 left edge
                self.setCornerC1(f, i=-f, j=0, val=False) #C0 left TE corner
            if right==0:
                self.setC1('surf', f, j=-1, v=-1, val=False) #C0 right edge
                self.setC1('edge', f, j=-1, v=-1, val=True) #C0 right edge
                self.setCornerC1(f, i=-f, j=-1, val=False) #C0 right TE corner

    def initializeVariables(self):
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        zeros = numpy.zeros
        ones = numpy.ones
        self.variables = {
            'offset':zeros(3),
            'chord':ones(nj),
            'pos':zeros((nj,3),order='F'),
            'rot':zeros((nj,3),order='F'),
            'nor':zeros((nj,3),order='F'),
            'shape':zeros((2,ni,nj,3),order='F')
            }
        
    def propagateQs(self):
        r = numpy.array([0.25,0,0])
        ni = self.Qs[0].shape[0]
        nj = self.Qs[0].shape[1]
        v = self.variables
        nQ = self.oml0.nQ
        dQ_dv = scipy.sparse.csr_matrix((nQ*3,nj*(5+6*ni)))
        for f in range(2):
            rot0, Da, Di, Dj = PAMlib.computewingrotations(nj, 9*(nj*3-2), v['pos'])
            drot0_dpos = scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(nj*3,nj*3))
            rot = v['rot']*numpy.pi/180.0 + rot0*v['nor']
            self.Qs[f][:,:,:], Da, Di, Dj = PAMlib.computewingsections(f, ni, nj, ni*nj*24, nQ, r, v['offset'], v['chord'], v['pos'], rot, v['shape'][f,:,:,:], self.Ns[f][:,:,:])
            dQ_dv = dQ_dv + scipy.sparse.csr_matrix((Da,(Di,Dj)),shape=(nQ*3,nj*(5+6*ni)))

        dPdc0 = self.oml0.JM.dot(dQ_dv[:nQ,4*nj+2])
        dPdc1 = self.oml0.JM.dot(dQ_dv[nQ:2*nQ,4*nj+2])
        dPdc2 = self.oml0.JM.dot(dQ_dv[2*nQ:,4*nj+2])
        
        self.updateQs()
        self.oml0.computePoints()
        export = PUBS.PUBSexport(self.oml0)
        export.write2Tec('derTest',['chord0','chord1','chord2'],[dPdc0,dPdc1,dPdc2])
            
        if 1:
            #k = 5
            #k = nj+5
            k = 4*nj+5
            h = 1e-5
            Q0 = numpy.zeros((nQ,3))
            Q = numpy.zeros((nQ,3))
            for f in range(2):
                Qf, Da, Di, Dj = PAMlib.computewingsections(f, ni, nj, ni*nj*24, nQ, r, v['offset'], v['chord'], v['pos'], rot, v['shape'][f,:,:,:], self.Ns[f][:,:,:])
                for i in range(ni):
                    for j in range(nj):
                        if self.Ns[f][i,j,0] != -1:
                            Q0[self.Ns[f][i,j,0],:] = Qf[i,j,:]
            #v['chord'][k] += h
            #v['pos'][5,0] += h
            rot[5,0] += h
            for f in range(2):
                Qf, Da, Di, Dj = PAMlib.computewingsections(f, ni, nj, ni*nj*24, nQ, r, v['offset'], v['chord'], v['pos'], rot, v['shape'][f,:,:,:], self.Ns[f][:,:,:])
                for i in range(ni):
                    for j in range(nj):
                        if self.Ns[f][i,j,0] != -1:
                            Q[self.Ns[f][i,j,0],:] = Qf[i,j,:]
            #v['chord'][k] -= h
            #v['pos'][5,0] -= h
            rot[5,0] -= h
            print (Q-Q0)[:,0]/h
            print dQ_dv[:nQ,k].todense().T
            print numpy.linalg.norm((Q-Q0)[:,0]/h - dQ_dv[:nQ,k].todense().T)
            exit()
                        

    def setAirfoil(self,filename):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.variables['shape'][f,:,j,:2] = Ps[f][:,:]


if __name__ == '__main__':
    h = 1e-5
    P = numpy.array([0.1,0.2])
    t0,dt_dP = PAMlib.arctan2pi(P)
    for k in range(2):
        P[k] += h
        t,dt_dP0 = PAMlib.arctan2pi(P)
        P[k] -= h
        print (t-t0)/h
        print dt_dP[k]
        print '-------'

    rot = numpy.array([0.1,0.2,0.3])
    T0,dT_drot = PAMlib.computertnmtx(rot)
    for k in range(3):
        rot[k] += h
        T,dT_drot0 = PAMlib.computertnmtx(rot)
        rot[k] -= h
        print (T-T0)/h
        print dT_drot[:,:,k]
        print (T-T0)/h - dT_drot[:,:,k]
        print '------'
    #exit()

    w = Wing(nx=2,nz=2)#,left=0)
    import PUBS
    from mayavi import mlab
    w.oml0 = PUBS.PUBS(w.Ps)
    w.setDOFs()
    w.oml0.updateBsplines()
    w.computems()
    w.initializeDOFmappings()
    w.initializeVariables()
    w.variables['pos'][:,2] = numpy.linspace(0,1,w.Qs[0].shape[1])
    for j in range(w.Qs[0].shape[1]):
        w.variables['shape'][0,:,j,0] = 1 - numpy.linspace(0,1,w.Qs[0].shape[0])
        w.variables['shape'][1,:,j,0] = numpy.linspace(0,1,w.Qs[0].shape[0])
    w.variables['pos'][:,0] = numpy.linspace(0,1,w.Qs[0].shape[1])
    w.variables['pos'][:,1] = numpy.linspace(0,1,w.Qs[0].shape[1])
    #w.variables['rot'][:,2] = 20
    w.variables['nor'][:,:] = 1.0
    w.setAirfoil("naca0012.dat")
    w.propagateQs()
    w.updateQs()
    w.oml0.computePoints()
    w.oml0.plot(pylab.figure(),False)
    pylab.show()


class Wing2(Component):

    def __init__(self, nb, nc, half=False, opentip=False):
        if half:
            self.faces = numpy.zeros((1,2),int)
            self.faces[0,:] = [1,2]
        else:
            self.faces = numpy.zeros((2,2),int)
            self.faces[0,:] = [1,2]
            self.faces[1,:] = [-1,2]

        Ps = []
        Ks = []

        if not half:
            P, K = self.createSurfaces(Ks, nc[::-1], nb, -1, 3, 0)
            for k in range(len(P)):
                for v in range(P[k].shape[1]):
                    for u in range(P[k].shape[0]):
                        if (opentip or P[k][u,v,2]!=1) and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                            P[k][u,v,1] = 1
            Ps.extend(P)
            Ks.append(K)

        P, K = self.createSurfaces(Ks, nc, nb, 1, 3, 0)
        for k in range(len(P)):
            for v in range(P[k].shape[1]):
                for u in range(P[k].shape[0]):
                    if (opentip or P[k][u,v,2]!=1) and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                        P[k][u,v,1] = -1
        Ps.extend(P)
        Ks.append(K)

        self.nb = nb
        self.nc = nc
        self.Ps = Ps
        self.Ks = Ks
        self.half = half
        self.opentip = opentip

        self.oml0 = []

    def setDOFs(self):        
        half = self.half
        opentip = self.opentip
        nf = len(self.Ks)
        for f in range(nf):
            self.setC1('surf', f)
            self.setC1('surf', f, j=0, v=0, val=False)
            self.setC1('surf', f, i=-f, u=-f, val=False)
            if opentip or half:
                self.setC1('surf', f, j=-1, v=-1, val=False)
            if half:
                self.setC1('surf', f, i=f-1, u=f-1, val=False)
        for f in range(nf):
            self.setC1('edge', f, j=0, v=0)
            self.setC1('edge', f, i=-f, u=-f)
            if opentip or half:
                self.setC1('edge', f, j=-1, v=-1)
            if half:
                self.setC1('edge', f, i=f-1, u=f-1)
        for f in range(nf):
            self.setCornerC1(f, i=-f, j=0, val=False)
            if opentip or half:
                self.setCornerC1(f, i=-f, j=-1, val=False)
            if half:
                self.setCornerC1(f, i=f-1, j=0, val=False)
                self.setCornerC1(f, i=f-1, j=-1, val=False)

    def isExteriorDOF(self, f, uType, vType, i, j):
        check = self.check
        half = self.half
        opentip = self.opentip
        value = check(uType,vType,v=0) or check(uType,vType,u=-f) or check(uType,vType,u=-f,v=0)
        if opentip or half:
            value = value or check(uType,vType,v=-1) or check(uType,vType,u=-f,v=-1)
        if half:
            value = value or check(uType,vType,u=f-1) or check(uType,vType,u=f-1,v=0) or check(uType,vType,u=f-1,v=-1)
        return value

    def initializeParameters(self):
        Ns = self.Ns
        self.offset = numpy.zeros(3)
        self.SECTshape = numpy.zeros((len(self.Ks),Ns[0].shape[0],Ns[0].shape[1],3),order='F')
        self.SECTrot0 = numpy.zeros((Ns[0].shape[1],3),order='F')
        self.props = {
            'chord':Property(Ns[0].shape[1]),
            'posx':Property(Ns[0].shape[1]),
            'posy':Property(Ns[0].shape[1]),
            'posz':Property(Ns[0].shape[1]),
            'rotx':Property(Ns[0].shape[1]),
            'roty':Property(Ns[0].shape[1]),
            'rotz':Property(Ns[0].shape[1]),
            'prpx':Property(Ns[0].shape[1]),
            'prpy':Property(Ns[0].shape[1])
            }
        self.setAirfoil("naca0012.dat")

    def setAirfoil(self,filename):
        airfoil = airfoils.getAirfoil(filename)
        if self.half:
            airfoil[0][:,:] = airfoil[1][:,:]
        Ps = airfoils.fitAirfoil(self,airfoil,rev=self.half)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.SECTshape[f,:,j,:2] = Ps[f][:,:]
        if self.half:
            self.SECTshape[0,:,-1,1] = 0
        
    def propagateQs(self):
        a = 0.25
        b = 0.0
        Ns = self.Ns
        Qs = self.Qs
        for f in range(len(self.Ks)):
            p = self.props
            Qs[f][:,:,:] = PAMlib.computewingsections(Ns[f].shape[0], Ns[f].shape[1], a, b, p['posx'].data, p['posy'].data, p['posz'].data, p['rotx'].data, p['roty'].data, p['rotz'].data, p['prpx'].data, p['prpy'].data, p['chord'].data, self.SECTshape[f])
            Qs[f][:,:,0] += self.offset[0]
            Qs[f][:,:,1] += self.offset[1]
            Qs[f][:,:,2] += self.offset[2]
        if self.half:
            Qs[0][:,-1,2] = 0
            Qs[0][0,:,2] = 0
            Qs[0][-1,:,2] = 0

    def getFlattenedC(self, f, i, j, ni, nj):
        ii = i/(ni-1)
        jj = j/(nj-1)
        return [jj,self.SECTshape[f][i,j,0],0]
        #if f==0:
        #    return [jj,1-ii,0]
        #else:
        #    return [jj,ii,0]

    def getAR(self):
        return 5

    def getSkinIndices(self):
        return [[0],[1]]
