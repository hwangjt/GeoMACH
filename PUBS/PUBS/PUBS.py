from __future__ import division
import PUBSlib
import export2EGADS
import numpy, scipy, pylab, copy, time
import math
import scipy.sparse, scipy.sparse.linalg
import mpl_toolkits.mplot3d.axes3d as p3
from mayavi import mlab


class PUBS(object):
    """
    A surface model represented as a watertight union of B-spline surfaces. 
    """

    def importSurfaces(self, P, ratio=3.0):
        self.symmPlane = 2
        self.initializeTopology(P)
        self.initializeBsplines(ratio)
        self.computeQindices()
        self.computeCindices()
        self.computeDindices()
        self.computePindices()
        self.initializePoints(P)
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()
        self.computeControlPts()
        self.computePoints()

    def importCGNSsurf(self, filename):
        n = PUBSlib.nsurfaces2(filename)  
        z,sizes = PUBSlib.surfacesizes2(filename, n) 
        P0 = []
        for i in range(n):
            P0.append(PUBSlib.getsurface2(filename,z[i],sizes[i,0],sizes[i,1]))
        self.importSurfaces(P0)

    def updateBsplines(self):
        self.computeQindices()
        self.computeCindices()
        self.computeDindices()
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()
        self.computeControlPts()
        self.computePoints()

    def updateEvaluation(self):
        self.computePindices()
        self.computeParameters()
        self.computeJacobian()
        self.computePoints()

    def initializeTopology(self, P):
        self.nsurf = len(P)
        print '# Surfaces =',self.nsurf
        Ps = numpy.zeros((self.nsurf,3,3,3),order='F')
        for k in range(self.nsurf):
            n = P[k].shape[0:2]
            for i in range(2):
                for j in range(2):
                    Ps[k,i*2,j*2] = P[k][-i,-j]
            left = 0.5*(P[k][0,int(numpy.ceil((n[1]-1)/2.0))] + P[k][0,int(numpy.floor((n[1]-1)/2.0))])
            right = 0.5*(P[k][-1,int(numpy.ceil((n[1]-1)/2.0))] + P[k][-1,int(numpy.floor((n[1]-1)/2.0))])
            bottom = 0.5*(P[k][int(numpy.ceil((n[0]-1)/2.0)),0] + P[k][int(numpy.floor((n[0]-1)/2.0)),0])
            top = 0.5*(P[k][int(numpy.ceil((n[0]-1)/2.0)),-1] + P[k][int(numpy.floor((n[0]-1)/2.0)),-1])
            Ps[k,0,1] = left
            Ps[k,2,1] = right
            Ps[k,1,0] = bottom
            Ps[k,1,2] = top
        self.ns = numpy.zeros((self.nsurf,2),order='F')
        for k in range(self.nsurf):
            self.ns[k,:] = P[k].shape[0:2]
        self.nvert,self.nedge,self.surf_vert,self.surf_edge = PUBSlib.computetopology(self.nsurf,1e-15,1e-5,Ps)
        print '# Vertices =',self.nvert
        print '# Edges =',self.nedge
        self.vert_count,self.edge_count = PUBSlib.countveptrs(self.nsurf,self.nvert,self.nedge,self.surf_vert,self.surf_edge)
        self.ngroup,self.edge_group = PUBSlib.computegroups(self.nsurf,self.nedge,self.surf_edge)
        print '# Groups =',self.ngroup
        self.surf_c1 = numpy.zeros((self.nsurf,3,3),bool,order='F')
        self.edge_c1 = numpy.zeros((self.nedge,2),bool,order='F')

    def initializeBsplines(self, ratio):
        k = 4
        self.group_k, self.group_m, self.group_n = PUBSlib.getkmn(k, self.nsurf, self.nedge, self.ngroup, ratio, self.ns, self.surf_edge, self.edge_group)

    def initializePoints(self, P):
        for k in range(self.nsurf):
            n1 = P[k].shape[0]
            n2 = P[k].shape[1]
            PUBSlib.populatep(self.nP, n1, n2, k+1, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, P[k], self.P)
        PUBSlib.avgboundaries(self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.edge_index_P, self.P)
        self.vert_symm, self.edge_symm = PUBSlib.determinesymm(self.symmPlane+1, 1e-10, 1e-8, self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.P)

    def computeQindices(self):
        self.surf_index_Q = PUBSlib.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_Q = PUBSlib.getedgeindicesq(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m, self.surf_c1)
        self.vert_index_Q = PUBSlib.getvertindicesq(self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_c1, self.edge_c1)
        self.nQ = 0
        self.nQ += max(self.vert_index_Q)
        self.nQ += max(self.edge_index_Q[:,1])
        self.nQ += self.surf_index_Q[-1,1]
        print '# Degrees of freedom =',self.nQ

    def computeCindices(self):
        self.surf_index_C = PUBSlib.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_C = PUBSlib.getedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_m)
        self.nD = 0
        for i in range(self.ngroup):
            self.nD += self.group_m[i] + self.group_k[i]
        self.nC = self.nvert
        self.nC += self.edge_index_C[-1,1]
        self.nC += self.surf_index_C[-1,1]
        self.C = numpy.zeros((self.nC,3),order='F')
        print '# Control points =',self.nC

    def computePindices(self):
        self.surf_index_P = PUBSlib.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        self.edge_index_P = PUBSlib.getedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_n)
        self.nT = 0
        self.nT += self.edge_index_P[-1,1]
        self.nT += 2*self.surf_index_P[-1,1]
        self.nP = self.nvert
        self.nP += self.edge_index_P[-1,1]
        self.nP += self.surf_index_P[-1,1]
        self.P = numpy.zeros((self.nP,3),order='F')
        print '# Points =',self.nP

    def computeDindices(self):
        self.knot_index = PUBSlib.getknotindices(self.ngroup, self.group_k, self.group_m)

    def computeKnots(self):
        self.group_d = PUBSlib.getd(self.ngroup,self.nD,self.group_k,self.group_m)

    def computeParameters(self):
        self.T = PUBSlib.initializet(self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.knot_index)
        
    def computeJacobian(self):
        self.nJ = PUBSlib.getjnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_k,self.group_n,self.edge_count)
        Ja, Ji, Jj = PUBSlib.getjacobian(self.nJ, self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.surf_index_P, self.edge_index_P, self.surf_index_C, self.edge_index_C, self.knot_index, self.edge_count, self.T)
        self.J = scipy.sparse.csc_matrix((Ja,(Ji,Jj)))
        print '# Jacobian non-zeros =',self.J.nnz

        self.nM = PUBSlib.getmnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C, self.edge_index_Q, self.vert_index_Q, self.edge_count, self.surf_c1, self.edge_c1)
        Ma, Mi, Mj = PUBSlib.getdofmapping(self.nM,self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_vert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C,self.edge_index_C,self.edge_index_Q,self.vert_index_Q, self.edge_count,self.surf_c1,self.edge_c1)
        self.M = scipy.sparse.csc_matrix((Ma,(Mi,Mj)))
        self.JM = self.J.dot(self.M)

    def computeControlPts(self):  
        AT = self.JM.transpose()
        ATA = AT.dot(self.JM)
        ATB = AT.dot(self.P)

        self.Q = numpy.zeros((self.JM.shape[1],3),order='F')
        for i in range(3):
            self.Q[:,i] = scipy.sparse.linalg.gmres(ATA,ATB[:,i])[0]
        self.C = self.M.dot(self.Q)

        #self.Q = numpy.zeros((self.JM.shape[1],3),order='F')
        #solve = scipy.sparse.linalg.dsolve.factorized(ATA)
        #for i in range(3):
        #    self.Q[:,i] = solve(ATB[:,i])
        #self.C = self.M.dot(self.Q)

    def computePoints(self):
        for i in range(3):
            self.P[:,i] = self.JM.dot(self.Q[:,i])
        self.C = self.M.dot(self.Q)

    def computeIndex(self, surf, u, v, quantity):
        ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
        vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
        if quantity==0:
            surf_index = self.surf_index_P
            edge_index = self.edge_index_P
            nvert = self.nvert
            group_m = self.group_n
        elif quantity==1:
            surf_index = self.surf_index_C
            edge_index = self.edge_index_C
            nvert = self.nvert
            group_m = self.group_m
        elif quantity==2:
            surf_index = self.surf_index_Q
            edge_index = self.edge_index_Q
            nvert = max(self.vert_index_Q)
            group_m = self.group_m            
        mu = group_m[ugroup-1]
        mv = group_m[vgroup-1]
        if u < 0:
            u += mu
        if v < 0:
            v += mv
        if (u==0 or u==mu-1) and (v==0 or v==mv-1) and quantity==2:
            vert = self.surf_vert[surf,int(u/(mu-1)),int(v/(mv-1))] - 1
            return self.vert_index_Q[vert] - 1
        else:
            return PUBSlib.computeindex(surf+1,u+1,v+1,mu,mv,self.nsurf,self.nedge,nvert,self.surf_vert,self.surf_edge,surf_index,edge_index) - 1

    def computePt(self,surf,u,v,uder=0,vder=0):
        ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
        vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
        ku = self.group_k[ugroup-1]
        kv = self.group_k[vgroup-1]   
        mu = self.group_m[ugroup-1]
        mv = self.group_m[vgroup-1]   
        nB = ku*kv
        P = PUBSlib.computept(surf+1,uder,vder,ku,kv,mu,mv,nB,self.nD,self.nC,self.nsurf,self.nedge,self.ngroup,self.nvert,u,v,self.surf_vert,self.surf_edge,self.edge_group,self.group_d,self.surf_index_C,self.edge_index_C,self.knot_index,self.C)
        return P

    def computeBases(self, s, u, v):
        nB = PUBSlib.computebnnz(s.shape[0], self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, s)
        Ba, Bi, Bj = PUBSlib.computebases(0, 0, s.shape[0], nB, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_d, self.surf_index_C, self.edge_index_C, self.knot_index, s, u, v)
        B = scipy.sparse.csc_matrix((Ba,(Bi,Bj)),shape=(s.shape[0],self.C.shape[0]))
        return B

    def computeProjection(self, P0, surfs=None, Q=None):
        if surfs==None:
            surfs = numpy.linspace(1,self.nsurf,self.nsurf)
        else:
            surfs = numpy.array(surfs) + 1
        if Q==None:
            s,u,v = PUBSlib.computeprojection(P0.shape[0],surfs.shape[0],self.nD,self.nT,self.nC,self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,surfs,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.knot_index,self.T,self.C,self.P,P0)
        else:
            s,u,v = PUBSlib.computepjtnalongq(P0.shape[0],surfs.shape[0],self.nD,self.nT,self.nC,self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,surfs,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.knot_index,self.T,self.C,self.P,P0,Q)
        P = numpy.zeros((P0.shape[0],3))
        for i in range(P0.shape[0]):
            P[i] = self.computePt(s[i]-1,u[i],v[i])
        return P,s,u,v

    def write2Tec(self,filename):
        f = open(filename+'.dat','w')
        f.write('title = "PUBSlib output"\n')
        f.write('variables = "x", "y", "z"\n')
        for surf in range(self.nsurf):
            ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
            vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
            nu = self.group_n[ugroup-1]
            nv = self.group_n[vgroup-1]      
            f.write('zone i='+str(nu)+', j='+str(nv)+', DATAPACKING=POINT\n')
            P = PUBSlib.getsurfacep(surf+1, self.nP, nu, nv, self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_index_P, self.edge_index_P, self.P)
            for v in range(nv):
                for u in range(nu):
                    f.write(str(P[u,v,0]) + ' ' + str(P[u,v,1]) + ' ' + str(P[u,v,2]) + '\n')
        f.close()

    def write2TecC(self,filename):
        f = open(filename+'.dat','w')
        f.write('title = "PUBSlib output"\n')
        f.write('variables = "x", "y", "z"\n')        
        f.write('zone i='+str(self.nC)+', DATAPACKING=POINT\n')
        for i in range(self.nC):
            f.write(str(self.C[i,0]) + ' ' + str(self.C[i,1]) + ' ' + str(self.C[i,2]) + '\n')
        f.close()

    def write2EGADS(self,filename):
        Ps = []
        for surf in range(self.nsurf):
            ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
            vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
            nu = self.group_n[ugroup-1]
            nv = self.group_n[vgroup-1]      
            P = PUBSlib.getsurfacep(surf+1, self.nP, nu, nv, self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_index_P, self.edge_index_P, self.P)
            Ps.append(P[:,:,:])
            Ps.append(copy.copy(P[::-1,:,:]))
            Ps[-1][:,:,2] *= -1
        export2EGADS.export(Ps, filename)

    def write2IGES(self, filename):
        def getProps(surf):
            ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
            vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
            ku = self.group_k[ugroup-1]
            kv = self.group_k[vgroup-1]      
            mu = self.group_m[ugroup-1]
            mv = self.group_m[vgroup-1]      
            du = self.group_d[self.knot_index[ugroup-1,0]:self.knot_index[ugroup-1,1]]
            dv = self.group_d[self.knot_index[vgroup-1,0]:self.knot_index[vgroup-1,1]]
            return ku,kv,mu,mv,du,dv            

        def write(f, val, dirID, parID, field, last=False):
            if last:
                f.write('%20.12e;' %(val.real))
            else:
                f.write('%20.12e,' %(val.real))
            field += 1
            if field==3:
                field = 0
                f.write('%9i' %(dirID))
                f.write('P')
                f.write('%7i\n' %(parID))
                parID += 1
            return parID, field

        f = open(filename+'.igs','w')
        f.write('                                                                        S      1\n')
        f.write('1H,,1H;,4HSLOT,37H$1$DUA2:[IGESLIB.BDRAFT.B2I]SLOT.IGS;,                G      1\n')
        f.write('17HBravo3 BravoDRAFT,31HBravo3->IGES V3.002 (02-Oct-87),32,38,6,38,15,  G      2\n')
        f.write('4HSLOT,1.,1,4HINCH,8,0.08,13H871006.192927,1.E-06,6.,                   G      3\n')
        f.write('31HD. A. Harrod, Tel. 313/995-6333,24HAPPLICON - Ann Arbor, MI,4,0;     G      4\n')

        dirID = 1
        parID = 1    
        for surf in range(self.nsurf):
            ku,kv,mu,mv,du,dv = getProps(surf)
            numFields = 4 + du.shape[0] + dv.shape[0] + 4*mu*mv
            numLines = 2 + numpy.ceil(numFields/3.0)
            for val in [128, parID, 0, 0, 1, 0, 0, 0]:
                f.write('%8i' %(val))
            f.write('00000001')
            f.write('D')
            f.write('%7i\n' %(dirID))
            dirID += 1
            for val in [128, 0, 2, numLines, 0]:
                f.write('%8i' %(val))
            f.write('%32i' %(0))
            f.write('D')
            f.write('%7i\n' %(dirID))
            dirID += 1
            parID += numLines
        nDir = dirID - 1

        dirID = 1    
        parID = 1
        for surf in range(self.nsurf):
            ku,kv,mu,mv,du,dv = getProps(surf)

            for val in [128, mu-1, mv-1, ku-1, kv-1]:
                f.write('%12i,' %(val))
            f.write('%7i' %(dirID))   
            f.write('P')
            f.write('%7i\n' %(parID))
            parID += 1

            for val in [0, 0, 1, 0, 0]:
                f.write('%12i,' %(val))
            f.write('%7i' %(dirID))   
            f.write('P')
            f.write('%7i\n' %(parID))
            parID += 1

            field = 0
            for i in range(du.shape[0]):
                parID,field = write(f, du[i], dirID, parID, field)
            for i in range(dv.shape[0]):
                parID,field = write(f, dv[i], dirID, parID, field)
            for i in range(mu*mv):
                parID,field = write(f, 1.0, dirID, parID, field)
            for j in range(mv):
                for i in range(mu):
                    C = self.C[self.computeIndex(surf,i,j,1)]
                    for k in range(3):
                        parID,field = write(f, C[k].real, dirID, parID, field)
            parID,field = write(f, 0, dirID, parID, field)
            parID,field = write(f, 1, dirID, parID, field)
            parID,field = write(f, 0, dirID, parID, field)
            parID,field = write(f, 1, dirID, parID, field, last=True)
            if not field==0:
                for i in range(3-field):
                    f.write('%21s' %(' '))
                f.write('%9i' %(dirID))
                f.write('P')
                f.write('%7i\n' %(parID))
                parID += 1

            dirID += 2

        nPar = parID - 1

        f.write('S%7i' %(1))   
        f.write('G%7i' %(4))   
        f.write('D%7i' %(nDir))   
        f.write('P%7i' %(nPar))   
        f.write('%40s' %(' '))   
        f.write('T')
        f.write('%7i\n' %(1))       
        f.close()

    def plotm(self,fig,mirror=True):
        mlab.figure(fig,bgcolor=(1,1,1))
        m = PUBSlib.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        n = PUBSlib.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        for i in range(self.nsurf):
            if 1:
                P = PUBSlib.getsurfacep(i+1, self.nP, n[i,0], n[i,1], self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_index_P, self.edge_index_P, self.P)
                mlab.mesh(P[:,:,0],P[:,:,1],P[:,:,2],color=(65/256,105/256,225/256))
                if mirror:
                    mlab.mesh(P[:,:,0],P[:,:,1],-P[:,:,2],color=(65/256,105/256,225/256))

    def plot(self,fig,mirror=True):
        ax = p3.Axes3D(fig)
        m = PUBSlib.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        n = PUBSlib.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        for i in range(self.nsurf):
            if 0:
                C = PUBSlib.getsurfacep(i+1, self.nC, m[i,0], m[i,1], self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_index_C, self.edge_index_C, self.C)
                for j in range(C.shape[0]):
                    ax.scatter(C[j,:,0],C[j,:,1],C[j,:,2])
                if mirror:
                    for j in range(C.shape[0]):
                        ax.scatter(C[j,:,0],-C[j,:,1],C[j,:,2])
            if 0:
                ax.scatter(self.C[:,0],self.C[:,1],self.C[:,2])
            if 1:
                P = PUBSlib.getsurfacep(i+1, self.nP, n[i,0], n[i,1], self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_index_P, self.edge_index_P, self.P)
                ax.plot_wireframe(P[:,:,0],P[:,:,1],P[:,:,2])
                if mirror:
                    ax.plot_wireframe(P[:,:,0],-P[:,:,1],P[:,:,2])
        mins = numpy.zeros(3)
        maxs = numpy.zeros(3)
        for i in range(3):
            mins[i] = min(self.P[:,i])
            maxs[i] = max(self.P[:,i])
        length = max(maxs-mins)
        limits = [-length/2.0,length/2.0]
        ax.scatter(limits,limits,limits,s=0.0)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        return ax
