from __future__ import division
import GeoMACH
import export2EGADS
import numpy, scipy, pylab, copy, time
import math
import scipy.sparse, scipy.sparse.linalg
import mpl_toolkits.mplot3d.axes3d as p3
from mayavi import mlab


class oml:

    def importSurfaces(self, P, ratio=3.0):
        self.symmPlane = 2
        self.initializeTopology(P)
        self.initializeBsplines(ratio)
        self.computeQindices()
        self.computeCindices()
        self.computePindices()
        self.initializePoints(P)
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()
        self.computeControlPts()
        self.computePoints()

    def importCGNSsurf(self, filename):
        n = GeoMACH.nsurfaces2(filename)  
        z,sizes = GeoMACH.surfacesizes2(filename, n) 
        P0 = []
        for i in range(n):
            P0.append(GeoMACH.getsurface2(filename,z[i],sizes[i,0],sizes[i,1]))
        self.importSurfaces(P0)

    def updateBsplines(self):
        self.computeQindices()
        self.computeCindices()
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
        self.nvert,self.nedge,self.surf_vert,self.surf_edge = GeoMACH.computetopology(self.nsurf,1e-15,1e-5,Ps)
        print '# Vertices =',self.nvert
        print '# Edges =',self.nedge
        self.vert_count,self.edge_count = GeoMACH.countveptrs(self.nsurf,self.nvert,self.nedge,self.surf_vert,self.surf_edge)
        self.ngroup,self.edge_group = GeoMACH.computegroups(self.nsurf,self.nedge,self.surf_edge)
        print '# Groups =',self.ngroup
        self.surf_c1 = numpy.zeros((self.nsurf,3,3),bool,order='F')
        self.edge_c1 = numpy.zeros((self.nedge,2),bool,order='F')

    def initializeBsplines(self, ratio):
        k = 4
        self.group_k, self.group_m, self.group_n = GeoMACH.getkmn(k, self.nsurf, self.nedge, self.ngroup, ratio, self.ns, self.surf_edge, self.edge_group)

    def initializePoints(self, P):
        for k in range(self.nsurf):
            n1 = P[k].shape[0]
            n2 = P[k].shape[1]
            GeoMACH.populatep(self.nP, n1, n2, k+1, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.surf_index_P, self.edge_index_P, P[k], self.P)
        GeoMACH.avgboundaries(self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.edge_index_P, self.P)
        self.vert_symm, self.edge_symm = GeoMACH.determinesymm(self.symmPlane+1, 1e-10, 1e-8, self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.P)

    def computeQindices(self):
        self.surf_index_Q = GeoMACH.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_Q = GeoMACH.getedgeindicesq(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m, self.surf_c1)
        self.vert_index_Q = GeoMACH.getvertindicesq(self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_c1, self.edge_c1)
        self.nQ = 0
        self.nQ += max(self.vert_index_Q)
        self.nQ += max(self.edge_index_Q[:,1])
        self.nQ += self.surf_index_Q[-1,1]
        print '# Degrees of freedom =',self.nQ

    def computeCindices(self):
        self.surf_index_C = GeoMACH.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_C = GeoMACH.getedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_m)
        self.nD = 0
        for i in range(self.ngroup):
            self.nD += self.group_m[i] + self.group_k[i]
        self.nC = self.nvert
        self.nC += self.edge_index_C[-1,1]
        self.nC += self.surf_index_C[-1,1]
        self.C = numpy.zeros((self.nC,3),order='F')
        print '# Control points =',self.nC

    def computePindices(self):
        self.surf_index_P = GeoMACH.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        self.edge_index_P = GeoMACH.getedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_n)
        self.nT = 0
        self.nT += self.edge_index_P[-1,1]
        self.nT += 2*self.surf_index_P[-1,1]
        self.nP = self.nvert
        self.nP += self.edge_index_P[-1,1]
        self.nP += self.surf_index_P[-1,1]
        self.P = numpy.zeros((self.nP,3),order='F')
        print '# Points =',self.nP

    def computeKnots(self):
        self.group_d = GeoMACH.getd(self.ngroup,self.nD,self.group_k,self.group_m)

    def computeParameters(self):
        self.T = GeoMACH.initializet(self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d)
        
    def computeJacobian(self):
        self.nJ = GeoMACH.getjnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_k,self.group_n,self.vert_count,self.edge_count)
        Ja, Ji, Jj = GeoMACH.getjacobian(self.nP, self.nJ, self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.surf_index_P, self.edge_index_P, self.surf_index_C, self.edge_index_C, self.edge_count, self.T)
        self.J = scipy.sparse.csr_matrix((Ja,(Ji,Jj)))
        print '# Jacobian non-zeros =',self.J.nnz

        self.nM = GeoMACH.getmnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C, self.edge_index_Q, self.vert_index_Q, self.edge_count, self.surf_c1, self.edge_c1)
        Ma, Mi, Mj = GeoMACH.getdofmapping(self.nM,self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_vert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C,self.edge_index_C,self.edge_index_Q,self.vert_index_Q, self.edge_count,self.surf_c1,self.edge_c1)
        self.M = scipy.sparse.csr_matrix((Ma,(Mi,Mj)))
        self.JM = self.J.dot(self.M)

    def computeControlPts(self):  
        AT = self.JM.transpose()
        ATA = AT.dot(self.JM)
        ATB = AT.dot(self.P)

        self.Q = numpy.zeros((self.JM.shape[1],3),order='F')
        for i in range(3):
            self.Q[:,i] = scipy.sparse.linalg.gmres(ATA,ATB[:,i])[0]
        self.C = self.M.dot(self.Q)

    def computePoints(self):
        for i in range(3):
            self.P[:,i] = self.JM.dot(self.Q[:,i])
        self.C = self.M.dot(self.Q)

    def computeIndex(self,surf,u,v,quantity):
        ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
        vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
        if quantity==0:
            surf_index = self.surf_index_P
            edge_index = self.edge_index_P
            nvert = self.nvert
            mu = self.group_n[ugroup-1]
            mv = self.group_n[vgroup-1]
        elif quantity==1:
            surf_index = self.surf_index_C
            edge_index = self.edge_index_C
            nvert = self.nvert
            mu = self.group_m[ugroup-1]
            mv = self.group_m[vgroup-1]   
        elif quantity==2:
            surf_index = self.surf_index_Q
            edge_index = self.edge_index_Q
            nvert = max(self.vert_index_Q)
            mu = self.group_m[ugroup-1]
            mv = self.group_m[vgroup-1]
        if u < 0:
            u += mu
        if v < 0:
            v += mv
        if (u==0 or u==mu-1) and (v==0 or v==mv-1):
            vert = self.surf_vert[surf,int(u/(mu-1)),int(v/(mv-1))] - 1
            if quantity==2:
                vert = self.vert_index_Q[vert] - 1
            return vert
        elif (v==0 or v==mv-1):
            edge = self.surf_edge[surf,0,int(v/(mv-1))]
            if edge < 0:
                u = mu-1 - u
            edge = abs(edge) - 1
            return nvert + edge_index[edge,0] + u - 1
        elif (u==0 or u==mu-1):
            edge = self.surf_edge[surf,1,int(u/(mu-1))]
            if edge < 0:
                v = mv-1 - v
            edge = abs(edge) - 1
            return nvert + edge_index[edge,0] + v - 1
        else:
            return nvert + max(edge_index[:,1]) + surf_index[surf,0] + (v-1)*(mu-2) + (u-1)

    def computePt(self,surf,u,v,uder=0,vder=0):
        P = GeoMACH.computept(surf+1,uder,vder,self.nD,self.nC,self.nsurf,self.nedge,self.ngroup,self.nvert,u,v,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.C)
        return P

    def computeProjection(self,P0,surf=None):
        if len(P0.shape)==1:
            P0 = numpy.array([[P0]],order='F')   
        if surf==None:
            projections = []
            for surf in range(self.nsurf):
                ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
                vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
                nu = self.group_n[ugroup-1]
                nv = self.group_n[vgroup-1]   
                minP,minu,minv = GeoMACH.computeprojection(surf+1,nu,nv,self.nD,self.nT,self.nC,P0.shape[0],P0.shape[1],self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.T,self.C,P0,self.P)
                projections.append([minP,minu,minv])
            minP = numpy.zeros((P0.shape[0],P0.shape[1],3))
            mins = numpy.zeros((P0.shape[0],P0.shape[1]),int)
            minu = numpy.zeros((P0.shape[0],P0.shape[1]))
            minv = numpy.zeros((P0.shape[0],P0.shape[1]))
            for u in range(P0.shape[0]):
                for v in range(P0.shape[1]):
                    mind = 1e10
                    for i in range(self.nsurf):
                        d = numpy.linalg.norm(projections[i][0][u,v]-P0)
                        if d < mind:
                            mind = d
                            mins[u,v] = i
                            minu[u,v] = projections[i][1][u,v]
                            minv[u,v] = projections[i][2][u,v]
                            minP[u,v] = projections[i][0][u,v]
            return minP,mins,minu,minv
        else:
            ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
            vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
            nu = self.group_n[ugroup-1]
            nv = self.group_n[vgroup-1]      
            minP,minu,minv = GeoMACH.computeprojection(surf+1,nu,nv,self.nD,self.nT,self.nC,P0.shape[0],P0.shape[1],self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.T,self.C,P0,self.P)
            return minP,surf,minu,minv

    def write2Tec(self,filename):
        f = open(filename+'.dat','w')
        f.write('title = "GeoMACH output"\n')
        f.write('variables = "x", "y", "z"\n')
        for surf in range(self.nsurf):
            ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
            vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
            nu = self.group_n[ugroup-1]
            nv = self.group_n[vgroup-1]      
            f.write('zone i='+str(nu)+', j='+str(nv)+', DATAPACKING=POINT\n')
            P = GeoMACH.getsurfacep(surf+1, self.nP, nu, nv, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, self.P)
            for v in range(nv):
                for u in range(nu):
                    f.write(str(P[u,v,0]) + ' ' + str(P[u,v,1]) + ' ' + str(P[u,v,2]) + '\n')
        f.close()

    def write2EGADS(self,filename):
        Ps = []
        for surf in range(self.nsurf):
            ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
            vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
            nu = self.group_n[ugroup-1]
            nv = self.group_n[vgroup-1]      
            P = GeoMACH.getsurfacep(surf+1, self.nP, nu, nv, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, self.P)
            Ps.append(P[:,:,:])
            Ps.append(copy.copy(P[::-1,:,:]))
            Ps[-1][:,:,2] *= -1
        export2EGADS.export(Ps, filename)

    def plotm(self,fig,mirror=True):
        mlab.figure(fig)
        m = GeoMACH.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        n = GeoMACH.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        for i in range(self.nsurf):
            if 1:
                P = GeoMACH.getsurfacep(i+1, self.nP, n[i,0], n[i,1], self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, self.P)
                mlab.mesh(P[:,:,0],P[:,:,1],P[:,:,2],color=(65/256,105/256,225/256))
                if mirror:
                    mlab.mesh(P[:,:,0],P[:,:,1],-P[:,:,2],color=(65/256,105/256,225/256))

    def plot(self,fig,mirror=True):
        ax = p3.Axes3D(fig)
        m = GeoMACH.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        n = GeoMACH.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        for i in range(self.nsurf):
            if 0:
                C = GeoMACH.getsurfacep(i+1, self.nC, m[i,0], m[i,1], self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_m, self.surf_index_C, self.edge_index_C, self.C)
                for j in range(C.shape[0]):
                    ax.scatter(C[j,:,0],C[j,:,1],C[j,:,2])
                if mirror:
                    for j in range(C.shape[0]):
                        ax.scatter(C[j,:,0],-C[j,:,1],C[j,:,2])
            if 1:
                ax.scatter(self.C[:,0],self.C[:,1],self.C[:,2])
            if 1:
                P = GeoMACH.getsurfacep(i+1, self.nP, n[i,0], n[i,1], self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, self.P)
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
