from __future__ import division
import GeoMACH
import numpy, scipy, pylab, copy, time
import math
import scipy.sparse, scipy.sparse.linalg
import mpl_toolkits.mplot3d.axes3d as p3
from collections import deque


class oml:

    def importSurfaces(self, P, ratio=4.0):
        self.initializeTopology(P)
        self.initializeBsplines(ratio)
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
        self.vert_fillet = numpy.ones(self.nvert,bool)
        self.edge_fillet = numpy.ones(self.nedge,bool)
        self.ngroup,self.edge_group = GeoMACH.computegroups(self.nsurf,self.nedge,self.surf_edge)
        print '# Groups =',self.ngroup

    def initializeBsplines(self, ratio):
        k = 4
        self.group_k, self.group_m, self.group_n = GeoMACH.getkmn(k, self.nsurf, self.nedge, self.ngroup, ratio, self.ns, self.surf_edge, self.edge_group)

    def initializePoints(self, P):
        for k in range(self.nsurf):
            n1 = P[k].shape[0]
            n2 = P[k].shape[1]
            GeoMACH.populatep(self.nP, n1, n2, k+1, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.surf_index_P, self.edge_index_P, P[k], self.P)
        GeoMACH.avgboundaries(self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.edge_index_P, self.P)
        self.vert_symm, self.edge_symm = GeoMACH.determinesymm(2, 1e-10, 1e-8, self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.P)

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

    def computeJacobianF(self, symm=False):
        Ja, Ji, Jj = GeoMACH.getjacobian(self.nP, self.nJ, self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.surf_index_P, self.edge_index_P, self.surf_index_C, self.edge_index_C, self.edge_count, self.T)
        self.nJf = GeoMACH.getjfnnz(self.nJ,self.nC,self.nedge,self.ngroup,self.nvert,self.edge_group,self.group_m,self.vert_count,self.edge_count,self.edge_index_C,Jj)
        if symm:
            vert_symm = self.vert_symm
            edge_symm = self.edge_symm
        else:
            vert_symm = numpy.zeros(self.nvert,bool)
            edge_symm = numpy.zeros(self.nedge,bool)
        Ja, Ji, Jj = GeoMACH.getjacobianf(self.nJf, self.nJ, self.nC, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_m, self.vert_count, self.edge_count, self.surf_index_C, self.edge_index_C, vert_symm, edge_symm, Ja, Ji+1, Jj+1)
        if symm:
            self.Jfs = scipy.sparse.csr_matrix((Ja,(Ji,Jj)))
        else:
            self.Jf = scipy.sparse.csr_matrix((Ja,(Ji,Jj)))
        self.nCf = self.Jf.shape[1]
        print '# Filleted control points =',self.nCf
        print '# Filleted jacobian non-zeros =',self.Jf.nnz

    def computeControlPts(self):  
        JT = self.J.transpose()
        JTJ = JT.dot(self.J)
        JTB = JT.dot(self.P)
        for i in range(3):
            self.C[:,i] = scipy.sparse.linalg.gmres(JTJ,JTB[:,i])[0]

    def computePoints(self):
        for i in range(3):
            self.P[:,i] = self.J.dot(self.C[:,i])

    def computeFilletedC(self):
        GeoMACH.computefilletedc(2, self.nC, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_m, self.vert_count, self.edge_count, self.vert_symm, self.edge_symm, self.surf_index_C, self.edge_index_C, self.C)

    def computePt(self,i,u,v,uder=0,vder=0):
        P = GeoMACH.computept(i+1,uder,vder,self.nD,self.nC,self.nsurf,self.nedge,self.ngroup,self.nvert,u,v,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.C)
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

    def plot(self,fig,mirror=True):
        ax = p3.Axes3D(fig)
        m = GeoMACH.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        n = GeoMACH.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        for i in range(self.nsurf):
            if 1:
                C = GeoMACH.getsurfacep(i+1, self.nC, m[i,0], m[i,1], self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_m, self.surf_index_C, self.edge_index_C, self.C)
                for j in range(C.shape[0]):
                    ax.scatter(C[j,:,0],C[j,:,1],C[j,:,2])
                if mirror:
                    for j in range(C.shape[0]):
                        ax.scatter(C[j,:,0],-C[j,:,1],C[j,:,2])
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
        ax.scatter(limits+length/2.0,limits,limits)
