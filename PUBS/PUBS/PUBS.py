from __future__ import division
import PUBSlib, PUBSexport
import numpy, scipy, pylab, copy, time
import math
import scipy.sparse, scipy.sparse.linalg


class PUBS(object):
    """ A surface model represented as a watertight union of B-spline surfaces. 

    group_d: double(nD)
        Global list of knots
    T: double(nT)
        Global list of parameter values

    Q: double(nQ,3)
        Vector of nQ DOFs
    C: double(nC,3)
        Vector of nC control points
    P: double(nP,3)
        Vector of nP points

    J: sparse double(nP,nC)
        Jacobian of points with respect to control points
    M: sparse double(nC,nQ)
        Jacobian of control points with respect to DOFs
    JM: sparse double(nP,nQ)
        Jacobian of points with respect to DOFs

    surf_index_Q: integer(nsurf,2)
        First-1 and last index of the interior Qs for each surface in the global list of Qs
    edge_index_Q: integer(nedge,2)
        First-1 and last index of the interior Qs for each edge in the global list of Qs
    vert_index_Q: integer(nvert)
        0 if not a DOF; otherwise, the index of the vertex Q in the global list of Qs
    surf_index_C: integer(nsurf,2)
        First-1 and last index of the interior Cs for each surface in the global list of Cs
    edge_index_C: integer(nedge,2)
        First-1 and last index of the interior Cs for each edge in the global list of Cs
    surf_index_P: integer(nsurf,2)
        First-1 and last index of the interior Ps for each surface in the global list of Ps
    edge_index_P: integer(nedge,2)
        First-1 and last index of the interior Ps for each edge in the global list of Ps
    knot_index: integer(ngroup,2)
        First-1 and last index of the knot vector in the global list of knots

    surf_vert: integer(nsurf,2,2)
        Mapping from surfaces to their corner vertices
    surf_edge: integer(nsurf,2,2)
        Mapping from surfaces to their 4 edges
    edge_group: integer(nedge)
        Mapping from edges to their associated group
    vert_count: integer(nvert)
        Number of surfaces each vertex touches
    edge_count: integer(nedge)
        Number of edges each vertex touches

    """

    def __init__(self, P_arrays, ratio=3.0, printInfo=False):
        """ Create an instance by specifying a list of surfaces

        * Input *
        P_arrays: list of ndarrays (nu,nv,3)
            Each element of the list is an nu x nv array of x-y-z coordinates
        ratio: integer
            Target ratio of points to control points for all edges  
        printInfo: boolean
            Whether to print output

        """

        self.printInfo = printInfo
        self.symmPlane = 2
        self.initializeTopology(P_arrays, ratio)
        self.computeQindices()
        self.computeCindices()
        self.computePindices()
        self.computeDindices()
        self.initializePoints(P_arrays)
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()
        self.computeControlPts()
        self.computePoints()

    def updateBsplines(self):
        """ Method to call after B-spline order, number of control points, or DOFs
            has been changed along any edge """

        self.computeQindices()
        self.computeCindices()
        self.computeDindices()
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()
        self.computeControlPts()
        self.computePoints()

    def updateEvaluation(self):
        """ Method to call after number of points has been changed along any edge """

        self.computePindices()
        self.computeParameters()
        self.computeJacobian()
        self.computePoints()

    def initializeTopology(self, P_arrays, ratio):
        """ Determine connectivities - mappings from surfaces to vertices and edges; from edges to groups
            Initialize order (k), # control points (m), and # points (n) for each edge 

        * Input *
        P_arrays: list of doubles(nu,nv,3)
            Each element of the list is an nu x nv array of x-y-z coordinates

        """

        self.nsurf = len(P_arrays)
        Ps = numpy.zeros((self.nsurf,3,3,3),order='F')
        for k in range(self.nsurf):
            n = P_arrays[k].shape[0:2]
            for i in range(2):
                for j in range(2):
                    Ps[k,i*2,j*2] = P_arrays[k][-i,-j]
            left = 0.5*(P_arrays[k][0,int(numpy.ceil((n[1]-1)/2.0))] + P_arrays[k][0,int(numpy.floor((n[1]-1)/2.0))])
            right = 0.5*(P_arrays[k][-1,int(numpy.ceil((n[1]-1)/2.0))] + P_arrays[k][-1,int(numpy.floor((n[1]-1)/2.0))])
            bottom = 0.5*(P_arrays[k][int(numpy.ceil((n[0]-1)/2.0)),0] + P_arrays[k][int(numpy.floor((n[0]-1)/2.0)),0])
            top = 0.5*(P_arrays[k][int(numpy.ceil((n[0]-1)/2.0)),-1] + P_arrays[k][int(numpy.floor((n[0]-1)/2.0)),-1])
            Ps[k,0,1] = left
            Ps[k,2,1] = right
            Ps[k,1,0] = bottom
            Ps[k,1,2] = top
        ns = numpy.zeros((self.nsurf,2),order='F')
        for k in range(self.nsurf):
            ns[k,:] = P_arrays[k].shape[0:2]

        self.nvert,self.nedge,self.surf_vert,self.surf_edge = PUBSlib.computetopology(self.nsurf,1e-13,1e-5,Ps)
        self.vert_count,self.edge_count = PUBSlib.countveptrs(self.nsurf,self.nvert,self.nedge,self.surf_vert,self.surf_edge)
        self.ngroup,self.edge_group = PUBSlib.computegroups(self.nsurf,self.nedge,self.surf_edge)
        self.surf_c1 = numpy.zeros((self.nsurf,3,3),bool,order='F')
        self.edge_c1 = numpy.zeros((self.nedge,2),bool,order='F')
        k = 4
        self.group_k, self.group_m, self.group_n = PUBSlib.getkmn(k, self.nsurf, self.nedge, self.ngroup, ratio, ns, self.surf_edge, self.edge_group)

        if self.printInfo:
            print '# Surfaces =',self.nsurf
            print '# Vertices =',self.nvert
            print '# Edges =',self.nedge
            print '# Groups =',self.ngroup

    def initializePoints(self, P_arrays):
        """ Rearrange the list of surfaces into a single vector of unique points 

        * Input *
        P_arrays: list of ndarrays (nu,nv,3)
            Each element of the list is an nu x nv array of x-y-z coordinates

        """

        self.P = numpy.zeros((self.nP,3),order='F')
        for k in range(self.nsurf):
            n1 = P_arrays[k].shape[0]
            n2 = P_arrays[k].shape[1]
            PUBSlib.populatep(self.nP, n1, n2, k+1, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, P_arrays[k], self.P)
        PUBSlib.avgboundaries(self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.edge_index_P, self.P)
        self.vert_symm, self.edge_symm = PUBSlib.determinesymm(self.symmPlane+1, 1e-10, 1e-8, self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.P)

    def computeQindices(self):
        """ Compute where each vertex, edge, and surface is located in the global list of Qs """

        self.surf_index_Q = PUBSlib.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_Q = PUBSlib.getedgeindicesq(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m, self.surf_c1)
        self.vert_index_Q = PUBSlib.getvertindicesq(self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_c1, self.edge_c1)
        self.nQ = 0
        self.nQ += max(self.vert_index_Q)
        self.nQ += max(self.edge_index_Q[:,1])
        self.nQ += self.surf_index_Q[-1,1]

        if self.printInfo:
            print '# Degrees of freedom =',self.nQ

    def computeCindices(self):
        """ Compute where each vertex, edge, and surface is located in the global list of Cs """

        self.surf_index_C = PUBSlib.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_C = PUBSlib.getedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_m)
        self.nC = self.nvert
        self.nC += self.edge_index_C[-1,1]
        self.nC += self.surf_index_C[-1,1]

        if self.printInfo:
            print '# Control points =',self.nC

    def computePindices(self):
        """ Compute where each vertex, edge, and surface is located in the global list of Ps """

        self.surf_index_P = PUBSlib.getsurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        self.edge_index_P = PUBSlib.getedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_n)
        self.nT = 0
        self.nT += self.edge_index_P[-1,1]
        self.nT += 2*self.surf_index_P[-1,1]
        self.nP = self.nvert
        self.nP += self.edge_index_P[-1,1]
        self.nP += self.surf_index_P[-1,1]

        if self.printInfo:
            print '# Points =',self.nP

    def computeDindices(self):
        """ Compute where each knot vector is located in the global list of knots """

        self.nD = 0
        for i in range(self.ngroup):
            self.nD += self.group_m[i] + self.group_k[i]
        self.knot_index = PUBSlib.getknotindices(self.ngroup, self.group_k, self.group_m)

    def computeKnots(self):
        """ Compute the global list of knots """

        self.group_d = PUBSlib.getd(self.ngroup,self.nD,self.group_k,self.group_m)

    def computeParameters(self):
        """ Compute the global list of parameter values """

        self.T = PUBSlib.initializet(self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.knot_index)
        
    def computeJacobian(self):
        """ Compute the global Jacobian """

        self.nJ = PUBSlib.getjnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_k,self.group_n,self.edge_count)
        Ja, Ji, Jj = PUBSlib.getjacobian(self.nJ, self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.surf_index_P, self.edge_index_P, self.surf_index_C, self.edge_index_C, self.knot_index, self.edge_count, self.T)
        self.J = scipy.sparse.csc_matrix((Ja,(Ji,Jj)))

        self.nM = PUBSlib.getmnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C, self.edge_index_Q, self.vert_index_Q, self.edge_count, self.surf_c1, self.edge_c1)
        Ma, Mi, Mj = PUBSlib.getdofmapping(self.nM,self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_vert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C,self.edge_index_C,self.edge_index_Q,self.vert_index_Q, self.edge_count,self.surf_c1,self.edge_c1)
        self.M = scipy.sparse.csc_matrix((Ma,(Mi,Mj)))
        self.JM = self.J.dot(self.M)

        if self.printInfo:
            print '# Jacobian non-zeros =',self.JM.nnz

    def computeControlPts(self):
        """ Perform fit to compute DOFs """
        
        AT = self.JM.transpose()
        ATA = AT.dot(self.JM)
        ATB = AT.dot(self.P)

        self.Q = numpy.zeros((self.JM.shape[1],3),order='F')
        
        solver = 1
        if solver==1:
            for i in range(3):
                self.Q[:,i] = scipy.sparse.linalg.cg(ATA,ATB[:,i])[0]
        elif solver==2:
            for i in range(3):
                self.Q[:,i] = scipy.sparse.linalg.gmres(ATA,ATB[:,i])[0]
        else:
            solve = scipy.sparse.linalg.dsolve.factorized(ATA)
            for i in range(3):
                self.Q[:,i] = solve(ATB[:,i])
        self.C = self.M.dot(self.Q)

    def computePoints(self):
        """ Compute matrix-vector product to find P from Q """

        self.P = numpy.zeros((self.nP,3),order='F')
        for i in range(3):
            self.P[:,i] = self.JM.dot(self.Q[:,i])
        self.C = self.M.dot(self.Q)

    def computePointsC(self):
        """ Compute matrix-vector product to find P from C """

        for i in range(3):
            self.P[:,i] = self.J.dot(self.C[:,i])

    def computeIndex(self, surf, u, v, quantity):
        """ Return the index of a Q, C, or P entry in the global list 

        * Input *
        surf: integer
            0-based surface index
        u: double
            Parametric coordinate [0,1]
        v: double
            Parametric coordinate [0,1]
        quantity: integer
            0 for P; 1 for C; 2 for Q

        * Output *
        return: integer
            0-based index in the Q, C, or P vector corresponding to (surf,u,v)

        """

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
        """ Return point or 1st or 2nd parametric derivative 

        * Input *
        surf: integer
            0-based surface index
        u: double
            Parametric coordinate [0,1]
        v: double
            Parametric coordinate [0,1]
        uder: integer
            Order of the desired derivative in the u direction
        vder: integer
            Order of the desired derivative in the v direction

        * Output *
        return: double(3)
            x-y-z values of the point or derivative requested

        """

        ugroup = self.edge_group[abs(self.surf_edge[surf,0,0])-1]
        vgroup = self.edge_group[abs(self.surf_edge[surf,1,0])-1]
        ku = self.group_k[ugroup-1]
        kv = self.group_k[vgroup-1]   
        mu = self.group_m[ugroup-1]
        mv = self.group_m[vgroup-1]   
        nB = ku*kv
        P = PUBSlib.computept(surf+1,uder,vder,ku,kv,mu,mv,nB,self.nD,self.nC,self.nsurf,self.nedge,self.ngroup,self.nvert,u,v,self.surf_vert,self.surf_edge,self.edge_group,self.group_d,self.surf_index_C,self.edge_index_C,self.knot_index,self.C)
        return P

    def computeBases(self, surf, u, v):
        """ Return matrix that multiples with C to give n points corresponding to (s,u,v)

        * Input *
        surf: integer(n)
            0-based surface index
        u: double(n)
            Parametric coordinate [0,1]
        v: double(n)
            Parametric coordinate [0,1]

        * Output *
        B: double(n,nC)
            Matrix whose rows correspond to the requested points

        """

        nB = PUBSlib.computebnnz(surf.shape[0], self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, surf)
        Ba, Bi, Bj = PUBSlib.computebases(0, 0, surf.shape[0], nB, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_d, self.surf_index_C, self.edge_index_C, self.knot_index, surf+1, u, v)
        B = scipy.sparse.csc_matrix((Ba,(Bi,Bj)),shape=(surf.shape[0],self.C.shape[0]))
        return B

    def computeProjection(self, P0, surfs=None, Q=None):
        """ Computes projections from P0 to the supplied list of surfaces
            and returns the parametric coordinates of the closest point

        * Input *
        P0: double(n,3)
            List of n points to project onto the B-spline surface model
        surfs: integer(n)
            List of surfaces (0-based) to check
        Q: double(n,3)
            Optional list of directions along which to compute projection

        * Output *
        surf: integer(n)
            0-based list of surfaces of the projected points
        u: double(n)
            List of parametric coordinates in u of the projected points
        v: double(n)
            List of parametric coordinates in v of the projected points

        """

        if surfs==None:
            surfs = numpy.linspace(1,self.nsurf,self.nsurf)
        else:
            surfs = numpy.array(surfs) + 1
        if Q==None:
            surf,u,v = PUBSlib.computeprojection(P0.shape[0],surfs.shape[0],self.nD,self.nT,self.nC,self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,surfs,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.knot_index,self.T,self.C,self.P,P0)
        else:
            surf,u,v = PUBSlib.computepjtnalongq(P0.shape[0],surfs.shape[0],self.nD,self.nT,self.nC,self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,surfs,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.knot_index,self.T,self.C,self.P,P0,Q)
        surf -= 1
        return surf,u,v

    def vstackSparse(self, Bs):
        """ Vertically stack a list of sparse matrices

        * Input *
        Bs: list of sparse doubles(mi,n)
            List of sparse matrices that all have n columns
        
        * Output *
        return: sparse double(m,n)
            Sparse matrix where m = sum(mi)
        
        """

        n0 = 0
        n1 = Bs[0].shape[1]
        data = []
        indices = []
        indptr = []
        offset = 0
        for k in range(len(Bs)):
            if not scipy.sparse.isspmatrix_csr(Bs[k]):
                Bs[k] = scipy.sparse.csr_matrix(Bs[k])
            n0 += Bs[k].shape[0]
            data.append(Bs[k].data)
            indices.append(Bs[k].indices)
            indptr.append(Bs[k].indptr[:-1]+offset)
            offset += Bs[k].indptr[-1]
        indptr.append(numpy.array([offset]))
        data = numpy.hstack(data)
        indices = numpy.hstack(indices)
        indptr = numpy.hstack(indptr)
        return scipy.sparse.csr_matrix((data, indices, indptr), shape=(n0, n1))

    def plotm(self,fig,mirror=True):
        from mayavi import mlab
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
        import mpl_toolkits.mplot3d.axes3d as p3
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
