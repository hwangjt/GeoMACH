from __future__ import division
import numpy, time
import scipy.sparse, scipy.sparse.linalg

from GeoMACH.PUBS import PUBSlib, PUBSexport


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

        self.export = PUBSexport.PUBSexport()
        self.printInfo = printInfo
        self.nvar = 6
        self.var = ['x','y','z','nx','ny','nz']
        self.symmPlane = 2
        self.__initializeTopology(P_arrays, ratio)
        self.update()
        self.__initializePoints(P_arrays)
        self.computeControlPts()
        self.computePoints()

    def addVars(self, var):
        self.var.extend(var)
        self.nvar += len(var)
        self.Q = numpy.hstack((self.Q, numpy.zeros((self.nQ,len(var)),order='F')))
        self.C = numpy.hstack((self.C, numpy.zeros((self.nC,len(var)),order='F')))
        self.P = numpy.hstack((self.P, numpy.zeros((self.nP,len(var)),order='F')))
        self.P0 = numpy.hstack((self.P0, numpy.zeros((self.Np0[-1],len(var)),order='F')))

    def update(self):
        self.computeQindices()
        self.computeCindices()
        self.computePindices()
        self.computeDindices()
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()

    def updateBsplines(self, refit=False):
        """ Method to call after B-spline order, number of control points, or DOFs
            has been changed along any edge """

        self.computeQindices()
        self.computeCindices()
        self.computeDindices()
        self.computeKnots()
        self.computeParameters()
        self.computeJacobian()
        if refit:
            self.computeControlPts()
            self.computePoints()

    def updateEvaluation(self):
        """ Method to call after number of points has been changed along any edge """

        self.computePindices()
        self.computeParameters()
        self.computeJacobian()
        self.computePoints()

    def __initializeTopology(self, P_arrays, ratio):
        """ Determine connectivities - mappings from surfaces to vertices and edges; from edges to groups
            Initialize order (k), # control points (m), and # points (n) for each edge 

        * Input *
        P_arrays: list of doubles(nu,nv,3)
            Each element of the list is an nu x nv array of x-y-z coordinates

        """

        self.nsurf = len(P_arrays)
        Ps = numpy.zeros((self.nsurf,3,3,3),order='F')
        for k in range(self.nsurf):
            nu = P_arrays[k].shape[0]
            nv = P_arrays[k].shape[1]
            for i in range(2):
                for j in range(2):
                    Ps[k,-i,-j] = P_arrays[k][-i,-j]
            for i in range(2):
                Ps[k,-i,1] += 0.5*P_arrays[k][-i,int(numpy.ceil((nv-1)/2.0))]
                Ps[k,-i,1] += 0.5*P_arrays[k][-i,int(numpy.floor((nv-1)/2.0))]
            for j in range(2):
                Ps[k,1,-j] += 0.5*P_arrays[k][int(numpy.ceil((nu-1)/2.0)),-j]
                Ps[k,1,-j] += 0.5*P_arrays[k][int(numpy.floor((nu-1)/2.0)),-j]

        self.nvert,self.nedge,self.surf_vert,self.surf_edge = PUBSlib.initializeconnectivities(self.nsurf,1e-13,1e-5,Ps)
        self.vert_count,self.edge_count = PUBSlib.initializevecounts(self.nsurf,self.nvert,self.nedge,self.surf_vert,self.surf_edge)
        self.ngroup,self.edge_group = PUBSlib.initializegroups(self.nsurf,self.nedge,self.surf_edge)
        self.surf_c1 = numpy.zeros((self.nsurf,3,3),bool,order='F')
        self.edge_c1 = numpy.zeros((self.nedge,2),bool,order='F')

        ns = numpy.zeros((self.nsurf,2),order='F')
        for k in range(self.nsurf):
            ns[k,:] = P_arrays[k].shape[0:2]
        k = 4
        self.group_k, self.group_m, self.group_n = PUBSlib.initializekmn(k, self.nsurf, self.nedge, self.ngroup, ratio, ns, self.surf_edge, self.edge_group)

        if self.printInfo:
            print '# Surfaces =',self.nsurf
            print '# Vertices =',self.nvert
            print '# Edges =',self.nedge
            print '# Groups =',self.ngroup

    def __initializePoints(self, P_arrays):
        """ Rearrange the list of surfaces into a single vector of unique points 

        * Input *
        P_arrays: list of ndarrays (nu,nv,3)
            Each element of the list is an nu x nv array of x-y-z coordinates

        """

        self.P = numpy.zeros((self.nP,3),order='F')
        for k in range(self.nsurf):
            nu = P_arrays[k].shape[0]
            nv = P_arrays[k].shape[1]
            PUBSlib.populatep(self.nP, nu, nv, k+1, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_n, self.surf_index_P, self.edge_index_P, P_arrays[k], self.P)
        PUBSlib.avgboundaries(self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.vert_count, self.edge_count, self.edge_index_P, self.P)
        self.vert_symm, self.edge_symm = PUBSlib.determinesymm(self.symmPlane+1, 1e-10, 1e-8, self.nP, self.nedge, self.ngroup, self.nvert, self.edge_group, self.group_n, self.P)
        self.P = numpy.hstack([self.P,numpy.zeros((self.nP,3),order='F')])

    def computeQindices(self):
        """ Compute where each vertex, edge, and surface is located in the global list of Qs """

        self.surf_index_Q = PUBSlib.computesurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_Q = PUBSlib.computeedgeindicesq(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m, self.surf_c1)
        self.vert_index_Q = PUBSlib.computevertindicesq(self.nsurf, self.nedge, self.nvert, self.surf_vert, self.surf_edge, self.surf_c1, self.edge_c1)
        self.nQ = 0
        self.nQ += max(self.vert_index_Q)
        self.nQ += max(self.edge_index_Q[:,1])
        self.nQ += self.surf_index_Q[-1,1]

        self.Q = numpy.zeros((self.nQ,self.nvar),order='F')  
        if self.printInfo:
            print '# Degrees of freedom =',self.nQ

    def computeCindices(self):
        """ Compute where each vertex, edge, and surface is located in the global list of Cs """

        self.surf_index_C = PUBSlib.computesurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_m)
        self.edge_index_C = PUBSlib.computeedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_m)
        self.nC = self.nvert
        self.nC += self.edge_index_C[-1,1]
        self.nC += self.surf_index_C[-1,1]

        if self.printInfo:
            print '# Control points =',self.nC

    def computePindices(self):
        """ Compute where each vertex, edge, and surface is located in the global list of Ps """

        self.surf_index_P = PUBSlib.computesurfindices(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        self.edge_index_P = PUBSlib.computeedgeindices(self.nedge, self.ngroup, self.edge_group, self.group_n)
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
        self.knot_index = PUBSlib.computeknotindices(self.ngroup, self.group_k, self.group_m)

    def computeKnots(self):
        """ Compute the global list of knots """

        self.group_d = PUBSlib.computeknots(self.ngroup,self.nD,self.group_k,self.group_m)

    def computeParameters(self):
        """ Compute the global list of parameter values """

        self.T = PUBSlib.computeparameters(self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.knot_index)
        
    def computeJacobian(self):
        """ Compute the global Jacobian """

        self.nJ = PUBSlib.computejnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_k,self.group_n,self.edge_count)
        Ja, Ji, Jj = PUBSlib.computejacobian(self.nJ, self.nT, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_n, self.group_d, self.surf_index_P, self.edge_index_P, self.surf_index_C, self.edge_index_C, self.knot_index, self.edge_count, self.T)
        self.J = scipy.sparse.csc_matrix((Ja,(Ji,Jj)))

        self.Nuv = PUBSlib.getsurfacesizes(self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_n)
        self.Np0 = numpy.zeros(self.nsurf+1)
        for s in range(self.nsurf+1):
            self.Np0[s] = int(sum(self.Nuv[:s,0]*self.Nuv[:s,1]))
        surf = numpy.zeros(self.Np0[-1],int)
        u = numpy.zeros(self.Np0[-1])
        v = numpy.zeros(self.Np0[-1])
        for s in range(self.nsurf):
            T = PUBSlib.getsurfacet(s+1, self.Nuv[s,0], self.Nuv[s,1], self.nT, self.nsurf, self.nedge, self.surf_edge, self.surf_index_P, self.edge_index_P, self.T)
            surf[self.Np0[s]:self.Np0[s+1]] = s
            u[self.Np0[s]:self.Np0[s+1]] = T[:,:,0].flatten(order='F')
            v[self.Np0[s]:self.Np0[s+1]] = T[:,:,1].flatten(order='F')
        self.J0 = self.evaluateBases(surf, u, v, 0, 0)
        self.Ju = self.evaluateBases(surf, u, v, 1, 0)
        self.Jv = self.evaluateBases(surf, u, v, 0, 1)

        self.nM = PUBSlib.computemnnz(self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C, self.edge_index_Q, self.vert_index_Q, self.edge_count, self.surf_c1, self.edge_c1)
        Ma, Mi, Mj = PUBSlib.computedofmapping(self.nM,self.nsurf,self.nedge,self.ngroup,self.nvert,self.surf_vert,self.surf_edge,self.edge_group,self.group_m,self.surf_index_C,self.edge_index_C,self.edge_index_Q,self.vert_index_Q, self.edge_count,self.surf_c1,self.edge_c1)
        self.M = scipy.sparse.csc_matrix((Ma,(Mi,Mj)))
        self.JM = self.J.dot(self.M)

        if self.printInfo:
            print '# Jacobian non-zeros =',self.JM.nnz

    def computeControlPts(self):
        """ Perform fit to compute DOFs """
        
        nvar = self.nvar
        AT = self.JM.transpose()
        ATA = AT.dot(self.JM)
        ATB = AT.dot(self.P)
      
        solver = 1
        Q = numpy.zeros((self.nQ,nvar),order='F')  
        if solver==1:
            for i in range(nvar):
                Q[:,i] = scipy.sparse.linalg.cg(ATA,ATB[:,i])[0]
        elif solver==2:
            for i in range(nvar):
                Q[:,i] = scipy.sparse.linalg.gmres(ATA,ATB[:,i])[0]
        else:
            solve = scipy.sparse.linalg.dsolve.factorized(ATA)
            for i in range(nvar):
                Q[:,i] = solve(ATB[:,i])
        self.Q = Q
        self.computePoints()

    def computePoints(self):
        """ Compute matrix-vector product to find P from Q """

        self.C = self.M.dot(self.Q)
        self.P = self.J.dot(self.C)
        self.P0 = self.J0.dot(self.C)
        self.computeNormals()

    def computeNormals(self):
        self.Pu = self.Ju.dot(self.C[:,:3])
        self.Pv = self.Jv.dot(self.C[:,:3])
        nor = numpy.cross(self.Pu,self.Pv)
        norms = numpy.sum(nor**2,axis=1)**0.5
        for k in range(3):
            self.P0[:,3+k] = nor[:,k]/norms

    def computePointsC(self):
        """ Compute matrix-vector product to find P from C """

        self.P = self.J.dot(self.C)

    def getIndex(self, surf, u, v, quantity):
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
            return PUBSlib.getindex(surf+1,u+1,v+1,mu,mv,self.nsurf,self.nedge,nvert,self.surf_vert,self.surf_edge,surf_index,edge_index) - 1

    def evaluatePoint(self,surf,u,v,uder=0,vder=0):
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
        P = PUBSlib.evaluatepoint(surf+1,uder,vder,ku,kv,mu,mv,self.nvar,self.nD,self.nC,self.nsurf,self.nedge,self.ngroup,self.nvert,u,v,self.surf_vert,self.surf_edge,self.edge_group,self.group_d,self.surf_index_C,self.edge_index_C,self.knot_index,self.C)
        return P

    def evaluateBases(self, surf, u, v, uder=0, vder=0):
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

        nB = PUBSlib.evaluatebnnz(surf.shape[0], self.nsurf, self.nedge, self.ngroup, self.surf_edge, self.edge_group, self.group_k, surf+1)
        Ba, Bi, Bj = PUBSlib.evaluatebases(uder, vder, surf.shape[0], nB, self.nD, self.nsurf, self.nedge, self.ngroup, self.nvert, self.surf_vert, self.surf_edge, self.edge_group, self.group_k, self.group_m, self.group_d, self.surf_index_C, self.edge_index_C, self.knot_index, surf+1, u, v)
        B = scipy.sparse.csc_matrix((Ba,(Bi,Bj)),shape=(surf.shape[0],self.nC))
        return B

    def evaluateProjection(self, P0, surfs=None, Q=None):
        """ Computes projections from P0 to the supplied list of surfaces
            and returns the parametric coordinates of the closest point

        * Input *
        P0: double(n,3)
            List of n points to project onto the B-spline surface model
        surfs: integer(n)
            List of surfaces (0-based) to check
        Q: double(n,3)
            Optional list of directions along which to evaluate projection

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
            surf,u,v = PUBSlib.evaluateprojection(P0.shape[0],surfs.shape[0],self.nD,self.nT,self.nC,self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,surfs,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.knot_index,self.T,self.C[:,:3],self.P[:,:3],P0)
        else:
            surf,u,v = PUBSlib.evaluatepjtnalongq(P0.shape[0],surfs.shape[0],self.nD,self.nT,self.nC,self.nP,self.nsurf,self.nedge,self.ngroup,self.nvert,surfs,self.surf_vert,self.surf_edge,self.edge_group,self.group_k,self.group_m,self.group_n,self.group_d,self.surf_index_P,self.edge_index_P,self.surf_index_C,self.edge_index_C,self.knot_index,self.T,self.C[:,:3],self.P[:,:3],P0,Q)
        surf -= 1
        return surf,u,v
    
    def edgeProperty(self, surf, p, d=None, val=None):
        """ Get/set the edge property for the u and v edges
        p: (integer)
          0: k
          1: m
          2: n
        """
        group = [self.edge_group[abs(self.surf_edge[surf,i,0])-1] - 1 for i in range(2)]

        if p==0:
            prop = self.group_k
        elif p==1:
            prop = self.group_m
        elif p==2:
            prop = self.group_n

        if not d==None:
            prop[group[d]] = val

        return [prop[group[i]] for i in range(2)]

    def exportPjtn(self, Q):
        return numpy.sum(self.P0[:,3:6]*self.J0.dot(self.M.dot(Q)),1)

    def exportPstr(self, surfs=None):
        if surfs==None:
            surfs = range(self.nsurf)

        Ps = []
        for s in surfs:
            nu, nv = self.Nuv[s,:]
            P = PUBSlib.inflatevector(nu, nv, self.nvar, nu*nv, self.P0[self.Np0[s]:self.Np0[s+1],:])
            Ps.append(P)

        return Ps

    def exportPtri(self, surfs=None):
        if surfs==None:
            surfs = range(self.nsurf)

        index = lambda i,j,ni: ni*j + i

        Tris = []
        for s in surfs:
            nu, nv = self.Nuv[s,:]
            ntri = 2*(nu-1)*(nv-1)
            Tri = self.Np0[s]*numpy.ones((ntri,3),int,order='F')
            iT = 0
            for j in range(nv-1):
                for i in range(nu-1):
                    Tri[iT,0] += index(i,j+1,nu)
                    Tri[iT,1] += index(i,j,nu)
                    Tri[iT,2] += index(i+1,j,nu)
                    iT += 1
                    Tri[iT,0] += index(i+1,j,nu)
                    Tri[iT,1] += index(i+1,j+1,nu)
                    Tri[iT,2] += index(i,j+1,nu)
                    iT += 1
            Tris.append(Tri)

        return Tris

    def exportSurfs(self):
        nsurf = self.nsurf
        ks = numpy.zeros((nsurf,2),int,order='F')
        ms = numpy.zeros((nsurf,2),int,order='F')
        ds = [[],[]]
        Cs = []
        for s in range(nsurf):
            for d in range(2):
                group = self.edge_group[abs(self.surf_edge[s,d,0])-1]
                ks[s,d] = self.group_k[group-1]
                ms[s,d] = self.group_m[group-1]
                ds[d].append(self.group_d[self.knot_index[group-1,0]:self.knot_index[group-1,1]])
            C = numpy.zeros((ms[s,0],ms[s,1],3),order='F')
            for j in range(ms[s,1]):
                for i in range(ms[s,0]):
                    C[i,j,:] = self.C[self.getIndex(s,i,j,1),:3]
            Cs.append(C)
        return ks, ms, ds, Cs

    def plot(self):
        Ps = self.exportPstr()
        self.export.plot(Ps)

    def write2Tec(self, filename):
        if not filename[-4:]=='.dat':
            filename = filename + '.dat'
        Ps = self.exportPstr()
        self.export.write2TecStruct(filename, Ps, self.var)

    def write2TecP(self, filename):
        if not filename[-4:]=='_P.dat':
            filename = filename + '_P.dat'
        self.export.write2TecScatter(filename, self.P, self.var)

    def write2TecC(self, filename):
        if not filename[-4:]=='_C.dat':
            filename = filename + '_C.dat'
        self.export.write2TecScatter(filename, self.C, self.var)

    def write2STL(self, filename):
        if not filename[-4:]=='.stl':
            filename = filename + '.stl'
        self.export.write2STL(filename, self.P0, self.exportPtri())

    def write2IGES(self, filename):
        if not filename[-4:]=='.igs':
            filename = filename + '.igs'
        ks, ms, ds, Cs = self.exportSurfs()
        self.export.write2IGES(filename, ks, ms, ds, Cs)
