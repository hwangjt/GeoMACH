from __future__ import division
import numpy
import time
import scipy.sparse
from collections import OrderedDict

from GeoMACH.PSM import PSMlib, QUADlib, CDTlib
from QUAD import QUAD


class Airframe(object):

    def __init__(self, geometry, maxL):
        self.geometry = geometry
        self.maxL = maxL
        self.quad = QUAD()
        self.members = []
        self.memberNames = []

    def addMember(self, name, member):
        self.members.append(member)
        self.memberNames.append(name)

    def addVertFlip(self, name, comp, p1, p2, w=[1,0], i=[0,1]):
        index = self.geometry.comps.keys().index
        self.addMember(name,[
                [[index(comp),i[0],0], [w[0],p1[0],p1[1]], [w[0],p2[0],p2[1]], [w[1],p1[0],p1[1]], [w[1],p2[0],p2[1]]],
                [[index(comp),i[1],0], [1-w[0],1-p1[0],p1[1]], [1-w[0],1-p2[0],p2[1]], [1-w[1],1-p1[0],p1[1]], [1-w[1],1-p2[0],p2[1]]],
                ])

    def addVert(self, name, comp, p1, p2, w=[1,0], i=[0,1]):
        index = self.geometry.comps.keys().index
        self.addMember(name,[
                [[index(comp),i[0],0], [w[0],p1[0],p1[1]], [w[0],p2[0],p2[1]], [w[1],p1[0],p1[1]], [w[1],p2[0],p2[1]]],
                [[index(comp),i[1],0], [1-w[0],p1[0],p1[1]], [1-w[0],p2[0],p2[1]], [1-w[1],p1[0],p1[1]], [1-w[1],p2[0],p2[1]]],
                ])

    def addCtrVert(self, name, comp1, comp2, p, w=[1,0]):
        index = self.geometry.comps.keys().index
        self.addMember(name,[
                [[index(comp1), 0,0], [w[0],p,0.0], [w[1],p,0.0], [0.0,p,0.0], [0.0,p,0.0]],
                [[index(comp1), 1,0], [1-w[0],1-p,0.0], [1-w[1],1-p,0.0], [0.0,1-p,0.0], [0.0,1-p,0.0]],
                [[index(comp2), 0,0], [0.0,p,1.0], [0.0,p,1.0], [w[0],p,1.0], [w[1],p,1.0]],
                [[index(comp2), 1,0], [0.0,1-p,1.0], [0.0,1-p,1.0], [1-w[0],1-p,1.0], [1-w[1],1-p,1.0]],
                ])

    def addCtr(self, name, comp1, comp2, i, p):
        index = self.geometry.comps.keys().index
        self.addMember(name,[
                [[index(comp1), i,0], [1,p[0],0.0], [1,p[1],0.0], [0,p[0],0.0], [0,p[1],0.0]],
                [[index(comp2), i,0], [0,p[0],1.0], [0,p[1],1.0], [1,p[0],1.0], [1,p[1],1.0]],
                ])

    def writeTecHeader(self, filename, title, variables):
        f = open(filename,'w')
        f.write('title = ' + title + '\n')
        f.write('variables = ')
        for i in range(len(variables)):
            f.write(variables[i] + ',')
        f.write('\n')
        return f

    def writeLine(self, f, data, label=''):
        f.write(label)
        for k in range(data.shape[0]):
            if data[k] == data[k]:
                f.write(str(data[k]) + ' ')
            else:
                f.write(str(0.0) + ' ')
        f.write('\n')

    def write2TecFEquads(self, filename, zones, variables=['x','y','z'], title="GeoMACH PSM output"):
        f = self.writeTecHeader(filename,title,variables) 
        for z in range(len(zones)):
            name, nodes, quads = zones[z]
            f.write('ZONE T=\"' + name +'\",')
            f.write('N=' + str(nodes.shape[0]) + ',')
            f.write('E=' + str(quads.shape[0]) + ',')
            f.write('DATAPACKING=POINT, ZONETYPE = FEQUADRILATERAL\n')
            for i in range(nodes.shape[0]):
                self.writeLine(f, nodes[i,:])
            for i in range(quads.shape[0]):
                self.writeLine(f, quads[i,:])
        f.close()  

    def preview(self, filename):
        self.preview = []

        self.computePreviewSurfaces()
        self.computeFaceDimensions()
        self.importMembers()
        self.computePreviewMembers()

        self.write2TecFEquads(filename,self.preview)

    def mesh(self):
        self.meshS = []
        self.meshM = []

        self.computeTopology()
        self.computeAdjoiningEdges()
        self.computeGroupIntersections()
        self.computeFaces()
        self.computeSurfaces()
        self.computeMembers()

    def computeMesh(self, filename):
        bse = self.geometry._bse

        B1, quads1, nnode1, mem1, ucoord1, vcoord1 = self.meshS
        B2, quads2, nnode2, mem2, ucoord2, vcoord2 = self.meshM
        nodes1 = B1 * bse.vec['cp_str'].array
        nodes2 = B2 * bse.vec['cp_str'].array

        new_nodes = numpy.vstack([nodes1, nodes2])
        new_ucoord = numpy.concatenate(ucoord1+ucoord2)
        new_vcoord = numpy.concatenate(vcoord1+vcoord2)

        mem1 = numpy.concatenate(mem1).astype(int)
        mem2 = numpy.concatenate(mem2).astype(int)
        mem2 += numpy.max(mem1) + 1
        new_mem = numpy.concatenate([mem1, mem2]).astype(int)

        mesh = []
        nodes = []
        quads = []
        quad_groups = []
        group_names = []
        nCount = 0
        nGroup = 0
        for i in range(len(self.surfaceNames)):
            mesh.append([self.surfaceNames[i], nodes1[nnode1[i]:nnode1[i+1]], quads1[i]])
            nodes.append(nodes1[nnode1[i]:nnode1[i+1]])
            quads.append(quads1[i]+nCount)
            nCount = nCount + nnode1[i+1] - nnode1[i]
            quad_groups.append(nGroup*numpy.ones(quads1[i].shape[0], int))
            group_names.append(self.surfaceNames[i])
            nGroup += 1
        for i in range(self.nmem):
            mesh.append([self.memberNames[i], nodes2[nnode2[i]:nnode2[i+1]], quads2[i]])
            nodes.append(nodes2[nnode2[i]:nnode2[i+1]])
            quads.append(quads2[i]+nCount)
            nCount = nCount + nnode2[i+1] - nnode2[i]
            quad_groups.append(nGroup*numpy.ones(quads2[i].shape[0], int))
            group_names.append(self.memberNames[i])
            nGroup += 1

        self.write2TecFEquads(filename,mesh)

        quad_groups = numpy.concatenate(quad_groups)

        nodes = numpy.array(numpy.concatenate(nodes,axis=0)[:,:3], order='F')
        quads = numpy.array(numpy.concatenate(quads,axis=0), order='F')

        nnode = nodes.shape[0]
        nquad = quads.shape[0]
        nid, ids = PSMlib.computeuniquenodes(nnode, nodes, 1e-7)
        nodes, quads = PSMlib.removeduplicatenodes(nnode, nid, nquad, ids, nodes, quads)

        nnode = nodes.shape[0]
        nquad = quads.shape[0]
        right = PSMlib.countrightquads(nnode, nquad, nodes, quads)
        quads = PSMlib.removerightquads(nquad, nquad-numpy.sum(right), quads, right)

        temp = numpy.zeros((nquad,4), int)
        temp[:,0] = quad_groups[:]
        temp = PSMlib.removerightquads(nquad, nquad-numpy.sum(right), temp, right)
        quad_groups = temp[:,0]

        nnode = nodes.shape[0]
        symm = PSMlib.identifysymmnodes(nnode, nodes)

        self.write2TecFEquads('test.dat',[['test',nodes,quads]])

        import BDFwriter
        BDFwriter.writeBDF(filename+'.bdf',nodes,quads,symm,quad_groups,group_names,
                           new_mem, new_nodes, new_ucoord, new_vcoord)

    def computePreviewSurfaces(self):
        bse = self.geometry._bse
        nsurf = bse._num['surf']

        quads,s,u,v = PSMlib.computepreviewsurfaces(4*nsurf,nsurf)
        bse.add_jacobian('temp', s, u, v, ndim=3)
        bse.apply_jacobian('temp', 'd(temp)/d(cp_str)', 'cp_str')
        nodes = bse.vec['temp'].array
        self.surfEdgeLengths = PSMlib.computeedgelengths(nodes.shape[1],nsurf,nodes,quads)

        self.preview.append(['surfs', nodes, quads])

    def computeSurfEdge(self):
        bse = self.geometry._bse
        nsurf = bse._num['surf']
        surf_ptrs = bse._topo['surf_ptrs']

        surf_edge = numpy.zeros((nsurf, 2, 2), int, order='F')
        surf_edge[:, 0, 0] = surf_ptrs[:, 1, 0]
        surf_edge[:, 0, 1] = surf_ptrs[:, 1, 2]
        surf_edge[:, 1, 0] = surf_ptrs[:, 0, 1]
        surf_edge[:, 1, 1] = surf_ptrs[:, 2, 1]

        return surf_edge
        
    def computeFaceDimensions(self):
        geometry = self.geometry
        bse = geometry._bse
        nsurf = bse._num['surf']
        nedge = bse._num['edge']
        ngroup = bse._num['group']
        surf_edge = self.computeSurfEdge()
        edge_group = bse._topo['edge_group']

        groupLengths = numpy.zeros(ngroup)
        groupCount = numpy.zeros(ngroup,int)
        for comp in geometry.comps.values():
            for face in comp.faces.values():
                ni, nj = face._num_surf['u'], face._num_surf['v']
                surf_indices = face._surf_indices
                groupLengths[:], groupCount[:] = PSMlib.addgrouplengths(ni, nj, nsurf, nedge, ngroup, surf_indices+1, surf_edge, edge_group, self.surfEdgeLengths, groupLengths, groupCount)

        groupLengths = groupLengths / groupCount
                 
        faceDims = OrderedDict()
        for comp in geometry.comps.values():
            faceDimsComp = OrderedDict()
            for face in comp.faces.values():
                ni, nj = face._num_surf['u'], face._num_surf['v']
                surf_indices = face._surf_indices
                idims, jdims = PSMlib.computefacedimensions(ni, nj, nsurf, nedge, ngroup, surf_indices+1, surf_edge, edge_group, groupLengths)
                faceDimsComp[face._name] = [idims,jdims]
            faceDims[comp._name] = faceDimsComp

        self.faceDims = faceDims

    def computeFaceDimensions0(self):
        geometry = self.geometry
        nsurf = bse._num['surf']

        faceDims = {}    
        for comp in geometry.comps.values():
            faceDimsComp = []
            jdim0 = numpy.zeros(comp.faces.values()[0]._num_surf['v']+1)
            for face in comp.faces.values():
                ni, nj = face._num_surf['u'], face._num_surf['v']
                surf_indices = face._surf_indices
                idims, jdims = PSMlib.computefacedimensions0(ni,nj,nsurf,surf_indices+1,self.surfEdgeLengths)
                jdim0 += jdims
                faceDimsComp[face._name] = [idims,jdims]
            jdim0 /= len(comp.faces)
            for face in comp.faces.values():
                faceDimsComp[face._name][1][:] = jdim0[:]
            faceDims[comp._name] = faceDimsComp

        self.faceDims = faceDims

    def importMembers(self):
        self.nmem = len(self.members)

        geometry = self.geometry
        nmem = self.nmem

        for imem in range(nmem):
            if len(self.members[imem]) is 2:
                self.members[imem].extend([[[-1,-1,-1] for j in range(5)] for i in range(2)])
        members = numpy.array(self.members, order='F')
        membersInt, membersFlt = PSMlib.importmembers(nmem, members)

        for comp in geometry.comps.values():
            for face in comp.faces.values():
                ni, nj = face._num_surf['u'], face._num_surf['v']
                surf_indices = face._surf_indices
                idims, jdims = self.faceDims[comp._name][face._name]
                PSMlib.computepreviewmembercoords(comp._num+1,face._num+1,ni,nj,nmem,surf_indices+1,idims,jdims,membersInt,membersFlt)

        self.membersInt = membersInt
        self.membersFlt = membersFlt

    def computePreviewMembers(self):
        nmem = self.nmem
        bse = self.geometry._bse
        nsurf = bse._num['surf']
        ncp = bse._size['cp_str']

        quads, Wa = PSMlib.computepreviewmemberweights(nmem, 4*nmem, self.membersFlt)
        linW = numpy.linspace(0,4*nmem-1,4*nmem)

        B0 = scipy.sparse.csr_matrix((4*nmem,ncp))
        for src in range(4):
            W = scipy.sparse.csr_matrix((Wa[:,src],(linW,linW)))
            for surf in range(nsurf):
                npts = PSMlib.countpreviewmembers(surf+1, src+1, nmem, self.membersFlt)
                if npts is not 0:
                    inds, P, Q = PSMlib.computepreviewmemberproj(surf+1, src+1, nmem, npts, self.membersFlt)
                    Ta = numpy.ones(npts)
                    Ti = inds - 1
                    Tj = numpy.linspace(0,npts-1,npts)
                    T = scipy.sparse.csr_matrix((Ta,(Ti,Tj)),shape=(4*nmem,npts))

                    mu = bse.get_bspline_option('num_cp', surf, 'u')
                    mv = bse.get_bspline_option('num_cp', surf, 'v')
                    nu = bse.get_bspline_option('num_pt', surf, 'u')
                    nv = bse.get_bspline_option('num_pt', surf, 'v')
                    
                    for u in range(mu):
                        for v in range(mv):
                            bse.vec['cp_str'](surf)[u, v, :] = [u/(mu-1), v/(mv-1), 0]
                    for u in range(nu):
                        for v in range(nv):
                            bse.vec['pt_str'](surf)[u, v, :] = [u/(nu-1), v/(nv-1), 0]

                    bse.compute_projection('temp', P, [surf], ndim=3)
                    B = bse.jac['d(temp)/d(cp_str)']
                    B0 = B0 + W * T * B

        bse.apply_jacobian('cp_str', 'd(cp_str)/d(cp)', 'cp')
        nodes = B0 * bse.vec['cp_str'].array
        self.memEdgeLengths = PSMlib.computeedgelengths(nodes.shape[1],nmem,nodes,quads)

        self.preview.append(['members', nodes, quads])

    def computeTopology(self):
        nmem = self.nmem
        bse = self.geometry._bse
        nsurf = bse._num['surf']

        Ps = numpy.zeros((nsurf,3,3,3),order='F')
        for s in range(nsurf):
            for i in range(3):
                for j in range(3):
                    bse.add_jacobian('temp', [s], [i/2.0], [j/2.0], ndim=3)
                    Ps[s,i,j] = bse.jac['d(temp)/d(cp_str)'] * bse.vec['cp_str'].array
        nvertS,ngroupS,surf_vert,surf_group = PSMlib.initializeconnectivities(nsurf,1e-13,1e-5,Ps)
        nvertM,ngroupM,mem_vert,mem_group = PSMlib.computemembertopology(nmem, self.membersInt, self.membersFlt)
        mem_group[:,:,:] += ngroupS

        self.surf_group = surf_group
        self.mem_group = mem_group
        self.ngroupS = ngroupS
        self.ngroupM = ngroupM

    def computeAdjoiningEdges(self):
        nmem = self.nmem

        nadj = PSMlib.countadjoiningedges(nmem, self.membersFlt)
        adjoiningInt, adjoiningFlt = PSMlib.computeadjoiningedges(nmem, nadj, self.membersInt, self.membersFlt)

        self.adjoiningInt = adjoiningInt
        self.adjoiningFlt = adjoiningFlt

    def computeGroupIntersections(self):
        geometry = self.geometry
        quad = self.quad
        nsurf = geometry._bse._num['surf']
        nmem = self.nmem
        ngroup = self.ngroupS + self.ngroupM
        nadj = self.adjoiningInt.shape[0]

        premeshFaces = []
        groupIntCount = numpy.zeros(ngroup,int)
        for comp in geometry.comps.values():
            for face in comp.faces.values():
                ni, nj = face._num_surf['u'], face._num_surf['v']
                surf_indices = face._surf_indices
                idims, jdims = self.faceDims[comp._name][face._name]
                nedge = PSMlib.countfaceedges(comp._num+1, face._num+1, ni, nj, nadj, self.adjoiningInt)
                edge_group, edgeLengths, edges = PSMlib.computefaceedges(comp._num+1, face._num+1, ni, nj, nsurf, nmem, nadj, nedge, idims, jdims, surf_indices+1, self.surf_group, self.mem_group, self.adjoiningInt, self.adjoiningFlt, self.surfEdgeLengths, self.memEdgeLengths)
                quad.importEdges(edges)
                quad.addIntersectionPts()
                quad.removeDuplicateVerts()
                verts, edges = quad.verts, quad.edges

                groupIntCount = PSMlib.countgroupintersections(verts.shape[0], edges.shape[0], ngroup, verts, edges, edge_group, groupIntCount)

                premeshFaces.append([verts,edges,edge_group,edgeLengths])

        groupIntPtr = PSMlib.computegroupintptr(ngroup, groupIntCount)
        nint = groupIntPtr[-1,-1]

        iList = 0
        groupInts = numpy.zeros(nint)
        for comp in geometry.comps.values():
            for face in comp.faces.values():
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                iList += 1
                groupInts = PSMlib.computegroupintersections(verts.shape[0], edges.shape[0], ngroup, nint, verts, edges, edge_group, groupIntPtr, groupInts)

        groupSplitCount = PSMlib.countgroupsplits(nsurf, nmem, ngroup, nint, self.maxL, self.surf_group, self.mem_group, self.surfEdgeLengths, self.memEdgeLengths, groupIntPtr, groupInts)

        groupSplitPtr = PSMlib.computegroupintptr(ngroup, groupSplitCount)
        nsplit = groupSplitPtr[-1,-1]

        groupSplits = PSMlib.computegroupsplits(nsurf, nmem, ngroup, nint, nsplit, self.maxL, self.surf_group, self.mem_group, self.surfEdgeLengths, self.memEdgeLengths, groupIntPtr, groupInts, groupSplitPtr)

        self.premeshFaces = premeshFaces
        self.groupIntPtr = groupIntPtr
        self.groupInts = groupInts
        self.groupSplitPtr = groupSplitPtr
        self.groupSplits = groupSplits

    def computeFaces(self):
        geometry = self.geometry
        quad = self.quad
        ngroup = self.ngroupS + self.ngroupM
        premeshFaces = self.premeshFaces
        groupIntPtr = self.groupIntPtr
        groupInts = self.groupInts
        groupSplitPtr = self.groupSplitPtr
        groupSplits = self.groupSplits
        nint = groupIntPtr[-1,-1]
        nsplit = groupSplitPtr[-1,-1]

        iList = 0
        for comp in geometry.comps.values():
            for face in comp.faces.values():
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                nvert = PSMlib.countintersectionverts(edges.shape[0], ngroup, edge_group, groupIntPtr, groupSplitPtr)
                verts = PSMlib.computeintersectionverts(verts.shape[0], edges.shape[0], ngroup, nint, nsplit, nvert + verts.shape[0], verts, edges, edge_group, groupIntPtr, groupInts, groupSplitPtr, groupSplits)
                quad.importVertsNEdges(verts, edges)
                quad.removeDuplicateVerts()
                quad.splitEdges()
                quad.removeDuplicateEdges()
                premeshFaces[iList][0] = quad.verts
                premeshFaces[iList][1] = quad.edges
                iList += 1

    def computeSurfaces(self):
        geometry = self.geometry
        bse = geometry._bse
        quad = self.quad
        premeshFaces = self.premeshFaces

        self.surfaceNames = []
        B0 = []
        quads0 = []
        nnode0 = [0]
        mem0 = []
        ucoord0 = []
        vcoord0 = []
        iList = 0
        imem = 0
        for comp in geometry.comps.values():
            for face in comp.faces.values():
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                iList += 1
                nvert = verts.shape[0]
                nedge = edges.shape[0]
                ni, nj = face._num_surf['u'], face._num_surf['v']
                idims, jdims = self.faceDims[comp._name][face._name]
                print 'Computing skin elements:', comp._name, face._name
                for i in range(ni):
                    for j in range(nj):
                        surf = face._surf_indices[i,j]
                        if surf >= 0:
                            nedge1 = PSMlib.countsurfaceedges(nvert, nedge, idims[i], idims[i+1], jdims[j], jdims[j+1], verts, edges)
                            edges1 = PSMlib.computesurfaceedges(nvert, nedge, nedge1, idims[i], idims[i+1], jdims[j], jdims[j+1], verts, edges)

                            print comp._name, face._name, i, j
                            quad.importEdges(edges1)

                            output = False
                            if comp._name=='lw' and face._name=='upp' and i==0 and j==2:
                                output = True
                            nodes, quads = quad.mesh(self.maxL, self.surfEdgeLengths[surf,:,:], output, output)
                            
                            mu = bse.get_bspline_option('num_cp', surf, 'u')
                            mv = bse.get_bspline_option('num_cp', surf, 'v')
                            nu = bse.get_bspline_option('num_pt', surf, 'u')
                            nv = bse.get_bspline_option('num_pt', surf, 'v')

                            for u in range(mu):
                                for v in range(mv):
                                    bse.vec['cp_str'](surf)[u, v, :] = [u/(mu-1), v/(mv-1), 0]
                            for u in range(nu):
                                for v in range(nv):
                                    bse.vec['pt_str'](surf)[u, v, :] = [u/(nu-1), v/(nv-1), 0]

                            P0, Q = PSMlib.computesurfaceprojections(nodes.shape[0], nodes)

                            bse.compute_projection('temp', P0, [surf], ndim=3)
                            B = bse.jac['d(temp)/d(cp_str)']

                            name = comp._name + ':' + str(face._name) + ':' + str(i) + ':' + str(j)

                            self.surfaceNames.append(name)
                            B0.append(B)
                            nnode0.append(nnode0[-1] + P0.shape[0])
                            quads0.append(quads)
                            mem0.append(imem*numpy.ones(P0.shape[0]))
                            imem += 1
                            ucoord0.append(P0[:,0])
                            vcoord0.append(P0[:,1])

        bse.apply_jacobian('cp_str', 'd(cp_str)/d(cp)', 'cp')

        B0 = scipy.sparse.vstack(B0)

        self.meshS = [B0, quads0, nnode0, mem0, ucoord0, vcoord0]

    def computeMembers(self):
        nmem = self.nmem
        geometry = self.geometry
        bse = geometry._bse
        groupIntPtr = self.groupIntPtr
        groupInts = self.groupInts
        groupSplitPtr = self.groupSplitPtr
        groupSplits = self.groupSplits
        quad = self.quad
        ngroup = self.ngroupS + self.ngroupM
        nint = groupIntPtr[-1,-1]
        nsplit = groupSplitPtr[-1,-1]
        nsurf = bse._num['surf']
        ncp = bse._size['cp_str']

        nodesInt0 = []
        nodesFlt0 = []
        quads0 = []
        nnode0 = [0]
        mem0 = []
        ucoord0 = []
        vcoord0 = []
        for imem in range(nmem):
            print 'Computing internal members:', self.memberNames[imem]
            edges, edge_group = PSMlib.computememberedges(imem+1, nmem, self.mem_group)
            quad.importEdges(edges)
            verts, edges = quad.verts, quad.edges
            nvert = PSMlib.countintersectionverts(edges.shape[0], ngroup, edge_group, groupIntPtr, groupSplitPtr)
            verts = PSMlib.computeintersectionverts(verts.shape[0], edges.shape[0], ngroup, nint, nsplit, nvert + verts.shape[0], verts, edges, edge_group, groupIntPtr, groupInts, groupSplitPtr, groupSplits)
            quad.importVertsNEdges(verts, edges)
            nodes, quads = quad.mesh(self.maxL, self.memEdgeLengths[imem,:,:])
            nodesInt, nodesFlt = PSMlib.computemembernodes(imem+1, nmem, nodes.shape[0], self.membersInt, self.membersFlt, nodes)
            nodesInt0.append(nodesInt)
            nodesFlt0.append(nodesFlt)
            quads0.append(quads)
            nnode0.append(nnode0[-1] + nodes.shape[0])
            P0, Q = PSMlib.computesurfaceprojections(nodes.shape[0], nodes)
            mem0.append(imem*numpy.ones(P0.shape[0]))
            ucoord0.append(P0[:,0])
            vcoord0.append(P0[:,1])

        nodesInt = numpy.array(numpy.vstack(nodesInt0),order='F')
        nodesFlt = numpy.array(numpy.vstack(nodesFlt0),order='F')
        nnode = nodesInt.shape[0]

        for comp in geometry.comps.values():
            for face in comp.faces.values():
                ni, nj = face._num_surf['u'], face._num_surf['v']
                surf_indices = face._surf_indices
                idims, jdims = self.faceDims[comp._name][face._name]
                PSMlib.computememberlocalcoords(comp._num+1, face._num+1, ni, nj, nnode, idims, jdims, surf_indices+1, nodesInt, nodesFlt)

        linW = numpy.linspace(0,nnode-1,nnode)
        B0 = scipy.sparse.csr_matrix((nnode,ncp))
        for src in range(4):
            W = scipy.sparse.csr_matrix((nodesFlt[:,src,0],(linW,linW)))
            for surf in range(nsurf):
                npts = PSMlib.countmembers(surf+1, src+1, nnode, nodesInt)
                if npts is not 0:
                    inds, P, Q = PSMlib.computememberproj(surf+1, src+1, nnode, npts, nodesInt, nodesFlt)
                    Ta = numpy.ones(npts)
                    Ti = inds - 1
                    Tj = numpy.linspace(0,npts-1,npts)
                    T = scipy.sparse.csr_matrix((Ta,(Ti,Tj)),shape=(nnode,npts))

                    mu = bse.get_bspline_option('num_cp', surf, 'u')
                    mv = bse.get_bspline_option('num_cp', surf, 'v')
                    nu = bse.get_bspline_option('num_pt', surf, 'u')
                    nv = bse.get_bspline_option('num_pt', surf, 'v')
                    
                    for u in range(mu):
                        for v in range(mv):
                            bse.vec['cp_str'](surf)[u, v, :] = [u/(mu-1), v/(mv-1), 0]
                    for u in range(nu):
                        for v in range(nv):
                            bse.vec['pt_str'](surf)[u, v, :] = [u/(nu-1), v/(nv-1), 0]

                    bse.compute_projection('temp', P, [surf], ndim=3)
                    B = bse.jac['d(temp)/d(cp_str)']
                    B0 = B0 + W * T * B

        bse.apply_jacobian('cp_str', 'd(cp_str)/d(cp)', 'cp')

        self.meshM = [B0, quads0, nnode0, mem0, ucoord0, vcoord0]
