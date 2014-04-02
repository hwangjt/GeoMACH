from __future__ import division
import numpy
import time
import scipy.sparse
import pylab

from GeoMACH.PUBS import PUBSlib
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

    def preview(self, filename):
        self.preview = []

        self.computePreviewSurfaces()
        self.computeFaceDimensions()
        self.importMembers()
        self.computePreviewMembers()

        self.geometry.oml0.export.write2TecFEquads(filename,self.preview,self.geometry.oml0.var)

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
        oml0 = self.geometry.oml0

        oml0.computePoints()

        B1, quads1, nnode1 = self.meshS
        B2, quads2, nnode2 = self.meshM
        nodes1 = B1.dot(oml0.C)
        nodes2 = B2.dot(oml0.C)

        mesh = []
        nodes = []
        quads = []
        nCount = 0
        for i in range(len(self.surfaceNames)):
            mesh.append([self.surfaceNames[i], nodes1[nnode1[i]:nnode1[i+1]], quads1[i]])
            nodes.append(nodes1[nnode1[i]:nnode1[i+1]])
            quads.append(quads1[i]+nCount)
            nCount = nCount + nnode1[i+1] - nnode1[i]
        for i in range(self.nmem):
            mesh.append([self.memberNames[i], nodes2[nnode2[i]:nnode2[i+1]], quads2[i]])
            nodes.append(nodes2[nnode2[i]:nnode2[i+1]])
            quads.append(quads2[i]+nCount)
            nCount = nCount + nnode2[i+1] - nnode2[i]

        self.geometry.oml0.export.write2TecFEquads(filename,mesh,self.geometry.oml0.var)

        nodes = numpy.array(numpy.concatenate(nodes,axis=0)[:,:3], order='F')
        quads = numpy.array(numpy.concatenate(quads,axis=0), order='F')

        nnode = nodes.shape[0]
        nquad = quads.shape[0]
        nid, ids = PSMlib.computeuniquenodes(nnode, nodes)
        nodes, quads = PSMlib.removeduplicatenodes(nnode, nid, nquad, ids, nodes, quads)

        #nnode = nodes.shape[0]
        #nquad = quads.shape[0]
        #right = PSMlib.countrightquads(nnode, nquad, nodes, quads)
        #quads = PSMlib.removerightquads(nquad, nquad-numpy.sum(right), quads, right)

        nnode = nodes.shape[0]
        symm = PSMlib.identifysymmnodes(nnode, nodes)

        self.geometry.oml0.export.write2TecFEquads('test.dat',[['test',nodes,quads]],self.geometry.oml0.var[:3])
        import BDFwriter
        BDFwriter.writeBDF('test.bdf',nodes,quads,symm)

        for i in range(3):
            nodes[:,i] *= symm
        self.geometry.oml0.export.write2TecScatter('symm.dat', nodes, ['x','y','z'])

    def computePreviewSurfaces(self):
        oml0 = self.geometry.oml0
        nsurf = oml0.nsurf

        quads,s,u,v = PSMlib.computepreviewsurfaces(4*nsurf,nsurf)
        B = oml0.evaluateBases(s,u,v)
        nodes = B.dot(oml0.C)
        self.surfEdgeLengths = PSMlib.computeedgelengths(nodes.shape[1],nsurf,nodes,quads)

        self.preview.append(['surfs', nodes, quads])

    def computeFaceDimensions(self):
        geometry = self.geometry
        oml0 = self.geometry.oml0

        groupLengths = numpy.zeros(oml0.ngroup)
        groupCount = numpy.zeros(oml0.ngroup,int)
        for comp in geometry.comps.values():
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                groupLengths, groupCount = PSMlib.addgrouplengths(ni, nj, oml0.nsurf, oml0.nedge, oml0.ngroup, comp.Ks[f]+1, oml0.surf_edge, oml0.edge_group, self.surfEdgeLengths, groupLengths, groupCount)

        groupLengths = groupLengths / groupCount
                 
        faceDims = {}     
        for comp in geometry.comps.values():
            faceDimsComp = []
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = PSMlib.computefacedimensions(ni, nj, oml0.nsurf, oml0.nedge, oml0.ngroup, comp.Ks[f]+1, oml0.surf_edge, oml0.edge_group, groupLengths)
                faceDimsComp.append([idims,jdims])
            faceDims[comp.name] = faceDimsComp

        self.faceDims = faceDims

    def computeFaceDimensions0(self):
        geometry = self.geometry
        nsurf = self.geometry.oml0.nsurf

        faceDims = {}    
        for comp in geometry.comps.values():
            faceDimsComp = []
            jdim0 = numpy.zeros(comp.Ks[0].shape[1]+1)
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = PSMlib.computefacedimensions(ni,nj,nsurf,comp.Ks[f]+1,self.surfEdgeLengths)
                jdim0 += jdims
                faceDimsComp.append([idims,jdims])
            jdim0 /= len(comp.Ks)
            for f in range(len(comp.Ks)):
                faceDimsComp[f][1][:] = jdim0[:]
            faceDims[comp.name] = faceDimsComp

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

        for k in range(len(geometry.comps)):
            comp = geometry.comps.values()[k]            
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[comp.name][f]
                PSMlib.computepreviewmembercoords(k+1,f+1,ni,nj,nmem,comp.Ks[f]+1,idims,jdims,membersInt,membersFlt)

        self.membersInt = membersInt
        self.membersFlt = membersFlt

    def computePreviewMembers(self):
        nmem = self.nmem
        oml0 = self.geometry.oml0

        quads, Wa = PSMlib.computepreviewmemberweights(nmem, 4*nmem, self.membersFlt)
        linW = numpy.linspace(0,4*nmem-1,4*nmem)

        B0 = scipy.sparse.csr_matrix((4*nmem,oml0.C.shape[0]))
        for src in range(4):
            W = scipy.sparse.csr_matrix((Wa[:,src],(linW,linW)))
            for surf in range(oml0.nsurf):
                npts = PSMlib.countpreviewmembers(surf+1, src+1, nmem, self.membersFlt)
                if npts is not 0:
                    inds, P, Q = PSMlib.computepreviewmemberproj(surf+1, src+1, nmem, npts, self.membersFlt)
                    Ta = numpy.ones(npts)
                    Ti = inds - 1
                    Tj = numpy.linspace(0,npts-1,npts)
                    T = scipy.sparse.csr_matrix((Ta,(Ti,Tj)),shape=(4*nmem,npts))

                    mu, mv = oml0.edgeProperty(surf,1)
                    for u in range(mu):
                        for v in range(mv):
                            oml0.C[oml0.getIndex(surf,u,v,1),:3] = [u/(mu-1), v/(mv-1), 0]
                    oml0.computePointsC()

                    s,u,v = oml0.evaluateProjection(P, [surf], Q)
                    B = oml0.evaluateBases(s, u, v)
                    B0 = B0 + W.dot(T.dot(B))

        oml0.computePoints()
        nodes = B0.dot(oml0.C)
        self.memEdgeLengths = PSMlib.computeedgelengths(nodes.shape[1],nmem,nodes,quads)

        self.preview.append(['members', nodes, quads])

    def computeTopology(self):
        nmem = self.nmem
        oml0 = self.geometry.oml0

        Ps = numpy.zeros((oml0.nsurf,3,3,3),order='F')
        for s in range(oml0.nsurf):
            for i in range(3):
                for j in range(3):
                    Ps[s,i,j] = oml0.evaluatePoint(s,i/2.0,j/2.0)[:3]
        nvertS,ngroupS,surf_vert,surf_group = PUBSlib.initializeconnectivities(oml0.nsurf,1e-13,1e-5,Ps)
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
        nsurf = geometry.oml0.nsurf
        nmem = self.nmem
        ngroup = self.ngroupS + self.ngroupM
        nadj = self.adjoiningInt.shape[0]

        premeshFaces = []
        groupIntCount = numpy.zeros(ngroup,int)
        for k in range(len(geometry.comps)):
            comp = geometry.comps.values()[k]     
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[comp.name][f]
                nedge = PSMlib.countfaceedges(k+1, f+1, ni, nj, nadj, self.adjoiningInt)
                edge_group, edgeLengths, edges = PSMlib.computefaceedges(k+1, f+1, ni, nj, nsurf, nmem, nadj, nedge, idims, jdims, comp.Ks[f]+1, self.surf_group, self.mem_group, self.adjoiningInt, self.adjoiningFlt, self.surfEdgeLengths, self.memEdgeLengths)
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
            for f in range(len(comp.Ks)):
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
            for f in range(len(comp.Ks)):
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
        oml0 = geometry.oml0
        quad = self.quad
        premeshFaces = self.premeshFaces

        self.surfaceNames = []
        B0 = []
        quads0 = []
        nnode0 = [0]
        iList = 0
        for comp in geometry.comps.values():
            for f in range(len(comp.Ks)):
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                iList += 1
                nvert = verts.shape[0]
                nedge = edges.shape[0]
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[comp.name][f]
                print 'Computing skin elements:', comp.name, f
                for i in range(ni):
                    for j in range(nj):
                        surf = comp.Ks[f][i,j]
                        if oml0.visible[surf]:
                            nedge1 = PSMlib.countsurfaceedges(nvert, nedge, idims[i], idims[i+1], jdims[j], jdims[j+1], verts, edges)
                            edges1 = PSMlib.computesurfaceedges(nvert, nedge, nedge1, idims[i], idims[i+1], jdims[j], jdims[j+1], verts, edges)

                            print comp.name, f, i, j
                            quad.importEdges(edges1)
                            if comp.name=='fu' and f==1 and i==0 and j==0:
                                nodes, quads = quad.mesh(self.maxL, self.surfEdgeLengths[surf,:,:])#,True,True)
                            else:
                                nodes, quads = quad.mesh(self.maxL, self.surfEdgeLengths[surf,:,:])
                            
                            mu, mv = oml0.edgeProperty(surf,1)
                            for u in range(mu):
                                for v in range(mv):
                                    oml0.C[oml0.getIndex(surf,u,v,1),:3] = [u/(mu-1), v/(mv-1), 0]
                            oml0.computePointsC()

                            P0, Q = PSMlib.computesurfaceprojections(nodes.shape[0], nodes)
                            s,u,v = oml0.evaluateProjection(P0, [surf], Q)

                            B = oml0.evaluateBases(s,u,v)

                            name = comp.name + ':' + str(f) + '-' + str(i) + '-' + str(j)

                            self.surfaceNames.append(name)
                            B0.append(B)
                            nnode0.append(nnode0[-1] + P0.shape[0])
                            quads0.append(quads)

        B0 = scipy.sparse.vstack(B0)

        self.meshS = [B0, quads0, nnode0]

    def computeMembers(self):
        nmem = self.nmem
        geometry = self.geometry
        oml0 = geometry.oml0
        groupIntPtr = self.groupIntPtr
        groupInts = self.groupInts
        groupSplitPtr = self.groupSplitPtr
        groupSplits = self.groupSplits
        quad = self.quad
        ngroup = self.ngroupS + self.ngroupM
        nint = groupIntPtr[-1,-1]
        nsplit = groupSplitPtr[-1,-1]

        nodesInt0 = []
        nodesFlt0 = []
        quads0 = []
        nnode0 = [0]
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

        nodesInt = numpy.array(numpy.vstack(nodesInt0),order='F')
        nodesFlt = numpy.array(numpy.vstack(nodesFlt0),order='F')
        nnode = nodesInt.shape[0]

        for k in range(len(geometry.comps)):
            comp = geometry.comps.values()[k]  
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[comp.name][f]
                PSMlib.computememberlocalcoords(k+1, f+1, ni, nj, nnode, idims, jdims, comp.Ks[f]+1, nodesInt, nodesFlt)

        linW = numpy.linspace(0,nnode-1,nnode)
        B0 = scipy.sparse.csr_matrix((nnode,oml0.C.shape[0]))
        for src in range(4):
            W = scipy.sparse.csr_matrix((nodesFlt[:,src,0],(linW,linW)))
            for surf in range(oml0.nsurf):
                npts = PSMlib.countmembers(surf+1, src+1, nnode, nodesInt)
                if npts is not 0:
                    inds, P, Q = PSMlib.computememberproj(surf+1, src+1, nnode, npts, nodesInt, nodesFlt)
                    Ta = numpy.ones(npts)
                    Ti = inds - 1
                    Tj = numpy.linspace(0,npts-1,npts)
                    T = scipy.sparse.csr_matrix((Ta,(Ti,Tj)),shape=(nnode,npts))

                    mu, mv = oml0.edgeProperty(surf,1)
                    for u in range(mu):
                        for v in range(mv):
                            oml0.C[oml0.getIndex(surf,u,v,1),:3] = [u/(mu-1), v/(mv-1), 0]
                    oml0.computePointsC()

                    s,u,v = oml0.evaluateProjection(P, [surf], Q)
                    B = oml0.evaluateBases(s, u, v)
                    B0 = B0 + W.dot(T.dot(B))

        self.meshM = [B0, quads0, nnode0]
