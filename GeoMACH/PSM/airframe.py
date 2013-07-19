from __future__ import division
import numpy
import time
import scipy.sparse
import pylab

from GeoMACH.PUBS import PUBSlib
from GeoMACH.PSM import PSMlib, QUADlib, CDTlib


class Airframe(object):

    def __init__(self, geometry, members, maxL):
        self.geometry = geometry
        self.maxL = maxL
        self.nmem = len(members)
        self.quad = QUAD()

        self.computePreviewSurfaces()
        self.computeFaceDimensions()
        self.importMembers(members)
        self.computePreviewMembers()

        self.geometry.oml0.export.write2TecQuads('previewSurfs.dat',self.nodesS,self.quadsS-1)
        self.geometry.oml0.export.write2TecQuads('previewMembers.dat',self.nodesM,self.quadsM-1)

    def mesh(self):
        self.computeTopology()
        self.computeAdjoiningEdges()
        self.computeGroupIntersections()
        self.computeFaces()
        self.computeSurfaces()

    def computePreviewSurfaces(self):
        oml0 = self.geometry.oml0
        nsurf = oml0.nsurf

        quadsS,s,u,v = PSMlib.computepreviewsurfaces(4*nsurf,nsurf)
        B = oml0.evaluateBases(s,u,v)

        self.nodesS = B.dot(oml0.C[:,:3])
        self.quadsS = quadsS
        self.surfEdgeLengths = PSMlib.computeedgelengths(4*nsurf,nsurf,self.nodesS,self.quadsS)

    def computeFaceDimensions(self):
        geometry = self.geometry
        nsurf = self.geometry.oml0.nsurf

        faceDims = []
        for k in range(len(geometry.comps)):
            comp = geometry.comps[geometry.keys[k]]
            faceDimsComp = []
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = PSMlib.computefacedimensions(ni,nj,nsurf,comp.Ks[f]+1,self.surfEdgeLengths)
                faceDimsComp.append([idims,jdims])
            faceDims.append(faceDimsComp)

        self.faceDims = faceDims

    def importMembers(self, members):
        geometry = self.geometry
        nmem = self.nmem

        for imem in range(nmem):
            if len(members[imem]) is 2:
                members[imem].extend([[[-1,-1,-1] for j in range(5)] for i in range(2)])
        members = numpy.array(members, order='F')
        membersInt, membersFlt = PSMlib.importmembers(nmem, members)

        for k in range(len(geometry.comps)):
            comp = geometry.comps[geometry.keys[k]]
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[k][f]
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

        self.geometry.computePoints()

        self.nodesM = B0.dot(oml0.C[:,:3])
        self.quadsM = quads
        self.memEdgeLengths = PSMlib.computeedgelengths(4*nmem,nmem,self.nodesM,self.quadsM)

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
        nmem = self.nmem
        ngroup = self.ngroupS + self.ngroupM
        nadj = self.adjoiningInt.shape[0]

        premeshFaces = []
        groupIntCount = numpy.zeros(ngroup,int)
        for k in range(len(geometry.comps)):
            comp = geometry.comps[geometry.keys[k]]
            for f in range(len(comp.Ks)):
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[k][f]
                nedge = PSMlib.countfaceedges(k+1, f+1, ni, nj, nadj, self.adjoiningInt)
                edge_group, edgeLengths, edges = PSMlib.computefaceedges(k+1, f+1, ni, nj, geometry.oml0.nsurf, nmem, nadj, nedge, idims, jdims, comp.Ks[f]+1, self.surf_group, self.mem_group, self.adjoiningInt, self.adjoiningFlt, self.surfEdgeLengths, self.memEdgeLengths)
                verts, edges = quad.importEdges(edges)
                verts, edges = quad.computeIntersections(verts, edges)
                verts, edges = quad.deleteDuplicateVerts(verts, edges)

                groupIntCount = PSMlib.countgroupintersections(verts.shape[0], edges.shape[0], ngroup, verts, edges, edge_group, groupIntCount)

                premeshFaces.append([verts,edges,edge_group,edgeLengths])

        groupIntPtr = PSMlib.computegroupintptr(ngroup, groupIntCount)
        nint = groupIntPtr[-1,-1]

        iList = 0
        groupInts = numpy.zeros(nint)
        for k in range(len(geometry.comps)):
            comp = geometry.comps[geometry.keys[k]]
            for f in range(len(comp.Ks)):
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                iList += 1
                groupInts = PSMlib.computegroupintersections(verts.shape[0], edges.shape[0], ngroup, nint, verts, edges, edge_group, groupIntPtr, groupInts)

        self.premeshFaces = premeshFaces
        self.groupIntPtr = groupIntPtr
        self.groupInts = groupInts

    def computeFaces(self):
        geometry = self.geometry
        quad = self.quad
        ngroup = self.ngroupS + self.ngroupM
        premeshFaces = self.premeshFaces
        groupIntPtr = self.groupIntPtr
        groupInts = self.groupInts
        nint = groupIntPtr[-1,-1]

        iList = 0
        for k in range(len(geometry.comps)):
            comp = geometry.comps[geometry.keys[k]]
            for f in range(len(comp.Ks)):
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                nvert = PSMlib.countintersectionverts(edges.shape[0], ngroup, edge_group, groupIntPtr)
                verts = PSMlib.computeintersectionverts(verts.shape[0], edges.shape[0], ngroup, nint, nvert + verts.shape[0], verts, edges, edge_group, groupIntPtr, groupInts)
                verts, edges = quad.deleteDuplicateVerts(verts, edges)
                verts, edges = quad.splitEdges(verts, edges)
                verts, edges = quad.deleteDuplicateEdges(verts, edges)
                premeshFaces[iList][0] = verts
                premeshFaces[iList][1] = edges
                iList += 1

                if k==1 and f==0 and 0:
                    quad.plot(verts, edges, 111, pt=False, pq=False)
                    pylab.show()

    def computeSurfaces(self):
        geometry = self.geometry
        oml0 = geometry.oml0
        quad = self.quad
        premeshFaces = self.premeshFaces

        B0 = []
        quads0 = []
        nquad0 = 0
        iList = 0
        for k in range(len(geometry.comps)):
            comp = geometry.comps[geometry.keys[k]]
            for f in range(len(comp.Ks)):
                verts,edges,edge_group,edgeLengths = premeshFaces[iList]
                iList += 1
                nvert = verts.shape[0]
                nedge = edges.shape[0]
                ni, nj = comp.Ks[f].shape
                idims, jdims = self.faceDims[k][f]
                for i in range(ni):
                    for j in range(nj):
                        surf = comp.Ks[f][i,j]
                        nedge1 = PSMlib.countsurfaceedges(nvert, nedge, idims[i], idims[i+1], jdims[j], jdims[j+1], verts, edges)
                        edges1 = PSMlib.computesurfaceedges(nvert, nedge, nedge1, idims[i], idims[i+1], jdims[j], jdims[j+1], verts, edges)
                        nodes, quads = quad.mesh(edges1, self.maxL, self.surfEdgeLengths[surf,:,:])

                        mu, mv = oml0.edgeProperty(surf,1)
                        for u in range(mu):
                            for v in range(mv):
                                oml0.C[oml0.getIndex(surf,u,v,1),:3] = [u/(mu-1), v/(mv-1), 0]
                        oml0.computePointsC()

                        P0, Q = PSMlib.computesurfaceprojections(nodes.shape[0], nodes)
                        s,u,v = oml0.evaluateProjection(P0, [surf], Q)
                        B0.append(oml0.evaluateBases(s,u,v))
                        quads0.append(quads + nquad0)
                        nquad0 += P0.shape[0]

                        if k==1 and f==0 and i==0 and j==3 and 0:
                            quad.plot(verts1, edges1, 111, pt=False, pq=False)
                            pylab.show()
                            exit()

        geometry.computePoints()
        P = scipy.sparse.vstack(B0).dot(oml0.C[:,:3])
        quads0 = numpy.vstack(quads0)

        self.geometry.oml0.export.write2TecQuads('structure.dat', P, quads0-1, loop=True)
        print P.shape[0]


class QUAD(object):

    def mesh(self, lines, maxL, lengths):
        verts, edges = self.importEdges(lines)
        verts, edges = self.computeDivisions(verts, edges, maxL, lengths)
        verts, edges = self.computeGrid(verts, edges, maxL, lengths)
        verts, edges = self.deleteDuplicateVerts(verts, edges)
        verts, edges = self.splitEdges(verts, edges)
        verts[:,0] *= 0.5*numpy.sum(lengths[0,:])
        verts[:,1] *= 0.5*numpy.sum(lengths[1,:])
        verts, edges = self.computeCDT(verts, edges)
        verts, edges = self.splitEdges(verts, edges)
        verts, edges = self.deleteDuplicateEdges(verts, edges)
        adjPtr, adjMap = self.computeAdjMap(verts, edges)
        tris = self.computeTriangles(verts, edges, adjPtr, adjMap)
        tris = self.deleteDuplicateTriangles(verts, edges, tris)
        verts, edges = self.computeTri2Quad(verts, edges, tris)
        verts, edges = self.deleteDuplicateVerts(verts, edges)
        verts, edges = self.splitEdges(verts, edges)
        adjPtr, adjMap = self.computeAdjMap(verts, edges)
        quads = self.computeQuads(verts, edges, adjPtr, adjMap)
        quads = self.deleteDuplicateQuads(verts, edges, quads)
        return verts, quads

    def importEdges(self, lines):
        verts, edges = QUADlib.importedges(2*lines.shape[0], lines.shape[0], lines)
        return verts, edges

    def computeIntersections(self, verts, edges):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        nint = QUADlib.countintersections(nvert, nedge, verts, edges)
        verts = QUADlib.computeintersections(nvert, nedge, nvert+nint, verts, edges)
        return verts, edges

    def deleteDuplicateVerts(self, verts, edges):
        nvert0 = verts.shape[0]
        nedge = edges.shape[0]
        nvert, ids = QUADlib.computeuniquevertids(nvert0, verts)
        verts, edges = QUADlib.deleteduplicateverts(nvert, nvert0, nedge, ids, verts, edges)
        return verts, edges

    def splitEdges(self, verts, edges):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        nsplit = QUADlib.countsplits(nvert, nedge, verts, edges)
        edges = QUADlib.splitedges(nvert, nedge, nedge+nsplit, verts, edges)
        return verts, edges

    def deleteDuplicateEdges(self, verts, edges):
        nedge0 = edges.shape[0]
        nedge, ids = QUADlib.computeuniqueedgeids(nedge0, edges)
        edges = QUADlib.deleteduplicateedges(nedge0, nedge, ids, edges)
        return verts, edges

    def computeDivisions(self, verts, edges, maxL, lengths):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        ndiv = QUADlib.countdivisions(nvert, nedge, maxL, lengths, verts, edges)
        verts = QUADlib.computedivisions(nvert, nedge, nvert+ndiv, maxL, lengths, verts, edges)
        return verts, edges

    def computeGrid(self, verts, edges, maxL, lengths):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        ngp = QUADlib.countgridpoints(maxL, lengths)
        verts = QUADlib.computegridpoints(nvert, nvert+ngp, maxL, lengths, verts)
        return verts, edges

    def computeCDT(self, verts, edges):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        verts, edges = QUADlib.reordercollinear(nvert, nedge, verts, edges)
        ntri, triangles = CDTlib.computecdt(nvert, nedge, 2*nvert-5, verts, edges)
        edges = QUADlib.trianglestoedges(2*nvert-5, ntri, 3*ntri, triangles)
        return verts, edges

    def computeAdjMap(self, verts, edges):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        adjPtr, adjMap = QUADlib.computeadjmap(nvert, nedge, 2*nedge, edges)
        return adjPtr, adjMap

    def computeTriangles(self, verts, edges, adjPtr, adjMap):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        nadj = adjMap.shape[0]
        ntri = nedge - nvert + 1
        triangles = QUADlib.computetriangles(nvert, nadj, 6*ntri, adjPtr, adjMap)
        return triangles

    def deleteDuplicateTriangles(self, verts, edges, triangles):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        ntri = nedge - nvert + 1
        triangles = QUADlib.deleteduplicatetriangles(6*ntri, ntri, triangles)
        return triangles

    def computeTri2Quad(self, verts, edges, triangles):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        ntri = nedge - nvert + 1
        verts, edges = QUADlib.computetri2quad(nvert, nedge, ntri, nvert+4*ntri, nedge+3*ntri, verts, edges, triangles)
        return verts, edges

    def computeQuads(self, verts, edges, adjPtr, adjMap):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        nadj = adjMap.shape[0]
        nquad = nedge - nvert + 1
        quads = QUADlib.computequads(nvert, nadj, 8*nquad, adjPtr, adjMap)
        return quads

    def deleteDuplicateQuads(self, verts, edges, quads):
        nvert = verts.shape[0]
        nedge = edges.shape[0]
        nquad = nedge - nvert + 1
        quads = QUADlib.deleteduplicatequads(8*nquad, nquad, quads)
        return quads

    def plot(self, verts, edges, plot, pv=True, pe=True, pt=True, pq=True):

        pylab.subplot(plot)
        pylab.axis('equal')
        if pv:
            for v in range(verts.shape[0]):
                pylab.plot([verts[v,0]],[verts[v,1]],'ok')
        if pe:
            for e in range(edges.shape[0]):
                pylab.plot(
                    [verts[edges[e,i]-1,0] for i in range(2)],
                    [verts[edges[e,i]-1,1] for i in range(2)],
                    )
        if pt:
            triangles = self.triangles
            for t in range(triangles.shape[0]):
                pylab.plot(
                    [verts[triangles[t,i]-1,0] for i in range(2)],
                    [verts[triangles[t,i]-1,1] for i in range(2)],
                    )
                pylab.plot(
                    [verts[triangles[t,i+1]-1,0] for i in range(2)],
                    [verts[triangles[t,i+1]-1,1] for i in range(2)],
                    )
                pylab.plot(
                    [verts[triangles[t,2*i]-1,0] for i in range(2)],
                    [verts[triangles[t,2*i]-1,1] for i in range(2)],
                    )
        if pq:
            quads = self.quads
            for q in range(quads.shape[0]):
                pylab.plot(
                    [verts[quads[q,i]-1,0] for i in range(2)],
                    [verts[quads[q,i]-1,1] for i in range(2)],
                    )
                pylab.plot(
                    [verts[quads[q,i+1]-1,0] for i in range(2)],
                    [verts[quads[q,i+1]-1,1] for i in range(2)],
                    )
                pylab.plot(
                    [verts[quads[q,i+2]-1,0] for i in range(2)],
                    [verts[quads[q,i+2]-1,1] for i in range(2)],
                    )
                pylab.plot(
                    [verts[quads[q,3*i]-1,0] for i in range(2)],
                    [verts[quads[q,3*i]-1,1] for i in range(2)],
                    )
