from __future__ import division
import numpy
import time
import scipy.sparse
#import pylab

from GeoMACH.PSM import PSMlib, QUADlib, CDTlib, BLSlib


class QUAD(object):

    def __init__(self):
        self.output = False

    def importEdges(self, lines):
        self.verts, self.edges = QUADlib.importedges(2*lines.shape[0], lines.shape[0], lines)
        self.edgeCon = numpy.ones(self.edges.shape[0], bool)
        if self.output: print 'Done: importEdges'

    def importVertsNEdges(self, verts, edges):
        self.verts, self.edges = verts, edges
        self.edgeCon = numpy.ones(self.edges.shape[0], bool)
        if self.output: print 'Done: importVertsNEdges'

    def mesh(self, maxL, lengths, plot=False, output=False):
        self.maxL = maxL
        self.lengths = lengths
        self.output = output

        self.addIntersectionPts()
        self.removeDuplicateVerts()
        self.splitEdges()
        self.removeDuplicateEdges()

        #self.addEdgePts()
        #self.splitEdges()
        self.addInteriorPts()
        self.removeDuplicateVerts()

        self.reorderCollinear()
        self.edges0 = self.edges
        self.verts[:,0] *= 0.5*numpy.sum(self.lengths[0,:])
        self.verts[:,1] *= 0.5*numpy.sum(self.lengths[1,:])
        self.removeDegenerateEdges()
        self.computeCDT()
        self.splitEdges()
        self.removeDuplicateEdges()
        self.computeConstrainedEdges()

        if plot:
            print self.verts.shape[0]
            print self.edges.shape[0]
            self.plot(111,pv=True,pe=True)
            import pylab
            pylab.show()
        self.computeAdjMap()
        self.computeTriangles()
        self.computeQuadDominant()
        self.splitTrisNQuads()
        self.removeDuplicateVerts()
        self.splitEdges()
        self.computeAdjMap()
        self.computeQuads()
        self.computeQuad2Edge()
        self.computeConstrainedVerts()
        self.smooth1()

        return self.verts, self.quads

    def addIntersectionPts(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nint = QUADlib.countintersectionpts(nvert, nedge, self.verts, self.edges)
        self.verts = QUADlib.addintersectionpts(nvert, nedge, nvert+nint, self.verts, self.edges)
        if self.output: print 'Done: addIntersectionPts'

    def removeDegenerateEdges(self):
        nedge = self.edges.shape[0]
        ndeg = QUADlib.countdegenerateedges(nedge, self.edges)
        self.edges = QUADlib.removedegenerateedges(nedge, nedge-ndeg, self.edges)
        if self.output: print 'Done: removeDegenerateEdges'

    def removeDuplicateVerts(self):
        nvert0 = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nvert, ids = QUADlib.computeuniqueverts(nvert0, self.verts)
        self.verts, self.edges = QUADlib.removeduplicateverts(nvert0, nvert, nedge, ids, self.verts, self.edges)
        if self.output: print 'Done: removeDuplicateVerts'

    def splitEdges(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nsplit = QUADlib.countsplits(nvert, nedge, self.verts, self.edges)
        self.edges, self.edgeCon = QUADlib.splitedges(nvert, nedge, nedge+nsplit, self.verts, self.edges, self.edgeCon)
        if self.output: print 'Done: splitEdges'

    def removeDuplicateEdges(self):
        nedge0 = self.edges.shape[0]
        nedge, ids = QUADlib.computeuniqueedges(nedge0, self.edges)
        self.edges, self.edgeCon = QUADlib.removeduplicateedges(nedge0, nedge, ids, self.edges, self.edgeCon)
        if self.output: print 'Done: removeDuplicateEdges'

    def addEdgePts(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        npts = QUADlib.countedgepts(nvert, nedge, self.maxL, self.lengths, self.verts, self.edges)
        self.verts = QUADlib.addedgepts(nvert, nedge, nvert+npts, self.maxL, self.lengths, self.verts, self.edges)
        if self.output: print 'Done: addEdgePts'

    def addInteriorPts(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        npts = QUADlib.countinteriorpts(self.maxL, self.lengths)
        self.verts = QUADlib.addinteriorpts(nvert, nvert+npts, self.maxL, self.lengths, self.verts)
        if self.output: print 'Done: addInteriorPts'

    def reorderCollinear(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        self.verts, self.edges = QUADlib.reordercollinear(nvert, nedge, self.verts, self.edges)
        if self.output: print 'Done: reorderCollinear'

    def computeCDT(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ntri, triangles = CDTlib.computecdt(nvert, nedge, 2*nvert-5, self.verts, self.edges)
        self.edges = QUADlib.importtriangles(2*nvert-5, ntri, 3*ntri, triangles)
        self.edgeCon = numpy.ones(self.edges.shape[0], bool)
        if self.output: print 'Done: computeCDT'

    def computeConstrainedEdges(self):
        nedge = self.edges.shape[0]
        nedge0 = self.edges0.shape[0]
        self.edgeCon = QUADlib.computeconstrainededges(nedge0, nedge, self.edges0, self.edges)
        if self.output: print 'Done: computeConstrainedEdges'

    def computeAdjMap(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        self.adjPtr, self.adjMap = QUADlib.computeadjmap(nvert, nedge, 2*nedge, self.edges)
        if self.output: print 'Done: computeAdjMap'

    def computeTriangles(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nadj = self.adjMap.shape[0]
        ntri = nedge - nvert + 1
        ntri = int(QUADlib.counttriangles(nvert, nadj, self.adjPtr, self.adjMap)/6)
        self.triangles = QUADlib.computetriangles(nvert, nadj, 6*ntri, self.adjPtr, self.adjMap)
        self.triangles = QUADlib.removeduplicatetriangles(6*ntri, ntri, self.triangles)
        QUADlib.rotatetriangles(nvert, ntri, self.verts, self.triangles)
        self.edge2tri, self.tri2edge = QUADlib.computetri2edge(nedge, ntri, self.edges, self.triangles)
        if self.output: print 'Done: computeTriangles'

    def computeQuadDominant(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ntri = self.triangles.shape[0]
        nrem, removeEdge, removeTri = QUADlib.computequaddominant(nvert, nedge, ntri, self.verts, self.edges, self.edgeCon, self.triangles, self.tri2edge, self.edge2tri)
        self.edgeCon, self.edges, self.triangles, self.quads = QUADlib.mergetriangles(nedge, ntri, nedge-nrem, ntri-2*nrem+1, nrem+1, self.edgeCon, self.edges, self.triangles, self.edge2tri, removeEdge, removeTri)
        if self.output: print 'Done: computeQuadDominant'

    def splitTrisNQuads(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ntri = self.triangles.shape[0] - 1
        nquad = self.quads.shape[0] - 1
        self.verts, self.edges, self.edgeCon = QUADlib.splittrisnquads(nvert, nedge, ntri+1, nquad+1, nvert+4*ntri+5*nquad, nedge+3*ntri+4*nquad, self.verts, self.edges, self.edgeCon, self.triangles, self.quads)
        if self.output: print 'Done: splitTrisNQuads'

    def computeQuads(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nadj = self.adjMap.shape[0]
        nquad = nedge - nvert + 1
        nquad = int(QUADlib.countquads(nvert, nadj, self.adjPtr, self.adjMap)/8)
        self.quads = QUADlib.computequads(nvert, nadj, 8*nquad, self.adjPtr, self.adjMap)
        self.quads = QUADlib.removeduplicatequads(8*nquad, nquad, self.quads)
        QUADlib.rotatequads(nvert, nquad, self.verts, self.quads)
        if self.output: print 'Done: computeQuads'

    def computeQuad2Edge(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nquad = self.quads.shape[0]
        self.quad2edge = QUADlib.computequad2edge(nedge, nquad, self.edges, self.quads)
        if self.output: print 'Done: computeQuad2Edge'

    def computeConstrainedVerts(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        self.vertCon = QUADlib.computeconstrainedverts(nvert, nedge, self.edges, self.edgeCon)
        if self.output: print 'Done: computeConstrainedVerts'

    def smooth1(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nquad = self.quads.shape[0]
        nM = nquad*16
        nC = 2 * numpy.sum(self.vertCon)
        nB = nvert + numpy.sum(self.vertCon)
        stack = numpy.hstack
        Ma, Mi, Mj = BLSlib.assemblemtx1(nquad, nM, self.quads)
        Ca, Ci, Cj = BLSlib.assembleconmtx1(nvert, nC, self.vertCon)
        Aa, Ai, Aj = stack([Ma,Ca]), stack([Mi,Ci]), stack([Mj,Cj])
        A = scipy.sparse.csc_matrix((Aa,(Ai-1,Aj-1)))
        B = BLSlib.assemblerhs1(nvert, nB, self.verts, self.vertCon)
        sol = scipy.sparse.linalg.factorized(A)
        for i in range(2):
            self.verts[:,i] = sol(B[:,i])[:nvert]
        if self.output: print 'Done: smooth1'

    def smooth2(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nquad = self.quads.shape[0]
        nM = nquad*81
        nC = 2 * (numpy.sum(self.vertCon) + numpy.sum(self.edgeCon))
        nB = nvert + nedge + nquad + numpy.sum(self.vertCon) + numpy.sum(self.edgeCon)
        stack = numpy.hstack
        Ma, Mi, Mj = BLSlib.assemblemtx2(nvert, nedge, nquad, nM, self.quads, self.quad2edge)
        Ca, Ci, Cj = BLSlib.assembleconmtx2(nvert, nedge, nquad, nC, self.vertCon, self.edgeCon)
        Aa, Ai, Aj = stack([Ma,Ca]), stack([Mi,Ci]), stack([Mj,Cj])
        A = scipy.sparse.csc_matrix((Aa,(Ai-1,Aj-1)))
        B = BLSlib.assemblerhs2(nvert, nedge, nB, self.verts, self.edges, self.vertCon, self.edgeCon)
        sol = scipy.sparse.linalg.factorized(A)
        for i in range(2):
            self.verts[:,i] = sol(B[:,i])[:nvert]
        if self.output: print 'Done: smooth2'

    def plot(self, plot, pv=False, pe=False, pt=False, pq=False, pv2=False, pe2=False):
        import pylab

        verts = self.verts
        def draw(v1, v2, clr='k'):
            pylab.plot([verts[v1-1,0], verts[v2-1,0]],
                       [verts[v1-1,1], verts[v2-1,1]],
                       clr)

        pylab.subplot(plot)
        pylab.axis('equal')
        if pe2:
            edges = self.edges
            for e in range(edges.shape[0]):
                draw(edges[e,0], edges[e,1], 'r' if self.edgeCon[e] else 'k')
        if pv2:
            for v in range(verts.shape[0]):
                pylab.plot([verts[v,0]],[verts[v,1]], 'or' if self.vertCon[v] else 'ok')
        if pe:
            edges = self.edges
            for e in range(edges.shape[0]):
                draw(edges[e,0], edges[e,1], 'k')
        if pv:
            for v in range(verts.shape[0]):
                pylab.plot([verts[v,0]],[verts[v,1]], 'ok')
        if pt:
            triangles = self.triangles
            for t in range(triangles.shape[0]):
                draw(triangles[t,0], triangles[t,1])
                draw(triangles[t,1], triangles[t,2])
                draw(triangles[t,0], triangles[t,2])
        if pq:
            quads = self.quads
            for q in range(quads.shape[0]):
                draw(quads[q,0], quads[q,1])
                draw(quads[q,1], quads[q,2])
                draw(quads[q,2], quads[q,3])
                draw(quads[q,0], quads[q,3])
