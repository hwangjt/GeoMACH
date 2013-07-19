from __future__ import division
import numpy, pylab, time, cProfile
from GeoMACH.PSM import PSMlib
from GeoMACH.CDT import CDTlib


class QuadMesh(object):

    def __init__(self, maxL, lines, limits=None):
        if limits==None:
            self.limits = numpy.array([
                    [numpy.amin(lines[:,:,0]), numpy.amin(lines[:,:,1])],
                    [numpy.amax(lines[:,:,0]), numpy.amax(lines[:,:,1])],
                    ], order='F')
        else:
            self.limits = limits
        self.maxL = maxL
        self.verts, self.edges = PSMlib.importedges(2*lines.shape[0], lines.shape[0], lines)

    def mesh(self):
        self.computeIntersections()
        self.deleteDuplicateVerts()
        self.splitEdges()
        self.deleteDuplicateEdges()
        self.computeDivisions()
        self.splitEdges()
        self.computeGrid()
        self.computeCDT()
        self.splitEdges()
        self.deleteDuplicateEdges()
        self.computeAdjMap()
        self.computeTriangles()
        self.deleteDuplicateTriangles()
        self.computeTri2Quad()
        self.deleteDuplicateVerts()
        self.splitEdges()
        self.computeAdjMap()
        self.computeQuads()
        self.deleteDuplicateQuads()

    def computeIntersections(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nint = PSMlib.countintersections(nvert, nedge, self.verts, self.edges)
        self.verts = PSMlib.computeintersections(nvert, nedge, nvert+nint, self.verts, self.edges)

    def deleteDuplicateVerts(self):
        nvert0 = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nvert, ids = PSMlib.computeuniquevertids(nvert0, self.verts)
        self.verts, self.edges = PSMlib.deleteduplicateverts(nvert, nvert0, nedge, ids, self.verts, self.edges)

    def splitEdges(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nsplit = PSMlib.countsplits(nvert, nedge, self.verts, self.edges)
        self.edges = PSMlib.splitedges(nvert, nedge, nedge+nsplit, self.verts, self.edges)

    def deleteDuplicateEdges(self):
        nedge0 = self.edges.shape[0]
        nedge, ids = PSMlib.computeuniqueedgeids(nedge0, self.edges)
        self.edges = PSMlib.deleteduplicateedges(nedge0, nedge, ids, self.edges)

    def computeDivisions(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ndiv = PSMlib.countdivisions(nvert, nedge, self.maxL, self.verts, self.edges)
        self.verts = PSMlib.computedivisions(nvert, nedge, nvert+ndiv, self.maxL, self.verts, self.edges)

    def computeGrid(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ngp = PSMlib.countgridpoints(nvert, self.maxL, self.limits, self.verts)
        self.verts = PSMlib.computegridpoints(nvert, nvert+ngp, self.maxL, self.limits, self.verts)

    def computeCDT(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        self.verts, self.edges = PSMlib.reordercollinear(nvert, nedge, self.verts, self.edges)
        ntri, triangles = CDTlib.computecdt(nvert, nedge, 2*nvert-5, self.verts, self.edges)
        self.edges = PSMlib.trianglestoedges(2*nvert-5, ntri, 3*ntri, triangles)

    def computeAdjMap(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        self.adjPtr, self.adjMap = PSMlib.computeadjmap(nvert, nedge, 2*nedge, self.edges)

    def computeTriangles(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nadj = self.adjMap.shape[0]
        ntri = nedge - nvert + 1
        self.triangles = PSMlib.computetriangles(nvert, nadj, 6*ntri, self.adjPtr, self.adjMap)

    def deleteDuplicateTriangles(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ntri = nedge - nvert + 1
        self.triangles = PSMlib.deleteduplicatetriangles(6*ntri, ntri, self.triangles)

    def computeTri2Quad(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        ntri = nedge - nvert + 1
        self.verts, self.edges = PSMlib.computetri2quad(nvert, nedge, ntri, nvert+4*ntri, nedge+3*ntri, self.verts, self.edges, self.triangles)

    def computeQuads(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nadj = self.adjMap.shape[0]
        nquad = nedge - nvert + 1
        self.quads = PSMlib.computequads(nvert, nadj, 8*nquad, self.adjPtr, self.adjMap)

    def deleteDuplicateQuads(self):
        nvert = self.verts.shape[0]
        nedge = self.edges.shape[0]
        nquad = nedge - nvert + 1
        self.quads = PSMlib.deleteduplicatequads(8*nquad, nquad, self.quads)

    def plot(self, plot, pv=True, pe=True, pt=True, pq=True):
        verts = self.verts
        edges = self.edges

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



if __name__ == "__main__":
    maxL = 0.5
    limits = numpy.array([
            [0,0],
            [3,3],
            ], order='F')

    lines = numpy.zeros((8,2,2),order='F')
    lines[0] = [[0,0],[3,0]]
    lines[1] = [[3,0],[3,3]]
    lines[2] = [[3,3],[0,3]]
    lines[3] = [[0,3],[0,0]]
    lines[4] = [[0,1],[3,0]]
    lines[5] = [[0,1],[3,3]]
    lines[6] = [[1,0],[1,3]]
    lines[7] = [[2,1],[2,2]]

    cProfile.run('QuadMesh(maxL, lines, limits)'); exit
    mesh = QuadMesh(maxL, lines, limits)
    mesh.plot(111,pv=False,pe=False,pt=False)
    pylab.show()
