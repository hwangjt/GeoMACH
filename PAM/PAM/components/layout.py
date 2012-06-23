from __future__ import division
import numpy, pylab


class Layout(object):

    def __init__(self):
        self.Lx = 1
        self.Ly = 1

    def structured(self, nu, nv, u1=0, u2=1, v1=0, v2=1, edges=[]):
        u = numpy.linspace(u1,u2,nu)
        v = numpy.linspace(v1,v2,nv)
        for i in range(nu):
            edges.append([[u[i],v1],[u[i],v2]])
        for j in range(nv):
            edges.append([[u1,v[j]],[u2,v[j]]])
        self.initialize(edges)

    def initialize(self, edges):
        self.edge_vert = []
        self.verts = []

        self.importEdges(numpy.array(edges))
        self.cleanup()

        done = False
        while not done:
            self.addConnectors()
            self.cleanup()
            self.computePolygons()
            done = True
            for p in range(len(self.poly_edge)):
                if len(self.poly_edge[p]) > 4:
                    done = False

        self.computeGroups()
        self.splitTriangles()
        self.propagateSplits()
        self.cleanup()
        self.computeSurfaces()

    def cleanup(self):
        self.computeIntersections()
        self.splitEdges()
        self.deleteDuplicates()

    def plot(self):
        self.computePolygons()
        print len(self.verts)
        print len(self.edge_vert)
        print len(self.poly_edge)
        for e in range(len(self.edge_vert)):
            P0 = self.verts[self.edge_vert[e][0]]
            P1 = self.verts[self.edge_vert[e][1]]
            pylab.plot([P0[0],P1[0]],[P0[1],P1[1]],'k')
        for v in range(len(self.verts)):
            P = self.verts[v]
            pylab.plot([P[0]],[P[1]],'ko')
        pylab.show()

    def importEdges(self, edges):
        self.addEdge([0,0],[0,1])
        self.addEdge([0,0],[1,0])
        self.addEdge([1,1],[0,1])
        self.addEdge([1,1],[1,0])
        for e in range(edges.shape[0]):
            self.addEdge(edges[e,0], edges[e,1]) 

    def addVert(self, P0):
        verts = self.verts
        P = numpy.array(P0)
        for v in range(len(verts)):
            if numpy.linalg.norm(verts[v]-P) < 1e-14:
                return v
        verts.append(P)
        return len(verts) - 1 

    def addEdge(self, P0, P1):      
        v0 = self.addVert(P0)
        v1 = self.addVert(P1)
        self.edge_vert.append([v0,v1])

    def computeIntersections(self):
        def contained(t):
            return 0 < t and t < 1

        verts = self.verts
        edge_vert = self.edge_vert
        nedge = len(edge_vert)
        for e1 in range(nedge):
            for e2 in range(e1+1,nedge):
                v0,v1 = edge_vert[e1]
                a1,b1 = verts[v0]
                m1,n1 = verts[v1] - verts[v0]
                v0,v1 = edge_vert[e2]
                a2,b2 = verts[v0]
                m2,n2 = verts[v1] - verts[v0]
                det = n1*m2 - n2*m1
                if not det==0:
                    da = a2 - a1
                    db = b2 - b1
                    t1 = 1/det*(-n2*da + m2*db)
                    t2 = 1/det*(-n1*da + m1*db)
                    if contained(t1) and contained(t2):
                        self.addVert([a1+t1*m1,b1+t1*n1])

    def splitEdges(self):
        verts = self.verts
        edge_vert = self.edge_vert
        e = 0
        while e < len(edge_vert):
            v = 0
            while v < len(verts):
                v0 = edge_vert[e][0]
                v1 = edge_vert[e][1]
                if (not v==v0) and (not v==v1):
                    r0 = (verts[v]-verts[v0])/numpy.linalg.norm(verts[v]-verts[v0])
                    r1 = (verts[v]-verts[v1])/numpy.linalg.norm(verts[v]-verts[v1])
                    if abs(numpy.dot(r0,r1)+1) < 1e-14:
                        v0,v1 = self.edge_vert[e]
                        edge_vert.append([v0,v])
                        edge_vert.append([v,v1])
                        edge_vert.pop(e)
                        v = -1
                v += 1
            e += 1

    def deleteDuplicates(self):
        edge_vert = self.edge_vert
        e1 = 0
        while e1 < len(edge_vert):
            e2 = e1 + 1
            while e2 < len(edge_vert):
                edge1 = edge_vert[e1]
                edge2 = edge_vert[e2]
                if edge1==edge2 or edge1==edge2[::-1]:
                    edge_vert.pop(e2)
                else:
                    e2 += 1
            e1 += 1

    def arctan(self, P):
        x,y = P
        x *= self.Lx
        y *= self.Ly
        if x==0:
            if y > 0:
                t = numpy.pi/2.0
            elif y < 0:
                t = 3*numpy.pi/2.0
        elif y==0:
            if x > 0:
                t = 0
            elif x < 0:
                t = numpy.pi
        elif x<0:
            t = numpy.arctan(y/x) + numpy.pi
        elif y<0:
            t = numpy.arctan(y/x) + 2*numpy.pi
        elif y>0:
            t = numpy.arctan(y/x)
        else:
            t = 0
        return t 

    def addConnectors(self):
        def isContained(t, t1, t2):
            return t1 <= t and t <= t2

        verts = self.verts
        edge_vert = self.edge_vert
        quadrants = numpy.zeros((len(verts),4),bool)
        for e in range(len(edge_vert)):
            for d in range(2):
                v = edge_vert[e][d]
                P1 = verts[edge_vert[e][-d]]
                P2 = verts[edge_vert[e][1-d]]
                t = self.arctan(P2-P1)/numpy.pi
                if isContained(t,0,0.25):
                    quadrants[v,0] = True
                elif isContained(t,0.25,0.75):
                    quadrants[v,1] = True
                elif isContained(t,0.75,1.25):
                    quadrants[v,2] = True
                elif isContained(t,1.25,1.75):
                    quadrants[v,3] = True
                elif isContained(t,1.75,2):
                    quadrants[v,0] = True
        for v in range(len(verts)):
            v0 = verts[v][0]
            v1 = verts[v][1]
            if not (v0==0 or v0==1 or v1==0 or v1==1):
                if not quadrants[v,0]:
                    self.addEdge(verts[v],[1,v1])
                if not quadrants[v,1]:
                    self.addEdge(verts[v],[v0,1])
                if not quadrants[v,2]:
                    self.addEdge(verts[v],[0,v1])
                if not quadrants[v,3]:
                    self.addEdge(verts[v],[v0,0])

    def computePolygons(self):
        def getOtherV(edge_vert, v, e):
            if edge_vert[e][0] == v:
                return edge_vert[e][1]
            else:
                return edge_vert[e][0]

        def turnRight(self, vert_edge, v, e):
            verts = self.verts
            edge_vert = self.edge_vert
            v.append(getOtherV(edge_vert,v[-1],e[-1]))

            t0 = self.arctan(verts[v[-2]]-verts[v[-1]])/numpy.pi
            t = []
            for ee in vert_edge[v[-1]]:
                vv = getOtherV(edge_vert, v[-1], ee)
                t.append(self.arctan(verts[vv]-verts[v[-1]])/numpy.pi)
                if t[-1] <= t0:
                    t[-1] += 2
            index = t.index(min(t))
            if t[index] - t0  < 1:
                e.append(vert_edge[v[-1]][index])
                return False
            else:
                return True                     

        verts = self.verts
        edge_vert = self.edge_vert

        vert_edge = []
        for v in range(len(verts)):
            vert_edge.append([])
        
        for e in range(len(edge_vert)):
            for d in range(2):
                vert_edge[edge_vert[e][d]].append(e)

        self.poly_vert = []
        self.poly_edge = []
        for v0 in range(len(verts)):
            for e0 in vert_edge[v0]:
                e = [e0]
                v = [v0]
                done = False
                while not done:
                    done = turnRight(self, vert_edge, v, e)
                    if (not done) and getOtherV(edge_vert, v[-1], e[-1]) == v0:
                        self.addPoly(v,e)
                        done = True

    def addPoly(self, v, e):
        poly_vert = self.poly_vert
        poly_edge = self.poly_edge
        found = False
        for p in range(len(poly_edge)):
            if len(poly_edge[p])==len(e):
                found = True
                for i in e:
                    if not poly_edge[p].count(i)==1:
                        found = False
                if found:
                    break
        if not found:
            poly_vert.append(v)
            poly_edge.append(e)

    def computeGroups(self):
        def setAll(edge_group, old, new):
            found = False
            for e in range(len(edge_group)):
                if edge_group[e] == old:
                    edge_group[e] = new
                    found = True
            return found

        nedge = len(self.edge_vert)
        edge_group = numpy.linspace(0,nedge-1,nedge)
        poly_edge = self.poly_edge
        for p in range(len(poly_edge)):
            if len(poly_edge[p])==4:
                setAll(edge_group, edge_group[poly_edge[p][0]], edge_group[poly_edge[p][2]])
                setAll(edge_group, edge_group[poly_edge[p][1]], edge_group[poly_edge[p][3]])
        group = 0
        for e in range(len(edge_group)):
            found = setAll(edge_group, e, group)
            if found:
                group += 1
        self.edge_group = edge_group
        self.group_split = numpy.zeros(max(self.edge_group)+1,bool)

    def splitPentagons(self):
        group_split = self.group_split
        verts = self.verts
        edge_group = self.edge_group
        edge_vert = self.edge_vert
        poly_vert = self.poly_vert
        poly_edge = self.poly_edge
        for p in range(len(poly_edge)):
            if len(poly_edge[p])==5:
                L = numpy.zeros(5)
                for e in range(5):
                    C0 = verts[edge_vert[poly_edge[p][e]][0]]
                    C1 = verts[edge_vert[poly_edge[p][e]][1]]
                    L[e] = numpy.linalg.norm(C1-C0)
                e = numpy.argmax(L)
                C0 = verts[edge_vert[poly_edge[p][e]][0]]
                C1 = verts[edge_vert[poly_edge[p][e]][1]]
                C = verts[poly_vert[p][e-2]]
                group_split[edge_group[poly_edge[p][e]]] = True
                self.addEdge(C,0.5*C0+0.5*C1)

    def splitTriangles(self):
        group_split = self.group_split
        verts = self.verts
        edge_group = self.edge_group
        edge_vert = self.edge_vert
        poly_vert = self.poly_vert
        poly_edge = self.poly_edge
        for p in range(len(poly_edge)):
            if len(poly_edge[p])==3:
                C = []
                for k in range(3):
                    C.append(verts[poly_vert[p][k]])
                ctd = 1/3.0*C[0] + 1/3.0*C[1] + 1/3.0*C[2]
                mid = []
                mid.append(1/2.0*C[1] + 1/2.0*C[2])
                mid.append(1/2.0*C[0] + 1/2.0*C[2])
                mid.append(1/2.0*C[0] + 1/2.0*C[1])
                for k in range(3):
                    group_split[edge_group[poly_edge[p][k]]] = True
                    self.addEdge(ctd,mid[k])

    def propagateSplits(self):
        group_split = self.group_split
        verts = self.verts
        edge_group = self.edge_group
        edge_vert = self.edge_vert
        poly_vert = self.poly_vert
        poly_edge = self.poly_edge
        for p in range(len(poly_edge)):
            if len(poly_edge[p])==4:
                C = []
                for k in range(4):
                    C.append(verts[poly_vert[p][k]])
                if group_split[edge_group[poly_edge[p][0]]]:
                    self.addEdge(0.5*C[0]+0.5*C[1],0.5*C[2]+0.5*C[3])
                if group_split[edge_group[poly_edge[p][1]]]:
                    self.addEdge(0.5*C[0]+0.5*C[3],0.5*C[2]+0.5*C[1])

    def computeSurfaces(self):
        def getSurface(verts):
            v0 = 1

        poly_vert = self.poly_vert
        surfs = []
        for p in range(len(poly_vert)):
            surfs.append(getSurface(poly_vert[p]))
        self.surfs = surfs
                    


if __name__ == '__main__':

    edges = []
    if 0:
        edges.append([[0.1,0.1],[0.9,0.9]])
        edges.append([[0.9,0.1],[0.1,0.9]])
    if 0:
        edges.append([[0.1,0.2],[0.9,0.9]])
        edges.append([[0.6,0.3],[0.2,0.4]])
        edges.append([[0.3,0.5],[0.3,0.2]])
        edges.append([[0.8,0.8],[0.1,0.4]])
    if 1:
        edges.append([[0,0],[0.3,0.3]])
        edges.append([[0,1],[0.3,0.7]])
        edges.append([[0,0.1],[0.1,0.1]])
        edges.append([[0,0.2],[0.2,0.2]])
        edges.append([[0,0.3],[0.3,0.3]])
        edges.append([[0,0.7],[0.1,0.7]])
        edges.append([[0,0.8],[0.2,0.8]])
        edges.append([[0,0.9],[0.3,0.9]])
        edges.append([[0.9,0.3],[0.3,0.3]])
        edges.append([[0.9,0.7],[0.3,0.7]])
        edges.append([[0,0.1],[0.1,0.1]])
        u = numpy.linspace(0,1,8)
        for i in range(8):
            edges.append([[u[i],0],[u[i],1]])
        
    l = Layout()
    l.initialize(edges)
#    l.structured(11,11)
    l.plot()
