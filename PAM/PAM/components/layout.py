from __future__ import division
import numpy, pylab, time
import PAM.PAMlib as PAMlib


class Member(object):

    def __init__(self, mat, nedge, SP1, EP1, SP2, EP2):
        self.mat = mat
        self.nedge = nedge
        self.SP1 = SP1
        self.EP1 = EP1
        self.SP2 = SP2
        self.EP2 = EP2
        self.types = [[0.0,'None'],[1.0]]

    def setTypes(self, mtype, limits=[0,1], n=1):
        if self.types == [[0.0,'None'],[1.0]]:
            self.types = [[1.0]]
        pos = numpy.linspace(limits[0],limits[1],n+1)
        for i in range(n):
            if mtype=='None':
                self.types.insert(-1,[pos[i],mtype])
            elif mtype=='Solid':
                self.types.insert(-1,[pos[i],mtype])
            elif mtype=='Stringer':
                self.types.insert(-1,[pos[i],mtype,0.1])
            elif mtype=='vStringer':
                self.types.insert(-1,[pos[i],mtype,0.1])
            elif mtype=='hStringer':
                self.types.insert(-1,[pos[i],mtype,0.1])
            elif mtype=='LErib':
                self.types.insert(-1,[pos[i],mtype,0.1])
            elif mtype=='TErib':
                self.types.insert(-1,[pos[i],mtype,0.1])
            elif mtype=='Hole':
                self.types.insert(-1,[pos[i],mtype,0.2,0.2,0.1,1.0])


class Layout(object):

    def __init__(self, Lx=1, Ly=1):
        self.Lx = Lx
        self.Ly = Ly
        self.members = {}
        self.keys = []

    def addMembers(self, name, mat, nedge, SP1=[0,0], EP1=[0,0], SP2=[0,0], EP2=[0,0]):
        self.members[name] = Member(mat, nedge, SP1, EP1, SP2, EP2)
        self.keys.append(name)

    def initialize(self):
        self.importEdges()
        self.computeIntersections()
        self.addConnectors()
        self.computeIntersections()
        print self.checkForPentagons()
        while self.checkForPentagons():
            self.addConnectors()
            self.computeIntersections()
            print 'a'
            

    def importEdges(self):
        def getv(r, v1, v2):
            return numpy.array(v1)*(1-r) + numpy.array(v2)*r

        self.addMembers('hWall', 1, 2, SP1=[0,0], EP1=[0,1], SP2=[1,0], EP2=[1,1])
        self.addMembers('vWall', 1, 2, SP1=[0,0], EP1=[1,0], SP2=[0,1], EP2=[1,1])

        nTotal = 0
        for m in self.keys:
            member = self.members[m]
            nedge = member.nedge
            ntypes = len(member.types)-1
            nTotal += nedge*ntypes
        edges = numpy.zeros((nTotal,4))
        segments = numpy.zeros((nTotal,3))
        index = 0
        for m in range(len(self.keys)):
            member = self.members[self.keys[m]]
            nedge = member.nedge
            ntypes = len(member.types)-1
            for i in range(nedge):
                SP = getv(i/(nedge-1), member.SP1, member.SP2)
                EP = getv(i/(nedge-1), member.EP1, member.EP2)
                for j in range(ntypes):
                    edges[index,0:2] = getv(member.types[j][0], SP, EP)
                    edges[index,2:4] = getv(member.types[j+1][0], SP, EP)
                    segments[index,:] = [m,i,j]
                    index += 1
        
        self.segments = segments
        self.edges = PAMlib.getedges(nTotal, edges)
        self.verts = PAMlib.getverts(numpy.max(self.edges[:,:2]), nTotal, edges, self.edges)
        self.nedge = self.edges.shape[0]
        self.nvert = self.verts.shape[0]

    def computeIntersections(self):
        nint = PAMlib.numintersections(self.nvert, self.nedge, self.verts, self.edges)
        self.verts = PAMlib.addintersections(self.nvert + nint, self.nvert, self.nedge, self.verts, self.edges)
        self.nvert = self.verts.shape[0]

        ndup = PAMlib.countduplicateverts(self.nvert, self.verts)
        self.verts, self.edges = PAMlib.deleteduplicateverts(self.nvert - ndup, self.nvert, self.nedge, self.verts, self.edges)
        self.nvert = self.verts.shape[0]
        self.nedge = self.edges.shape[0]

        nsplit = PAMlib.numsplits(self.nvert, self.nedge, self.verts, self.edges)
        self.edges = PAMlib.splitedges(self.nedge + nsplit, self.nvert, self.nedge, self.verts, self.edges)
        self.nedge = self.edges.shape[0]

        ndup = PAMlib.countduplicateedges(self.nedge, self.edges)
        self.edges = PAMlib.deleteduplicateedges(self.nedge - ndup, self.nedge, self.edges)
        self.nedge = self.edges.shape[0]

    def addConnectors(self):
        ncon, quadrants = PAMlib.countconnectors(self.nvert, self.nedge, self.Lx, self.Ly, self.verts, self.edges)
        self.verts, self.edges = PAMlib.addconnectors(self.nvert + ncon, self.nedge + ncon, self.nvert, self.nedge, self.verts, self.edges, quadrants)
        self.nvert = self.verts.shape[0]
        self.nedge = self.edges.shape[0]

    def checkForPentagons(self):
        val = PAMlib.haspentagons(self.nvert, self.nedge, self.verts, self.edges)
        return val

    def plot(self):
        print '# verts:', self.nvert
        print '# edges:', self.nedge
        v = self.verts
        for e in range(self.edges.shape[0]):
            v0,v1 = self.edges[e,:2]
            v0 -= 1
            v1 -= 1
            pylab.plot([v[v0,0],v[v1,0]],[v[v0,1],v[v1,1]],'k')            
        pylab.plot(self.verts[:,0],self.verts[:,1],'ok')
        pylab.show()


if __name__ == '__main__':
        
    l = Layout()
    l.addMembers('Spars', 1, 5, SP1=[0.1,0.1], EP1=[0.8,0.2], SP2=[0.1,0.9], EP2=[1,1])
    l.addMembers('Ribs', 1, 5, SP1=[0.1,0.1], EP1=[0,1], SP2=[1,0], EP2=[1,1])
    #l.addMembers('Ribs', 1, 5, SP1=[0.1,0], EP1=[0,1], SP2=[1,0], EP2=[1,1])
    #l.addMembers('Spars', 1, 5, SP1=[0,0], EP1=[1,0], SP2=[0,1], EP2=[1,1])
    l.initialize()
    l.plot()
    
