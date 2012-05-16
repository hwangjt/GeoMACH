from __future__ import division
import component
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class fulljunction(component.component):

    def __init__(self, toComp, toFace, toDir, toCorner1, toCorner2, fromComp, fromFace, fromDir, fromSide):
        self.toComp = toComp
        self.toFace = toFace
        self.toDir = toDir
        self.toCorner1 = toCorner1
        self.toCorner2 = toCorner2
        self.fromComp = fromComp
        self.fromFace = fromFace
        self.fromDir = fromDir
        self.fromSide = fromSide

        u1 = toCorner1[0]
        u2 = toCorner2[0]
        v1 = toCorner1[1]
        v2 = toCorner2[1]

        Ps = []
        Ks = []
        
        for i in range(fromComp.Ks[0].shape[0]):
            if toDir==0:
                n = toComp.Ps[toComp.Ks[toFace][u1,v1]].shape[0]
                edge1 = toComp.Ps[toComp.Ks[toFace][u1,v1+1+i]][0,:,:]
            else:
                n = toComp.Ps[toComp.Ks[toFace][u2,v1]].shape[1]
                edge1 = toComp.Ps[toComp.Ks[toFace][u2-1-i,v1]][::-1,0,:]
            if fromDir==0:
                edge2 = fromComp.Ps[fromComp.Ks[fromFace][-1-i,-fromSide]][::-1,-fromSide,:]
            else:
                edge2 = fromComp.Ps[fromComp.Ks[fromFace][i,-fromSide]][:,-fromSide,:]
            Ps.append(self.createInterface(n, edge1, edge2))

            if toDir==0:
                n = toComp.Ps[toComp.Ks[toFace][u2,v1]].shape[0]
                edge2 = toComp.Ps[toComp.Ks[toFace][u2,v1+1+i]][-1,:,:]
            else:
                n = toComp.Ps[toComp.Ks[toFace][u2,v2]].shape[1]
                edge2 = toComp.Ps[toComp.Ks[toFace][u2-1-i,v2]][::-1,-1,:]
            if fromDir==0:
                edge1 = fromComp.Ps[fromComp.Ks[1-fromFace][i,-fromSide]][:,-fromSide,:]
            else:
                edge1 = fromComp.Ps[fromComp.Ks[1-fromFace][-1-i,-fromSide]][::-1,-fromSide,:]
            Ps.append(self.createInterface(n, edge1, edge2))

        if toDir==0:
            n = toComp.Ps[toComp.Ks[toFace][u1,v1]].shape[1]
            edge1 = toComp.Ps[toComp.Ks[toFace][u1,v1]][:,0,:]
            edge2 = Ps[0][:,0,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))
        else:
            n = toComp.Ps[toComp.Ks[toFace][u2,v1]].shape[0]
            edge1 = toComp.Ps[toComp.Ks[toFace][u2,v1]][-1,:,:]
            edge2 = Ps[0][:,0,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))

        if toDir==0:
            n = toComp.Ps[toComp.Ks[toFace][u2,v1]].shape[1]
            edge1 = toComp.Ps[toComp.Ks[toFace][u2,v1]][:,0,:]
            edge2 = Ps[1][:,0,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))
        else:
            n = toComp.Ps[toComp.Ks[toFace][u2,v2]].shape[0]
            edge1 = toComp.Ps[toComp.Ks[toFace][u2,v2]][-1,:,:]
            edge2 = Ps[1][:,0,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))

        if toDir==0:
            n = toComp.Ps[toComp.Ks[toFace][u1,v2]].shape[1]
            edge2 = toComp.Ps[toComp.Ks[toFace][u1,v2]][:,-1,:]
            edge1 = Ps[2*fromComp.Ks[0].shape[0]-2][:,-1,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))
        else:
            n = toComp.Ps[toComp.Ks[toFace][u1,v1]].shape[0]
            edge2 = toComp.Ps[toComp.Ks[toFace][u1,v1]][0,:,:]
            edge1 = Ps[2*fromComp.Ks[0].shape[0]-2][:,-1,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))

        if toDir==0:
            n = toComp.Ps[toComp.Ks[toFace][u2,v2]].shape[1]
            edge2 = toComp.Ps[toComp.Ks[toFace][u2,v2]][:,-1,:]
            edge1 = Ps[2*fromComp.Ks[0].shape[0]-1][:,-1,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))
        else:
            n = toComp.Ps[toComp.Ks[toFace][u1,v2]].shape[0]
            edge2 = toComp.Ps[toComp.Ks[toFace][u1,v2]][0,:,:]
            edge1 = Ps[2*fromComp.Ks[0].shape[0]-1][:,-1,:]
            Ps.append(self.createInterface(n, edge1, edge2, True))

        K = numpy.zeros((2,2+fromComp.Ks[0].shape[0]),int)
        counter = 0
        for j in range(fromComp.Ks[0].shape[0]):
            for i in range(2):
                K[i,j+1] = counter
                counter += 1
        for j in range(2):
            for i in range(2):
                K[i,-j] = counter
                counter += 1
        Ks.append(K)

        for j in range(v2,v1-1,-1):
            for i in range(u2,u1-1,-1):
                toComp.Ps.pop(toComp.Ks[toFace][i,j])
                for f in range(len(toComp.Ks)):
                    for v in range(toComp.Ks[f].shape[1]):
                        for u in range(toComp.Ks[f].shape[0]):
                            if toComp.Ks[f][u,v] > toComp.Ks[toFace][i,j]:
                                toComp.Ks[f][u,v] -= 1
                toComp.Ks[toFace][i,j] = -1
            
        self.Ps = Ps
        self.Ks = Ks

        self.oml0 = []

    def setDOFs(self):
        oml0 = self.oml0
        Ks = self.Ks
        for f in [0]:
            for j in range(Ks[f].shape[1]):
                for i in range(2):
                    oml0.surf_c1[Ks[f][i,j],:,:] = True
            for j in range(1,Ks[f].shape[1]-1):
                for i in range(2):
                    oml0.surf_c1[Ks[f][i,j],i-1,:] = False
            oml0.surf_c1[Ks[f][0,0],-1,-1] = False
            oml0.surf_c1[Ks[f][1,0],0,-1] = False
            oml0.surf_c1[Ks[f][0,-1],-1,0] = False
            oml0.surf_c1[Ks[f][1,-1],0,0] = False

    def isExteriorDOF(self, f, uType, vType):
        return False

    def initializeParameters(self):
        Ns = self.Ns

    def propagateQs(self):
        toComp = self.toComp
        toFace = self.toFace
        toDir = self.toDir
        toCorner1 = self.toCorner1
        toCorner2 = self.toCorner2
        fromComp = self.fromComp
        fromFace = self.fromFace
        fromDir = self.fromDir
        fromSide = self.fromSide

        u1 = toCorner1[0]
        u2 = toCorner2[0]
        v1 = toCorner1[1]
        v2 = toCorner2[1]

        Ns = self.Ns
        Qs = self.Qs

        Qs[0][:,:,:] = 0

        for f in range(1):
            si = self.getni(f,0)
            sj = self.getni(f,1)
            ti = toComp.getni(toFace,0)
            tj = toComp.getni(toFace,1)
            for j in range(int(Ns[f].shape[1]-sj[0]-sj[-1])):
                if toDir==0:
                    Q1 = toComp.Qs[toFace][sum(ti[:u1])-1,sum(tj[:v1+1])+j]
                else:
                    Q1 = toComp.Qs[toFace][sum(ti[:u2])-j,sum(tj[:v1])-1]
                if fromDir==0:
                    Q2 = fromComp.Qs[fromFace][-1-j,-fromSide]
                else:
                    Q2 = fromComp.Qs[fromFace][j,-fromSide]
                Qs[f][:si[0]+1,j+sj[0],:] += numpy.outer(numpy.linspace(1,0,si[0]+2)[1:],Q1)
                Qs[f][:si[0]+1,j+sj[0],:] += numpy.outer(numpy.linspace(0,1,si[0]+2)[1:],Q2)

                if toDir==0:
                    Q2 = toComp.Qs[toFace][sum(ti[:u2+1])+1,sum(tj[:v1+1])+j]
                else:
                    Q2 = toComp.Qs[toFace][sum(ti[:u2])-j,sum(tj[:v1+1])+1]
                if fromDir==0:
                    Q1 = fromComp.Qs[1-fromFace][j,-fromSide]
                else:
                    Q1 = fromComp.Qs[1-fromFace][-1-j,-fromSide]
                Qs[f][si[0]:,j+sj[0],:] += numpy.outer(numpy.linspace(1,0,si[1]+2)[:-1],Q1)
                Qs[f][si[0]:,j+sj[0],:] += numpy.outer(numpy.linspace(0,1,si[1]+2)[:-1],Q2)                 

            for i in range(si[0]):
                if toDir==0:
                    Q1 = toComp.Qs[toFace][sum(ti[:u1])+i,sum(tj[:v1])-1]
                else:
                    Q1 = toComp.Qs[toFace][sum(ti[:u2+1])+1,sum(tj[:v1])+i]
                Q2 = Qs[f][i,sj[0]+1]
                Qs[f][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(1,0,sj[0]+3)[1:-1],Q1)
                Qs[f][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(0,1,sj[0]+3)[1:-1],Q2)

                if toDir==0:
                    Q2 = toComp.Qs[toFace][sum(ti[:u1])+i,sum(tj[:v2+1])+1]
                else:
                    Q2 = toComp.Qs[toFace][sum(ti[:u1])-1,sum(tj[:v1])+i]
                Q1 = Qs[f][i,sum(sj[:-1])-1]
                Qs[f][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(1,0,sj[-1]+3)[1:-1],Q1)
                Qs[f][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(0,1,sj[-1]+3)[1:-1],Q2)

            for i in range(si[0],sum(si)):
                if toDir==0:
                    Q1 = toComp.Qs[toFace][sum(ti[:u2])+i-si[0],sum(tj[:v1])-1]
                else:
                    Q1 = toComp.Qs[toFace][sum(ti[:u2+1])+1,sum(tj[:v2])+i-si[0]]
                Q2 = Qs[f][i,sj[0]+1]
                Qs[f][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(1,0,sj[0]+3)[1:-1],Q1)
                Qs[f][i,:sj[0]+1,:] += numpy.outer(numpy.linspace(0,1,sj[0]+3)[1:-1],Q2)

                if toDir==0:
                    Q2 = toComp.Qs[toFace][sum(ti[:u2])+i-si[0],sum(tj[:v2+1])+1]
                else:
                    Q2 = toComp.Qs[toFace][sum(ti[:u1])-1,sum(tj[:v2])+i-si[0]]
                Q1 = Qs[f][i,sum(sj[:-1])-1]
                Qs[f][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(1,0,sj[-1]+3)[1:-1],Q1)
                Qs[f][i,sum(sj[:-1]):,:] += numpy.outer(numpy.linspace(0,1,sj[-1]+3)[1:-1],Q2)
