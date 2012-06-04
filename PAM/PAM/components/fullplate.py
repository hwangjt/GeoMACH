from __future__ import division
from PAM.components import component, airfoils
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class fullplate(component):

    def __init__(self, nb, nc, half=False):
        Ps = []
        Ks = []

        P, K = self.createSurfaces(Ks, nc[::-1], nb, -1, 3, 0)
        for k in range(len(P)):
            for v in range(P[k].shape[1]):
                for u in range(P[k].shape[0]):
                    if P[k][u,v,2]!=1 and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                        P[k][u,v,1] = 1
        Ps.extend(P)
        Ks.append(K)

        if not half:
            P, K = self.createSurfaces(Ks, nc, nb, 1, 3, 0)
            for k in range(len(P)):
                for v in range(P[k].shape[1]):
                    for u in range(P[k].shape[0]):
                        if P[k][u,v,2]!=1 and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                            P[k][u,v,1] = -1
            Ps.extend(P)
            Ks.append(K)

        self.nb = nb
        self.nc = nc
        self.Ps = Ps
        self.Ks = Ks

        self.oml0 = []

    def setDOFs(self):
        oml0 = self.oml0
        nf = len(self.Ks)
        for f in range(nf):
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,:] = True
            for i in range(self.Ks[f].shape[0]):
                oml0.surf_c1[self.Ks[f][i,0],:,0] = False
            for j in range(self.Ks[f].shape[1]):
                oml0.surf_c1[self.Ks[f][-f,j],-f,:] = False
        for f in range(nf):
            for i in range(self.Ks[f].shape[0]):
                edge = oml0.surf_edge[self.Ks[f][i,0],0,0]
                edge = abs(edge) - 1
                oml0.edge_c1[edge,:] = True
            edge = oml0.surf_edge[self.Ks[f][-f,0],0,0]
            if edge > 0:
                oml0.edge_c1[abs(edge)-1,-f] = False
            else:
                oml0.edge_c1[abs(edge)-1,f-1] = False                
        for f in range(nf):
            for j in range(self.Ks[f].shape[1]):
                edge = oml0.surf_edge[self.Ks[f][-f,j],1,-f]
                edge = abs(edge) - 1
                oml0.edge_c1[edge,:] = True
            edge = oml0.surf_edge[self.Ks[f][-f,0],1,-f]
            if edge > 0:
                oml0.edge_c1[abs(edge)-1,0] = False
            else:
                oml0.edge_c1[abs(edge)-1,-1] = False

    def isExteriorDOF(self, f, uType, vType):
        value = False
        if f==0:
            if uType==2 and vType==0:
                value = True
            elif uType==0 and (vType==2 or vType==0):
                value = True
        elif f==1:
            if uType==2 and vType==0:
                value = True
            elif uType==-1 and (vType==2 or vType==0):
                value = True
        return value

    def initializeParameters(self):
        Ns = self.Ns
        self.T = numpy.zeros(3)
        self.A = numpy.zeros((Ns[0].shape[0],Ns[0].shape[1],len(self.Ks),2))
        self.S = numpy.zeros((Ns[0].shape[1],3))
        self.R = numpy.zeros((Ns[0].shape[1],3))
        self.B = numpy.zeros((Ns[0].shape[1],2))
        self.C = numpy.zeros(Ns[0].shape[1])
        self.setAirfoil("naca0012.dat")

    def setSpan(self, span):
        self.S[:-1,2] = numpy.linspace(0,span,self.S.shape[0]-1)

    def setTaper(self, root, tip):
        self.C[:-1] = numpy.linspace(root,tip,self.C.shape[0]-1)

    def setSweep(self, sweep):
        self.S[:-1,0] = numpy.linspace(0,sweep,self.S.shape[0]-1)
        
    def propagateQs(self):
        Ns = self.Ns
        Qs = self.Qs
        for f in range(len(self.Ks)):
            Qs[f][:,:,:] = 0
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    Qs[f][i,j,:2] = self.A[i,j,f,:]*self.C[j]
                    Qs[f][i,j,:] += self.T + self.S[j]

    def setAirfoil(self,filename):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.A[1:-1,j,f,:] = Ps[f][:,:]
            self.A[-f,:,f,0] = 1


if __name__ == '__main__':

    f = fullplate([7,8],[9,10])
    P = f.Ps
    
    ax = p3.Axes3D(pylab.figure())
    for k in range(len(P)):
        ax.plot_wireframe(P[k][:,:,0],P[k][:,:,1],P[k][:,:,2])
    pylab.show()
