from __future__ import division
import component
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class fullplate(component.component):

    def __init__(self, nb, nc):
        Ps = []
        Ks = []

        P, K = self.createSurfaces(Ks, nc, nb, -1, 3, 0)
        for k in range(len(P)):
            for v in range(P[k].shape[1]):
                for u in range(P[k].shape[0]):
                    if P[k][u,v,2]!=1 and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                        P[k][u,v,1] = 1
        Ps.extend(P)
        Ks.append(K)

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
        for f in range(2):
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,:] = True
            for i in range(self.Ks[f].shape[0]):
                oml0.surf_c1[self.Ks[f][i,0],:,0] = False
            for j in range(self.Ks[f].shape[1]):
                oml0.surf_c1[self.Ks[f][-f,j],-f,:] = False
        for f in range(2):
            for i in range(self.Ks[f].shape[0]):
                edge = oml0.surf_edge[self.Ks[f][i,0],0,0]
                edge = abs(edge) - 1
                oml0.edge_c1[edge,:] = True
            edge = oml0.surf_edge[self.Ks[f][-f,0],0,0]
            if edge > 0:
                oml0.edge_c1[abs(edge)-1,-f] = False
            else:
                oml0.edge_c1[abs(edge)-1,f-1] = False                
        for f in range(2):
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
        self.offset = numpy.zeros(3)
        self.span = 4
        self.chord = numpy.ones(Ns[0].shape[1])
        self.sweep = numpy.zeros(Ns[0].shape[1])
        self.dihedral = numpy.zeros(Ns[0].shape[1])
        self.twist = numpy.zeros(Ns[0].shape[1])
        self.thickness = 1

    def propagateQs(self):
        Ns = self.Ns
        Qs = self.Qs
        for f in range(2):
            for j in range(Ns[f].shape[1]):
                for i in range(Ns[f].shape[0]):
                    Qs[f][i,j,:] = self.offset
                    Qs[f][i,j,0] += self.sweep[j]
                    Qs[f][i,j,1] += self.dihedral[j]
                    Qs[f][i,j,2] += j/(Ns[f].shape[1]-2)*self.span
                    if f==0:
                        Qs[f][i,j,0] += (1-i/(Ns[f].shape[0]-2))*self.chord[j]
                        if i>0:
                            Qs[f][i,j,1] += self.thickness/2.0
                    else:
                        Qs[f][i,j,0] += (i-1)/(Ns[f].shape[0]-2)*self.chord[j]
                        if i<Ns[f].shape[0]-1:
                            Qs[f][i,j,1] -= self.thickness/2.0
                            
                                     


if __name__ == '__main__':

    f = fullplate([2,3],[3,4])
    P = f.Ps
    print f.Ks
    
    ax = p3.Axes3D(pylab.figure())
    for k in range(len(P)):
        ax.plot_wireframe(P[k][:,:,0],P[k][:,:,1],P[k][:,:,2])
    pylab.show()
