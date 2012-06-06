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
        self.offset = numpy.zeros(3)
        self.SECTshape = numpy.zeros((Ns[0].shape[0],Ns[0].shape[1],len(self.Ks),3))
        self.SECTrot0 = numpy.zeros((Ns[0].shape[1],3))
        self.props = {
            'chord':Property(Ns[0].shape[1]),
            'posx':Property(Ns[0].shape[1]),
            'posy':Property(Ns[0].shape[1]),
            'posz':Property(Ns[0].shape[1]),
            'rotx':Property(Ns[0].shape[1]),
            'roty':Property(Ns[0].shape[1]),
            'rotz':Property(Ns[0].shape[1]),
            'prpx':Property(Ns[0].shape[1]),
            'prpy':Property(Ns[0].shape[1])
}
        self.SECTchord = numpy.zeros((Ns[0].shape[1],1))
        self.SECTpos = numpy.zeros((Ns[0].shape[1],3))
        self.SECTrot = numpy.zeros((Ns[0].shape[1],3))
        self.SECTbend = numpy.ones((Ns[0].shape[1],3))
        self.setAirfoil("naca0012.dat")

    def setBend(self,a=0,b=1,p=1):
        for k in range(3):
            self.SECTbend[:-1,k] = numpy.linspace(a,b,self.SECTpos.shape[0]-1)**p

    def setSpan(self, span):
        self.SECTpos[:-1,2] = numpy.linspace(0,span,self.SECTpos.shape[0]-1)

    def setTaper(self, root, tip):
        self.SECTchord[:-1,0] = numpy.linspace(root,tip,self.SECTchord.shape[0]-1)

    def setSweep(self, sweep):
        self.SECTpos[:-1,0] = numpy.linspace(0,sweep,self.SECTpos.shape[0]-1)

    def setTaper2(self, root, tip):
        self.SECTchord[:-1,0] = root + (tip-root)*numpy.linspace(0,1,self.SECTchord.shape[0]-1)

    def setSweep2(self, sweep):
        w = 0.4
        self.SECTpos[:-1,0] = w*sweep*numpy.linspace(0,1,self.SECTpos.shape[0]-1)**2 + (1-w)*sweep*numpy.linspace(0,1,self.SECTpos.shape[0]-1)

    def setAirfoil(self,filename):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.SECTshape[1:-1,j,f,:2] = Ps[f][:,:]
            self.SECTshape[-f,:,f,0] = 1
        
    def propagateQs(self):
        a = 0.25
        b = 0.0
        Ns = self.Ns
        Qs = self.Qs
        self.computeRotations()
        for f in range(len(self.Ks)):
            Qs[f][:,:,:] = 0
            for j in range(Ns[f].shape[1]):
                pos = [self.props['posx'].data[j], self.props['posy'].data[j], self.props['posz'].data[j]]
                rot = [self.props['rotx'].data[j], self.props['roty'].data[j], self.props['rotz'].data[j]]
                prp = [self.props['prpx'].data[j], self.props['prpy'].data[j], 0]
                pos = numpy.array(pos)
                rot = numpy.array(rot)
                prp = numpy.array(prp)
                T = self.computeRtnMtx(rot+self.SECTrot0[j,:]*prp)
                for i in range(Ns[f].shape[0]):
                    Qs[f][i,j,:] = numpy.dot(T,self.SECTshape[i,j,f,:]-[a,b,0])*self.props['chord'].data[j]
                    Qs[f][i,j,:] += self.offset + pos

    def computeRotations(self):
        for j in range(self.Ns[0].shape[1]-1):
            pos0 = [self.props['posx'].data[j-1], self.props['posy'].data[j-1], self.props['posz'].data[j-1]]
            pos1 = [self.props['posx'].data[j], self.props['posy'].data[j], self.props['posz'].data[j]]
            pos2 = [self.props['posx'].data[j+1], self.props['posy'].data[j+1], self.props['posz'].data[j+1]]
            pos0 = numpy.array(pos0)
            pos1 = numpy.array(pos1)
            pos2 = numpy.array(pos2)
            if j==0:
                tangent = pos2 - pos1
            elif j==self.Ns[0].shape[1]-2:
                tangent = pos1 - pos0
            else:
                t1 = pos2 - pos1
                t2 = pos1 - pos0
                tangent = t1/numpy.linalg.norm(t1) + t2/numpy.linalg.norm(t2)
            x,y,z = tangent
            if y==0 and z==0:
                q = numpy.pi/2.0
            else:
                q = numpy.arctan(x/(y**2+z**2)**0.5)
            if z==0:
                p = numpy.pi/2.0
            else:
                p = numpy.arctan(y/z)
            self.SECTrot0[j,:2] = [p,q]
            self.SECTrot0[j,:2] *= 180.0/numpy.pi

    def computeCproj(self):
        ax = p3.Axes3D(pylab.figure())
        Ns = self.Ns
        iLE = numpy.zeros((Ns[0].shape[1],3))
        iTE = numpy.zeros((Ns[0].shape[1],3))
        for j in range(Ns[0].shape[1]):
            v = Ns[0][0,j,4]
            jj = Ns[0][0,j,3]
            surf = self.Ks[0][-1,jj]
            col = self.oml0.J.getcol(self.oml0.computeIndex(surf,-1,v,1)).tocsc()
            index = col.indices[numpy.argmax(col.data)]
            iLE[j] = self.oml0.computeProjection(self.oml0.P[index],surf)[1:]
            surf = self.Ks[0][0,jj]
            col = self.oml0.J.getcol(self.oml0.computeIndex(surf,0,v,1)).tocsc()
            index = col.indices[numpy.argmax(col.data)]
            iTE[j] = self.oml0.computeProjection(self.oml0.P[index],surf)[1:]

    def computeRtnMtx(self, rot):
        p,q,r = rot*numpy.pi/180.0
        cos = numpy.cos
        sin = numpy.sin
        T0 = numpy.eye(3)
        T = numpy.zeros((3,3))
        T[0,:] = [   1   ,   0   ,   0   ]
        T[1,:] = [   0   , cos(p), sin(p)]
        T[2,:] = [   0   ,-sin(p), cos(p)]
        T0 = numpy.dot(T,T0)
        T[0,:] = [ cos(q),   0   , sin(q)]
        T[1,:] = [   0   ,   1   ,   0   ]
        T[2,:] = [-sin(q),   0   , cos(q)]
        T0 = numpy.dot(T,T0)
        T[0,:] = [ cos(r), sin(r),   0   ]
        T[1,:] = [-sin(r), cos(r),   0   ]
        T[2,:] = [   0   ,   0   ,   1   ]
        T0 = numpy.dot(T,T0)
        return T0

class Property(object):

    def __init__(self, n0):
        self.data = numpy.zeros(n0)
        self.set([0.0,1.0],[0,1])

    def set(self, val, ind, p=None, w=None, d=None):
        self.G = numpy.zeros((numpy.array(val).shape[0],5))
        self.G[:,0] = val
        self.G[:,1] = ind
        if p==None:
            self.G[:,2] = 2
        else:
            self.G[:,2] = p
        if w==None:
            self.G[:,3] = 0
        else:
            self.G[:,3] = w
        if d==None:
            self.G[:,4] = 0
        else:
            self.G[:,4] = d
        self.evaluate()

    def evaluate(self):
        indices = numpy.round((self.data.shape[0]-1)*self.G[:,1])
        for i in range(self.G.shape[0]-1):
            v1 = self.G[i,0]
            v2 = self.G[i+1,0]
            i1 = int(indices[i])
            i2 = int(indices[i+1])
            p = self.G[i,2]
            w = self.G[i,3]
            x = numpy.linspace(0,1,i2-i1+1)
            if self.G[i,4]==0:
                self.data[i1:i2+1] = v1 + (v2-v1)*(1-w)*x + (v2-v1)*w*x**p
            else:
                self.data[i1:i2+1] = v2 + (v1-v2)*(1-w)*x[::-1] + (v1-v2)*w*x[::-1]**p


if __name__ == '__main__':

    f = fullplate([7,8],[9,10])
    P = f.Ps
    
    ax = p3.Axes3D(pylab.figure())
    for k in range(len(P)):
        ax.plot_wireframe(P[k][:,:,0],P[k][:,:,1],P[k][:,:,2])
    pylab.show()
