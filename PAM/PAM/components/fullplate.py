from __future__ import division
from PAM.components import component, Property, airfoils
import numpy, pylab
import mpl_toolkits.mplot3d.axes3d as p3


class fullplate(component):

    def __init__(self, nb, nc, half=False, opentip=False):
        self.faces = numpy.zeros((2,2),int)
        self.faces[0,:] = [1,2]
        self.faces[1,:] = [-1,2]

        Ps = []
        Ks = []

        P, K = self.createSurfaces(Ks, nc[::-1], nb, -1, 3, 0)
        for k in range(len(P)):
            for v in range(P[k].shape[1]):
                for u in range(P[k].shape[0]):
                    if (opentip or P[k][u,v,2]!=1) and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                        P[k][u,v,1] = 1
        Ps.extend(P)
        Ks.append(K)

        if not half:
            P, K = self.createSurfaces(Ks, nc, nb, 1, 3, 0)
            for k in range(len(P)):
                for v in range(P[k].shape[1]):
                    for u in range(P[k].shape[0]):
                        if (opentip or P[k][u,v,2]!=1) and P[k][u,v,0]!=0 and P[k][u,v,0]!=1:
                            P[k][u,v,1] = -1
            Ps.extend(P)
            Ks.append(K)

        self.nb = nb
        self.nc = nc
        self.Ps = Ps
        self.Ks = Ks
        self.opentip = opentip

        self.oml0 = []

    def setDOFs(self):
        opentip = self.opentip
        oml0 = self.oml0
        nf = len(self.Ks)
        for f in range(nf):
            for j in range(self.Ks[f].shape[1]):
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,j],:,:] = True
            for i in range(self.Ks[f].shape[0]):
                oml0.surf_c1[self.Ks[f][i,0],:,0] = False
            if opentip:
                for i in range(self.Ks[f].shape[0]):
                    oml0.surf_c1[self.Ks[f][i,-1],:,-1] = False
            for j in range(self.Ks[f].shape[1]):
                oml0.surf_c1[self.Ks[f][-f,j],-f,:] = False
        for f in range(nf):
            for i in range(self.Ks[f].shape[0]):
                edge = oml0.surf_edge[self.Ks[f][i,0],0,0]
                edge = abs(edge) - 1
                oml0.edge_c1[edge,:] = True
                if opentip:
                    edge = oml0.surf_edge[self.Ks[f][i,-1],0,1]
                    edge = abs(edge) - 1
                    oml0.edge_c1[edge,:] = True                    
            edge = oml0.surf_edge[self.Ks[f][-f,0],0,0]
            if edge > 0:
                oml0.edge_c1[abs(edge)-1,-f] = False
            else:
                oml0.edge_c1[abs(edge)-1,f-1] = False  
            if opentip:            
                edge = oml0.surf_edge[self.Ks[f][-f,-1],0,1]
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
            if opentip:
                edge = oml0.surf_edge[self.Ks[f][-f,-1],1,-f]
                if edge > 0:
                    oml0.edge_c1[abs(edge)-1,-1] = False
                else:
                    oml0.edge_c1[abs(edge)-1,0] = False

    def isExteriorDOF(self, f, uType, vType):
        opentip = self.opentip
        value = False
        if f==0:
            if uType==2 and vType==0:
                value = True
            elif uType==0 and (vType==2 or vType==0):
                value = True
            elif (uType==2 or uType==0) and vType==-1 and opentip:
                value = True
        elif f==1:
            if uType==2 and vType==0:
                value = True
            elif uType==-1 and (vType==2 or vType==0):
                value = True
            elif (uType==2 or uType==-1) and vType==-1 and opentip:
                value = True
        return value

    def initializeParameters(self):
        Ns = self.Ns
        self.offset = numpy.zeros(3)
        self.SECTshape = numpy.zeros((len(self.Ks),Ns[0].shape[0],Ns[0].shape[1],3))
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
        self.setAirfoil("naca0012.dat")

    def setAirfoil(self,filename):
        Ps = airfoils.fitAirfoil(self,filename)
        for f in range(len(self.Ks)):
            for j in range(self.Ns[f].shape[1]):
                self.SECTshape[f,:,j,:2] = Ps[f][:,:]
            self.SECTshape[f,-f,:,0] = 1
        
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
                    Qs[f][i,j,:] = numpy.dot(T,self.SECTshape[f,i,j,:]-[a,b,0])*self.props['chord'].data[j]
                    Qs[f][i,j,:] += self.offset + pos

    def computeRotations(self):
        def arctan(x,y):
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

        pos = numpy.zeros((self.Ns[0].shape[1],3))
        for j in range(self.Ns[0].shape[1]):
            pos[j,:] = [self.props['posx'].data[j], self.props['posy'].data[j], self.props['posz'].data[j]]
        for j in range(self.Ns[0].shape[1]):
            if j==0:
                tangent = pos[j+1] - pos[j]
            elif j==self.Ns[0].shape[1]-1:
                tangent = pos[j] - pos[j-1]
            else:
                t1 = pos[j+1] - pos[j]
                t2 = pos[j] - pos[j-1]
                tangent = t1/numpy.linalg.norm(t1) + t2/numpy.linalg.norm(t2)
            x,y,z = tangent
            p = arctan(z,y)
            q = arctan((y**2+z**2)**0.5,x) 
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


if __name__ == '__main__':

    f = fullplate([7,8],[9,10])
    P = f.Ps
    
    ax = p3.Axes3D(pylab.figure())
    for k in range(len(P)):
        ax.plot_wireframe(P[k][:,:,0],P[k][:,:,1],P[k][:,:,2])
    pylab.show()
