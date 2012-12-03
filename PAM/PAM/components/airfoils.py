from __future__ import division
import numpy, pylab, copy
import PUBS



def getAirfoil(filename):
    if filename[:4]=='naca' or filename[:4]=='NACA':
        n = 50
        m = int(filename[4])/100.0
        p = int(filename[5])/10.0
        t = int(filename[6:8])/100.0
        x = numpy.linspace(0,1,n)**2
        ys = t/0.2*(0.2969*x**0.5 - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1036*x**4)
        yc = numpy.zeros(n)
        if p != 0:
            y1 = m*x/p**2*(2*p-x)
            y2 = m*(1-x)/(1-p)**2*(1+x-2*p)
            for i in range(n):
                if x[i] < p:
                    yc[i] = y1[i]
                else:
                    yc[i] = y2[i]
        upper = numpy.zeros((n,2))
        lower = numpy.zeros((n,2))
        upper[:,0] = x[::-1]
        lower[:,0] = x
        upper[:,1] = ys[::-1] + yc[::-1]
        lower[:,1] = -ys + yc
    else:
        data = numpy.genfromtxt('../components/airfoils/'+filename)    
        if data[0,0] > data[1,0]:
            mark = numpy.argmin(data,0)[0]
            upper = data[:mark+1,:]
            lower = data[mark:,:]
        else:
            for i in range(data.shape[0]-1):
                if abs(data[i+1,0]-data[i,0]) > 0.8:
                    mark = i
            upper = data[mark::-1,:]
            lower = data[mark+1:,:]
    return [upper, lower]

def getP(nP, airfoil):
    P = numpy.zeros((airfoil.shape[0],4,3))
    for i in range(4):
        P[:,i,:2] = airfoil[:,:]
    oml0 = PUBS.PUBS([P])
    oml0.group_n[1] = nP
    oml0.updateEvaluation()
    P = numpy.zeros((nP,2))
    for i in range(nP):
        P[i,:] = oml0.P[oml0.getIndex(0,i,0,0),:2]
    return P

def getQ(ms, ns, P0):
    Ps = []
    for i in range(ns.shape[0]):
        P = numpy.zeros((4,ns[i]+1,3))
        for j in range(4):
            P[j,:,:2] = P0[sum(ns[:i]):sum(ns[:i+1])+1,:]
        Ps.append(P)
    oml0 = PUBS.PUBS(Ps)
    for i in range(ms.shape[0]):
        oml0.group_m[1+i] = ms[i]+1
    oml0.updateBsplines()
    Q = numpy.zeros((sum(ms) + 1,2))
    for i in range(ns.shape[0]):
        for j in range(ms[i]+1):
            Q[sum(ms[:i])+j] = oml0.Q[oml0.getIndex(i,0,j,2),:2]
    return Q

def fitAirfoil(wing,filename):
    airfoil = getAirfoil(filename)
    oml0 = wing.oml0
    Qs = []
    for f in range(2):
        nsurf = wing.Ks[f].shape[0]
        ms = numpy.zeros(nsurf,int)
        ns = numpy.zeros(nsurf,int)
        for i in range(nsurf):
            surf = wing.Ks[f][i,0]
            edge = abs(oml0.surf_edge[surf,0,0]) - 1
            group = oml0.edge_group[edge] - 1
            ms[i] = oml0.group_m[group] - 1
            ns[i] = oml0.group_n[group] - 1
        nP = sum(ns) + 1

        P = getP(nP, airfoil[f])
        Q = getQ(ms, ns, P)
        Qs.append(Q)

    return Qs


