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

def fitAirfoil(wing,airfoil):
    oml0 = wing.oml0
    P1 = []
    P2 = []
    for f in range(len(wing.Ks)):
        ms = []
        ns = []
        for i in range(wing.Ks[f].shape[0]):
            edge = oml0.surf_edge[wing.Ks[f][i,0],0,0]
            edge = abs(edge) - 1
            group = oml0.edge_group[edge] - 1
            ms.append(oml0.group_m[group])
            ns.append(oml0.group_n[group])
        n = sum(ns) - len(ns) + 1

        P = numpy.zeros((airfoil[f].shape[0],4,3))
        for i in range(4):
            P[:,i,:2] = airfoil[f][:,:]
        oml1 = PUBS.PUBS()
        oml1.importSurfaces([P])
        oml1.group_n[1] = n
        oml1.updateEvaluation()
        P1.append(numpy.zeros((n,2)))
        for i in range(n):
            P1[-1][i,:] = oml1.P[oml1.computeIndex(0,i,0,0),:2]

        jQ = numpy.zeros(sum(ms)-2*len(ms)+1)
        counter = 0
        for i in range(wing.Ks[f].shape[0]):
            for u in range(ms[i]-2):
                jQ[counter] = oml0.computeIndex(wing.Ks[f][i,0],u+1,0,2)
                counter += 1
        jQ[counter] = oml0.computeIndex(wing.Ks[f][-f,0],-f,0,2)
        iP = numpy.zeros(sum(ns)-len(ns)-1)
        counter = 0
        for i in range(wing.Ks[f].shape[0]):
            for u in range(ns[i]-1):
                if not (i==0 and u==0):
                    iP[counter] = oml0.computeIndex(wing.Ks[f][i,0],u,0,0)
                    counter += 1

        B = numpy.zeros((iP.shape[0],jQ.shape[0]))
        for j in range(jQ.shape[0]):
            for i in range(iP.shape[0]):
                B[i,j] = oml0.JM[iP[i],jQ[j]]

        R = P1[-1][1:-1]
        dR = B[:,-1]
        B = numpy.delete(B,-1,1)
        R[:,0] -= dR 
        BTB = numpy.dot(B.T,B)
        BTR = numpy.dot(B.T,R)
        sol = numpy.linalg.solve(BTB,BTR)
        sol = numpy.insert(sol, 0, 0, axis=0)
        for i in range(wing.Ks[f].shape[0]):
            sol = numpy.insert(sol, sum(ms[:i+1])-(i+1), 0, axis=0)
        sol[-f,0] = 1
        P2.append(sol)
#        Ps = numpy.dot(B,sol)
#        Ps[:,0] += dR
#        pylab.plot(Ps[:,0],Ps[:,1])
#    pylab.show()

    return P2


