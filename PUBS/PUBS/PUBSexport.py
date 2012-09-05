from __future__ import division
import numpy, pylab, copy
import PUBS, PUBSlib

def export(Ps, filename):
    oml0 = PUBS.PUBS()
    oml0.importSurfaces(Ps)

    file = open(filename+'.dat',"w")

    file.write(str(oml0.nvert)+'\n')
    file.write(str(oml0.nedge)+'\n')
    file.write(str(oml0.nsurf)+'\n')

    for i in range(oml0.nvert):
        file.write(str(oml0.C[i,0])+' ')
        file.write(str(oml0.C[i,1])+' ')
        file.write(str(oml0.C[i,2])+'\n')

    edgePtr = numpy.zeros((oml0.nedge,3),int)
    for i in range(oml0.nsurf):
        for u in range(2):
            for v in range(2):
                edgePtr[abs(oml0.surf_edge[i,u,v])-1,:] = [i,u,v]
    ms = PUBSlib.getsurfacesizes(oml0.nsurf, oml0.nedge, oml0.ngroup, oml0.surf_edge, oml0.edge_group, oml0.group_m)
    for i in range(oml0.nedge):
        surf = edgePtr[i,0]
        edge0 = edgePtr[i,1]
        edge1 = edgePtr[i,2]
        if edge0==0:
            if edge1==0:
                start = [0,0]
                end = [1,0]
            else:
                start = [0,1]
                end = [1,1]
        else:
            if edge1==0:
                start = [0,0]
                end = [0,1]
            else:
                start = [1,0]
                end = [1,1]
        if oml0.surf_edge[surf,edge0,edge1] < 0:
            temp = start
            start = end
            end = temp
        file.write(str(oml0.surf_vert[surf,start[0],start[1]])+' ')
        file.write(str(oml0.surf_vert[surf,end[0],end[1]])+'\n')
        group = oml0.edge_group[i] - 1
        k = oml0.group_k[group]
        m = oml0.group_m[group]
        file.write(str(k)+' ')
        file.write(str(m)+'\n')
        d = oml0.group_d[oml0.knot_index[group,0]:oml0.knot_index[group,1]]
        for j in range(d.shape[0]):
            file.write(str(d[j])+' ')
        file.write('\n')
        C = PUBSlib.getsurfacep(surf+1, oml0.nC, ms[surf,0], ms[surf,1], oml0.nsurf, oml0.nedge, oml0.nvert, oml0.surf_vert, oml0.surf_edge, oml0.surf_index_C, oml0.edge_index_C, oml0.C)
        if edge0==0:
            if edge1==0:
                edgeCs = copy.copy(C[:,0,:])
            else:
                edgeCs = copy.copy(C[:,-1,:])
        else:
            if edge1==0:
                edgeCs = copy.copy(C[0,:,:])
            else:
                edgeCs = copy.copy(C[-1,:,:])
        if oml0.surf_edge[surf,edge0,edge1] < 0:
            edgeCs[:,:] = copy.copy(edgeCs[::-1,:])
        for j in range(edgeCs.shape[0]):
            file.write(str(edgeCs[j,0])+' ')
            file.write(str(edgeCs[j,1])+' ')
            file.write(str(edgeCs[j,2])+'\n')

    for i in range(oml0.nsurf):
        for u in range(2):
            for v in range(2):
                file.write(str(oml0.surf_edge[i,u,v])+' ')
        file.write('\n')
        ugroup = oml0.edge_group[abs(oml0.surf_edge[i,0,0])-1] - 1
        vgroup = oml0.edge_group[abs(oml0.surf_edge[i,1,0])-1] - 1
        ku = oml0.group_k[ugroup]
        mu = oml0.group_m[ugroup]
        kv = oml0.group_k[vgroup]
        mv = oml0.group_m[vgroup]
        file.write(str(ku)+' ')
        file.write(str(mu)+' ')
        file.write(str(kv)+' ')
        file.write(str(mv)+' ')
        file.write('\n')    
        d = oml0.group_d[oml0.knot_index[ugroup,0]:oml0.knot_index[ugroup,1]]
        for j in range(d.shape[0]):
            file.write(str(d[j])+' ')
        file.write('\n') 
        d = oml0.group_d[oml0.knot_index[vgroup,0]:oml0.knot_index[vgroup,1]]
        for j in range(d.shape[0]):
            file.write(str(d[j])+' ')
        file.write('\n')
        C = PUBSlib.getsurfacep(i+1, oml0.nC, ms[i,0], ms[i,1], oml0.nsurf, oml0.nedge, oml0.nvert, oml0.surf_vert, oml0.surf_edge, oml0.surf_index_C, oml0.edge_index_C, oml0.C)
        for v in range(C.shape[1]):
            for u in range(C.shape[0]):
                file.write(str(C[u,v,0])+' ')
                file.write(str(C[u,v,1])+' ')
                file.write(str(C[u,v,2])+'\n')

    file.close()
