"""
B-spline surface-modeling engine (BSE) vec class
"""

# pylint: disable=E1101

from __future__ import division
import numpy



class BSEvec(object):

    def __init__(self, name, size, ndim, hidden):
        self.array = numpy.zeros((size, ndim))
        self.name = name
        self.size = size
        self.ndim = ndim
        self._hidden = hidden
        self._file = None
        self._default_var_names = ['v' + str(idim) 
                                   for idim in xrange(self.ndim)]

    def _open_file(self, filename):
        self._file = open(filename, 'w')

    def _write_tec_header(self, title, variables):
        self._file.write('title = ' + title + '\n')
        self._file.write('variables = ')
        for ivar in range(len(variables)):
            self._file.write(variables[ivar] + ',')
        self._file.write('\n')

    def _write(self, text):
        self._file.write(text)

    def _write_line(self, data, label=''):
        self._file.write(label)
        for ind in range(data.shape[0]):
            if data[ind] == data[ind]:
                self._file.write(str(data[ind]) + ' ')
            else:
                self._file.write(str(0.0) + ' ')
        self._file.write('\n')

    def _close_file(self):
        self._file.close()

    def export_tec_scatter(self, filename=None, var_names=None):
        if filename is None:
            filename = self.name + '_scatter.dat'
        if var_names is None:
            var_names = self._default_var_names

        self._open_file(filename)
        self._write_tec_header('BSE output', var_names)
        for ind in xrange(self.size):
            self._write_line(self.array[ind, :])
        self._close_file()



class BSEvecUns(BSEvec):

    pass



class BSEvecStr(BSEvec):

    def __init__(self, name, size, ndim, surf_sizes, hidden,
                 bse=None):
        super(BSEvecStr, self).__init__(name, size, ndim, hidden)

        self._bse = bse
        self.surfs = []
        if surf_sizes is not None:
            ind1, ind2 = 0, 0
            for isurf in xrange(surf_sizes.shape[0]):
                num_u, num_v = surf_sizes[isurf, :]
                ind2 += num_u * num_v
                surf = self.array[ind1:ind2]
                surf = surf.reshape((num_u, num_v, ndim), 
                                    order='F')
                ind1 += num_u * num_v
                self.surfs.append(surf)

    def __call__(self, isurf):
        return self.surfs[isurf]

    def export_tec_str(self, filename=None, var_names=None):
        if filename is None:
            filename = self.name + '_surf.dat'
        if var_names is None:
            var_names = self._default_var_names

        self._open_file(filename)
        self._write_tec_header('BSE output', var_names)
        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                surf = self.surfs[isurf]
                num_u, num_v = surf.shape[:2]
                self._write('zone i='+str(num_u) + \
                            ', j=' + str(num_v) + \
                            ', DATAPACKING=POINT\n')
                for ind_v in xrange(num_v):
                    for ind_u in xrange(num_u):
                        self._write_line(surf[ind_u, ind_v, :])
        self._close_file()

    def export_STL(self, filename=None, var_names=None):
        if filename is None:
            filename = self.name + '.stl'
        if var_names is None:
            var_names = self._default_var_names

        self._open_file(filename)
        self._write('solid model\n')
        for surf in self.surfs:
            num_u, num_v = surf.shape[:2]
            for ind_v in xrange(num_v - 1):
                for ind_u in xrange(num_u - 1):
                    pt1 = surf[ind_u, ind_v+1, :]
                    pt2 = surf[ind_u, ind_v, :]
                    pt3 = surf[ind_u+1, ind_v, :]
                    
                    for pts in [[surf[ind_u, ind_v+1, :],
                                 surf[ind_u, ind_v, :],
                                 surf[ind_u+1, ind_v, :]],
                                [surf[ind_u+1, ind_v, :],
                                 surf[ind_u+1, ind_v+1, :],
                                 surf[ind_u, ind_v+1, :]]]:
                        nor = numpy.cross(pts[2] - pts[1], 
                                          pts[0] - pts[1])
                        nor /= numpy.linalg.norm(nor)
                        self._write_line(nor, 'facet normal')
                        self._write('outer loop\n')
                        self._write_line(pts[0], 'vertex ')
                        self._write_line(pts[1], 'vertex ')
                        self._write_line(pts[2], 'vertex ')
                        self._write('endloop\n')
                        self._write('endfacet\n')
        self._write('endsolid model')
        self._close_file()

    def export_IGES(self, filename=None):
        ks = []
        ms = []
        ds = [[],[]]
        Cs = []

        for isurf in xrange(len(self.surfs)):
            if not self._hidden[isurf]:
                ku = self._bse.get_bspline_option('order', isurf, 'u')
                kv = self._bse.get_bspline_option('order', isurf, 'v')
                mu = self._bse.get_bspline_option('num_cp', isurf, 'u')
                mv = self._bse.get_bspline_option('num_cp', isurf, 'v')
                du = numpy.zeros(ku + mu)
                dv = numpy.zeros(kv + mv)
                du[mu:] = 1.0
                dv[mv:] = 1.0
                du[ku-1:mu+1] = numpy.linspace(0, 1, mu-ku+2)
                dv[kv-1:mv+1] = numpy.linspace(0, 1, mv-kv+2)

                ks.append([ku,kv])
                ms.append([mu,mv])
                ds[0].append(du)
                ds[1].append(dv)
                Cs.append(self.surfs[isurf])
        ks = numpy.array(ks)
        ms = numpy.array(ms)

        def write(val, dirID, parID, field, last=False):
            if last:
                self._write('%20.12e;' %(val.real))
            else:
                self._write('%20.12e,' %(val.real))
            field += 1
            if field==3:
                field = 0
                self._write('%9i' %(dirID))
                self._write('P')
                self._write('%7i\n' %(parID))
                parID += 1
            return parID, field

        if filename is None:
            filename = self.name + '.igs'

        self._open_file(filename)
        self._write('                                                                        S      1\n')
        self._write('1H,,1H;,4HSLOT,37H$1$DUA2:[IGESLIB.BDRAFT.B2I]SLOT.IGS;,                G      1\n')
        self._write('17HBravo3 BravoDRAFT,31HBravo3->IGES V3.002 (02-Oct-87),32,38,6,38,15,  G      2\n')
        self._write('4HSLOT,1.,1,4HINCH,8,0.08,13H871006.192927,1.E-06,6.,                   G      3\n')
        self._write('31HD. A. Harrod, Tel. 313/995-6333,24HAPPLICON - Ann Arbor, MI,4,0;     G      4\n')

        dirID = 1
        parID = 1
        for s in range(ks.shape[0]):
            numFields = 4 + ds[0][s].shape[0] + ds[1][s].shape[0] + 4*ms[s,0]*ms[s,1]
            numLines = 2 + numpy.ceil(numFields/3.0)
            for val in [128, parID, 0, 0, 1, 0, 0, 0]:
                self._write('%8i' %(val))
            self._write('00000001')
            self._write('D')
            self._write('%7i\n' %(dirID))
            dirID += 1
            for val in [128, 0, 2, numLines, 0]:
                self._write('%8i' %(val))
            self._write('%32i' %(0))
            self._write('D')
            self._write('%7i\n' %(dirID))
            dirID += 1
            parID += numLines
        nDir = dirID - 1

        dirID = 1
        parID = 1
        for s in range(ks.shape[0]):
            ku = ks[s,0]
            kv = ks[s,1]
            mu = ms[s,0]
            mv = ms[s,1]
            du = ds[0][s]
            dv = ds[1][s]

            for val in [128, mu-1, mv-1, ku-1, kv-1]:
                self._write('%10i,' %(val))  
            self._write('          ')
            self._write('%7i' %(dirID))   
            self._write('P')
            self._write('%7i\n' %(parID))
            parID += 1

            for val in [0, 0, 1, 0, 0]:
                self._write('%10i,' %(val))
            self._write('          ')
            self._write('%7i' %(dirID))   
            self._write('P')
            self._write('%7i\n' %(parID))
            parID += 1

            field = 0
            for i in range(du.shape[0]):
                parID,field = write(du[i], dirID, parID, field)
            for i in range(dv.shape[0]):
                parID,field = write(dv[i], dirID, parID, field)
            for i in range(mu*mv):
                parID,field = write(1.0, dirID, parID, field)
            for j in range(mv):
                for i in range(mu):
                    for k in range(3):
                        parID,field = write(Cs[s][i,j,k], dirID, parID, field)
            parID,field = write(0, dirID, parID, field)
            parID,field = write(1, dirID, parID, field)
            parID,field = write(0, dirID, parID, field)
            parID,field = write(1, dirID, parID, field, last=True)
            if not field==0:
                for i in range(3-field):
                    self._write('%21s' %(' '))
                self._write('%9i' %(dirID))
                self._write('P')
                self._write('%7i\n' %(parID))
                parID += 1

            dirID += 2

        nPar = parID - 1

        self._write('S%7i' %(1))   
        self._write('G%7i' %(4))   
        self._write('D%7i' %(nDir))   
        self._write('P%7i' %(nPar))   
        self._write('%40s' %(' '))   
        self._write('T')
        self._write('%7i\n' %(1))
        self._close_file()
