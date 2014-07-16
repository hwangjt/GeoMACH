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



class BSEvecUns(BSEvec):

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



class BSEvecStr(BSEvec):

    def __init__(self, name, size, ndim, surf_sizes, hidden):
        super(BSEvecStr, self).__init__(name, size, ndim, hidden)

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
