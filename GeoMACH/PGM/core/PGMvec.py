"""
GeoMACH PGM vec class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
import scipy.sparse


class PGMvec(object):
    """ PGM vector class """

    def __init__(self, name, storing_objects, computing_objects):
        size = 0
        for obj in storing_objects:
            size += numpy.prod(obj.get_vec_shape(name))

        data = numpy.zeros(size)
        inds = numpy.array(numpy.linspace(0, size-1, size), int)

        start, end = 0, 0
        for obj in storing_objects:
            shape = obj.get_vec_shape(name)

            end += numpy.prod(shape)
            obj_data = data[start:end].reshape(shape, order='F')
            obj_inds = inds[start:end].reshape(shape, order='F')
            obj.initialize_vec(name, obj_data, obj_inds)
            start += numpy.prod(shape)

        self._name = name
        self._objects = computing_objects
        self._size = size
        self._data = data
        self._inds = inds
        self._col_size = 0

    def __call__(self):
        return self._data

    def get_size(self):
        return self._size

    def set_col_size(self, size):
        self._col_size = size

    def compute(self):
        list_vals = [numpy.array([])]
        list_rows = [numpy.array([], int)]
        list_cols = [numpy.array([], int)]
        for obj in self._objects:
            vals, rows, cols = obj.compute(self._name)
            list_vals.extend(vals)
            list_rows.extend(rows)
            list_cols.extend(cols)

        vals = numpy.concatenate(list_vals)
        rows = numpy.concatenate(list_rows)
        cols = numpy.concatenate(list_cols)

        jac = scipy.sparse.csr_matrix((vals, (rows, cols)),
                                      shape=(self._size, 
                                             self._col_size))

        return jac
