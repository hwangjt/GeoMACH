"""
GeoMACH component class
John Hwang, July 2014
"""
# pylint: disable=E1101
from __future__ import division
import numpy
from collections import OrderedDict

from GeoMACH.PGM.core.PGMproperty import PGMproperty


class PGMcomponent(object):
    """ Base class for aircraft components """

    def __init__(self):
        """
        Attributes
        ----------
        **faces** : OrderedDict
           Dictionary of ``PGMface`` objects
        **props** : OrderedDict
           Dictionary of ``PGMprop`` objects
        **_name** : str
           Name of component assigned by user
        **_num** : str
           Global index of this component object
        **_bse** : BSEmodel
           Pointer to configuration's BSE object
        **_num_surf** : dictionary of ``int``s
           User-chosen number of surfaces in each direction
           with keys ['x', 'y', 'z']
        **_shapes** : dictionary with same keys as ``faces``
           Normalized shape control points in section frame           
        """
        self.faces = OrderedDict()
        self.props = OrderedDict()

        self._name = None
        self._num = None
        self._bse = None
        self._num_surf = {}
        self._shapes = {}

    def initialize_props(self):
        """ Adds the *X*, *Y*, and *Z* shape properties """
        props = self.props
        for name in self.faces:
            props['shX', name] = PGMproperty()
            props['shY', name] = PGMproperty()
            props['shZ', name] = PGMproperty()

    def assemble_sizes(self, bse):
        """ 
        Calls on faces to store the numbers of control points
        and passes this information to properties.

        Arguments
        ---------
        **bse** : ``BSEmodel``
           BSE model instance of parent configuration object.
           The value is ``None`` the first time this method is called. 
        """
        for face in self.faces.values():
            face.assemble_sizes(bse)

        props = self.props
        for name in self.faces:
            num_u = self.faces[name]._num_cp_total['u']
            num_v = self.faces[name]._num_cp_total['v']
            props['shX', name].assemble_sizes(num_u, num_v)
            props['shY', name].assemble_sizes(num_u, num_v)
            props['shZ', name].assemble_sizes(num_u, num_v)
            self._shapes[name] = numpy.zeros((num_u, num_v, 3), 
                                             order='F')

    def initialize_bse(self, surfs_list, surf_index):
        """
        Collects initial list of surfaces from each face
        and initializes surface indices for each face.

        Arguments
        ---------
        **surfs_list** : ``list`` of ``numpy.ndarray``s, ``float``[ni, nj, 3]
           List of initial surfaces computed by PGM (df_surf)
        **surf_index** : ``int``
           Smallest available surface index

        Returns
        -------
        **surf_index** : ``int``
           Next available after indices have been assigned
        """
        for face in self.faces.values():
            surf_index = face.initialize_bse(surfs_list, surf_index)
        return surf_index

    def set_info(self, name, num, bse):
        """
        Assigns component name and index as well as BSE object
        and does the same for faces
        """
        self._name = name
        self._num = num
        self._bse = bse

        face_num = 0
        for face_name in self.faces:
            self.faces[face_name].set_info(face_name, face_num, bse)
            face_num += 1

    def set_diff(self):
        """ Sets differentiability options for surface boundaries """
        pass

    def set_hidden_surfaces(self):
        """ Hides surfaces replaced by junctions """
        pass
