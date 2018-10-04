# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:38:09 2016

@author: Ben
"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import numpy as np
#import scipy as sp

class Node(object):

    def __init__(self, location, ID, coord=0):
        self.location = location
        self.coord = coord
        self.id = ID
        self.boundary_conditions = [None, None]
        self.force = np.array([0,0])

    @property
    def location(self):
        return self.__location

    @location.setter
    def location(self, array):
        self.__location = np.array(array)


class CartesianNode(Node):

    def __init__(self, location, ID, coord=0):
        super(CartesianNode, self).__init__(location, ID, coord=coord)

    @property
    def x(self):
        return self.location[0]

    @property
    def y(self):
        return self.location[1]

    @property
    def z(self):
        return self.location[2]

class PolarNode(Node):

    def __init__(self, location, ID, coord=0):
        super(PolarNode, self).__init__(location, ID, coord=coord)

    @property
    def location(self):
        return self.__location

    @location.setter
    def location(self, array):
        """
        x = r cos(theta)
        y = r sin(theta)
        z = z
        """
        self.__location = np.array([array[0]*np.cos(array[1]),
                                    array[0]*np.sin(array[1]), array[2]])

    @property
    def r(self):
        return np.sqrt(self.x**2 + self.y**2)

    @property
    def theta(self):
        return np.arctan(self.y/self.x)

    @property
    def x(self):
        return self.location[0]

    @property
    def y(self):
        return self.location[1]

    @property
    def z(self):
        return self.location[2]