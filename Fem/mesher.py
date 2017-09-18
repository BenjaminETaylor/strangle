# -*- coding: utf-8 -*-
"""

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from abc import ABCMeta, abstractmethod#, abstractproperty
import numpy as np
from StrEngL.Fem import nodes, elements

def isSquare(aspect_ratio):
    return 0.9 < aspect_ratio < 1.1

def logb(x, base=10):
    return np.log(x)/np.log(base)

class Mesh(object):
    """
    Base class / Interface for all meshing schemes.
    """
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def _build_nodes(self):
        raise NotImplementedError
        
    @abstractmethod
    def _build_elements(self):
        raise NotImplementedError
    
    def generate_mesh(self):
        self._build_nodes()
        self._build_elements()

    def mirror_nodes(self):
        raise NotImplementedError

class BoltedJoint(Mesh):

    def __init__(self, pitch, width, bolt_diameter, sample_distance,
                 material_matrix=None, num_samples=8*4, num_radial_nodes=8):
        """
        @inputs: num_samples must be a multiple of 8
        """
        self.pitch = pitch
        self.width = width
        self.aspect_ratio = pitch/width
        self.bolt_diameter = bolt_diameter
        self.sample_distance = sample_distance
        self.material_matrix = material_matrix
        self.num_samples = num_samples
        self.num_radial_nodes = num_radial_nodes
        self.nodes = []
        self.elements = []

    def _build_nodes(self):
        """
        Generates nodes using polar coordinate system around fastener hole.
        """
        ray_angles = self._generate_ray_angles()

        node_locations = []
        for ray_angle in ray_angles:
            ray_length = self._calculate_ray_length(ray_angle)
            node_locations.extend(self._distribute_nodes(ray_angle,
                                                         ray_length))

#        if not isSquare(self.aspect_ratio):
#            node_locations.extend(self._generate_extension_nodes())

        self._generate_node_objects(node_locations)

    def _build_elements(self):

        rows = self.num_samples
        columns = self.num_radial_nodes
        num_nodes = len(self.nodes)

        for row in range(rows):
            for column in range(columns - 1):

                node1 = row*columns + column
                node2 = node1 + 1
                node3 = node2 + columns
                node4 = node1 + columns

                node3 = node3 if node3 <= num_nodes else node3 - num_nodes - 1
                node4 = node4 if node4 <= num_nodes else node4 - num_nodes - 1

                node_list = [self.nodes[node1],
                             self.nodes[node2],
                             self.nodes[node3],
                             self.nodes[node4]]

                self.elements.append(elements.Quad2D(node_list,
                                                     self.material_matrix))

    def _generate_ray_angles(self):
        """
        @todo: include capability to handle rectangles
        """
        assert self.num_samples % 8 == 0, "num_samples must be multiple of 8!"
        if self.aspect_ratio == 1.0: #square       
            step = 2.0*np.pi/self.num_samples
            return np.arange(0, 2.0*np.pi + step, step)
        elif self.aspect_ratio < 1.0: #pitch < width
            raise NotImplementedError
        elif self.aspect_ratio > 1.0: #pitch > width
            raise NotImplementedError

    def _calculate_ray_length(self, ray_angle):
        """
        @todo: include capability to handle rectangles
        """
        if (0 <= ray_angle <= np.pi/4) or (3*np.pi/4 <= ray_angle <= 5*np.pi/4) \
        or (7*np.pi/4 <= ray_angle <= 2*np.pi):

            return self.width/abs(np.cos(ray_angle))/2

        elif (np.pi/4 <= ray_angle <= 3*np.pi/4) \
        or (5*np.pi/4 <= ray_angle <= 7*np.pi/4):

            return self.pitch/abs(np.sin(ray_angle))/2

    def _distribute_nodes(self, ray_angle, ray_length):
        """
        @returns: list of polar (radius, angle) pairs for each node on ray

        Uses the following equation to solve for step length used in the
        numpy.logspace() function.

        endpoint = base**stop = base**((num_nodes - 1)*step)

        #setting base = 10,

        np.log10(endpoint) = (num_nodes - 1)*step

        #and solving for step gives,

        step = np.log10(endpoint)/(num_nodes - 1)

        #then use step to solve for stop

        stop = (num_nodes - 1)*step
        
        @todo: fix 1st row element (at bolt hole radius) aspect ratio when
        high values of num_samples are chosen.
        """
        bolt_radius = self.bolt_diameter/2.0
        first_element_width = 2.0*self.sample_distance
        second_row_radius = bolt_radius + first_element_width
        node_list = [(bolt_radius, ray_angle, 0)]

        #Use logarihmic distribution for remaining nodes
        offset = 1.0 - second_row_radius #1.0 is used because start=0 below
        endpoint = ray_length + offset
        num_nodes_remaining = self.num_radial_nodes - 1 #one node already defined above
        step = np.log10(endpoint)/(num_nodes_remaining - 1)
        stop = (num_nodes_remaining - 1)*step

        distribution = np.logspace(0, stop, num=num_nodes_remaining)
        distribution -= offset
        node_list.extend([(node_radius, ray_angle, 0) for 
                         node_radius in distribution])

        return node_list

    def _generate_node_objects(self, node_locations):
        self.nodes = list(map(nodes.PolarNode,
                              node_locations,
                              list(range(len(node_locations)))))
                              
class RectangularPlate(Mesh):
    
    def __init__(self, length, width, material_matrix=None, num_elements=10):
        self.length = length
        self.width = width
        self.material_matrix = material_matrix
        self.num_elements = num_elements
        self.nodes = []
        self.elements = []
    
    @property
    def _xpts(self):
        return np.linspace(0, self.length, self.num_elements + 1)
        
    @property
    def _ypts(self):
        return np.linspace(0, self.width, self.num_elements + 1)
        
    def _build_nodes(self):
        for ii,x in enumerate(self._xpts):
            for jj,y in enumerate(self._ypts):
                nID = ii*len(self._xpts) + jj
                self.nodes.append(nodes.CartesianNode(np.array([x,y,0]), nID))
                
    def _build_elements(self):
        for ii in range(len(self._xpts) - 1):
            for jj in range(len(self._ypts) - 1):
                ID = ii*len(self._xpts) + jj
                n1 = self.nodes[ID]
                n2 = self.nodes[ID + len(self._xpts)]
                n3 = self.nodes[ID + len(self._xpts) + 1]
                n4 = self.nodes[ID + 1]
                nodeList = [n1, n2, n3, n4]
                self.elements.append(elements.Quad2D(nodeList,
                                                     self.material_matrix))