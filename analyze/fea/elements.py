# -*- coding: utf-8 -*-
"""

"""
#import os, sys
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
from abc import ABCMeta, abstractmethod#, abstractproperty
import numpy as np
from scipy import linalg
from StrEngL.Geometry.Primitives import vector

class Element(object):
    """
    Base class / Interface for all element types.
    """
    __metaclass__ = ABCMeta
    
    #Class variables
    NUM_INT_PTS = None
    NODAL_DOF = None

    def __init__(self, nodes, material_matrix):
        """
        Intialize an Element object and define required attributes.
        
        Arguments:
        nodes -- List of Node objects (or objects that inherit from Node)
        material_matrix -- 3x3 material properties <np.array>
        """
        self.nodes = nodes
        self.material_matrix = material_matrix
#        self.property = propertyObj
        self.num_nodes = len(nodes)
        self.k = None
        self.f = None
        self._body_force = np.zeros(self.num_nodes*self.NODAL_DOF)
        self._traction = np.zeros(self.num_nodes*self.NODAL_DOF)

    @property
    def xLocs(self):
        """
        Return:
        np.array of nodal x-locations
        """
        return np.array([n.x for n in self.nodes])

    @property
    def yLocs(self):
        """
        Return:
        np.array of nodal y-locations
        """
        return np.array([n.y for n in self.nodes])
        
    @property
    def edges(self):
        """
        Create an edge vector object for each element edge.
        
        Return:
        List of vector objects for each element edge
        """
        edges = []
        for i in range(self.num_nodes):
            if i == (self.num_nodes - 1):
                edges.append(vector(
                    self.nodes[0].location - self.nodes[i].location))
            else:
                edges.append(vector(
                    self.nodes[i + 1].location - self.nodes[i].location))
        return edges
        
    @abstractmethod
    def _N(self, xi, eta):
        """
        Return: 
        np.array of shape function results. Output length should match
        number of nodes.
        (must be overidden by subclass)
        """
        raise NotImplementedError

    @abstractmethod    
    def J(self, xi, eta):
        """
        Return: 
        Jacobian matrix.
        (must be overridden by subclass)
        """
        raise NotImplementedError
        
    @abstractmethod
    def B(self, xi, eta):
        """
        Return: 
        Strain-displacement 'B' matrix.
        (must be overriden by subclass)
        """
        raise NotImplementedError
        
#    @abstractproperty
#    def k(self):
#        """
#        Returns to element's stiffness matrix.
#        (must be overriden by subclass)
#        """
#        raise NotImplementedError

    def detJ(self, xi, eta):
        """
        Get the determinate of the jacobian.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        Determinate of the jacobian matrix
        """
        return linalg.det(self.J(xi, eta))

    def x(self, xi, eta):
        """
        Get the x location in global coordinate system.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        """
        return self._N(xi, eta).dot(self.xLocs)

    def y(self, xi, eta):
        return self._N(xi, eta).dot(self.yLocs)
        

class Quad2D(Element):
    """
    A 2-dimensional quadrilateral element with 2x2 gauss integration.
    
    Edge and local node numbering:
    
    3         2
    +---------+
    |    2    |
    |        1|
    |3        |
    |    0    |
    +---------+
    0         1

    Ref. 1: Hughes, T., "The Finite Element Method: Linear Static and Dynamic 
    Finite Element Analysis". 2000. ISBN 0-486-41181-8
    """
    
    #Class variables
    NUM_INT_PTS = 2 #1D, used in each direction
    NODAL_DOF = 2

    def __init__(self, nodes, material_matrix):
        """
        Initialize a Quad2D element.
        
        Arguments:
        nodes -- List of Node objects (or objects that inherit from Node)
        material_matrix -- 3x3 material properties <np.array>
        """
        super(Quad2D, self).__init__(nodes, material_matrix)

    def _N(self, xi, eta):
        """
        Get shape function results.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        1x4 <np.array> shape function results at coordinates
        """
        return 1/4*np.array([(1 - xi)*(1 - eta), (1 + xi)*(1 - eta),
                             (1 + xi)*(1 + eta), (1 - xi)*(1 + eta)])

    def _dNdxi(self, xi, eta):
        """
        Get shape function derivatives wrt xi.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        1x4 <np.array> shape function derivatives at coordinates
        """
        return 1/4*np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])

    def _dNdeta(self, xi, eta):
        """
        Get shape function derivatives wrt eta.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        1x4 <np.array> shape function derivatives at coordinates
        """
        return 1/4*np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])

    def J(self, xi, eta):
        """
        Get Jacobian matrix.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        Jacobian, 2x2 <np.array> at coordinates
        """
        J11 = self._dNdxi(xi, eta).dot(self.xLocs)
        J12 = self._dNdxi(xi, eta).dot(self.yLocs)
        J21 = self._dNdeta(xi, eta).dot(self.xLocs)
        J22 = self._dNdeta(xi, eta).dot(self.yLocs)
        
        return np.array([[J11, J12], [J21, J22]])

    def _dN(self, xi, eta):
        """
        Get shape function derivatives wrt x and y.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        2x1 <np.array> shape function derivatives at coordinates
        """
        J = self.J(xi, eta)
        Jinv = linalg.inv(J)
        return Jinv.dot(np.array([self._dNdxi(xi, eta), 
                                  self._dNdeta(xi, eta)]))

    def B(self, xi, eta):
        """
        Get strain-displacement matrix.
        
        Arguments:
        xi -- x-location in natural (element) coordinates
        eta -- y-location in natural (element) coordinates
        
        Return:
        3x8 <np.array> strian-displacement matrix
        """
        Brow1 = []
        Brow2 = []
        Brow3 = []
        
        dNdx = self._dN(xi, eta)[0]
        dNdy = self._dN(xi, eta)[1]
        
        for i in range(self.num_nodes):
            Brow1.extend([dNdx[i], 0])
            Brow2.extend([0, dNdy[i]])
            Brow3.extend([dNdy[i], dNdx[i]])
            
        return np.array([Brow1,Brow2,Brow3])
                
    def generate_stiffness_matrix(self):
        """
        Overwrites the element stiffness matrix with current node locations
        and material matrix.
        """
        from numpy.polynomial.legendre import leggauss
        
        dof = self.NODAL_DOF
        Nn = self.num_nodes

        self.k = np.zeros((dof*Nn, dof*Nn))

        # Points and weights for Gauss quadrature
        xi, Wi = leggauss(self.NUM_INT_PTS)
        eta, Wj = xi, Wi

        D = self.material_matrix
        B = self.B
        
        # Integrate over the the element (natural coordinates)
        for ii in range(len(xi)):
            for jj in range(len(eta)):
                self.k += (B(xi[ii],eta[jj]).T.dot(D).dot(B(xi[ii],eta[jj]))
                           *Wi[ii]*Wj[jj]*self.detJ(xi[ii], eta[jj]))
                               
    def generate_force_vector(self):
        """
        Builds the element force vector given the current state of applied
        loads.
        
        Return: 
        8x1 <np.array> Element force vector
        """
        
        self.f = self._body_force + self._traction
        
        bc_array = self._get_bc_array()
        
        # Include boundary condition terms ref. [1] pp. 99
        for A in range(self.num_nodes):
            for I in range(self.NODAL_DOF):
                p = self.NODAL_DOF*A + I
                self.f[p] -= self.k[p,:].dot(bc_array)
                
        return self.f + self._nodal_forces()
        
    def _get_bc_array(self):
        """
        Get boundary condition array from element nodes. Each index (dof)
        contains displacement if given, otherwise it is set to zero.
        
        Return: 
        1x8 <np.array> boundary condition array
        """
        bc_array = []
        for node in self.nodes:
            for bc in node.boundary_conditions:
                if bc != None:
                    bc_array.append(bc)
                else:
                    bc_array.append(0)
        return np.array(bc_array)
        
    def _nodal_forces(self):
        """
        Nodal force array where each index (dof) contains scalar force
        value if given, otherwise it is set to zero.
        
        Return: 
        1x8 <np.array> nodal forces
        """
        nodal_forces = np.zeros(self.num_nodes*self.NODAL_DOF)
        for A, node in enumerate(self.nodes):
            #node.force is a 1x2 <np.array>
            nodal_forces[A:A+2] += node.force
        return nodal_forces
        
    def apply_body_force(self, force_vector):
        """
        Overwrites the element body force array with given force vector.
        
        Arguments:
        force_vector -- 1x2 <np.array> body force vector
        """
        from numpy.polynomial.legendre import leggauss
        
        #erase incase of re-write
        self._body_force = np.zeros(self.num_nodes*self.NODAL_DOF)
        
        xi, Wi = leggauss(self.NUM_INT_PTS)
        eta, Wj = xi, Wi
        
        for A in range(self.num_nodes):
            for I in range(self.NODAL_DOF):
            
                p = self.NODAL_DOF*A + I
                
                for ii in range(len(xi)):
                    for jj in range(len(eta)):
                        self._body_force[p] += \
                                (self._N(xi[ii], eta[jj])[A]*force_vector[I]
                                 *Wi[ii]*Wj[jj]*self.detJ(xi[ii], eta[jj]))
                                 
    def apply_traction(self, traction_vector, edge=0):
        """
        Adds the given traction force vector, applied to the given edge, 
        to the element traction force array.
        
        Arguments:
        traction_vector -- 1x2 <np.array> traction force vector
        edge -- edge number (1 : node_0 --> node_1
                             2 : node_1 --> node_2
                             3 : node_2 --> node_3
                             4 : node_3 --> node_0)
        """
        from numpy.polynomial.legendre import leggauss
        
        #erase incase of re-write
#        self._traction = np.zeros(self.num_nodes*self.NODAL_DOF)
        
        if edge == 0: #between nodes 1 & 2
            xi, Wi = leggauss(self.NUM_INT_PTS)
            eta, Wj = [(-1,-1), (1,1)]
        elif edge == 1: #between nodes 2 & 3
            xi, Wi = [(1,1), (1,1)]
            eta, Wj = leggauss(self.NUM_INT_PTS)
        elif edge == 2: #between nodes 3 & 4
            xi, Wi = leggauss(self.NUM_INT_PTS)
            eta, Wj = [(1,1), (1,1)]
        elif edge == 3: #between nodes 4 & 1
            xi, Wi = [(-1,-1), (1,1)]
            eta, Wj = leggauss(self.NUM_INT_PTS)
            
        j = self.edges[edge].magnitude/2.0
                
        for A in range(self.num_nodes):
            for I in range(self.NODAL_DOF):
            
                p = self.NODAL_DOF*A + I
                
                for jj in range(2):
                    self._traction[p] += \
                            (self._N(xi[jj], eta[jj])[A]*traction_vector[I]
                             *Wi[jj]*Wj[jj]*j)
