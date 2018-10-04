# -*- coding: utf-8 -*-
"""

"""
import numpy as np

def removeRow(arr, row):
    """
    Removes <row(index)> in <arr(2D array)>.
    """
    top, row, bot = np.split(arr, [row, row + 1])
    newarr = np.vstack((top, bot))
    return newarr, row

class Assembly(object):
    """
    Create global stiffness matrix and force vector.
    
    Ref. 1: Hughes, T., "The Finite Element Method: Linear Static and Dynamic 
    Finite Element Analysis". pp. 92-98
    """
    
    def __init__(self, mesh):
        """
        Initialize an Assembly object.
        
        Arguments:
        mesh -- instance of subclass from Mesh interface
        """
        self.mesh = mesh
        self.K = None
        self.F = None
        self._global_eqns = {}
        
    def assemble(self):
        """
        Assemble the stiffness matrix and force vector.
        """
        self._build_assembly_array()
        self._assemble_stiffness_matrix()
        self._assemble_force_vector()
        
    def _build_assembly_array(self):
        """
        Populates ID array as defined in ref. [1] pp. 94. 
        """
        eqn_num = 0
        for node in self.mesh.nodes:
            self._global_eqns[node.id] = []
            for bc in node.boundary_conditions:
                if bc == None:
                    self._global_eqns[node.id].append(eqn_num)
                    eqn_num += 1
                else:
                    self._global_eqns[node.id].append(None)
                    
        # assign space for stiffness matrix and force vector
        self.K = np.zeros((eqn_num,eqn_num))
        self.F = np.zeros(eqn_num)
                    
    def _get_eqn_numbers(self, element):
        """
        Find global equation numbers for each dof at nodes of element.
        
        Arguments:
        element -- ojbect conforming to the Element interface
        
        Return:
        List of equation numbers
        """
        eqn_numbers = []
        for node in element.nodes:
            eqn_numbers.extend(self._global_eqns[node.id])
        return eqn_numbers
        
    def _assemble_stiffness_matrix(self): 
        """
        Assemble the global stiffness matrix.
        """
        for element in self.mesh.elements:
            eqn_numbers = self._get_eqn_numbers(element)
            element.generate_stiffness_matrix()
            
            for I, eqnI in enumerate(eqn_numbers):
                for J, eqnJ in enumerate(eqn_numbers):
                    if (eqnI != None) and (eqnJ != None):
                        self.K[eqnI, eqnJ] += element.k[I, J]
                        
    def _assemble_force_vector(self):
        """
        Assemble the global force vector.
        """
        for element in self.mesh.elements:
            eqn_numbers = self._get_eqn_numbers(element)
            element.generate_force_vector()
            
            for I, eqnI in enumerate(eqn_numbers):
                if (eqnI != None):
                    self.F[eqnI] += element.f[I]