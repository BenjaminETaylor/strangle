"""
===============================================================================
    Copyright (C) Benjamin E. Taylor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
===============================================================================
"""
import numpy as np


class Coord:

    def __init__(self, origin, cosines):
        """
        :param origin: <np.array> (3x1)
        :param cosines: <np.array> (3x3) direction cosines (rotation matrix)
        """
        self.origin = np.array(origin)
        self.cosines = np.array(cosines)

    def convert_vector(self, vector, def_coord):
        """
        Converts a vector to this coordinate system.
        :param vector: <np.array> (3x1)
        :param def_coord: <strengl.coord>
        :return: <np.array> (3x1) new vector
        """
        vec = np.array(vector)
        rotation_matrix = self.cosines.dot(def_coord.cosines)
        return vec.dot(rotation_matrix)

    def convert_point(self, point, def_coord):
        """
        Converts a point to this coordinate system.
        :param point: <np.array> (3x1) point location
        :param def_coord: <strengl.coord>
        :return: new point location
        """
        translation = def_coord.origin - self.origin
        return translation + point


class Tensor:

    def __init__(self, tensor):
        """
        :param tensor: <np.array> (3x3)
        Tensor form: [[t11, t12, t13],
                      [t21, t22, t23],
                      [t31, t32, t33]]
         """
        self.tensor = np.array(tensor)
        assert self.tensor.shape == (3, 3), "Tensor must be 3x3."

    @property
    def is_symmetric(self):
        return (self.tensor.transpose() == self.tensor).all()

    def rotate(self, rotation_matrix):
        """
        :param rotation_matrix: <np.array> (3x3)
        :return: <np.array> (3x3) rotated strain tensor
        """
        return np.dot(rotation_matrix, self.tensor.dot(rotation_matrix.transpose()))

    def to_vector(self):
        if self.is_symmetric:
            vec = self.tensor.flatten()
            return np.array([vec[0], vec[4], vec[7], vec[5], vec[2], vec[2]]).transpose()
        else:
            raise ValueError("Tensor cannot become a vector. It is not symmetric.")


class Vector:

    @staticmethod
    def magnitude(vec):
        return np.linalg.norm(vec)
    
    @staticmethod
    def unit(vec):
        return vec/np.linalg.norm(vec)

    @staticmethod
    def component(vec1, vec2):
        unit_vec2 = vec2/np.linalg.norm(vec2)
        return np.dot(vec1, unit_vec2)

    @staticmethod
    def project(vec1, vec2):
        unit_vec2 = vec2/np.linalg.norm(vec2)
        comp_vec = np.dot(vec1, unit_vec2)
        return comp_vec*unit_vec2

    @staticmethod
    def angle(vec1, vec2):
        vec1_mag = np.linalg.norm(vec1)
        vec2_mag = np.linalg.norm(vec2)
        return np.arccos(np.dot(vec1, vec2)/(vec1_mag*vec2_mag))


class Line:
    def __init__(self, m, x0):
        self.m = m.unit()
        self.x0 = x0
        
    def __repr__(self):
        return 'x=mt+x0 = %st + %s' % (self.m, self.x0)
        
    def __call__(self, t):
        return self.m*t + self.x0
        
    def intersect_plane(self, plane):
        t = plane.n*(plane.x0 - self.x0)/(plane.n*self.m)
        return self(t)


class Plane:
    def __init__(self, n, x0):
        #n = normal vector()
        #x0 = location vector()
        self.n = n.unit()
        self.x0 = x0
        self.a = None
        self.b = None
        
    def __repr__(self):
        return '0=(x-x0)*n = (x - %s)*%s' % (self.x0, self.n)
        
    def paramaterize(self, vec, p0):
        #x = x0+sa+tb
        #where: s,t => scalars
        self.a = self.projVector(vec).unit
        self.b = self.n.cross(self.a)
        
    def projVector(self, v):
        #try: proj = v-self.n.dot(v)*self.n
        try: proj = v - v.proj(self.n)
        except: proj = v
        return proj
        
    def projPoint(self, p0):
        return Line(self.n, p0).intersect_plane(self)
