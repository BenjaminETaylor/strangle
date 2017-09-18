"""
===============================================================================
    StrEngL.Geometry.Components
    Copyright (C) 2015  Benjamin E. Taylor

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

class matrix(np.matrix):
    pass

class vector():
    def __init__(self, components):
        assert len(components) == 3, "Vector must be 3D!"
        self._array = np.array(map(float,components))
        
    def __add__(self, v):
        return self._array + v._array
        
    def __sub__(self, v):
        return self._array - v._array
        
    def __mul__(self, x):
        if isinstance(x, vector):
            return self._array.dot(x._array)
        else: #assume scalar
            return self._array*x
            
    def __truediv__(self, s):
        return self._array/s
        
    def __floordiv__(self, s):
        return self._array//s
        
    def __neg__(self):
        return self._array*-1
        
    def __lt__(self, vec):
        return self.mag() < vec.mag()
        
    def __eq__(self, vec):
        test = self.unit() == vec.unit()
        return all([self.mag() == vec.mag(), test.all()])
        
    def __ne__(self, vec):
        test = self.unit() == vec.unit()
        return not all([self.mag() == vec.mag(), test.all()])
        
    def __gt__(self, vec):
        return self.mag() > vec.mag()
        
    def __repr__(self):
        return repr(self._array).replace('array', 'vector')
        
    def __getitem__(self, key):
        return self._array.__getitem__(key)
        
    @property
    def T(self):
        return self._array.T
    
    @property
    def magnitude(self):
        return np.linalg.norm(self._array)
    
    @property    
    def unit(self):
        return self._array/self.mag()
        
    def component(self, v):
        return np.dot(self._array,v._array)/v.mag()
        
    def project(self, v):
        return self.comp(v)*v.unit()
        
    def cross(self, v):
        return  vector(np.cross(self._array,v._array))
        
    def angle(self, v):
        return np.arccos(np.dot(self._array,v._array)/(self.mag()*v.mag()))

class line():
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

class plane():
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
        return line(self.n, p0).intersect_plane(self)
