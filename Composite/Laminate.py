# -*- coding: utf-8 -*-
"""
===============================================================================
    StrEngL.Composite.Laminate
    Copyright (C) 2016  Benjamin E. Taylor

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

#from __future__ import division
import numpy as np

# Private Functions --------
def _c(theta): return np.cos(np.deg2rad(theta))
def _s(theta): return np.sin(np.deg2rad(theta))

def _T_eps(a):
    return np.array([[_c(a)**2, _s(a)**2, _s(a)*_c(a)],
                    [_s(a)**2, _c(a)**2, -_s(a)*_c(a)],
                    [-2*_s(a)*_c(a),2*_s(a)*_c(a),_c(a)**2-_s(a)**2]])

def _T_sig(a):
    return np.array([[_c(a)**2, _s(a)**2, 2*_s(a)*_c(a)],
                    [_s(a)**2, _c(a)**2, -2*_s(a)*_c(a)],
                    [-_s(a)*_c(a),_s(a)*_c(a),_c(a)**2-_s(a)**2]])
#-----------------------
    
# Public Functions -------
def buildLam(matObj, thickness, angles, offset=0):
    """generates a laminate with minimal input"""
    plies = []
    for angle in angles:
        plies.append(Ply(matObj,thickness,angle))
    return Laminate(plies, offset=offset)
# ----------------------

class Ply():
    """defines a ply object"""
    def __init__(self, matObj, thickness, angle):
        self.mat = matObj
        self.t = thickness
        self.angle = angle

    def Qbar(self):
        """returns a rotated (3x3) plane stress elastic constants matrix."""
        Q = _T_eps(self.angle).T.dot(self.mat.Qij())
        return Q.dot(_T_eps(self.angle))

    def sbar(self):
        """returns a rotated (3x3) plane stress compliance matrix."""
        s = _T_sig(self.angle).T.dot(self.mat.sij())
        return s.dot(_T_sig(self.angle))

    def alphaBar(self):
        """returns global thermal expansion coefficient vector (alpha bar)"""
        return _T_sig(self.angle).T.dot(self.mat.alpha())
        
    def rotate(self, angle):
        """rotates ply by angle"""
        newAngle = self.angle + angle
        while newAngle >= 180:
            newAngle -= 180
        while newAngle <= -180:
            newAngle += 180
        self.angle = newAngle

class Laminate():
    """defines a laminate object"""
    def __init__(self, plies, offset=0):
        """plies = List of ply objects"""
        self.plies = plies
        self.z = self.__plyCoordinates(offset=offset)
        self.__A = None
        self.__B = None
        self.__D = None

    def __plyCoordinates(self, offset=0):
        """returns list of ply coordinates"""
        self.H = sum([ply.t for ply in self.plies])
        z = [-self.H/2 + offset]
        for k,ply in enumerate(self.plies):
            z.append(z[k] + ply.t)
        return z

    def rotate(self, angle):
        """rotates plies in laminate by angle"""
        for ply in set(self.plies):
            ply.rotate(angle)
        self.__A, self.__B, self.__D = None, None, None

    @property    
    def A(self):
        """returns A matrix"""
        if self.__A is None:        
            A = np.zeros((3,3))
            for k,ply in enumerate(self.plies):
                A += ply.Qbar() * (self.z[k+1] - self.z[k])
            self.__A = A
        return self.__A

    @property    
    def B(self):
        """returns B matrix"""
        if self.__B is None:
            B = np.zeros((3,3))
            for k,ply in enumerate(self.plies):
                B += 1/2 * ply.Qbar() * (self.z[k+1]**2 - self.z[k]**2)
            self.__B = B
        return self.__B

    @property    
    def D(self):
        """returns D matrix"""
        if self.__D is None:
            D = np.zeros((3,3))
            for k,ply in enumerate(self.plies):
                D += 1/3 * ply.Qbar() * (self.z[k+1]**3 - self.z[k]**3)
            self.__D = D
        return self.__D

    @property    
    def ABD(self):
        """returns ABD matrix"""
        return np.array(np.bmat([[self.A, self.B],
                                 [self.B, self.D]]))

    @property    
    def abd(self):
        """returns abd matrix"""
        return np.linalg.inv(self.ABD)

    @property    
    def a(self):
        """returns a matrix"""
        return self.abd[:3,:3]

    @property    
    def b(self):
        """returns b matrix"""
        return self.abd[3:,:3]

    @property    
    def d(self):
        """returns d matrix"""
        return self.abd[3:,3:]

    @property    
    def Ex(self):
        """returns laminate Ex average modulus (assumes symmetry)"""
        return 1/(self.H*self.a[0,0])

    @property    
    def Ey(self):
        """returns laminate Ey average modulus (assumes symmetry)"""
        return 1/(self.H*self.a[1,1])

    @property    
    def Gxy(self):
        """returns laminate Gxy average modulus (assumes symmetry)"""
        return 1/(self.H*self.a[2,2])

    @property    
    def nuxy(self):
        """returns laminate nuxy effective Poisson's ratio
        (assumes symmetry)"""
        return -self.a[0,1]/self.a[0,0]

    @property    
    def nuyx(self):
        """returns laminate nuyx effective Poisson's ratio
        (assumes symmetry)"""
        return -self.a[0,1]/self.a[1,1]

    def Nt(self, deltaT):
        """returns thermal in-plane loads"""
        Nt = np.zeros((3,1))
        for ply in self.plies:
            Nt += deltaT * ply.Qbar().dot(ply.alphaBar()) * ply.t
        return Nt

    def Mt(self, deltaT):
        """returns thermal moments"""
        Mt = np.zeros((3,1))
        for k,ply in enumerate(self.plies):
            Mt += 1/2 * deltaT * ply.Qbar().dot(ply.alphaBar()) \
                * (self.z[k+1]**2 - self.z[k]**2)
        return Mt

    def alpha(self):
        """returns laminate average thermal expansion coefficients"""
        alpha = np.zeros((3,1))
        for ply in self.plies:
            alpha += self.a().dot(ply.Qbar()).dot(ply.alphaBar()) * ply.t
        return alpha