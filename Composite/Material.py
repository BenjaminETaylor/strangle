# -*- coding: utf-8 -*-
"""
===============================================================================
    StrEngL.Composite.Material
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

class Material():
    pass

class Orthotropic(Material):
    """ defines an orthotropic material. \n
    (E1, E2, E3, nu12, nu23, nu13, G12, G23, G13, alpha1=0, alpha2=0, alpha3=0,
    Xt=0, Yt=0, Zt=0, Xc=0, Yc=0, Zc=0, S12=0, S23=0, S13=0)"""
    def __init__(self, E1=0, E2=0, E3=0, nu12=0, nu23=0, nu13=0, 
                 G12=0, G23=0, G13=0, alpha1=0, alpha2=0, alpha3=0, 
                 Xt=0, Yt=0, Zt=0, Xc=0, Yc=0, Zc=0, S12=0, S23=0, S13=0, 
                 e1t=0, e2t=0, e3t=0, e1c=0, e2c=0, e3c=0, e12=0, e23=0, e13=0):
        """enter compression strengths as negative values!"""
        # material stiffness properties
        self.E1 = E1
        self.E2 = E2
        self.E3 = E3
        self.nu12 = nu12
        self.nu23 = nu23
        self.nu13 = nu13
        self.G12 = G12
        self.G23 = G23
        self.G13 = G13
        self.nu21 = E2/E1 * nu12
        self.nu31 = E3/E1 * nu13
        self.nu32 = E3/E2 * nu23
        # Thermal expansion coefficients
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.alpha3 = alpha3
        # axial strength properties (t=tension, c=compression)
        self.Xt = Xt
        self.Yt = Yt
        self.Zt = Zt
        self.Xc = Xc
        self.Yc = Yc
        self.Zc = Zc
        # shear strength properties
        self.S12 = S12
        self.S23 = S23
        self.S13 = S13
        # axial strain allowables (t=tension, c=compression)
        self.e1t = e1t
        self.e2t = e2t
        self.e3t = e3t
        self.e1c = e1c
        self.e2c = e2c
        self.e3c = e3c
        # shear strain allowables
        self.e12 = e12
        self.e23 = e23
        self.e13 = e13

    def Sij(self):
        """ returns full (6x6) compliance matrix."""
        Sij = np.zeros((6,6))
        Sij[0,0] = 1/self.E1
        Sij[0,1] = Sij[1,0] = -self.nu21/self.E2
        Sij[0,2] = Sij[2,0] = -self.nu31/self.E3
        Sij[1,2] = Sij[2,1] = -self.nu32/self.E3
        Sij[1,1] = 1/self.E2
        Sij[2,2] = 1/self.E3
        Sij[3,3] = 1/self.G23
        Sij[4,4] = 1/self.G13
        Sij[5,5] = 1/self.G12
        return Sij

    def Cij(self):
        """ returns full (6x6) elastic constants matrix."""
        return np.linalg.inv(self.Sij())

    def sij(self):
        """ returns (3x3) 2D plane stress compliance matrix."""
        Sij = self.Sij()
        return np.array([[Sij[0,0], Sij[0,1], Sij[0,5]],
                         [Sij[1,0], Sij[1,1], Sij[1,5]],
                         [Sij[5,0], Sij[5,1], Sij[5,5]]])

    def Qij(self):
        """ returns (3x3) 2D plane stress elastic constants matrix."""
        return np.linalg.inv(self.sij())

    def alpha(self):
        """ returns (3x1) 2D plane stress thermal expansion coefficient
        vector."""
        return np.array([[self.alpha1],[self.alpha2],[0]])


class TransverseIso(Orthotropic):
    """ defines a transversely isotropic material. \n
    (E1, E2, G12, nu12, nu23, alpha1=0, alpha2=0,
    Xt=0, Yt=0, Xc=0, Yc=0, S12=0, S23=0,
    e1t=0, e2t=0, e1c=0, e2c=0, e12=0, e23=0,)"""
    def __init__(self, E1=0, E2=0, G12=0, nu12=0.3, nu23=0.3, 
                 alpha1=0, alpha2=0, Xt=0, Yt=0, Xc=0, Yc=0, S12=0, S23=0,
                 e1t=0, e2t=0, e1c=0, e2c=0, e12=0, e23=0):
        """enter compression strengths as negative values!"""
        G23 = E2 / (2*(1 + nu23))
        Orthotropic.__init__(self, E1=E1, E2=E2, E3=E2, nu12=nu12, nu23=nu23, 
                             nu13=nu12, G12=G12, G23=G23, G13=G12,
                             alpha1=alpha1, alpha2=alpha2, alpha3=alpha2, 
                             Xt=Xt, Yt=Yt, Zt=Yt, Xc=Xc, Yc=Yc, Zc=Yc,
                             S12=S12, S23=S23, S13=S12, e1t=e1t, e2t=e2t, 
                             e3t=e2t, e1c=e1c, e2c=e2c, e3c=e2c, e12=e12, 
                             e23=e23, e13=e12)


class Isotropic(TransverseIso):
    """ defines an isotropic material. \n
    (E, nu, alpha=0, Ft=0, Fc=0, Fs=0, et=0, ec=0, es=0)"""
    def __init__(self, E=0, nu=0, alpha=0, Ft=0, Fc=0, Fs=0, et=0, ec=0, es=0):
        """enter compression strengths as negative values!"""
        G = E / (2*(1 + nu))
        TransverseIso.__init__(self, E1=E, E2=E, G12=G, nu12=nu, nu23=nu, 
                               alpha1=alpha, alpha2=alpha, Xt=Ft, Yt=Ft,
                               Xc=Fc, Yc=Fc, S12=Fs, S23=Fs, e1t=et, e2t=et, 
                               e1c=ec, e2c=ec, e12=es, e23=es)
                               
                               
