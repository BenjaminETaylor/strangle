# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 07:56:23 2018

@author: Ben
"""
#import abc
import numpy as np

class Hole:
    """
    Class for defining a hole in an anisotropic plate for stress
    calculations.
        
    References:
        1. "Stress distribution and strength prediction of composite 
        laminates with multiple holes" by Esp, Brian, Ph.D., 
        The University of Texas at Arlington, 2007, 177; 3295219
        2. Lekhnitskii, S.G. (1968), "Anisotropic Plates", 2nd Ed., 
        Translated from the 2nd Russian Ed. by S.W. Tsai and Cheron, 
        Bordon and Breach.
        3. Garbo, S.P. and Ogonowski, J.M., "Effect of Variances and
        Manufacturing Tolerances on the Design Strength and Life of
        Mechanically Fastened Composite Joints", AFWAL-TR-81-3041,
        Volumes 1, 2 and 3, April 1981
        4. J. P. Waszczak and T. A. Cruse, "A Synthesis Procedure for
        Mechanically Fastened Joints in Advanced Composite Materials",
        AFML-TR-73-145, Volume II, 1973.

    """
    def __init__(self, r, lam):
        """
        Class contructor.
        
        :param r: hole radius
        :param lam: <strengl.composite.laminate>
        """
        self.r = r
        self.lam = lam
        self.mu1, self.mu2, self.mu1_bar, self.mu2_bar = self.roots()
        
    def roots(self):
        """
        Finds the roots to the characteristic equation [Eq. A.2, Ref. 1]
        """
        a11 = self.lam.a[0,0]
        a16 = self.lam.a[0,2]
        a12 = self.lam.a[0,1]
        a66 = self.lam.a[2,2]
        a26 = self.lam.a[1,2]
        a22 = self.lam.a[1,1]
        
        roots = np.roots([a11, -2*a16, (2*a12 + a66), -2*a26, a22])
        
        if np.imag(roots[0]) >= 0.0:
            mu1 = roots[0]
            mu1_bar = roots[1]
        elif np.imag(roots[1]) >= 0.0:
            mu1 = roots[1]
            mu1_bar = roots[0]
        else:
            raise ValueError("mu1 cannot be solved")
            
        if np.imag(roots[2]) >= 0.0:
            mu2 = roots[2]
            mu2_bar = roots[3]
        elif np.imag(roots[3]) >= 0.0:
            mu2 = roots[3]
            mu2_bar = roots[2]
        else:
            raise ValueError("mu2 cannot be solved")
            
        return mu1, mu2, mu1_bar, mu2_bar
        
    def ksi_1(self, z1):
        """
        Calculates the first mapping parameter [Eq. A.4 & Eq. A.5, Ref. 1]
        """
        mu1 = self.mu1
        a = self.r
        b = self.r
        
        ksi_1_pos = (z1 + np.sqrt(z1*z1 - a*a - mu1*mu1*b*b))/(a - 1j*mu1*b)
        ksi_1_neg = (z1 - np.sqrt(z1*z1 - a*a - mu1*mu1*b*b))/(a - 1j*mu1*b)
        
        #ksi_1 = []
        #sign_1 = []
        
        #for ii in range(len(ksi_1_pos)):
        if np.abs(ksi_1_pos) >= 1.0 - 1.0e-15:
            ksi_1 = ksi_1_pos
            sign_1 = 1
        elif np.abs(ksi_1_neg) >= 1.0 - 1.0e-15:
            ksi_1 = ksi_1_neg
            sign_1 = -1
        else:
            raise ValueError(
            "ksi_1 unsolvable:\n ksi_1_pos={0}, ksi_1_neg={1}".format(
            ksi_1_pos, ksi_1_neg))
        
        #return np.array(ksi_1), np.array(sign_1)
        return ksi_1, sign_1

    def ksi_2(self, z2):
        """
        Calculates the second mapping parameter [Eq. A.4 & Eq. A.5, Ref. 1]
        """
        mu2 = self.mu2
        a = self.r
        b = self.r
        
        ksi_2_pos = (z2 + np.sqrt(z2*z2 - a*a - mu2*mu2*b*b))/(a - 1j*mu2*b)
        ksi_2_neg = (z2 - np.sqrt(z2*z2 - a*a - mu2*mu2*b*b))/(a - 1j*mu2*b)
        
        #ksi_2 = []
        #sign_2 = []        
        
        #for ii in range(len(ksi_2_pos)):
        if np.abs(ksi_2_pos) >= 1.0 - 1.0e-15:
            ksi_2 = ksi_2_pos
            sign_2 = 1
        elif np.abs(ksi_2_neg) >= 1.0 - 1.0e-15:
            ksi_2 = ksi_2_neg
            sign_2 = -1
        else:
            raise ValueError(
            "ksi_1 unsolvable:\n ksi_1_pos={0}, ksi_1_neg={1}".format(
            ksi_2_pos, ksi_2_neg))
        
        #return np.array(ksi_2), np.array(sign_2)
        return ksi_2, sign_2
    
    #@abc.abstractclassmethod()
    def phi_1_prime(self, z1):
        raise NotImplementedError("You must implement this function.")
    
    #@abc.abstractclassmethod()
    def phi_2_prime(self, z2):
        raise NotImplementedError("You must implement this function.")
        
    def stress(self, x, y):
        """
        Calculates the stress at point (x, y) in the plate. 
        [Eq. 8.2, Ref. 2]
        """
        mu1 = self.mu1
        mu2 = self.mu2
        
        z1 = x + mu1*y
        z2 = x + mu2*y
        
        phi_1_prime = self.phi_1_prime(z1)
        phi_2_prime = self.phi_2_prime(z2)
        
        sx = 2.0*np.real(mu1*mu1*phi_1_prime + mu2*mu2*phi_2_prime)
        sy = 2.0*np.real(phi_1_prime + phi_2_prime)
        sxy = -2.0*np.real(mu1*phi_1_prime + mu2*phi_2_prime)
        
        return np.array([sx, sy, sxy])

    

class Unloaded(Hole):
    """
    Class for defining an unloaded hole in an anisotropic homogeneous plate
    with farfield forces (Nx, Ny, Nxy) [force / unit length] applied.
    """
    
    def __init__(self, r, lam, bypass):
        """
        Class constructor.
        
        :param r: hole radius
        :param lam: <strengl.composite.Laminate>
        :param bypass: [Nx, Ny, Nxy] force / unit length
        """
        super().__init__(r, lam)
        self.applied_stress = np.array(bypass)/lam.H
        
    def alpha(self):
        """
        Calculates the alpha loading term. [Eq. A.7, Ref. 1]
        """
        sy = self.applied_stress[1]
        sxy = self.applied_stress[2]
        r = self.r

        return 1j*sxy*r/2 - sy*r/2 
        
    def beta(self):
        """
        Calculates the beta loading term. [Eq. A.7, Ref. 1]
        """
        sx = self.applied_stress[0]
        sxy = self.applied_stress[2]
        r = self.r
        
        return sxy*r/2 - 1j*sx*r/2
        
    def phi_1_prime(self, z1):
        """
        Calculates derivative of the stress function. [Eq. A.8, Ref. 1]
        """
        a = self.r
        b = self.r
        mu1 = self.mu1
        mu2 = self.mu2
        alpha = self.alpha()
        beta = self.beta()
        ksi_1, sign_1 = self.ksi_1(z1)
        
        C1 = (beta - mu2*alpha)/(mu1 - mu2)
        eta1 = sign_1*np.sqrt(z1*z1 - a*a - mu1*mu1*b*b)
        kappa1 = 1/(a - 1j*mu1*b)
        
        return -C1/(ksi_1**2)*(1 + z1/eta1)*kappa1
                
    def phi_2_prime(self, z2):
        """
        Calculates derivative of the stress function. [Eq. A.8, Ref. 1]
        """
        a = self.r
        b = self.r
        mu1 = self.mu1
        mu2 = self.mu2
        alpha = self.alpha()
        beta = self.beta()
        ksi_2, sign_2 = self.ksi_2(z2)
        
        C2 = -(beta - mu1*alpha)/(mu1 - mu2)
        eta2 = sign_2*np.sqrt(z2*z2 - a*a - mu2*mu2*b*b)
        kappa2 = 1/(a - 1j*mu2*b)
        
        return -C2/(ksi_2**2)*(1 + z2/eta2)*kappa2
        
    def stress(self, x, y):
        """
        Calculates the laminate average stress at a point (x, y). 
        [Eq. A.9, Ref. 1]
        """
        sx, sy, sxy = super().stress(x, y)
        
        sx_app = self.applied_stress[0]
        sy_app = self.applied_stress[1]
        sxy_app = self.applied_stress[2]
        
        return np.array([sx + sx_app, sy + sy_app, sxy + sxy_app])


class Loaded(Hole):
    
    def __init__(self, r, lam, bearing):
        super().__init__(r, lam)
        self.bearing = np.array(bearing)
        self.px = self.bearing[0]
        self.py = self.bearing[1]
        self.p = np.sum(self.bearing**2)#/(2*r*lam.H) #stress of load?
        self.alpha = np.arctan(self.py/self.px)
        self.A, self.B, self.A_bar, self.B_bar = self.solve_constants()
        
    def solve_constants(self):
        mu1 = self.mu1
        mu2 = self.mu2
        mu1_bar = self.mu1_bar
        mu2_bar = self. mu2_bar
        px = self.bearing[0]
        py = self.bearing[1]
        h = self.lam.H
        a11 = self.lam.a[0,0]
        a22 = self.lam.a[1,1]
        a12 = self.lam.a[0,1]
        a16 = self.lam.a[0,2]
        a26 = self.lam.a[1,2]
        pi = np.pi
        
        mu_mat = np.array([[1, 1, -1, -1],
                           [mu1, mu2, -mu1_bar, -mu2_bar],
                           [mu1**2, mu2**2, -mu1_bar**2, -mu2_bar**2],
                           [1/mu1, 1/mu2, -1/mu1_bar, -1/mu2_bar]])
        
        load_vec = np.array([[py/(2*pi*h*1j)],
                             [-px/(2*pi*h*1j)],
                             [-a16/a11*px/(2*pi*h*1j) - a12/a11*py/(2*pi*h*1j)],
                             [a12/a22*px/(2*pi*h*1j) + a26/a22*py/(2*pi*h*1j)]])
        
        return np.dot(np.linalg.inv(mu_mat), load_vec)
        
    def series_term(self, k, m):
        a = self.r
        p = self.p
        mu1 = self.mu1
        mu2 = self.mu2
        pi = np.pi
        
        if k == 1:
            if m == 2:
                return a*p*1j*(1 + 1j*mu2)/(16*(mu2 - mu1))
            elif m % 2 == 0:
                # m is even                
                return 0
            elif m % 2 == 1:
                # m is odd
                return (-a*p*1j*(-1)**((m - 1)/2)*(2 + 1j*m*mu2)/
                        (pi*m**2*(m**2 - 4)*(mu2 - mu1)))
        elif k == 2:
            if m == 2:
                return -a*p*1j*(1 + 1j*mu1)/(16*(mu2 - mu1))
            elif m % 2 == 0:
                # m is even                
                return 0
            elif m % 2 == 1:
                # m is odd
                return (a*p*1j*(-1)**((m - 1)/2)*(2 + 1j*m*mu1)/
                        (pi*m**2*(m**2 - 4)*(mu2 - mu1)))
        
    def phi_1_prime(self, z1):
        mu1 = self.mu1
        a = self.r
        A = self.A
        A_bar = self.A_bar
        ksi_1, sign_1 = self.ksi_1(z1)
        
        eta_1 = sign_1*np.sqrt(z1*z1 - a*a - a*a*mu1*mu1)
        
        phi_1_series = 0.
        for m in range(1, 45 + 1):
            A1m = self.series_term(1, m)
            phi_1_series += m*A1m/(ksi_1**m)
        
        sai_1 = A + 1j*A_bar        
        
        return (sai_1 - phi_1_series)/eta_1
        
    def phi_2_prime(self, z2):
        mu2 = self.mu2
        a = self.r
        B = self.B
        B_bar = self.B_bar
        ksi_2, sign_2 = self.ksi_2(z2)
        
        eta_2 = sign_2*np.sqrt(z2*z2 - a*a - a*a*mu2*mu2)
        
        phi_2_series = 0.
        for m in range(1, 45 + 1):
            A2m = self.series_term(2, m)
            phi_2_series += m*A2m/(ksi_2**m)
            
        sai_2 = B + 1j*B_bar
            
        return (sai_2 - phi_2_series)/eta_2














































        