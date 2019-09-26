# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 07:56:23 2018

@author: Ben
"""
#import abc
import numpy as np


def fourier_series_coefficients(f, T, N, return_complex=True):
    """
    Ref: https://stackoverflow.com/questions/4258106/how-to-calculate-a-fourier-series-in-numpy

    Calculates the first 2*N+1 Fourier series coeff. of a periodic function.

    Given a periodic, function f(t) with period T, this function returns the
    complex coefficients {c0,c1,c2,...}
    such that:

    f(t) ~= sum_{k=-N}^{N} c_k * exp(i*2*pi*k*t/T)

    where we define c_{-n} = complex_conjugate(c_{n})

    Refer to wikipedia for the relation between the real-valued and complex
    valued coeffs at http://en.wikipedia.org/wiki/Fourier_series.

    Parameters
    ----------
    f : the periodic function, a callable like f(t)
    T : the period of the function f, so that f(0)==f(T)
    N : the function will return the first N + 1 Fourier coeff.

    Returns
    -------
    c : numpy 1-dimensional complex-valued array of size N+1

    """
    # From Shanon theoreom we must use a sampling freq. larger than the maximum
    # frequency you want to catch in the signal.
    f_sample = 2 * N

    t = np.linspace(-T/2, T/2, f_sample, endpoint=False)

    y = np.fft.rfft(f(t)) / t.size

    # multiply odd terms by -1 to match SageMath
    y[1::2] *= -1

    if return_complex:
        return y
    else:
        y *= 2
        return y[0].real, y[1:].real, -y[1:].imag


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
    def __init__(self, radius, thickness, a_matrix):
        """
        Class contructor.
        
        :param radius: hole radius
        :param thickness: laminate thickness
        :param a_matrix: <np.array> inverser laminate A-matrix
        """
        self.r = radius
        self.a = a_matrix
        self.h = thickness
        self.mu1, self.mu2, self.mu1_bar, self.mu2_bar = self.roots()

    def roots(self):
        """
        Finds the roots to the characteristic equation [Eq. A.2, Ref. 1]
        """
        a11 = self.a[0, 0]
        a12 = self.a[0, 1]
        a16 = self.a[0, 2]
        a22 = self.a[1, 1]
        a26 = self.a[1, 2]
        a66 = self.a[2, 2]

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

    def __init__(self, radius, thickness, a_matrix, bypass):
        """
        Class constructor.
        
        :param radius: hole radius
        :param thickness: laminate thickness
        :param a_matrix: inverse laminate A-matrix
        :param bypass: [Nx, Ny, Nxy] force / unit length
        """
        super().__init__(radius, thickness, a_matrix)
        self.applied_stress = np.array(bypass)/self.h

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
    """
    Class for defining an loaded hole in an anisotropic homogeneous plate
    with bearing force (Px) applied [Fig. 10, Ref. 4].
    """
    def __init__(self, radius, thickness, a_matrix, bearing):
        """
        :param radius: Hole radius
        :param thickness: laminate thickness
        :param a_matrix: inverse laminate A-matrix
        :param bearing: bearing force (x-dir only)
        """
        super().__init__(radius, thickness, a_matrix)
        self.p = bearing
        self.A, self.A_bar, self.B, self.B_bar = self.solve_constants()

    def alphas(self, N):

        h = self.h
        p = self.p

        def brg_load_x_component(thetas):
            """
            Cosine load distribution [Fig. 10, Ref. 4]
            :param thetas: <np.array> angles
            """
            new_array = np.zeros(len(thetas))
            for i, theta in enumerate(thetas):
                # assume bearing load is only applied in x-dir
                if -np.pi/2 <= theta <= np.pi/2:
                    new_array[i] = np.cos(theta)**2
            return new_array

        coefficients = fourier_series_coefficients(brg_load_x_component, 2*np.pi, N, return_complex=False)[1]

        return -p/(2*np.pi*h)*coefficients

    def betas(self, N):

        h = self.h
        p = self.p

        def brg_load_y_component(thetas):
            """
            Cosine load distribution [Fig. 10, Ref. 4]
            :param thetas: <np.array> angles
            """
            new_array = np.zeros(len(thetas))
            for i, theta in enumerate(thetas):
                # assume bearing load is only applied in x-dir
                if -np.pi/2 <= theta <= np.pi/2:
                    new_array[i] = np.cos(theta)*np.sin(theta)
            return new_array

        m = np.arange(1, N + 1)
        coefficients = fourier_series_coefficients(brg_load_y_component, 2*np.pi, N, return_complex=False)[2]

        return 2*p/(np.pi*h*m**2*1j)*coefficients

    def solve_constants(self):
        """
        Eq. 37.5 [ref. 2] expanding complex terms and resolving for A, A_bar, B and B_bar (setting Px equal to zero)
        """
        R1, R2 = np.real(self.mu1), np.imag(self.mu1)
        R3, R4 = np.real(self.mu2), np.imag(self.mu2)
        p = self.p
        h = self.h
        a11 = self.a[0, 0]
        a12 = self.a[0, 1]
        a22 = self.a[1, 1]
        a16 = self.a[0, 2]
        pi = np.pi

        mu_mat = np.array([[0., 1., 0., 1.],
                           [R2, R1, R4, R3],
                           [2*R1*R2, R1**2 - R2**2, 2*R3*R4, R3**2 - R4**2],
                           [R2/(R1**2 + R2**2), -R1/(R1**2 + R2**2), R4/(R3**2 + R4**2), -R3/(R3**2 + R4**2)]])

        load_vec = p/(4.*pi*h)*np.array([[0.],
                                        [1.],
                                        [a16/a11],
                                        [a12/a22]])

        return np.dot(np.linalg.inv(mu_mat), load_vec)

    def phi_1_prime(self, z1):
        """
        Calculates derivative of the stress function. [Eq. 37.6, Ref. 2]
        """
        mu1 = self.mu1
        mu2 = self.mu2
        a = self.r
        b = self.r
        A = self.A+1j*self.A_bar
        ksi_1, sign_1 = self.ksi_1(z1)

        eta_1 = sign_1*np.sqrt(z1*z1 - a*a - b*b*mu1*mu1)

        N = 45
        m = np.arange(1, N + 1)
        alphas = self.alphas(N)
        betas = self.betas(N)

        return 1/eta_1*(A - np.sum(m*(betas - mu2*alphas)/(mu1 - mu2)/ksi_1**m))

    def phi_2_prime(self, z2):
        """
        Calculates derivative of the stress function. [Eq. 37.6, Ref. 2]
        """
        mu1 = self.mu1
        mu2 = self.mu2
        a = self.r
        b = self.r
        B = self.B+1j*self.B_bar
        ksi_2, sign_2 = self.ksi_2(z2)

        eta_2 = sign_2*np.sqrt(z2*z2 - a*a - b*b*mu2*mu2)

        N = 45
        m = np.arange(1, N + 1)
        alphas = self.alphas(N)
        betas = self.betas(N)

        return 1/eta_2*(B + np.sum(m*(betas - mu1*alphas)/(mu1 - mu2)/ksi_2**m))














































