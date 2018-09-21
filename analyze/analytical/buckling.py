# -*- coding: utf-8 -*-
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

from math import pi
import numpy as np
from scipy.linalg import eigvals


def absBottom(i, iterable):
    """
    Returns the absolute bottom i values from sorted iterable while maintaining
    sign.
    """
    return sorted(iterable,key=abs)[:i]


class CurvedPlate:
    """
    This class defines an anisotropic curved plate buckling solution.
    Rerference:
         Ward, S.W., McClesky, S.F., "Buckling of Laminated Shells Including
         Transverse Shear Flexibility" 28th SDM Cont. April 8, 1987.
         AIAA-87-0728-CP
    """
    #class variables
    NUMEIGS = 3 #int
    I = J = 8 #int
    MAX_FLAT_ASPECT_RATIO = 6.0 #float
    MAX_R10_to_FLAT_ASPECT_RATIO = 6.0 #float
    MAX_R4_to_R10_ASPECT_RATIO = 3.0 #float
    MAX_LESS_THAN_R4_ASPECT_RATIO = 2.0 #float

    def __init__(self, a=0., b=0., r=1.e8, lam=None, bc='simple'):
        self.a = a
        self.b = b
        self.r = r
        self.bc = bc
        self.__negShear = False
        self.__lam = lam
        self.__checkLongPlate()
        self.__T = np.array([])
        self.__generateStiffnessMatrix()

    def setLaminate(self, lam):
        self.__lam = lam
        self.__generateStiffnessMatrix()

    def getLaminate(self):
        return self.__lam

    def __checkLongPlate(self):
        a = self.a
        b = self.b
        r = self.r
        lam = self.getLaminate()

        flat = self.MAX_FLAT_ASPECT_RATIO
        r10 = self.MAX_R10_to_FLAT_ASPECT_RATIO
        r4to10 = self.MAX_R4_to_R10_ASPECT_RATIO
        ltR4 = self.MAX_LESS_THAN_R4_ASPECT_RATIO

        # Check for long plate
        if r/b >= 50.0:
            # assume flat plate
            if not 1/flat < a/b < flat:
                # confirm long plate assumption applies and modify dimensions
                if a > b and a > 3*b*(lam.D[0,0]/lam.D[1,1])**(1/4):
                    self.a = b*flat
                elif b > a and b > 3*a*(lam.D[1,1]/lam.D[0,0])**(1/4):
                    self.b = a*flat
                else:
                    raise ValueError("Long plate assumption does not apply! Try "
                                     "removing anisotropy from the laminate.")
        elif 10.0 <= r/b < 50.0:
            # shallow curved plate
            if not 1/r10 < a/b < r10 and a > b:
                self.a = b*r10
            elif not 1/r10 < a/b < r10 and b > a:
                self.b = a*r10
        elif 4.0 < r/b < 10.0:
            # medium curved plate
            if not 1/r4to10 < a/b < r4to10 and a > b:
                self.a = b*r4to10
            elif not 1/r4to10 < a/b < r4to10 and b > a:
                self.b = a*r4to10
        elif r/b <= 4.0:
            # deep curved plate
            if not 1/ltR4 < a/b < ltR4 and a > b:
                self.a = b*ltR4
            elif not 1/ltR4 < a/b < ltR4 and b > a:
                self.b = a*ltR4

    def __L1L4L6(self, m, n, p, q, a, b):
        if m==p and n==q:
            L1 = L4 = L6 = a*b/4.0
        else:
            L1 = L4 = L6 = 0
        return L1, L4, L6

    def __L2L3L5(self, m, n, p, q, a, b):
        if (n+q)%2==0 or (m+p)%2==0:
            L2 = L3 = L5 = 0
        else:
            denom = pi**2*(p**2-m**2)*(q**2-n**2)
            L2 = -4.0*a*b*m*q/denom
            L3 = -4.0*a*b*n*p/denom
            L5 = 4.0*a*b*p*q/denom
        return L2, L3, L5

    def __generateStiffnessMatrix(self):

        lam = self.getLaminate()

        if not lam:
            print("No buckling laminate defined!")
            return False

        A11, A12, A16 = lam.A[0,0], lam.A[0,1], lam.A[0,2]
        A22, A26 = lam.A[1,1], lam.A[1,2]
        A66 = lam.A[2,2]

        B11, B12, B16 = lam.B[0,0], lam.B[0,1], lam.B[0,2]
        B22, B26 = lam.B[1,1], lam.B[1,2]
        B66 = lam.B[2,2]

        D11, D12, D16 = lam.D[0,0], lam.D[0,1], lam.D[0,2]
        D22, D26 = lam.D[1,1], lam.D[1,2]
        D66 = lam.D[2,2]

        if self.__negShear:
            A16 *= -1
            A26 *= -1
            B16 *= -1
            B26 *= -1
            D16 *= -1
            D26 *= -1

        a, b, r = self.a, self.b, self.r
        I, J = self.I, self.J

        Tmnpq = np.zeros((3,3,(I*J),(I*J)))

        for m in range(1,I+1):
            for n in range(1,J+1):
                for p in range(1,I+1):
                    for q in range(1,J+1):

                        k = (p-1)*J + q - 1
                        l = (m-1)*J + n - 1

                        M = m*pi/a
                        N = n*pi/b

                        L1, L4, L6 = self.__L1L4L6(m,n,p,q,a,b)
                        L2, L3, L5 = self.__L2L3L5(m,n,p,q,a,b)

                        if (n+q) % 2 == 0:
                            lam1 = 0
                        else:
                            lam1 = 2*b*q*(1-(-1)**(m+p))/(q**2 - n**2)/pi

                        if (m+p) % 2 == 0:
                            lam2 = 0
                        else:
                            lam2 = 2*a*p*(1-(-1)**(n+q))/(p**2 - m**2)/pi

                        # w
                        Tmnpq[0,0,k,l] += ((2*B12*M**2*L6 - 4*B26*M*N*L5 + 2*B22*N**2*L6)/r
                                          + L6*(D11*M**4 + D22*N**4
                                          + (2*D12 + 4*D66)*M**2*N**2 + A22/r**2)
                                          - 4*L5*M*N*(D16*M**2 + D26*N**2)
                                          - 2*m*n*pi**3*(D16*lam1*p + D26*lam2*q)/(a**2*b))

                        Tmnpq[0,1,k,l] += (L6*M*(B11*M**2 + (B12 + 2*B66)*N**2 + A12/r)
                                          - L5*N*(3*B16*M**2 + B26*N**2 + A26/r)
                                          + n*pi**2*(B16*lam1*p + B26*lam2*q)/(a*b))

                        Tmnpq[0,2,k,l] += (L6*N*(B22*N**2 + (B12 + 2*B66)*M**2 + A22/r)
                                          - L5*M*(B16*M**2 + 3*B26*N**2 + A26/r)
                                          + m*pi**2*(B16*lam1*p + B26*lam2*q)/a**2)

                        # u
                        Tmnpq[1,0,k,l] += ((A12*M*L1 + A26*N*L2)/r + L1*M*(B11*M**2
                                          + (B12 + 2*B66)*N**2)
                                          + L2*N*(3*B16*M**2 + B26*N**2))

                        Tmnpq[1,1,k,l] += (L1*(A11*M**2 + A66*N**2) + 2*L2*M*N*A16)

                        Tmnpq[1,2,k,l] += (L1*M*N*(A12 + A66) + L2*(A26*N**2 + A16*M**2))

                        # v
                        Tmnpq[2,0,k,l] += ((A22*N*L4 + A26*M*L3)/r + L4*N*(B22*N**2
                                          + (B12 + 2*B66)*M**2)
                                          + L3*M*(B16*M**2 + 3*B26*N**2))

                        Tmnpq[2,1,k,l] += (L4*M*N*(A12 + A66) + L3*(A26*N**2 + A16*M**2))

                        Tmnpq[2,2,k,l] += (L4*(A22*N**2 + A66*M**2) + 2*L3*M*N*A26)

        T11 = Tmnpq[0,0]

        T1 = Tmnpq[0,1:].reshape((2*(I*J),(I*J))).T

        T2 = Tmnpq[1:,0].reshape((2*(I*J),(I*J)))

        T3 = np.zeros((2*(I*J),2*(I*J)))
        T3[:I*J,:I*J] += Tmnpq[1,1]
        T3[:I*J,I*J:] += Tmnpq[1,2]
        T3[I*J:,:I*J] += Tmnpq[2,1]
        T3[I*J:,I*J:] += Tmnpq[2,2]

        self.__T = T11 - T1.dot(np.linalg.pinv(T3).dot(T2))

        return True

    def run(self, Nx=0., Ny=0., Nxy=0.):

        if len(self.__T) == 0:
            return [None]

        a, b = self.a, self.b
        I, J = self.I, self.J

        if Nxy < 0 and not self.__negShear:
            self.__negShear = True
            self.__generateStiffnessMatrix()
        elif Nxy > 0 and self.__negShear:
            self.__negShear = False
            self.__generateStiffnessMatrix()
        else:
            # do nothing
            pass
        
        Nxy = abs(Nxy)
        
        N1 = np.zeros(((I*J),(I*J)))

        for m in range(1,I+1):
            for n in range(1,J+1):
                for p in range(1,I+1):
                    for q in range(1,J+1):

                        k = (p-1)*J + q - 1
                        l = (m-1)*J + n - 1

                        M = m*pi/a
                        N = n*pi/b

                        L6 = self.__L1L4L6(m,n,p,q,a,b)[-1]
                        L5 = self.__L2L3L5(m,n,p,q,a,b)[-1]

                        N1[k,l] = (2*Nxy*M*N*L5 - L6*(Nx*M**2 + Ny*N**2))

        return absBottom(self.NUMEIGS, np.real(eigvals(self.__T,N1)))