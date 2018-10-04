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
from strengl.geometry.primitives import Tensor


class Stress(Tensor):

    def __init__(self, stresses):
        """

        :param stresses: <np.array> (6x1) or (3x3),
            Voight notation: [s11, s22, s33, s23, s13, s12]
            or
            Tensor form: [[s11, s12, s13],
                          [s21, s22, s23],
                          [s31, s32, s33]]
        """
        if self.tensor.shape == (6,):
            s11, s22, s33, s23, s13, s12 = stresses
            tensor = np.array([[s11, s12, s13],
                               [s12, s22, s23],
                               [s13, s23, s33]])
        else:
            tensor = np.array(stresses)
        super().__init__(tensor)

