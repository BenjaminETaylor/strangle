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


class Strain(Tensor):

    def __init__(self, strain_tensor):
        """
        :param strain_tensor: <np.array> (3x3)
            Tensor form: [[e11, e12, e13],
                          [e21, e22, e23],
                          [e31, e32, e33]]
        """
        super().__init__(strain_tensor)

    @property
    def engineering(self):
        """
        Converts strain tensor to engineering strain
        :return: <np.array> (6x1)
            [e11, e22, e33, 2*e23, 2*e13, 2*e12]
        """
        e = self.tensor
        e11, e22, e33, e23, e13, e12 = e[0, 0], e[1, 1], e[2, 2], e[1, 2], e[0, 2], e[0, 1]
        return np.array([e11, e22, e33, 2*e23, 2*e13, 2*e12]).transpose()


