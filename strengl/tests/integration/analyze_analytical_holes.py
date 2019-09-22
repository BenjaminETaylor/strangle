import unittest
import numpy as np


class LoadedHoleTests(unittest.TestCase):
    def test_stress_at_x_equals_r(self):
        a_matrix = np.array([[2.65646e-6, -8.91007e-7, 0.], [-8.91007e-7, 2.65646e-6, 0.], [0., 0., 7.09494e-6]])
        r = 0.125
        h = 0.058
        p = 100.
