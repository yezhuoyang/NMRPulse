import unittest
import unittest
from util import *
from NMRdensity import *
import numpy as np

CNOT_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)

CNOT_rev_matrix = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]], dtype=complex)

P1_matrix = np.matmul(CNOT_rev_matrix,CNOT_matrix)



class P1test(unittest.TestCase):
    def test_P1(self):
        NMRsample = chloroform()
        NMRsample.add_p1_perm_pulse()
        NMRsample.evolve_all_pulse()
        pulse_density = np.abs(NMRsample.get_pulse_unitary())
        distance = hilbert_schmidt_distance(P1_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("P1 Test passed!")