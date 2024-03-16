import unittest
from Simulator.util import *
from Simulator.NMRdensity import *
import numpy as np

CNOT_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)

CNOT_rev_matrix = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]], dtype=complex)

P2_matrix = np.matmul(CNOT_matrix,CNOT_rev_matrix)



class P2test(unittest.TestCase):
    def test_P2(self):
        NMRsample = chloroform()
        NMRsample.add_p2_perm_pulse()
        NMRsample.evolve_all_pulse()
        pulse_density = np.abs(NMRsample.get_pulse_unitary())
        distance = hilbert_schmidt_distance(P2_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("P2 Test passed!")