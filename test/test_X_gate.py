import unittest
from Simulator.util import *
from Simulator.NMRdensity import *
import numpy as np

Imatrix = np.array([[1, 0], [0, 1]], dtype=complex)
Xmatrix = np.array([[0, 1], [1, 0]], dtype=complex)
X2_matrix = np.kron(Imatrix, Xmatrix)
X1_matrix = np.kron(Xmatrix, Imatrix)


class Xtest(unittest.TestCase):
    def test_X1(self):
        NMRsample = chloroform()
        NMRsample.add_X_gate_first_pulse()
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(X1_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("X1 Test passed!")

    def test_X2(self):
        NMRsample = chloroform()
        NMRsample.add_X_gate_second_pulse()
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(X2_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("X2 Test passed!")
