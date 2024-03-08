import unittest
import unittest
from util import *
from NMRdensity import *
import numpy as np


Imatrix = np.array([[1, 0], [0, 1]], dtype=complex)
Hmatrix = 1/np.sqrt(2)*np.array([[1, 1], [1, -1]], dtype=complex)
H2_matrix = np.kron(Imatrix, Hmatrix)
H1_matrix = np.kron(Hmatrix, Imatrix)



class Hadamardtest(unittest.TestCase):
    def test_H1(self):
        NMRsample = chloroform()
        NMRsample.add_H_gate_first_pulse()
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(H1_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("H1 Test passed!")
    def test_H2(self):
        NMRsample = chloroform()
        NMRsample.add_H_gate_second_pulse()
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(H2_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("H2 Test passed!")