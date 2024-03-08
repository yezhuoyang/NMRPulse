import unittest
import unittest
from util import *
from NMRdensity import *
import numpy as np

CZ_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]], dtype=complex)

CZ_rev_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 1], [0, 0, 1, 0], [0, 0, 0, -1]], dtype=complex)





class CZtest(unittest.TestCase):
    def test_Hcontrol(self):
        NMRsample = chloroform()
        NMRsample.add_CZ_pulse(Hcontrol=True, approximate=False)
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(CZ_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("CZ(H,C) Test passed!")
    def test_Ccontrol(self):
        NMRsample = chloroform()
        NMRsample.add_CZ_pulse(Hcontrol=False, approximate=False)
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(CZ_rev_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("CZ(C,H) Test passed!")