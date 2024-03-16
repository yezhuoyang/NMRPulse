import unittest
from Simulator.util import *
from Simulator.NMRdensity import *
import numpy as np

CNOT_matrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)

CNOT_rev_matrix = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]], dtype=complex)


class CNOTtest(unittest.TestCase):
    def test_Hcontrol(self):
        NMRsample = chloroform()
        NMRsample.add_CNOT_pulse(Hcontrol=True, approximate=False)
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(CNOT_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("CNOT(H,C) Test passed!")
    def test_Ccontrol(self):
        NMRsample = chloroform()
        NMRsample.add_CNOT_pulse(Hcontrol=False, approximate=False)
        NMRsample.evolve_all_pulse()
        pulse_density = NMRsample.get_pulse_unitary()
        distance = hilbert_schmidt_distance(CNOT_rev_matrix, pulse_density)
        self.assertAlmostEqual(distance, 0)
        print("CNOT(C,H) Test passed!")


