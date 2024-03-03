import qiskit
import functools
from qiskit_aer import AerSimulator
from NMRdensity import *
from Pulses import *




class NMRalgorithm:
    def __init__(self,
                 *params):
        pass

    def set_function(self,
                     *params):
        raise NotImplementedError

    def construct_pulse(self):
        raise NotImplementedError


    def construct_circuit(self):
        raise NotImplementedError


    def calculate_result(self):
        raise NotImplementedError


class Djalgorithm:

    def __init__(self):
        pass

    def construct_pulse(self):
        pass

    def calculate_result(self):
        pass


    def construct_circuit(self):
        pass


    def set_function(self,
                     *params):
        pass


class Grover:

    def __init__(self):
        pass

    def construct_pulse(self):
        pass

    def calculate_result(self):
        pass


    def construct_circuit(self):
        pass


    def set_function(self,
                     *params):
        pass
