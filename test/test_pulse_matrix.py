from NMRdensity import *
from Pulses import *
from params import *
from util import *


def test_pulse_two_case1():
    NMRsample1 = chloroform()
    NMRsample1.set_density(np.array([[0.5, 0, 0, 0],
                                     [0, 0.3, 0, 0],
                                     [0, 0, -0.3, 0],
                                     [0, 0, 0, -0.5]], dtype=complex))

    NMRsample1.add_pulse(pulseTwo(1, 0.3 * pl90H, wH, 0, 0.4 * pl90C, wC))
    NMRsample1.evolve_all_pulse()
    matrix1 = NMRsample1.get_pulse_unitary()

    print(matrix1)

    NMRsample2 = chloroform()
    NMRsample2.set_density(np.array([[0.5, 0, 0, 0],
                                     [0, 0.3, 0, 0],
                                     [0, 0, -0.3, 0],
                                     [0, 0, 0, -0.5]], dtype=complex))

    NMRsample2.add_pulse(pulseSingle(1, 0.3 * pl90H, wH))
    NMRsample2.add_pulse(pulseSingle(0, 0.4 * pl90C, wC))
    NMRsample2.evolve_all_pulse()
    matrix2 = NMRsample2.get_pulse_unitary()

    print(matrix2)

    assert all_close(matrix1, matrix2)



def test_pulse_two_case2():
    NMRsample1 = chloroform()
    NMRsample1.set_density(np.array([[0.5, 0, 0, 0],
                                     [0, 0.3, 0, 0],
                                     [0, 0, -0.3, 0],
                                     [0, 0, 0, -0.5]], dtype=complex))

    NMRsample1.add_pulse(pulseTwo(1, 0.3 * pl90C, wC, 0, 0.4 * pl90H, wH))
    NMRsample1.evolve_all_pulse()
    matrix1 = NMRsample1.get_pulse_unitary()

    print(matrix1)

    NMRsample2 = chloroform()
    NMRsample2.set_density(np.array([[0.5, 0, 0, 0],
                                     [0, 0.3, 0, 0],
                                     [0, 0, -0.3, 0],
                                     [0, 0, 0, -0.5]], dtype=complex))

    NMRsample2.add_pulse(pulseSingle(1, 0.3 * pl90C, wC))
    NMRsample2.add_pulse(pulseSingle(0, 0.4 * pl90H, wH))
    NMRsample2.evolve_all_pulse()
    matrix2 = NMRsample2.get_pulse_unitary()

    print(matrix2)

    assert all_close(matrix1, matrix2)





if __name__ == "__main__":
    test_pulse_two_case1()
    test_pulse_two_case2()