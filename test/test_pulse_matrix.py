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


def test_pulse_two_case3():
    NMRsample1 = chloroform()
    NMRsample1.set_density(np.array([[0.5, 0, 0, 0],
                                     [0, 0.4, 0, 0],
                                     [0, 0, 0.05, 0],
                                     [0, 0, 0, 0.05]], dtype=complex))

    NMRsample1.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
    NMRsample1.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample1.add_pulse(pulseTwo(1, 0.5 * pl90H, wH, 0, 0.5 * pl90C, wC))
    NMRsample1.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample1.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
    NMRsample1.evolve_all_pulse()
    matrix1 = NMRsample1.get_pulse_unitary()

    NMRsample1.read_proton_time()
    NMRsample1.read_proton_spectrum()
    NMRsample1.show_proton_spectrum_real(-5, 15)
    NMRsample1.show_proton_fid_real(maxtime=0.01, miny=-0.000001, maxy=0.000001)

    print(NMRsample1._proton_time_domain)
    print(matrix1)

    NMRsample2 = chloroform()
    NMRsample2.set_density(np.array([[0.5, 0, 0, 0],
                                     [0, 0.3, 0, 0],
                                     [0, 0, -0.3, 0],
                                     [0, 0, 0, -0.5]], dtype=complex))

    NMRsample2.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
    NMRsample2.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample2.add_pulse(pulseSingle(1, 0.5 * pl90H, wH))
    NMRsample2.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
    NMRsample2.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample2.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))

    NMRsample2.evolve_all_pulse()
    matrix2 = NMRsample2.get_pulse_unitary()

    NMRsample2.read_proton_time()
    NMRsample2.read_proton_spectrum()
    NMRsample2.show_proton_spectrum_real(-5, 15)
    NMRsample2.show_proton_fid_real(maxtime=0.01, miny=-0.000001, maxy=0.000001)

    print(NMRsample2._proton_time_domain)

    print(matrix2)

    print(norm_2_distance(matrix1, matrix2))

    assert all_close(matrix1, matrix2)

    assert all_close(NMRsample1.get_density(), NMRsample2.get_density())


if __name__ == "__main__":
    # test_pulse_two_case1()
    # test_pulse_two_case2()
    test_pulse_two_case3()
