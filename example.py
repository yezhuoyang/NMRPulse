import numpy as np

from NMRdensity import *
from params import *
from Pulses import *
import matplotlib.pyplot as plt

'''
Rz(theta)=Rx(-\pi/2)Ry(theta)Rx(\pi/2)
'''


def IZ_pulse(theta):
    NMRsample = chloroform()
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))

    Rz_theta = Rz_matrix(theta * np.pi)

    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.add_pulse(pulseSingle(2, 1/2*pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, theta * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1/2*pl90H, wH))


    NMRsample.evolve_all_pulse()
    matrix = NMRsample.get_pulse_unitary()
    print("Ix(-\pi/2)Iy(theta)Ix(\pi/2):")
    print(matrix)

    true_matrix = np.kron(Rz_theta, I_matrix())
    print("True value")
    print(true_matrix)
    return


def SZ_pulse(theta):
    NMRsample = chloroform()
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))

    Rz_theta = Rz_matrix(theta * np.pi)

    '''
    Clear the pulses
    '''
    NMRsample.add_pulse(pulseSingle(2, 1/2*pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, theta * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1/2*pl90C, wC))


    NMRsample.evolve_all_pulse()
    matrix = NMRsample.get_pulse_unitary()
    print("Sx(-\pi/2)Sy(theta)Sx(\pi/2):")
    print(matrix)

    true_matrix = np.kron(I_matrix(), Rz_theta)
    print("True value")
    print(true_matrix)
    return


def approx_hadamard_carbon_pulse():
    NMRsample = chloroform()
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))

    '''
    Clear the pulses
    '''
    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))

    NMRsample.evolve_all_pulse()
    matrix = NMRsample.get_pulse_unitary()
    print("(pi/4)Sy--(pi)Sx--(-pi/4)Sy")
    print(matrix)
    return


def approx_hadamard_proton_pulse():
    NMRsample = chloroform()
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))

    '''
    Clear the pulses
    '''
    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))

    NMRsample.evolve_all_pulse()
    matrix = NMRsample.get_pulse_unitary()
    print("(pi/4)Iy--(pi)Ix--(-pi/4)Iy")
    print(matrix)
    return






def pulse_length_change():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))
    '''
    Add a single pulse on proton.
    '''
    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_fid_real(maxtime=0.1, store=True, path="Figure/45pulsesFID.png")
    '''
    Read the data signal in the frequency
    '''
    NMRsample.read_proton_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/45pulsespec.png")


def approx_CNOT():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))

    '''
    Add approximate CNOT pulse sequence
    (pi/2)Ix2---(2Iz1Iz2)---(pi/2)Iy2
    Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
    '''
    NMRsample.add_pulse(pulseSingle(0, 0.5*pl90C, wC))
    NMRsample.add_pulse(delayTime(1 / Jfreq))
    NMRsample.add_pulse(pulseSingle(1, 0.5*pl90C, wC))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Print the unitary of all pulses:
    '''
    print(NMRsample.get_pulse_unitary())

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    NMRsample.read_carbon_time()
    '''
    Read the spectrum
    '''
    NMRsample.read_proton_spectrum()
    NMRsample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/CNOTapproxproton.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/CNOTapproxcarbon.png")


def exact_CNOT():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))
    CNOTmatrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
    '''
    Directly evolve the density matrix by CNOT matrix
    '''
    NMRsample.evolve_density(CNOTmatrix)
    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    NMRsample.read_carbon_time()
    '''
    Read the spectrum
    '''
    NMRsample.read_proton_spectrum()
    NMRsample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/CNOTExactproton.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/CNOTExactcarbon.png")


def exact_CZ_pulse():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))

    '''
    Add pulse sequence for exact CZ gate
    (pi/2)Iz1---(pi/2)Iz2---(-2)Iz1Iz2
    Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
    '''
    NMRsample.add_pulse(pulseSingle(2, 1/2*pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, 1/2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1/2*pl90H, wH))

    NMRsample.add_pulse(pulseSingle(2, 1/2*pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, 1/2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1/2*pl90C, wC))


    NMRsample.add_pulse(delayTime(1 / Jfreq))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Print the unitary of all pulses:
    '''
    print(NMRsample.get_pulse_unitary())

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    NMRsample.read_carbon_time()
    '''
    Read the spectrum
    '''
    NMRsample.read_proton_spectrum()
    NMRsample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/CNOTapproxproton.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/CNOTapproxcarbon.png")


def exact_CNOT_pulse():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))


def P1_pulse():
    return


def P2_pulse():
    return


if __name__ == "__main__":
    # pulse_length_change()
    #approx_hadamard_carbon_pulse()
    #IZ_pulse(0.5)
    #SZ_pulse(0.3)
    #approx_CNOT()
    exact_CZ_pulse()