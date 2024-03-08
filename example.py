import numpy as np

from NMRdensity import *
from params import *
from Pulses import *
import matplotlib.pyplot as plt


def pulse_length_calib_proton():
    pulses_length_list = np.linspace(0, 1, 20)
    integra_list = []
    NMRsample = chloroform()
    for pulse in pulses_length_list:
        NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                        [0, 0.3, 0, 0],
                                        [0, 0, -0.3, 0],
                                        [0, 0, 0, -0.5]], dtype=complex))
        NMRsample.set_pulses([])
        '''
        The first 1/2 pi pulse is added to cancel 
        the sigmax in the measurement operator.
        '''
        NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
        '''
        This is the actual varying Ix pulse we add in the pulse
        length calibration 
        '''
        NMRsample.add_pulse(pulseSingle(0, pulse * pl90H, wH))
        NMRsample.evolve_all_pulse()
        NMRsample.read_proton_time()
        NMRsample.read_proton_spectrum(normalize=False)
        integra_list.append(NMRsample.integral_proton_spectrum_real())

    pulses_length_list = [2 * x * pl90H * 10 ** 6 for x in pulses_length_list]
    plt.scatter(pulses_length_list, integra_list, label="Integral value of proton spectrum")
    plt.axvline(x=pl90H * 10 ** 6, color="red", linestyle="--", label="Measured 90-x pulse for proton")
    plt.xlabel("Pulse length Time/ microsecond")
    plt.ylabel("Integral value")
    plt.legend(fontsize=8)
    plt.savefig("Figure/protoncalib.png")
    plt.show()


def pulse_length_calib_carbon():
    pulses_length_list = np.linspace(0, 1, 20)
    integra_list = []
    NMRsample = chloroform()
    for pulse in pulses_length_list:
        NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                        [0, 0.3, 0, 0],
                                        [0, 0, -0.3, 0],
                                        [0, 0, 0, -0.5]], dtype=complex))
        NMRsample.set_pulses([])
        '''
        The first 1/2 pi pulse is added to cancel 
        the sigmax in the measurement operator.
        '''
        NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
        '''
        This is the actual varying Ix pulse we add in the pulse
        length calibration 
        '''
        NMRsample.add_pulse(pulseSingle(0, pulse * pl90C, wC))
        NMRsample.evolve_all_pulse()
        NMRsample.read_carbon_time()
        NMRsample.read_carbon_spectrum(normalize=False)
        integra_list.append(NMRsample.integral_carbon_spectrum_real())

    pulses_length_list = [2 * x * pl90C * 10 ** 6 for x in pulses_length_list]
    plt.scatter(pulses_length_list, integra_list, label="Integral value of carbon spectrum")
    plt.axvline(x=pl90C * 10 ** 6, color="red", linestyle="--", label="Measured 90-x pulse for carbon")
    plt.xlabel("Pulse length Time/ microsecond")
    plt.ylabel("Integral value")
    plt.legend(fontsize=8)
    plt.savefig("Figure/carboncalib.png")
    plt.show()


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
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, theta * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))

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
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, theta * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))

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
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))
    '''

    NMRsample.set_thermal_equilibrium()

    '''
    Add approximate CNOT pulse sequence
    (pi/2)Ix2---(2Iz1Iz2)---(pi/2)Iy2
    Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
    '''
    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))
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
    NMRsample.read_and_plot("Figure/CNOTapprox")


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
    NMRsample.read_and_plot("Figure/CNOTexact")


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
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))

    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))

    '''
    The pulse of (-2)Iz1Iz2. The angle theta is actually (-\pi/2). However, since we cannot rotate 
    a minus angle, we should plus another (4\pi)
    '''
    NMRsample.add_pulse(delayTime((4 - 0.5) / Jfreq))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Print the unitary of all pulses:
    '''
    print(np.exp(np.pi / 4 * 1j) * NMRsample.get_pulse_unitary())

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_and_plot("Figure/CNOTapprox")



def exact_CNOT_pulse_Hcontrol():
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
    Add pulse sequence for approximate h gate on carbon
    '''
    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))

    '''
    Add pulse sequence for exact CZ gate
    '''
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(delayTime((4 - 0.5) / Jfreq))

    '''
    Add pulse sequence for approximate h gate on carbon
    '''

    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))

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
    NMRsample.read_and_plot("Figure/CNOTExact-Hcontrol")


def exact_CNOT_pulse_Ccontrol():
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
    Add pulse sequence for approximate h gate on carbon
    '''
    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))

    '''
    Add pulse sequence for exact CZ gate
    '''
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(delayTime((4 - 0.5) / Jfreq))

    '''
    Add pulse sequence for approximate h gate on carbon
    '''

    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))

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
    NMRsample.read_and_plot("Figure/CNOTExact-Ccontrol")


def Xgate_proton():
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
    Add a single pulse on proton.
    '''
    NMRsample.add_pulse(pulseSingle(0, pl90H, wH))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()

    matrix = NMRsample.get_pulse_unitary()
    print(matrix)

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_and_plot("Figure/Xgateproton")
    return


def Xgate_carbon():
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
    Add a single pulse on proton.
    '''
    NMRsample.add_pulse(pulseSingle(0, pl90C, wC))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()

    matrix = NMRsample.get_pulse_unitary()
    print(matrix)

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_and_plot("Figure/Xgatecarbon")
    return


def P1_pulse():
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
    Add two approximate CNOT pulse sequence
    (pi/2)Ix2---(2Iz1Iz2)---(pi/2)Iy2
    Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
    '''

    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseTwo(3, 0.5 * pl90C, wC, 0, 0.5 * pl90H, wH))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90H, wH))

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
    NMRsample.read_and_plot("Figure/P1")

    NMRsample.print_pulses()
    return


def P2_pulse():
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
    Add two approximate CNOT pulse sequence
    (pi/2)Ix2---(2Iz1Iz2)---(pi/2)Iy2
    Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
    '''

    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseTwo(3, 0.5 * pl90H, wH, 0, 0.5 * pl90C, wC))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    print(NMRsample.get_pulse_unitary())

    NMRsample.read_and_plot("Figure/P2")
    return


'''
Simulate the result of pseudo pure state
uf: A list of the state input
'''


def pseudo_pure_state(uf):
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()

    NMRsample.set_thermal_equilibrium()

    '''
    Initialize the initial state 00,01,10 or 11
    '''
    if uf[0] == 1:
        NMRsample.add_pulse(pulseSingle(0, pl90H, wH))
    if uf[1] == 1:
        NMRsample.add_pulse(pulseSingle(1, pl90C, wC))

    NMRsample.evolve_all_pulse()

    NMRsample.read_and_plot("Figure/pseudoP0")

    density0 = NMRsample.get_density()

    NMRsample.set_thermal_equilibrium()
    NMRsample.set_pulses([])

    '''
    Add P1 permutation
    '''

    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseTwo(3, 0.5 * pl90C, wC, 0, 0.5 * pl90H, wH))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90H, wH))

    '''
    Initialize the initial state 00,01,10 or 11
    '''
    if uf[0] == 1:
        NMRsample.add_pulse(pulseSingle(0, pl90H, wH))
    if uf[1] == 1:
        NMRsample.add_pulse(pulseSingle(1, pl90C, wC))

    NMRsample.evolve_all_pulse()

    NMRsample.read_and_plot("Figure/pseudoP1")

    # np.set_printoptions(precision=2)
    # print("P1 matrix")
    # print(NMRsample.get_pulse_unitary())

    density1 = NMRsample.get_density()

    NMRsample.set_thermal_equilibrium()
    NMRsample.set_pulses([])

    '''
    Add P2 permutation
    '''

    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseTwo(3, 0.5 * pl90H, wH, 0, 0.5 * pl90C, wC))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))

    '''
    Initialize the initial state 00,01,10 or 11
    '''
    if uf[0] == 1:
        NMRsample.add_pulse(pulseSingle(0, pl90H, wH))
    if uf[1] == 1:
        NMRsample.add_pulse(pulseSingle(1, pl90C, wC))

    NMRsample.evolve_all_pulse()

    # print("P2 matrix")
    # print(NMRsample.get_pulse_unitary())

    NMRsample.read_and_plot("Figure/pseudoP1")

    density2 = NMRsample.get_density()

    pseudo_pure_density = (density0 + density1 + density2) / 3

    new_NMRsample = chloroform()
    new_NMRsample.set_density(pseudo_pure_density)

    new_NMRsample.read_and_plot("Figure/pseudoaverage")

    print(pseudo_pure_density)


def spectrum_only_a():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[1, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0]], dtype=complex))

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
                                        path="Figure/proton-spectruma.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/carbon-spectruma.png")
    return


def spectrum_only_b():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0, 0, 0, 0],
                                    [0, 1, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0]], dtype=complex))

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
                                        path="Figure/proton-spectrumb.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/carbon-spectrumb.png")
    return


def spectrum_only_c():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 1, 0],
                                    [0, 0, 0, 0]], dtype=complex))

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
                                        path="Figure/proton-spectrumc.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/carbon-spectrumc.png")
    return


def spectrum_only_d():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 0],
                                    [0, 0, 0, 1]], dtype=complex))

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
                                        path="Figure/proton-spectrumd.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/carbon-spectrumd.png")
    return


if __name__ == "__main__":
    # pulse_length_change()
    # approx_hadamard_carbon_pulse()
    # IZ_pulse(0.5)
    # SZ_pulse(0.3)
    # approx_CNOT()
    # exact_CZ_pulse()
    # exact_CNOT_pulse_Ccontrol()
    # pulse_length_calib_carbon()
    # pulse_length_calib_proton()
    # P1_pulse()

    # Xgate_proton()
    # Xgate_carbon()

    # approx_CNOT()
    # pseudo_pure_state()
    # spectrum_only_a()
    # spectrum_only_b()
    # spectrum_only_c()
    # spectrum_only_d()

    P1_pulse()
