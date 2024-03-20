from Simulator.NMRdensity import *
from Simulator.Pulses import *
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

'''
Rz(-theta)=Rx(-\pi/2)Ry(-theta)Rx(\pi/2)
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


def new_CNOT_pulse():
    NMRsample = chloroform()
    NMRsample.set_thermal_equilibrium()
    # Rz_theta = Rz_matrix(theta * np.pi)

    '''
    Clear the pulses
    '''
    NMRsample.add_CZ_pulse(exact=True)
    #NMRsample.add_CNOT_pulse(exact=True)

    NMRsample.evolve_all_pulse()
    matrix = NMRsample.get_pulse_unitary()
    print("New CZ")
    print(matrix)

    return


def approx_hadamard_carbon_pulse():
    NMRsample = chloroform()
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))

    '''
    Add the pulses
    '''
    NMRsample.add_H_gate_first_pulse(approximate=True)

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
    Add the pulses
    '''
    NMRsample.add_H_gate_second_pulse(approximate=True)

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
    NMRsample.add_CNOT_pulse(Hcontrol=True, approximate=True)
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

    NMRsample.add_CZ_pulse(Hcontrol=True, approximate=False)
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
    NMRsample.add_CNOT_pulse(approximate=False, Hcontrol=True)

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

    NMRsample.add_CNOT_pulse(approximate=False, Hcontrol=False)

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


def pseudo_pure_state(uf, add_CNOT=True, add_approx_CNOT=False, generateprogram=False):
    assert not (add_CNOT and add_approx_CNOT)

    ufstring = str(uf[0]) + str(uf[1])

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

    if add_CNOT:
        NMRsample.add_CNOT_pulse(approximate=False, Hcontrol=True)
    elif add_approx_CNOT:
        NMRsample.add_CNOT_pulse(approximate=True, Hcontrol=True)

    NMRsample.evolve_all_pulse()

    if add_CNOT:
        NMRsample.read_and_plot("Figure/Pseudo/{}/CNOTpseudoP0f{}".format(ufstring, ufstring))
    elif add_approx_CNOT:
        NMRsample.read_and_plot("Figure/Pseudo/{}/approxCNOTpseudoP0f{}".format(ufstring, ufstring))
    else:
        NMRsample.read_and_plot("Figure/Pseudo/{}/pseudoP0f{}".format(ufstring, ufstring))

    if generateprogram:
        NMRsample.insert_delay()
        pulse_program = NMRsample.print_pulses()
        if add_CNOT:
            store_string_to_file(pulse_program, "data/CNOT/{}/Code_P0f{}.txt".format(ufstring, ufstring))
        elif add_approx_CNOT:
            store_string_to_file(pulse_program,
                                 "data/ApproxCNOT/{}/Code_P0f{}.txt".format(ufstring, ufstring))
        else:
            store_string_to_file(pulse_program,
                                 "data/Pseudo/{}/Code_P0f{}.txt".format(ufstring, ufstring))

    density0 = NMRsample.get_density()

    NMRsample.set_thermal_equilibrium()
    NMRsample.set_pulses([])

    '''
    Add P1 permutation
    '''

    NMRsample.add_p1_perm_pulse()

    '''
    Initialize the initial state 00,01,10 or 11
    '''
    if uf[0] == 1:
        NMRsample.add_pulse(pulseSingle(0, pl90H, wH))
    if uf[1] == 1:
        NMRsample.add_pulse(pulseSingle(1, pl90C, wC))

    if add_CNOT:
        NMRsample.add_CNOT_pulse(approximate=False, Hcontrol=True)
    if add_approx_CNOT:
        NMRsample.add_CNOT_pulse(approximate=True, Hcontrol=True)

    NMRsample.evolve_all_pulse()

    if add_CNOT:
        NMRsample.read_and_plot("Figure/Pseudo/{}/CNOTpseudoP1f{}".format(ufstring, ufstring))
    elif add_approx_CNOT:
        NMRsample.read_and_plot("Figure/Pseudo/{}/approxCNOTpseudoP1f{}".format(ufstring, ufstring))
    else:
        NMRsample.read_and_plot("Figure/Pseudo/{}/pseudoP1f{}".format(ufstring, ufstring))

    if generateprogram:
        NMRsample.insert_delay()
        pulse_program = NMRsample.print_pulses()
        if add_CNOT:
            store_string_to_file(pulse_program, "data/CNOT/{}/Code_P1f{}.txt".format(ufstring, ufstring))
        elif add_approx_CNOT:
            store_string_to_file(pulse_program,
                                 "data/ApproxCNOT/{}/Code_P1f{}.txt".format(ufstring, ufstring))
        else:
            store_string_to_file(pulse_program,
                                 "data/Pseudo/{}/Code_P1f{}.txt".format(ufstring, ufstring))

    # np.set_printoptions(precision=2)
    # print("P1 matrix")
    # print(NMRsample.get_pulse_unitary())

    density1 = NMRsample.get_density()

    NMRsample.set_thermal_equilibrium()
    NMRsample.set_pulses([])

    '''
    Add P2 permutation
    '''

    NMRsample.add_p2_perm_pulse()

    '''
    Initialize the initial state 00,01,10 or 11
    '''
    if uf[0] == 1:
        NMRsample.add_pulse(pulseSingle(0, pl90H, wH))
    if uf[1] == 1:
        NMRsample.add_pulse(pulseSingle(1, pl90C, wC))

    if add_CNOT:
        NMRsample.add_CNOT_pulse(approximate=False, Hcontrol=True)
    if add_approx_CNOT:
        NMRsample.add_CNOT_pulse(approximate=True, Hcontrol=True)

    NMRsample.evolve_all_pulse()

    # print("P2 matrix")
    # print(NMRsample.get_pulse_unitary())

    if add_CNOT:
        NMRsample.read_and_plot("Figure/Pseudo/{}/CNOTpseudoP2f{}".format(ufstring, ufstring))
    elif add_approx_CNOT:
        NMRsample.read_and_plot("Figure/Pseudo/{}/approxCNOTpseudoP2f{}".format(ufstring, ufstring))
    else:
        NMRsample.read_and_plot("Figure/Pseudo/{}/pseudoP2f{}".format(ufstring, ufstring))

    if generateprogram:
        NMRsample.insert_delay()
        pulse_program = NMRsample.print_pulses()
        if add_CNOT:
            store_string_to_file(pulse_program, "data/CNOT/{}/Code_P2f{}.txt".format(ufstring, ufstring))
        elif add_approx_CNOT:
            store_string_to_file(pulse_program,
                                 "data/ApproxCNOT/{}/Code_P2f{}.txt".format(ufstring, ufstring))
        else:
            store_string_to_file(pulse_program,
                                 "data/Pseudo/{}/Code_P2f{}.txt".format(ufstring, ufstring))

    density2 = NMRsample.get_density()

    pseudo_pure_density = (density0 + density1 + density2) / 3

    new_NMRsample = chloroform()
    new_NMRsample.set_density(pseudo_pure_density)

    if add_CNOT:
        new_NMRsample.read_and_plot("Figure/Pseudo/{}/CNOTpseudoaveragef{}".format(ufstring, ufstring))
    elif add_approx_CNOT:
        new_NMRsample.read_and_plot("Figure/Pseudo/{}/approxCNOTpseudoaveragef{}".format(ufstring, ufstring))
    else:
        new_NMRsample.read_and_plot("Figure/Pseudo/{}/pseudoaveragef{}".format(ufstring, ufstring))

    print(pseudo_pure_density)


def pseudo_pure_all_cases():
    pseudo_pure_state([0, 0], add_CNOT=False, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([0, 1], add_CNOT=False, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([1, 0], add_CNOT=False, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([1, 1], add_CNOT=False, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)


def CNOT_all_cases():
    pseudo_pure_state([0, 0], add_CNOT=True, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([0, 1], add_CNOT=True, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([1, 0], add_CNOT=True, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([1, 1], add_CNOT=True, add_approx_CNOT=False, generateprogram=True)
    time.sleep(2)


def ApproxCNOT_all_cases():
    pseudo_pure_state([0, 0], add_CNOT=False, add_approx_CNOT=True, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([0, 1], add_CNOT=False, add_approx_CNOT=True, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([1, 0], add_CNOT=False, add_approx_CNOT=True, generateprogram=True)
    time.sleep(2)
    pseudo_pure_state([1, 1], add_CNOT=False, add_approx_CNOT=True, generateprogram=True)
    time.sleep(2)


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
    NMRsample.read_and_plot("Figure/spectruma-")
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
    NMRsample.read_and_plot("Figure/spectrumb-")
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
    NMRsample.read_and_plot("Figure/spectrumc-")


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
    NMRsample.read_and_plot("Figure/spectrumd-")


def plot_data():
    NMRsample = chloroform()
    NMRsample.load_data_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_H_P0.csv")
    # NMRsample.read_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")


from Simulator.algorithm import *
import time

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

    # pseudo_pure_state([1, 1], add_CNOT=True)
    # spectrum_only_a()
    # spectrum_only_b()
    # spectrum_only_c()
    # spectrum_only_d()

    # P1_pulse()
    # CNOT_all_cases()
    # ApproxCNOT_all_cases()
    # Grover_print_pulse([1, 0, 0, 0])
    # plot_data()

    # pseudo_pure_state([0,0], add_CNOT=False, add_approx_CNOT=True)
    # pseudo_pure_state([0,1], add_CNOT=False, add_approx_CNOT=True)
    # pseudo_pure_state([1,0], add_CNOT=False, add_approx_CNOT=True)
    # pseudo_pure_state([1,1], add_CNOT=False, add_approx_CNOT=True)

    #permute_DJ([0, 0])
    #permute_DJ([0, 1])
    #permute_DJ([1, 0])
    #permute_DJ([1, 1])
    permute_grover([0, 0])
    permute_grover([0, 1])
    permute_grover([1, 0])
    permute_grover([1, 1])
    # strings = Grover_print_pulse([1, 0, 0, 0])

    # print("SSS")
    # print(strings)

    generate_grover_program("00")
    generate_grover_program("01")
    generate_grover_program("10")
    generate_grover_program("11")

    #generate_DJ_program(1)
    #generate_DJ_program(2)
    #generate_DJ_program(3)
    #generate_DJ_program(4)

    #pseudo_pure_all_cases()

    # CNOT_all_cases()
    # ApproxCNOT_all_cases()
    # time.sleep(2)
    #CNOT_all_cases()
    #new_CNOT_pulse()
