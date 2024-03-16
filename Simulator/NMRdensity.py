import numpy as np
from Simulator.Pulses import pulse, pulseSingle, pulseTwo, delayTime, BarrierPulse
from Simulator.params import *
import pandas as pd
from Simulator.util import *


def thermal_equilibrium_density():
    thermal_density = 1 / 4 * np.identity(4, dtype=complex) + 10 ** (-4) * np.array(
        [[5, 0, 0, 0], [0, 3, 0, 0], [0, 0, -3, 0], [0, 0, 0, -5]], dtype=complex)

    return thermal_density


class chloroform:
    '''
    Initialize a chloroform instance
    :param
            wTMS: The TMS frequency, in unit of MHZ
            dwell: Dwell time, in unit of second
            Npoints: Number of points
            T1p: T1 time of proton, in unit of second
            T2p: T2 time of proton, in unit of second
            T2starp: T2* time of proton, in unit of second
            T1C: T1 time of carbon, in unit of second
            T2C: T2 time of carbon, in unit of second
            T2starC: T2* time of carbon, in unit of second
            wH: Lamor frequency of proton, in unit of MHz
            cfH: Center frequency of proton, in unit of PPM
            pl90H: The calibrated 90 degree pulse for proton
            wC: Lamor frequency of carbon, in unit of MHz
            cfC: Center frequency of carbon, in unit of PPM
            pl90C: The calibrated 90 degree pulse for carbon
            Jfreq: J coupling frequency, in unit of Hz
    '''

    def __init__(self,
                 wTMS=wTMS,
                 dwell=dwell,
                 Npoints=Npoints,
                 T1p=T1p,
                 T2p=T2p,
                 T2starp=T2starp,
                 T1C=T1C,
                 T2C=T2C,
                 T2starC=T2starC,
                 wH=wH,
                 cfH=cfH,
                 pl90H=pl90H,
                 wC=wC,
                 cfC=cfC,
                 pl90C=pl90C,
                 Jfreq=Jfreq,
                 ):
        self._wTMS = wTMS
        self._dwell = dwell
        self._Npoints = Npoints

        self._times = np.linspace(start=0, stop=(self._Npoints - 1) * self._dwell, num=self._Npoints)
        self._proton_time_domain = []
        self._proton_freq_domain = []
        self._proton_freq_ppm = []

        '''
        Store the real ppm data and the 
        spectrum data for proton
        '''
        self._data_proton_ppm = []
        self._data_proton_freq_domain_real = []
        self._data_proton_peaks_real = []
        self._data_proton_peaks_pos_real = []
        self._data_proton_freq_domain_imag = []
        self._data_proton_peaks_imag = []
        self._data_proton_peaks_pos_imag = []

        self._carbon_time_domain = []
        self._carbon_freq_domain = []
        self._carbon_freq_ppm = []

        '''
        Store the real ppm data and the 
        spectrum data for proton
        '''
        self._data_carbon_ppm = []
        self._data_carbon_freq_domain_real = []
        self._data_carbon_peaks_real = []
        self._data_carbon_peaks_pos_real = []
        self._data_carbon_freq_domain_imag = []
        self._data_carbon_peaks_imag = []
        self._data_carbon_peaks_pos_imag = []

        self._density = np.zeros((4, 4), dtype=complex)
        self._density[0][0] = 0
        self._density[1][1] = 0
        self._density[2][2] = 1
        self._density[3][3] = 0

        self._pulses = []
        self._pulses_evolved = False
        self._pulse_unitary = None

        self._T1p = T1p
        self._T2p = T2p
        self._T2starp = T2starp
        self._T1C = T1C
        self._T2C = T2C
        self._T2starC = T2starC
        self._wH = wH
        self._cfH = cfH
        self._pl90H = pl90H
        self._wC = wC
        self._cfC = cfC
        self._pl90C = pl90C
        self._Jfreq = Jfreq
        pass

    def get_times(self):
        return self._times

    '''
    Calculate the thermal equilibrium density density matrix 
    '''

    def set_thermal_equilibrium(self):
        self._density = thermal_equilibrium_density()

    def reset_proton_params(self,
                            T1p,
                            T2p,
                            T2starp,
                            wH,
                            cfH,
                            pl90H
                            ):
        self._T1p = T1p
        self._T2p = T2p
        self._T2starp = T2starp
        self._wH = wH
        self._cfH = cfH
        self._pl90H = pl90H

    def reset_carbon_params(self,
                            T1C,
                            T2C,
                            T2starC,
                            wC,
                            cfC,
                            pl90C
                            ):
        self._T1C = T1C
        self._T2C = T2C
        self._T2starC = T2starC
        self._wC = wC
        self._cfC = cfC
        self._pl90C = pl90C

    def add_pulse(self, pl: pulse):
        self._pulses.append(pl)
        self._pulses_evolved = False

    def set_pulses(self, pulses):
        self._pulses = pulses

    def get_density(self):
        return self._density

    def print_density(self):
        print(self._density)

    def set_density(self, density: np.ndarray):
        assert density.shape == self._density.shape
        self._density = density

    def evolve_density(self, matrix: np.ndarray):
        assert matrix.shape == self._density.shape
        new_rho = self._density
        new_rho = np.matmul(matrix, new_rho)
        matrix_dag = np.conj(matrix)
        matrix_dag = np.transpose(matrix_dag)
        new_rho = np.matmul(new_rho, matrix_dag)
        self._density = new_rho
        # print(self._density)

    def evolve_all_pulse(self):
        if self._pulses_evolved:
            return
        self._pulse_unitary = np.identity(4, dtype=complex)
        for pulse in self._pulses:
            if isinstance(pulse, BarrierPulse):
                continue

            if isinstance(pulse, pulseSingle):
                matrix = pulse.get_matrix(self._pl90H, self._pl90C)
                self._pulse_unitary = np.matmul(matrix, self._pulse_unitary)
            elif isinstance(pulse, pulseTwo):
                matrix = pulse.get_matrix(self._pl90H, self._pl90C)
                self._pulse_unitary = np.matmul(matrix, self._pulse_unitary)
            elif isinstance(pulse, delayTime):
                matrix = pulse.get_matrix(self._Jfreq)
                self._pulse_unitary = np.matmul(matrix, self._pulse_unitary)
        self.evolve_density(self._pulse_unitary)
        self._pulses_evolved = True

    def get_pulse_unitary(self):
        return self._pulse_unitary

    def read_proton_time(self):
        self._proton_time_domain = []
        for t in self._times:
            self._proton_time_domain.append(self.measure_proton(t))
        return self._proton_time_domain

    def add_p1_perm_pulse(self):
        self.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
        self.add_pulse(delayTime(0.5 / Jfreq))
        self.add_pulse(pulseTwo(3, 0.5 * pl90C, wC, 0, 0.5 * pl90H, wH))
        self.add_pulse(delayTime(0.5 / Jfreq))
        self.add_pulse(pulseSingle(3, 0.5 * pl90H, wH))

    def add_p2_perm_pulse(self):
        self.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
        self.add_pulse(delayTime(0.5 / Jfreq))
        self.add_pulse(pulseTwo(3, 0.5 * pl90H, wH, 0, 0.5 * pl90C, wC))
        self.add_pulse(delayTime(0.5 / Jfreq))
        self.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))

    def add_X_gate_first_pulse(self):
        self.add_pulse(pulseSingle(0, pl90H, wH))

    def add_X_gate_second_pulse(self):
        self.add_pulse(pulseSingle(0, pl90C, wC))

    def add_H_gate_first_pulse(self, approximate=False):
        if not approximate:
            self.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
            self.add_pulse(pulseSingle(0, 1 * pl90H, wH))
            self.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))
        else:
            self.add_pulse(pulseSingle(3, 1 / 2 * pl90H, wH))

    def add_H_gate_second_pulse(self, approximate=False):
        if not approximate:
            self.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
            self.add_pulse(pulseSingle(0, 1 * pl90C, wC))
            self.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))
        else:
            self.add_pulse(pulseSingle(3, 1 / 2 * pl90C, wC))

    def add_CZ_pulse(self, approximate=False, Hcontrol=True):
        if not approximate:
            if Hcontrol:
                '''
                Add pulse sequence for exact CZ gate
                (pi/2)Iz1---(pi/2)Iz2---(-2)Iz1Iz2
                Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
                '''
                self.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
                self.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
                self.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))

                self.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
                self.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
                self.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))

                '''
                The pulse of (-2)Iz1Iz2. The angle theta is actually (-\pi/2). However, since we cannot rotate 
                a minus angle, we should plus another (4\pi)
                '''
                self.add_pulse(delayTime((4 - 0.5) / Jfreq))
            else:
                '''
                Change H and C
                '''
                self.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
                self.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
                self.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))

                self.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
                self.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
                self.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))

                self.add_pulse(delayTime((4 - 0.5) / Jfreq))
        else:
            pass

    def add_CNOT_pulse(self, approximate=False, Hcontrol=True):
        if not approximate:
            if Hcontrol:
                '''
                Add pulse sequence for approximate h gate on carbon
                '''
                self.add_H_gate_second_pulse(approximate=False)

                '''
                Add pulse sequence for exact CZ gate
                '''
                self.add_CZ_pulse(Hcontrol=True, approximate=False)

                '''
                Add pulse sequence for approximate h gate on carbon
                '''
                self.add_H_gate_second_pulse(approximate=False)
            else:
                '''
                Add pulse sequence for approximate h gate on proton
                '''
                self.add_H_gate_first_pulse(approximate=False)

                '''
                Add pulse sequence for exact CZ gate
                '''
                self.add_CZ_pulse(Hcontrol=False, approximate=False)

                '''
                Add pulse sequence for approximate h gate on proton
                '''
                self.add_H_gate_first_pulse()
        else:
            if Hcontrol:
                self.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
                self.add_pulse(delayTime(0.5 / Jfreq))
                self.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))
            else:
                self.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
                self.add_pulse(delayTime(0.5 / Jfreq))
                self.add_pulse(pulseSingle(3, 0.5 * pl90H, wH))

    def read_proton_spectrum(self, normalize=True):

        fft_result = np.fft.fft(self._proton_time_domain)
        fft_freq = np.fft.fftfreq(self._Npoints, self._dwell)

        # Combine the lists into a list of tuples and sort them by frequency
        combined_list = sorted(zip(list(fft_freq), list(fft_result)), key=lambda x: x[0])

        # Unzip the combined list back into two lists
        sorted_fft_freq, sorted_fft_result = zip(*combined_list)

        # If you need the results as lists (zip returns tuples)
        fft_freq = list(sorted_fft_freq)
        fft_result = list(sorted_fft_result)

        fft_freq = [x * (10 ** 6) / wTMS + cfH for x in fft_freq]

        self._proton_freq_domain = fft_result
        self._proton_freq_ppm = fft_freq
        '''
        Normalize the spectrum by integral value
        '''
        if normalize:
            integral_real = abs(np.trapz(self._proton_freq_ppm, abs(np.real(self._proton_freq_domain))))
            integral_imag = abs(np.trapz(self._proton_freq_ppm, abs(np.imag(self._proton_freq_domain))))
            self._proton_freq_domain = np.real(self._proton_freq_domain) / integral_real + 1j * np.imag(
                self._proton_freq_domain) / integral_imag

        return self._proton_freq_ppm, self._proton_freq_domain

    def read_carbon_time(self):
        self._carbon_time_domain = []
        for t in self._times:
            self._carbon_time_domain.append(self.measure_Carbon(t))
        return self._carbon_time_domain

    def read_carbon_spectrum(self, normalize=True):

        fft_result = np.fft.fft(self._carbon_time_domain)
        fft_freq = np.fft.fftfreq(self._Npoints, self._dwell)

        # Combine the lists into a list of tuples and sort them by frequency
        combined_list = sorted(zip(list(fft_freq), list(fft_result)), key=lambda x: x[0])

        # Unzip the combined list back into two lists
        sorted_fft_freq, sorted_fft_result = zip(*combined_list)

        # If you need the results as lists (zip returns tuples)
        fft_freq = list(sorted_fft_freq)
        fft_result = list(sorted_fft_result)

        fft_freq = [x * (10 ** 6) / wTMS + cfC for x in fft_freq]

        self._carbon_freq_domain = fft_result
        self._carbon_freq_ppm = fft_freq
        '''
        Normalize the spectrum by integral value
        '''
        if normalize:
            integral_real = abs(np.trapz(self._carbon_freq_ppm, abs(np.real(self._carbon_freq_domain))))
            integral_imag = abs(np.trapz(self._carbon_freq_ppm, abs(np.imag(self._carbon_freq_domain))))
            self._carbon_freq_domain = np.real(self._carbon_freq_domain) / integral_real + 1j * np.imag(
                self._carbon_freq_domain) / integral_imag

        return self._carbon_freq_ppm, self._carbon_freq_domain

    '''
    Measure the density matrix 
    '''

    def measure_proton(self, t):
        Mp = self.Mp_matrix(t)
        return np.trace(np.matmul(self._density, Mp)) * np.exp(-t / self._T2starp)

    def integral_proton_spectrum_real(self):
        return abs(np.trapz(self._proton_freq_ppm, abs(np.real(self._proton_freq_domain))))

    def integral_proton_spectrum_imag(self):
        return abs(np.trapz(self._proton_freq_ppm, abs(np.imag(self._proton_freq_domain))))

    def integral_carbon_spectrum_real(self):
        return abs(np.trapz(self._carbon_freq_ppm, abs(np.real(self._carbon_freq_domain))))

    def integral_carbon_spectrum_imag(self):
        return abs(np.trapz(self._carbon_freq_ppm, abs(np.imag(self._carbon_freq_domain))))

    def measure_Carbon(self, t):
        MC = self.MC_matrix(t)
        return np.trace(np.matmul(self._density, MC)) * np.exp(-t / self._T2starC)

    def Mp_matrix(self, t):
        J = self._Jfreq * 2 * np.pi
        Mp = np.array([[np.exp(-1j * J * t / 2), 0, -1j * np.exp(-1j * J * t / 2), 0],
                       [0, np.exp(1j * J * t / 2), 0, -1j * np.exp(1j * J * t / 2)],
                       [-1j * np.exp(-1j * J * t / 2), 0, -np.exp(-1j * J * t / 2), 0],
                       [0, -1j * np.exp(1j * J * t / 2), 0, -np.exp(1j * J * t / 2)]]
                      )
        return Mp

    def MC_matrix(self, t):
        J = self._Jfreq * 2 * np.pi
        MC = np.array([[np.exp(-1j * J * t / 2), -1j * np.exp(-1j * J * t / 2), 0, 0],
                       [-1j * np.exp(-1j * J * t / 2), -np.exp(-1j * J * t / 2), 0, 0],
                       [0, 0, np.exp(1j * J * t / 2), -1j * np.exp(1j * J * t / 2)],
                       [0, 0, -1j * np.exp(1j * J * t / 2), -np.exp(1j * J * t / 2)]]
                      )
        return MC

    def read_and_plot(self, path):
        '''
        Read the data signal in the time domain
        '''
        self.read_proton_time()
        self.read_carbon_time()
        '''
        Read the spectrum
        '''
        self.read_proton_spectrum()
        self.read_carbon_spectrum()
        '''
        Simulate what is shown on the screen
        '''
        self.show_proton_spectrum_real(-5, 15, store=True,
                                       path=path + "proton.png")

        self.show_carbon_spectrum_real(74, 80, store=True,
                                       path=path + "carbon.png")

    def determine_peaks(self, isproton=True, isreal=True):
        if isproton:
            if isreal:
                positions, values = find_two_largest_peaks(self._data_proton_ppm, self._data_proton_freq_domain_real)
                self._data_proton_peaks_real = values
                self._data_proton_peaks_pos_real = positions
            else:
                positions, values = find_two_largest_peaks(self._data_proton_ppm, self._data_proton_freq_domain_imag)
                self._data_proton_peaks_imag = values
                self._data_proton_peaks_pos_imag = positions
        else:
            if isreal:
                positions, values = find_two_largest_peaks(self._data_carbon_ppm, self._data_carbon_freq_domain_real)
                self._data_carbon_peaks_real = values
                self._data_carbon_peaks_pos_real = positions
            else:
                positions, values = find_two_largest_peaks(self._data_carbon_ppm, self._data_carbon_freq_domain_imag)
                self._data_carbon_peaks_imag = values
                self._data_carbon_peaks_pos_imag = positions

    '''
    Load the spectrum data from a give path
    '''

    def load_data_and_plot(self, path, isproton=True, store=False,storePath=None):
        # Load the CSV file into a DataFrame
        data = pd.read_csv(path, header=None)
        data.columns = ['Frequency (ppm)', 'Real Part', 'Imaginary Part']

        # Plotting
        plt.figure(figsize=(10, 6))

        if isproton:
            self._data_proton_ppm = list(data['Frequency (ppm)'])
            self._data_proton_freq_domain_real = list(data['Real Part'])
            self._data_proton_freq_domain_imag = list(data['Imaginary Part'])
            # Find the peaks for real and imaginary data
            self.determine_peaks(isproton=True, isreal=True)
            self.determine_peaks(isproton=True, isreal=False)
        else:
            self._data_carbon_ppm = list(data['Frequency (ppm)'])
            self._data_carbon_freq_domain_real = list(data['Real Part'])
            self._data_carbon_freq_domain_imag = list(data['Imaginary Part'])
            # Find the peaks for real and imaginary data
            self.determine_peaks(isproton=False, isreal=True)
            self.determine_peaks(isproton=False, isreal=False)
        # Plot the real part of the spectrum
        plt.plot(data['Frequency (ppm)'], data['Real Part'], label='Real Part', color='blue')

        if isproton:
            plt.scatter(self._data_proton_peaks_pos_real[0], self._data_proton_peaks_real[0], color="red",
                        label="First peak f={}, p={}".format(self._data_proton_peaks_pos_real[0],
                                                             self._data_proton_peaks_real[0]))
            plt.scatter(self._data_proton_peaks_pos_real[1], self._data_proton_peaks_real[1], color="red",
                        label="First peak f={}, p={}".format(self._data_proton_peaks_pos_real[1],
                                                             self._data_proton_peaks_real[1]))
        else:
            plt.scatter(self._data_carbon_peaks_pos_real[0], self._data_carbon_peaks_real[0], color="red",
                        label="First peak f={}, p={}".format(self._data_carbon_peaks_pos_real[0],
                                                             self._data_carbon_peaks_real[0]))
            plt.scatter(self._data_carbon_peaks_pos_real[1], self._data_carbon_peaks_real[1], color="red",
                        label="First peak f={}, p={}".format(self._data_carbon_peaks_pos_real[1],
                                                             self._data_carbonpeaks_real[1]))
        # Optionally, plot the imaginary part of the spectrum on the same plot
        # Uncomment the next line if you want to include the imaginary part in the plot
        # plt.plot(data['Frequency (ppm)'], data['Imaginary Part'], label='Imaginary Part', color='red')

        # Adding plot title and labels
        if isproton:
            plt.title('Real part of proton NMR Spectrum')
        else:
            plt.title('Real part of carbon NMR Spectrum')
        plt.xlabel('Frequency (ppm)')
        plt.ylabel('Spectrum Intensity')
        plt.legend()

        if store:
            plt.savefig(storePath)

        # Display the plot
        plt.show()

    def show_proton_fid_real(self, maxtime, store=False, path=None, title="Proton FID real", miny=-1, maxy=1):
        plt.plot(self._times, np.real(self._proton_time_domain), label="Proton test_fid.py(Real part)")
        plt.xlim(0, maxtime)
        plt.ylim(miny, maxy)
        plt.xlabel("Time/second")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_proton_fid_imag(self, maxtime, store=False, path=None, title="Proton FID imaginary", miny=-1, maxy=1):
        plt.plot(self._times, np.imag(self._proton_time_domain), label="Proton test_fid.py(Imaginary part)")
        plt.xlim(0, maxtime)
        plt.ylim(miny, maxy)
        plt.xlabel("Time/second")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_proton_spectrum_real(self, minppm, maxppm, store=False, path=None, title="Proton spectrum real"):
        plt.plot(self._proton_freq_ppm, np.real(self._proton_freq_domain), label="Proton spectrum(Real part)")
        plt.xlim(minppm, maxppm)
        plt.axvline(x=cfH, color="red", linestyle="--", label="Center frequency of proton")
        plt.xlabel("Frequency/ppm")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_proton_spectrum_imag(self, minppm, maxppm, store=False, path=None, title="Proton spectrum imaginary"):
        plt.plot(self._proton_freq_ppm, np.imag(self._proton_freq_domain), label="Proton spectrum(Imaginary part)")
        plt.xlim(minppm, maxppm)
        plt.axvline(x=cfH, color="red", linestyle="--", label="Center frequency of proton")
        plt.xlabel("Frequency/ppm")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_carbon_fid_real(self, maxtime, store=False, path=None, miny=-1, maxy=1, title="Carbon FID real"):
        plt.plot(self._times, np.real(self._carbon_time_domain), label="Carbon test_fid.py(Real part)")
        plt.xlim(0, maxtime)
        plt.ylim(miny, maxy)
        plt.xlabel("Time/second")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_carbon_fid_imag(self, maxtime, store=False, path=None, miny=-1, maxy=1, title="Carbon FID imaginary"):
        plt.plot(self._times, np.real(self._carbon_freq_domain), label="Carbon test_fid.py(Imaginary part)")
        plt.xlim(0, maxtime)
        plt.ylim(miny, maxy)
        plt.xlabel("Time/second")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_carbon_spectrum_real(self, minppm, maxppm, store=False, path=None, title="Carbon spectrum real"):
        plt.plot(self._carbon_freq_ppm, np.real(self._carbon_freq_domain), label="Carbon spectrum(Real part)")
        plt.xlim(minppm, maxppm)
        plt.axvline(x=cfC, color="red", linestyle="--", label="Center frequency of carbon")
        plt.xlabel("Frequency/ppm")
        plt.title(title)
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def show_carbon_spectrum_imag(self, minppm, maxppm, store=False, path=None, title="Carbon spectrum imaginary"):
        plt.plot(self._carbon_freq_ppm, np.imag(self._carbon_freq_domain), label="Carbon spectrum(Imaginary part)")
        plt.xlim(minppm, maxppm)
        plt.axvline(x=cfC, color="red", linestyle="--", label="Center frequency of carbon")
        plt.xlabel("Frequency/ppm")
        plt.legend()
        if store:
            plt.savefig(path)
        plt.show()

    def print_pulses(self):
        print("initpp(dir)")
        for pulse in self._pulses:
            print(pulse)
        print("--------------------------------")
        print("parList=endpp()")


import matplotlib.pyplot as plt

if __name__ == "__main__":
    NMRsample = chloroform()

    proton_time_domain = NMRsample.read_proton_time()

    print(proton_time_domain)
