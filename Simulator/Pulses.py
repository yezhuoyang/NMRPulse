import numpy as np
from Simulator.params import *


def I_matrix():
    return np.array([[1, 0], [0, 1]], dtype=complex)


def Rx_matrix(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex)


def Rz_matrix(theta):
    return np.array([[np.exp(-1j * theta / 2), 0],
                     [0, np.exp(1j * theta / 2)]], dtype=complex)


def Rx_matrix_first(theta):
    return np.kron(
        np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                  [-1j * np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex),
        I_matrix())


def Rx_matrix_second(theta):
    return np.kron(I_matrix(), np.array(
        [[np.cos(theta / 2), -1j * np.sin(theta / 2)],
         [-1j * np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex))


def Ry_matrix(theta):
    return np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)],
                     [np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex)


def Ry_matrix_first(theta):
    return np.kron(np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)],
                             [np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex),
                   I_matrix())


def Ry_matrix_second(theta):
    return np.kron(I_matrix(),
                   np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)],
                             [np.sin(theta / 2), np.cos(theta / 2)]], dtype=complex))


def Rzz_matrix(theta):
    return np.array([[np.exp(-1j * theta / 2), 0, 0, 0],
                     [0, np.exp(1j * theta / 2), 0, 0],
                     [0, 0, np.exp(1j * theta / 2), 0],
                     [0, 0, 0, np.exp(-1j * theta / 2)]], dtype=complex)


class pulse:

    def __init__(self):
        pass

    def get_matrix(self, *params):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError


class BarrierPulse(pulse):

    def __init__(self, name, endpulse=False):
        self.name = name
        self.endpulse = endpulse

    def get_matrix(self, *params):
        raise NotImplementedError

    def __str__(self):
        if not self.endpulse:
            return "############Pulse for {} gate:###########".format(self.name)
        else:
            return "###########Pulse for {} gate end Here####".format(self.name)


'''
Single qubit NMR pulse
'''


class pulseSingle(pulse):
    '''
    Initialize a single pulse instance
    :param:
          channel: Which channel are you adding the pulse, 0 for +x, 1 for +y, 2 for -x, 3 for -y
          length: The length of the pulse
          freq: The frequency of the pulse.
    '''

    def __init__(self,
                 channel,
                 length,
                 freq,
                 ):
        super().__init__()
        self._channel = channel
        self._length = length
        self._freq = freq

        assert freq == wC or freq == wH
        self._is_carbon = False
        if freq == wC:
            self._is_carbon = True

    def get_matrix(self, pl90lengthp, pl90lengthC):
        if self._is_carbon:
            theta = (self._length / pl90lengthC) * np.pi
        else:
            theta = (self._length / pl90lengthp) * np.pi
        if self._channel == 0:
            if self._is_carbon:
                return Rx_matrix_second(theta)
            else:
                return Rx_matrix_first(theta)
        elif self._channel == 1:
            if self._is_carbon:
                return Ry_matrix_second(theta)
            else:
                return Ry_matrix_first(theta)
        elif self._channel == 2:
            if self._is_carbon:
                return Rx_matrix_second(-theta)
            else:
                return Rx_matrix_first(-theta)
        elif self._channel == 3:
            if self._is_carbon:
                return Ry_matrix_second(-theta)
            else:
                return Ry_matrix_first(-theta)

    '''
    Print the string format of the pulse sequence
    '''

    def __str__(self):
        if self._is_carbon:
            prop = self._length / pl90C
            if prop == 0.5:
                return "pulse(2,a90C,{},d90C)".format(self._channel)
            elif prop == 1:
                return "pulse(2,a90C,{},d180C)".format(self._channel)
            else:
                return "pulse(2,a90C,{},{}d90C)".format(self._channel, prop)
        else:
            prop = self._length / pl90H
            if prop == 0.5:
                return "pulse(1,a90H,{},d90H)".format(self._channel)
            elif prop == 1:
                return "pulse(1,a90H,{},d180H)".format(self._channel)
            else:
                return "pulse(1,a90H,{},{}d90H)".format(self._channel, prop)


'''
Two qubits NMR pulse
'''


class pulseTwo(pulse):
    '''
    Initialize a single pulse instance
    :param:
          channel1: Which channel are you adding the first pulse, 0 for +x, 1 for +y, 2 for -x, 3 for -y
          length1: The length of the first pulse
          freq1: The frequency of the first pulse
          channel2: Which channel are you adding the second pulse, 0 for +x, 1 for +y, 2 for -x, 3 for -y
          length2: The length of the second pulse
          freq2: The frequency of the second pulse
    '''

    def __init__(self, channel1, length1, freq1, channel2, length2, freq2):
        super().__init__()
        self._channel1 = channel1
        self._length1 = length1
        self._freq1 = freq1
        assert freq1 == wC or freq1 == wH
        self._is_1_carbon = False
        if freq1 == wC:
            self._is_1_carbon = True

        self._channel2 = channel2
        self._length2 = length2
        self._freq2 = freq2

        assert freq1 != freq2 and (freq2 == wC or freq2 == wH)
        self._is_2_carbon = False
        if freq2 == wC:
            self._is_2_carbon = True

    def get_matrix(self, pl90lengthp, pl90lengthC):

        if self._is_1_carbon:
            theta1 = (self._length1 / pl90lengthC) * np.pi
            theta2 = (self._length2 / pl90lengthp) * np.pi
        else:
            theta1 = (self._length1 / pl90lengthp) * np.pi
            theta2 = (self._length2 / pl90lengthC) * np.pi

        matrix1 = None
        matrix2 = None

        if self._channel1 == 0:
            matrix1 = Rx_matrix(theta1)
        elif self._channel1 == 1:
            matrix1 = Ry_matrix(theta1)
        elif self._channel1 == 2:
            matrix1 = Rx_matrix(-theta1)
        elif self._channel1 == 3:
            matrix1 = Ry_matrix(-theta1)

        if self._channel2 == 0:
            matrix2 = Rx_matrix(theta2)
        elif self._channel2 == 1:
            matrix2 = Ry_matrix(theta2)
        elif self._channel2 == 2:
            matrix2 = Rx_matrix(-theta2)
        elif self._channel2 == 3:
            matrix2 = Ry_matrix(-theta2)

        if not self._is_1_carbon:
            return np.kron(matrix1, matrix2)
        else:
            return np.kron(matrix2, matrix1)

    def __str__(self):
        if self._is_1_carbon:
            return "pulse(1,a90HC,{},freq1H,2,a90C,{},freq13C,d90C)".format(self._channel2, self._channel1)
        else:
            return "pulse(1,a90HC,{},freq1H,2,a90C,{},freq13C,d90C)".format(self._channel1, self._channel2)


class delayTime(pulse):
    '''
    Initialize a delay time instance
    :param:
          delaytime: The delaytime
    '''

    def __init__(self, delaytime, isgapdelay=False):
        super().__init__()
        self._delaytime = delaytime
        self._isgapdelay = isgapdelay

    def get_matrix(self, Jfreq):
        theta = np.pi * Jfreq * self._delaytime
        return Rzz_matrix(theta)

    def __str__(self):
        if not self._isgapdelay:
            if self._delaytime == ((4 - 0.5) / Jfreq):
                return "delay(dCZEvolution)"
            else:
                return "delay(dEvolution)"
        else:
            return "delay(0.25)"
        # return "delay({:.2f})".format(self._delaytime * 10 ** (3))
