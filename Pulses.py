import numpy as np
from params import *


def I_matrix():
    return np.array([[1, 0], [0, 1]])


def Rx_matrix(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                     [-1j * np.sin(theta / 2), np.cos(theta / 2)]])


def Rx_matrix_first(theta):
    return np.kron(
        np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)],
                  [-1j * np.sin(theta / 2), np.cos(theta / 2)]]),
        I_matrix())


def Rx_matrix_second(theta):
    return np.kron(I_matrix(), np.array(
        [[np.cos(theta / 2), -1j * np.sin(theta / 2)],
         [-1j * np.sin(theta / 2), np.cos(theta / 2)]]))


def Ry_matrix(theta):
    return np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)],
                     [np.sin(theta / 2), np.cos(theta / 2)]])


def Ry_matrix_first(theta):
    return np.kron(np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)],
                             [np.sin(theta / 2), np.cos(theta / 2)]]),
                   I_matrix())


def Ry_matrix_second(theta):
    return np.kron(I_matrix(),
                   np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)],
                             [np.sin(theta / 2), np.cos(theta / 2)]]))


def Rzz_matrix(theta):
    return np.array([[np.exp(-1j * theta / 2), 0, 0, 0],
                     [0, np.exp(1j * theta / 2), 0, 0],
                     [0, 0, np.exp(1j * theta / 2), 0],
                     [0, 0, 0, np.exp(-1j * theta / 2)]])


class pulse:

    def __init__(self):
        pass

    def get_matrix(self, *params):
        raise NotImplementedError


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

    def get_matrix(self,pl90lengthp,pl90lengthC):
        if self._is_carbon:
            theta = (self._length / pl90lengthC) * np.pi / 2
        else:
            theta = (self._length / pl90lengthp) * np.pi / 2
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
        self._is_carbon = False
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
            theta1 = (self._length1 / pl90lengthC) * np.pi / 2
            theta2 = (self._length2 / pl90lengthp) * np.pi / 2
        else:
            theta1 = (self._length1 / pl90lengthp) * np.pi / 2
            theta2 = (self._length2 / pl90lengthC) * np.pi / 2

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


class delayTime(pulse):
    '''
    Initialize a delay time instance
    :param:
          delaytime: The delaytime
    '''

    def __init__(self, delaytime):
        super().__init__()
        self._delaytime = delaytime

    def get_matrix(self, Jfreq):
        theta = np.pi * Jfreq * self._delaytime
        return Rzz_matrix(theta)
