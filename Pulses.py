import numpy as np


def I_matrix():
    return np.array([[1, 0], [0, 1]])


def Rx_matrix(theta):
    return np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)], [-1j * np.sin(theta / 2), np.cos(theta / 2)]])


def Rx_matrix_first(theta):
    return np.kron(np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)], [-1j * np.sin(theta / 2), np.cos(theta / 2)]]),I_matrix())


def Rx_matrix_second(theta):
    return np.kron(I_matrix(),np.array([[np.cos(theta / 2), -1j * np.sin(theta / 2)], [-1j * np.sin(theta / 2), np.cos(theta / 2)]]))


def Ry_matrix(theta):
    return np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)], [np.sin(theta / 2), np.cos(theta / 2)]])


def Ry_matrix_first(theta):
    return np.kron(np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)], [np.sin(theta / 2), np.cos(theta / 2)]]),I_matrix())


def Ry_matrix_second(theta):
    return np.kron(I_matrix(),np.array([[np.cos(theta / 2), -1 * np.sin(theta / 2)], [np.sin(theta / 2), np.cos(theta / 2)]]))


def Rzz_matrix(theta):
    return np.array([[np.exp(-1j * theta / 2), 0, 0, 0],
                     [0, np.exp(1j * theta / 2), 0, 0],
                     [0, 0, np.exp(1j * theta / 2), 0],
                     [0, 0, 0, np.exp(-1j * theta / 2)]])


class pulse:

    def __init__(self):
        pass

    def get_matrix(self, pl90length):
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
          freq: The frequency of the pulse
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
        pass

    def get_matrix(self, pl90length):
        raise NotImplementedError


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
        self._channel2 = channel2
        self._length2 = length2
        self._freq2 = freq2

    def get_matrix(self, pl90length):
        raise NotImplementedError


class delayTime(pulse):
    '''
    Initialize a delay time instance
    :param:
          delaytime: The delaytime
    '''

    def __init__(self, delaytime):
        super().__init__()
        self._delaytime = delaytime

    def get_matrix(self, pl90length):
        raise NotImplementedError
