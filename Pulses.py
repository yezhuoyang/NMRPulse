class pulse:

    def __init__(self):
        pass


'''
Single qubit NMR pulse
'''


class pulseSingle(pulse):
    '''
    Initialize a single pulse instance
    :param:
          channel: Which channel are you adding the pulse
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


'''
Two qubits NMR pulse
'''


class pulseTwo(pulse):
    '''
    Initialize a single pulse instance
    :param:
          channel: Which channel are you adding the pulse
          length: The length of the pulse
          freq: The frequency of the pulse
    '''

    def __init__(self,
                 channel1,
                 length1,
                 freq1,
                 channel2,
                 length2,
                 freq2,
                 ):
        self._channel1 = channel1
        self._length1 = length1
        self._freq1 = freq1
        self._channel2 = channel2
        self._length2 = length2
        self._freq2 = freq2
        pass


class delayTime(pulse):
    '''
    Initialize a delay time instance
    :param:
          channel: Which channel are you adding the pulse
          length: The length of the pulse
          freq: The frequency of the pulse
    '''

    def __init__(self,
                 delaytime
                 ):
        self._delaytime = delaytime
        pass
