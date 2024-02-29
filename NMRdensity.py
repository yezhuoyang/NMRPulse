import numpy as np
from .Pulses import pulse, pulseSingle, pulseTwo, delayTime


class chloroform:
    '''
    Initialize a chloroform instance
    :param
            T1p: T1 time of proton
            T2p: T2 time of proton
            T2starp: T2* time of proton
            T1C: T1 time of carbon
            T2C: T2 time of carbon
            T2starC: T2* time of carbon
            wH: Lamor frequency of proton
            cfH: Center frequency of proton
            pl90H: The calibrated 90 degree pulse for proton
            wC: Lamor frequency of carbon
            cfC: Center frequency of carbon
            pl90C: The calibrated 90 degree pulse for carbon
    '''

    def __init__(self,
                 T1p,
                 T2p,
                 T2starp,
                 T1C,
                 T2C,
                 T2starC,
                 wH,
                 cfH,
                 pl90H,
                 wC,
                 cfC,
                 pl90C,
                 ):
        self._density = np.zeros((2, 2), dtype=complex)
        self._density[0][0] = 1
        self._pulses = []
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
        pass

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

    def get_density(self):
        return self._density

    def print_density(self):
        print(self._density)

    def set_density(self, density: np.ndarray):
        assert density.shape == self._density.shape
        self._density = density

    def evolve_all_pulse(self):
        for pulse in self._pulses:
            if isinstance(pulse, pulseSingle):
                pass
            elif isinstance(pulse, pulseTwo):
                pass
            elif isinstance(pulse, delayTime):
                pass

    def evolve(self, acquiretime):
        self.evolve_all_pulse()
        pass

    '''
    Measure the density matrix 
    
    '''

    def measure(self, channel):
        pass
