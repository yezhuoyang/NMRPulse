import numpy as np
from .Pulses import pulse, pulseSingle, pulseTwo, delayTime

wTMS = 62
T1p = 10
T2p = 2
T2starp = 200
T1C = 10
T2C = 2
T2starC = 50
wH = 62.3750106999999970
cfH = 5
pl90H = 13
wC = 15.6858599000000010
cfC = 77
pl90C = 77
Jfreq=215
















class chloroform:
    '''
    Initialize a chloroform instance
    :param
            wTMS: The TMS frequency, in unit of MHZ
            T1p: T1 time of proton, in unit of second
            T2p: T2 time of proton, in unit of second
            T2starp: T2* time of proton, in unit of ms
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
                 pl90C=pl90H,
                 Jfreq=Jfreq,
                 ):
        self._wTMS = wTMS
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
        self._Jfreq=Jfreq
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



    def evolve_density(self,matrix:np.ndarray):
        assert matrix.shape==self._density.shape
        new_rho=self._density
        new_rho=np.matmul(matrix,new_rho)
        matrix_dag=np.conj(matrix)
        matrix_dag=np.transpose(matrix_dag)
        new_rho=np.matmul(new_rho,matrix_dag)
        return new_rho



    def evolve_all_pulse(self):
        for pulse in self._pulses:
            if isinstance(pulse, pulseSingle):
                pass
            elif isinstance(pulse, pulseTwo):
                pass
            elif isinstance(pulse, delayTime):
                pass


    '''
    Measure the density matrix 
    
    '''

    def measure(self, channel):
        pass



    def Mp_matrix(self):
        return




    def Mc_matrix(self):
        return





if __name__ == "__main__":
    pass