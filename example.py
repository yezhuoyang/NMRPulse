from NMRdensity import *
from params import *
from Pulses import *
import matplotlib.pyplot as plt

if __name__ == "__main__":
    NMRsample = chloroform()

    #NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))

    #NMRsample.evolve_all_pulse()
    NMRsample.read_proton_time()
    NMRsample.read_proton_spectrum()
    NMRsample.show_proton_spectrum_real(-5, 15)
