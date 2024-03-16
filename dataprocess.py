from Simulator.NMRdensity import *
from Simulator.Pulses import *
import matplotlib.pyplot as plt

'''
Store the real ppm data and the 
spectrum data for proton
'''


def average_spectrum(sample0: chloroform, sample1: chloroform, sample2: chloroform):
    data_proton_ppm0 = sample0._data_proton_ppm
    data_proton_freq_domain_real0 = np.array(sample0._data_proton_freq_domain_real)
    data_proton_peaks_pos_real0 = sample0._data_proton_peaks_pos_real

    data_proton_ppm1 = sample1._data_proton_ppm
    data_proton_freq_domain_real1 = np.array(sample1._data_proton_freq_domain_real)
    data_proton_peaks_pos_real1 = sample1._data_proton_peaks_pos_real

    data_proton_ppm2 = sample2._data_proton_ppm
    data_proton_freq_domain_real2 = np.array(sample2._data_proton_freq_domain_real)
    data_proton_peaks_pos_real2 = sample2._data_proton_peaks_pos_real

    assert (data_proton_ppm0 == data_proton_ppm1 == data_proton_ppm2)

    avg_data_freq_real = (data_proton_freq_domain_real0 + data_proton_freq_domain_real1 + data_proton_freq_domain_real2) / 3

    plt.plot(data_proton_ppm0, avg_data_freq_real)
    plt.show()

def plot_data():
    NMRsample = chloroform()
    NMRsample.load_data_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")
    # NMRsample.read_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")


def process_CNOT_data():
    return


def process_approx_CNOT_data():
    return


def process_pseudo_state00_data():
    path00proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P0.csv"

    sample00_P0 = chloroform()
    sample00_P0.load_data_and_plot(path00proton_P0, isproton=True)

    path00proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P1.csv"

    sample00_P1 = chloroform()
    sample00_P1.load_data_and_plot(path00proton_P1, isproton=True)

    path00proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P2.csv"

    sample00_P2 = chloroform()
    sample00_P2.load_data_and_plot(path00proton_P2, isproton=True)

    average_spectrum(sample00_P0, sample00_P1, sample00_P2)

    return


def process_pseudo_state01_data():
    path00proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P0.csv"

    sample00_P0 = chloroform()
    sample00_P0.load_data_and_plot(path00proton_P0, isproton=True)

    path00proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P1.csv"

    sample00_P1 = chloroform()
    sample00_P1.load_data_and_plot(path00proton_P1, isproton=True)

    path00proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P2.csv"

    sample00_P2 = chloroform()
    sample00_P2.load_data_and_plot(path00proton_P2, isproton=True)

    return


def process_DJ_data():
    return


def process_Grover_data():
    return


from Simulator.algorithm import *

if __name__ == "__main__":
    process_pseudo_state00_data()
