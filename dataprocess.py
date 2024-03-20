from Simulator.NMRdensity import *
from Simulator.Pulses import *
import matplotlib.pyplot as plt

'''
Store the real ppm data and the 
spectrum data for proton
'''


def average_proton_spectrum(sample0: chloroform, sample1: chloroform, sample2: chloroform, store=False, storetitle=None,
                            storePath=None):
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

    avg_data_freq_real = (data_proton_freq_domain_real0 +
                          data_proton_freq_domain_real1 +
                          data_proton_freq_domain_real2) / 3

    avg_sample = chloroform()
    avg_sample._data_proton_ppm = data_proton_ppm0
    avg_sample._data_proton_freq_domain_real = avg_data_freq_real
    avg_sample.determine_peaks(isproton=True, isreal=True)

    plt.plot(data_proton_ppm0, avg_data_freq_real, label="Averaged proton spectrum result")
    plt.xlabel("Frequency/ppm")

    plt.axvline(x=avg_sample._data_proton_peaks_pos_real[0] - proton_int_range, color="red", linestyle="--")
    plt.axvline(x=avg_sample._data_proton_peaks_pos_real[0] + proton_int_range, color="red", linestyle="--")
    plt.axvline(x=avg_sample._data_proton_peaks_pos_real[1] - proton_int_range, color="green", linestyle="--")
    plt.axvline(x=avg_sample._data_proton_peaks_pos_real[1] + proton_int_range, color="green", linestyle="--")
    plt.scatter(avg_sample._data_proton_peaks_pos_real[0], avg_sample._data_proton_peaks_real[0], color="red",
                label="First peak f={:.3f}, p={:.3f}, integral={:.3f}".format(avg_sample._data_proton_peaks_pos_real[0],
                                                                              avg_sample._data_proton_peaks_real[0],
                                                                              avg_sample._data_proton_peaks_integral_real[
                                                                                  0]))
    plt.scatter(avg_sample._data_proton_peaks_pos_real[1], avg_sample._data_proton_peaks_real[1], color="red",
                label="Second peak f={:.3f}, p={:.3f}, integral={:.3f}".format(
                    avg_sample._data_proton_peaks_pos_real[1],
                    avg_sample._data_proton_peaks_real[1],
                    avg_sample._data_proton_peaks_integral_real[1]))

    if store:
        plt.title(storetitle)
        plt.legend()
        plt.savefig(storePath)


def average_carbon_spectrum(sample0: chloroform, sample1: chloroform, sample2: chloroform, store=False, storetitle=None,
                            storePath=None):
    data_carbon_ppm0 = sample0._data_carbon_ppm
    data_carbon_freq_domain_real0 = np.array(sample0._data_carbon_freq_domain_real)
    data_carbon_peaks_pos_real0 = sample0._data_carbon_peaks_pos_real

    data_carbon_ppm1 = sample1._data_carbon_ppm
    data_carbon_freq_domain_real1 = np.array(sample1._data_carbon_freq_domain_real)
    data_carbon_peaks_pos_real1 = sample1._data_carbon_peaks_pos_real

    data_carbon_ppm2 = sample2._data_carbon_ppm
    data_carbon_freq_domain_real2 = np.array(sample2._data_carbon_freq_domain_real)
    data_carbon_peaks_pos_real2 = sample2._data_carbon_peaks_pos_real

    assert (data_carbon_ppm0 == data_carbon_ppm1 == data_carbon_ppm2)

    avg_data_freq_real = (data_carbon_freq_domain_real0 +
                          data_carbon_freq_domain_real1 +
                          data_carbon_freq_domain_real2) / 3

    avg_sample = chloroform()
    avg_sample._data_carbon_ppm = data_carbon_ppm0
    avg_sample._data_carbon_freq_domain_real = avg_data_freq_real
    avg_sample.determine_peaks(isproton=False, isreal=True)

    plt.plot(data_carbon_ppm0, avg_data_freq_real, label="Averaged carbon spectrum result")
    plt.xlabel("Frequency/ppm")

    plt.axvline(x=avg_sample._data_carbon_peaks_pos_real[0] - carbon_int_range, color="red", linestyle="--")
    plt.axvline(x=avg_sample._data_carbon_peaks_pos_real[0] + carbon_int_range, color="red", linestyle="--")
    plt.axvline(x=avg_sample._data_carbon_peaks_pos_real[1] - carbon_int_range, color="green", linestyle="--")
    plt.axvline(x=avg_sample._data_carbon_peaks_pos_real[1] + carbon_int_range, color="green", linestyle="--")
    plt.scatter(avg_sample._data_carbon_peaks_pos_real[0], avg_sample._data_carbon_peaks_real[0], color="red",
                label="First peak f={:.3f}, p={:.3f}, integral={:.3f}".format(avg_sample._data_carbon_peaks_pos_real[0],
                                                                              avg_sample._data_carbon_peaks_real[0],
                                                                              avg_sample._data_carbon_peaks_integral_real[
                                                                                  0]))
    plt.scatter(avg_sample._data_carbon_peaks_pos_real[1], avg_sample._data_carbon_peaks_real[1], color="red",
                label="Second peak f={:.3f}, p={:.3f}, integral={:.3f}".format(
                    avg_sample._data_carbon_peaks_pos_real[1],
                    avg_sample._data_carbon_peaks_real[1],
                    avg_sample._data_carbon_peaks_integral_real[1]))

    plt.plot(data_carbon_ppm0, avg_data_freq_real, label="Averaged carbon spectrum result")
    plt.xlabel("Frequency/ppm")
    if store:
        plt.title(storetitle)
        plt.legend()
        plt.savefig(storePath)
    print("Plot succeed!")


def plot_data():
    NMRsample = chloroform()
    NMRsample.load_data_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")
    # NMRsample.read_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")


def process_CNOT_data():
    return

def process_approx_CNOT_data_H_00():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/00/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/00/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/00/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/00/00_H_avg.png",
                            storetitle="Averaged result in H channel for DJ(ApproxCNOT 00 State).")
    return

def process_approx_CNOT_data_H_10():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/10_H_avg.png",
                            storetitle="Averaged result in H channel for DJ(ApproxCNOT 10 State).")
    return


def process_approx_CNOT_data_H_01():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/01/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/01/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/01/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/01/01_H_avg.png",
                            storetitle="Averaged result in H channel for DJ(ApproxCNOT 01 State).")
    return


def process_exact_CNOT_data_H_11():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/11/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/11/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/11/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/11/11_H_avg.png",
                            storetitle="Averaged result in H channel(Exact CNOT 11 State).")
    return


def process_exact_CNOT_data_H_00():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/00/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/00/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/00/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/00/00_H_avg.png",
                            storetitle="Averaged result in H channel(Exact CNOT 00 State).")
    return


def process_exact_CNOT_data_H_01():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/01/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/01/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/01/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/01/01_H_avg.png",
                            storetitle="Averaged result in H channel(Exact CNOT 01 State).")
    return

def process_exact_CNOT_data_H_10():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/10/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/10/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/10/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/CNOT/10/10_H_avg.png",
                            storetitle="Averaged result in H channel(Exact CNOT 10 State).")
    return


def process_approx_CNOT_data_H_11():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/P0_H.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/P1_H.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/P2_H.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_H_avg.png",
                            storetitle="Averaged result in H channel(ApproxCNOT 11 State).")
    return



def process_approx_CNOT_data_C_10():
    path10carbon_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/10_C_P0_apprxCNOT.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10carbon_P0, minppm=60, maxppm=85, isproton=False)

    path10carbon_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/10_C_P1_apprxCNOT.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10carbon_P1, minppm=60, maxppm=85, isproton=False)

    path10carbon_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/10_C_P2_apprxCNOT.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10carbon_P2, minppm=60, maxppm=85, isproton=False)

    average_carbon_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/10/10_C_avg.png",
                            storetitle="Averaged result in C channel(ApproxCNOT 10 State).")
    return


def process_approx_CNOT_data_H_11():
    path11proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_H_P0_apprxCNOT.csv"

    sample11_P0 = chloroform()
    sample11_P0.load_data_and_plot(path11proton_P0, isproton=True)

    path11proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_H_P1_apprxCNOT.csv"

    sample11_P1 = chloroform()
    sample11_P1.load_data_and_plot(path11proton_P1, isproton=True)

    path11proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_H_P2_apprxCNOT.csv"

    sample11_P2 = chloroform()
    sample11_P2.load_data_and_plot(path11proton_P2, isproton=True)

    average_proton_spectrum(sample11_P0, sample11_P1, sample11_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_H_avg.png",
                            storetitle="Averaged result in H channel(ApproxCNOT 11 State).")
    return





def process_approx_CNOT_data_C_11():
    path11carbon_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_C_P1_apprxCNOT.csv"

    sample11_P0 = chloroform()
    sample11_P0.load_data_and_plot(path11carbon_P0, minppm=60, maxppm=85, isproton=False)

    path11carbon_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_C_P1_apprxCNOT.csv"

    sample11_P1 = chloroform()
    sample11_P1.load_data_and_plot(path11carbon_P1, minppm=60, maxppm=85, isproton=False)

    path11carbon_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_C_P2_apprxCNOT.csv"

    sample11_P2 = chloroform()
    sample11_P2.load_data_and_plot(path11carbon_P2, minppm=60, maxppm=85, isproton=False)

    average_carbon_spectrum(sample11_P0, sample11_P1, sample11_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/ApproxCNOT/11/11_C_avg.png",
                            storetitle="Averaged result in C channel(ApproxCNOT 11 State).")
    return


def process_pseudo_state00_proton_data():
    path00proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P0.csv"

    sample00_P0 = chloroform()
    sample00_P0.load_data_and_plot(path00proton_P0, isproton=True)

    path00proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P1.csv"

    sample00_P1 = chloroform()
    sample00_P1.load_data_and_plot(path00proton_P1, isproton=True)

    path00proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_P2.csv"

    sample00_P2 = chloroform()
    sample00_P2.load_data_and_plot(path00proton_P2, isproton=True)

    average_proton_spectrum(sample00_P0, sample00_P1, sample00_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_H_avg.png",
                            storetitle="Averaged result in Proton channel")
    return


def process_pseudo_state00_carbon_data():
    path00carbon_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_C_P0.csv"

    sample00_P0 = chloroform()
    sample00_P0.load_data_and_plot(minppm=60, maxppm=85, path=path00carbon_P0, isproton=False)

    path00carbon_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_C_P1.csv"

    sample00_P1 = chloroform()
    sample00_P1.load_data_and_plot(minppm=60, maxppm=85, path=path00carbon_P1, isproton=False)

    path00carbon_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_C_P2.csv"

    sample00_P2 = chloroform()
    sample00_P2.load_data_and_plot(minppm=60, maxppm=85, path=path00carbon_P2, isproton=False)

    average_carbon_spectrum(sample00_P0, sample00_P1, sample00_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/00/00_C_avg.png",
                            storetitle="Averaged result in Carbon channel for 00 state.")
    return


def process_pseudo_state01_proton_data():
    path01proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_H_P0.csv"

    sample01_P0 = chloroform()
    sample01_P0.load_data_and_plot(path01proton_P0, isproton=True)

    path01proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_H_P1.csv"

    sample01_P1 = chloroform()
    sample01_P1.load_data_and_plot(path01proton_P1, isproton=True)

    path01proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_H_P2.csv"

    sample01_P2 = chloroform()
    sample01_P2.load_data_and_plot(path01proton_P2, isproton=True)

    average_proton_spectrum(sample01_P0, sample01_P1, sample01_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_H_avg.png",
                            storetitle="Averaged result in Proton channel")
    return


def process_pseudo_state01_carbon_data():
    path01carbon_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_C_P0.csv"

    sample01_P0 = chloroform()
    sample01_P0.load_data_and_plot(minppm=60, maxppm=85, path=path01carbon_P0, isproton=False)

    path01carbon_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_C_P1.csv"

    sample01_P1 = chloroform()
    sample01_P1.load_data_and_plot(minppm=60, maxppm=85, path=path01carbon_P1, isproton=False)

    path01carbon_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_C_P2.csv"

    sample01_P2 = chloroform()
    sample01_P2.load_data_and_plot(minppm=60, maxppm=85, path=path01carbon_P2, isproton=False)

    average_carbon_spectrum(sample01_P0, sample01_P1, sample01_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/01/01_C_avg.png",
                            storetitle="Averaged result in Carbon channel")
    return


def process_pseudo_state10_proton_data():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_H_P0.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_H_P1.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_H_P2.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_H_avg.png",
                            storetitle="Averaged result in Proton channel")
    return


def process_pseudo_state10_carbon_data():
    path10carbon_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_C_P0.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(minppm=60, maxppm=85, path=path10carbon_P0, isproton=False)

    path10carbon_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_C_P1.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(minppm=60, maxppm=85, path=path10carbon_P1, isproton=False)

    path10carbon_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_C_P2.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(minppm=60, maxppm=85, path=path10carbon_P2, isproton=False)

    average_carbon_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/10/10_C_avg.png",
                            storetitle="Averaged result in Carbon channel")
    return


def process_pseudo_state11_proton_data():
    path11proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_H_P0.csv"

    sample11_P0 = chloroform()
    sample11_P0.load_data_and_plot(path11proton_P0, isproton=True)

    path11proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_H_P1.csv"

    sample11_P1 = chloroform()
    sample11_P1.load_data_and_plot(path11proton_P1, isproton=True)

    path11proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_H_P2.csv"

    sample11_P2 = chloroform()
    sample11_P2.load_data_and_plot(path11proton_P2, isproton=True)

    average_proton_spectrum(sample11_P0, sample11_P1, sample11_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_H_avg.png",
                            storetitle="Averaged result in Proton channel")
    return


def process_pseudo_state11_carbon_data():
    path11carbon_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_C_P0.csv"

    sample11_P0 = chloroform()
    sample11_P0.load_data_and_plot(minppm=60, maxppm=85, path=path11carbon_P0, isproton=False)

    path11carbon_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_C_P1.csv"

    sample11_P1 = chloroform()
    sample11_P1.load_data_and_plot(minppm=60, maxppm=85, path=path11carbon_P1, isproton=False)

    path11carbon_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_C_P2.csv"

    sample11_P2 = chloroform()
    sample11_P2.load_data_and_plot(minppm=60, maxppm=85, path=path11carbon_P2, isproton=False)

    average_carbon_spectrum(sample11_P0, sample11_P1, sample11_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Pseudo/11/11_C_avg.png",
                            storetitle="Averaged result in Carbon channel")
    return


def process_DJ_00_data():
    path00proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/P0.csv"

    sample00_P0 = chloroform()
    sample00_P0.load_data_and_plot(path00proton_P0, isproton=True)

    path00proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/P1.csv"

    sample00_P1 = chloroform()
    sample00_P1.load_data_and_plot(path00proton_P1, isproton=True)

    path00proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/P2.csv"

    sample00_P2 = chloroform()
    sample00_P2.load_data_and_plot(path00proton_P2, isproton=True)

    average_proton_spectrum(sample00_P0, sample00_P1, sample00_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_H_avg.png",
                            storetitle="Averaged result in Proton channel for DJ. f1")
    return


def process_DJ_01_data():
    path01proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f2/P0.csv"

    sample01_P0 = chloroform()
    sample01_P0.load_data_and_plot(path01proton_P0, isproton=True)

    path01proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f2/P1.csv"

    sample01_P1 = chloroform()
    sample01_P1.load_data_and_plot(path01proton_P1, isproton=True)

    path01proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f2/P2.csv"

    sample01_P2 = chloroform()
    sample01_P2.load_data_and_plot(path01proton_P2, isproton=True)

    average_proton_spectrum(sample01_P0, sample01_P1, sample01_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f2/f2_H_avg.png",
                            storetitle="Averaged result in Proton channel for DJ. f2")

    return


def process_DJ_10_data():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f3/P0.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f3/P1.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f3/P2.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f3/f3_H_avg.png",
                            storetitle="Averaged result in Proton channel for DJ. f3")

    return


def process_DJ_11_data():
    path11proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f4/P0.csv"

    sample11_P0 = chloroform()
    sample11_P0.load_data_and_plot(path11proton_P0, isproton=True)

    path11proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f4/P1.csv"

    sample11_P1 = chloroform()
    sample11_P1.load_data_and_plot(path11proton_P1, isproton=True)

    path11proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f4/P2.csv"

    sample11_P2 = chloroform()
    sample11_P2.load_data_and_plot(path11proton_P2, isproton=True)

    average_proton_spectrum(sample11_P0, sample11_P1, sample11_P2)

    average_proton_spectrum(sample11_P0, sample11_P1, sample11_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f4/f4_H_avg.png",
                            storetitle="Averaged result in Proton channel for DJ. f4")

    return


def process_Grover_00_data():
    path00proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/00/P0.csv"

    sample00_P0 = chloroform()
    sample00_P0.load_data_and_plot(path00proton_P0, isproton=True)

    path00proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/00/P1.csv"

    sample00_P1 = chloroform()
    sample00_P1.load_data_and_plot(path00proton_P1, isproton=True)

    path00proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/00/P2.csv"

    sample00_P2 = chloroform()
    sample00_P2.load_data_and_plot(path00proton_P2, isproton=True)

    average_proton_spectrum(sample00_P0, sample00_P1, sample00_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/00/f00_H_avg.png",
                            storetitle="Averaged result in Proton channel for Grover. f00")
    return




def process_Grover_01_data():
    path01proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/01/P0.csv"

    sample01_P0 = chloroform()
    sample01_P0.load_data_and_plot(path01proton_P0, isproton=True)

    path01proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/01/P1.csv"

    sample01_P1 = chloroform()
    sample01_P1.load_data_and_plot(path01proton_P1, isproton=True)

    path01proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/01/P2.csv"

    sample01_P2 = chloroform()
    sample01_P2.load_data_and_plot(path01proton_P2, isproton=True)

    average_proton_spectrum(sample01_P0, sample01_P1, sample01_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/01/f01_H_avg.png",
                            storetitle="Averaged result in Proton channel for Grover. f01")
    return


def process_Grover_10_data():
    path10proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/10/P0.csv"

    sample10_P0 = chloroform()
    sample10_P0.load_data_and_plot(path10proton_P0, isproton=True)

    path10proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/10/P1.csv"

    sample10_P1 = chloroform()
    sample10_P1.load_data_and_plot(path10proton_P1, isproton=True)

    path10proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/10/P2.csv"

    sample10_P2 = chloroform()
    sample10_P2.load_data_and_plot(path10proton_P2, isproton=True)

    average_proton_spectrum(sample10_P0, sample10_P1, sample10_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/10/f10_H_avg.png",
                            storetitle="Averaged result in Proton channel for Grover. f10")
    return


def process_Grover_11_data():
    path11proton_P0 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/11/P0.csv"

    sample11_P0 = chloroform()
    sample11_P0.load_data_and_plot(path11proton_P0, isproton=True)

    path11proton_P1 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/11/P1.csv"

    sample11_P1 = chloroform()
    sample11_P1.load_data_and_plot(path11proton_P1, isproton=True)

    path11proton_P2 = "C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/11/P2.csv"

    sample11_P2 = chloroform()
    sample11_P2.load_data_and_plot(path11proton_P2, isproton=True)

    average_proton_spectrum(sample11_P0, sample11_P1, sample11_P2, store=True,
                            storePath="C:/Users/73747/PycharmProjects/NMRPulse/data/Grover/11/f11_H_avg.png",
                            storetitle="Averaged result in Proton channel for Grover. f11")
    return


from Simulator.algorithm import *

if __name__ == "__main__":
    # process_pseudo_state00_data()
    process_DJ_00_data()
    process_DJ_01_data()
    process_DJ_10_data()
    process_DJ_11_data()
    #process_pseudo_state00_proton_data()
    #process_pseudo_state01_proton_data()
    #process_pseudo_state10_proton_data()
    #process_pseudo_state11_proton_data()

    # process_pseudo_state00_carbon_data()
    # process_pseudo_state01_carbon_data()
    # process_pseudo_state10_carbon_data()
    # process_pseudo_state11_carbon_data()


    process_exact_CNOT_data_H_00()
    process_exact_CNOT_data_H_01()
    process_exact_CNOT_data_H_10()
    process_exact_CNOT_data_H_11()

    #process_approx_CNOT_data_H_00()
    #process_approx_CNOT_data_H_01()
    #process_approx_CNOT_data_H_10()
    #process_approx_CNOT_data_H_11()


    #process_approx_CNOT_data_C_10()
    #process_approx_CNOT_data_C_11()
    #process_approx_CNOT_data_H_10()
    #process_approx_CNOT_data_H_11()

    process_Grover_00_data()
    process_Grover_01_data()
    process_Grover_10_data()
    process_Grover_11_data()