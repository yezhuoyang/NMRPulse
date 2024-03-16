from Simulator.NMRdensity import *
from Simulator.Pulses import *
import matplotlib.pyplot as plt


def plot_data():
    NMRsample = chloroform()
    NMRsample.load_data_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")
    # NMRsample.read_and_plot("C:/Users/73747/PycharmProjects/NMRPulse/data/DJ/f1/f1_C_P0.csv")


def process_CNOT_data():
    return


def process_approx_CNOT_data():
    return


def process_pseudo_state_data():
    return


def process_DJ_data():
    return


def process_Grover_data():
    return


from Simulator.algorithm import *

if __name__ == "__main__":
    plot_data()
