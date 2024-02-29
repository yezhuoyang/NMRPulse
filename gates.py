import numpy as np
from Pulses import pulseTwo, pulseSingle, delayTime




CNOT_approx_HC = []
CNOT_approx_CH = []
Hadamard_approx = []
Hadamard_approx_inverse=[]



def Xgate():
    return


def Hadamard():
    return


def ControlledZ():
    return


def CNOT():
    return


if __name__ == "__main__":
    qc = QubitCircuit(2)
    qc.add_gate("X", targets=0)
    processor = LinearSpinChain(2)
    processor.load_circuit(qc)
    fig, axis = processor.plot_pulses()
    fig.show()
