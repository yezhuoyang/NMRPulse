import numpy as np
from qutip import sigmaz
from qutip.qip.device import Processor
from qutip.qip.pulse import Pulse
from qutip.qip.device import LinearSpinChain
from qutip.qip.circuit import QubitCircuit





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