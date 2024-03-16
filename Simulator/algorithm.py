from typing import List

import qiskit
from qiskit_aer import AerSimulator

from Simulator.NMRdensity import *
from Simulator.Pulses import *


class NMRalgorithm:
    def __init__(self):
        self.circuit = qiskit.QuantumCircuit(2, 1)
        self.UF = []
        self.computed = False
        self.balance = False
        self.simulator = AerSimulator()
        self.NMRsample = chloroform()
        '''
        Initialize the state to be thermal equilibrium state
        '''
        self.NMRsample.set_thermal_equilibrium()
        self.approximate = False

        self.perm_index_value = 0

    def init_sample(self):
        self.NMRsample.set_thermal_equilibrium()
        self.NMRsample.set_pulses([])

    def set_prem_value(self, value):
        self.perm_index_value = value

    def add_p1_perm_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="P1"))
        self.NMRsample.add_p1_perm_pulse()
        self.NMRsample.add_pulse(BarrierPulse(name="End P1"))

    def add_p2_perm_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="P2"))
        self.NMRsample.add_p2_perm_pulse()
        self.NMRsample.add_pulse(BarrierPulse(name="P2", endpulse=True))

    def add_X_gate_first_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="X1"))
        self.NMRsample.add_X_gate_first_pulse()
        self.NMRsample.add_pulse(BarrierPulse(name="X1", endpulse=True))

    def add_X_gate_second_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="X2"))
        self.NMRsample.add_X_gate_second_pulse()
        self.NMRsample.add_pulse(BarrierPulse(name="X2", endpulse=True))

    def add_H_gate_first_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="H1"))
        self.NMRsample.add_H_gate_first_pulse(approximate=self.approximate)
        self.NMRsample.add_pulse(BarrierPulse(name="H1", endpulse=True))

    def add_H_gate_second_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="H2"))
        self.NMRsample.add_H_gate_second_pulse(approximate=self.approximate)
        self.NMRsample.add_pulse(BarrierPulse(name="H2", endpulse=True))

    def add_CZ_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="CZ(H,C)"))
        self.NMRsample.add_CZ_pulse(approximate=self.approximate, Hcontrol=True)
        self.NMRsample.add_pulse(BarrierPulse(name="CZ(H,C)", endpulse=True))

    def add_CNOT_pulse(self):
        self.NMRsample.add_pulse(BarrierPulse(name="CNOT(H,C)"))
        self.NMRsample.add_CNOT_pulse(approximate=self.approximate, Hcontrol=True)
        self.NMRsample.add_pulse(BarrierPulse(name="CNOT(H,C)", endpulse=True))

    def plot_measure_all(self):
        raise NotImplementedError

    def get_final_density(self):
        return self.NMRsample.get_density()

    def set_approxmiate(self, value):
        self.approximate = value

    def set_function(self,
                     *params):
        raise NotImplementedError

    def construct_pulse(self):
        raise NotImplementedError

    def construct_circuit(self):
        raise NotImplementedError

    def calculate_result(self):
        raise NotImplementedError

    def print_pulses(self):
        self.NMRsample.print_pulses()

    def show_spectrum(self, path: str, title: str):
        '''
        Read the data signal in the time domain
        '''
        self.NMRsample.read_proton_time()
        self.NMRsample.read_carbon_time()
        '''
        Read the spectrum
        '''
        self.NMRsample.read_proton_spectrum()
        self.NMRsample.read_carbon_spectrum()
        '''
        Simulate what is shown on the screen
        '''
        self.NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                                 path=path + "proton.png", title=title)

        self.NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                                 path=path + "carbon.png", title=title)


class Djalgorithm(NMRalgorithm):

    def __init__(self):
        super().__init__()

    def construct_pulse(self):
        if self.perm_index_value == 1:
            self.add_p1_perm_pulse()
        elif self.perm_index_value == 2:
            self.add_p2_perm_pulse()

        # self.circuit.x(1)
        self.add_X_gate_second_pulse()
        # self.circuit.h(list(range(0, 2)))
        self.add_H_gate_first_pulse()
        self.add_H_gate_second_pulse()
        self.compile_uf_pulse()
        self.add_H_gate_first_pulse()
        self.add_H_gate_second_pulse()

    def calculate_result_circuit(self):
        compiled_circuit = qiskit.transpile(self.circuit, self.simulator)

        # Execute the circuit on the aer simulator
        job = self.simulator.run(compiled_circuit, shots=1)

        # Grab results from the job
        result = job.result()
        # print(result)
        # Returns counts
        counts = result.get_counts(compiled_circuit)
        result = list(counts.keys())[0]
        if result == '0':
            self.balance = False
            print("The function is constant")
        else:
            self.balance = True
            print("The function is balanced")

    def calculate_result_pulse(self):
        self.NMRsample.evolve_all_pulse()
        '''
        Read the data signal in the time domain
        '''
        self.NMRsample.read_proton_time()
        self.NMRsample.read_carbon_time()
        '''
        Read the spectrum
        '''
        self.NMRsample.read_proton_spectrum()
        self.NMRsample.read_carbon_spectrum()
        '''
        Simulate what is shown on the screen
        '''
        self.NMRsample.show_proton_spectrum_real(-5, 15, store=False)

        self.NMRsample.show_carbon_spectrum_real(74, 80, store=False)

    def construct_circuit(self):
        '''
        The first layer of Hadmard
        '''
        self.circuit.x(1)
        self.circuit.h(list(range(0, 2)))
        self.compile_uf_circuit()
        self.circuit.h(list(range(0, 2)))
        self.circuit.measure(list(range(0, 2 - 1)), list(range(0, 2 - 1)))

    def set_function(self,
                     uf: List):
        assert len(uf) == 2
        self.UF = uf

    def compile_uf_circuit(self):
        for i in range(0, 1 << (2 - 1)):
            if self.UF[i] == 1:
                if i == 0:
                    self.circuit.x(i)
                self.circuit.cx(0, 1)
                if i == 0:
                    self.circuit.x(i)
        return

    def compile_uf_pulse(self):
        for i in range(0, 1 << (2 - 1)):
            if self.UF[i] == 1:
                if i == 0:
                    self.add_H_gate_first_pulse()
                self.add_CNOT_pulse()
                if i == 0:
                    self.add_H_gate_first_pulse()
        return

    def plot_measure_all(self):
        self.circuit.clear()
        self.circuit.x(1)
        self.circuit.h(list(range(0, 2)))
        self.compile_uf_circuit()
        self.circuit.h(list(range(0, 2)))
        self.circuit.measure_all()

        # Use Aer's qasm_simulator
        simulator = self.simulator

        # Execute the circuit on the qasm simulator
        job = qiskit.execute(self.circuit, simulator, shots=5000)

        # Grab results from the job
        result = job.result()

        # Returns counts with suffix in keys
        counts_with_suffix = result.get_counts(self.circuit)

        # Remove suffix and sum counts if necessary
        counts = {}
        for key_with_suffix, count in counts_with_suffix.items():
            key = key_with_suffix.split()[0]  # Assumes suffix is after a space and we only want the first part
            if key in counts:
                counts[key] += count
            else:
                counts[key] = count

        # Calculate probabilities
        probabilities = {state: count / 5000 for state, count in counts.items()}

        print("Prob!")
        print(probabilities)

        # Ensure all possible outcomes are present in the probabilities dictionary
        for state in ['00', '01', '10', '11']:
            if state not in probabilities:
                probabilities[state] = 0

        # Sorting states to ensure consistent order
        sorted_states = sorted(probabilities.keys())
        sorted_probs = [probabilities[state] for state in sorted_states]

        # Plotting using matplotlib
        plt.bar(sorted_states, sorted_probs)
        plt.xlabel('State')
        plt.ylabel('Probability')
        plt.title('Measurement result of DJ for f{}{}'.format(self.UF[0], self.UF[1]))
        plt.savefig("Figure/DJ/{}{}/DJSimu{}{}".format(self.UF[0], self.UF[1], self.UF[0], self.UF[1]))
        plt.show()


class Grover(NMRalgorithm):

    def __init__(self):
        super().__init__()
        self.computed = False
        self.circuit = qiskit.QuantumCircuit(2, 2)
        self._simulator = AerSimulator()
        '''
        How many time we may need to call the grover operator,
        initially set to be 1
        '''
        self.grover_step = 1
        self.max_step = 1
        self._database = []
        self._solution = -1
        self._succeed = False

    def set_grover_step(self, value):
        self.grover_step = value

    def construct_circuit(self):
        self.circuit.x(1)
        self.circuit.h(list(range(0, 2)))
        '''
        Construct grover circuit many times
        '''
        for i in range(self.grover_step):
            self.circuit.barrier([0, 1])
            self.construct_grover_op_circuit()
        self.circuit.x(0)
        self.circuit.measure(list(range(0, 2)), list(range(0, 2)))
        return

    def construct_pulse(self):
        if self.perm_index_value == 1:
            self.add_p1_perm_pulse()
        elif self.perm_index_value == 2:
            self.add_p2_perm_pulse()
        self.add_X_gate_second_pulse()
        self.add_H_gate_first_pulse()
        self.add_H_gate_second_pulse()
        for i in range(self.grover_step):
            self.construct_grover_op_pulse()
        self.add_X_gate_first_pulse()
        return

    def calculate_result_circuit(self):

        self.construct_circuit()
        compiled_circuit = qiskit.transpile(self.circuit, self._simulator)
        # Execute the circuit on the aer simulator
        job = self._simulator.run(compiled_circuit, shots=1)

        # Grab results from the job
        result = job.result()

        print(result)

        # Returns counts
        counts = result.get_counts(compiled_circuit)
        print(counts)
        result = list(counts.keys())[0]
        print("Measure:")
        print(result)

    def calculate_result_pulse(self):
        self.NMRsample.evolve_all_pulse()
        '''
        Read the data signal in the time domain
        '''
        self.NMRsample.read_proton_time()
        self.NMRsample.read_carbon_time()
        '''
        Read the spectrum
        '''
        self.NMRsample.read_proton_spectrum()
        self.NMRsample.read_carbon_spectrum()
        '''
        Simulate what is shown on the screen
        '''
        self.NMRsample.show_proton_spectrum_real(-5, 15, store=False)

        self.NMRsample.show_carbon_spectrum_real(74, 80, store=False)

    def construct_zf_circuit(self):
        for i in range(0, len(self._database)):
            if self._database[i] == 1:
                if i == 0:
                    self.circuit.x(0)
                    self.circuit.x(1)
                    self.circuit.cz(0, 1)
                    self.circuit.x(0)
                    self.circuit.x(1)
                elif i == 1:
                    self.circuit.x(0)
                    self.circuit.cz(0, 1)
                    self.circuit.x(0)
                elif i == 2:
                    self.circuit.x(1)
                    self.circuit.cz(0, 1)
                    self.circuit.x(1)
                elif i == 3:
                    self.circuit.cz(0, 1)
        return

    def construct_zf_pulse(self):
        for i in range(0, len(self._database)):
            if self._database[i] == 1:
                if i == 0:
                    self.add_X_gate_first_pulse()
                    self.add_X_gate_second_pulse()
                    self.add_CZ_pulse()
                    self.add_X_gate_first_pulse()
                    self.add_X_gate_second_pulse()
                elif i == 1:
                    self.add_X_gate_first_pulse()
                    self.add_CZ_pulse()
                    self.add_X_gate_first_pulse()
                elif i == 2:
                    self.add_X_gate_second_pulse()
                    self.add_CZ_pulse()
                    self.add_X_gate_second_pulse()
                elif i == 3:
                    self.add_CZ_pulse()
        return

    def construct_zo_circuit(self):
        self.circuit.x(0)
        self.circuit.x(1)
        self.circuit.cz(0, 1)
        self.circuit.x(0)
        self.circuit.x(1)

    def construct_zo_pulse(self):
        self.add_X_gate_first_pulse()
        self.add_X_gate_second_pulse()
        self.add_CZ_pulse()
        self.add_X_gate_first_pulse()
        self.add_X_gate_second_pulse()

    def construct_grover_op_circuit(self):
        self.circuit.barrier([0, 1])
        self.construct_zf_circuit()
        self.circuit.barrier([0, 1])
        self.circuit.h(list(range(0, 2)))
        self.circuit.barrier([0, 1])
        self.construct_zo_circuit()
        self.circuit.barrier([0, 1])
        self.circuit.h(list(range(0, 2)))
        self.circuit.barrier([0, 1])
        return

    def construct_grover_op_pulse(self):
        self.construct_zf_pulse()
        self.add_H_gate_first_pulse()
        self.add_H_gate_second_pulse()
        self.construct_zo_pulse()
        self.add_H_gate_first_pulse()
        self.add_H_gate_second_pulse()
        return

    def set_function(self,
                     database: List):
        # assert len(database) == 2
        self._database = database

    def plot_measure_all(self):
        self.circuit.clear()
        self.circuit.x(1)
        self.circuit.barrier([0, 1])
        self.circuit.h(list(range(0, 2)))
        self.circuit.barrier([0, 1])
        '''
        Construct grover circuit many times
        '''
        for i in range(self.grover_step):
            self.construct_grover_op_circuit()

        self.circuit.x(0)
        self.circuit.measure_all()

        # Use Aer's qasm_simulator
        simulator = self.simulator

        # Execute the circuit on the qasm simulator
        job = qiskit.execute(self.circuit, simulator, shots=5000)

        # Grab results from the job
        result = job.result()

        # Returns counts with suffix in keys
        counts_with_suffix = result.get_counts(self.circuit)

        # Remove suffix and sum counts if necessary
        counts = {}
        for key_with_suffix, count in counts_with_suffix.items():
            key = key_with_suffix.split()[0]  # Assumes suffix is after a space and we only want the first part
            if key in counts:
                counts[key] += count
            else:
                counts[key] = count

        # Calculate probabilities
        probabilities = {state: count / 5000 for state, count in counts.items()}

        print("Prob!")
        print(probabilities)

        # Ensure all possible outcomes are present in the probabilities dictionary
        for state in ['00', '01', '10', '11']:
            if state not in probabilities:
                probabilities[state] = 0

        # Sorting states to ensure consistent order
        sorted_states = sorted(probabilities.keys())
        sorted_probs = [probabilities[state] for state in sorted_states]

        index = -1
        bitstr = ""
        for i in range(0, 4):
            if self._database[i] == 1:
                if i == 0:
                    bitstr = "00"
                elif i == 1:
                    bitstr = "01"
                elif i == 2:
                    bitstr = "10"
                elif i == 3:
                    bitstr = "11"

        # Plotting using matplotlib
        plt.bar(sorted_states, sorted_probs)
        plt.xlabel('State')
        plt.ylabel('Probability')
        plt.title('Measurement result of Grover for f{}(Grover Step:{})'.format(bitstr,
                                                                                self.grover_step))

        print("bitstr")
        print(bitstr)
        plt.savefig("Figure/Grover/{}/GroverSimu{}".format(bitstr, bitstr))
        plt.show()


def permute_DJ(uf):
    DJ = Djalgorithm()
    DJ.set_function(uf)

    DJ.construct_circuit()
    DJ.calculate_result_circuit()

    DJ.plot_measure_all()

    '''
    First, calculate the result without permutation
    '''
    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJ/{}{}/DJP0f{}{}".format(uf[0], uf[1], uf[0], uf[1]),
                     title="Result of DJ algorithm after P0 for f{}{}".format(uf[0], uf[1]))

    density0 = DJ.get_final_density()

    '''
    Reinitialize the sample, add a P1 permutation:
    '''
    DJ.init_sample()
    DJ.set_prem_value(1)

    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJ/{}{}/DJP1f{}{}".format(uf[0], uf[1], uf[0], uf[1]),
                     title="Result of DJ algorithm after P1 for f{}{}".format(uf[0], uf[1]))

    density1 = DJ.get_final_density()

    '''
    Reinitialize the sample, add a P2 permutation:
    '''
    DJ.init_sample()
    DJ.set_prem_value(2)

    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJ/{}{}/DJP2f{}{}".format(uf[0], uf[1], uf[0], uf[1]),
                     title="Result of DJ algorithm after P2 for f{}{}".format(uf[0], uf[1]))

    density2 = DJ.get_final_density()

    final_density = 1 / 3 * (density0 + density1 + density2)

    pseudo_sample = chloroform()
    pseudo_sample.set_density(final_density)

    pseudo_sample.read_proton_time()
    pseudo_sample.read_carbon_time()
    '''
    Read the spectrum
    '''
    pseudo_sample.read_proton_spectrum()
    pseudo_sample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    pseudo_sample.show_proton_spectrum_real(-5, 15, store=True,
                                            path="Figure/DJ/{}{}/DJproton{}{}.png".format(uf[0], uf[1], uf[0], uf[1]))

    pseudo_sample.show_carbon_spectrum_real(74, 80, store=True,
                                            path="Figure/DJ/{}{}/DJcarbon{}{}.png".format(uf[0], uf[1], uf[0], uf[1]))


def permute_grover(db):
    none_zero_index = db[0] * 2 + db[1]
    print(none_zero_index)

    database = [0, 0, 0, 0]
    database[none_zero_index] = 1

    grover = Grover()
    grover.set_grover_step(1)
    grover.set_function(database)

    grover.construct_circuit()
    grover.calculate_result_circuit()

    grover.plot_measure_all()

    '''
    First, calculate the result without permutation
    '''
    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/Grover/{}{}/GroverP0f{}{}".format(db[0], db[1], db[0], db[1]),
                         title="Result of Grover algorithm after P0 for f{}{}".format(db[0], db[1]))

    density0 = grover.get_final_density()

    '''
    Reinitialize the sample, add a P1 permutation:
    '''
    grover.init_sample()
    grover.set_prem_value(1)

    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/Grover/{}{}/GroverP1f{}{}".format(db[0], db[1], db[0], db[1]),
                         title="Result of Grover algorithm after P1 for f{}{}".format(db[0], db[1]))

    density1 = grover.get_final_density()

    '''
    Reinitialize the sample, add a P2 permutation:
    '''
    grover.init_sample()
    grover.set_prem_value(2)

    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/Grover/{}{}/GroverP2f{}{}".format(db[0], db[1], db[0], db[1]),
                         title="Result of Grover algorithm after P2 for f{}{}".format(db[0], db[1]))

    density2 = grover.get_final_density()

    final_density = 1 / 3 * (density0 + density1 + density2)

    pseudo_sample = chloroform()
    pseudo_sample.set_density(final_density)

    pseudo_sample.read_proton_time()
    pseudo_sample.read_carbon_time()
    '''
    Read the spectrum
    '''
    pseudo_sample.read_proton_spectrum()
    pseudo_sample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    pseudo_sample.show_proton_spectrum_real(-5, 15, store=True,
                                            path="Figure/Grover/{}{}/Groverproton{}{}.png".format(
                                                db[0], db[1], db[0], db[1]))

    pseudo_sample.show_carbon_spectrum_real(74, 80, store=True,
                                            path="Figure/Grover/{}{}/Grovercarbon{}{}.png".format(
                                                db[0], db[1], db[0], db[1]))


def DJ_print_pulse(uf):
    DJ = Djalgorithm()
    DJ.set_prem_value(1)
    '''
    Initialize the input function
    f1:uf=[0,0]
    f2:uf=[0,1]
    f3:uf=[1,0]
    f4:uf=[1,1]
    '''
    DJ.set_function(uf)
    DJ.construct_pulse()
    '''
    Print the real Spinsolve pulses
    '''
    DJ.print_pulses()


def Grover_print_pulse(uf):
    GV = Grover()
    '''
    Initialize the input function
    f1:uf=[1,0,0,0]
    f2:uf=[0,1,0,0]
    f3:uf=[0,0,1,0]
    f4:uf=[0,0,0,1]
    '''
    GV.set_function(uf)
    GV.construct_pulse()
    '''
    Print the real Spinsolve pulses
    '''
    GV.print_pulses()


if __name__ == "__main__":
    # DJ_print_pulse([0, 0])

    # permute_DJ([1,1])
    # permute_grover([1, 0])
    Grover_print_pulse([1, 0, 0, 0])
