import qiskit
import functools
from qiskit_aer import AerSimulator
from NMRdensity import *
from Pulses import *
from typing import List, Union, Any
from qiskit.visualization import plot_histogram


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
        self.NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
        self.NMRsample.add_pulse(delayTime(0.5 / Jfreq))
        self.NMRsample.add_pulse(pulseTwo(3, 0.5 * pl90C, wC, 0, 0.5 * pl90H, wH))
        self.NMRsample.add_pulse(delayTime(0.5 / Jfreq))
        self.NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90H, wH))

    def add_p2_perm_pulse(self):
        self.NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
        self.NMRsample.add_pulse(delayTime(0.5 / Jfreq))
        self.NMRsample.add_pulse(pulseTwo(3, 0.5 * pl90H, wH, 0, 0.5 * pl90C, wC))
        self.NMRsample.add_pulse(delayTime(0.5 / Jfreq))
        self.NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))

    def add_X_gate_first_pulse(self):
        self.NMRsample.add_pulse(pulseSingle(0, pl90H, wH))

    def add_X_gate_second_pulse(self):
        self.NMRsample.add_pulse(pulseSingle(0, pl90C, wC))

    def add_H_gate_first_pulse(self):
        self.NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
        self.NMRsample.add_pulse(pulseSingle(0, 1 * pl90H, wH))
        self.NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))

    def add_H_gate_second_pulse(self):
        self.NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
        self.NMRsample.add_pulse(pulseSingle(0, 1 * pl90C, wC))
        self.NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))

    def add_CZ_pulse(self):
        if not self.approximate:
            '''
            Add pulse sequence for exact CZ gate
            (pi/2)Iz1---(pi/2)Iz2---(-2)Iz1Iz2
            Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
            '''
            self.NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))

            self.NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
            self.NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
            self.NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))

            '''
            The pulse of (-2)Iz1Iz2. The angle theta is actually (-\pi/2). However, since we cannot rotate 
            a minus angle, we should plus another (4\pi)
            '''
            self.NMRsample.add_pulse(delayTime((4 - 0.5) / Jfreq))
        else:
            pass

    def add_CNOT_pulse(self):
        if not self.approximate:
            '''
            Add pulse sequence for approximate h gate on carbon
            '''
            self.NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(0, 1 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))

            '''
            Add pulse sequence for exact CZ gate
            '''
            self.NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
            self.NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
            self.NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))
            self.NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))
            self.NMRsample.add_pulse(delayTime((4 - 0.5) / Jfreq))

            '''
            Add pulse sequence for approximate h gate on carbon
            '''

            self.NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(0, 1 * pl90H, wH))
            self.NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90H, wH))
        else:
            self.NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
            self.NMRsample.add_pulse(delayTime(0.5 / Jfreq))
            self.NMRsample.add_pulse(pulseSingle(3, 0.5 * pl90C, wC))

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
        plt.savefig("Figure/DJSimu{}{}".format(self.UF[0], self.UF[1]))
        plt.show()


class Grover(NMRalgorithm):

    def __init__(self):
        super().__init__()
        self.computed = False
        self.circuit = qiskit.QuantumCircuit(2, 1)
        self._simulator = AerSimulator()
        '''
        How many time we may need to call the grover operator,
        initially set to be 1
        '''
        self.grover_step = 2
        self.max_step = 2
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
            self.construct_grover_op_circuit()
        self.circuit.measure(list(range(0, 1)), list(range(0, 1)))
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
        return

    def calculate_result_circuit(self):
        # print(self._database)
        while not self.computed:
            self.construct_circuit()
            compiled_circuit = qiskit.transpile(self.circuit, self._simulator)
            # Execute the circuit on the aer simulator
            job = self._simulator.run(compiled_circuit, shots=1)

            # Grab results from the job
            result = job.result()
            # Returns counts
            counts = result.get_counts(compiled_circuit)
            print(counts)
            result = list(counts.keys())[0]

            self.computed = False
            result = int(result[::-1], 2)
            if self._database[result] == 1:
                self.computed = True
                self._succeed = True
                self._solution = result
                print(f"Found a solution {result}, Grover step={self.grover_step}")
                return
            else:
                self.grover_step = self.grover_step + 1
                self.circuit.clear()
                if self.grover_step >= self.max_step:
                    break
        if not self.computed:
            print("No solution!")

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
                if i == 1:
                    self.circuit.cz(0, 1)
                elif i == 0:
                    self.circuit.x(0)
                    self.circuit.cz(0, 1)
                    self.circuit.x(0)
        return

    def construct_zf_pulse(self):
        for i in range(0, len(self._database)):
            if self._database[i] == 1:
                if i == 1:
                    self.add_CZ_pulse()
                elif i == 0:
                    self.add_X_gate_first_pulse()
                    self.add_CZ_pulse()
                    self.add_X_gate_first_pulse()
        return

    def construct_zo_circuit(self):
        self.circuit.cz(0, 1)

    def construct_zo_pulse(self):
        self.add_CZ_pulse()

    def construct_grover_op_circuit(self):
        self.construct_zf_circuit()
        self.circuit.h(list(range(0, 2)))
        self.construct_zo_circuit()
        self.circuit.h(list(range(0, 2)))
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
        assert len(database) == 2
        self._database = database

    def plot_measure_all(self):
        self.circuit.clear()
        self.circuit.x(1)
        self.circuit.h(list(range(0, 2)))
        '''
        Construct grover circuit many times
        '''
        for i in range(self.grover_step):
            self.construct_grover_op_circuit()

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
        plt.title('Measurement result of Grover for f{}{}(Grover Step:{})'.format(self._database[0], self._database[1],
                                                                                  self.grover_step))
        plt.savefig("Figure/GroverSimu{}{}".format(self._database[0], self._database[1]))
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

    DJ.show_spectrum("Figure/DJP0f{}{}".format(uf[0], uf[1]), title="Result of DJ algorithm after P0 for f{}{}".format(uf[0], uf[1]))

    density0 = DJ.get_final_density()

    '''
    Reinitialize the sample, add a P1 permutation:
    '''
    DJ.init_sample()
    DJ.set_prem_value(1)

    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJP1f{}{}".format(uf[0], uf[1]), title="Result of DJ algorithm after P1 for f{}{}".format(uf[0], uf[1]))

    density1 = DJ.get_final_density()

    '''
    Reinitialize the sample, add a P2 permutation:
    '''
    DJ.init_sample()
    DJ.set_prem_value(2)

    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJP2f{}{}".format(uf[0], uf[1]), title="Result of DJ algorithm after P2 for f{}{}".format(uf[0], uf[1]))


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
                                            path="Figure/DJproton%d%d.png" % (uf[0], uf[1]))

    pseudo_sample.show_carbon_spectrum_real(74, 80, store=True,
                                            path="Figure/DJcarbon%d%d.png" % (uf[0], uf[1]))


def permute_grover(db):
    grover = Grover()
    grover.set_function(db)

    grover.construct_circuit()
    grover.calculate_result_circuit()

    '''
    First, calculate the result without permutation
    '''
    grover.construct_pulse()
    grover.calculate_result_pulse()


    grover.show_spectrum("Figure/GroverP0f{}{}".format(db[0], db[1]), title="Result of Grover algorithm after P0 for f{}{}".format(db[0], db[1]))

    density0 = grover.get_final_density()

    '''
    Reinitialize the sample, add a P1 permutation:
    '''
    grover.init_sample()
    grover.set_prem_value(1)

    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/GroverP1f{}{}".format(db[0], db[1]), title="Result of Grover algorithm after P1 for f{}{}".format(db[0], db[1]))

    density1 = grover.get_final_density()

    '''
    Reinitialize the sample, add a P2 permutation:
    '''
    grover.init_sample()
    grover.set_prem_value(2)

    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/GroverP2f{}{}".format(db[0], db[1]), title="Result of Grover algorithm after P2 for f{}{}".format(db[0], db[1]))

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
                                            path="Figure/Groverproton%d%d.png" % (db[0], db[1]))

    pseudo_sample.show_carbon_spectrum_real(74, 80, store=True,
                                            path="Figure/Grovercarbon%d%d.png" % (db[0], db[1]))


def DJ_print_pulse(uf):
    DJ = Djalgorithm()
    DJ.set_function(uf)

    DJ.construct_pulse()

    DJ.print_pulses()


if __name__ == "__main__":
    #permute_DJ([0, 0])
    #permute_DJ([0, 1])
    #permute_DJ([1, 0])
    #permute_DJ([1, 1])
    '''
    grover = Grover()

    grover.set_grover_step(1)

    grover.set_function([0, 0])

    grover.plot_measure_all()
    '''
    #permute_grover([0,0])
    #permute_DJ([0, 1])
    permute_DJ([1, 0])
    #permute_DJ([1, 1])



    '''
    DJ = Djalgorithm()
    DJ.set_function([1, 0])

    DJ.construct_circuit()
    DJ.calculate_result_circuit()

    DJ.plot_measure_all()
    '''
