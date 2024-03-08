# NMRPulse

Simulation of all NMR quantum gates and pulses.
I want to highlight the following property of my simulator:


__1. Close to reality simulation and visualization__

   The proton spectrum and carbon spectrum are shown in real ppm scale.

   
__2. Handy interface and many examples__
   
   User can easily understand and call the interface for simulation.
   Many examples are provided in example.py for users to learn and get familiar.

   
__3. Correctness guaranteed by testcases__
   
   All pulses sequence and evolution function has passed tests cases in test/ folder

   
__4. Print Spinsolve Pulse sequences for Pulse Debug__

   User can print the real pulse sequence directly for further check and debug in Spinsolve machine
   
__4. Four lines of code for Dj algorithm and Grover algorithm__

   User only need four lines of code to calculate the pulse level DJ algorithm and Grover algorithm




# Example of Chloroform samples


## Add a 45 degree pulse

### Code

```python
def pulse_length_change():
    '''
    Initialize the
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.4, 0, 0, 0],
                                    [0, 0.4, 0, 0],
                                    [0, 0, 0.1, 0],
                                    [0, 0, 0, 0.1]], dtype=complex))
    '''
    Add a single pulse on proton.
    '''
    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_fid_real(maxtime=0.1, store=True, path="Figure/45pulsesFID.png")
    '''
    Read the data signal in the frequency
    '''
    NMRsample.read_proton_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/45pulsespec.png")
```

### Result of proton FID
![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/45pulsesFID.png)

### Result of proton spectrum

![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/45pulsespec.png)





## Proton pulse length calibration

### Code

```python
def pulse_length_calib_proton():
    pulses_length_list = np.linspace(0, 1, 20)
    integra_list = []
    NMRsample = chloroform()
    for pulse in pulses_length_list:
        NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                        [0, 0.3, 0, 0],
                                        [0, 0, -0.3, 0],
                                        [0, 0, 0, -0.5]], dtype=complex))
        NMRsample.set_pulses([])
        '''
        The first 1/2 pi pulse is added to cancel 
        the sigmax in the measurement operator.
        '''
        NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90H, wH))
        '''
        This is the actual varying Ix pulse we add in the pulse
        length calibration 
        '''
        NMRsample.add_pulse(pulseSingle(0, pulse * pl90H, wH))
        NMRsample.evolve_all_pulse()
        NMRsample.read_proton_time()
        NMRsample.read_proton_spectrum(normalize=False)
        integra_list.append(NMRsample.integral_proton_spectrum_real())

    pulses_length_list = [2 * x * pl90H * 10 ** 6 for x in pulses_length_list]
    plt.scatter(pulses_length_list, integra_list, label="Integral value of proton spectrum")
    plt.axvline(x=pl90H*10**6, color="red", linestyle="--", label="Measured 90-x pulse for proton")
    plt.xlabel("Pulse length Time/ microsecond")
    plt.ylabel("Integral value")
    plt.legend(fontsize=8)
    plt.savefig("Figure/protoncalib.png")
    plt.show()
```

### Result of proton pulse length calibration
![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/protoncalib.png)




## Carbon pulse length calibration

### Code

```python
def pulse_length_calib_carbon():
    pulses_length_list = np.linspace(0, 1, 20)
    integra_list = []
    NMRsample = chloroform()
    for pulse in pulses_length_list:
        NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                        [0, 0.3, 0, 0],
                                        [0, 0, -0.3, 0],
                                        [0, 0, 0, -0.5]], dtype=complex))
        NMRsample.set_pulses([])
        '''
        The first 1/2 pi pulse is added to cancel 
        the sigmax in the measurement operator.
        '''
        NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
        '''
        This is the actual varying Ix pulse we add in the pulse
        length calibration 
        '''
        NMRsample.add_pulse(pulseSingle(0, pulse * pl90C, wC))
        NMRsample.evolve_all_pulse()
        NMRsample.read_carbon_time()
        NMRsample.read_carbon_spectrum(normalize=False)
        integra_list.append(NMRsample.integral_carbon_spectrum_real())

    pulses_length_list = [2 * x * pl90C * 10 ** 6 for x in pulses_length_list]
    plt.scatter(pulses_length_list, integra_list, label="Integral value of carbon spectrum")
    plt.axvline(x=pl90C*10**6, color="red", linestyle="--", label="Measured 90-x pulse for carbon")
    plt.xlabel("Pulse length Time/ microsecond")
    plt.ylabel("Integral value")
    plt.legend(fontsize=8)
    plt.savefig("Figure/carboncalib.png")
    plt.show()
```

### Result of carbon pulse length calibration
![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/carboncalib.png)






## Exact CNOT gate matrix

### Code

```python
def exact_CNOT():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))
    CNOTmatrix = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
    '''
    Directly evolve the density matrix by CNOT matrix
    '''
    NMRsample.evolve_density(CNOTmatrix)
    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    NMRsample.read_carbon_time()
    '''
    Read the spectrum
    '''
    NMRsample.read_proton_spectrum()
    NMRsample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/CNOTExactproton.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/CNOTExactcarbon.png")
```


### Result of proton spectrum after CNOT

![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/CNOTExactproton.png)

### Result of carbon spectrum after CNOT

![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/CNOTExactcarbon.png)



## Approximate CNOT gate controlled by proton

### Code

```python
def approx_CNOT():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))

    '''
    Add approximate CNOT pulse sequence
    (pi/2)Ix2---(2Iz1Iz2)---(pi/2)Iy2
    Recall that channel 0 for +x, 1 for +y, 2 for -x, 3 for -y
    '''
    NMRsample.add_pulse(pulseSingle(0, 0.5 * pl90C, wC))
    NMRsample.add_pulse(delayTime(0.5 / Jfreq))
    NMRsample.add_pulse(pulseSingle(1, 0.5 * pl90C, wC))
    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Print the unitary of all pulses:
    '''
    print(NMRsample.get_pulse_unitary())

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    NMRsample.read_carbon_time()
    '''
    Read the spectrum
    '''
    NMRsample.read_proton_spectrum()
    NMRsample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/CNOTapproxproton.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/CNOTapproxcarbon.png")
```


### Result of proton spectrum after approximate CNOT



![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/CNOTapproxproton.png)

### Result of carbon spectrum after approximate CNOT

![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/CNOTapproxcarbon.png)





## Exact CNOT gate controlled by proton

### Code

```python
def exact_CNOT_pulse_Hcontrol():
    '''
    Initialize the chloroform instance
    '''
    NMRsample = chloroform()
    '''
    Set the initial density matrix
    '''
    NMRsample.set_density(np.array([[0.5, 0, 0, 0],
                                    [0, 0.3, 0, 0],
                                    [0, 0, -0.3, 0],
                                    [0, 0, 0, -0.5]], dtype=complex))

    '''
    Add pulse sequence for approximate h gate on carbon
    '''
    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))

    '''
    Add pulse sequence for exact CZ gate
    '''
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90H, wH))
    NMRsample.add_pulse(pulseSingle(2, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(1, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 / 2 * pl90C, wC))
    NMRsample.add_pulse(delayTime((4 - 0.5) / Jfreq))

    '''
    Add pulse sequence for approximate h gate on carbon
    '''

    NMRsample.add_pulse(pulseSingle(1, 1 / 4 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(0, 1 * pl90C, wC))
    NMRsample.add_pulse(pulseSingle(3, 1 / 4 * pl90C, wC))

    '''
    Evolve the density matrix with all pulses
    '''
    NMRsample.evolve_all_pulse()
    '''
    Print the unitary of all pulses:
    '''
    print(NMRsample.get_pulse_unitary())

    '''
    Read the data signal in the time domain
    '''
    NMRsample.read_proton_time()
    NMRsample.read_carbon_time()
    '''
    Read the spectrum
    '''
    NMRsample.read_proton_spectrum()
    NMRsample.read_carbon_spectrum()
    '''
    Simulate what is shown on the screen
    '''
    NMRsample.show_proton_spectrum_real(-5, 15, store=True,
                                        path="Figure/CNOTExactproton-Hcontrol.png")

    NMRsample.show_carbon_spectrum_real(74, 80, store=True,
                                        path="Figure/CNOTExactcarbon-Hcontrol.png")
```


### Result of proton spectrum after approximate CNOT


![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/CNOTExactproton-Hcontrol.png)

### Result of carbon spectrum after approximate CNOT

![alt text](https://github.com/yezhuoyang/NMRPulse/blob/main/Figure/CNOTExactcarbon-Hcontrol.png)



## Deutsch jozsa algorithm

### Code


```python
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

    DJ.show_spectrum("Figure/DJP0f{}{}".format(uf[0], uf[1]),
                     title="Result of DJ algorithm after P0 for f{}{}".format(uf[0], uf[1]))

    density0 = DJ.get_final_density()

    '''
    Reinitialize the sample, add a P1 permutation:
    '''
    DJ.init_sample()
    DJ.set_prem_value(1)

    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJP1f{}{}".format(uf[0], uf[1]),
                     title="Result of DJ algorithm after P1 for f{}{}".format(uf[0], uf[1]))

    density1 = DJ.get_final_density()

    '''
    Reinitialize the sample, add a P2 permutation:
    '''
    DJ.init_sample()
    DJ.set_prem_value(2)

    DJ.construct_pulse()
    DJ.calculate_result_pulse()

    DJ.show_spectrum("Figure/DJP2f{}{}".format(uf[0], uf[1]),
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
                                            path="Figure/DJproton%d%d.png" % (uf[0], uf[1]))

    pseudo_sample.show_carbon_spectrum_real(74, 80, store=True,
                                            path="Figure/DJcarbon%d%d.png" % (uf[0], uf[1]))
```






## Grover's algorithm

### Code


```python
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

    grover.show_spectrum("Figure/GroverP0f{}{}".format(db[0], db[1]),
                         title="Result of Grover algorithm after P0 for f{}{}".format(db[0], db[1]))

    density0 = grover.get_final_density()

    '''
    Reinitialize the sample, add a P1 permutation:
    '''
    grover.init_sample()
    grover.set_prem_value(1)

    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/GroverP1f{}{}".format(db[0], db[1]),
                         title="Result of Grover algorithm after P1 for f{}{}".format(db[0], db[1]))

    density1 = grover.get_final_density()

    '''
    Reinitialize the sample, add a P2 permutation:
    '''
    grover.init_sample()
    grover.set_prem_value(2)

    grover.construct_pulse()
    grover.calculate_result_pulse()

    grover.show_spectrum("Figure/GroverP2f{}{}".format(db[0], db[1]),
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
                                            path="Figure/Groverproton%d%d.png" % (db[0], db[1]))

    pseudo_sample.show_carbon_spectrum_real(74, 80, store=True,
                                            path="Figure/Grovercarbon%d%d.png" % (db[0], db[1]))
```



## Print pulse sequences:


```python
def DJ_print_pulse(uf):
    DJ = Djalgorithm()
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


if __name__ == "__main__":
    DJ_print_pulse([0,0])
```


