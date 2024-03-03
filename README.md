# NMRPulse

Simulation of all NMR quantum gates and pulses.




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



