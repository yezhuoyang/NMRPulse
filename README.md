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


## Exact CNOT gate

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

