# blazar_model
Attempt to create a time-dependent radiative model for AGN

# SSC_test guide 
This is a small guide to use the macros contained in this repository

## what you need
To use the code you need python (code is written with python 2.7), the following packages are also required
These are the ones I'm actually using on Red Hat Enterprise Linux 6.8
  numpy==1.10.2
  scipy==0.16.0
  astropy==1.1.2
  matplotlib==1.4.3
  naima==0.8
last one - naima - is the most important.
Naima and all of its dependencies can be installed in an Anaconda distribution through the Astropy conda channel:
    $ conda config --add channels astropy
    $ conda install naima
For more info http://naima.readthedocs.io/en/latest/index.html.


## structure
The main reference for all the calculations is  [arxiv of the paper](https://arxiv.org/abs/astro-ph/9810263) - pdf version available in the references.  
We are using the great python package [naima](http://naima.readthedocs.io/en/latest/index.html)!

### the ssc_model class
It contains three files:  

1. `model.py`  
in which all the detils of the SSC object you are creating are contained:  
energy (Lorentz Factor) grid, time grid, dimension and magnetic field of the emitting region, injected spectrum.

2. `numerics.py`  
in which all the calculations are performed. A parent population of electrons is evolved, cooling on a magnetic and radiation field.  
naima is used for evaluating the emissivities (synchrotron and SSC photon spectra)

3. `constants.py`    
In which we define some variables (speed of light, Thomson cross section) with module-wise scope.


### ssc_test.py
All you need is here, you have to define four python dictionaries with the details of the model you want to create:  
1. a time grid  
2. a energy grid  
3. details of the emitting region (magnetic field, radius, escape time)  
4. details of the injected spectrum
5. distance [in astropy units] (used for evaluating the SEDs, TODO: include it in one of the dictionaries - emitting region)


```python
time_grid = dict(time_min = 0., time_max = 3., time_bins = 50)
gamma_grid = dict(gamma_min = 2., gamma_max = 2e5, gamma_bins = 20)
emission_region = dict(R = 1e16, B = 1., t_esc = 1.5)
injected_spectrum = dict(norm = 1e-3, alpha = -2., t_inj = 1.5)
dist = 2*u.Mpc
```

after this you define the SSC object with this dictionaries

```python
SSC = model(time_grid, gamma_grid, emission_region, injected_spectrum)
```

the numerics associated to it
```python
num = numerics(SSC)
```

Let it evolve!
```python
N_e = num.evolve()
```
and you get an array `N_e` with the evolved spectrum (over the energy grid defined).

Then we invoke naima's Synchrotron and IC objects for building the photon SED.
```python
SYN = num.synchrotron(N_e)
IC = num.inverse_compton(N_e)
```

the rest of the macro is for nicely plotting the results, after invoking naima Sycnhrotron and IC objects

