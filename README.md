# tweezepy
This is [tweezepy](https://github.com/ianlmorgan/tweezepy), a package of single-molecule pulling experiment related analysis code.

# How to install
Using pip:

    pip install tweezepy
Using setuptools:

Clone repository onto local machine. 
Navigate to directory: 
    
    cd path/to/tweezepy

    python setup.py install
    
# Contents
The `tweezepy` package includes the following modules:
* 'smmcalibration' - tool to determine force of probe trajectories by calculating and performing maximum likelihood estimation fits to power spectral density and allan variance
* 'simulations' - tool to simulate single-molecule probe trajectories

# Example use:
Simulate data:
```python
>>> import matplotlib.pyplot as plt
>>> from tweezepy.simulations import downsampled_trace
>>> alpha,kappa,fsample,N = 1e-5,.001,400,10240
>>> xtrace = downsampled_trace(alpha,kappa,fsample,N)
>>> plt.plot(xtrace)
>>> plt.show()
```
<img src="examples/example_trace.png" width=600>

Power spectral density:
```python
>>> from tweezepy.smmcalibration import PSD
>>> psd = PSD(xtrace,fsample)
>>> pars,errs,covs = psd.mlefit()
>>> psd.plot()
```

<img src="examples/example1.png" width="600">

Allan variance:
```python
>>> from tweezepy.smmcalibration import AV
>>> av = AV(xtrace,fsample)
>>> pars,errs,covs = av.mlefit()
>>> av.plot()
```
<img src="examples/example2.png" width="600">

