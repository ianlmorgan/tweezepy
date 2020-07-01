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
* 'smmcalibration' - tools to determine force of probe trajectories by calculating and performing maximum likelihood estimation fits to power spectral density and allan variance

# Example use:

```python
>>> from tweezepy.smmcalibration import PSD
>>> psd = PSD(xtrace,fsample)
>>> psd.plot()
```
<img src="examples/example1.png" width="600">

```python
>>> from tweezepy.smmcalibration import AV
>>> psd = PSD(xtrace,fsample)
>>> psd.plot()
```
<img src="examples/example2.png" width="600">

