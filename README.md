# Tweezepy
[![DOI](https://zenodo.org/badge/261266475.svg)](https://zenodo.org/badge/latestdoi/261266475)

This is [Tweezepy](https://github.com/ianlmorgan/tweezepy), a Python package for calibrating forces in single-molecule force spectroscopy video-tracking experiments using the power spectral density (PSD) and Allan variance (AV).

## Documentation
Read the documentation for `Tweezepy` [here](https://tweezepy.readthedocs.io/).

## How to install
The simplest method of installing the `Tweezepy` package is via the [Python Package Index](https://packaging.python.org/glossary/#term-python-package-index-pypi) (PyPI). To install from PyPI, you will need to be able to run python from the command line and make sure you have [pip](https://packaging.python.org/key_projects/#pip) available.

Install from PyPI:

    pip install tweezepy
An alternative method to install `Tweezepy` is with setuptools.  Clone the repository onto a local machine, then navigate to the directory.

Using setuptools:
    
    cd path/to/tweezepy

    python setup.py install
    
## Contents
The `Tweezepy` package includes the following modules:
* 'smmcalibration' - classes for calibration methods using the PSD and AV
* 'expressions' - functions with closed-form expressions for thermal motion in the PSD and AV
* 'MLE' - classes for maximum likelihood estimation (MLE) and Monte Carlo Markov chain (MCMC) sampling
* 'allanvar' - tools for calculating the AV and equivalent degrees of freedom
* 'simulations' - tools to simulate bead thermal motion

## Example use:
Simulate data:
```python
>>> import matplotlib.pyplot as plt
>>> from tweezepy.simulations import downsampled_trace
>>> xtrace = downsampled_trace()
>>> plt.plot(xtrace)
>>> plt.show()
```
<img src="docs/example_trace.png" width=600>

Power spectral density:
```python
>>> from tweezepy.smmcalibration import PSD
>>> psd = PSD(xtrace,fsample)
>>> psd.mlefit()
>>> print(psd.results)
>>> psd.plot()
```

<img src="docs/example_PSD.png" width="600">

Allan variance:
```python
>>> from tweezepy.smmcalibration import AV
>>> av = AV(xtrace,fsample)
>>> av.mlefit()
>>> print(av.results)
>>> av.plot()
```
<img src="docs/example_AV.png" width="600">

# Jupyter notebooks with examples
Jupyter notebooks are interactive Python scripts, embedded in a browser, allowing you to manipulate data and display plots like easily. For guidance on installing jupyter, please refer to https://jupyter.org/install.

See /docs for some examples in notebook format.

github formats the notebooks into nice web-pages, for example
* [Basic usage](https://github.com/ianlmorgan/tweezepy/tree/master/docs/basic_usage.ipynb)
* [MCMC](https://github.com/ianlmorgan/tweezepy/tree/master/docs/MCMC.ipynb)