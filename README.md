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