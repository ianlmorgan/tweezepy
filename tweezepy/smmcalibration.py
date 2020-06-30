# -*- coding: utf-8 -*-
"""
The allanvar module.

Example usage:
    from tweezepy.allanvar import AV
    trace = np.arange(0,5,0.1)
    av = AV(trace,0.1)
    av.plot()
    
Created on Mon May  4 11:51:23 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""

import matplotlib.pyplot as plt
import autograd.numpy as np 

import pandas as pd

from autograd import hessian
from autograd.scipy.stats import gamma
from scipy.optimize import curve_fit,minimize
from scipy.signal import welch
from scipy import stats

class calibration:
    def __init__(self, trace, fsample, calibfunc,**kwargs):
        self.trace = trace
        self.fsample = fsample
        x,n,y = calibfunc(trace,fsample,**kwargs)
        yerr = stats.gamma.std(n,scale = y/n)
        self.results = pd.DataFrame({'x':x,'n':n,'y':y,'yerr':yerr})
        
    def plot(self, ax = None,**kwargs):
        if ax == None:
            fig,ax = plt.subplots()
        ax.errorbar('x','y',yerr='yerr',data=self.results,
                    fmt = 'o', **kwargs,zorder = 0)
        if 'yfit' in self.results:
            ax.plot('x','yfit',data=self.results,lw = 2, label='', zorder = 1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax
    def mlefit(self,func,pedantic = True,**kwargs):
        if pedantic == False:
            np.seterr('ignore')
        elif pedantic == True:
            np.seterr('warn')
        results = self.results
        params, se, cov = MLEfit(func,
                                 results.x,
                                 results.n,
                                 results.y,
                                 **kwargs)
        self.params = params
        self.se = se
        self.cov = cov
        self.results['yfit'] = func(results.x,*params)
        return params,se,cov
    
class AV(calibration):
    def __init__(self, trace, fsample,**kwargs):
        calibration.__init__(self, trace, fsample, allanvar, **kwargs)
    def plot(self,**kwargs):
        ax = calibration.plot(self,**kwargs)
        ax.set_xlabel(r'$\tau$ (s)')
        ax.set_ylabel(r'$\sigma^2$ [nm$^2$]')
        return ax
    def mlefit(self, fitfunc = None, **kwargs):
        if fitfunc == None:
            fitfunc = SMMAV2
        return calibration.mlefit(self, fitfunc, **kwargs)

class PSD(calibration):
    def __init__(self, trace, fsample,**kwargs):
        calibration.__init__(self, trace, fsample, psd, **kwargs)
    def plot(self,**kwargs):
        ax = calibration.plot(self,**kwargs)
        ax.set_xlabel(r'f [Hz]')
        ax.set_ylabel(r'PSD [nm$^2$/Hz]')
        return ax
    def mlefit(self, fitfunc = 'lansdorpPSD', **kwargs):
        if fitfunc == 'lansdorpPSD':
            fitfunc = lansdorpPSD
        elif fitfunc == 'aliasPSD':
            fitfunc = aliasPSD        
        mlefunc = lambda x,a,k: fitfunc(x,self.fsample,a,k)
        return calibration.mlefit(self, mlefunc, **kwargs)
        
def allanvar(xtrace,freq,taus = 'octave'):
    """
    Estimate overlapping allan variance 
    Takes an array of numbers and returns the overlapping allan variance.
    Returns the taus, etas, and oavs.

    Parameters
    ----------
    xtrace : array, series, or list
        1-D array of numbers.
    freq : float
        Frequency of acquisition.

    Returns
    -------
    taus : array
        Taus.
    etas : array
        Etas.
    oavs : array
        Overlapping allan variance.

    """
    xtrace = np.asarray(xtrace) # convert to numpy array
    dt = 1.0/freq # shutter speed in s tau_s
    N = len(xtrace)
    if taus == 'all':
        # all-tau sampling not particularly useful but why not?
        m = np.linspace(1.0,N,N,dtype='int')
    elif taus == 'octave':
        # octave sampling break bin sizes, m, into powers of 2^n
        maxn = int(np.floor(np.log2(N/2.))) # m =< N/2
        m = np.logspace(0,maxn,maxn+1,base=2,dtype='int')  #bin sizes
    taus = m*dt # tau = m*tau_c
    etas = np.floor(N/m-1)/2 # shape factors
    # See eratum for Eq. 18b
    phasedata = np.cumsum(xtrace) * dt # convert frequency data to phase data
    phasedata = np.insert(phasedata, 0, 0) # phase data should start at 0
    oavs = np.empty_like(taus) # creates array with the same shape as taus
    # Calculate OAV as in Eq. 18a (in erratum)
    for i,mj in enumerate(m):
        # calculate overlapping allan variance for each m
        d2,d1,d0 = phasedata[2*mj:],phasedata[mj:],phasedata[:]
        n = min(len(d2),len(d1),len(d0))
        v_arr = d2[:n]-2.*d1[:n]+d0[:n]
        s = np.sum(v_arr*v_arr)
        oavs[i] = s/(2.*(N-2.*mj+1.)*(mj*dt)**2.)
    return taus,etas,oavs
def psd(xtrace,freq,nperseg = None,return_onesided=False,**kwargs):
    """Takes 1-D array and returns frequency,etas, and psd."""
    N = len(xtrace)
    if nperseg is None:
        nperseg = N//27
    f, dens = welch(xtrace, freq, nperseg=nperseg,return_onesided=return_onesided,**kwargs)
    msk = f>0
    f, dens = f[msk], dens[msk]
    b = 2*N//nperseg - 1
    etas = np.full_like(f,b)
    return f,etas,dens
def SMMAV(t,a,k):
    """
    Modified Eq. 17 from Lansdorp et al. (2012) 
    for the single-molecule allan variance.
    ## FIX-ME
    ## Still unclear whether I want to use this or the other equation
    ## The issue is that minimization can be hard with parameters that are 
    ## different by several orders of magnitude.

    Parameters
    ----------
    t : array
        taus.
    a : float
        alpha/kappa.
    k : float
        kappa.

    Returns
    -------
    oav : array
        predicted allan variance.

    """
    kT = 4.1 # in pN*nm
    oav = 2.*kT*a/(k*t) * (1 + 2*a*np.exp(-t/a)/t - 
                           a*np.exp(-2*t/a)/(2*t) - 3*a/(2*t))
    return oav
def SMMAV2(t,a,k):
    """
    Eq. 17 from Lansdorp et al. (2012) for the single-molecule allan variance.

    Parameters
    ----------
    t : array
        taus.
    alpha : float
        alpha.
    kappa : float
        kappa.

    Returns
    -------
    oav : array
        predicted allan variance.

    """
    kT = 4.1 # in pN*nm
    oav = 2.*kT*a/(k**2.*t) * (1. +
                               2.*a/(k*t) * np.exp(-k*t/a) -
                               a/(2.*k*t) * np.exp(-2.*k*t/a) -
                               3*a/(2.*k*t)
                              )
    return oav
def lansdorpPSD(f,fs,a,k):
    """
    Eq. 7 in Lansdorp et al. (2012) for the single-molecule power spectral density that accounts for btoh aliasing and boxcar filtering. Takes frequency and gives the PSD.

    Parameters
    ----------
    f : array-like
        frequency.
    fs : float
        Acquisition frequency.
    a : float
        alpha.
    k : float
        kappa.

    Returns
    -------
    PSD : array
        power spectral density.

    """
    kT = 4.1 # in pN*nm
    PSD = 2.*kT*a/k**3. * (k + 
                           (2.*a*fs*np.sin(np.pi*f/fs)**2. * np.sinh(k/(a*fs))) / (np.cos(2.*np.pi*f/fs) - 
                                                                                   np.cosh(k/(a*fs)))
                          )
    return PSD
def aliasPSD(f,fs,a,k):
    kT = 4.1 # thermal energy in pNnm
    return kT/(k*fs) * (np.sinh(k/(a*fs))/(np.cosh(k/(a*fs))-np.cos(2*np.pi*f/fs)))
   

def negLL(p,func,x,e,y):
    """Takes params, taus,etas, and oavs and returns the negative log likelihood. """
    yhat = func(x,*p)
    scale = yhat/e
    likelihood = gamma.pdf(y/scale,e)/scale
    negloglikelihood = -np.sum(np.log(likelihood))
    return negloglikelihood 
def MLEfit(func,x,e,y,guess = None,**kwargs):
    """
    Performs a basic maximum likelihood estimation on xtrace.
    Returns alpha and kappa.

    Parameters
    ----------
    taus : array_like
        overlapping time bins.
    etas : array_like
        shape factors.
    oavs : array_like
        overlapping allan variances
    guess : list, optional
        alpha and kappa guesses. The default is [1.15E-5,0.001].

    Returns
    -------
    params : array
        best-fit parameters 
    se : array
        error associated with parameters

    Notes
    -----
    Need to update so it does fixed alpha.
    """  
    x = np.asarray(x)
    e = np.asarray(e)
    y = np.asarray(y)
    
    # Bounds on biased nonlinear least squares fit
    # ((alpha lower,kappa lower),(alpha upper,kappa upper))
    bounds = ((0,0),(np.inf,np.inf)) 
    popt,_ = curve_fit(func,x,y,p0 = guess, bounds = bounds)
    # Minimize the negative log likelihood
    _negLL = lambda p: negLL(p,func,x,e,y)
    hess = hessian(_negLL)
    fit = minimize(_negLL,x0=popt,method='Nelder-Mead',**kwargs)
    pars = fit['x']
    var = np.linalg.inv(hess(fit['x']))
    errors = np.sqrt(np.diag(var))
    return pars,errors,var

