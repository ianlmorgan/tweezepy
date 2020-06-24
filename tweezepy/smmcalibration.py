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


class AV:
    """
    The AV object holds allan variance related functions. 
    
    Parameters
    ----------
    trace : array
        1-D array of data collected with acquisition freq
    freq : float or int
        acquisition frequency in Hz
        
    Attributes
    ----------
    trace : array
        stored trace
    freq : float or int
        stored acquisition frequency
    results : dataframe
        allan variance data
    """
    def __init__(self, trace, freq):
        self.trace = trace
        self.freq = freq
        taus,etas,oavs = allanvar(trace, freq = freq)
        self.results = pd.DataFrame({'taus':taus,'oavs': oavs,'etas':etas})
    def plot(self):
        """
        Plots the overlapping allan variance on a loglog plot.

        Returns
        -------
        ax : ax.
            The matplotlib axis instance.

        """
        fig,ax = plt.subplots()
        taus = self.results['taus']
        oavs = self.results['oavs']
        #Error bars not implemented before fit
        #var_l,var_h = confidence_interval(oavs,etas)
        #ax.errorbar(taus,oavs,yerr=[var_l,var_h],fmt='o')
        ax.plot(taus,oavs,'o')
        ax.set_xlabel(r'$\tau$ (s)')
        ax.set_ylabel(r'$\sigma^2$ [nm$^2$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax
    def mlefit(self,fixed_alpha = None,guess = [1.15E-5,0.001],**kwargs):   
        """
        Parameters
        ----------
        guess : list, optional
            parameter guesses. The default is [1.15E-5,0.001].

        Returns
        -------
        alpha : float
            alpha
        kappa : float
            kappa in pNnm

        """          

        #assert type(guess) is list, "guess must be a list: %r" %guess
            
        if fixed_alpha == None:                
            func = lambda t,a,k: SMMAV(t,a,k)
        else:
            assert type(fixed_alpha) is float, "fixed_alpha must be a float: %r" %fixed_alpha
            guess = [.001]
            func = lambda t,k: SMMAV(t,fixed_alpha,k)
        params, se = MLEfit(func,
                            self.results.taus,
                            self.results.etas,
                            self.results.oavs,
                            guess = guess,
                            **kwargs)
        self.params = params
        self.se = se
        self.results['yhat'] = func(self.results.taus, *self.params)
        scale = self.results.yhat/self.results.etas
        self.results['ostd'] = gamma.std(self.results.oavs,
                                         self.results.etas,
                                         scale=scale)
            
        def plot():
            """
            Reassigns the plotting function after the mle fit to plot errors (standard deviations) and fits.  
            
            Returns
            -------
            ax : ax.
                The matplotlib axis instance.

            """
            fig, ax = plt.subplots()
            # Plot the points with errorbars
            ax.errorbar('taus','oavs',yerr='ostd',
                         data=self.results,fmt = 'o',label = None, zorder=0)
            # Plot the function with best-fit values from the mle fit
            ax.plot('taus','yhat',
                     data=self.results,label = 'MLE fit', zorder=1)
            ax.set_xlabel(r'$\tau$ [s]')
            ax.set_ylabel(r'$\sigma^2$ [nm$^2$]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend()
            return ax
        self.plot = plot
        return self.params,self.se
class PSD:
    def __init__(self, trace, freq,**kwargs):
        self.trace = trace
        self.freq = freq
         
        #f, Pxx_den = welch(trace, freq, nperseg=1024)
        # By default, scipy returns the one-sided PSD
        # Lansdorp et al. (2012) uses the positive part of the two-sided PSD
        # Not going to lie, figuring this out drove me crazy
        # Divide by 2 to get the positive part of the two-sided PSD
        f, etas, dens = psd(trace,freq,**kwargs)
        #Pxx_den /= 2
        stds = stats.gamma.std(etas,scale = dens/etas)
        self.results = pd.DataFrame({'f':f,'etas':etas,'psd': dens,'psd_stds':stds})
    def plot(self):
        """
        Plots the overlapping allan variance on a loglog plot.

        Returns
        -------
        ax : ax.
            The matplotlib axis instance.

        """
        fig,ax = plt.subplots()
        ax.errorbar('f','psd',yerr='psd_stds',data=self.results)
        ax.set_xlabel(r'f [Hz]')
        ax.set_ylabel(r'PSD [nm$^2$/Hz]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax
    
    def mlefit(self,fixed_alpha = None,guess = [1.15E-5,0.001],**kwargs):   
        """
        Parameters
        ----------
        guess : list, optional
            parameter guesses. The default is [1.15E-5,0.001].

        Returns
        -------
        alpha : float
            alpha
        kappa : float
            kappa in pNnm

        """
        assert type(guess) is list, "guess must be a list: %r" %guess
        if fixed_alpha == None:
            func = lambda t,a,k: SMMPSD(t,self.freq,a,k)
        else:
            assert type(fixed_alpha) is float, "fixed_alpha must be a float: %r" %fixed_alpha
            func = lambda t,k: SMMPSD(t,self.freq,fixed_alpha,k)
        params,se = MLEfit(func,
                           self.results.f,
                           self.results.etas,
                           self.results.psd,
                           guess = guess,
                           **kwargs)
        self.params = params
        self.se = se
        self.results['yhat'] = func(self.results.f, *params)
        #scale = self.results.yhat/self.results.etas
        #self.results['psd_std'] = gamma.std(self.results.psd,
        #                                    self.results.etas,
        #                                    scale=scale)
            
        def plot():
            """
            Reassigns the plotting function after the mle fit to plot errors (standard deviations) and fits.  
            
            Returns
            -------
            ax : ax.
                The matplotlib axis instance.

            """
            fig, ax = plt.subplots()
            # Plot the points with errorbars
            ax.errorbar('f','psd',yerr='psd_stds',
                         data=self.results,fmt = 'o',label = None,zorder=0)
            # Plot the function with best-fit values from the mle fit
            ax.plot('f','yhat',
                     data=self.results,label = 'MLE fit',zorder=1)
            ax.set_xlabel(r'f [Hz]')
            ax.set_ylabel(r'PSD [nm$^2$/Hz]')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend()
            return ax
        self.plot = plot
        return params,se
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
    etas = (N/m-1)/2 # shape factors
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
def psd(xtrace,freq, nperseg = 1024):
    """Takes 1-D array and returns frequency,etas, and psd."""
    f, dens = welch(xtrace, freq, nperseg=nperseg,return_onesided=False)
    msk = f>0
    f, dens = f[msk], dens[msk]
    b = 2*np.floor(len(xtrace)/nperseg) - 1
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
def SMMPSD(f,fs,a,k):
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
def negLL(p,func,x,e,y):
    """Takes params, taus,etas, and oavs and returns the negative log likelihood. """
    yhat = func(x,*p)
    scale = yhat/e
    likelihood = gamma.pdf(y/scale,e)/scale
    negloglikelihood = -np.sum(np.log(likelihood))
    return negloglikelihood 
def MLEfit(func,x,e,y,guess = [1.15E-5,0.001],**kwargs):
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

    popt,_ = curve_fit(func,x,y,p0 = guess)
    # Minimize the negative log likelihood
    _negLL = lambda p: negLL(p,func,x,e,y)
    hess = hessian(_negLL)
    fit = minimize(_negLL,x0=popt,method='Nelder-Mead')
    pars = fit['x']
    var = np.linalg.inv(hess(fit['x']))
    errors = np.sqrt(np.diag(var))
    return pars,errors

