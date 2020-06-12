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
import numpy as np

import pandas as pd


from scipy.optimize import minimize,curve_fit
from scipy.signal import welch
from scipy.stats import gamma,chi2
from scipy.special import erf
from iminuit import Minuit
from statsmodels.tools.numdiff import approx_hess
from numpy import pi,exp,sin,sinh,cos,cosh



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
        etas = self.results['etas']
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
    def __init__(self, trace, freq, nperseg = 1024):
        self.trace = trace
        self.freq = freq
         
        f, Pxx_den = welch(trace, freq, nperseg=nperseg)
        # By default, scipy returns the one-sided PSD
        # Lansdorp et al. (2012) uses the positive part of the two-sided PSD
        # Not going to lie, figuring this out drove me crazy
        # Divide by 2 to get the positive part of the two-sided PSD
        Pxx_den /= 2
        b = (2*len(trace)/nperseg) - 1
        etas = np.full_like(f,b)
        self.results = pd.DataFrame({'f':f,'psd': Pxx_den,'etas':etas})
    def plot(self):
        """
        Plots the overlapping allan variance on a loglog plot.

        Returns
        -------
        ax : ax.
            The matplotlib axis instance.

        """
        fig,ax = plt.subplots()
        ax.plot('f','psd',data=self.results)
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
        scale = self.results.yhat/self.results.etas
        self.results['psd_std'] = gamma.std(self.results.psd,
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
            ax.errorbar('f','psd',yerr='psd_std',
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
    f, Pxx_den = welch(xtrace, freq, nperseg=nperseg)
    Pxx_den /= 2
    b = (2*len(xtrace)/nperseg) - 1
    etas = np.full_like(f,b)
    return f,etas,Pxx_den
def SMMAV(t,a,k):
    """
    Modified Eq. 17 from Lansdorp et al. (2012) 
    for the single-molecule allan variance.

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
                               2.*a/(k*t) * exp(-k*t/a) -
                               a/(2.*k*t) * exp(-2.*k*t/a) -
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
    PSD = 2.*kT*a/k**3. * (k + (2.*a*fs*sin(pi*f/fs)**2. * sinh(k/(a*fs))) / (cos(2.*pi*f/fs) - 
                                                                              cosh(k/(a*fs)))
                          )
    return PSD
def negLL(p,func,x,e,y):
    """Takes params, taus,etas, and oavs and returns the negative log likelihood. """
    yhat = func(x,*p)
    loglikelihood = gamma.logpdf(y,e,scale=yhat/e).sum()
    negloglikelihood = -1 * loglikelihood
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
    popt,pcov = curve_fit(func,x,y,p0 = guess)
    # Minimize the negative log likelihood
    
    negLL = lambda a,k: -gamma.logpdf(y,e,scale=func(x,a,k)/e).sum()
    m = Minuit(negLL,
               pedantic = False,
               a=popt[0],k=popt[1],
               error_a=0.1*popt[0], error_k=0.1*popt[1],
               errordef = 0.5)
    m.migrad()
    params = m.np_values()
    se = m.np_errors()
    print(m.get_param_states())
    return params,se


# Confidence Intervals
ONE_SIGMA_CI = erf(1/np.sqrt(2))
#    = 0.68268949213708585
def confidence_interval(oav, eta, ci=ONE_SIGMA_CI):
    ci_l = min(np.abs(ci), np.abs((ci-1))) / 2
    ci_h = 1 - ci_l

    # function from scipy, works OK, but scipy is large and slow to build
    chi2_l = chi2.ppf(ci_l, eta)
    chi2_h = chi2.ppf(ci_h, eta)

    var_l = eta * oav / chi2_h  # NIST SP1065 eqn (45)
    var_h = eta * oav / chi2_l
    return var_l,var_h