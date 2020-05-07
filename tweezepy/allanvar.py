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
from scipy.stats import gamma

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
        taus,etas,oavs = allanvar(trace, freq = 400)
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
        ax.plot('taus','oavs',data=self.results,marker='o',lw=0)
        ax.set_xlabel(r'$\tau$ (s)')
        ax.set_ylabel(r'$\sigma^2$ (s)')
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax
    def mlefit(self, guess = [1.15E-5,0.001]):   
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
        alpha, kappa = MLEfit(self.results.taus,
                                        self.results.etas,
                                        self.results.oavs,
                                        guess = guess)
        self.alpha = alpha
        self.kappa = kappa
        self.results['yhat'] = SMMAV(self.results.taus,
                                     self.alpha,
                                     self.kappa)
        scale = self.results.yhat/self.results.etas
        self.results['ostd'] = gamma.std(self.results.oavs,
                                         self.results.etas,
                                         scale=scale)
        
        def plot():
            """
            Reassigns the plotting function after the mle fit to plot errors and fits.  
            
            Returns
            -------
            ax : ax.
                The matplotlib axis instance.

            """
            fig, ax = plt.subplots()
            # Plot the points with errorbars
            ax.errorbar('taus','oavs',yerr='ostd',
                         data=self.results,fmt = 'o',label = None)
            # Plot the function with best-fit values from the mle fit
            ax.plot('taus','yhat',
                     data=self.results,label = 'MLE fit')
            ax.set_xlabel(r'$\tau$ (s)')
            ax.set_ylabel(r'$\sigma^2$ (s)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend()
            return ax
        self.plot = plot
        return alpha,kappa
    
def allanvar(xtrace,freq):
    """
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
    # octave sampling break bin sizes, m, into powers of 2^n
    maxn = int(np.floor(np.log2(N/2.))) # m =< N/2
    m = np.logspace(0,maxn,maxn+1,base=2,dtype='int')  #bin sizes
    taus = m*dt # tau = m*tau_c
    etas = (N/m-1)/2 # shape factors
    phasedata = np.cumsum(xtrace) * dt # convert frequency data to phase data
    phasedata = np.insert(phasedata, 0, 0) # phase data should start at 0
    oavs = np.zeros_like(taus) # creates array of zeros with the same length as taus
    for i,mj in enumerate(m):
        # calculate overlapping allan variance for each m
        d2,d1,d0 = phasedata[2*mj:],phasedata[mj:],phasedata[:]
        n = min(len(d2),len(d1),len(d0))
        v_arr = d2[:n]-2.*d1[:n]+d0[:n]
        s = np.sum(v_arr*v_arr)
        oavs[i] = s/(2.*(N-2.*mj+1.)*(mj*dt)**2.)
    return taus,etas,oavs

def SMMAV(t,alpha,kappa):
    """
    Function from Lansdorp et al. (2012) for the single-molecule allan variance.

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
        allan variance.

    """
    kT = 4.1 # in pN*nm
    oav = 2*kT*alpha/(kappa**2.*t) * (1. +
                                       2.*alpha/(kappa*t) * np.exp(-kappa*t/alpha) -
                                       alpha/(2.*kappa*t) * np.exp(-2.*kappa*t/alpha) -
                                       3*alpha/(2.*kappa*t))
    return oav
def negLL(p,t,e,o):
    """Takes params, taus,etas, and oavs and returns the negative log likelihood. """
    a,k = p[0],p[1]
    yhat = SMMAV(t,a,k)
    return -np.sum(gamma.logpdf(o,e,scale=yhat/e)) 
def MLEfit(taus,etas,oavs,guess = [1.15E-5,0.001]):
    """
    Performs a basic maximum likelihood estimation on xtrace.
    Returns alpha and kappa.

    Parameters
    ----------
    taus : array
        overlapping time bins.
    etas : array
        shape factors.
    oavs : array
        overlapping allan variances
    guess : list, optional
        alpha and kappa guesses. The default is [1.15E-5,0.001].

    Returns
    -------
    alpha : float
        alpha - drag on bead [input units].
    kappa : float
        kappa - spring constant in pN/nm.

    Notes
    -----
    Need to update so it does fixed alpha.
    """  
    popt,pcov = curve_fit(SMMAV,taus[:-4],oavs[:-4],p0 = guess)
    results = minimize(negLL, popt, args = (taus[:-4],etas[:-4],oavs[:-4]), 
                       method = 'Nelder-Mead')
    alpha,kappa = results['x']
    return alpha,kappa

