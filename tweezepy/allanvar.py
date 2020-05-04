# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:51:23 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""

import numpy as np

from scipy.optimize import minimize,curve_fit
from scipy.stats import gamma

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


def MLEfit(xtrace,freq,guess = [1.15E-5,0.001]):
    """
    Performs a basic maximum likelihood estimation on xtrace.
    Returns alpha and kappa.

    Parameters
    ----------
    xtrace : series,array, or list?
        The xtrace data.
    freq : float or int
        data acquisition frequency.
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
    taus,etas,oavs = allanvar(xtrace,freq)    
    popt,pcov = curve_fit(SMMAV,taus[:-4],oavs[:-4],p0 = guess)
    # Define the negative log likelihood
    negLL = lambda p,t,e,o: -np.sum(gamma.logpdf(oavs,etas,scale=SMMAV(t,p[0],p[1])/etas)) 
    results = minimize(negLL, # minimize the negative log likelihood 
                       popt, # initial guess
                       args = (taus[:-4],etas[:-4],oavs[:-4]), 
                       method = 'Nelder-Mead')
    alpha,kappa = results['x']
    return alpha,kappa

