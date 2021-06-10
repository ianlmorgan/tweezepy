import autograd.numpy as np
"""
Analytical models
"""

def SMMAV(t,g,k,kT = 4.1):
    """
    Analytical function for the AV of a trapped bead.
    Eq. 17 from Lansdorp et al. (2012) for the single-molecule allan variance.

    Parameters
    ----------
    t : array
        taus.
    g : float
        drag coefficient.
    k : float
        spring constant.
    kT : float
        Thermal energy in pN nm. Default value is 4.1

    Returns
    -------
    oav : array
        theoretical allan variance.

    """
    tc = np.true_divide(g,k)
    oav = 2.*kT*tc/(k*t) * (1. + 2. * (tc/t)*np.exp(-t/tc) - (tc/(2.*t))*np.exp(-2.*t/tc) - 3.*tc/(2.*t))
    return oav

def lansdorpPSD(f,fs,g,k,kT = 4.1):
    """
    Analytical function for the PSD of a trapped bead with aliasing and lowpass filtering.
    Eq. 7 in Lansdorp et al. (2012).

    Parameters
    ----------
    f : array-like
        frequency.
    fs : float
        Acquisition frequency.
    g : float
        drag coefficient.
    k : float
        spring constant.
    kT : float
        Thermal energy in pN nm, defaults to 4.1

    Returns
    -------
    PSD : array
        theoretical power spectral density.
    """
    tc = g/k
    fc = 1./tc
    PSD = 2.*kT*tc/k * (1. + 2.*tc*fs*np.sin(np.pi*f/fs)**2 * np.sinh(fc/fs)/(np.cos(2.*np.pi*f/fs) - np.cosh(fc/fs)))
    return PSD
    

def aliasPSD(f,fs,a,k, kT = 4.1):
    """
    Analytical function for the PSD of a trapped bead with aliasing.
    Eq. 8 in Lansdorp et al. (2012).

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
    kT : float
        Thermal energy in pN nm. Default value is 4.1

    Returns
    -------
    PSD : array
        theoretical power spectral density.
    """
    kT = 4.1 # thermal energy in pNnm
    return kT/(k*fs) * (np.sinh(k/(a*fs))/(np.cosh(k/(a*fs))-np.cos(2*np.pi*f/fs)))
