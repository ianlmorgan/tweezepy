# -*- coding: utf-8 -*-
"""
The smmcalibration module.

Example usage:
    from tweezepy.simulations import simulate_trace
    from tweezepy.smmcalibration import AV
    trace = simulate_trace()
    av = AV(trace,100)
    av.plot()
    
Created on Mon May  4 11:51:23 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
import corner
import emcee
import matplotlib.pyplot as plt
import autograd.numpy as np 
import pandas as pd

from autograd import hessian
from autograd.scipy.stats import gamma
from scipy.optimize import curve_fit,minimize
from scipy.signal import welch
from scipy import stats
class MLEfit:
    """
    Perform maximum likelihood estimation and uncertainty calculation.
    """
    def __init__(self,func,x,e,y,guess = [1.15E-5,0.002],*args,**kwargs):
        """
        Minimize negative log likelihood and find optimal parameters
        Args:
            func (func): Function to fit to.
            x (list or array): x data
            e (list or array): shape parameter for gamma function
            y (list or array): y data
            guess (list, optional): initial parameter guesses. Defaults to [1.15E-5,0.002].
        """
        
        self.func = func
        x = np.asarray(x); self.x=x
        e = np.asarray(e); self.e=e
        y = np.asarray(y); self.y=y
        # Bounds on biased nonlinear least squares fit
        # ((alpha lower,kappa lower),(alpha upper,kappa upper))
        bounds = ((0,0),(np.inf,np.inf)) 
        popt,_ = curve_fit(func,x,y,p0 = guess, bounds = bounds)
        # Minimize the negative log likelihood
        self.LL = lambda p: loglikelihood(x,e,y,func(x,*p))
        self.negLL =  lambda p: -loglikelihood(x,e,y,func(x,*p)) #(p,func,x,e,y)
        self.fit = minimize(self.negLL,x0=popt,method='Nelder-Mead',*args,**kwargs)
    def mcmc(self, walkers = 32, steps = 1500):
        """Run a Markov Chain Monte Carlo Sampler to determine possible asymmetric fitting uncertainty.

        Args:
            scale (list, optional): Parameter scale. Defaults to [1e-5,1e-3].
            N (int, optional): Number of samples. Defaults to 1500.
        """
        params = self.params
        scale = np.power(10,np.floor(np.log10(params)))
        pos = params + 1e-4* np.random.randn(walkers,len(params))*scale
        nwalkers,ndims = pos.shape
        self.sampler = emcee.EnsembleSampler(nwalkers,ndims,self.LL)
        self.sampler.run_mcmc(pos, steps,progress = True)
        

    def corner_plot(self,labels = [r'$\alpha$',r'$\kappa$'], quantiles = [0.16, 0.5, 0.84],
                    truths=None, levels = None,**kwargs):
        """
        Make a corner plot with parameter distributions.

        Args:
            labels (list, optional): Labels for axis. Defaults to [r'$\alpha$',r'$\kappa$'].
            quantiles (list, optional): Parameter distribution quantiles. Defaults to 1 standard deviation [0.16, 0.5, 0.84].
            truths ([list], optional): True or fitted values. Defaults to None.
            levels ([type], optional): . Defaults to None.
        """
        if truths is None:
            truths = self.params
        if levels is None:
            levels = (1-np.exp(-0.5),)
        corner.corner(self.samples,
                      labels = labels,
                      truths = truths,
                      quantiles = quantiles,
                      levels = levels,
                      **kwargs)
    @property
    def params(self):
        return self.fit['x']
    @property
    def cov(self):
        hess = hessian(self.negLL)
        return np.linalg.inv(hess(self.params))
    @property
    def errors(self):
        cov = self.cov
        return np.sqrt(2*np.diag(cov))
    @property
    def redchi2(self):
        """
        Calculates the reduced chi-squared from the likelihood ratio. 
        Returns:
            rechi2 (float): Reduced chi-squared
        """
        func = self.func
        x = self.x
        e = self.e
        y = self.y
        params = self.params
        chi2 = -2*(loglikelihood(x,e,y,func(x,*params))-loglikelihood(x,e,y,y))
        redchi2 = chi2/(len(y)-len(params)-1) # needs to be checked
        return redchi2
    @property
    def mc_errors(self,discard = 100, thin = 10, percentiles = [16,50,84]):
        """
        Returns 16th and 84th percentiles from monte carlo samples.

        Args:
            discard (int): how many steps to discard
            thin (int): how many steps to thin by
            percentiles (list): Percentiles to calculate. Defaults to 16th and 84 pecentiles [16,50,84].
        Returns:
            q (array): Pecentiles for each parameter.
        """
        samples = self.sampler.get_chain(discard=discard, thin=thin, flat=True)
        q = []
        for i,p in enumerate(self.params):
            mcmc = np.percentile(samples[:, i], percentiles)
            q.append(np.diff(mcmc))
        q = np.array(q)*np.sqrt(2) # Still a bit confused why I have to multiply the errors by sqrt(2)
        return q

class calibration(MLEfit):
    """
    Base class for PSD and AV calibration.

    Args:
        MLEfit (class): A base class for maximum likelihood fitting.
    """
    def __init__(self, trace, fsample):
        """
        Loads in trace and fsample as class attributes.

        Args:
            trace (array): Bead trajectory positions.
            fsample (array): Camera sample rate
        """
        self.trace = trace
        self.fsample = fsample

    def plot(self, ax = None,**kwargs):
        """
        Generate plot of data and fit to data (if applicable).

        Args:
            ax (axes, optional): If applicable, define axes to plot on. Defaults to None.

        Returns:
            ax [axes]: Returns axes for plotting.
        """
        if ax == None:
            fig,ax = plt.subplots()
        ax.errorbar('x','y',yerr='yerr',data=self.results,
                    fmt = 'o', zorder = 0, **kwargs)
        if 'yfit' in self.results:
            ax.plot('x','yfit',data=self.results,
                    lw = 2, label='', zorder = 1,c='black')
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax

    def mlefit(self,func,pedantic = True, guess = [1.15E-5,0.002],*args,**kwargs):
        """
        Perform maximum likelihood estimation for optimal parameters.

        Args:
            func (func): Function to fit to.
            pedantic (bool, optional): Avoid annoying numpy warnings. Defaults to True.
            guess (list, optional): Initial parameter guesses. Defaults to [1.15E-5,0.002].
        """
        if pedantic == False:
            np.seterr('ignore')
        elif pedantic == True:
            np.seterr('warn')
        results = self.results
        x = results['x']
        e = results['shapes']
        y = results['y']
        MLEfit.__init__(self,func,x,e,y,guess=guess,*args,**kwargs)
        self.results['yfit'] = self.func(x,*self.params)

class AV(calibration):
    """
    Class for calibration with allan variance method.

    Args:
        calibration (class): Base class with shared calibration methods
    """
    def __init__(self, trace, fsample,taus = 'octave'):
        """
        Load trace and calculate allan variance.

        Args:
            trace (array): [description]
            fsample (float): Sample 
            taus (str, optional): Choose sampling type. Defaults to octave.
        """
        calibration.__init__(self, trace, fsample)
        if taus is None:
            taus = 'octave'
        taus,shapes,oavs = oavar(trace, rate = fsample, taus = taus)
        yerr = stats.gamma.std(shapes, scale = oavs/shapes)
        self.results = pd.DataFrame({'x':taus,'shapes':shapes,'y':oavs,'yerr':yerr})
    def plot(self,**kwargs):
        """
        Plot av data and fit (if applicable).

        Returns:
            **kwargs: Extra named parameters to send to errorbar.
        """
        ax = calibration.plot(self,**kwargs)
        ax.set_xlabel(r'$\tau$ (s)')
        ax.set_ylabel(r'$\sigma^2$ (nm$^2$)')
        plt.show()
        return ax
    def mlefit(self, fitfunc = 'SMMAV', *args, **kwargs):
        """
        Fit av data to analytical function.

        Args:
            fitfunc (func, optional): Funciton to fit to. Defaults to SMMAV.
        """
        if fitfunc is 'SMMAV':
            fitfunc = SMMAV
        calibration.mlefit(self, fitfunc,*args, **kwargs)

class PSD(calibration):
    """
    Class for PSD calibration.

    Args:
        calibration (base class): Base class with fitting methods
    """
    def __init__(self, trace, fsample,blocks=8,**kwargs):
        """
        Load trace and calculate PSD.

        Args:
            trace (array): Bead trajectory positions
            fsample (float): Camera sampling frequency
            blocks (int, optional): Number of blocks. Defaults to 8.
        """
        calibration.__init__(self, trace, fsample)
        f, dens = psd(trace, fsample, blocks = blocks)
        
        msk = f>0
        f, dens = f[msk], dens[msk]
        shapes = np.full_like(f,blocks)
        yerr = stats.gamma.std(shapes,scale = dens/shapes)
        self.results = pd.DataFrame({'x':f,'shapes':shapes,'y':dens,'yerr':yerr})
    def plot(self,**kwargs):
        """
        Plot PSD data and fit (if applicable).

        Returns:
            **kwargs: Extra named parameters to send to errorbar.
        """
        ax = calibration.plot(self,**kwargs)
        ax.set_xlabel(r'$f$ (Hz)')
        ax.set_ylabel(r'PSD (nm$^2$/Hz)')
        plt.show()
        return ax
    def mlefit(self, fitfunc = 'lansdorpPSD', *args, **kwargs):
        """
        Fit psd data to analytical function.

        Args:
            fitfunc (func, optional): Funciton to fit to. Defaults to lansdorpPSD.
        """
        if fitfunc == 'lansdorpPSD':
            fitfunc = lansdorpPSD
        elif fitfunc == 'aliasPSD':
            fitfunc = aliasPSD
        mlefunc = lambda x,a,k: fitfunc(x,self.fsample,a,k)
        calibration.mlefit(self, mlefunc,*args, **kwargs)

def m_generator(N,taus = None):
    if taus is None:
        taus = 'octave'
    if taus is 'all':
        # all-tau sampling not particularly useful but why not?
        m = np.linspace(1.0,N,N,dtype='int')
    elif taus is 'octave':
        # octave sampling break bin sizes, m, into powers of 2^n
        maxn = int(np.floor(np.log2(N//2))) # m =< N/2
        m = np.logspace(0,maxn-1,maxn,base=2,dtype='int')  #bin sizes
    elif taus is 'decade':
        # again not particularly useful, but why not?
        maxn = int(np.floor(np.log10(N)))
        m = [np.array([1,2,4])*k for k in np.logspace(0,maxn,maxn+1,base=10)]
        m = np.ravel(m)
    return m
def calc_totvar(x,rate,mj,mid,N):
    d0 = x[mid + 1:]
    d1 = x[mid + mj + 1:]
    d1n = x[mid - mj + 1:]
    e = min(len(d0), len(d1), len(d1n))
    v_arr = d1n[:e] - 2.0 * d0[:e] + d1[:e]
    var = np.sum(v_arr[:mid] * v_arr[:mid])
    var /= float(2 * pow(mj / rate, 2) * (N-2))
    return var
def calc_avar(phase,rate,mj,stride = 1):
    """
    Calculates the overlapping allan variance from
    Eq. 18b in erratum of Lansdorp et al. (2012)
    

    Args:
        phase ([type]): [description]
        rate ([type]): [description]
        mj ([type]): [description]
        stride (int, optional): [description]. Defaults to 1.

    Returns:
        [type]: [description]
    """
    stride = int(stride)
    d2 = phase[2 * mj:]
    d1 = phase[1 * mj:]
    d0 = phase
    # n = N - 2*m + 1
    # this handles all tau or octave
    n = min(len(d0), len(d1), len(d2))
    if n == 0:
        RuntimeWarning("Data array length is too small: %i" % len(phase))
        n = 1
    v_arr = d2[:n] - 2 * d1[:n] + d0[:n]
    s = np.sum(v_arr * v_arr)
    var =  s/(2.*n*(mj/rate)**2.)
    return var
def totvar(data, rate=1.0, taus=None):
    """ Total variance.
        Better confidence at long averages for Allan variance.
    .. math::
        \\sigma^2_{TOTDEV}( m\\tau_0 ) = { 1 \\over 2 (m\\tau_0)^2 (N-2) }
            \\sum_{i=2}^{N-1} ( {x}^*_{i-m} - 2x^*_{i} + x^*_{i+m} )^2
    where :math:`x^*_i` is a new time-series of length :math:`3N-4`
    derived from the original phase time-series :math:`x_n` of
    length :math:`N` by reflection at both ends.
    Parameters
    ----------
    data: np.array
        Bead positions.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement
    """
    data = np.asarray(data)
    N = len(data)
    m = m_generator(N,taus=taus)
    taus = m/rate
    shapes = (N/m-1)//2
    
    phase = np.cumsum(data)/rate  # convert frequency data to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    N = len(phase)
    
    # totdev requires a new dataset
    # Begin by adding reflected data before dataset
    x1 = 2.0 * phase[0] * np.ones((N - 2,))
    x1 = x1 - phase[1:-1]
    x1 = x1[::-1]

    # Reflected data at end of dataset
    x2 = 2.0 * phase[-1] * np.ones((N - 2,))
    x2 = x2 - phase[1:-1][::-1]

    # check length of new dataset
    assert len(x1)+len(phase)+len(x2) == 3*N - 4
    # Combine into a single array
    x = np.zeros((3*N - 4))
    x[0:N-2] = x1
    x[N-2:2*(N-2)+2] = phase  # original data in the middle
    x[2*(N-2)+2:] = x2
    
    mid = len(x1)
    calc_var = lambda mj: calc_totvar(x,rate,mj,mid,N)
    tvars = np.array([calc_var(mj) for mj in m])
    
    return taus,shapes,tvars
def oavar(data,rate = 1.0,taus = None):
    """
    Calculate overlapping allan variance 
    Takes an array of bead positions.
    Returns the taus, shapes, and oavs.

    .. math::
        \\z_j = \\tau_s\\sum_{i=1}^{j-1}{x_i}

        \\sigma^2_(m) = { 1 \\over 2 (m \\tau_c )^2 (N-2m+1) }
        \\sum_{n=1}^{N-2m+1} ( {z}_{k+2m} - 2z_{k+m} + z_{k} )^2\\
        z_j = \\tau_c \\sum_{i=1}^{j-1}x_i
        

    Parameters
    ----------
    data : array, series, or list
        1-D array of numbers.
    freq : float
        Frequency of acquisition.

    Returns
    -------
    taus : array
        taus.
    shape : array
        shapes.
    oavs : array
        Overlapping allan variance.

    """
    data = np.asarray(data) # convert to numpy array
    N = len(data)
    m = m_generator(N,taus = taus)
    taus = m/rate # tau = m*tau_c
    shapes = (N/m-1)//2 # shape factors
    # Calculate phasedata from Eq. 18b (in erratum)
    phase = np.cumsum(data)/rate  # integrate positions, converting frequency to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    # Calculate oav from Eq. 18a (in erratum)
    oavs = np.array([calc_avar(phase,rate,mj) for mj in m])
    return taus,shapes,oavs

def psd(trace,fsample,blocks):
    """
    Calculates the windowed and blocked experimental power spectral densities.

    Args:
        trace (array): Bead positions
        fsample (float): Sampling rate in Hz
        blocks (in): Number of blocks for binning

    Returns:
        f (array): Frequencies in Hz
        dens (array): Power spectral densities in nm^2/Hz
    """
    N = len(trace)
    nperseg = (2*N)/(blocks+1)
    f, dens = welch(trace, fsample, nperseg=nperseg,return_onesided=False)
    return f, dens

def SMMAV(t,a,k):
    """
    Analytical function for the AV of a trapped bead.
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
    Analytical function for the PSD of a trapped bead with aliasing and lowpass filtering.
    Eq. 7 in Lansdorp et al. (2012).

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

def loglikelihood(x,e,y,yhat):
    scale = yhat/e
    likelihood = gamma.pdf(y/scale,e)/scale
    return np.sum(np.log(likelihood))