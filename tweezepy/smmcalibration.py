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

import autograd.numpy as np
import corner
import emcee
import matplotlib.pyplot as plt
import warnings

from autograd import hessian
from autograd.scipy.stats import gamma
from inspect import signature
from scipy.optimize import curve_fit,minimize
from scipy.signal import welch
from scipy import stats

class MCMC:
    def __init__(self, walkers = 32, steps = 1600, progress = True, plot = True,**kwargs):
        """
        Run a Markov Chain Monte Carlo Sampler to determine possible asymmetric fitting uncertainty.

        Args:
            scale (list, optional): Parameter scale. Defaults to [1e-5,1e-3].
            N (int, optional): Number of samples. Defaults to 1500.
            **kwargs: Pass kwargs to the mcmc sampler
        """
        #self.params = params
        #self.LL = LL
        #self.walkers = walkers
        #self.steps = steps
        N = len(self.params)

        scale = np.power(10,np.floor(np.log10(self.params)))
        pos = self.params + 1e-4 * np.random.randn(walkers,N) * scale
        nwalkers,ndims = pos.shape
        self.sampler = emcee.EnsembleSampler(nwalkers,ndims,self.LL,**kwargs)
        self.sampler.run_mcmc(pos, steps,progress = progress)
        self.samples = self.sampler.get_chain()
        self.autocorr_time = self.sampler.get_autocorr_time()
        if plot:
            fig, axes = plt.subplots(N, figsize=(10, 3*N), sharex=True,squeeze=0)
            for i in range(N):
                ax = axes[i][0]
                ax.plot(self.samples[:, :, i], "k", alpha=0.3)
                ax.axvline(2*self.autocorr_time[i],c="k",lw=2)
                ax.set_xlim(0, len(self.samples))
                ax.set_ylabel(self.names[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)
            axes[-1][0].set_xlabel("step number")

    def calc_mc_errors(self, percentiles = [16,84], discard = 100, thin = 10):
        """
        Returns quantiles from monte carlo samples.

        Args:
            quantiles (list): Quantiles to calculate. Defaults to 16th and 84 pecentiles [16,50,84].
        """
        for tau in self.autocorr_time:
            if any(discard < 2*tau for tau in self.autocorr_time):
                warnings.warn('discard should be greater than twice the autocorrelation time %s.'%self.autocorr_time)
            if any(thin < tau%2 for tau in self.autocorr_time):
                warnings.warn('thin should be greater than half the autocorrelation time %s.'%self.autocorr_time)
        self.flat_samples = self.sampler.get_chain(discard=discard,thin = thin,flat = True)
        percentiles.insert(1,50)
        for i,name in enumerate(self.names):
            for p in percentiles:
                perc = np.percentile(self.flat_samples[:, i], (50,p))
                q = abs(np.diff(perc))
                self.results['%s|mcerr|%s'%(name,p)] = q

    def corner_plot(self,quantiles = (0.16,0.84),**kwargs):
        fig = corner.corner(self.flat_samples,  
                            truths=self.params,
                            quantiles=quantiles,
                            levels=(1-np.exp(-0.5)),
                            show_titles = True,
                            names = self.names,
                            **kwargs)

class MLEfit(MCMC):
    """
    Perform maximum likelihood estimation and uncertainty calculations.
    """
    def __init__(self,func,x,shape,y,guess = None,**kwargs):
        """
        Minimize negative log likelihood and find optimal parameters
        Args:
            func (func): Function to fit to.
            x (list or array): x data
            shape (list or array): shape parameter for gamma function
            y (list or array): y data
            guess (list, optional): initial parameter guesses.
        """
        self.func = func
        # Fancy way of determining fit param names        
        names = signature(func).parameters # inspect fit function parameters
        names = list(names.keys())  # make list of parameter names
        names = names[1:]; self.names = names # only take fit param names
        self.nparams = len(names)
        
        # Data
        self.x = np.asarray(x)
        self.shape = np.asarray(shape)
        self.y = np.asarray(y)
        
        self.ndata = len(y)
        logpdfs = lambda p: Gamma_Distribution(shape,func(x,*p)).logpdf(y)
        # Log likelihood
        #self.LL = lambda p: 0.5*loglikelihood(self.shape,self.y,func(self.x,*p))
        self.LL = lambda p: 0.5*np.sum(logpdfs(p))
        # Negative log likelihood
        self.negLL =  lambda p: -self.LL(p)
        # Use automatic differentiation to calculate hessian
        self.hess = hessian(self.negLL)
        # Minimize negative log likelihood
        self.fit = minimize(self.negLL,x0=guess,method='Nelder-Mead',**kwargs)
        # Save minimizer fit results
        self._params = self.fit['x']
        self.success = self.fit['success']
        # Collect results into dictionary
        self.results = {}
        self.results['success'] = self.success
        # Calculate information matrix from Hessian
        self.information_matrix = self.hess(self.params)
        if not self.success:
            warnings.warn('MLE fitting failed. %s'%self.fit['message'])
            self._params = [float('nan') for i in range(self.nparams)]
            self._std_errors = [float('nan') for i in range(self.nparams)]
        elif np.isnan(self.information_matrix).all():
            warnings.warn('Hessian method failed. Switching to MCMC uncertainty estimation.')
            self.mcmc()
            self._std_errors = [np.std(self.flat_samples[:,i]) for i in range(self.nparams)]
        else:
            # Covariance matrix
            inv_information_matrix = np.linalg.inv(self.information_matrix)
            #cov = 2 * inv_information_matrix*information_matrix*inv_information_matrix; self._cov = cov
            self._cov = inv_information_matrix
            # Calculate errors based on hessian
            self._std_errors = np.sqrt(np.diag(self.cov))
        for i,p in enumerate(names):
            self.results['%s'%p] = self.params[i]
            self.results['%s|err'%p] = self.std_errors[i]
        # Fit values
        self.yfit = func(self.x,*self.params)
        yerr = stats.gamma.std(self.shape, scale = self.yfit/self.shape)
        # Normalized residuals
        self.residuals = (self.y-self.yfit)/yerr
        # Calculate chi-squared based on likelihood ratio        
        LR = (loglikelihood(self.shape,self.y,self.yfit)-loglikelihood(self.shape,self.y,self.y))
        self._redchi2 = -2*LR/(self.ndata-self.nparams)
        self.results['redchi2'] = self.redchi2
        #results['AIC'] = self.AIC
        #results['BIC'] = self.BIC      

    def mcmc(self, walkers = 32, steps = 2000, percentiles = [16,84],**kwargs):
        MCMC.__init__(self,walkers=walkers,steps=steps,**kwargs)
        mc_errors = self.calc_mc_errors(percentiles=percentiles)

    @property
    def params(self):
        return self._params
    
    @params.setter
    def params(self, value):
        self._params = value

    @property
    def cov(self):
        return self._cov
    
    @property
    def std_errors(self):
        return self._std_errors

    @property
    def redchi2(self):
        return self._redchi2

    @property
    def AIC(self):
        return 2*(self.nparams-self.LL(self.params))

    @property
    def BIC(self):
        return self.nparams*np.log(self.ndata)-2*self.LL(self.params)
import math
from scipy import special
class Gamma_Distribution:
    def __init__(self,shape,yhat):
        self.yhat = yhat
        self.shape = shape
        self.scale = yhat/shape
    def pdf(self,y):
        shape = self.shape
        scale = self.scale
        return gamma.pdf(y/scale,shape)/scale
    def logpdf(self,y):
        shape = self.shape
        scale = self.scale
        return gamma.logpdf(y/scale,shape) - np.log(scale)
    def cdf(self,y):
        self.y = y
        return gamma.cdf(y,self.shape,scale=self.scale)
    def logcdf(self, y):
        return gamma.logcdf(y,self.shape,scale = self.scale)

def loglikelihood(shape,y,yhat):
    scale = yhat/shape
    likelihood = gamma.pdf(y/scale,shape)/scale
    loglikelihood = np.log(likelihood).sum()
    return loglikelihood
    #loglikelihood = gamma.logpdf(y,shape,0,scale).sum()
    #return loglikelihood

class calibration(MLEfit):
    """
    Base class for PSD and AV calibration.
    """
    def __init__(self,trace,fsample):
        """
        Loads in trace and fsample as class attributes.

        Args:
            trace (array): Bead trajectory positions in nm.
            fsample (array): Camera sample rate in Hz.
        """
        self.trace = trace
        self.fsample = fsample # in Hz

    def plot(self, ax = None,**kwargs):
        """
        Generate plot of data and fit to data (if applicable).

        Args:
            ax (axes, optional): If applicable, define axes to plot on. Defaults to None.

        Returns:
            ax [axes]: Returns axes for plotting.
        """
        if not ax:
            fig,ax = plt.subplots()
        data = self.data
        x = data['x']
        y = data['y']
        yerr = data['yerr']
        eb = ax.errorbar(x,y,yerr=yerr,
                        fmt = 'o', zorder = 0,
                        ms=3,capsize = 3,
                        **kwargs)
        if 'yfit' in data:
            yfit = data['yfit']
            c = eb[0].get_color()
            ax.plot(x,yfit,
                    lw = 2, label='',c=c, zorder = 1)
        ax.set_xscale('log')
        ax.set_yscale('log')
        return ax
    def guess(self,kT = 4.1, viscosity = 8.9e-4, radius = 0.5e-6):
        # Make initial guess for kappa based on equipartition theorem
        kappa = kT/np.var(self.trace)
        # Make initial guess for alpha based on myone bead in water 
        _alpha = 6 * np.pi * viscosity * radius * 1000
        guess = [_alpha,kappa]    
        return guess
    def mlefit(self, pedantic = True, **kwargs):
        """
        Perform maximum likelihood estimation for optimal parameters.

        Args:
            func (func): Function to fit to.
            guess (list, optional): Initial parameter guesses. Defaults to None.
            pedantic (bool, optional): Avoid annoying numpy warnings. Defaults to True.
            viscosity (float): Viscosity of water in Pa*s. Defaults to 8.9e-4 Pa*s.
            radius (float): Radius of bead in meters. Defaults to 0.5e-6 m.
            
        """
        if pedantic == False:
            np.seterr('warn')
        elif pedantic == True:
            np.seterr('ignore')

        # Biased nonlinear least squares fit to get good starting values
        popt,_ = curve_fit(self.func,self.data['x'],self.data['y'], sigma = self.data['yerr'],
                           absolute_sigma=True, p0 = self.guess, bounds = (1e-10,np.inf))
        MLEfit.__init__(self,self.func,self.data['x'],self.data['shape'],self.data['y'], guess=popt,**kwargs)
        self.data['yfit'] = self.yfit

class AV(calibration):
    """
    Class for calibration with allan variance method.

    Args:
        calibration (class): Base class with shared calibration methods
    """
    def __init__(self, trace, fsample,taus = 'octave',mode = 'oavar',**kwargs):
        """
        Load trace and calculate allan variance.

        Args:
            trace (array): [description]
            fsample (float): Sample 
            taus (str, optional): Choose sampling type. Defaults to octave.
        """
        calibration.__init__(self,trace,fsample)
        if mode == 'oavar':
            taus, shapes, oavs = oavar(trace, rate = fsample, taus = taus,**kwargs)
        elif mode == 'avar':
            taus, shapes, oavs = avar(trace, rate = fsample, taus = taus,**kwargs)
        elif mode == 'totvar':
            taus, shapes, oavs = totvar(trace, rate = fsample, taus = taus,**kwargs)
        yerr = stats.gamma.std(shapes, scale = oavs/shapes)
        self.data = {'x':taus,'shape':shapes,'y':oavs,'yerr':yerr}
        
    
    def plot(self,**kwargs):
        """
        Plot av data and fit (if applicable).

        Returns:
            **kwargs: Extra named parameters to send to errorbar.
        """
        ax = calibration.plot(self,**kwargs)
        ax.set_xlabel(r'$\tau$ (s)')
        ax.set_ylabel(r'$\sigma^2$ (nm$^2$)')
        return ax
    
    def mlefit(self, fitfunc = 'SMMAV', guess = None, alpha = None, pedantic = True,
               kT = 4.1, viscosity = 8.9e-4, radius = 0.5e-6,  **kwargs):
        """
        Fit av data to analytical function.

        Args:
            fitfunc (func, optional): Function to fit to. Defaults to SMMAV.
            guess (list, optional): Initial guess for parameters. Defaults to None.
            alpha (float, optional): Value for alpha, if alpha is fixed. Defaults to None.
            kT (float, optional): Value for thermal energy in pN nm. Defaults to 4.1. Defaults to
            viscosity (float, optional): Value for solution viscosity in Pa s. Only used if guess is None.
                                         Defaults to viscosity of water at 25 C, 8.9e-4.
            radius (float, optional): Value for bead radius in m. Only used if guess is None.
                                      Defaults to myone bead radius, 0.5e-6.
            pedantic (bool, optional): If True, turns off annoying numpy warnings. Defaults to False. 
            **kwargs
        """
                
        if not guess:
            guess = calibration.guess(self,kT=kT,viscosity = viscosity,radius = radius)
        if alpha:
            guess.pop(0)
        if fitfunc == 'SMMAV':
            func = SMMAV     
            if alpha:
                func = lambda x,k: SMMAV(x,alpha,k,kT)
            else:
                func = lambda x,a,k: SMMAV(x,a,k,kT)
                
        elif fitfunc == 'SMMAV2':
            func = SMMAV2
            if not guess:
                guess.append(1)
            if alpha:
                func = lambda x,k,b: SMMAV2(x,alpha,k,b,kT)
            else:
                func = lambda x,a,k,b: SMMAV2(x,a,k,b,kT)
        else:
            func = fitfunc
        self.guess = guess
        self.func = func
        calibration.mlefit(self, pedantic = pedantic, **kwargs)
        return self
        
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
        #N = len(trace)
        #blocks = (2*N/nperseg) - 1
        f, dens = psd(trace, fsample, blocks = blocks)
        msk = f>0
        f, dens = f[msk], dens[msk]
        shapes = np.full_like(f,blocks)
        yerr = stats.gamma.std(shapes,scale = dens/shapes)
        self.f,self.shapes,self.dens,self.yerr = f,shapes,dens, yerr
        self.data = {'x':f,'shape':shapes,'y':dens,'yerr':yerr}
    def plot(self,**kwargs):
        """
        Plot PSD data and fit (if applicable).

        Returns:
            **kwargs: Extra named parameters to send to errorbar.
        """
        ax = calibration.plot(self,**kwargs)
        ax.set_xlabel(r'$f$ (Hz)')
        ax.set_ylabel(r'PSD (nm$^2$/Hz)')
        return ax
    def mlefit(self, fitfunc = 'lansdorpPSD', guess = None, alpha = None, kT = 4.1, 
               viscosity = 8.9e-4, radius = 0.5e-6, pedantic = True, **kwargs):
        """
        Fit psd data to analytical function.

        Args:
            fitfunc (func, optional): Function to fit to. Defaults to lansdorpPSD.
            guess (list, optional): Initial guess for parameters. Defaults to None.
            alpha (float, optional): Value for alpha, if alpha is fixed. Defaults to None.
            kT (float, optional): Value for thermal energy in pN nm. Defaults to 4.1. Defaults to
            viscosity (float, optional): Value for solution viscosity in Pa s. Only used if guess is None.
                                         Defaults to viscosity of water at 25 C, 8.9e-4.
            radius (float, optional): Value for bead radius in m. Only used if guess is None.
                                      Defaults to myone bead radius, 0.5e-6.
            pedantic (bool, optional): If True, turns off annoying numpy warnings. Defaults to False. 
            **kwargs
        """
        if not guess:
            guess = calibration.guess(self,kT=kT,viscosity = viscosity,radius = radius)
        if alpha:
            guess.pop(0)
        # Determine fit function. Remove fsample variable for convenience.
        if fitfunc == 'lansdorpPSD':
            func = lambda x,a,k: lansdorpPSD(x,self.fsample,a,k,kT)
            if alpha:
                func = lambda x,k: func(x,alpha,k)
        elif fitfunc == 'lansdorpPSD2':
            func = lambda x,a,k,b: lansdorpPSD2(x,self.fsample,a,k,b,kT)
            if not guess:
                guess.append(1)
            if alpha:
                func = lambda x,k,b: func(x,alpha,k,b)
        elif fitfunc == 'aliasPSD':
            func = lambda x,a,k: aliasPSD(x,self.fsample,a,k,kT)
            if alpha:
                func = lambda x,k: aliasPSD(x,alpha,k)
        else:
            func = fitfunc
        self.guess = guess
        self.func = func
        calibration.mlefit(self, pedantic = pedantic, **kwargs)
        return self

        
def m_generator(N,taus = 'octave'):
    assert type(taus) is str, "taus must be a string"
    if taus == 'all':
        # all-tau sampling not particularly useful but why not?
        m = np.linspace(1.0,N,N,dtype='int')
        m = m[m<=N//2]
    elif taus == 'octave':
        # octave sampling break bin sizes, m, into powers of 2^n
        maxn = int(np.floor(np.log2(N/2))) # m =< N/2
        m = np.logspace(0,maxn-1,maxn,base=2,dtype='int')  #bin sizes
    elif taus == 'decade':
        # again not particularly useful, but why not?
        maxn = int(np.floor(np.log10(N)))
        m = [np.array([1,2,4])*k for k in np.logspace(0,maxn,maxn+1,base=10,dtype='int')]
        m = np.ravel(m)
        m = m[m<=N//2]
    return m
def avar(data,rate = 1.0,taus = 'octave'):
    """
    Calculate standard allan variance 
    Takes an array of bead positions.
    Returns the taus, shapes, and oavs.

    Parameters
    ----------
    data : array, series, or list
        1-D array of numbers.
    rate : float
        Frequency of acquisition.

    Returns
    -------
    taus : array
        taus.
    shape : array
        shapes.
    oavs : array
        Allan variance.

    """
    rate = float(rate)
    data = np.asarray(data) # convert to numpy array
    N = len(data)
    m = m_generator(N,taus = taus)
    n = N - 2*m + 1 
    m = m[n>=2]
    taus = m/rate # tau = m*tau_c

    # Calculate phasedata from Eq. 18b (in erratum)
    phase = np.cumsum(data)/rate  # integrate positions, converting frequency to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    assert len(phase) > 0, "Data array length is too small: %i" % len(phase)
    # Calculate oav from Eq. 18a (in erratum)
    oavs = np.zeros_like(taus)
    edfs = np.zeros_like(taus)
    for idx, mj in enumerate(m):
        d2 = phase[2 * mj::mj]
        d1 = phase[1 * mj::mj]
        d0 = phase[::mj]
        # n = N - 2*m + 1
        # this handles all tau or octave
        nmin = min(len(d0), len(d1), len(d2))
        edfs[idx] = approxedf(N,mj)
        v_arr = d2[:nmin] - 2 * d1[:nmin] + d0[:nmin]
        s = np.sum(v_arr * v_arr)
        var = np.divide(s,(2.*nmin*np.power(mj/rate,2)))
        oavs[idx] =  var
    shapes = edfs/2
    return taus,shapes,oavs

def oavar(data,rate = 1.0,taus = 'octave',edf = 'approx'):
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
    rate : float
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
    rate = float(rate)
    data = np.asarray(data) # convert to numpy array
    N = len(data)
    m = m_generator(N,taus = taus)
    n = N - 2*m + 1 
    m = m[n>=2]
    taus = m/rate # tau = m*tau_c

    # Calculate phasedata from Eq. 18b (in erratum)
    phase = np.cumsum(data)/rate  # integrate positions, converting frequency to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    assert len(phase) > 0, "Data array length is too small: %i" % len(phase)
    # Calculate oav from Eq. 18a (in erratum)
    oavs = np.zeros_like(taus)
    edfs = np.zeros_like(taus)
    for idx, mj in enumerate(m):
        d2 = phase[2 * mj:]
        d1 = phase[1 * mj:]
        d0 = phase
        # n = N - 2*m + 1
        # this handles all tau or octave
        nmin = min(len(d0), len(d1), len(d2))
        if edf == 'real':
            # autocorrelation is only good for N/mj >30
            if N/mj >= 32:
                alpha = noise_id(data,mj)[0]
            # use approximate edf for N/mj < 30
            else:
                alpha = 0
            edfs[idx] = edf_simple(N,mj,alpha)
        elif edf == 'approx':
            edfs[idx] = approxedf(N,mj)
        v_arr = d2[:nmin] - 2 * d1[:nmin] + d0[:nmin]
        s = np.sum(v_arr * v_arr)
        var = np.divide(s,(2.*nmin*np.power(mj/rate,2)))
        oavs[idx] =  var
    shapes = edfs/2
    return taus,shapes,oavs

def totvar(data, rate=1.0, taus='octave'):
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
    rate = float(rate)
    data = np.asarray(data) # make sure data is an array, not a series or list
    N = len(data)
    m = m_generator(N,taus=taus)
    #edfs = 3*N//(2*m)
    #shapes = edfs/2
    taus = m/rate
    #shapes = (N/m-1)//2
    
    phase = np.cumsum(data)/rate  # convert frequency data to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    N = len(phase)
    assert N > 0, "Data array length is too small: %i" % len(phase)
    
    # totvar requires a new dataset
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
    
    tvars = np.zeros_like(taus)
    shapes = np.zeros_like(taus)
    mid = len(x1)
    for idx,mj in enumerate(m):
        d0 = phase[mid + 1:]
        d1 = phase[mid + mj + 1:]
        d1n = phase[mid - mj + 1:]
        e = min(len(d0), len(d1), len(d1n))
        v_arr = d1n[:e] - 2.0 * d0[:e] + d1[:e]
        var = np.sum(v_arr[:mid] * v_arr[:mid])
        var /= float(2 * pow(mj / rate, 2) * (N-2))
        shapes[idx] = 3*N/(4*mj)
        tvars[idx] = var
    return taus,shapes,tvars

def noise_id(x,af, dmin = 0, dmax = 2):
    N = len(x)
    y_cut = np.array(x[:N-(N % af)])  # cut to length
    assert len(y_cut) % af == 0
    y_shaped = y_cut.reshape((int(len(y_cut)/af), af))
    x = np.average(y_shaped, axis=1)  # average
    
    
    # require minimum length for time-series
    if len(x) < 32:
        print(("noise_id() Can't determine noise-ID for"
               " time-series length= %d") % len(x))
        raise NotImplementedError
    d = 0 # number of differentiations
    while True:
        #c = np.corrcoef(x[:-af],x[af:])
        #r1 = c[0,1]
        r1 = autocorrelation(x)
        rho = r1/(1.0+r1)
        if d >= dmin and (rho < 0.25 or d >= dmax):
            p = -2.*(rho+d)
            alpha = p
            alpha_int = int(-1.0*np.round(2*rho) - 2.0*d)
            return alpha_int,alpha
        else:
            x = np.diff(x)
            d = d + 1
def edf_simple(N, m, alpha, pedantic = False):
    """Equivalent degrees of freedom.
    Simple approximate formulae.
    Parameters
    ----------
    N : int
        the number of phase samples
    m : int
        averaging factor, tau = m * tau0
    alpha: int
        exponent of f for the frequency PSD:
        'wp' returns white phase noise.             alpha=+2
        'wf' returns white frequency noise.         alpha= 0
        'fp' returns flicker phase noise.           alpha=+1
        'ff' returns flicker frequency noise.       alpha=-1
        'rf' returns random walk frequency noise.   alpha=-2
        If the input is not recognized, it defaults to idealized, uncorrelated
        noise with (N-1) degrees of freedom.
    Notes
    -----
       S. Stein, Frequency and Time - Their Measurement and
       Characterization. Precision Frequency Control Vol 2, 1985, pp 191-416.
       http://tf.boulder.nist.gov/general/pdf/666.pdf
       
       Modified from allantools.
    Returns
    -------
    edf : float
        Equivalent degrees of freedom
    """

    N = float(N)
    m = float(m)
    if alpha in [2, 1, 0, -1, -2]:
        # NIST SP 1065, Table 5
        if alpha == +2:
            edf = (N + 1) * (N - 2*m) / (2 * (N - m))

        if alpha == 0:
            edf = (((3 * (N - 1) / (2 * m)) - (2 * (N - 2) / N)) *
                   ((4*pow(m, 2)) / ((4*pow(m, 2)) + 5)))

        if alpha == 1:
            a = (N - 1)/(2 * m)
            b = (2 * m + 1) * (N - 1) / 4
            edf = np.exp(np.sqrt(np.log(a) * np.log(b)))

        if alpha == -1:
            if m == 1:
                edf = 2 * (N - 2)/(2.3 * N - 4.9)
            if m >= 2:
                edf = 5 * N**2 / (4 * m * (N + (3 * m)))

        if alpha == -2:
            a = (N - 2) / (m * (N - 3)**2)
            b = (N - 1)**2
            c = 3 * m * (N - 1)
            d = 4 * m**2
            edf = a * (b - c + d)

    else:
        edf = (N/m - 1) # assume correlated noise
        if pedantic == True:
            print("Noise type not recognized."
                  " Defaulting to N/m - 1 degrees of freedom.")

    return edf
def autocorrelation(x):
    n = len(x)
    variance = x.var()
    x = x-x.mean()
    r = np.correlate(x, x, mode = 'full')[-n:]
    result = r/(variance*(np.arange(n, 0, -1)))
    return result[1]


def approxedf(N,mj):
    N = float(N)
    mj = float(mj)
    edf = (N/mj - 1)
    return edf

def psd(trace,fsample,blocks=8):
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

def SMMAV(t,a,k,kT = 4.1):
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
    kT : float
        Thermal energy in pN nm. Default value is 4.1

    Returns
    -------
    oav : array
        predicted allan variance.

    """
    oav = 2.*kT*a/(np.power(k,2)*t) * (1. +
                               2.*a/(k*t) * np.exp(-k*t/a) -
                               a/(2.*k*t) * np.exp(-2.*k*t/a) -
                               3*a/(2.*k*t)
                              )
    return oav

def SMMAV2(t,a,k,b,kT = 4.1):
    """
    Same as SMMAV but with an extra white noise term.

    Parameters
    ----------
    t : array
        taus.
    alpha : float
        alpha.
    kappa : float
        kappa.
    kT : float
        Thermal energy in pN nm. Default value is 4.1

    Returns
    -------
    oav : array
        predicted allan variance.

    """
    kT = 4.1 # in pN*nm
    oav = SMMAV(t,a,k)+ b*a/t
    return oav

def lansdorpPSD(f,fs,a,k,kT = 4.1):
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
    kT : float
        Thermal energy in pN nm. Default value is 4.1

    Returns
    -------
    PSD : array
        power spectral density.
    """
    kT = 4.1 # themal energy in pN*nm
    PSD = 2.*kT*a/np.power(k,3) * (k + 
         (2.*a*fs*np.power(np.sin(np.pi*f/fs),2) * np.sinh(k/(a*fs))) / (np.cos(2.*np.pi*f/fs) - 
                                                                                   np.cosh(k/(a*fs)))
                          )
    return PSD
def lansdorpPSD2(f,fs,a,k,b,kT = 4.1):
    """
    Same as lansdorpPSD but with an extra white noise term.

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
    b : float
        white noise parameter
    kT : float
        Thermal energy in pN nm. Default value is 4.1

    Returns
    -------
    PSD : array
        power spectral density.
    """
    PSD = lansdorpPSD(f,fs,a,k) + b
    return PSD

def aliasPSD(f,fs,a,k, kT):
    """
    Alias only PSD.

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
        power spectral density.
    """
    kT = 4.1 # thermal energy in pNnm
    return kT/(k*fs) * (np.sinh(k/(a*fs))/(np.cosh(k/(a*fs))-np.cos(2*np.pi*f/fs)))



def mlefit(loglikelihood):
    negLL = lambda p: -loglikelihood(shape,y,yfit)
    fit = minimize(negLL,x0=popt,method='Nelder-Mead',**kwargs)
