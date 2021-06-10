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

from tweezepy.allanvar import avar,totvar
from tweezepy.expressions import aliasPSD,lansdorpPSD,SMMAV

class MCMC:
    def __init__(self, walkers = 32, steps = 1600, progress = True,**kwargs):
        """
        Monte Carlo sampler

        Parameters
        ----------
        walkers : int, optional
            Number of walkers, by default 32
        steps : int, optional
            Number of steps, by default 1600
        progress : bool, optional
            Print progress bar, by default True
        """
        self.walkers = walkers
        self.steps = steps

        scale = np.power(10,np.floor(np.log10(self.params)))
        pos = self.params + 1e-4 * np.random.randn(walkers,self.nparams) * scale
        nwalkers,ndims = pos.shape
        self.sampler = emcee.EnsembleSampler(nwalkers,ndims,self.logL,**kwargs)
        self.sampler.run_mcmc(pos, steps,progress = progress)
        self.samples = self.sampler.get_chain()
        self.autocorr_time = self.sampler.get_autocorr_time()

    def sample_plot(self,fig=None,labels = [],fig_kwgs={},ax_kwgs={}):
        if not isinstance(fig, plt.Figure):
            fig = plt.figure(figsize = (10,3*self.nparams),**fig_kwgs)
        axes = fig.get_axes()
        gs = plt.GridSpec(nrows=self.nparams,ncols =1)
        if len(labels) == 0:
            labels = self.names
        if len(axes) == 0:
            for i in range(self.nparams):
                axes.append(fig.add_subplot(gs[i]))
        #fig, axes = plt.subplots(self.nparams, figsize=(10, 3*self.nparams),
        #                         sharex=True,squeeze=0)
        for i,p in enumerate(labels):
            ax = axes[i]
            ax.plot(self.samples[:, :, i], c="k", alpha=0.3,**ax_kwgs)
            ax.axvline(2*self.autocorr_time[i],c="k",lw=2)
            ax.set_xlim(0, len(self.samples))
            ax.set_ylabel(p)
            ax.yaxis.set_label_coords(-0.1, 0.5)
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        axes[-1].set_xlabel("step number")
        return fig,axes

    def calc_mc_errors(self, percentiles = [15.87,50,84.13], discard = 100, thin = 10):
        """
        Computes percentiles from Monte Carlo samples.

        Parameters
        ----------
        percentiles : list, optional
            Percentiles for each parameter, by default [15.87,50,84.13]
        discard : int, optional
            Number of "burn-in" steps to discard, by default 100
        thin : int, optional
            N, by default 10
        """
        for tau in self.autocorr_time:
            if any(discard < 2*tau for tau in self.autocorr_time):
                warnings.warn('discard should be greater than twice the autocorrelation time %s.'%self.autocorr_time)
            if any(thin < tau%2 for tau in self.autocorr_time):
                warnings.warn('thin should be greater than half the autocorrelation time %s.'%self.autocorr_time)
        self.flat_samples = self.sampler.get_chain(discard=discard,thin = thin,flat = True)
        return np.percentile(self.flat_samples,percentiles,axis=0).T


    def corner_plot(self,quantiles = [0.16,0.84],labels = None,**kwargs):
        """
        Utility function for generating corner plots.

        Parameters
        ----------
        quantiles : list, optional
            Quantiles to annotate, by default (0.16,0.84)
        labels : list, optional
            Parameter labels, by default None

        Returns
        -------
        fig,ax
            Figure and axes objects.
        """
        if not labels:
            labels = self.names
        fig = corner.corner(self.flat_samples,  
                            truths=self.params,
                            quantiles=quantiles,
                            #levels=(1-np.exp(-0.5),),
                            #title_fmt = '.2e',
                            #show_titles = True,
                            labels = labels,
                            title_kwargs={'fontdict':{'fontsize':12}},
                            **kwargs)
        ax = fig.get_axes()
        for a in ax:
            if a.get_xlabel():
                a.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
            if a.get_ylabel():
                a.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                a.get_yaxis().get_offset_text().set_position((-.5,0.9))
        return fig,ax

class MLEfit(MCMC):
    """
    Perform maximum likelihood estimation and uncertainty calculations.
    """
    def __init__(self,guess = None, scale_covar = False,fit_kwargs = {'method':'Nelder-Mead'}):
        """
        Parameters
        ----------
        guess : list, optional
            Initial parameter guesses, by default None
        scale_covar : bool, optional
            Whether to scale standard errors by reduced chi-squared, by default False
        fit_kwargs : dict, optional
            Fitting keywork arguments to send to scipy.minimize, by default {'method':'Nelder-Mead'}
        """
        data = self.data
        shape,y,yerr = data['shape'],data['y'],data['yerr']
        self.nparams = len(self.names)
        self.ndata = len(y)
        self.nfree = self.ndata-self.nparams
        # Probability distribution function
        self.gd = Gamma_Distribution(shape,y)
        # Log likelihood
        self.logL = lambda p: self.gd.logpdf(self.func(*p)).sum()
        # Negative log likelihood
        self.negLL =  lambda p: -self.logL(p)
        # Use automatic differentiation to calculate hessian
        hess = hessian(self.negLL)
        # Minimize negative log likelihood
        self.fit = minimize(self.negLL,x0=guess,**fit_kwargs)
        # Save minimizer fit results
        self.params = self.fit['x']
        self.success = self.fit['success']
        # Collect results into dictionary
        self.results = {}
        # Throw warning if fit fails
        if not self.success:
            warnings.warn('MLE fitting failed. %s'%self.fit['message'])
            self.params = np.array([float('nan') for i in range(self.nparams)])
            self.std_errors = np.array([float('nan') for i in range(self.nparams)])
        # Compute standard errors by inverting the expected Hessian
        else:
            # Covariance matrix
            inv_hessian = np.linalg.inv(hess(self.params))
            self.cov = 2. * inv_hessian
            # Calculate errors from diagonals of covariance matrix
            self.std_errors = np.sqrt(np.diag(self.cov))
        for i,p in enumerate(self.names):
            self.results['%s'%p] = self.params[i]
            self.results['%s_error'%p] = self.std_errors[i]
        # Calculate loglikelihood
        self.loglikelihood = self.logL(self.params)
        # Calculate fit values
        yfit = self.func(*self.params); self.data['yfit'] = yfit
        # Calculate residuals, chi2, and reduced chi2
        residuals = (y-yfit)/yerr; self.data['residuals'] = residuals
        self.chi2 = np.power(residuals,2).sum(); self.results['chi2'] = self.chi2
        self.redchi2 = self.chi2/self.nfree; self.results['redchi2'] = self.redchi2
        # Scale errors by reduched chi-squared value
        if scale_covar:
            self.std_errors *= self.redchi2
        # Calculate fit support and p-value
        ks = stats.kstest(residuals,'chi2',args=(self.nfree,))
        self.support,self.p = ks
        self.results['support'] = self.support
        self.results['p-value'] = self.p
        # Calculate AIC, AICc, and BIC
        self.AIC = 2.*(self.nparams-self.loglikelihood); self.results['AIC'] = self.AIC
        self.AICc = self.AIC + (2*(pow(self.nparams,2)+self.nparams))/(self.ndata-self.nparams-1)
        #self.BIC = self.nparams*np.log(self.ndata)-2.*self.loglikelihood; self.results['BIC'] = self.BIC

    def mcmc(self, walkers = 32, steps = 2000, discard = 100, thin = 10):
        """
        Runs Monte Carlo sampler and computes standard errors as 0.5*(std_u - std_l)

        Parameters
        ----------
        walkers : int, optional
            Number of walkers, by default 32
        steps : int, optional
            Number of steps to take, by default 2000
        discard : int, optional
            Number of initial steps to discard, by default 100
        thin : int, optional
            Distance between independent steps, by default 10
        """
        MCMC.__init__(self,walkers=walkers,steps=steps)
        percentiles = self.calc_mc_errors(percentiles = [15.87,50,84.13],discard=discard,thin=thin)
        for i,name in enumerate(self.names):
            std_l, median, std_u = percentiles[i]
            self.results['%s'%name] = median
            self.std_errors[i] = 0.5 * (std_u - std_l)
            self.results['%s_error'%name] = 0.5 * (std_u - std_l)


class Gamma_Distribution:
    """
    Gamma probability distribution for AV and PSD.
    """
    def __init__(self,shape,yhat):
        """
        Parameters
        ----------
        shape : numeric
            [description]
        yhat : numeric
            Measured 
        """
        self.yhat = yhat
        self.shape = shape
    def scale(self,y):
        return y/self.shape 
    def std(self,y):
        scale = self.scale(y)
        return stats.gamma.std(self.shape, scale = scale)
    def pdf(self,y):
        scale = self.scale(y)
        return gamma.pdf(self.yhat/scale,self.shape)/scale
    def logpdf(self,y):
        scale = self.scale(y)
        return gamma.logpdf(self.yhat/scale,self.shape) - np.log(scale)
    def cdf(self,y):
        scale = self.scale(y)
        return stats.gamma.cdf(self.yhat,self.shape,scale=scale)
    def logcdf(self, y):
        scale = self.scale(y)
        return gamma.logcdf(self.yhat,self.shape,scale = scale)
    def interval(self,y,alpha = 0.95):
        scale = self.scale(y)
        return stats.gamma.interval(alpha,self.shape,scale = scale)

class calibration(MLEfit):
    """
    Base class for PSD and AV calibration.
    """
    def __init__(self,trace,fsample):
        """
        Parameters
        ----------
        trace : numpy.array
            1-D array of bead positions
        fsample : float
            Sampling frequency in Hz
        """
        self.trace = trace
        self.fsample = fsample

    def plot(self, 
             fig = None,
             fig_kwgs={},
             ax_fit_kwgs = {},
             ax_res_kwgs = {},
             label = None):
        """
        Utility function for plotting 

        Parameters
        ----------
        fig : Figure, optional
            Matplotlibe Figure object, by default None
        fig_kwgs : dict, optional
            Figure keyword arguments, by default {}
        ax_fit_kwgs : dict, optional
            Axis keyword arguments, by default {}
        ax_res_kwgs : dict, optional
            Axis keyword arguments, by default {}
        label : str, optional
            Legend label for data

        Returns
        -------
        fig,ax : Figure, Axes
            Figure and axes objects.
        """
        data = self.data
        if not isinstance(fig, plt.Figure):
            fig = plt.figure(figsize = (5,5),**fig_kwgs)
        ax = fig.get_axes()
        gs = plt.GridSpec(nrows=2, ncols=1, height_ratios=[4, 1],hspace = 0)
        if len(ax) == 0:
            ax.append(fig.add_subplot(gs[0], **ax_fit_kwgs))    
            if 'yfit' in data.keys():
                ax.append(fig.add_subplot(gs[1], **ax_res_kwgs))
                ax[1].axhline(0,c='k',**ax_res_kwgs)
                ax[1].set_ylabel(r'$\Delta$')
        data = self.data
        # Plot data with errorbars    
        eb = ax[0].errorbar(data['x'], data['y'], yerr=data['yerr'],
                            fmt = 'o', zorder = 0,
                            ms=3, capsize = 3, label = label)
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        # if there is a fit plot it
        if ('yfit' in data.keys()):
            # Plot fit
            c = eb[0].get_color() # get last color
            ax[0].plot(data['x'], data['yfit'],
                       lw = 2, label='',c=c, zorder = 1) 
            if (len(ax) == 2):
                ax[1].plot(data['x'],data['residuals'])
                ax[1].set_xscale('log')
                ax[0].get_xaxis().set_ticklabels([])       
        return fig,ax

    def make_guess(self, kT = 4.1, viscosity = 8.94e-10, radius = 530):
        """
        Estimate gamma based on Stokes drag and kappa based on equipartition theorem.

        Parameters
        ----------
        trace : np.array
            Bead positions in nm
        kT : float, optional
            Thermal energy in pN nm, by default 4.1
        viscosity : float, optional
            Dynamic viscosity of solution in pN s/nm^2, by default 8.94e-10
        radius : float, optional
            Bead radius in nm, by default 530

        Returns
        -------
        guess : list
            gamma and kappa guesses
        """
        viscosity = viscosity
        radius = radius
        # Make initial guess for kappa based on equipartition theorem
        kappa = kT/np.var(self.trace)
        # Make initial guess for alpha based on myone bead in water 
        gamma = 6 * np.pi * viscosity * radius 
        guess = [gamma,kappa]    
        return guess

    def mlefit(self, pedantic = True, **kwargs):
        """
        Estimate parameters and uncertainties using maximum likelihood estimation.

        Parameters
        ----------
        pedantic : bool, optional
            Ignore unhelpful warnings, by default True
        """
        if pedantic == False:
            np.seterr('warn')
        elif pedantic == True:
            np.seterr('ignore')
        # Fancy way of determining fit param names        
        names = signature(self.func).parameters # inspect fit function parameters
        names = list(names.keys())  # make list of parameter names
        self.names = names
        # Biased nonlinear least squares fit to get good starting values
        #try:
        #    popt,_ = curve_fit(self.func,self.x,self.y, sigma = self.yerr,
        #                       absolute_sigma=True, p0 = self.guess, bounds = (1e-10,np.inf))
        #except:
        #    popt = self.guess
        #    warnings.warn('Initial least-squares fit failed to find starting values.')
        MLEfit.__init__(self, guess=self.guess,**kwargs)

class AV(calibration):
    """
    Class for calibration with allan variance method.

    Args:
        calibration (class): Base class with shared calibration methods
    """
    def __init__(self, trace, fsample,taus = 'octave',mode = 'oavar',edf='approx'):
        """
        Load trace and calculate allan variance.

        Args:
            trace (array): [description]
            fsample (float): Sample 
            taus (str, optional): Choose sampling type. Defaults to octave.
        """
        calibration.__init__(self,trace,fsample)
        if mode == 'avar':
            taus, edfs, oavs = avar(trace, rate = fsample, taus = taus,edf=edf)
        elif mode == 'oavar':
            taus, edfs, oavs = avar(trace, rate = fsample, taus = taus,edf=edf,overlapping = True)
        elif mode == 'totvar':
            taus, edfs, oavs = totvar(trace, rate = fsample, taus = taus,edf=edf)
        else:
            raise ValueError('%s is not a valid mode.'%mode)
        shapes = edfs/2
        yerr = stats.gamma.std(shapes, scale = oavs/shapes)
        self.x = taus
        self.y = oavs
        self.yerr = yerr
        self.shape = shapes
        self.data = {'x':taus,'shape':shapes,'y':oavs,'yerr':yerr}
    
    def plot(self,
             fig = None,
             fig_kwgs = {},
             ax_fit_kwgs = {},
             ax_res_kwgs = {},
             label = None,
             **kwargs):
        """
        Utility function for plotting 

        Parameters
        ----------
        fig : Figure, optional
            Matplotlibe Figure object, by default None
        fig_kwgs : dict, optional
            Figure keyword arguments, by default {}
        ax_fit_kwgs : dict, optional
            Axis keyword arguments, by default {}
        ax_res_kwgs : dict, optional
            Axis keyword arguments, by default {}
        label : str or None, optional
            Data label for plotting multiple traces

        Returns
        -------
        fig,ax : Figure, Axes
            Figure and axes objects.
        """
        """
        Plot av data and fit (if applicable).

        Returns:
            **kwargs: Extra named parameters to send to errorbar.
        """
        fig,ax = calibration.plot(self,
                                  fig = fig,
                                  fig_kwgs=fig_kwgs,
                                  ax_fit_kwgs=ax_fit_kwgs,
                                  ax_res_kwgs=ax_res_kwgs,
                                  label = label)
        ax[-1].set_xlabel(r'$\tau$ (s)')
        ax[0].set_ylabel(r'$\sigma_{AV}^2$ (nm$^2$)')
        return fig,ax
    
    def mlefit(self, 
               fitfunc = 'SMMAV', 
               tracking_error = False, 
               guess = None, 
               gamma = None, 
               kappa = None, 
               epsilon = None,
               pedantic = True, 
               kT = 4.1, 
               viscosity = 8.94e-10, 
               radius = 530.,  
               **kwargs):
        """
        Estimate parameters and uncertainties via maximum likelihood estimation.

        Parameters
        ----------
        fitfunc : str or function, optional
            String or user defined function, by default 'SMMAV'
        tracking_error : bool, optional
            Incorporate tracking errors, by default False
        guess : list, optional
            Initial parameter guesses, by default None
        gamma : float, optional
            Fixed value for gamma parameter, by default None
        kappa : float, optional
            Fixed value for kappa parameter, by default None
        epsilon : float, optional
            Fixed value for epsilon parameter, by default None
        pedantic : bool, optional
            Ignore unhelpful warning messages, by default True
        kT : float, optional
            Thermal energy used in functions and initial parameter guesses, by default 4.1
        viscosity : float, optional
            Dynamic viscosity for initial parameter guesses in pN s/nm^2, by default 8.94e-10
        radius : float, optional
            Bead radius for inital parameter guesses, by default 530.

        Returns
        -------
        MLE fit results 
            Returns AV class object with MLE fit results
        """
        x = self.data['x']
        ts = 1/self.fsample
        # By default, use lansdorp SMM AV function
        if fitfunc == 'SMMAV':
            # Select appropriate function
            if gamma and kappa and tracking_error:
                func = lambda e: e/(2.*t)
            elif gamma and tracking_error:
                func = lambda k,e: SMMAV(x,gamma,k,kT=kT) + pow(e,2)*ts/x
            elif kappa and tracking_error:
                func = lambda g,e: SMMAV(x,g,kappa,kT=kT) + pow(e,2)*ts/x
            elif epsilon and tracking_error:
                func = lambda g,k: SMMAV(x,g,k,kT=kT) + pow(epsilon,2)*ts/x
            elif gamma:
                func = lambda k: SMMAV(x,gamma,k,kT=kT)
            elif kappa:
                func = lambda g: SMMAV(x,g,kappa,kT=kT)
            elif tracking_error:
                func = lambda g,k,e: SMMAV(x,g,k,kT=kT) + pow(e,2)*ts/x
            else:
                func = lambda g,k: SMMAV(x,g,k,kT=kT)
        else:
            # User provided function
            func = fitfunc
            if not guess:
                warnings.warn("User-provided function requires initial parameter guesses.")
        self.func = func
        # If no guess provided, make default guess based on Stoke's drag and equipartition theorem.
        if not guess:
            guess = calibration.make_guess(self,kT=kT, viscosity = viscosity, radius = radius)
            if gamma:
                guess.pop(0)
            if kappa:
                guess.pop(1)
            if tracking_error:
                guess.append(5)
        self.guess = guess
        calibration.mlefit(self, pedantic = pedantic, **kwargs)
        return self
        
class PSD(calibration):
    """
    Class for PSD calibration.

    Args:
        calibration (base class): Base class with fitting methods
    """
    def __init__(self, trace, fsample,bins=3,**kwargs):
        """
        Parameters
        ----------
        trace : np.array
            Bead positions in nm.
        fsample : float
            Sampling frequency in Hz
        bins : int, optional
            Number of bins, by default 3
        """
        calibration.__init__(self, trace, fsample)
        N = len(trace)
        nperseg = (2.*N)/(bins+1)
        f, dens = welch(trace, fsample, nperseg=nperseg,return_onesided=False)
        msk = f>0
        f, dens = f[msk], dens[msk]
        shapes = np.full_like(f,bins)
        yerr = stats.gamma.std(shapes,scale = dens/shapes)
        self.x = f
        self.shape = shapes
        self.y = dens
        self.yerr = yerr
        self.data = {'x':f,'shape':shapes,'y':dens,'yerr':yerr}

    def plot(self,
             fig = None,
             fig_kwgs = {},
             ax_fit_kwgs = {},
             ax_res_kwgs = {},
             label = None,
             **kwargs):
        """
        Utility function for plotting 

        Parameters
        ----------
        fig : Figure, optional
            Matplotlibe Figure object, by default None
        fig_kwgs : dict, optional
            Figure keyword arguments, by default {}
        ax_fit_kwgs : dict, optional
            Axis keyword arguments, by default {}
        ax_res_kwgs : dict, optional
            Axis keyword arguments, by default {}
        label : str or None, optional
            Data label for plotting multiple traces

        Returns
        -------
        fig,ax : Figure, Axes
            Figure and axes objects.
        """
        """
        Plot av data and fit (if applicable).

        Returns:
            **kwargs: Extra named parameters to send to errorbar.
        """
        fig,ax = calibration.plot(self,
                                  fig = fig,
                                  fig_kwgs=fig_kwgs,
                                  ax_fit_kwgs=ax_fit_kwgs,
                                  ax_res_kwgs=ax_res_kwgs,
                                  label = label)     
        ax[-1].set_xlabel(r'$f$ (Hz)')
        ax[0].set_ylabel(r'PSD (nm$^2$/Hz)')
        return fig,ax
    def mlefit(self, 
               fitfunc = 'lansdorpPSD', 
               tracking_error = False, 
               guess = None, 
               gamma = None, 
               kappa = None,
               epsilon = None, 
               kT = 4.1, 
               viscosity = 8.94e-10, 
               radius = 530,
               pedantic = True, **kwargs):
        """
        Estimate parameters and uncertainties via maximum likelihood estimation.

        Parameters
        ----------
        fitfunc : str or function, optional
            String or user defined function, by default 'lansdorpPSD'
        tracking_error : bool, optional
            Incorporate tracking errors, by default False
        guess : list, optional
            Initial parameter guesses, by default None
        gamma : float, optional
            Fixed value for gamma parameter, by default None
        kappa : float, optional
            Fixed value for kappa parameter, by default None
        epsilon : float, optional
            Fixed value for epsilon parameter, by default None
        pedantic : bool, optional
            Ignore unhelpful warning messages, by default True
        kT : float, optional
            Thermal energy used in functions and initial parameter guesses, by default 4.1
        viscosity : float, optional
            Dynamic viscosity for initial parameter guesses in pN s/nm^2, by default 8.94e-10
        radius : float, optional
            Bead radius for inital parameter guesses, by default 530.

        Returns
        -------
        MLE fit results 
            Returns AV class object with MLE fit results
        """
        fs = self.fsample
        x = self.data['x']
        # Select appropriate function
        if fitfunc == 'lansdorpPSD':
            if gamma and kappa and tracking_error:
                func = lambda e: e/fs
            elif gamma and tracking_error:
                func = lambda k,e: lansdorpPSD(x,fs,gamma,k,kT=kT) + pow(e,2)/fs
            elif kappa and tracking_error:
                func = lambda g,e: lansdorpPSD(x,fs,g,kappa,kT=kT) + pow(e,2)/fs
            elif epsilon and tracking_error:
                func = lambda g,k: lansdorpPSD(x,fs,g,k,kT=kT) + pow(epsilon,2)/fs
            elif gamma:
                func = lambda k: lansdorpPSD(x,fs,gamma,k,kT=kT)
            elif kappa:
                func = lambda g: lansdorpPSD(x,fs,g,kappa,kT=kT)
            elif tracking_error:
                func = lambda g,k,e: lansdorpPSD(x,fs,g,k,kT=kT) + pow(e,2)/fs
            else:
                func = lambda g,k: lansdorpPSD(x,fs,g,k,kT=kT)
        elif fitfunc == 'aliasPSD':
            if gamma and kappa and tracking_error:
                func = lambda e: e/fs
            elif gamma and tracking_error:
                func = lambda k,e: aliasPSD(x,fs,gamma,k,kT=kT) + pow(e,2)/fs
            elif kappa and tracking_error:
                func = lambda g,e: aliasPSD(x,fs,g,kappa,kT=kT) + pow(e,2)/fs
            elif epsilon and tracking_error:
                func = lambda g,k: aliasPSD(x,fs,g,k,kT=kT) + pow(epsilon,2)/fs
            elif gamma:
                func = lambda k: aliasPSD(x,fs,gamma,k,kT=kT)
            elif kappa:
                func = lambda g: aliasPSD(x,fs,g,kappa,kT=kT)
            elif tracking_error:
                func = lambda g,k,e: aliasPSD(x,fs,g,k,kT=kT) + pow(e,2)/(2.*t)
            else:
                func = lambda g,k: aliasPSD(x,fs,g,k,kT=kT)
        else:
            func = fitfunc
            if not guess:
                warnings.warn("User-provided functions require initial parameter guesses.")
        # If no guess provided, make default guess
        if not guess:
            guess = calibration.make_guess(self,kT=kT,viscosity = viscosity,radius = radius)
            if gamma:
                guess.pop(0)
            if kappa:
                guess.pop(1)
            if tracking_error:
                guess.append(10)
        self.guess = guess
        self.func = func
        calibration.mlefit(self, pedantic = pedantic, **kwargs)
        
        return self

