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
import numpy as np

from scipy.signal import welch
from tweezepy.allanvar import avar,totvar
from tweezepy.MLE import MLEfit, Gamma_Distribution
from tweezepy.expressions import aliasPSD,lansdorpPSD,SMMAV

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
             data_label = None,
             fit_label = None,
             data_color = None,
             fit_color = 'k'):
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
        data_label : str, optional
            Legend label for data.
        fit_label : str, optional
            Legend label for fit line.

        Returns
        -------
        fig,ax : Figure, Axes
            Figure and axes objects.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise RuntimeError("Matplotlib is required for plotting.")
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
        # Plot data with errorbars    
        eb = ax[0].errorbar(data['x'], data['y'], yerr=data['yerr'],
                            fmt = 'o', zorder = 0, c = data_color,
                            ms=3, capsize = 3, label = data_label)
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        # if there is a fit plot it
        if ('yfit' in data.keys()):
            # Plot fit
            c = eb[0].get_color() # get last color
            ax[0].plot(data['x'], data['yfit'],
                       lw = 2, label=fit_label,c=fit_color, zorder = 1) 
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
    def predefined(func,
                   gamma = None,
                   kappa = None,
                   tracking_error = False,
                   epsilon = None):
        if gamma and kappa and tracking_error:
             func2 = lambda e: func(gamma, kappa, e)
        elif gamma and tracking_error:
            func2 = lambda k,e: func(gamma, k, e)
        elif kappa and tracking_error:
            func2 = lambda g,e: func(g, kappa, e)
        elif epsilon and tracking_error:
            func2 = lambda g,k: func(g, k, epsilon)
        elif gamma:
            func2 = lambda k: func(gamma,k,0)
        elif kappa:
            func2 = lambda g: func(g,kappa,0)
        elif tracking_error:
            func2 = lambda g,k,e: func(g,k,e)
        else:
            func2 = lambda g,k: func(g,k,0)
        return func2
    def mlefit(self,**kwargs):
        # If no guess provided, make default guess
        if not self.guess:
            #make default guess based on Stoke's drag and equipartition theorem
            guess = self.make_guess()
            if self.gamma:
                guess.pop(0)
            if self.kappa:
                guess.pop(1)
            if self.tracking_error:
                guess.append(10)
        self.guess = guess
        MLEfit.__init__(self, pedantic = self.pedantic, scale_covar = self.scale_covar, **kwargs)


class AV(calibration):
    """
    A class for computing and fitting the Allan variance using MLE.

    Parameters
    ----------
    calibration : class
        Base class with utility functions.
    """
    def __init__(self, trace, fsample,taus = 'octave',mode = 'oavar',edf='real'):
        """
        Computes AV values.

        Parameters
        ----------
        trace : array-like
            Bead trajectory.
        fsample : float
            Sampling frequency.
        taus : str, optional
            Tau interval sampling, by default 'octave'
        mode : str, optional
            Allan variance type, either 'avar', 'oavar', and 'totvar', by default 'oavar'
        edf : str, optional
            Equivalent degrees of freedom for AV, by default 'real'

        Raises
        ------
        ValueError
            [description]
        """
        calibration.__init__(self,trace,fsample)
        if not mode in ['avar','oavar','totvar']:
            raise RuntimeError('taus should be avar, oavar, or totvar.')
        if mode == 'avar':
            taus, edfs, oavs = avar(trace, rate = fsample, taus = taus,edf=edf)
        elif mode == 'oavar':
            taus, edfs, oavs = avar(trace, rate = fsample, taus = taus,edf=edf,overlapping = True)
        elif mode == 'totvar':
            taus, edfs, oavs = totvar(trace, rate = fsample, taus = taus,edf=edf)
        else:
            raise ValueError('%s is not a valid mode.'%mode)
        self.x = taus
        self.y = oavs
        self.shape = edfs/2
        # Probability distribution function
        self.gd = Gamma_Distribution(self.shape,self.y)
        self.yerr = self.gd.std(self.y)
        self.data = {'x':self.x,
                     'shape':self.shape,
                     'y':self.y,
                     'yerr':self.yerr}
        self.data_init = self.data.copy()
    
    def plot(self,
             fig = None,
             fig_kwgs = {},
             ax_fit_kwgs = {},
             ax_res_kwgs = {},
             data_label = None,
             fit_label = None,
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
        data_label : str or None, optional
            Data label for plotting multiple traces
        fit_label : str or None, optional
            Fit label for plotting multiple traces

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
        fig,ax = calibration.plot(self)
        ax[-1].set_xlabel(r'$\tau$ (s)')
        ax[0].set_ylabel(r'$\sigma_{AV}^2$ (nm$^2$)')
        return fig,ax
    
    def mlefit(self, 
               fitfunc = 'SMMAV', 
               cutoffs = [-np.inf,np.inf],
               tracking_error = False, 
               guess = None, 
               gamma = None, 
               kappa = None, 
               epsilon = None,
               kT = 4.1, 
               viscosity = 8.94e-10, 
               radius = 530.,  
               pedantic = True, 
               scale_covar = False,
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
        kT : float, optional
            Thermal energy used in functions and initial parameter guesses, by default 4.1
        viscosity : float, optional
            Dynamic viscosity for initial parameter guesses in pN s/nm^2, by default 8.94e-10
        radius : float, optional
            Bead radius for inital parameter guesses, by default 530.
        pedantic : bool, optional
            Ignore unhelpful warning messages, by default True
        scale_covar: bool, optional
            Whether to scale covariance by reduced chi-squared, by default False
        **kwargs
            keyword arguments passed to scipy's minimizer

        Returns
        -------
        MLE fit results 
            Returns AV class object with MLE fit results
        """
        self.fitfunc = fitfunc
        self.cutoffs = cutoffs
        self.tracking_error = tracking_error
        self.guess = guess
        self.gamma = gamma
        self.kappa = kappa
        self.epsilon = epsilon
        self.kT = kT
        self.viscosity = viscosity
        self.radius = radius
        self.pedantic = pedantic
        self.scale_covar = scale_covar
        if cutoffs:
            assert len(cutoffs) == 2, "cutoffs should have a lower an upper cutoff in the form of [lower,upper]."
            msk = (self.data_init['x']>=cutoffs[0])&(self.data_init['x']<=cutoffs[1])
            for key,d in self.data_init.items():
                self.data[key] = d[msk]        
        x = self.data['x']
        # Sampling time
        ts = 1/self.fsample
        # By default, use lansdorp SMM AV function
        if fitfunc == 'SMMAV':
            func = lambda g,k,e: SMMAV(x,ts,g,k,e,kT=kT)
            func = calibration.predefined(func,
                                          gamma=gamma,
                                          kappa = kappa,
                                          tracking_error=tracking_error,
                                          epsilon = epsilon)
        else:
            # User provided function
            func = fitfunc
            if not guess:
                raise RuntimeError("User-provided function requires initial parameter guesses.")
        self.func = func
        calibration.mlefit(self,**kwargs)
        
class PSD(calibration,MLEfit):
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
        self.x = f
        self.shape = shapes
        self.y = dens
        # Probability distribution function
        self.gd = Gamma_Distribution(self.shape,self.y)
        self.yerr = self.gd.std(self.y)
        self.data = {'x':self.x,
                     'shape':self.shape,
                     'y':self.y,
                     'yerr':self.yerr}
        self.data_init = self.data.copy()

    def plot(self,
             fig = None,
             fig_kwgs = {},
             ax_fit_kwgs = {},
             ax_res_kwgs = {},
             data_label = None,
             fit_label = None,
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
        fig,ax = calibration.plot(self)     
        ax[-1].set_xlabel(r'$f$ (Hz)')
        ax[0].set_ylabel(r'PSD (nm$^2$/Hz)')
        return fig,ax
    def mlefit(self, 
               fitfunc = 'lansdorpPSD', 
               cutoffs = [-np.inf,np.inf],
               tracking_error = False, 
               guess = None, 
               gamma = None, 
               kappa = None, 
               epsilon = None,
               kT = 4.1, 
               viscosity = 8.94e-10, 
               radius = 530.,  
               pedantic = True, 
               scale_covar = False,
               **kwargs):
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
        self.cutoffs = cutoffs
        self.tracking_error = tracking_error
        self.guess = guess
        self.gamma = gamma
        self.kappa = kappa
        self.epsilon = epsilon
        self.kT = kT
        self.viscosity = viscosity
        self.radius = radius
        self.pedantic = pedantic
        self.scale_covar = scale_covar
        if cutoffs:
            assert len(cutoffs) == 2, "cutoffs should have a lower an upper cutoff in the form of [lower,upper]."
            msk = (self.data_init['x']>=cutoffs[0])&(self.data_init['x']<=cutoffs[1])
            for key,d in self.data_init.items():
                self.data[key] = d[msk]        
        x = self.data['x']
        # Select appropriate function
        if fitfunc == 'lansdorpPSD':
            func = lambda g,k,e: lansdorpPSD(x,self.fsample,g,k,e,kT=kT)
            func = calibration.predefined(func,
                                          gamma=gamma,
                                          kappa = kappa,
                                          tracking_error=tracking_error,
                                          epsilon = epsilon)
        elif fitfunc == 'aliasPSD':
            func = lambda g,k,e: aliasPSD(x,self.fsample,g,k,e,kT=kT)
            func = calibration.predefined(func,
                                          gamma=gamma,
                                          kappa = kappa,
                                          tracking_error=tracking_error,
                                          epsilon = epsilon)
        else:
            func = fitfunc
            if not guess:
                raise RuntimeError("User-provided function requires initial parameter guesses.")
        self.func = func
        calibration.mlefit(self, **kwargs)

