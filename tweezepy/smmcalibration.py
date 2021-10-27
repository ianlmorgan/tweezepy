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
import pkg_resources

from scipy.signal import welch
from tweezepy.allanvar import avar,totvar
from tweezepy.MLE import MLEfit, Gamma_Distribution
from tweezepy.expressions import aliasPSD,lansdorpPSD,SMMAV

def load_trajectory():
    """
    Loads test bead trajectory sampled at 400 Hz.

    Returns
    -------
    data : array-like
        Bead trajectory data in nm
    """
    fname = pkg_resources.resource_stream(__name__, 'data/trajectory.csv')
    data = np.loadtxt(fname)
    return data

class calibration(MLEfit):
    """
    Base class for PSD and AV calibration.
    
    Parameters
    ----------
    trace : array
        1-D array of bead positions in nm
    fsample : float
        Sampling frequency in Hz

    Attributes
    ----------
    trace : array
        1-D array of bead positions in nm
    fsample : float
        Sampling frequency in Hz
    """
    def __init__(self,trace,fsample):
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
        Utility function for plotting.

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
        fig, ax : Figure, Axes
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
        if isinstance(self,AV):
            ax[-1].set_xlabel(r'$\tau$ (s)')
            ax[0].set_ylabel(r'$\sigma_{AV}^2$ (nm$^2$)')
        if isinstance(self,PSD):
            ax[-1].set_xlabel(r'$f$ (Hz)')
            ax[0].set_ylabel(r'PSD (nm$^2$/Hz)')   
        return fig, ax
    def _make_guess(self, 
                    kT = 4.1,
                    viscosity = 8.94e-10,
                    radius = 530):
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
        # Make initial guess for kappa based on equipartition theorem
        kappa = kT/np.var(self.trace)
        # Make initial guess for alpha based on myone bead in water 
        gamma = 6 * np.pi * viscosity * radius 
        guess = [gamma,kappa]    
        return guess
    def _predefined(self,
                    func,
                    gamma = None,
                    kappa = None,
                    tracking_error = False,
                    epsilon = None,
                    ):
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
    def mlefit(self,
               fitfunc = None, 
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
            String or user defined function, 
            by default 'lansdorpPSD' for PSD and 'SMMAV' for AV
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
            Thermal energy in pN nm, by default 4.1
        viscosity : float, optional
            Dynamic viscosity for initial parameter guesses in pN s/nm^2, by default 8.94e-10
        radius : float, optional
            Bead radius for inital parameter guesses in nm, by default 530.

        """
        if cutoffs:
            assert len(cutoffs) == 2, "cutoffs should have a lower an upper cutoff in the form of [lower,upper]."
        msk = (self.data_init['x']>=cutoffs[0])&(self.data_init['x']<=cutoffs[1])
        for key,d in self.data_init.items():
            self.data[key] = d[msk]        
        x = self.data['x']

        if isinstance(self,AV):
            if not fitfunc:
                fitfunc = 'SMMAV'
            if fitfunc == 'SMMAV':
                ts = 1/self.fsample
                func = lambda g,k,e: SMMAV(x,ts,g,k,e,kT=kT)
                func = self._predefined(func,
                                        gamma=gamma,
                                        kappa = kappa,
                                        tracking_error=tracking_error,
                                        epsilon = epsilon,
                                        )
            else:
                func = fitfunc
                if not guess:
                    raise RuntimeError("User-provided function requires initial parameter guesses.")
            
        if isinstance(self,PSD):
            # Select appropriate function
            if not fitfunc:
                fitfunc = 'lansdorpPSD'

            if fitfunc == 'lansdorpPSD':
                func = lambda g,k,e: lansdorpPSD(x,self.fsample,g,k,e,kT=kT)
                func = self._predefined(func,
                                        gamma=gamma,
                                        kappa = kappa,
                                        tracking_error=tracking_error,
                                        epsilon = epsilon,
                                        )
            elif fitfunc == 'aliasPSD':
                func = lambda g,k,e: aliasPSD(x,self.fsample,g,k,e,kT=kT)
                func = self._predefined(func,
                                        gamma=gamma,
                                        kappa = kappa,
                                        tracking_error=tracking_error,
                                        epsilon = epsilon,
                                        )
            else:
                func = fitfunc
                if not guess:
                    raise RuntimeError("User-provided function requires initial parameter guesses.")
        self.func = func
        # If no guess provided, make default guess
        if not guess:
            #make default guess based on Stoke's drag and equipartition theorem
            guess = self._make_guess(kT=kT,
                                     viscosity=viscosity,
                                     radius=radius)
            if gamma:
                guess.pop(0)
            if kappa:
                guess.pop(1)
            if tracking_error:
                guess.append(10)
        self.guess = guess
        MLEfit.__init__(self, pedantic = pedantic, scale_covar = scale_covar, **kwargs)
    @property
    def data(self):
        """
        Returns dictionary of data.

        Returns
        -------
        dict
            Dictionary of data.
        """
        return self._data
    @data.setter
    def data(self, value):
        self._data = value
class AV(calibration):
    """
    A class for computing and fitting the Allan variance using maximum likelihood estimation.

    Parameters
    ----------
    trace : array
        Bead trajectory in nm.
    fsample : float
        Sampling frequency in Hz.
    taus : str, optional
        Tau interval sampling. Options are 'octave', 'decade', and 'all', by default 'octave'
    mode : str, optional
        Allan variance type, either 'avar', 'oavar', and 'totvar', by default 'oavar'. 
        See allanvar module for more information.
    edf : str, optional
        Equivalent degrees of freedom for AV. Options are 'real' and 'approx', by default 'real'

    Raises
    ------
    AssertionError
        If taus, mode, or edf is not recognized.

    Examples
    --------
    >>> from tweezepy import load_trajectory, AV
    >>> trace = load_trajectory() # Load trajectory in nm
    >>> av = AV(trace, fsample = 400) # Compute overlapping AV
    >>> av.mlefit() # Perform MLE fit
    >>> print(av.results)
    """

    def __init__(self, trace, fsample,taus = 'octave',mode = 'oavar',edf='real'):
        calibration.__init__(self,trace,fsample)
        assert taus in ['all','octave','decade'], "taus must be either all, octave, or decade."
        assert mode in ['avar','oavar','totvar'], "mode must be either avar, oavar, or totvar."
        assert edf in ['real','approx'], "edf must be either real or approx."
        if mode == 'avar':
            taus, edfs, oavs = avar(trace, rate = fsample, taus = taus,edf=edf)
        elif mode == 'oavar':
            taus, edfs, oavs = avar(trace, rate = fsample, taus = taus,edf=edf,overlapping = True)
        elif mode == 'totvar':
            taus, edfs, oavs = totvar(trace, rate = fsample, taus = taus,edf=edf)
        self.x = taus
        self.y = oavs
        self.shape = edfs/2
        # Probability distribution function
        self.gd = Gamma_Distribution(self.shape,self.y)
        self.yerr = self.gd.std(self.y)
        self._data = {'x':self.x,
                      'shape':self.shape,
                      'y':self.y,
                      'yerr':self.yerr}
        self.data_init = self.data.copy()
        
class PSD(calibration,MLEfit):
    """
    A class for computing and fitting the power spectral density using maximum likelihood estimation.

    Parameters
    ----------
    trace : array
        Bead positions in nm.
    fsample : float
        Sampling frequency in Hz
    bins : int, optional
        Number of bins, by default 3

    Example
    -------
    >>> from tweezepy import load_trajectory, PSD            
    >>> trace = load_trajectory() # Load trajectory in nm
    >>> psd = PSD(trace, fsample = 400) # Compute PSD using Welch's method
    >>> psd.mlefit() # Perform MLE fit
    >>> print(psd.results)
    """

    def __init__(self, 
                 trace, # input bead trajectory 
                 fsample, # input sampling frequency in Hz
                 bins=3, # input number of half-overlapping bins
                 ):
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
        self._data = {'x':self.x,
                      'shape':self.shape,
                      'y':self.y,
                      'yerr':self.yerr}
        self.data_init = self.data.copy()
    
    #def mlefit(self, 
    #           fitfunc = 'lansdorpPSD', 
    #           cutoffs = [-np.inf,np.inf],
    #           tracking_error = False, 
    #           guess = None, 
    #           gamma = None, 
    #           kappa = None, 
    #           epsilon = None,
    #           kT = 4.1, 
    #           viscosity = 8.94e-10, 
    #           radius = 530.,  
    #           pedantic = True, 
    #           scale_covar = False,
    #           **kwargs):
        
    #    self.fitfunc = fitfunc
    #    self.cutoffs = cutoffs
    #    self.tracking_error = tracking_error
    #    self.guess = guess
    #    self.gamma = gamma
    #    self.kappa = kappa
    #    self.epsilon = epsilon
    #    self.kT = kT
    #    self.viscosity = viscosity
    #    self.radius = radius
    #    self.pedantic = pedantic
    #    self.scale_covar = scale_covar
    #    
    #    
    #    calibration.mlefit(self, **kwargs)

