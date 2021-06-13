"""
This file is part of tweezepy.py
"""

from autograd import hessian
from autograd.scipy.stats import gamma
from inspect import signature
from scipy.optimize import minimize
from scipy import stats

import autograd.numpy as np

class MCMC:
    """
    Monte Carlo sampler class.

    :Example:
        ::

            from tweezepy import AV,downsampled_trace
            trace = downsampled_trace
            av = AV(trace, 100)
            av.mlefit()
            av.mcmc()
            av.sample_plot()
            av.corner_plot()
    
    """
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
        try:
            import emcee
        except ImportError:
            RuntimeError("Monte Carlo sampling requires the emcee package.")
        self.walkers = walkers
        self.steps = steps

        scale = np.power(10,np.floor(np.log10(self.params)))
        pos = self.params + 1e-4 * np.random.randn(walkers,self.nparams) * scale
        nwalkers,ndims = pos.shape
        self.sampler = emcee.EnsembleSampler(nwalkers,ndims,self.logL,**kwargs)
        self.sampler.run_mcmc(pos, steps,progress = progress)
        self.samples = self.sampler.get_chain()
        self.autocorr_time = self.sampler.get_autocorr_time()

    def sample_plot(self,
                    fig=None,
                    labels = [],
                    fig_kwgs={},
                    ax_kwgs={}):
        """
        Plot the accepted Monte Carlo samples.

        Parameters
        ----------
        fig : object, optional
            Figure object, by default None
        labels : list, optional
            Plot labels, by default []
        fig_kwgs : dict, optional
            Figure keywords, by default {}
        ax_kwgs : dict, optional
            Axes keywords, by default {}

        Returns
        -------
        (fig, axes) : tuple
            Figure and axes objects.
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise RuntimeError("Matplotlib is required for plotting.")
        if not isinstance(fig, plt.Figure):
            fig = plt.figure(figsize = (10,3*self.nparams),**fig_kwgs)
        axes = fig.get_axes()
        gs = plt.GridSpec(nrows=self.nparams,ncols =1)
        if len(labels) == 0:
            labels = self.names
        if len(axes) == 0:
            for i in range(self.nparams):
                axes.append(fig.add_subplot(gs[i]))
        for i,p in enumerate(labels):
            ax = axes[i]
            ax.plot(self.samples[:, :, i], c="k", alpha=0.3,**ax_kwgs)
            ax.axvline(2*self.autocorr_time[i],c="k",lw=2)
            ax.set_xlim(0, len(self.samples))
            ax.set_ylabel(p)
            ax.yaxis.set_label_coords(-0.1, 0.5)
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
        axes[-1].set_xlabel("step number")
        return fig, axes

    def calc_mc_errors(self, percentiles = [15.87,50,84.13], discard = None, thin = None):
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
        tau = max(self.autocorr_time)

        if not discard:
            discard = np.ceil(2*tau)
        if not thin:
            thin = np.ceil(tau/2)
        if discard < 2*tau:
            raise RuntimeError('discard should be greater than twice the autocorrelation time %s.'%self.autocorr_time)
        if thin < tau%2:
            raise RuntimeError('thin should be greater than half the autocorrelation time %s.'%self.autocorr_time)
        self.flat_samples = self.sampler.get_chain(discard=discard,thin = thin,flat = True)
        return np.percentile(self.flat_samples,percentiles,axis=0).T

    def corner_plot(self,
                    quantiles = [0.16,0.84],
                    labels = None,
                    **kwargs):
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
        try:
            import corner
        except ImportError:
            raise RuntimeError("Corner plot requires the corner module.")
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
    def __init__(self, pedantic = True, scale_covar = False,fit_kwargs = {'method':'Nelder-Mead'}):
        """
        Parameters
        ----------
        pedantic : bool, optional
            Ignore unhelpful warnings, by default True
        scale_covar : bool, optional
            Whether to scale standard errors by reduced chi-squared, by default False
        fit_kwargs : dict, optional
            Fitting keywork arguments to send to scipy.minimize, by default {'method':'Nelder-Mead'}
        """
        if pedantic == False:
            np.seterr('warn')
        elif pedantic == True:
            np.seterr('ignore')
        # Fancy way of determining fit param names        
        names = signature(self.func).parameters # inspect fit function parameters
        self.names = list(names.keys()); # make list of parameter names
        # Data
        data = self.data
        shape,y,yerr = data['shape'],data['y'],data['yerr']
        self.ndata = len(y)
        # Log likelihood
        self.logL = lambda p: self.gd.logpdf(self.func(*p)).sum()
        # Negative log likelihood
        self.negLL =  lambda p: -self.logL(p)
        # Use automatic differentiation to calculate hessian
        hess = hessian(self.negLL)
        # Minimize negative log likelihood
        self.fit = minimize(self.negLL,x0=self.guess,**fit_kwargs)
        # Save minimizer fit results
        self.params = self.fit['x']
        self.nparams = len(self.fit['x'])
        self.success = self.fit['success']
        # Collect results into dictionary
        self.results = {}
        # Throw warning if fit fails
        if not self.success:
            print('MLE fitting failed. %s'%self.fit['message'])
            #self.params = np.array([float('nan') for i in range(self.nparams)])
            #self.std_errors = np.array([float('nan') for i in range(self.nparams)])
        # Compute standard errors by inverting the Hessian
        inv_hessian = np.linalg.inv(hess(self.params))
        # Covariance matrix
        self.cov = 2. * inv_hessian
        # Calculate errors from diagonals of covariance matrix
        self.std_errors = np.sqrt(np.diag(self.cov))
        # Calculate fit values
        yfit = self.func(*self.params); 
        self.fit_data = data.copy()
        self.data['yfit'] = yfit
        # Calculate residuals, chi2, and reduced chi2
        residuals = (y-yfit)/yerr; self.residuals = residuals
        self.data['residuals'] = residuals
        self.chi2 = np.power(residuals,2).sum(); self.results['chi2'] = self.chi2
        self.nfree = self.ndata-self.nparams
        self.redchi2 = self.chi2/self.nfree; self.results['redchi2'] = self.redchi2
        # Scale errors by reduched chi-squared value
        if scale_covar:
            self.std_errors *= self.redchi2
        # Collect errors into results
        for i,p in enumerate(self.names):
            self.results['%s'%p] = self.params[i]
            self.results['%s_error'%p] = self.std_errors[i]
        # Calculate fit support and p-value
        ks = stats.kstest(residuals,'chi2',args=(self.nfree,))
        self.support,self.p = ks
        self.results['support'] = self.support
        self.results['p-value'] = self.p
        # Calculate loglikelihood
        self.loglikelihood = self.logL(self.params)
        # Calculate AIC, AICc, and BIC
        self.AIC = 2.*(self.nparams-self.loglikelihood); self.results['AIC'] = self.AIC
        self.AICc = self.AIC + (2*(pow(self.nparams,2)+self.nparams))/(self.ndata-self.nparams-1)

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
        self.mcmc_params = np.zeros(self.nparams)
        self.mcmc_std_errors = np.zeros(self.nparams)
        for i,name in enumerate(self.names):
            std_l, median, std_u = percentiles[i]
            self.mcmc_params[i] = median
            self.results['%s_mcmc'%name] = median
            std_error = 0.5 * (std_u - std_l)
            self.mcmc_std_errors[i] = std_error
            self.results['%s_mcmc_error'%name] = 0.5 * (std_u - std_l)

class Gamma_Distribution:
    """
    Gamma probability distribution for AV and PSD.
    """
    def __init__(self,shape,yhat):
        """
        Parameters
        ----------
        shape : np.array
            Shape parameter.
        yhat : np.array
            Experimental values. 
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