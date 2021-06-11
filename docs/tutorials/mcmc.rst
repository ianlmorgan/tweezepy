.. module:: tweezepy

**Note:** This tutorial was generated from an IPython notebook that can be
downloaded `here <../../_static/notebooks/mcmc.ipynb>`_.

.. _mcmc:


Robust Uncertainties
--------------------

In most cases, the Hessian method is sufficient for calculating the
parameter uncertainties. However, in some cases (e.g., small sample
sizes and parasitic noise), the Gaussian approximation for the parameter
erros is not valid. In these cases, it may be necessary to run a Monte
Carlo sampler to estimate errors. This method is slower, but it has the
advantage of not assuming

.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from tweezepy import downsampled_trace, PSD, AV

.. code:: python

    fc = 10 # corner frequency
    gamma = 1e-5 # dissipation due to viscous drag, in pN s/nm
                 # 1e-5 is a typical value for an MT experiment
    kappa = gamma*2*np.pi*fc # kappa in pN/nm
    fsample = 400 # sampling frequency in Hz
    N  = 102400 # number of points in trajectory
    seed = 0 # random seed for reproducibility
    xtrace_ds = downsampled_trace(gamma,kappa,fsample,N, seed = seed)
    plt.plot(xtrace_ds,label = 'Downsampled')
    plt.xlabel('t (s)')
    plt.ylabel('x (nm)')
    plt.show()



.. image:: mcmc_files%5Cmcmc_3_0.svg


.. code:: python

    psd = PSD(xtrace_ds,fsample,bins = 15)
    psd.mlefit()
    psd.plot()
    plt.show()



.. image:: mcmc_files%5Cmcmc_4_0.svg


Now, that we have done an MLE fit, we can use the emcee package to
perform Monte Carlo sampling of the parameter space.

.. code:: python

    psd.mcmc()
    


.. parsed-literal::

    100%|██████████| 2000/2000 [01:12<00:00, 27.45it/s]
    

.. code:: python

    psd.calc_mc_errors()




.. parsed-literal::

    array([[9.95272903e-06, 9.98671257e-06, 1.00220528e-05],
           [6.28507478e-04, 6.33540959e-04, 6.38723243e-04]])



We can take a look at the accepted steps. We throw out twice the
autocorrelation time as ‘burn-in’ (black vertical line) and thin by half
the autocorrelation time to make sure the steps are uncorrelated.

.. code:: python

    psd.sample_plot()




.. parsed-literal::

    (<Figure size 720x432 with 2 Axes>,
     [<AxesSubplot:ylabel='g'>, <AxesSubplot:xlabel='step number', ylabel='k'>])




.. image:: mcmc_files%5Cmcmc_9_1.svg


We can use a corner plot to look at the 1 and 2D histograms. This uses
the corner package, which is separate from the emcee package. The MLE
best-fit parameter estimates are indicated by the blue lines. The dotted
black lines in the 1D histograms are the 16th and 84th% percentiles. The
black lines in the 2D histogram are the 1,2, and 3 standard deviations.

.. code:: python

    psd.corner_plot()
    plt.show()



.. image:: mcmc_files%5Cmcmc_11_0.svg


.. code:: python

    # The mcmc results are automatically added to the results dictionary
    psd.results




.. parsed-literal::

    {'chi2': 9541.523752905014,
     'redchi2': 1.4915622561990016,
     'g': 9.987274585927754e-06,
     'g_error': 4.939065221213263e-08,
     'k': 0.0006338066057220219,
     'k_error': 7.270886693644723e-06,
     'support': 1.0,
     'p-value': 0.0,
     'AIC': 13644.000903245327,
     'g_mcmc': 9.986712570927105e-06,
     'g_mcmc_error': 3.466189176357288e-08,
     'k_mcmc': 0.0006335409588759638,
     'k_mcmc_error': 5.107882519335026e-06}



.. code:: python

    # We can also access the mcmc median and standard error values
    print(psd.mcmc_params)
    print(psd.mcmc_std_errors)


.. parsed-literal::

    [9.98671257e-06 6.33540959e-04]
    [3.46618918e-08 5.10788252e-06]
    

