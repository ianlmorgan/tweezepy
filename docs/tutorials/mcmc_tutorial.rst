.. module:: tweezepy

**Note:** This tutorial was generated from an IPython notebook that can be
downloaded `here <../../_static/notebooks/mcmc_tutorial.ipynb>`_.

.. _mcmc_tutorial:


Monte carlo sampler for uncertainty
-----------------------------------

For most uses, the Hessian method should be sufficient for calculating
parameter uncertainties. However, it is possible that instrumental noise
could cause non-linear effects that disrupt the gaussian approximation
in the Hessian method. For these cases, tweezepy comes built in with a
Markov chain Monte Carlo (MCMC) sampler that can sample the parameter
space to estimate non-gaussian noise.

Here, we demonstrate how to use the MCMC feature.

To start, we’ll simulate some data.

.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from tweezepy.smmcalibration import PSD
    from tweezepy.simulations import downsampled_trace

.. code:: python

    fc = 10 # corner frequency
    alpha = 1e-5 # dissipation due to viscous drag, in pN s/nm
                 # 1e-5 is a typical value for an MT experiment
    kappa = alpha*2*np.pi*fc # kappas in pN/nm
    fsample = 400 # sampling frequency in Hz
    N  = 10240 # number of points in trace
    seed = 0 # gives the same random numbers each time
    time = np.arange(N)/fsample # generate time for plotting purposes
    xtrace = downsampled_trace(alpha,kappa,fsample,N,seed = seed) # simulated the position data
    plt.plot(time, xtrace) # plot the data
    plt.xlabel('Time (s)'); plt.ylabel('Position (nm)')
    plt.show()



.. image:: mcmc_tutorial_files%5Cmcmc_tutorial_4_0.svg


Now, we’ll calculate the power spectral densities and fit them using the
maximum likelihood estimation.

.. code:: python

    psd = PSD(xtrace,fsample)
    psd.mlefit()
    psd.plot();



.. image:: mcmc_tutorial_files%5Cmcmc_tutorial_6_0.svg


Next, we’ll run the MCMC sampler. The default values in tweezepy are 32
walkers and 1500 steps of the MCMC. On my computer, this usually takes
15 seconds. This is usually sufficient to get good estimates, but the
user can test out different values.

.. code:: python

    psd.mcmc(steps=1500) 

Let’s take a look at what the sampler has done. This usually happens
“under the hood,” but it’s useful to get a sense of how these things
work. The figure below shows the position of each walker as a function
of the number of steps in the chain.

.. code:: python

    sampler = psd.sampler
    samples = sampler.get_chain()
    fig, axes = plt.subplots(2, figsize=(8, 7), sharex=True)
    labels = [r"$\alpha$", r"$\kappa$"]
    for i,label in enumerate(labels):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(label)
    axes[-1].set_xlabel("step number");



.. image:: mcmc_tutorial_files%5Cmcmc_tutorial_10_0.svg


The walkers start in small distributions around the maximum likelihood
values and then quickly wander and start exploring the full
distribution. After fewer than 50 steps, the samples are pretty well
“burnt-in.” To check this statement quantitatively, we can look at the
integrated autocorrelation time.

.. code:: python

    tau = sampler.get_autocorr_time()
    print(tau)


.. parsed-literal::

    [25.07859788 22.33785708]
    

This suggests that only about 25 steps are required for the chain to
“forget” where it started. So, we’ll throw away 4x this amount as
“burn-in” and thin the samples by about half the autocorrelation time
(10 steps).

.. code:: python

    flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)
    print(flat_samples.shape)


.. parsed-literal::

    (4480, 2)
    

Results
-------

With out list of samples , we’ll take a look at a corner plot. To do so,
you’ll need the corner.py module.

.. code:: python

    import corner
    
    fig = corner.corner(
        flat_samples, labels=labels, truths=[alpha, kappa],
    );



.. image:: mcmc_tutorial_files%5Cmcmc_tutorial_16_0.svg


This plot shows us all the one and two dimensional projections of the
probability distributions of the parameters. This allows us to see if
their are any covariances between parameters. From this, we see that the
distributions are fairly symmetrical, suggesting that the gaussian
approximation would be appropriate. It also shows us that the noise in
the simulation throws off our estimation of alpha slightly, but does a
good job of capturing kappa.

Typically, parameter errors are given as 1\ :math:`\sigma`, so, in
analogy, we calculate the 16th, 50th, and 84th percentiles.

.. code:: python

    from IPython.display import display, Math
    labels = [r"\alpha",r"\kappa"]
    for i in range(2):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        txt = "{{{3}}} = {0:.3g}_{{-{1:.3g}}}^{{{2:.3g}}}"
        txt = txt.format(mcmc[1], q[0], q[1], labels[i])
        display(Math(txt))



.. math::

    \displaystyle {\alpha} = 9.8e-06_{-1.08e-07}^{1.07e-07}



.. math::

    \displaystyle {\kappa} = 0.000618_{-1.7e-05}^{1.6e-05}


