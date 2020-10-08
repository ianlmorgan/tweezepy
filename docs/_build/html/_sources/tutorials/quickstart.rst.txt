.. module:: tweezepy

**Note:** This tutorial was generated from an IPython notebook that can be
downloaded `here <../../_static/notebooks/quickstart.ipynb>`_.

.. _quickstart:


Quickstart
----------

To get you started, here’s an annotated example of using tweezpy with
simulated data.

To start, we’ll import some necessary modules and functions:

.. code:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from tweezepy.smmcalibration import PSD,AV
    from tweezepy.simulations import downsampled_trace

Let’s simulate some downsampled data to emulate the effects of lowpass
filtering from a camera.

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



.. image:: quickstart_files%5Cquickstart_4_0.svg


Now that we have some data, let’s calculate the power spectral densities
and allan variances.

.. code:: python

    psd = PSD(xtrace, fsample) # convert to power spectral densities 
    av = AV(xtrace,fsample)
    psd.plot();av.plot(); # plot power spectral densities and allan variances



.. image:: quickstart_files%5Cquickstart_6_0.svg



.. image:: quickstart_files%5Cquickstart_6_1.svg


Now, we can fit an analytical function to the data and extract the
best-fit parameters.

.. code:: python

    psd.mlefit(); av.mlefit()
    psd.plot(); av.plot();



.. image:: quickstart_files%5Cquickstart_8_0.svg



.. image:: quickstart_files%5Cquickstart_8_1.svg


We see that after performing a fit, when we call psd.plot() or
av.plot(), it automatically overlays the fit on the data.

Once we have the fits, it would be nice to know what best-fit
parameters, parameter uncertainties (using the Hessian method), and
goodness-of-fit.

.. code:: python

    print(psd.params,psd.errors,psd.redchi2)
    print(av.params,av.errors,av.redchi2)


.. parsed-literal::

    [9.79175568e-06 6.18752762e-04] [1.57251007e-07 2.32127388e-05] 1.1976827885649373
    [9.63570984e-06 6.26491847e-04] [1.79912941e-07 1.93802292e-05] 0.5814388323124055
    

Or to display them in a pretty way.

.. code:: python

    from IPython.display import display, Math
    labels = [r"\alpha", r"\kappa"] # labels
    txt = "{0} = {1:.3g}\pm{2:.3g}" # shared text format
    for l,p,e in zip(labels,psd.params,psd.errors):
        display(Math(txt.format(l,p,e))) # display parameters and errors
    display(Math(r"\chi^2 = {:.2f}".format(psd.redchi2))) # display reduced chi-squared



.. math::

    \displaystyle \alpha = 9.79e-06\pm1.57e-07



.. math::

    \displaystyle \kappa = 0.000619\pm2.32e-05



.. math::

    \displaystyle \chi^2 = 1.20


.. code:: python

    from IPython.display import display, Math
    labels = [r"\alpha", r"\kappa"] # labels
    txt = "{0} = {1:.3g}\pm{2:.3g}" # shared text format
    for l,p,e in zip(labels,av.params,av.errors):
        display(Math(txt.format(l,p,e))) # display parameters and errors
    display(Math(r"\chi^2 = {:.2f}".format(av.redchi2))) # display reduced chi-squared



.. math::

    \displaystyle \alpha = 9.64e-06\pm1.8e-07



.. math::

    \displaystyle \kappa = 0.000626\pm1.94e-05



.. math::

    \displaystyle \chi^2 = 0.58


