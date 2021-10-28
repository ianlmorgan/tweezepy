Welcome to Tweezepy's documentation!
====================================

``Tweezepy`` is a Python package for calibrating forces in single-molecule force spectroscopy experiments and these pages will show you how to use it.

Development of ``Tweezepy`` happens on `Github <https://github.com/ianlmorgan/tweezepy>`_ so you can raise any bugs, documentation improvements, or feature requests `there <https://github.com/ianlmorgan/tweezepy/issues>`_.

Installing Tweezepy
-------------------
To install ``Tweezepy``, follow the :ref:`install` guide. 

.. toctree::
    :maxdepth: 1
    :caption: User Guide

    install

Using Tweezepy
--------------
After you install ``Tweezepy``, you can learn how to use it with the tutorials listed below 
(start with quickstart and go from there). 

.. toctree::
    :maxdepth: 1
    :caption: Tutorials

    pages/quickstart
    pages/mcmc
    pages/simulations

API
---
A description of all the modules, classes, and functions in ``Tweezepy``, 
including expected inputs and outputs, can be found in the API section. 
The main user-facing classes ``AV`` and ``PSD`` are pat of the ``smmcalibration`` module.

.. toctree::
    :maxdepth: 1
    :caption: API

    api/smmcalibration
    api/MLE
    api/simulations
    api/allanvar

License & Attribution
---------------------
Built by `Ian L. Morgan <https://github.com/ianlmorgan>`_ and licensed under the GNU General Public License (see ``LICENSE``).

If you make use of ``Tweezepy`` in your work, please cite our code on `Zenodo <https://doi.org/10.5281/zenodo.4948230>`_.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4948230.svg
   :target: https://doi.org/10.5281/zenodo.4948230
