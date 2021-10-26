Welcome to Tweezepy's documentation!
====================================

``Tweezepy`` is a Python package for calibrating forces in single-molecule force spectroscopy experiments and these pages will show you how to use it.

Development of ``Tweezepy`` happens on `Github <https://github.com/ianlmorgan/tweezepy>`_ so you can raise any bugs, documentation improvements, or feature reqeusts `there <https://github.com/ianlmorgan/tweezepy/issues>`_.

How to use this guide
---------------------
To install ``Tweezepy``, you should follow the :ref:`install` guide. After you install it, you can learn about how to use it from the tutorials listed below (start with quickstart and go from there). A description of all the code, including expected inputs and outputs, can be found in the API section.

.. toctree::
    :maxdepth: 1
    :caption: User Guide

    install

.. toctree::
    :maxdepth: 1
    :caption: Tutorials

    pages/quickstart
    pages/mcmc
    pages/simulations

In the API, you can find a list of the classes and modules in ``Tweezepy`` including all their expected inputs and outputs. The main user-facing classes ``AV`` and ``PSD`` are pat of the ``smmcalibration`` module.

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
