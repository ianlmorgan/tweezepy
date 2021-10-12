.. _install:

Installation
============

Dependencies
------------

``Tweezepy`` depends on ``numpy``,  ``scipy``, and ``autograd``. Optionally, ``matplotlib`` is used for plotting and ``emcee`` is used for monte carlo sampling. You can install these using your favorite Python package manager.

Using pip
---------

The recommended way to install the latest stable verision of ``tweezepy`` is
with `pip <http://www.pip-installer.org/>`_:

.. code-block:: bash

    python -m pip install tweezepy

From source
-----------

Alternatively, ``tweezepy`` can be installed by cloning the source repository from `Github <https://github.com/ianlmorgan/tweezepy>`_:

.. code-block:: bash

    git clone https://github.com/ianlmorgan/tweezepy.git

Once you've downloaded the source, you can navigate into the root source directory and run:

.. code-block:: bash

    python -m pip install .

Tests
-----

If you installed from source, you should run the unit tests to make sure everything worked properly. From the root of the source directory, run:

.. code-block:: bash

    python -m pip install pytest -v tests

This might take a few minutes but you shouldn't get any errors if all went as planned.