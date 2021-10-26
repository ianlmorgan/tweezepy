.. _install:

============
Installation
============
To use ``Tweezepy``, you will need to `install Python <https://realpython.com/installing-python/>`_ on your system if you haven't already. 

Requirements for Installing Packages
====================================

This section describes the steps to follow before installing Python packagaes.

Ensure you can run Python from the command line
-----------------------------------------------

Before you can install Python packages, make sure you have Python and that the expected version is available from your command prompt. You can check this by running:

.. code-block:: bash

    python --version

You should get an output like ``Python 3.8.x``. If you do not have Python,
please install the latest 3.9.x from `python.org`_ or a `specialized distribution for scientific computing <https://www.anaconda.com/products/individual>`_.
Note that, as of 2021-10-25, many packages, including ones required by ``Tweezepy``, are not compatible with Python 3.10. 

.. Note:: If you get an error like this:

    .. code-block:: bash

        python --version
        Traceback (most recent call last):
          File "<stdin>", line 1, in <module>
        NameError: name 'python' is not defined

    You are not running the command on the command line.

Using Python
------------

Once Python is installed, it can be accessed by running:

.. code-block:: bash
    python

After Python opens, it will show some contextual information like this:

.. code-block:: python

    Python 3.8.x (default Jul 2 2020 17:30:36)
    Type "help", "copyright", "credits", or "license" for more information.
    >>>

.. Note:: If you're a newcomer to Python, see the Python for Beginners `getting started tutorial <https://opentechschool.github.io/python-beginners/en/getting_started.html#what-is-python-exactly>`_ for an introduction to using your operating system's terminal and interacting with Python.
    
Installing a Package
--------------------

Installing using pip
####################

The recommended way to install ``Tweezepy`` is using `pip <http://www.pip-installer.org/>`_, the `package installer for Python <https://packaging.python.org/guides/tool-recommendations/>`_, from the `Python Package Index <https://packaging.python.org/glossary/#term-python-package-index-pypi>`_. Usually, pip comes pre-installed with Python.

To check whether you have a working Python installation with pip installed, run the following commands from the command line.

.. code-block:: bash

    pip --version

If pip is installed properly, you should get something like this:

.. code-block:: bash

    pip 21.0.1 from ...

Once you've verified that you have a working Python installation with pip installed, you can run the following command from the command line to download and install ``Tweezepy`` from the Python package index.

.. code-block:: bash

    python -m pip install tweezepy

Installing from source
######################

Alternatively, ``Tweezepy`` can be installed by cloning the source repository from `Github <https://github.com/ianlmorgan/tweezepy>`_:

.. code-block:: bash

    git clone https://github.com/ianlmorgan/tweezepy.git

Once you've downloaded the source, you can navigate into the root source directory and run:

.. code-block:: bash

    python -m pip install .

Tests
*****

If you installed from source, you should run the unit tests to make sure everything worked properly. From the root of the source directory, run:

.. code-block:: bash

    python -m pip install pytest -v tests

This might take a few minutes but you shouldn't get any errors if all went as planned.