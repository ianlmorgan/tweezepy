.. _install:

Installation
============
To use ``Tweezepy``, you will need to `install Python <https://realpython.com/installing-python/>`_ on your system if you haven't already. 

Package manager
---------------
``Tweezepy`` can be found on the `Python Package Index <https://packaging.python.org/glossary/#term-python-package-index-pypi>`_ and installed using `pip <http://www.pip-installer.org/>`_, the `package installer for Python <https://packaging.python.org/guides/tool-recommendations/>`_. Usually, pip comes pre-installed with Python.

To check whether you have a working Python installation with pip installed, run the following commands from the command line.

.. code-block:: bash

    python --version
    pip --version

Once you've verified that you have a working Python installation with pip installed, you can run the following command to download and install ``Tweezepy`` from the Python package index.

.. code-block:: bash

    python -m pip install tweezepy

From source
-----------

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