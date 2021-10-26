.. _install:

============
Installation
============

This page contains installation information and links about installing ``Tweezepy``.

Prerequisite: Python
--------------------

**Tweezepy** requires a working Python installation (Python 3.3 or greater).
If you do not have a working installation, you will need to download and install one.

For new users, we recommend using the `Anaconda <https://www.anaconda.com/download>`_
distribution to install Python and Jupyter. Anaconda conveniently installs Python, Jupyter, and other commonly used packages for scientific computing and data science.

To install Python via Anaconda, use the following installation steps:

1. Download `Anaconda <https://www.anaconda.com/download>`_. We recommend 
   downloading Anaconda's latest Python 3 version (currently Python 3.8).
2. Install the version of Anaconda which you downloaded, following the
   instructions on the download page.
3. Congratulations, you have installed Python. 

Alternative Python installation: 
********************************
If you do not want to install Python via the Anaconda distribution,
there are a number of other ways to install Python. For example,
you can install it via the `official distribution <https://www.python.org/downloads>`_.
.. important::

    Currently, Python 3.10.x is incompatible with many standard scientific computing
    packages, including those used by **Tweezepy**. You should install the latest
    version of 3.9.x. until further notice.

If you install Python separately, you may need to install some additional packages.
These additional packages are included in the requirements.txt file.

Ensure you can run Python
-------------------------
Once you install Python, you can access the Python intepreter by going to the 
command prompt and typing

   .. code-block:: bash

       python

To exit the Python interpreter, type

    .. code-block:: python

        quit()

You can check the installed Python version by typing

    .. code-block:: bash
        
        python --version

If everything is installed correctly, you should get an output like this ``Python 3.8.x``.

Installing Tweezepy
-------------------

The recommended way to install **Tweezepy** is using `pip <http://www.pip-installer.org/>`_, 
the `package installer for Python <https://packaging.python.org/guides/tool-recommendations/>`_, 
from the `Python Package Index <https://packaging.python.org/glossary/#term-python-package-index-pypi>`_. 
Usually, pip comes pre-installed with Python.

To check whether you have pip installed, run the following commands from the command prompt.

.. code-block:: bash

    pip --version

If pip is installed properly, you should get something like this:

.. code-block:: bash

    pip 21.0.1 from ...

Once you've verified that you have a working Python installation with pip installed, 
you can run the following command from the command line to download and install **Tweezepy**:

.. code-block:: bash

    python -m pip install tweezepy

Congratulations, you have now installed **Tweezepy**.

Installing from source
----------------------
For advanced Python users who want the latest development version of **Tweezepy**, it can also be installed from source.
To install from source, clone the source repository from `Github <https://github.com/ianlmorgan/tweezepy>`_:

.. code-block:: bash

    git clone https://github.com/ianlmorgan/tweezepy.git

Once you've downloaded the source, you can navigate into the root source directory and run:

.. code-block:: bash

    python -m pip install .

Running tests
*************

If you installed from source, you should run the unit tests to make sure everything worked properly. 
From the root of the source directory, run:

.. code-block:: bash

    python -m pip install -U pytest
    python -m pytest -v tests

This will take a few seconds. You may get a few deprecation warnings, but you shouldn't get any errors if all went as planned.


Optional: Jupyter notebooks
---------------------------

You can run Python code directly in the Python interpreter or as a script in an integrated development editor (IDE), 
such as Spyder, Visual Studio Code, or Sublime text. 

Alternatively, it is often convenient to use Jupyter Notebooks,which is similar to notebook format used in Mathematica.

Jupyter notebooks
-----------------

If you installed Python via the Anaconda distribution, you should already have installed Jupyter.
If you installed Python in a differnt way, you may need to install Jupyter separately.
To install Jupyter via pip, type the following into the command prompt:

    .. code-block:: bash

        python -m pip install jupyter

Once you've installed Jupyter, you can launch a Jupyter notebook via the command prompt:

    .. code-block:: bash

        jupyter notebook

This will open a Jupyter notebook in your browser. For more information on working with
Jupyter notebooks go `here <https://jupyter.readthedocs.io/en/latest/running.html>`_. 