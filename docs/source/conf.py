# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
from pkg_resources import DistributionNotFound,get_distribution

try:
    __version__ = get_distribution("tweezepy").version
except DistributionNotFound:
    __version__ = "unknown version"

# -- Project information -----------------------------------------------------

project = 'Tweezepy'
copyright = '2021, Ian L. Morgan'
author = 'Ian L. Morgan'

# The full version, including alpha/beta/rc tags
version = __version__
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx.ext.autodoc",
              "sphinx.ext.intersphinx",
              "sphinx.ext.napoleon",
              "sphinx.ext.mathjax",
              "nbsphinx",
              ]

#source_suffix = {
#                ".rst":"restructuredtext",
#                ".ipynb":"nbsphinx",
#}
source_suffix = ".rst"
master_doc = "index"

# Add any paths that contain templates here, relative to this directory.
#templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    "_build", 
    "pages/**.ipynb_checkpoints",
    ]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'
html_copy_source = True
html_show_sourcelink = True
html_sourcelink_suffix = ""
html_title = 'Tweezepy'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
html_theme_options = {
    "path_to_docs": "docs",
    "repository_url": "https://github.com/ianlmorgan/tweezepy",
    #"repository_branch": "ian-local", # For testing
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "colab_url": "https://colab.research.google.com/",
        #"notebook_interface": "classic",
    },
    "use_edit_page_button": True,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_download_button": True,
}
#html_baseurl = "https://tweezepy.readthedocs.io/en/latest/"
jupyter_execute_notebooks = "off"