# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
import pathlib
import setuptools
# The directory containing this file
here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

# This call to setup() does all the work
setuptools.setup(
   name='tweezepy',
   version='1.0.0',
   author='Ian Morgan',
   author_email='ilmorgan@ucsb.edu',
   packages=setuptools.find_packages(),
   license='LICENSE.txt',
   description='Single-molecule video-tracking calibration package',
   long_description=long_description,
   long_description_content_type='text/markdown',
   url='https://github.com/pypa/tweezepy',
   install_requires = ['autograd','emcee'],
   project_urls={  # Optional
      'Source': 'https://github.com/ianlmorgan/tweezepy/',
    },
)