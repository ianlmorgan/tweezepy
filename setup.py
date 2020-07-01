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

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setuptools.setup(
   name='tweezepy',
   version='0.1.3',
   author='Ian Morgan',
   author_email='ilmorgan@ucsb.edu',
   packages=setuptools.find_packages(),
   license='LICENSE.txt',
   description='Single-molecule pulling analysis package',
   long_description=long_description,
   long_description_content_type='text/markdown',
   url='https://github.com/pypa/tweezepy',
   install_requires = ['autograd'],
   project_urls={  # Optional
      'Source': 'https://github.com/ianlmorgan/tweezepy/',
    },
)