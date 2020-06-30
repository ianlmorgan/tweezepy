# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
import pathlib
import setuptools
# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setuptools.setup(
   name='tweezepy',
   version='0.1.2',
   author='Ian Morgan',
   author_email='ilmorgan@ucsb.edu',
   packages=setuptools.find_packages(),
   license='LICENSE.txt',
   description='Single-molecule pulling analysis',
   long_description=open('README.md').read(),
   install_requires = ['autograd']
)