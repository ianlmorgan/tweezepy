# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""

import setuptools

setuptools.setup(
   name='tweezepy',
   version='0.1.2',
   author='Ian Morgan',
   author_email='ilmorgan@ucsb.edu',
   packages=setuptools.find_packages(),
   license='LICENSE.txt',
   description='Single-molecule pulling analysis',
   long_description=open('README.md').read()
)