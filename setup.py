# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""

from setuptools import setup

setup(
   name='tweezepy',
   version='0.1.1',
   author='Ian Morgan',
   author_email='ilmorgan@ucsb.edu',
   packages=['tweezepy'],
   license='LICENSE.txt',
   description='Saleh lab tweezer code',
   long_description=open('README.md').read(),
   install_requires=[
       "numpy",
       "pandas",
       "scipy",
       "autograd"
       ],
)