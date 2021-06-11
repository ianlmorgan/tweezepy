# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
from io import open
import setuptools
import json
import os

pkginfo_path = os.path.join(os.path.dirname(__file__),
                            'tweezepy',
                            'tweezepy_info.json')
pkginfo = json.load(open(pkginfo_path))

# This call to setup() does all the work
setuptools.setup(name=pkginfo['name'],
                 version=pkginfo['version'],
                 description=pkginfo['description'],
                 author=pkginfo['author'],
                 author_email=pkginfo['author_email'],
                 url=pkginfo['url'],
                 license=pkginfo['license'],
                 packages=setuptools.find_packages(),
                 long_description=open('README.md', 'r', encoding='utf8').read(),
                 install_requires = ['numpy',
                                     'scipy',
                                     'autograd',
                                     'emcee'],
                 tests_require=['pytest','numpy'],
                 ) 