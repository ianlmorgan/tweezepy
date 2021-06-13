# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
import setuptools
import json
import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

pkginfo_path = os.path.join(this_directory,
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
                 long_description=long_description,
                 long_description_content_type = 'text/markdown',
                 install_requires = ['numpy',
                                     'scipy',
                                     'autograd',
                                     'emcee'],
                 tests_require=['unittest',
                                'numpy'],
                 ) 