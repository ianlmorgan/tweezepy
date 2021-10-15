# -*- coding: utf-8 -*-
"""
Created on Tue May  5 16:31:29 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
import setuptools
#import json
import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

NAME = "Tweezepy"
VERSION = "1.1.3"
DESCRIPTION = "Single-molecule force spectroscopy calibration"
AUTHOR = "Ian L. Morgan"
AUTHOR_EMAIL = "ilmorgan@ucsb.edu"
URL = "https://github.com/ianlmorgan/tweezepy"
LICENSE = "LGPLv3+"
INSTALL_REQUIRES = [
                    "numpy>=1.15,<1.20",
                    "scipy",
                    "autograd",
                    "numba",
                    "emcee",
                    ]
SETUP_REQUIRES = INSTALL_REQUIRES + [
    "setuptools>=40.6.0",
]
TESTS_REQUIRE = [
    "unittest",
]
#pkginfo_path = os.path.join(this_directory,
#                            'tweezepy',
#                            'tweezepy_info.json')
#pkginfo = json.load(open(pkginfo_path))

# This call to setup() does all the work
setuptools.setup(name=NAME,
                 version=VERSION,
                 description=DESCRIPTION,
                 author=AUTHOR,
                 author_email=AUTHOR_EMAIL,
                 url=URL,
                 license=LICENSE,
                 packages=setuptools.find_packages(),
                 long_description=long_description,
                 long_description_content_type = 'text/markdown',
                 install_requires = INSTALL_REQUIRES,
                 setup_requires = SETUP_REQUIRES,
                 tests_require=TESTS_REQUIRE,
                 ) 

# This call to setup() does all the work
#setuptools.setup(name=pkginfo['name'],
#                 version=pkginfo['version'],
#                 description=pkginfo['description'],
#                 author=pkginfo['author'],
#                 author_email=pkginfo['author_email'],
#                 url=pkginfo['url'],
#                 license=pkginfo['license'],
#                 packages=setuptools.find_packages(),
#                 long_description=long_description,
#                 long_description_content_type = 'text/markdown',
#                 install_requires = ['numba',
#                                     'numpy',
#                                     'scipy',
#                                     'autograd',
#                                     'emcee'],
#                 tests_require=['unittest',
#                                'numpy'],
#                 ) 