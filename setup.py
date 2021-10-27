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
VERSION = "1.1.6"
DESCRIPTION = "Single-molecule force spectroscopy calibration"
AUTHOR = "Ian L. Morgan"
AUTHOR_EMAIL = "ilmorgan@ucsb.edu"
URL = "https://github.com/ianlmorgan/tweezepy"
LICENSE = "LGPLv3+"

INSTALL_REQUIRES = [
                    "matplotlib",
                    "numpy>=1.15,<1.20",
                    "scipy",
                    ]
PACKAGE_DATA = {
    # Include example trajectory data in the "data" subdirectory
    "trajectory": ['data/*.csv'],
}
PYTHON_REQUIRES = ">=3.7, !=3.10.*"
SETUP_REQUIRES = INSTALL_REQUIRES + [
    "setuptools>=40.6.0",
]
TESTS_REQUIRE = [
    "unittest",
]
CLASSIFIERS = [
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
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
                 include_package_data=True,
                 package_data = PACKAGE_DATA,
                 #data_files=[('trajectory',['data/trajectory.csv'])]
                 python_requires=PYTHON_REQUIRES,
                 long_description=long_description,
                 long_description_content_type = 'text/markdown',
                 install_requires = INSTALL_REQUIRES,
                 setup_requires = SETUP_REQUIRES,
                 tests_require=TESTS_REQUIRE,
                 classifiers=CLASSIFIERS,
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