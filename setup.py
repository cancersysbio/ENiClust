#!/usr/bin/env python

"""
This is a setup script for ENiClust.

@author: Aziz Khan
@email: azizk@stanford.edu

"""
import os
from distutils.core import setup
from setuptools import find_packages
from eniclust import __version__ as VERSION

CLASSIFIERS = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
]

install_requires = [
    'matplotlib==3.2.2',
    'pandas==1.1.0',
    'numpy==1.21.4',
    'seaborn==0.11.1',
    'scikit-learn==0.21.3',
    'imbalanced-learn==0.5'
]

#def readme():
#    with open('README.rst') as f:
#        return f.read()

def readme(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="ENiClust",
    description="Ensemble Integrative Clustering",
    version=VERSION,
    author="Curtis Lab",
    license='MIT',
    platforms='linux/unix',
    author_email="",
    url="https://github.com/cancersysbio/ENiClust",
    long_description=readme("README.md"),
    package_dir={'eniclust': 'eniclust'},

    packages=['eniclust',
        #'eniclust.ClassifyIC',
        #'eniclust.PreProc'
        ],

    scripts=['eniclust/ENiClust.py',
              'eniclust/ClassifyIC.py',
              'eniclust/PreProc.py',
              'eniclust/prepare_ENiClust.R',
                   ],
    package_data={'eniclust': ['data/*','data/*/*', 'data/*/*/*',]},
    include_package_data=True,
    install_requires = install_requires,
    classifiers=CLASSIFIERS,
)
