#!/usr/bin/env python

"""
setup.py file for pyBioLCCC
"""

from distutils.core import setup, Extension
import glob
#from setuptools import setup, Extension

sources = glob.glob("src/core/*.cpp") + ['./src/bindings/pyBioLCCC_wrap.cc']
#sources=[
#    './src/core/biolcccexception.cpp',
#    './src/core/chemicalgroup.cpp',
#    './src/core/chemicalbasis.cpp',
#    './src/core/gradientpoint.cpp',
#    './src/core/gradient.cpp',
#    './src/core/chromoconditions.cpp',
#    './src/core/biolccc.cpp',
#    './src/bindings/pyBioLCCC_wrap.cc',
#    #'./src/bindings/pyBioLCCC.i',
#    ]

pyBioLCCC_ext = Extension(
    '_pyBioLCCC',
    sources=sources,
    include_dirs=['./include'],
    )

version = open('./VERSION').readline().strip()
setup(name = 'pyBioLCCC',
    version = version,
    description      = """Bindings for the libBioLCCC""",
    long_description = ''.join(open('README').readlines()),
    author           = 'Anton Goloborodko',
    author_email     = 'goloborodko.anton@gmail.com',
    url              = 'http://theorchromo.ru',
    ext_modules      = [pyBioLCCC_ext],
    py_modules       = ['pyBioLCCC'],
    classifiers      = ['Intended Audience :: Science/Research',
                      'Topic :: Scientific/Engineering :: Bio-Informatics',
                      'Topic :: Scientific/Engineering :: Chemistry',
                      'Topic :: Scientific/Engineering :: Physics',],
    license          = 'License :: Free for non-commercial use',
    package_dir      = {'': './src/bindings'},
    )
