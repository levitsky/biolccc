#!/usr/bin/env python

"""
setup.py file for pyteomics.biolccc
"""

import glob
from setuptools import setup, Extension
#from distutils.core import setup, Extension

sources = glob.glob("src/core/*.cpp") + ['pyteomics/biolccc_wrap.cxx']
version = open('./VERSION').readline().strip()

biolccc_ext = Extension(
    'pyteomics._biolccc',
    sources=sources,
    include_dirs=['./include'],
    extra_compile_args=['-DVERSION=\"%s\"' % version, '-std=c++14']
    )

setup(
    name = 'pyteomics.biolccc',
    py_modules = ['pyteomics.biolccc'],
    namespace_packages = ['pyteomics'],
    version = version,
    description      = """Bindings for the libBioLCCC""",
    long_description = ''.join(open('README.rst').readlines()),
    author           = 'Anton Goloborodko',
    author_email     = 'goloborodko.anton@gmail.com',
    url              = 'http://theorchromo.ru',
    ext_modules      = [biolccc_ext],
    classifiers      = ['Intended Audience :: Science/Research',
                      'Topic :: Scientific/Engineering :: Bio-Informatics',
                      'Topic :: Scientific/Engineering :: Chemistry',
                      'Topic :: Scientific/Engineering :: Physics',],
    license          = 'License :: Free for non-commercial use',
    )
