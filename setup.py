#!/usr/bin/env python

"""
setup.py file for pyteomics.biolccc
"""

from distutils.core import setup, Extension
import glob
#from setuptools import setup, Extension

sources = glob.glob("src/core/*.cpp") + ['./src/pyteomics/biolccc_wrap.cc']
version = open('./VERSION').readline().strip()

biolccc_ext = Extension(
    'pyteomics._biolccc',
    sources=sources,
    include_dirs=['./include'],
    extra_compile_args=['-DVERSION=\"%s\"' % version,]
    )

setup(name = 'pyteomics.biolccc',
    version = version,
    description      = """Bindings for the libBioLCCC""",
    long_description = ''.join(open('README').readlines()),
    author           = 'Anton Goloborodko',
    author_email     = 'goloborodko.anton@gmail.com',
    url              = 'http://theorchromo.ru',
    ext_modules      = [biolccc_ext],
    packages         = ['pyteomics'],
    classifiers      = ['Intended Audience :: Science/Research',
                      'Topic :: Scientific/Engineering :: Bio-Informatics',
                      'Topic :: Scientific/Engineering :: Chemistry',
                      'Topic :: Scientific/Engineering :: Physics',],
    license          = 'License :: Free for non-commercial use',
    package_dir      = {'pyteomics': './src/pyteomics'},
    )
