#!/usr/bin/env python

"""
setup.py file for pyBioLCCC
"""

from distutils.core import setup, Extension

sources=[
   'aminoacid.cpp',
   'terminus.cpp',
   'chemicalbasis.cpp',
   'gradientpoint.cpp',
   'gradient.cpp',
   'chromoconditions.cpp',
   'BioLCCC.cpp',
   'pyBioLCCC_wrap.cxx',
   ]

pyBioLCCC_module = Extension('_pyBioLCCC',
                           sources=sources,
                           language='c++',
                           )

setup (name = 'pyBioLCCC',
       version = '0.1',
       author      = "Anton Goloborodko",
       description = """pyBioLCCC package""",
       ext_modules = [pyBioLCCC_module],
       py_modules = ['pyBioLCCC'],
       package_dir = {'pyBioLCCC':''},
       )
