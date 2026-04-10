#!/usr/bin/env python

import glob

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup


version = open("VERSION").readline().strip()
sources = ["src/pyteomics/_biolccc.cpp"] + glob.glob("src/core/*.cpp")

ext_modules = [
    Pybind11Extension(
        "pyteomics._biolccc",
        sources=sources,
        include_dirs=["./include"],
        define_macros=[("VERSION", f'\"{version}\"')],
        cxx_std=17,
    )
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
