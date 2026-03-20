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
        extra_compile_args=[f'-DVERSION=\"{version}\"', "-std=c++14"],
    )
]

setup(ext_modules=ext_modules, cmdclass={"build_ext": build_ext})
