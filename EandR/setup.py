#!/usr/bin/env python3

from distutils.core import setup, Extension
from Cython.Build import cythonize
import os
from glob import glob
import numpy as np

srcs = ["threelocus/pyexp.pyx", "threelocus/cexp.cpp"]

setup(
    ext_modules = cythonize(Extension(
        "pyexp",
        srcs,
        extra_compile_args=["-std=c++0x", "-Wno-unused-function", "-Wno-unused-variable"],
        include_dirs=[np.get_include()],
        language="c++"))
)
