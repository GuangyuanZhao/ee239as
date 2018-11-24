from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import os

#run setup.py build_ext --inplace

sources = ['csunpos.pyx']
cd = os.getcwd()

setup(
    ext_modules = cythonize([Extension('sunpos.csunpos',sources,include_dirs=[cd,np.get_include()],define_macros=[('_CRT_SECURE_NO_WARNINGS',None)])])
)