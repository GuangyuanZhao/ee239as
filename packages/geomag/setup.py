from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import os

sources = ['emm.pyx','GeomagnetismLibrary.c','Mesh_SubLibrary.c']
cd = os.getcwd()

setup(
    ext_modules = cythonize([Extension('geomag.emm',sources,include_dirs=[cd],define_macros=[('_CRT_SECURE_NO_WARNINGS',None)])])
)