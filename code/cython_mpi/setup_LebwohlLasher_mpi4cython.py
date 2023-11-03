try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

from Cython.Build import cythonize
from distutils.extension import Extension
import numpy


LebwohlLasher_mpi4cython = Extension(
    "LebwohlLasher_mpi4cython",
    sources=["LebwohlLasher_mpi4cython.pyx"],
    extra_compile_args=['-O3'],
    extra_link_args=['-O3']
)

setup(name="LebwohlLasher_mpi4cython",
    ext_modules = cythonize(LebwohlLasher_mpi4cython),
    include_dirs=[numpy.get_include()],  # Include NumPy header files
    packages=find_packages(),
)