try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

from Cython.Build import cythonize
from distutils.extension import Extension
import numpy


LebwohlLasher = Extension(
    "LebwohlLasher_cython",
    sources=["LebwohlLasher_cython.pyx"],
    extra_compile_args=['-O3'],
    extra_link_args=['-O3']
)

setup(name="LebwohlLasher",
    ext_modules = cythonize(LebwohlLasher),
    include_dirs=[numpy.get_include()],  # Include NumPy header files
    packages=find_packages(),
)