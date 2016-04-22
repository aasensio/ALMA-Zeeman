from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'make clean ; make'
print(fortran_mod_comp)
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         "pydelopar",
                         # source file:
                         ['formal.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_link_args=['formal.o', '-lgfortran'])]

setup(name = 'pydelopar',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

system('cp pydelopar.*.so ../')
