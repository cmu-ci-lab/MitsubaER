import os, sys

BUILDDIR       = '#build/release'
DISTDIR        = '#dist'
CXX            = 'g++'
CC             = 'gcc'
CXXFLAGS       = ['-O3', '-std=c++11', '-Wall', '-g', '-pipe', '-march=nocona', '-msse2', '-ftree-vectorize', '-mfpmath=sse', '-funsafe-math-optimizations', '-fno-rounding-math', '-fno-signaling-nans', '-fno-math-errno', '-fomit-frame-pointer', '-DMTS_DEBUG', '-DSINGLE_PRECISION', '-DSPECTRUM_SAMPLES=3', '-DMTS_SSE', '-DMTS_HAS_COHERENT_RT', '-DFLOATDEBUG', '-fopenmp', '-fvisibility=hidden', '-mtls-dialect=gnu2']
LINKFLAGS      = []
SHLINKFLAGS    = ['-rdynamic', '-shared', '-fPIC', '-lstdc++']
BASEINCLUDE    = ['#include', '/usr/include']
BASELIB        = ['dl', 'm', 'pthread', 'gomp']
EIGENINCLUDE   = ['/usr/include/eigen3']
OEXRINCLUDE    = ['/usr/include/OpenEXR']
OEXRLIB        = ['Half', 'IlmImf', 'z']
PNGLIB         = ['png']
JPEGLIB        = ['jpeg']
FASTWINDINGINCLUDE  = ['/home/ubuntu/fast-winding-number-soups/tbb/include', '/home/ubuntu/fast-winding-number-soups/libigl/include', '/home/ubuntu/fast-winding-number-soups/libigl/cmake/../external/eigen', '/home/ubuntu/fast-winding-number-soups']
FASTWINDINGFLAGS    = ['-fPIC', '-rdynamic', '-DNDEBUG']
FASTWINDINGLIBDIR   = ['/home/ubuntu/fast-winding-number-soups/build/CMakeFiles/fastwinding.dir/WindingNumber/', '/home/ubuntu/fast-winding-number-soups/build/tbb/']
FASTWINDINGLIB      = ['pthread', 'dl']
CERESINCLUDE   = ['/home/ubuntu/ceres-bin/config', '/home/ubuntu/ceres-solver-1.14.0/include']
CERESFLAGS     = ['-DNDEBUG', '-DCERES_GFLAGS_NAMESPACE=google']
CERESLIBDIR    = ['/home/ubuntu/ceres-bin', '/usr/lib/x86_64-linux-gnu', '/usr/lib']
CERESLIB      = ['ceres', 'glog', 'gflags', 'pthread', 'spqr', 'cholmod', 'ccolamd', 'camd', 'colamd', 'amd', 'lapack', 'f77blas', 'atlas', 'suitesparseconfig', 'rt', 'cxsparse', 'gomp', 'pthread', 'spqr', 'cholmod', 'ccolamd', 'camd', 'colamd', 'amd', 'lapack', 'f77blas', 'atlas', 'suitesparseconfig', 'rt', 'cxsparse', 'gomp']
XERCESINCLUDE  = []
XERCESLIB      = ['xerces-c']
GLLIB          = ['GL', 'GLU', 'GLEWmx', 'Xxf86vm', 'X11']
GLFLAGS        = ['-DGLEW_MX']
BOOSTLIB       = ['boost_system', 'boost_filesystem', 'boost_thread']
COLLADAINCLUDE = ['/usr/include/collada-dom', '/usr/include/collada-dom/1.4']
COLLADALIB     = ['collada14dom', 'xml2']
FFTWLIB        = ['fftw3_threads', 'fftw3']

# The following runs a helper script to search for installed Python
# packages that have a Boost Python library of matching version.
# A Mitsuba binding library will be compiled for each such pair.
# Alternatively, you could also specify the paths and libraries manually
# using the variables PYTHON27INCLUDE, PYTHON27LIB, PYTHON27LIBDIR etc.

import sys, os
sys.path.append(os.path.abspath('../data/scons'))
from detect_python import detect_python
locals().update(detect_python())

