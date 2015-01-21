Python and C/C++ Integration
--------------------

Cython and distutils are both included in the AnacondaCE python distribution. The PySmallK framework utilizes both of these packages to allow it's NMF library (high performance custom C/C++ code) to be compiled into shared Python libraries (.so files for OSX&trade; and Linux platforms). 

Distutils provides support for building and installing modules into a Python installation. The PySmallK framework relies on distutils to help build it's NMF library.

This document explains the details of using Cython to compile the .so or .dll for user's custom C/C++ code. A detailed example and unit test is also included examples directory.

Using Distutils and Cython to build nmflib.so
-------------------------
	# Change directory to the smallk/PySmallK python directory
		cd path/to/smallk/PySmallK/cython

	# Build libnmfpy.so/libnmfdistpy.so
	python setup.py build_ext --inplace

	# Execute an example that uses the libnmfpy.so library
	python examples/test_nmflib.py --k 32 --rows 256 --cols 256 --tol 0.005 --algorithm HALS

	# Test an example that runs the distributed version with 4 processes
	mpirun -np 4 python examples/test_nmflib_dist.py --k 32 --rows 256 --cols 256 --tol 0.005 --algorithm HALS --gridrows=2 --gridcols=2 --verbose

Distutils - setup.py
--------------------

Distutils uses a setup.py script to define the various options for building a module. Within the script, calling setup() supplies distutils with the information to assist in building, distributing and installing modules into Python. 

Below is an example of a minimal setup.py script:

	from distutils.core import setup
	from distutils.extension import Extension
	from Cython.Distutils import build_ext
	ext_modules = [Extension("file", ["somelib.py"])]
	setup(
	  name = 'name of the app',
	  cmdclass = {'build_ext': build_ext},
	  ext_modules = ext_modules
	)

Build using:
	python setup.py build_ext --inplace


Extension modules
-----------------

The file setup.py describes an extension module for building the nmflib software. In it's description, it provides the name for the extension module, the source for the interface file (interface/nmf_lib.pyx), include directories (which includes the numpy.get_include() function), libraries, library directories, and the language of the module (C++).

Environment variables for compiler paths can be included by defining the os.environ dictionary, e.g. os.environ['CPP'] = '/usr/local/bin/g++-4.7'

Cython & the interface file
---------------------------
The interface file (interface/nmf_lib.pyx) is a Python language extension which allows for the writing of fast Python extension modules and interfaces with C/C++ libraries. Cython code allows us to take advantage of the speed of the C/C++ NMF library and has support for distutils.

To generate C source for our Python extension, you may execute:

	cython nmf_lib.pyx

This will generate C code that compiles with all major C/C++ compilers, major platforms, and in Python 2.3 to 3.1. Currently Cython language syntax follows Python 2.6 and supports top-level classes and functions, control loops, object operations, arithmetic, etc, along with Python 3 features such as lists, set, dict comprehensions, keyword-only arguments and extended utterable unpacking.

Cython generates very efficient C code. Benchmarks on optimized code have been shown to run code significantly faster than in regular Python code.

Optimizing Cython code in the interface file
--------------------------------------------
There are several ways to optimize Cython code. We will briefly overview a few declarations that provide extreme optimization.

First, declaring Python types can provide a 30% speed increase

Rather than

	def f(d):
	  return ...

We declare our dict

	def f(dict d):
	  return ...

Using a <i>cdef</i> keyword provides a C signature and C call semantics. <i>cdef</i> can be used when C attribute types are used (e.g. wrapping C structs/pointers/etc) and when speed is more important than the convenience of Pythons lack of need for type declaration. Below are some examples where you can use cdef to optimize your Cython code to use C types.

Use variables with C or builtin types
	cdef double x
Use functions with C signatures
	cdef double f(double x):
	  return x*x
Use classes with 'builtin' extension types
	cdef class MyClass:
	  cdef int ivar

The <i>cpdef</i> keyword can be used to wrap Python around a <i>cdef</i> function. C will call the <i>cdef</i> function and Python will call the wrapper. 

Let's assume we have some Python code: 

	from math import sin
	def f(x):
      return sin(x**2)

    def integrate_f(a, b, N):
      dx = (b-a)/N
      s = 0
      for i in range(N):
        s += f(a+i*dx)
      return s * dx

This Python code can be optimized to Cython code:

	cdef extern from "math.h":
	  double sin(double x)
    
	cdef double f(double x):
	  return sin(x**2)

    cpdef double integrate_f(double a, double b, int N):
      cdef double dx, s
      cdef int i

      dx = (b-a)/N
      s = 0
      for i in range(N):
        s += f(a+i*dx)
      return s * dx

Disclaimer
----------

This software is a work in progress.  It will be updated throughout the course of the 
XDATA program with additional algorithms and examples.  Any given instance of this software may not
function as expected. Every effort will be made to ensure that the test cases produce the expected results as development proceeds.

References
------------

1. Markdown created

Contact Info
------------

For comments, questions, bug reports, suggestions, etc., contact:
    	Daniel Lee
    	Research Scientist
    	Georgia Tech - School of Computational Science & Engineering
    	dlee@cc.gatech.edu


