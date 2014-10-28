Documentation for SmallK (including build and installation instructions) 
can be found in doc/smallk_readme.pdf.

Setup instructions for first time use:

Create a folder for your repos, such as ~/repos/bitbucket.

Change directories into your repo folder:

       cd ~/repos/bitbucket

Instructions for cloning the 'xdata_data' and 'xdata3' repos can be found
on your dashboard at bitbucket.org.  You should have access to both projects,
and if you click on each, you can find the clone commands.

    Clone the xdata_data repo.
    Clone the xdata3 repo.

You should see the following directory structure:

    ~/repos/bitbucket/xdata_data
    ~/repos/bitbucket/xdata3

Follow the smallk_readme.pdf documentation and install Elemental 
and all dependencies.

In order to build pysmallk (the python library that provides simple access 
to SmallK), a number of scientific python libraries are required. These
are easiest acquired via Homebrew. The default installation directory
for brew is /usr/local/lib/. 

IMPORTANT: If you already have an installation of numpy in /usr/local/,
you should uninstall it (i.e. pip uninstall numpy). This will ensure
a clean environment to work in. You will also need a gfortran compiler
to build the numpy libraries. This can be acquired (as brew will 
recommend) by running:

	brew install gcc

This can take a long time to run and may be redundant if you already have
a recent GCC available on your system. You can also directly set the 
environment variable for the fortran compiler:

	export FC=/path/to/your/gcc/bin/gfortran-*

To install numpy, run:

	brew install python
	brew install numpy
	brew install scipy

To check your installation, run:

	brew test numpy

IMPORTANT: Check to see that your numpy installation has correctly linked to
the needed BLAS libraries.

Ensure that you are running the correct python:

	which python

This should print out /usr/local/bin/python.

Open a python terminal and run the following:

	import numpy as np
	np.__config__.show()

For Mac, you should see something similar to the following:

	lapack_opt_info:
	    extra_link_args = ['-Wl,-framework', '-Wl,Accelerate']
	    extra_compile_args = ['-msse3']
	    define_macros = [('NO_ATLAS_INFO', 3)]
	blas_opt_info:
	    extra_link_args = ['-Wl,-framework', '-Wl,Accelerate']
	    extra_compile_args = ['-msse3', '-I/System/Library/Frameworks/vecLib.framework/Headers']
	    define_macros = [('NO_ATLAS_INFO', 3)]

If you are using OpenBLAS, you should see that indicated as well. 

To install Cython:
	
	pip install cython

The Makefile assumes an installation path of /usr/local/lib/python2.7/site-packages
for the compiled library file. If you are not using brew to install your packages,
you will need to tell the Makefile where the appropriate site-packages directory is.
This can be done by setting the SITE_PACKAGES_DIR command line variable when
running ‘make’.

After installing the above dependencies, checkout the 'develop' branch of xdata3:

      cd xdata3
      git checkout -b develop origin/develop

Verify that you are on the 'develop' branch:

       git branch

You should see two branches (master and develop); the develop branch
should have an asterisk by it.

Build the xdata3 code as follows:

      make clean
      make all

Test the build as follows:

     make check

The test script should report that all tests passed.

If you place the 'xdata_data' repo in a different location, you will need to 
explicitly tell the build system where to find it.  In this case you will need
to run this command, instead of 'make check' (don't type the brackets):

    make check DATA_DIR=<path_to_xdata_data>

If you build a tarball, you will need to run either of these commands:

    make distcheck
    make distcheck DATA_DIR=<path_to_xdata_data>

If you test individual python files or write your own, you will need to ensure that
the pysmallk.so shared library file is accessible. This is most simply accomplished
by running ‘make install’, which will install the shared library in the appropriate
site-packages directory for the python binary being used. Alternatively, you could
set the PYTHONPATH variable to include the directory in which pysmallk.so is found 
(default is xdata3/pysmallk/) or copy pysmallk.so into the directory in which 
the test files are located. 

