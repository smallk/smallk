############
Installation
############

.. toctree::
   :maxdepth: 8

**********
Quickstart
**********

Vagrant Virtual Machine
=======================

Installing SmallK into a virtual machine (OSX, Linux, Windows) is intended for those who are not doing development and/or do not have a reason to do the full installation on Linux or OSX outlined in the sections to follow.

The complete stack of software dependencies for SmallK as well as SmallK itself can be rapidly set up and configured through use of Vagrant and VirtualBox and the files included in the repository. The Vagrant install has been tested on Linux Ubuntu 16.04, Mac OSX Sierra 10.12.6, and Windows 10. 

Note that the `smallk/vagrant/bootstrap.sh` file can be modified to perform various tasks when provisioning the vagrant session. Consider customizing `bootstrap.sh` to set up a custom install of libsmallk as required.

To deploy the SmallK VM:

**1.** Install `Vagrant <http://www.vagrantup.com/downloads.html>`_ and `VirtualBox <https://www.virtualbox.org/wiki/Downloads>`_.

.. tip::
   Note: For Windows, ensure that you have a VirtualBox version >= 4.3.12. After installing Vagrant, you may need to log out and log back in to ensure that you can run vagrant commands in the command prompt.

**Optional:** `git clone` the `smallk_data <https://github.com/smallk/smallk_data>`_ repository so that it is parallel with the smallk repository. This is an alternate way to test the installation and begin to work with SmallK. This directory can be synced with a directory of the same name in the VM by adding or uncommenting the following line in `smallk/vagrant/Vagrantfile`:

		config.vm.synced_folder "../../smallk_data", "/home/vagrant/smallk_data"	

**2.** From within the `smallk/vagrant` directory, run:
		
		vagrant up
		
This can take as long as an hour to build the VM, which will be based on a minimal Ubuntu 16.04 installation. The `smallk/vagrant/Vagrantfile` can be customized in many ways to change the specifications for the VM that is built. See more information `here <http://docs.vagrantup.com/v2/>`_. The default configuration provides the VM with 4 GB of memory and 3 CPUs. Increasing these allocations will improve the performance of the application. This can be done by modifying these lines in the Vagrantfile:

		vb.memory = 4096
		vb.cpus = 3

After ‘vagrant up’ has completed, the SmallK and pysmallk libraries will have been built and tested. Additionally, the smallk_data directory, if cloned as in the optional step above, will have been synced into the VM. For more details regarding what is being built and executed while provisioning the VM, please inspect `smallk/vagrant/bootstrap.sh`.

[--back to top--](#top)

**3.** Once the VM has been built, run:

		vagrant ssh

.. tip::
   Note: For Windows, you will need an ssh client in order to run the above command. This can be obtained via `CygWin <https://www.cygwin.com/>`_ `MinGW <http://sourceforge.net/projects/mingw/files/>`_, or `Git <http://git-scm.com/downloads>`_. If you would like to use PuTTY to connect to your virtual machine, follow `these <https://github.com/Varying-Vagrant-Vagrants/VVV/wiki/Connect-to-Your-Vagrant-Virtual-Machine-with-PuTTY>`_ instructions.

In case you need it, the username/password for the VM created will be vagrant/vagrant.

This will drop you into the command line of the VM that was just created, in a working directory at `/home/vagrant`. From there, you can navigate to `/home/vagrant/libsmallk-<version>`, (e.g., libsmallk-1.6.2), and run::

		make check PYSMALLK=1 ELEMVER=0.85 DATA_DIR=../smallk_data		
		
to verify your installation was successful. 

**4.** To test the installation at the command line, run::

		nmf

This will produce the help output for the nmf library function::

	Usage: nmf
	        --matrixfile <filename>  Filename of the matrix to be factored.
	                                 Either CSV format for dense or MatrixMarket format for sparse.
	        --k <integer value>      The common dimension for factors W and H.
	        [--algorithm  BPP]       NMF algorithms:
	                                     MU:    multiplicative updating
	                                     HALS:  hierarchical alternating least squares
	                                     RANK2: rank2 with optimal active set selection
	                                     BPP:   block principal pivoting
	        [--stopping  PG_RATIO]   Stopping criterion:
	                                     PG_RATIO: Ratio of projected gradients
	                                     DELTA:    Change in relative F-norm of W
	        [--tol  0.005]           Tolerance for the selected stopping criterion.
	        [--tolcount  1]          Tolerance count; declare convergence after this many
	                                 iterations with metric < tolerance; default is to
	                                 declare convergence on the first such iteration.
	        [--infile_W  (empty)]    Dense mxk matrix to initialize W; CSV file.
	                                 If unspecified, W will be randomly initialized.
	        [--infile_H  (empty)]    Dense kxn matrix to initialize H; CSV file.
	                                 If unspecified, H will be randomly initialized.
	        [--outfile_W  w.csv]     Filename for the W matrix result.
	        [--outfile_H  h.csv]     Filename for the H matrix result.
	        [--miniter  5]           Minimum number of iterations to perform.
	        [--maxiter  5000]        Maximum number of iterations to perform.
	        [--outprecision  6]      Write results with this many digits of precision.
	        [--maxthreads    3]      Upper limit to thread count.
	        [--normalize  1]         Whether to normalize W and scale H.
	                                     1 == yes, 0 == no
	        [--verbose  1]           Whether to print updates to the screen.
	                                     1 == print updates, 0 == silent

**5.** To test the installation of pysmallk, attempt to import numpy and pysmallk; numpy must be imported BEFORE pysmallk is imported. Running the following command from the command line should produce no output::

		python –c "import numpy; import pysmallk"
		
If there is no import error, pysmallk was installed correctly and is globally available.


**6.** When you are ready to shut down the VM, run `exit` from within the vagrant machine, then run one of the following from the command line of your host machine (wherever `vagrant up` was executed):

Save the current running state::

		vagrant suspend

Gracefully shut down the machine::

		vagrant halt

Remove the VM from your machine (this will require rebuilding the VM to restart it)::

		vagrant destroy

If you want to work with the VM again, from any of the above states you can run::

		vagrant up
		
again and the VM will be resumed or recreated.

Docker Instructions
===================

Running SmallK in a Docker container is intended for those who would like a fast, simple install that keeps their environment unmodified, in exchange for a loss in runtime performance. The basic process is to first build the Docker image, then run the Docker container to execute the desired command.

**1.** Install `Docker <https://docs.docker.com/engine/installation/>`_. If you are new to Docker, it may be worth exploring a `quick introduction <https://docs.docker.com/get-started/)>`_, or at least a `cheat-sheet <https://github.com/wsargent/docker-cheat-sheet>`_. There are `platform specific <https://docs.docker.com/manuals/>`_ installation, configuiration, and execution instructions for Mac, Windows, and Linux. The following instructions were tested on Ubuntu 16.04 with Docker version 17.06.0-ce.

**2.** Build the smallk Docker image.

First, make sure you have all submodules and their own submodules. From within the root of the smallk directory, run::

    	git submodule update --init --recursive

Now we can build the image. In the same (project root) directory, run this::

    	docker build -t smallk .

This will download all dependencies from the Ubuntu repositories, PyPI, GitHub, etc. Everything will be built including smallk itself. You will end up with a Docker image tagged "smallk". At the end of the build process you should see the following::

		Step 40/40 : CMD /bin/bash
		 ---> Running in 3fdb5e73afdc
		 ---> f8afa9f6a532
		Removing intermediate container 3fdb5e73afdc
		Successfully built f8afa9f6a532
		Successfully tagged smallk:latest

This can take as long as an hour to build the image, which is based on a minimal Ubuntu 16.04 installation. The `smallk/Dockerfile` can be customized in many ways to change the specifications for the image that is built. 

**3.** Run the Docker container.

The Docker container may be executed from any directory. Regardless of where you run it, you will need a volume for any input/output data. As an example, you may run the built-in PySmallk tests. The instructions below assume that your work directory is named `/home/ubuntu`. Replace it with the appropriate name. (The Docker daemon requires an absolute path for the local volume reference.)::

	    cd /home/ubuntu
	    git clone https://github.com/smallk/smallk_data.git smallk_data
	    docker run --volume /home/ubuntu/smallk_data:/data smallk make check PYSMALLK=1 ELEMVER=0.85 DATA_DIR=/data

Here is a breakdown of that Docker command to explain each part:

- `docker run`: Run a new container from an image

  - `--volume`: Add a volume (persistent storage area) to the container

    - `/home/ubuntu/smallk_data`: Local absolute path that will be exposed within the running container
    - `/data`: Internal path to use within the container

  - `smallk`: Image tag from which to spawn the new container
  - `make check PYSMALLK=1 ELEMVER=0.85`: Command to run within the container (run the smallk test suite)

     - `DATA_DIR=/data`: Tell the test suite where the local data is stored (from the perspective of the container)

If your execution of the PySmallk tests is successful, you should see a lot of output, ending with the following lines::

		assignment file test passed
		***** PysmallK: All tests passed. *****

***************************
Standard Build Instructions
***************************

Prerequisites
=============
* A modern C++ compiler that supports the C++11 standard, such as the latest release of the GNU or clang compilers
* `Elemental <http://libelemental.org/>`_, a high-performance library for dense, distributed linear algebra, which requires:

	* An MPI installation, such as `OpenMPI <http://www.open-mpi.org/software/ompi/v1.6/>`_ and `mpich <http://www.mpich.org/>`_
	* A BLAS implementation, preferably optimized/tuned for the local system
	* `libFLAME <http://www.cs.utexas.edu/~flame/web/libFLAME.html>`_: a high-performance library for dense linear algebra
	* `OpenMP <http://openmp.org/wp/>`_ (optional, see below)
	* CMake 

* Python 2.7, including the following libraries: 

	* numpy
	* scipy
	* cython version 0.22

Elemental can make use of MPI parallelization if available. This is generally advantageous for large problems. The SmallK code is also internally parallelized to take full advantage of multiple CPU cores for maximum performance. SmallK does not currently support distributed computation. However, future updates are planned that provide this capability. Please see the `About <http://smallk.github.io/about/>`_ page for information regarding distributed versions of many of the algorithms within SmallK.

We **strongly** recommend that users install both the HybridRelease and PureRelease builds of `Elemental <http://libelemental.org/>`_. OpenMP is enabled in the HybridRelease build and disabled in the PureRelease build. So why install both? For smaller problems the overhead of *MPI can actually cause code to run slower* than without it. Whereas for large problems MPI parallelization generally helps, but there is no clear transition point between where it helps and where it hurts. Thus, we encourage users to experiment with both builds to find the one that performs best for their typical problems.

We also recommend that users clearly separate the different build types as well as the versions of Elemental on their systems. Elemental is under active development, and new releases can introduce changes to the API that are not backwards-compatible with previous releases. To minimize build problems and overall hassle, we recommend that Elemental be installed so that the different versions and build types are cleanly separated.

Thus, two versions of Elemental need to be built. One is a hybrid release build with OpenMP parallelization, and the other is the pure release build without OpenMP parallelization. A separate build folder will be created for each build. The build that uses internal OpenMP parallelization is called a ‘HybridRelease’ build; the build that doesn’t is called a ‘PureRelease’ build. The debug build is called a ‘PureDebug’ build. The HybridRelease build is best for large problems, where the problem size is large enough to overcome the OpenMP parallel overhead. The following is for the 0.84 version of elemental. Set the version to that specified in the README.html file. Note that the files will be installed in `/usr/local/elemental/[version]/[build type]`.

**The SmallK software supports the latest stable release of Elemental, version 0.85 and above.**


.. CAUTION
   A note of caution: copying the command lines from this website and pasting them into a terminal may result in the commands not properly executing due to how characters are interpreted: the double dash --, “double quotes”, etc. For pasting the commands to a terminal, first copy the command lines to a text editor and copy/paste from there.

How to Install Elemental on MacOSX
----------------------------------

On MacOSX we recommend using `Homebrew <http://mxcl.github.io/homebrew/>`_ as the package manager. Homebrew does not require sudo privileges for package installation, unlike other package managers such as MacPorts. Thus the chances of corrupting vital system files are greatly reduced using Homebrew.

It is convenient to be able to view hidden files (like .file) in the MacOSX Finder. To do so run the following at the command line::

	defaults write com.apple.finder AppleShowAllFiles -bool YES

To revert back to hiding hidden files, set the Boolean flag to NO::

	defaults write com.apple.finder AppleShowAllFiles -bool NO

If you use Homebrew, ensure that your PATH is configured to search Homebrew’s installation directory first. Homebrew’s default installation location is `/usr/local/bin`, so that location needs to be first on your path. To check, run this command from a terminal window::

	cat /etc/paths

We also recommend running the following commands on a daily basis to refresh your brewed installations::

	brew update
	brew upgrade
	brew cleanup
	brew doctor

This will maintain your Homebrew installed software and diagnose any issues with the installations.

If the first entry is not `/usr/local/bin`, you will need to edit the `/etc/paths` file. This is a system file, so first create a backup. Move the line `/usr/local/bin` so that it is on the first line of the file. Save the file, then close the terminal session and start a new terminal session so that the path changes will take effect.

OSX:Install the latest GNU compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Elemental and SmallK both require a modern C++ compiler compliant with the C++11 standard. We recommend that you install the latest stable version of the clang and GNU C++ compilers. To do this, first install the XCode command line tools with this command::

		xcode-select --install

If this command produces an error, download and install XCode from the AppStore, then repeat the command. If that should still fail, install the command line tools from the XCode preferences menu. After the installation completes, run this command from a terminal window::

		clang++ --version

You should see output similar to this::

		Apple LLVM version 8.1.0 (clang-802.0.42)
		Target: x86_64-apple-darwin16.7.0
		Thread model: posix
		InstalledDir: /Library/Developer/CommandLineTools/usr/bin

The latest version of the GNU compiler at the time of writing is g++-7 (gcc 7.1.0), which is provided by the ‘gcc’ homebrew package. In addition to the gcc package, homebrew also provides a gcc49 package from the homebrew/versions tap. If this alternative gcc49 package is installed on your system it will prevent homebrew from symlinking the gcc package correctly. We recommend uninstalling the gcc49 versioned package and just using the gcc package instead. The Fortran compiler provided with the gcc package will also be configured to properly build numpy, which is required for the python interface to SmallK.

If you need to uninstall the gcc49 package, run the following commands::

		brew uninstall gcc49
		brew cleanup
		brew doctor

Then install the gcc package as follows::

		brew install gcc

The Apple-provided gcc and g++ will not be overwritten by this installation. The new compilers will be installed into /usr/local/bin as gcc-7, g++-7, and gfortran-6. The Fortran compiler is needed for the installation of MPI and for building the python interface to SmallK.

OSX:Install MPI Tools
^^^^^^^^^^^^^^^^^^^^^

Install the latest version of `mpich <http://www.mpich.org/>`_ with Homebrew as follows:

		brew install mpich

We recommend installing mpich rather than openMPI due to some superior features of mpich (prior versions of Elemental use openMPI, which can be installed using Homebrew as well). Also, Elemental 0.85 (discussed below) now uses mpich. Please see some discussion regarding openMPI vs mpich at:
http://stackoverflow.com/questions/2427399/mpich-vs-openmpi

OSX:Install libFlame
^^^^^^^^^^^^^^^^^^^^

Next we detail the installation of the high performance numerical library libflame. The library can be gotten from the libflame git repository on github.

It’s important to perform the git clone into a subdirectory NOT called ‘flame’ since this can cause name conflicts with the installation. Typically, a git clone is performed into a directory called ‘libflame’. However, other directory names will work as well. **Please do not use the directory name ‘flame’**.

To obtain the latest version of the FLAME library, clone the FLAME git repository with this command::

		git clone https://github.com/flame/libflame.git

Run the configure script in the top-level FLAME directory as follows (assuming the install path is `/usr/local/flame`)::

	./configure –-prefix=/usr/local/flame –-with-cc=/usr/local/bin/gcc-6 –-with-ranlib=/usr/local/bin/gcc-ranlib-6

A complete list of configuration options can be obtained by running `./configure –-help`.

After the configuration process completes, build the FLAME library as follows::

		make –j4

The –j4 option tells Make to use four processes to perform the build.  This number can be increased if you have a more capable system. Libflame will be installed with the following command::

		make install

The FLAME library is now installed.

OSX:Install Elemental
^^^^^^^^^^^^^^^^^^^^^

### Here is a recommended installation scheme for Elemental: ###

Choose a directory for the root of the Elemental installation.  For example, this may be::

		/usr/local/elemental

Download one of the SmallK-supported releases of Elemental, unzip and untar the distribution, and cd to the top-level directory of the unzipped distribution. This directory will be denoted by UNZIP_DIR in the following instructions.

**We now recommend using Elemental 0.85 or later. Earlier versions will no longer be supported**.

HybridRelease Build
"""""""""""""""""""

From the Elemental-0.85 directory, run the following command to create a local build directory for the HybridRelease build::

		mkdir build_hybrid
		cd build_hybrid

Use the following CMake command for the HybridRelease build, substituting 0.85 for <VERSION_STRING>::

	cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/<VERSION_STRING>/HybridRelease
	-D CMAKE_BUILD_TYPE=HybridRelease 
	-D CMAKE_CXX_COMPILER=/usr/local/bin/g++-7 
	-D CMAKE_C_COMPILER=/usr/local/bin/gcc-7 
	-D CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-7 
	-D MATH_LIBS="/usr/local/flame/lib/libflame.a;-framework Accelerate" 
	–D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..

Note that we have installed g++-7 into `/usr/local/bin` and libFLAME into `/usr/local/flame`. Alter these paths, if necessary, to match the installation location on your system.

Once the CMake configuration step completes, you can build Elemental from the generated Makefiles with the following command::

		make –j4

The –j4 option tells Make to use four processes to perform the build. This number can be increased if you have a more capable system.

After the build completes, install elemental as follows::

		make install

For Elemental version 0.85 and later, you need to setup your system to find the Elemental dynamic libraries. Method 2 below is preferred:

1. If your Mac OSX is earlier than Sierra, then, in your startup script (`~/.bash_profile`) or in a terminal window, enter the following command on a single line, replacing VERSION_STRING as above::

		export DYLD_LIBRARY_PATH=
			$DYLD_LIBRARY_PATH:/usr/local/elemental/VERSION_STRING/HybridRelease/lib/

2. If your Mac OSX is Sierra or higher Apple’s System Integrity Protection (SIP) will prevent using the `DYLD_LIBRARY_PATH` variable. We highly discourage disabling SIP as a workaround. Instead, in your startup script (`~/.bash_profile`) or in a terminal window, enter the following command on a single line, replacing `VERSION_STRING` as above::

		ln -s /usr/local/elemental/<VERSION_STRING>/HybridRelease/lib/*.dylib* /usr/local/lib

This will symlink the required Elemental libraries.


PureRelease Build
"""""""""""""""""

Run these commands to create a build directory for the PureRelease build::

		cd ..
		mkdir build_pure
		cd build_pure 

Then repeat the CMake configuration process, this time with the following command for the PureRelease build::

	cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/<VERSION_STRING>/PureRelease 
	-D CMAKE_BUILD_TYPE=PureRelease -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-7 
	-D CMAKE_C_COMPILER=/usr/local/bin/gcc-7 
	-D CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-7 
	-D MATH_LIBS="/usr/local/flame/lib/libflame.a;-framework Accelerate"  
	–D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON ..

Repeat the build commands and install this build of Elemental. 

For Elemental version 0.85 and later, you need to setup your system to find the Elemental dynamic libraries. Method 2 below is preferred:

1. If your Mac OSX is earlier than Sierra, then, in your startup script (`~/.bash_profile`) or in a terminal window, enter the following command on a single line, replacing VERSION_STRING as above::

		export DYLD_LIBRARY_PATH=
			$DYLD_LIBRARY_PATH:/usr/local/elemental/VERSION_STRING/HybridRelease/lib/

2. If your Mac OSX is Sierra or higher Apple’s System Integrity Protection (SIP) will prevent using the `DYLD_LIBRARY_PATH` variable. We highly discourage disabling SIP as a workaround. Instead, in your startup script (`~/.bash_profile`) or in a terminal window, enter the following command on a single line, replacing `VERSION_STRING` as above::

		ln -s /usr/local/elemental/<VERSION_STRING>/HybridRelease/lib/*.dylib* /usr/local/lib

This will symlink the required Elemental libraries.

The two builds of Elemental are now complete.

To test the installation, follow Elemental’s `test instructions <http://libelemental.org/documentation/0.85/build.html>`_  for the SVD test to verify that Elemental is working correctly.


How to Install Elemental on Linux
---------------------------------

We strongly recommend using a package manager for your Linux distribution for installation and configuration of the required dependencies. We cannot provide specific installation commands for every variant of Linux, so we specify the high-level steps below. The following was tested on a system with Ubuntu 16.04 installed.

Linux:Install the latest GNU compilers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend installation of the latest stable release of the GNU C++ compiler, which is g++-6 at the time of this writing. 

Also, install the latest version of GNU Fortran, which is needed for the installation of the Message Passing Interface (MPI) tools. 

Linux:Install MPI Tools
^^^^^^^^^^^^^^^^^^^^^^^

Elemental version 0.85 and higher uses `mpich <http://www.mpich.org/>`_ for its MPI implementation.::

		sudo apt-get update
		sudo apt-get install mpich

This completes the installation of the MPI tools. It should also be noted that the Open MP implementation of the MPI tools could also be used for the following installations.

Linux:Install libFlame
^^^^^^^^^^^^^^^^^^^^^^

Next we detail the installation of the high performance numerical library libflame. The library can be gotten from the libflame git repository on github.

It’s important to perform the git clone into a subdirectory NOT called ‘flame’ since this can cause name conflicts with the installation. We normally do a git clone into a directory called ‘libflame’. However, other directory names will work as well, but not ‘flame’.

To obtain the latest version of the FLAME library, clone the FLAME git repository with this command::

		git clone https://github.com/flame/libflame.git

Run the configure script in the top-level FLAME folder as follows (assuming you want to install to `/usr/local/flame`; if not, change the prefix path)::

		./configure --prefix=/usr/local/flame --with-cc=/usr/local/bin/gcc-6
			--with-ranlib=/usr/local/bin/gcc-ranlib-6

A complete list of configuration options can be obtained by running::

		./configure –-help

Then build and install the code as follows::

		make -j4
		make install

This completes the installation of the FLAME library.

Linux:Install an accelerated BLAS library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is essential to link Elemental with an accelerated BLAS library for maximum performance. Linking Elemental with a ‘reference’ BLAS implementation will cripple performance, since the reference implementations are designed for correctness not speed. 

If you do not have an accelerated BLAS on your system, you can download and build OpenBLAS. Download, unzip, and untar the tarball (version 0.2.19 as of this writing) and cd into the top-level folder. Build OpenBLAS with this command, assuming you have a 64-bit system::

		make BINARY=64 USE_OPENMP=1

Install with this command, assuming the installation directory is `/usr/local/openblas/0.2.19/`::

		make PREFIX=/usr/local/openblas/0.2.19/ install

This completes the installation of OpenBLAS.

Linux:Install Elemental
^^^^^^^^^^^^^^^^^^^^^^^

### Here is our suggested installation scheme for Elemental: ###

We strongly recommend that users install both the HybridRelease and PureRelease builds of Elemental. MPI tools are enabled in the HybridRelease build and disabled in the PureRelease build. So why install both? For smaller problems the overhead of MPI can actually cause code to run slower than without it. On the other hand, for large problems, MPI parallelization generally helps. However, there is no clear transition point between where it helps and where it hurts. Thus, we encourage users to experiment with both builds to find the one that performs best for their typical problems.

Another strong recommendation is that users clearly separate the different build types as well as the versions of Elemental on their systems. Elemental is under active development, and new releases can introduce changes to the API that are not backwards compatible with previous releases. To minimize build problems and overall hassle, we recommend that Elemental be installed so that the different versions and build types are cleanly separated.

Choose a directory for the root of the Elemental installation. A good choice is::

		/usr/local/elemental

Download one of the SmallK-supported releases of Elemental (see above), unzip and untar the distribution, and cd to the top-level folder of the unzipped distribution.  This directory will be denoted by UNZIP_DIR in the following instructions.

Note that Elemental version 0.85 or later is the version currently supported; earlier versions are not supported. If an earlier version is needed for Linux, use the following instructions.


For the first step of the installation, for Elemental versions prior to 0.85, we need to fix a few problems with the CMake configuration files.  Open the following file in a text editor::

		UNZIP_DIR/cmake/tests/OpenMP.cmake

On the first line of the file, change::

		if(HYBRID)

to this::

		if(ELEM_HYBRID)

Next, open this file in a text editor::

		UNZIP_DIR/cmake/tests/Math.cmake

Near the first line of the file, change::

		if(PURE)

to this::

		if(ELEM_PURE)

Save both files.

Run these commands to create the required directories for the build types::

		mkdir build_hybrid
		mkdir build_pure

HybridRelease build
"""""""""""""""""""

From the Elemental-<VERSION> folder, run the following command to change to the local build folder for the HybridRelease build::

		cd build_hybrid

For the first step of the installation, we need to fix a few problems with the CMake configuration files. Open the following file in a text editor::

		Elemental-<VERSION>/cmake/tests/OpenMP.cmake

On the first line of the file, change::

		if(HYBRID)

to this::

		if(ELEM_HYBRID)

Next, open this file in a text editor::

		Elemental-<version>/cmake/tests/Math.cmake

Near the first line of the file, change::

		if(PURE)

to this::

		if(ELEM_PURE)

Save both files.

Run the following command to create a local build folder for the HybridRelease build::

		cd build_hybrid

Use the following CMake command for the HybridRelease build::

	cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/<VERSION>/HybridRelease
	-D CMAKE_BUILD_TYPE=HybridRelease -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-6
	-D CMAKE_C_COMPILER=/usr/local/bin/gcc-6
	-D CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-6
	-D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/openblas/0.2.19/ –lopenblas –lm"
	–D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..

Note that we have installed g++-6 into `/usr/local/bin` and libFLAME into `/usr/local/flame`. Alter these paths, if necessary, to match the installation location on your system.

If this command does not work on your system, you may need to define the BLAS_LIBS and/or GFORTRAN_LIB config options.

Version 0.85 of Elemental has an error in one of its cmake files. The file is::

		Elemental-0.85/cmake/tests/CXX.cmake

Modify the first line of this file from::

		include(FindCXXFeatures)

to::

		include_directories(FindCXXFeatures)

since FindCXXFeatures is now a directory. After this change, Elemental should Make without errors.

Once the CMake configuration step completes, you can build Elemental from the generated Makefiles with the following command::

		make –j4

The –j4 option tells Make to use four processes to perform the build. This number can be increased if you have a more capable system.

After the build completes, install elemental as follows::

		make install


After installing Elemental version 0.85, setup the system to find the Elemental shared library.  Either in the startup script (~/.bashrc) or in a terminal window, enter the following command on a single line, replacing VERSION_STRING as above::

		export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/elemental/VERSION_STRING/HybridRelease/lib/


PureRelease build
"""""""""""""""""

After this, run these commands to create a build folder for the PureRelease build::

		cd ..
		cd build_pure

Then repeat the CMake configuration process, this time with the following command for the PureRelease build::

	cmake -D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.84-p1/PureRelease
	-D CMAKE_BUILD_TYPE=PureRelease -D CMAKE_CXX_COMPILER=/usr/local/bin/g++-6
	-D CMAKE_C_COMPILER=/usr/local/bin/gcc-6
	-D CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran-6
	-D MATH_LIBS="/usr/local/flame/lib/libflame.a;-L/usr/local/openblas/0.2.19/ –lopenblas –lm"
	–D ELEM_EXAMPLES=ON –D ELEM_TESTS=ON  ..

If this command does not work on your system, you may need to define the `BLAS_LIBS` and/or `GFORTRAN_LIB` config options.

Repeat the build commands and install this build of Elemental. Then, if you installed a version of Elemental **prior** to the 0.84 release, edit the `/usr/local/elemental/<version>/PureRelease/conf/ElemVars` file and replace the CXX line as indicated above.

Version 0.85 of Elemental has an error in one of its cmake files. The file is::

		Elemental-0.85/cmake/tests/CXX.cmake

Modify the first line of this file from::

		include(FindCXXFeatures)

to::

		include_directories(FindCXXFeatures)

since FindCXXFeatures is now a directory. After this change, Elemental should Make without errors.

If Elemental version 0.85 or later was installed, setup the system to find the Elemental shared library for the PureRelease build. Enter the following command in a terminal window on a single line, replacing `VERSION_STRING` as above::

	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/elemental/VERSION_STRING/PureRelease/lib/

Note: set this variable to point to either the HybridRelease or the PureRelease build of the Elemental shared library whenever you want to use SmallK.

This completes the two builds of Elemental.

To test the installation, follow Elemental’s `test instructions <http://libelemental.org/documentation/\<version\>/build.html>`_ for the SVD test to verify that Elemental is working correctly.


Installation of Python libraries
--------------------------------

**Note: the following section for installing the Python libraries can be skipped if not needed.**

OSX:Install Python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install Python scientific packages
""""""""""""""""""""""""""""""""""

Assuming that you have used brew to install gcc, as indicated earlier, you can run the following commands to install the necessary libraries::

		brew install python
		brew install numpy
		brew install scipy

To check your installation, run::

		brew test numpy

IMPORTANT: Check to see that your numpy installation has correctly linked to the needed BLAS libraries.

Ensure that you are running the correct python::

		which python

This should print out `/usr/local/bin/python`. Open a python terminal by typing `python` at the command line and run the following::

		import numpy as np
		np.__config__.show()

You should see something similar to the following::

		lapack_opt_info:
		extra_link_args = ['-Wl,-framework', '-Wl,Accelerate']
		extra_compile_args = ['-msse3']
		define_macros = [('NO_ATLAS_INFO', 3)]

		blas_opt_info:
		extra_link_args = ['-Wl,-framework', '-Wl,Accelerate']
		extra_compile_args = ['-msse3', '-I/System/Library/Frameworks/vecLib.framework/Header’]
		define_macros = [('NO_ATLAS_INFO', 3)]

If you are using OpenBLAS, you should see that indicated as well.

Install Cython: a Python interface to C/C++
"""""""""""""""""""""""""""""""""""""""""""

First install the Python Package Index utility, pip. Many Python packages are configured to use this package manager, Cython being one.::

		brew install pip

Only Cython 0.22 is supported at this time. To check which version is installed on your system use this commands::

		$ python
		>> import Cython
		>> Cython.__version__
		'0.22'
		>> 

To install Cython version 0.22 (if not already installed)::

		pip uninstall cython
		pip install cython==0.22

Check the version of cython as above to ensure that Cython version 0.22 is installed.


Linux:Install Python libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Python libraries can easily be installed via pip and apt-get with the following commands::

		apt-get install pip
		pip install numpy
		apt-get install python-scipy
		pip uninstall cython
		pip install cython==0.22

This also ensures that cython version 0.22 is installed, which is the currently supported version. The Makefile assumes an installation path of `/usr/local/lib/python2.7/site-packages` for the compiled library file. If you are not using apt-get to install your packages, you will need to tell the Makefile where the appropriate site-packages directory is located on your system. Setting the `SITE_PACKAGES_DIR` command line variable when running make accomplishes this. If this doesn’t work, an alternative way to set this up is to add a line to the `.bash_profile` file (always back up first)::

		export SITE_PACKAGES_DIR="<path to lib/python2.7>/site-packages/"

This allows for special installations of Python such as Continuum Analytics’ `Anaconda <https://www.continuum.io/>`_ distribution site-packages to be accessed.


Build and Installation of SmallK
================================

Obtain the Source Code
----------------------

The source code for the SmallK library can be downloaded from the `SmallK repository <https://github.com/smallk/smallk.github.io/tree/master/code>`_ on github.
Once downloaded uncompress the tar ball and follow the installation instructions below.

Build the SmallK library
------------------------

After downloading and unpacking the code tarball cd into the top-level libsmallk1_<version> directory, where version is MAJOR.MINOR.PATCH (for example 1.6.2). The makefiles assume that you followed our suggested installation plan for Elemental. If this is NOT the case you will need to do one of the following:

		1. Create an environment variable called ELEMENTAL_INSTALL_DIR which contains the 
			path to the root folder of your Elemental installation
		2. Define the variable ELEMENTAL_INSTALL_DIR on the make command line
		3. Edit the SmallK makefile so that it can find your Elemental installation

Assuming that the default install locations are acceptable, build the SmallK code by running this command from the root directory of the distribution::

		make all PYSMALLK=1 ELEMVER=0.85

or::

		make all PYSMALLK=0 ELEMVER=0.85

This will build the SmallK and pysmallk (optional; see section [Installation of Python libraries]) below for setup of the Python libraries) libraries and several command-line applications. These are:

	1. libsmallk.a, the SmallK library
	2. preprocess_tf, a command-line application for processing and scoring term-frequency matrices
	3. matrixgen, a command-line application for generating random matrices
	4. nmf, a command-line application for NMF
	5. hierclust, a command-line application for fast hierarchical clustering
	6. flatclust, a command-line application for flat clustering via NMF
	7. pysmallk.so, if PYSMALLK=1 (0: default), the Python-wrapped SmallK library, making SmallK available via Python

Install the SmallK library
--------------------------

To install the code, run this command to install to the default location, which is `/usr/local/smallk`::

		make install PYSMALLK=1 ELEMVER=0.85

or::

		make install PYSMALLK=0 ELEMVER=0.85

This will install the binary files listed above into the `/usr/local/smallk/bin` directory, which needs to be on your path to run the executables from anywhere on your system and avoid prepending with the entire path. To install the binary code to a different location, either create an environment variable called `SMALLK_INSTALL_DIR` and set it equal to the desired installation location prior to running the install command, or supply a prefix argument::

		make prefix=/path/to/smallk  install

If `PYSMALLK=1`, this will install pysmallk.so into the site-packages directory associated with the Python binary, which is determined by ‘brew install python’ as discussed above or wherever the python distribution is installed on the system, e.g., `Continuum’s Anaconda Python <https://www.continuum.io/>`_ distribution is installed in the user’s home directory. To install the Python library to a different location, create an environment variable called `SITE_PACKAGES_DIR` and set it equal to the desired installation location prior to running the install command, or supply this as an argument for make::

		make SITE_PACKAGES_DIR=/path/to/site-packages install

Or, as a last resort, you can edit the top-level SmallK makefile to conform to the installation scheme of your system.  You may need root privileges to do the installation, depending on where you choose to install it.

Before testing the installation, the test code needs to access data. The data is located in a separate github repository so that when cloning the code, the large amount of data is not included. The data repository is located on github at `smallk_data <https://github.com/smallk/smallk_data>`_:

Check the build and installation
--------------------------------

To test the build, run this command with DATA_DIR set to wherever the SmallK data repository was cloned::

		make check PYSMALLK=1 ELEMVER=0.85 DATA_DIR=../smallk_data

or::

		make check PYSMALLK=0 ELEMVER=0.85 DATA_DIR=../smallk_data

This will run a series of tests, none of which should report a failure.  Sample output from a run of these tests can be found in section `SmallK Test Results <http://smallk.github.io/documentation/tests/#smalk_tests>`_.

Note: if you installed Elemental version 0.85, you will need to configure your system to find the Elemental shared library.  See the Elemental installation instructions above for information on how to do this.

The command-line applications can be built individually by running the appropriate make command from the top-level SmallK directory.  These commands are::

		To build the smallk library only: 		make libsmallk
		To build the preprocessor only:			make preprocessor
		To build the matrix generator only:		make matrixgen
		To build the nmf only:					make nmf
		To build hierclust only:				make hierclust
		To build flatclust only:				make flatclust
		To build pysmallk only:					make pysmallk

This completes the SmallK NMF library installation.

.. NOTE
   Note: Pysmallk requires builds of libsmallk, preprocessor, matrixgen, hierclust, and flatclust.

Build and Installation of pysmallk shared library
=================================================

Before building pysmallk, you must ensure that you have already built the standard SmallK library and applications: libsmallk, preprocessor, matrixgen, hierclust, and flatclust.

All C++ and python libraries and applications can be built simultaneously by setting the PYSMALLK command line variable::

		make PYSMALLK=1

To build pysmallk individually from the pysmallk subdirectory (`<path to SmallK>/libsmallk-<version>/pysmallk`)::
		
		make pysmallk
		
To check the library installation::

		make pysmallk_check DATA_DIR=../smallk_data

This will run a series of tests, none of which should report a failure.

To install the shared library in a globally accessible location, enable the `PYSMALLK` command line variable and, if needed, specify an `INSTALLATION_DIR`.

The Makefile assumes an installation path of `/usr/local/lib/python2.7/site-packages` for the compiled library file. If you are not using brew to install your packages, you will need to tell the Makefile where the appropriate site-packages directory is located on your system. Setting the `INSTALLATION_DIR` command line variable when running make accomplishes this. Also, make sure that there is not another site-packages directory in your PATH before the site-packages you intend to use since `make install` will copy pysmallk.so to `/usr/local/lib/python2.7/site-packages` by default. Other Python distributions will probably interfere with the pysmallk installation.::

		make install PYSMALLK=1	INSTALLATION_DIR=/usr/local/lib/python2.7/site-packages/

To uninstall the libraries::

		make uninstall PYSMALLK=1	INSTALLATION_DIR=/usr/local/lib/python2.7/site-packages/

*******************
Matrix file formats
*******************

The SmallK software supports comma-separated value (CSV) files for dense matrices and `Matrix Market <http://math.nist.gov/MatrixMarket/formats.html>`_ files for sparse matrices.

For example, the 5x3 dense matrix::

		42	47	52
		43	48	53
		44	49	54
		45	50	55
		46	51	56

would be stored in a CSV file as follows::
	
		42,47,52
		43,48,53
		44,49,54
		45,50,55
		46,51,56

The matrix is loaded exactly as it appears in the file.  **Internally, SmallK stores dense matrices in column-major order**. Sparse matrices are stored in **compressed column format**.


**********
Disclaimer
**********

This software is a work in progress. It will be updated throughout the course of the XDATA program with additional algorithms and examples. The distributed NMF factorization routine uses sequential algorithms, but it replaces the matrices and matrix operations with distributed versions. The GA Tech research group is working on proper distributed NMF algorithms, and when such algorithms are available they will be added to the library. Thus, the performance of the distributed code should be viewed as being the baseline for our future distributed NMF implementations.

************
Contact Info
************

For comments, questions, bug reports, suggestions, etc., contact:

Barry Drake 
Research Scientist 
Information and Communications Laboratory (ICL) 
Information and Cyber Sciences Directorate (ICSD) 
Georgia Tech Research Institute (GTRI) 
75 5TH St. NW STE 900 
ATLANTA, GA 30308-1018
barry.drake@gtri.gatech.edu

Stephen Lee-Urban
Research Scientist
Information and Communications Laboratory (ICL)
Information and Cyber Sciences Directorate (ICSD)
Georgia Tech Research Institute (GTRI)
75 5TH St. NW STE 900
ATLANTA, GA 30308-1018
stephen.lee-urban@gtri.gatech.edu




