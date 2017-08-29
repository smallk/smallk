#########################
Quickstart - Installation
#########################

***********************
Vagrant Virtual Machine
***********************

Installing SmallK into a virtual machine (OSX, Linux, Windows) is intended for those who are not doing development and/or do not have a reason to do the full installation on Linux or OSX outlined in the sections to follow.

The complete stack of software dependencies for SmallK as well as SmallK itself can be rapidly set up and configured through use of Vagrant and VirtualBox and the files included in the repository. The Vagrant install has been tested on Linux Ubuntu 16.04, Mac OSX Sierra 10.12.6, and Windows 10. 

Note that the ``smallk/vagrant/bootstrap.sh`` file can be modified to perform various tasks when provisioning the vagrant session. Consider customizing ``bootstrap.sh`` to set up a custom install of libsmallk as required.

To deploy the SmallK VM:

**1.** Install `Vagrant <http://www.vagrantup.com/downloads.html>`_ and `VirtualBox <https://www.virtualbox.org/wiki/Downloads>`_.

.. tip::
   Note: For Windows, ensure that you have a VirtualBox version >= 4.3.12. After installing Vagrant, you may need to log out and log back in to ensure that you can run vagrant commands in the command prompt.

**Optional:** ``git clone`` the `smallk_data <https://github.com/smallk/smallk_data>`_ repository so that it is parallel with the smallk repository. This is an alternate way to test the installation and begin to work with SmallK. This directory can be synced with a directory of the same name in the VM by adding or uncommenting the following line in ``smallk/vagrant/Vagrantfile``::

		config.vm.synced_folder "../../smallk_data", "/home/vagrant/smallk_data"	

**2.** From within the ``smallk/vagrant`` directory, run::
		
		vagrant up
		
This can take as long as an hour to build the VM, which will be based on a minimal Ubuntu 16.04 installation. The ``smallk/vagrant/Vagrantfile`` can be customized in many ways to change the specifications for the VM that is built. See more information `here <http://docs.vagrantup.com/v2/>`_. The default configuration provides the VM with 4 GB of memory and 3 CPUs. Increasing these allocations will improve the performance of the application. This can be done by modifying these lines in the Vagrantfile::

		vb.memory = 4096
		vb.cpus = 3

After ``vagrant up`` has completed, the SmallK and pysmallk libraries will have been built and tested. Additionally, the smallk_data directory, if cloned as in the optional step above, will have been synced into the VM. For more details regarding what is being built and executed while provisioning the VM, please inspect ``smallk/vagrant/bootstrap.sh``.

**3.** Once the VM has been built, run::

		vagrant ssh

.. tip::
   Note: For Windows, you will need an ssh client in order to run the above command. This can be obtained via `CygWin <https://www.cygwin.com/>`_ `MinGW <http://sourceforge.net/projects/mingw/files/>`_, or `Git <http://git-scm.com/downloads>`_. If you would like to use PuTTY to connect to your virtual machine, follow `these <https://github.com/Varying-Vagrant-Vagrants/VVV/wiki/Connect-to-Your-Vagrant-Virtual-Machine-with-PuTTY>`_ instructions.

In case you need it, the username/password for the VM created will be vagrant/vagrant.

This will drop you into the command line of the VM that was just created, in a working directory at ``/home/vagrant``. From there, you can navigate to ``/home/vagrant/libsmallk-<version>``, (e.g., libsmallk-1.6.2), and run::

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

		python -c "import numpy; import pysmallk"
		
If there is no import error, pysmallk was installed correctly and is globally available.


**6.** When you are ready to shut down the VM, run ``exit`` from within the vagrant machine, then run one of the following from the command line of your host machine (wherever ``vagrant up`` was executed):

Save the current running state::

		vagrant suspend

Gracefully shut down the machine::

		vagrant halt

Remove the VM from your machine (this will require rebuilding the VM to restart it)::

		vagrant destroy

If you want to work with the VM again, from any of the above states you can run::

		vagrant up
		
again and the VM will be resumed or recreated.

*******************
Docker Instructions
*******************

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

This can take as long as an hour to build the image, which is based on a minimal Ubuntu 16.04 installation. The ``smallk/Dockerfile`` can be customized in many ways to change the specifications for the image that is built. 

**3.** Run the Docker container.

The Docker container may be executed from any directory. Regardless of where you run it, you will need a volume for any input/output data. As an example, you may run the built-in PySmallk tests. The instructions below assume that your work directory is named ``/home/ubuntu``. Replace it with the appropriate name. (The Docker daemon requires an absolute path for the local volume reference.)::

	    cd /home/ubuntu
	    git clone https://github.com/smallk/smallk_data.git smallk_data
	    docker run --volume /home/ubuntu/smallk_data:/data smallk make check PYSMALLK=1 ELEMVER=0.85 DATA_DIR=/data

Here is a breakdown of that Docker command to explain each part:

- ``docker run``: Run a new container from an image

  - ``--volume``: Add a volume (persistent storage area) to the container

    - ``/home/ubuntu/smallk_data``: Local absolute path that will be exposed within the running container
    - ``/data``: Internal path to use within the container

  - ``smallk``: Image tag from which to spawn the new container
  - ``make check PYSMALLK=1 ELEMVER=0.85``: Command to run within the container (run the smallk test suite)

     - ``DATA_DIR=/data``: Tell the test suite where the local data is stored (from the perspective of the container)

If your execution of the PySmallk tests is successful, you should see a lot of output, ending with the following lines::

		assignment file test passed
		***** PysmallK: All tests passed. *****

