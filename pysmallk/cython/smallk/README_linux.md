With some additional notes for Linux specifics

Building Open MPI
-----------------
After downloading the latest version of Open MPI, and untarring the directory, execute the following to build:

*Note the --disable-dlopen, this may/may not be needed.

	./configure --disable-dlopen --enable-mpi-thread-multiple --prefix=/usr/local  CC=/usr/local/bin/gcc-4.8 CXX=/usr/local/bin/g++-4.8 F77=/usr/local/bin/gfortran  FC=/usr/local/bin/gfortran

Without the --disable-dlopen, Python may not be able to use Open MPI.

Building the Elemental Library
------------------------------
Supposing you have all the dependencies of the Elemental library (Open MPI, BLAS, libflame), use the following to build the elemental library, making the appropriate changes to locate your compilers and libraries:

(on one line) 

	cmake 
	-D CMAKE_INSTALL_PREFIX=/usr/local/elemental/0.81/HybridRelease 
	-D CMAKE_BUILD_TYPE=HybridRelease 
	-D CMAKE_CXX_COMPILER=/usr/local/bin/g++-4.8 
	-D CMAKE_CXX_FLAGS=-fPIC
	-D CMAKE_C_COMPILER=/usr/local/bin/gcc-4.8
	-D CMAKE_C_fPIC=-FLAGS
	-D CXX_FLAGS=”-std=c++11 –O3 -fPIC” 
	-D CMAKE_Fortran_COMPILER=/usr/local/bin/gfortran 
	-D MATH_LIBS="/usr/local/flame/lib/libflame.a; -L/usr/local/lib –lopenblas -lm"
	-D ELEM_EXAMPLES=ON
	-D ELEM_TESTS=ON  ..
	
Then
 
	make
	make install

Edit the Elemental configuration file as follows:

	cd /usr/local/elemental-0.81-HybridRelease/conf/

Open the <i>elemvariables</i> file in a text editor and add the string -std=c++11 to the CXX macro line. For instance, if the CXX macro is

	CXX = /usr/local/bin/g++-4.8

change it to

	CXX = /usr/local/bin/g++-4.8 -std=c++11
	
Without the -fPIC flag you may have trouble building the Cython shared object

Building the Smallk Library
------------------------
From the main xdata2 folder, edit the Makefile to add the -fPIC flag into the CXXFLAGS, similar to the Elemental build, then run these commands:

	make
	make install

In your home directory, there should now be a 'smallk' folder.  This folder contains 'include',' conf', and 'lib' folders.  You will note that libsmallk.a is in the lib folder.  If you see this, you are looking at the 'installed' version of smallk.

If you run this command from the xdata2 folder:

	make check

It should build and run the smallk test code.  IF the tests pass you should see text at the end stating that the W and H tests passed.

To uninstall the code and clean up everything, do this from the xdata2 folder:

	make clean
	make uninstall

The items in ~/smallk should now be gone.  The folders will remain, though.

Building the Python shared object libsmallkpy.so
-------------------------------------------
In the xdata2/smallk/cython\_Linux folder, open the setup.py script, and check that <i>hpp\_path</i> points to the smallk header, and that <i>archive_path</i> points to the smallk library libsmallk.a, which if you followed the default build above, will be in $HOME/smallk/include and $HOME/smallk/lib, respectively. Then run

	python setup.py build_ext --inplace
	
Here's where you might see something like

	/path/to/smallk/lib/libsmallk.a(smallk.o): relocation R_X86_64_32 against `.rodata.str1.1' can not be used when making a shared object; recompile with -fPIC
	
Prompting the rebuild of the smallk library with the -fPIC flag.

After a successful build, you should now have a libsmallkpy.so file that you can import into any Python script. 

Testing
-------
To run the included test script, execute:

	python test_smallk.py --k 2 --algorithm MU --data_dir ../../data/
	
and check that your output is similar to:

	$ python test_smallk.py --k 2 --algorithm MU --data_dir ../../data/
	sys.argv =  ['test_smallk.py', '--k', '2', '--algorithm', 'MU', '--data_dir', '../../data/']
	Initialization oK
	Loading matrix...
	Running NMF on w (../../data/test/W.csv) and h (../../data/test/H.csv)
	Initializing matrix W...
	Initializing matrix H...

		             parameters: 

			      algorithm: Multiplicative Updating
		stopping criterion: Relative Change in the F-norm of W
			         height: 12411
			          width: 7984
			              k: 2
			        miniter: 5
			        maxiter: 5000
			            tol: 0.005
			  success_count: 1
			     matrixfile: ../../data/reuters.mtx
			     maxthreads: 4

	1:	progress metric: 	(min_iter)
	2:	progress metric: 	(min_iter)
	3:	progress metric: 	(min_iter)
	4:	progress metric: 	(min_iter)
	5:	progress metric: 	(min_iter)
	6:	progress metric:	0.763943
	7:	progress metric:	0.0844044
	8:	progress metric:	0.0614892
	9:	progress metric:	0.0433576
	10:	progress metric:	0.0291242
	11:	progress metric:	0.019805
	12:	progress metric:	0.0141457
	13:	progress metric:	0.0106286
	14:	progress metric:	0.00830373
	15:	progress metric:	0.00667069
	16:	progress metric:	0.00548692
	17:	progress metric:	0.00465811

	Solution converged after 17 iterations.

	Norms for matrix W: 
	|| W ||_max = 0.391874
	|| W ||_1   = 13.4866
	|| W ||_oo  = 0.576506
	|| W ||_F   = 1.41421

	Norms for matrix H: 
	|| H ||_max = 0.856171
	|| H ||_1   = 0.880084
	|| H ||_oo  = 989.695
	|| H ||_F   = 23.6168

	Elapsed wall clock time: 0.490 sec.

