#########
SmallkAPI
#########

********
Contents
********

.. toctree::
   :maxdepth: 8

Examples of API Usage
=====================

In the examples folder you will find a file called smallk_example.cpp. This file contains several examples of how to use the SmallK library.  Also included in the examples folder is a makefile that you can customize for your use.  Note that the SmallK library must first be installed before the example project can be built.

As an example of how to use the sample project, assume the smallk software has been installed into /usr/local/smallk.  Also assume that the user chose to create the recommended environment variable SMALLK_INSTALL_DIR that stores the location of the top-level install folder, i.e. the user’s .bashrc file contains this statement::

		export SMALLK_INSTALL_DIR=/usr/local/smallk 

To build the smallk example project, open a terminal window and cd to the smallk/examples folder and run this command:: 

		make

To run the example project, run this command::

		./bin/example ../../smallk_data

Note: the output will be *similar* to the following not identical since some problems are randomly initialized::

	Smallk major version: 1
	Smallk minor version: 0
	Smallk patch level:   0
	Smallk version string: 1.0.0
	Loading matrix...

	************************************************************
	*                                                          *
	*                Running NMF-BPP using k=32                *
	*                                                          *
	************************************************************
	Initializing matrix W...
	Initializing matrix H...

                parameters: 

             algorithm: Nonnegative Least Squares with Block Principal Pivoting
    stopping criterion: Ratio of Projected Gradients
                height: 12411
                 width: 7984
                     k: 32
               miniter: 5
               maxiter: 5000
                   tol: 0.005
            matrixfile: ../data/reuters.mtx
            maxthreads: 8

	1:  progress metric:    (min_iter)
	2:  progress metric:    (min_iter)
	3:  progress metric:    (min_iter)
	4:  progress metric:    (min_iter)
	5:  progress metric:    (min_iter)
	6:  progress metric:    0.0747031
	7:  progress metric:    0.0597987
	8:  progress metric:    0.0462878
	9:  progress metric:    0.0362883
	10: progress metric:    0.030665
	11: progress metric:    0.0281802
	12: progress metric:    0.0267987
	13: progress metric:    0.0236731
	14: progress metric:    0.0220778
	15: progress metric:    0.0227083
	16: progress metric:    0.0244029
	17: progress metric:    0.0247552
	18: progress metric:    0.0220007
	19: progress metric:    0.0173831
	20: progress metric:    0.0137033

	Solution converged after 39 iterations.

	Elapsed wall clock time: 4.354 sec.

	Writing output files...

	************************************************************
	*                                                          *
	*               Running NMF-HALS using k=16                *
	*                                                          *
	************************************************************
	Initializing matrix W...
	Initializing matrix H...

                parameters: 

             algorithm: HALS
    stopping criterion: Ratio of Projected Gradients
                height: 12411
                 width: 7984
                     k: 16
               miniter: 5
               maxiter: 5000
                   tol: 0.005
            matrixfile: ../data/reuters.mtx
            maxthreads: 8

	1:  progress metric:    (min_iter)
	2:  progress metric:    (min_iter)
	3:  progress metric:    (min_iter)
	4:  progress metric:    (min_iter)
	5:  progress metric:    (min_iter)
	6:  progress metric:    0.710219
	7:  progress metric:    0.580951
	8:  progress metric:    0.471557
	9:  progress metric:    0.491855
	10: progress metric:    0.531999
	11: progress metric:    0.353302
	12: progress metric:    0.201634
	13: progress metric:    0.1584
	14: progress metric:    0.142572
	15: progress metric:    0.12588
	16: progress metric:    0.113239
	17: progress metric:    0.0976934
	18: progress metric:    0.0821207
	19: progress metric:    0.0746089
	20: progress metric:    0.0720616
	40: progress metric:    0.0252854
	60: progress metric:    0.0142085
	80: progress metric:    0.0153269

	Solution converged after 88 iterations.

	Elapsed wall clock time: 1.560 sec.

	Writing output files...

	************************************************************
	*                                                          *
	*       Running NMF-RANK2 with W and H initializers        *
	*                                                          *
	************************************************************
	Initializing matrix W...
	Initializing matrix H...

                parameters: 

             algorithm: Rank 2
    stopping criterion: Ratio of Projected Gradients
                height: 12411
                 width: 7984
                     k: 2
               miniter: 5
               maxiter: 5000
                   tol: 0.005
            matrixfile: ../data/reuters.mtx
            maxthreads: 8

	1:  progress metric:    (min_iter)
	2:  progress metric:    (min_iter)
	3:  progress metric:    (min_iter)
	4:  progress metric:    (min_iter)
	5:  progress metric:    (min_iter)
	6:  progress metric:    0.0374741
	7:  progress metric:    0.0252389
	8:  progress metric:    0.0169805
	9:  progress metric:    0.0113837
	10: progress metric:    0.00761077
	11: progress metric:    0.0050782
	12: progress metric:    0.00338569

	Solution converged after 12 iterations.

	Elapsed wall clock time: 0.028 sec.

	Writing output files...

	************************************************************
	*                                                          *
	*       Repeating the previous run with tol = 1.0e-5       *
	*                                                          *
	************************************************************
	Initializing matrix W...
	Initializing matrix H...

                parameters: 

             algorithm: Rank 2
    stopping criterion: Ratio of Projected Gradients
                height: 12411
                 width: 7984
                     k: 2
               miniter: 5
               maxiter: 5000
                   tol: 1e-05
            matrixfile: ../data/reuters.mtx
            maxthreads: 8

	1:  progress metric:    (min_iter)
	2:  progress metric:    (min_iter)
	3:  progress metric:    (min_iter)
	4:  progress metric:    (min_iter)
	5:  progress metric:    (min_iter)
	6:  progress metric:    0.0374741
	7:  progress metric:    0.0252389
	8:  progress metric:    0.0169805
	9:  progress metric:    0.0113837
	10: progress metric:    0.00761077
	11: progress metric:    0.0050782
	12: progress metric:    0.00338569
	13: progress metric:    0.00225761
	14: progress metric:    0.00150429
	15: progress metric:    0.00100167
	16: progress metric:    0.000666691
	17: progress metric:    0.000443654
	18: progress metric:    0.000295213
	19: progress metric:    0.000196411
	20: progress metric:    0.000130604

	Solution converged after 27 iterations.

	Elapsed wall clock time: 0.061 sec.

	Writing output files...
	Minimum value in W matrix: 0.
	Maximum value in W matrix: 0.397027.


	************************************************************
	*                                                          *
	*      Running HierNMF2 with 5 clusters, JSON format       *
	*                                                          *
	************************************************************
	loading dictionary...
	creating random W initializers...
	creating random H initializers...

            parameters: 

                height: 12411
                 width: 7984
            matrixfile: ../data/reuters.mtx
              dictfile: ../data/reuters_dictionary.txt
                   tol: 0.0001
               miniter: 5
               maxiter: 5000
              maxterms: 5
            maxthreads: 8
	[1] [2] [3] [4] 

	Elapsed wall clock time: 391 ms.
	9/9 factorizations converged.

	Writing output files...

	************************************************************
	*                                                          *
	* Running HierNMF2 with 10 clusters, 12 terms, XML format  *
	*                                                          *
	************************************************************
	creating random W initializers...
	creating random H initializers...

            parameters: 

                height: 12411
                 width: 7984
            matrixfile: ../data/reuters.mtx
              dictfile: ../data/reuters_dictionary.txt
                   tol: 0.0001
               miniter: 5
               maxiter: 5000
              maxterms: 12
            maxthreads: 8
	[1] [2] [3] [4] [5] [6] dropping 20 items ...
	[7] [8] [9] 

	Elapsed wall clock time: 837 ms.
	21/21 factorizations converged.

	Writing output files...

	************************************************************
	*                                                          *
	*  Running HierNmf2 with 18 clusters, 8 terms, with flat   *
	*                                                          *
	************************************************************
	creating random W initializers...
	creating random H initializers...

            parameters: 

                height: 12411
                 width: 7984
            matrixfile: ../data/reuters.mtx
              dictfile: ../data/reuters_dictionary.txt
                   tol: 0.0001
               miniter: 5
               maxiter: 5000
              maxterms: 8
            maxthreads: 8
	[1] [2] [3] [4] [5] [6] dropping 20 items ...
	[7] [8] [9] dropping 25 items ...
	[10] [11] [12] [13] [14] [15] [16] [17] 

	Running NNLS solver...
	1:  progress metric:    1
	2:  progress metric:    0.264152
	3:  progress metric:    0.0760648
	4:  progress metric:    0.0226758
	5:  progress metric:    0.00743562
	6:  progress metric:    0.00280826
	7:  progress metric:    0.00103682
	8:  progress metric:    0.000361738
	9:  progress metric:    0.000133087
	10: progress metric:    5.84849e-05

	Elapsed wall clock time: 1.362 s.
	40/40 factorizations converged.

	Writing output files...

The output files are written to the default directory or the directory specified on the command line.

SmallK API
==========

The SmallK API is an extremely simplistic API for basic NMF and clustering.  Users who require more control over the factorization or clustering algorithms can instead run one of the command-line applications in the SmallK distribution.

The SmallK API is exposed by the file smallk.hpp, which can be found in this location:: 

		SMALLK_INSTALL_DIR/include/smallk.hpp.  

All API functions are contained within the smallk namespace. 

An example of how to use the API can be found in the file examples/smallk_example.cpp.

The smallk library maintains a set of state variables that are used to control the Nmf and clustering routines.  Once set, the state variables maintain their values until changed by an API function.  For instance, one state variable represents the matrix to be factored (or used for clustering).  The API provides a function to load this matrix; once loaded, it can be repeatedly factored without the need for reloading.  The state variables and their default values are documented below.

All computations with the smallk library are performed in double precision.

Enumerations
------------

The SmallK API provides two enumerated types, one for the supported NMF algorithms and one for the clustering file output format.  These are::

	enum Algorithm
	{
		MU,      // Multiplicative Updating, Lee & Seung
		BPP,     // Block Principal Pivoting, Kim and Park
		HALS,    // Hierarchical Alternating Least Squares, Cichocki & Pan
		RANK2    // Rank2, Kuang and Park
	};

The default NMF algorithm is BPP.  The Rank2 algorithm is optimized for two-column or two-row matrices and is the underlying factorization routine for the clustering code.

:: 

	enum OutputFormat
	{
		XML,  // Extensible Markup Language
		JSON  // JavaScript Object Notation
	};

API functions
-------------

Initialization and cleanup
^^^^^^^^^^^^^^^^^^^^^^^^^^
:: 

	void Initialize(int& argc,     // in
		char**& argv)  // in

Call this function first, before all others in the API; initializes Elemental and the smallk library.

::

	bool IsInitialized()
    
Returns true if the library has been initialized via a call to Initialize(), false otherwise.

Call this function last, after all others in the API; performs cleanup for Elemental and the smallk library::

	void Finalize()

Versioning
^^^^^^^^^^
:: 

	unsigned int GetMajorVersion()

Returns the major release version number of the library as an unsigned integer.
:: 

	unsigned int GetMinorVersion()

Returns the minor release version number of the library as an unsigned integer.
:: 

	unsigned int GetPatchLevel()

Returns the patch version number of the library as an unsigned integer.
:: 

	std::string GetVersionString()

Returns the version of the library as a string, formatted as major.minor.patch.

Common functions
^^^^^^^^^^^^^^^^
:: 

	unsigned int GetOutputPrecision()

Returns the floating point precision with which numerical output will be written (i.e., the computed W and H matrix factors from the Nmf routine).  The default precision is six digits. 
:: 

	void SetOutputPrecision(const unsigned int num_digits)

Sets the floating point precision with which numerical output will be written.  Input values should be within the range [1, precision(double)].  Any inputs outside of this range will be adjusted. 
:: 

	unsigned int GetMaxIter()

Returns the maximum number of iterations allowed for NMF computations.  The default value is 5000.
:: 

	void SetMaxIter(const unsigned int max_iterations = 5000)

Sets the maximum number of iterations allowed for NMF computations.  The default of 5000 should be more than sufficient for most computations. 
:: 

	unsigned int GetMinIter()

Returns the minimum number of NMF iterations. The default value is 5.
:: 

	void SetMinIter(const unsigned int min_iterations = 5)

Sets the minimum number of NMF iterations to perform before checking for convergence. The convergence and progress estimation routines are non-trivial calculations, so increasing this value may result in faster performance. 
:: 

	unsigned int GetMaxThreads()

Returns the maximum number of threads used for NMF or clustering computations. The default value is hardware-dependent, but is generally the maximum number allowed by the hardware.
:: 

	void SetMaxThreads(const unsigned int max_threads);

Sets an upper limit to the number of threads used for NMF and clustering computations.  Inputs that exceed the capabilities of the hardware will be adjusted. This function is provided for scaling and performance studies.  
:: 

	void Reset()

Resets all state variables to their default values. 
:: 

	void SeedRNG(const int seed)

Seeds the random number generator (RNG) within the smallk library. Normally this RNG is seeded from the system time whenever the library is initialized.  The RNG is the ‘19937’ Mersenne Twister implementation provided by the C++ standard library.
:: 

	void LoadMatrix(const std::string& filepath)

Loads a matrix contained in the given file.  The file must either be a comma-separated value (.CSV) file for a dense matrix, or a MatrixMarket-format file (.MTX) for a sparse matrix. If the matrix cannot be loaded the library throws a std::runtime_error exception.
:: 

	bool IsMatrixLoaded()

Returns true if a matrix is currently loaded, false if not.
::

	std::string GetOuputDir()

Returns a string indicating the directory into which output files will be written.  The default is the current directory.
::

	void SetOutputDir(const std::string& outdir)

Sets the directory into which output files should be written. The ‘outdir’ argument can either be an absolute or relative path.  The default is the current directory.

NMF functions
^^^^^^^^^^^^^
:: 

	void Nmf(const unsigned int k, 
		const Algorithm algorithm     = Algorithm::BPP,
		const std::string& initfile_w = std::string(“”),
		const std::string& initfile_h = std::string(“”))

This function factors the input matrix A of nonnegative elements into nonnegative factors such that: A &cong; WH.  If a matrix is not currently loaded a std::logic_error exception will be thrown.  The default algorithm is NMF-BPP; provide one of the enumerated algorithm values to use a different algorithm.

Where A is mxn, W is mxk, and H is kxn.  The value of k a user defined argument, e.g., for clustering applications, k is the number of clusters.

Optional initializer matrices can be provided for the W and H factors via the ‘initfile_w’ and ‘initfile_h’ arguments. These files must contain fully dense matrices in .CSV format.  The W matrix initializer must have dimension mxk, and the H matrix initializer must have dimension kxn. If the initializer matrices do not match these dimensions exactly a std::logic_error exception is thrown.  If initializers are not provided, matrices W and H will be randomly initialized.

The computed factors W and H will be written to the output directory in the files ‘w.csv’ and ‘h.csv’.
    
Exceptions will be thrown (either from Elemental or smallk) in case of error.
:: 

	const double* LockedBufferW(unsigned int& ldim, unsigned int& height, unsigned int& width)

This function returns a READONLY pointer to the buffer containing the W factor computed by the Nmf routine, along with buffer and matrix dimensions.  The ‘ldim’, ‘height’, and ‘width’ arguments are all out parameters.  The buffer has a height of ‘ldim’ and a width of ‘width’.  The matrix W has the same width but a height of ‘height’, which may differ from ldim.  The W matrix is stored in the buffer in column-major order.  See the examples/smallk_example.cpp file for an illustration of how to use this function. 
:: 

	const double* LockedBufferH(unsigned int& ldim, unsigned int& height, unsigned int& width)

Same as LockedBufferW, but for the H matrix.
:: 

	double GetNmfTolerance()

Returns the tolerance value used to determine NMF convergence. The default value is 0.005. 
:: 

	void SetNmfTolerance(const double tol=0.005)

Sets the tolerance value used to determine NMF convergence.  The NMF algorithms are iterative, and at each iteration a progress metric is computed and compared with the tolerance value.  When the metric falls below the tolerance value the iterations stop and convergence is declared.  The tolerance value should satisfy 0.0 < tolerance < 1.0.  Any inputs outside this range will cause a `std::logic_error` exception to be thrown.
Clustering Functions
:: 

	void LoadDictionary(const std::string& filepath)

Loads the dictionary used for clustering. The dictionary is an ASCII file of text strings as described in the preprocessor input files section below.  If the dictionary file cannot be loaded a `std::runtime_error` exception is thrown.

::

	unsigned int GetMaxTerms()

Returns the number of highest-probability dictionary terms to store per cluster. The default value is 5.
:: 

	void SetMaxTerms(const unsigned int max_terms = 5)

Sets the number of highest-probability dictionary terms to store per cluster.
:: 

	OutputFormat GetOutputFormat()

Returns a member of the OutputFormat enumerated type; this is the file format for the clustering results.  The default output format is JSON.
:: 

	void SetOutputFormat(const OutputFormat = OutputFormat::JSON)

Sets the output format for the clustering result file. The argument must be one of the values in the OutputFormat enumerated type.
:: 

	double GetHierNmf2Tolerance()

Returns the tolerance value used by the NMF-RANK2 algorithm for hierarchical clustering.  The default value is 1.0e-4.
:: 

	void SetHierNmf2Tolerance(const double tol=1.0e-4)

Sets the tolerance value used by the NMF-RANK2 algorithm for hierarchical clustering.  The tolerance value should satisfy 0.0 < tolerance < 1.0.  Any inputs outside this range will cause a `std::logic_error` exception to be thrown.
:: 

	void HierNmf2(const unsigned int num_clusters)

This function performs hierarchical clustering on the loaded matrix, generating the number of clusters specified by the ‘num_clusters’ argument.  For an overview of the hierarchical clustering process, see the description below for the hierclust command line application.

This function generates two output files in the output directory: `assignments_N.csv` and `tree_N.{json, xml}`.  Here N is the number of clusters specified as an argument, and the tree file can be in either JSON XML format.

The content of the files is described below in the section on the hierclust command line application.
:: 

	void HierNmf2WithFlat(const unsigned int num_clusters)

This function performs hierarchical clustering on the loaded matrix, exactly as described for HierNmf2. In addition, it also computes a flat clustering result.  Thus four output files are generated.  The flat clustering result files are ‘assignments_flat_N.csv’ and ‘clusters_N.{json, xml}’.  The cluster file contents are documented below in the section on the flatclust command line application.

