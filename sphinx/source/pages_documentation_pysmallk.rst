########
Pysmallk
########

.. contents:: :local:

************
Introduction
************

Why Python? Although it’s perfectly fine to run SmallK from the command line, Python provides a great deal more flexibility that augments the C++ code with other tasks that are much more easily accomplished with a very high level language. Python distributions can be easily extended with open source libraries from third party sources as well, two examples being numpy and scipy, well-known standards for scientific computing in the Python community. There are numerous packages available that extend these scientific libraries into the data analytics domain as well, such as `scikit-learn <http://scikit-learn.org/stable/index.html>`_.

For using scientific Python, we strongly recommend the Anaconda Python distribution provided by `Continuum Analytics <http://continuum.io/>`_. Download and installation instructions for all platforms can be found `here <https://store.continuum.io/cshop/anaconda/>`_. Anaconda includes many if not most of the commonly used scientific and data analytics packages available and a very easy to use package manager and updating system. After installing Anaconda there will be available at the command line both a standard Python interpreter (type ‘python’) and an iPython interpreter (type ‘ipython’). We recommend using the iPython interpreter. In addition to the command line interfaces to Python, Anaconda includes the Spyder visual development environment featuring a very well thought out interface that makes developing Python code almost “too easy”. Spyder has many features found in the Matlab™ editor and a similar look and feel.

Anaconda also includes the Cython package, which is used by SmallK to integrate the Python and C++ code. `Cython <http://cython.org/>`_ includes support for most of the C++ standard and supports the latest GNU C++ compilers. Most if not all the standard libraries are supported and the latest version (20.2) has support for the standard template library (STL) as well.

**************************
Examples of Pysmallk Usage
**************************

Pysmallk has five classes, each of which represents one of the SmallK tools: SmallkAPI (the simplistic Smallk API), Flatclust, Hierclust, Matrixgen, and Preprocessor. These tools can be strung together into various kind of applications. Examples of such applications can be found in `examples/pysmallk_example.py` and in the `pysmallk/tests/` subdirectory.

The smallk_data repository contains several files (articles_matrix.mtx, articles_documents.txt, articles_dictionary.txt) that contain the matrix and associated text files created from 2,424 news articles. 

First, we will need to import numpy and the shared libary::

	import numpy as np
	import pysmallk
	
We then should apply the preprocessor to our data::

	p = pysmallk.Preprocessor()
	p.load_matrix(filepath='smallk_data/articles_matrix.mtx')
	p.load_dictionary(filepath='smallk_data/articles_dictionary.txt')
	p.load_documents(filepath='smallk_data/articles_documents.txt')
	
We will begin with the default inputs and run preprocess::
	
	p.preprocess()
	
Instead of writing the results to files, we can get the outputs from the Preprocessor class and pass them directly as inputs to the SmallkAPI class.::

	reduced_docs = p.get_reduced_documents()
	reduced_dict = p.get_reduced_dictionary()
	reduced_scores = p.get_reduced_scores()
	reduced_row_indices = p.get_reduced_row_indices()
	reduced_col_offsets = p.get_reduced_col_offsets()
	reduced_height = len(reduced_dict)
	reduced_width = len(reduced_docs)

Now let's instantiate the SmallkAPI object that we will use to do further computations.::

	sk = pysmallk.SmallkAPI()

One of the options for matrix loading is to pass in the appropriate fields for a sparse matrix, as so::

	sk.load_matrix(buffer=reduced_scores, row_indices=reduced_row_indices, col_offsets=reduced_col_offsets, height=reduced_height, width=reduced_width, nz=len(reduced_scores))

The input matrix alone is sufficient to run NMF and compute the factor matricies.
::

	sk.nmf(5, 'BPP')

This will compute the W and H factor matrices and subsequently write them to the files w.csv and h.csv, respectively.

We can continue with further calcuations using the same input matrix. For example, we can extract topic models from the input matrix if we also provide a dictionary.
::

	sk.load_dictionary(dictionary=reduced_dict)
	sk.hiernmf2(5)

This will use Hierarchical NMF to determine the final leaf nodes to use for the topic models and will output assignments_5.csv (cluster labels) and tree_5.xml.

Now let's say we want to create our own random matrix and pass that as a numpy matrix into SmallK.
::	

	a = np.random.random((256, 256))

In order to run the Hierclust or Flatclust applications, we will need to provide a dictionary file from which to select the top terms.
::

	pathtodict = args.indir + 'reuters_dictionary.txt'
	with open(pathtodict) as dictionary:
    	terms = dictionary.read().split("\n")
	    
For illustration, let's use the Flatclust object and extract the resulting assignments from running NMF.
::

	f = pysmallk.Flatclust()

	f.load_matrix(matrix=a)
	f.load_dictionary(dictionary=terms)
	f.cluster(16, algorithm='HALS')
	a = f.get_assignments()

Now the variable 'a' holds a list of the computed assignment labels for each of the 256 elements in our original matrix.

When we are finished, we should clean up the environment before exiting::

	sk.finalize()
	f.finalize()


******************
Pysmallk Functions
******************

Pysmallk has five classes, each of which represents one of the SmallK tools: SmallkAPI (the simplistic Smallk API), Flatclust, Hierclust, Matrixgen, and Preprocessor. Each of these classes can be imported as follows::

	from pysmallk import SmallkAPI
	from pysmallk import Flatclust
	from pysmallk import Hierclust
	from pysmallk import Matrixgen
	from pysmallk import Preprocessor

Each class’s primary functions are documented in the sections below. The parameters are either marked [in] or [kwarg] which represent, respectively, positional and keyword arguments.

Preprocessor
============
:: 

	def parser()

Returns the parsed arguments for the default command line application. The command line arguments are the same as those for the C++ binary application preprocessor.
::

	def load_matrix(filepath=””, height=0, width=0, nz=0, buffer=[], row_indices=[], col_offsets=[])

Load an input matrix.

1. To load a matrix from a file::

    * filepath:      The path to the input matrix

2. To load a sparse matrix from Matrixgen::

    * height:        The height of the sparse matrix
    * width:         The width of the sparse matrix
    * sparse_matrix: The sparse matrix returned from Matrixgen

3. To load a sparse matrix from python::

    * height:        The height of the sparse matrix
    * width:         The width of the sparse matrix
    * nz:            The number of non-zeros in the sparse matrix
    * buffer:        List of doubles containing the non-zero elements of the sparse matrix
    * row_indices:   List of integers representing the row indices of the sparse matrix
    * col_offsets:   List of integers representing the column offsets of the sparse matrix

::

	def load_dictionary(filepath=None, dictionary=None)

Loads a dictionary from either a filepath or a list of dictionary strings.
::

	def load_documents(filepath=None, documents=None)

Loads a documents from either a filepath or a list of document strings.
::

	def get_reduced_documents()

Returns the reduced documents.
::

	def get_reduced_dictionary()

Returns the reduced dictionary.
::

	def get_reduced_scores()

Returns the non-zero scores from the reduced matrix.
::

	def get_reduced_row_indices ()
	
Returns the row indices for the reduced matrix.
::

	def get_reduced_col_offsets ()

Returns the column offsets for the reduced matrix.
::

	def get_reduced_field (filepath="", values=[])
	
Loads a field from either a filepath or a list of field strings. Returns the reduced fields.
::

	def preprocess(maxiter=1000, docsperterm=3,termsperdoc=5, boolean_mode=0)

Preprocesses the matrix.
    
* maxiter:      The maximum number of iterations (optional)
* docsperterm:  The number of documents required per term (optional)
* termsperdoc:  The number of terms requried per document (optional)
* boolean_mode: All nonzero matrix elements will be treated as if they had the value 1.0  (optional)

::

	def write_output(matrix_filepath, dict_filepath, docs_filepath, precision=4)

Writes the preprocessor results to files.

* matrix_filepath:     The filepath for writing the matrix
* dict_filepath:       The filepath for writing the dictionary
* docs_filepath:       The filepath for the documents
* precision:           The precision with which to write the outputs (optional)

Matrixgen
=========
::

	def parser()

Returns the parsed arguments for the default command line application. The command line arguments are the same as those for the C++ binary application matrixgen.
::

	def uniform(m, n, center=0.5, radius=0.5)

Generates a uniform matrix. Returns a tuple of the list of values, the height, and the width.

* m:       The desired height
* n:       The desired width
* center:  Center with which to initialize the RNG 
* radius:  Radius with which to initialize the RNG 

::

	def densediag(m, n, center=0.5, radius=0.5)

Generates a dense diagonal matrix. Returns a tuple of the list of values, the height, and the width.

* m:       The desired height
* n:       The desired width
* center:  Center with which to initialize the RNG 
* radius:  Radius with which to initialize the RNG 

::

	def identify(m, n)

Generates an identify matrix. Returns a tuple of the list of values, the height, and the width.

* m:       The desired height
* n:       The desired width

::

	def sparsediag(n, center=0.5, radius=0.5)

Generates a sparse diagonal matrix. Returns a tuple of the list of values, the height, and the width.

* n:       The desired width
* center:  Center with which to initialize the RNG 
* radius:  Radius with which to initialize the RNG 

::

	def ones(m, n)

Generates a matrix of ones. Returns a tuple of the list of values, the height, and the width.

* m:       The desired height
* n:       The desired width

::

	def zeros(m, n)

Generates a matrix of zeros. Returns a tuple of the list of values, the height, and the width.

* m:       The desired height
* n:       The desired width

::

	def sparse(m, n, nz)

Generates a random sparse matrix. Returns a tuple of the list of values, the height, and the width.

* m:       The desired height
* n:       The desired width
* nz:      The number of non zeros in the matrix

::

	def write_output(filename, precision=6)

Writes the generated matrix to file.

* filename:     The filepath for writing the matrix
* precision:    The precision with which to write the matrix

SmallkAPI
=========
::

	def parser()

Returns the parsed arguments for the default command line application. The dictionary containing the parsed arguments.
::

	def get_major_version()

Returns the major version of SmallK.
::

	def get_minor_version()

Returns the minor version of SmallK.
::

	def get_patch_level()

Returns the patch level of SmallK.
::

	def get_version_string()

Returns a string representation of the version of SmallK.
::

	def load_matrix(filepath="", height=0, width=0, delim="", buffer=[], matrix=[], 
        nz=0, row_indices=[], col_offsets=[], column_major=False, sparse_matrix=None):


Load an input matrix.

1. To load a matrix from a file:

	* filepath:      The path to the input matrix

2. To load a sparse matrix from python:

	* height:        The height of the sparse matrix
	* width:         The width of the sparse matrix
	* nz:            The number of non-zeros in the sparse matrix
	* buffer:        List of doubles containing the non-zero elements of the sparse matrix
	* row_indices:   List of integers representing the row indices of the sparse matrix
	* col_offsets:   List of integers representing the column offsets of the sparse matrix

3. To load a dense matrix from python:

	* height: The height of the dense matrix	
	* width:         The width of the dense matrix
	* buffer: List of doubles containing the elements of the dense matrix

4. To load a numpy matrix from python:

	* matrix:        The numpy matrix
	* column_major:  Boolean for whether or not the matrix is column major (optional)

**Note: Internal to SmallK, the matrix is stored in column-major order. When you are loading a numpy matrix, the assumption is that your matrix is in row-major order. If this is not the case, you can pass column_major=True in as a keyword argument. When directly loading a dense matrix, the assumption is that your buffer holds the data in column-major order as well.**

:: 

	def is_matrix_loaded()

Indicates whether or not a matrix has been loaded.

:: 

	def nmf(k, algorithm, infile_W="", infile_H="", precision=4, min_iter=5, max_iter=5000, tol=0.005, max_threads=8, outdir=".")

Runs NMF on the loaded matrix using the supplied algorithm and implementation details.

*    k:           The desired number of clusters
*    algorithm:   The desired NMF algorithm
* infile_W:    Initialization for W (optional)
* infile_H:    Initialization for H (optional)
* precision:   Precision for calcuations (optional)
* min_iter:    Minimum number of iterations (optional)
* max_iter:    Maximum number of iterations (optional)
* tol:         Tolerance for determing convergence (optional)
* max_threads: Maximum number of threads to use (optional)
* outdir:      Output directory for files (optional)

::

	def get_inputs()

Returns a dictionary of the supplied inputs to the nmf function.
::

	def get_H()

Returns the output H matrix.
::

	def get_W()

Returns the output W matrix.
::

	def load_dictionary (filepath="", dictionary=[])

Loads a dictionary from either a filepath or a list of dictionary strings.
::

	def hiernmf2(k, format="XML", maxterms=5, tol=0.0001)

Runs HierNMF2 on the loaded matrix.

*    k:           The desired number of clusters
* format:      Output format, XML or JSON (optional)
* maxterms:    Maximum number of terms (optional)
* tol:         Tolerance to use for determining convergence (optional)

::

	def finalize()

Cleans up the elemental and smallk environment.

Flatclust
=========
::

	def parser()

Returns the parsed arguments for the default command line application. The command line arguemnts are the same as those for the C++ binary application flatclust.
:: 

	def load_matrix(**kwargs)

Load an input matrix.

1. To load a matrix from a file:

	* filepath:      The path to the input matrix

2. To load a sparse matrix from python:

	* height:        The height of the sparse matrix
	* width:         The width of the sparse matrix
	* nz:            The number of non-zeros in the sparse matrix
	* buffer:        List of doubles containing the non-zero elements of the sparse matrix
	* row_indices:   List of integers representing the row indices of the sparse matrix
	* col_offsets:   List of integers representing the column offsets of the sparse matrix

3. To load a sparse matrix from Matrixgen:

  	* height:        The height of the sparse matrix
  	* width:         The width of the sparse matrix
  	* sparse_matrix: The sparse matrix returned from Matrixgen

4. To load a dense matrix from python:

	* height: The height of the dense matrix	
	* width:         The width of the dense matrix
	* buffer: List of doubles containing the elements of the dense matrix

5. To load a numpy matrix from python:

	* matrix:        The numpy matrix
	* column_major:  Boolean for whether or not the matrix is column major (optional)

**Note: Internal to SmallK, the matrix is stored in column-major order. When you are loading a numpy matrix, the assumption is that your matrix is in row-major order. If this is not the case, you can pass column_major=True in as a keyword argument. When directly loading a dense matrix, the assumption is that your buffer holds the data in column-major order as well.**

::

	def load_dictionary (filepath="", dictionary=[])

Loads a dictionary from either a filepath or a list of dictionary strings.

::

	def cluster(k, infile_W='', infile_H='', algorithm="BPP", maxterms=5, verbose=True, min_iter=5, max_iter=5000, max_threads=8, tol=0.0001)

Runs NMF on the loaded matrix using the supplied algorithm and implementation details.

* k:           The desired number of clusters
* infile_W:    Initialization for W (optional)
* infile_H:    Initialization for H (optional)
* algorithm:   The desired NMF algorithm (optional)
* maxterms:    Maximum number of terms per cluster (optional)
* verbose:     Boolean for whether or not to be verbose (optional)
* min_iter:    Minimum number of iterations (optional)
* max_iter:    Maximum number of iterations (optional)
* max_threads: Maximum number of threads to use (optional)
* tol:         Tolerance for determing convergence (optional)

::

	def get_top_indices()

Return the top term indices for each cluster. The length of the returned array is maxterms*k, with the first maxterms elements belonging to the first cluster, the second maxterms elements belonging to the second cluster, etc.
::

	def get_top_terms()

Return the top terms for each cluster.The length of the returned array is maxterms*k, with the first maxterms elements belonging to the first cluster, the second maxterms elements belonging to the second cluster, etc.
::

	def get_assignments()

Return the list of cluster assignments for each document.
::

	def write_output(assignfile, treefile, outdir='./', format='XML')

Writes the flatclust results to files.

*    assignfile:     The filepath for writing assignments
*    fuzzyfile:      The filepath for writing fuzzy assignments
*    treefile:       The filepath for the tree results
*    outdir:         The output directory for the output files (optional)
*    format:         The output format JSON or XML (optional)

::

	def finalize()

Cleans up the elemental and smallk environment.

Heirclust
=========
::

	def parser()

Returns the parsed arguments for the default command line application. The command line arguemnts are the same as those for the C++ binary application hierclust.

::

	def load_matrix(**kwargs)

Load an input matrix.

1. To load a matrix from a file:

	* filepath:      The path to the input matrix

2. To load a sparse matrix from python:

	* height:        The height of the sparse matrix
	* width:         The width of the sparse matrix
	* nz:            The number of non-zeros in the sparse matrix
	* buffer:        List of doubles containing the non-zero elements of the sparse matrix
	* row_indices:   List of integers representing the row indices of the sparse matrix
	* col_offsets:   List of integers representing the column offsets of the sparse matrix

3. To load a sparse matrix from Matrixgen:
  	* height:        The height of the sparse matrix
  	* width:         The width of the sparse matrix
  	* sparse_matrix: The sparse matrix returned from Matrixgen

4. To load a dense matrix from python:
	* height: The height of the dense matrix	
	* width:         The width of the dense matrix
	* buffer: List of doubles containing the elements of the dense matrix

5. To load a numpy matrix from python:
	* matrix:        The numpy matrix
	* column_major:  Boolean for whether or not the matrix is column major (optional)

**Note: Internal to SmallK, the matrix is stored in column-major order. When you are loading a numpy matrix, the assumption is that your matrix is in row-major order. If this is not the case, you can pass column_major=True in as a keyword argument. When directly loading a dense matrix, the assumption is that your buffer holds the data in column-major order as well.**

::

	def load_dictionary (filepath="", dictionary=[])

Loads a dictionary from either a filepath or a list of dictionary strings.

::

	def cluster(k, infile_W='', infile_H='', maxterms=5, unbalanced=0.1, trial_allowance=3,  verbose=True, flat=0, min_iter=5, max_iter=5000, max_threads=8, tol=0.0001)

Runs NMF on the loaded matrix using the supplied algorithm and implementation details.

* k:           The desired number of clusters
* infile_W:    Initialization for W (optional)
* infile_H:    Initialization for H (optional)
* maxterms:    Maximum number of terms per cluster (optional)
* unbalanced:      Unbalanced parameter (optional)
* trial_allowance: Number of trials to use (optional)
* verbose:     Boolean for whether or not to be verbose (optional)
* flat:            Whether or not to flatten the results (optional)
* min_iter:    Minimum number of iterations (optional)
* max_iter:    Maximum number of iterations (optional)
* max_threads: Maximum number of threads to use (optional)
* tol:         Tolerance for determing convergence (optional)

::

	def get_top_indices()

Return the top term indices for each cluster. The length of the returned array is maxterms*k, with the first maxterms elements belonging to the first cluster, the second maxterms elements belonging to the second cluster, etc.

::

	def get_assignments()

Return the list of cluster assignments for each document.
::

	def write_output(assignfile, fuzzyfile, treefile, outdir='./', format='XML')

Writes the flatclust results to files.

*    assignfile:     The filepath for writing assignments
*    fuzzyfile:      The filepath for writing fuzzy assignments
*    treefile:       The filepath for the tree results
*    outdir:         The output directory for the output files (optional)
*    format:         The output format JSON or XML (optional)

::

	def finalize()

Cleans up the elemental and smallk environment.


