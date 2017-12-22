##################
Command Line Tools
##################

.. toctree::
   :maxdepth: 8

************
Introduction
************

The SmallK library provides a number of algorithm implementations for low rank approximation of a matrix. These can be used for performing various data analytics tasks such as topic modeling, clustering, and dimension reduction. This section will provide more in-depth description of the tools available with examples that can be expanded/modified for other application domains.

Before diving into the various tools, it will be helpful to set up the command line environment to easily run the various executables that comprise the SmallK library. First the command line needs to know where to find the executable files to run the tools. Since while installing SmallK  ``make_install`` was run, the executables are located in ``/usr/local/smallk/bin``. Thus, this should be added to the ``$PATH`` system variable or added to the environment. The following command line performs the task of modifying the path avoiding the need to cd into directories were the tools are located:

.. code-block:: bash

   export PATH=/usr/local/smallk/bin:$PATH

This allows the tools to be executed from any directory.

A subset of these tools are also available from the pysmallk library: smallkapi (mirrors the NMF command line application), matrixgen, preprocessor, flatclust, and hierclust. The command line arguments are the same as those documented below. These tools are available within the ``/pysmallk/tests/ directory`` and can be executed as follows:

.. code-block:: none

   [python binary] [tool].py [command line arguments]

For example:

.. code-block:: bash

   python preprocessor.py --indir smallk_data

************
Preprocessor
************

Overview
========

The preprocessor prunes rows and columns from term-frequency matrices, attempting to generate a result matrix that is more suitable for clustering.  It also computes tf-idf weights for the remaining entries.  Therefore the input matrix consists of nonnegative integers, and the output matrix consists of floating point numbers between 0.0 and 1.0.  The MatrixMarket file format is used for the input and output matrices.

Rows (terms) are pruned if a given term appears in fewer than ``DOCS_PER_TERM`` documents.  The value of ``DOCS_PER_TERM`` is a command-line parameter with a default value of 3.  For a term-frequency
input matrix, in which the matrix elements represent occurrence counts for the terms, this parameter actually specifies the minimum row sum for each term.  Any rows whose rowsums are less than this value will be pruned.

Columns (docs) are pruned if a given document contains fewer than ``TERMS_PER_DOC`` terms.  The value of ``TERMS_PER_DOC`` is a command-line parameter with a default value of 5.

Whenever columns (documents) are pruned the preprocessor checks the remaining columns for uniqueness.  Any duplicate columns are identified and a representative column is chosen as the survivor. The code always selects the column with the largest column index in such groups as the survivor. The preprocessor continues to prune rows and columns until it finds no further candidates for pruning. It then computes new tf-idf scores for the resulting entries and writes out the result matrix in MatrixMarket format.

If the preprocessor should prune all rows or columns, it writes an error message to the screen and terminates without generating any output.

Input Files
===========

The preprocessor requires three input files: a matrix file, a dictionary file, and a document file.  The matrix file contains a sparse matrix in MatrixMarket format (.mtx).  This is a term-frequency matrix, and all entries should be positive integers. The preprocessor can also read in matrices containing floating-point inputs, but only if ``boolean mode`` is enabled; this will be described below. The preprocessor does not support dense matrices, since the typical matrices encountered in topic modeling problems are extremely sparse, with occupancies generally less than 1%.

The second file required by the preprocessor is a ``dictionary file``.  This is a simple ASCII text file containing one entry per line.  Entries represent keywords, bigrams, or other general text strings the user is interested in.  Each line of the file is treated as a ``keyword``, so multi-word keywords are supported as well.  The smallk/data folder contains a sample dictionary file called ``dictionary.txt``.  The first few entries are:

.. code-block:: none

   triumph
   dey
   canada
   finger
   circuit
   ...

The third file required by the preprocessor is a ``documents file``.  This is another simple ASCII text file containing one entry per line.  Entries represent document names or other unique identifiers.  The smallk/data folder also contains a sample documents file called ``documents.txt``.  The first few entries of this file are:

.. code-block:: none

   52828-11101.txt
   51820-10202.txt
   104595-959.txt
   60259-3040.txt
   ...

These are the unique document identifiers for the user who generated the file.  Your identifiers will likely have a different format.

Finally, the preprocessor requires these files to have the following names: matrix.mtx, dictionary.txt, and documents.txt.  The input folder containing these files can be specified on the command line (described below).  The output of the preprocessor is a new set of files called ``reduced_matrix.mtx``, ``reduced_dictionary.txt``, and ``reduced_documents.txt``.

Command Line Options
====================

The preprocessor binary is called ``preprocess_tf``, to emphasize the fact that it operates on term-frequency matrices.  If the binary is run with no arguments, it prints out the following information:

.. code-block:: none

   preprocess_tf
   	--indir <path>
   	[--outdir (defaults to current directory)]
   	[--docs_per_term 3]
   	[--terms_per_doc 5]
   	[--maxiter 1000]
   	[--precision 4]
   	[--boolean_mode 0]

Only the first parameter, ``--indir``, is required.  All remaining params are optional and have the default values indicated.

The meanings of the various options are as follows:

	1. ``--indir``: path to the folder containing the files ``matrix.mtx``, ``dictionary.txt``, and ``documents.txt``; may be in small_data for example
	2. ``--outdir``: path to the folder to into which results should be written
	3. ``--docs_per_term``: any rows whose entries sum to less than this value will be pruned
	4. ``--terms_per_doc``: any columns whose entries sum to less than this value will be pruned
	5. ``--maxiter``: perform no more than this many iterations
	6. ``--precision``: the number of digits of precision with which to write the output matrix
	7. ``--boolean_mode``:  all nonzero matrix elements will be treated as if they had the value 1.0.  In other words, the preprocessor will ignore the actual frequency counts and treat all nonzero entries as if they were 1.0.

Sample Runs
===========

Here is a sample run of the preprocessor using the data provided in the smallk distribution.  This run was performed from the top-level smallk folder after building the code:

.. code-block:: none

   preprocessor/bin/preprocess_tf --indir data
   
     Command line options: 
   
                 indir: data/
                outdir: current directory
         docs_per_term: 3
         terms_per_doc: 5
              max_iter: 1000
             precision: 4
          boolean_mode: 0

	Loading input matrix data/matrix.mtx
    		Input file load time: 1.176s.

	Starting iterations...
    		[1] height: 39771, width: 11237, nonzeros: 877453
	Iterations finished.
	New height: 39727
	New width: 11237
	New nonzero count: 877374
	Processing time: 0.074s.

	Writing output matrix reduced_matrix.mtx
	Output file write time: 2.424s.
	Writing dictionary file reduced_dictionary.txt
	Writing documents file reduced_documents.txt
	Dictionary + documents write time: 0.08s.

*********
Matrixgen
*********

Overview
========

The matrix generator application is a simple tool for generating simple matrices.  These matrices can be loaded by the NMF and clustering tools for various testing scenarios.  Use of the matrix generator is entirely optional.

Command Line Options
====================

Running the matrixgen binary with no options generates the following output:

.. code-block:: none

	matrixgen 

	Usage: matrixgen
         	--height <number of rows> 
         	--width  <number of cols> 
         	--filename <path> 
        	[--type  UNIFORM]  UNIFORM:     matrix with uniformly-distributed random entries
                           DENSE_DIAG:  dense diagonal matrix with uniform random entries
                           SPARSE_DIAG: sparse diagonal matrix with uniform random entries
                           IDENTITY:    identity matrix
                           ONES:        matrix of all ones
                           ZEROS:       matrix of all zeros
                           SPARSE:      sparse matrix with uniform random entries
                                        specify 'nz_per_col' to control occupancy

        	[--rng_center  0.5]   center of random numbers
        	[--rng_radius  0.5]   radius of random numbers
        	[--precision   6]     digits of precision
        	[--nz_per_col  1]     (SPARSE only) nonzeros per column

The ``--height``, ``--width``, and ``--filename`` options are required.  All others are optional and have the default values indicated.

The meanings of the various options are as follows:

	1. ``--height``: number of rows in the generated matrix
	2. ``--width``: number of columns in the generated matrix
	3. ``--filename``: name of the output file
	4. ``--type``: the type of matrix to be generated; the default is a uniformly-distributed random matrix
	5. ``--rng_center``: random number distribution will be centered on this value
	6. ``--rng_radius``: random numbers will span this distance to either side of the center value
	7. ``--precision``: the number of digits of precision with which to write the output matrix
	8. ``--nz_per_col``:  number of nonzero entries per sparse matrix column; valid only for SPARSE type

Sample Runs
===========

Suppose we want to generate a matrix of uniformly-distributed random numbers.  The matrix should have a height of 100 and a width of 16, and should be written to a file called ``w_init.csv``.  Use the matrix generator as follows:

.. code-block:: none

   matrixgen --height 100 --width 16 --filename w_init.csv

**************************************
Nonnegative Matrix Factorization (NMF)
**************************************

Overview
========

The NMF command line application performs nonnegative matrix factorization on dense or sparse matrices. If the input matrix is denoted by A, nonnegative matrix factors Wand H are computed such that :math:`\matr{A} \cong \matr{W} \matr{H}`.

Matrix :math:`\matr{A}` can be either dense or sparse; matrices :math:`\matr{W}` and :math:`\matr{H}` are always dense. Matrix :math:`\matr{A}` has m rows and n columns; matrix :math:`\matr{W}` has m rows and k columns; matrix :math:`\matr{H}` has k rows and n columns. Parameter k is a positive integer and is typically much less than either m or n.

Command Line Options
====================

Running the nmf application with no command line parameters will cause the application to display all params that it supports. These are:

.. code-block:: none

	Usage: nmf
        --matrixfile <filename>  Filename of the matrix to be factored.
                                 Either CSV format for dense or MatrixMarket format for sparse.
        --k <integer value>      Inner dimension for factors W and H.
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
        [--maxthreads    8]      Upper limit to thread count. 
        [--normalize  1]         Whether to normalize W and scale H.
                                     1 == yes, 0 == no 
        [--verbose  1]           Whether to print updates to the screen. 
                                     1 == print updates, 0 == silent

The --matrixfile and --k options are required; all others are optional and have the default values indicated.  The meanings of the various options are as follows:

	1.  ``--matrixfile``: Filename of the matrix to be factored.  CSV files are supported for dense matrices and MTX files for sparse matrices.
	2.  ``--k``: the width of the W matrix (inner dimension of the matrix factors)
	3.  ``--algorithm``: identifier for the factorization algorithm
	4.  ``--stopping``: the method used to terminate the iterations; use PG_RATIO unless you have a specific reason not to
	5.  ``--tol``: tolerance value used to terminate iterations; when the progress metric falls below this value iterations will stop; typical values are in the 1.0e-3 or 1.0e-4 range
	6.  ``--tolcount``: a positive integer representing the number of successive iterations for which the progress metric must have a value <= tolerance; default is 1, which means the iterations will terminate on the first iteration with: ``progress_metric <= tolerance`` 
	7.  ``--infile_W``: CSV file containing the mxk initial values for matrix W; if omitted, W is randomly initialized
	8.  ``--infile_H``:  CSV file containing the kxn initial values for matrix H; if omitted, H is randomly initialized
	9.  ``--outfile_W``: filename for the computed W factor; default is w.csv
	10. ``--outfile_H``: filename for the computed H factor; default is h.csv
	11. ``--miniter``: the minimum number of iterations to perform before checking progress; for smaller tolerance values, you may want to increase this number to avoid needless progress checks
	12. ``--maxiter``: the maximum number of iterations to perform
	13. ``--outprecision``: matrices W and H will be written to disk using this many digits of precision
	14. ``--maxthreads``: the maximum number of threads to use; the default is to use as many threads as the hardware can support (your number may differ from that shown) 
	15. ``--normalize``: whether to normalize the columns of the W matrix and correspondingly scale the rows of H after convergence
	16. ``--verbose``: whether to display updates to the screen as the iterations progress 

Sample Runs
===========

The smallk distribution utilizes another repository `smallk_data <https://github.com/smallk/smallk_data>`_ (clone this repository from github) with a matrix file ``reuters.mtx``.  This is a tf-idf weighted matrix derived from the popular Reuters data set used in machine learning experiments.  

Suppose we want to factor the Reuters matrix using a k value of 8.  We would do that as follows, assuming that we are in the top-level smallk folder after building the code and that the ``smallk_data`` repository was cloned into `data`:

.. code-block:: none

		nmf/bin/nmf --matrixfile data/reuters.mtx  --k 8

Note that if ``make install`` was run during installation and the $PATH variable or environment variable was set as above, this could also be called with:

.. code-block:: none

		usr/local/bin/nmf --matrixfile data/reuters.mtx  --k 8

If we want to instead use the HALS algorithm with k=16, a tolerance of 1.0e-4, and also perform 10 iterations prior to checking progress, we would use this command line:

.. code-block:: none

		nmf/bin/nmf --matrixfile data/reuters.mtx --k 16 --algorithm HALS --tol 1.0e-4 --miniter 10

To repeat the previous experiment but with new names for the output files, we would do this:

.. code-block:: none

		nmf/bin/nmf --matrixfile data/reuters.mtx --k 16 --algorithm HALS --tol 1.0e-4 
			--miniter 10 --outfile_W w_hals.csv -outfile_H h_hals.csv

*********
Hierclust
*********

Overview
========

First, we briefly describe the algorithm and the references section provides pointers to papers with detailed descriptions of the algorithms.  NMF-RANK2 for hierarchical clustering generates a binary tree of items. We refer to a node in the binary tree and the items associated with the node interchangeably. This method begins by placing all data items in the root node. The number of leaf nodes to generate is specified (user input). The algorithm proceeds with the following steps, repeated until the maximum number of leaf nodes, ``max_leaf_nodes``, is reached:

	1.  Pick the leaf node with the highest score (at the very beginning where only a root node is present, just pick the root node)
	2.  Apply NMF-RANK2 to the node selected in step 1, and generate two new leaf nodes
	3.  Compute a score for each of the two leaf nodes generated in step 2
	4.  Repeat until the desired number of leaf nodes has been generated

Step 2 implements the details of the node splitting into child nodes. Outlier detection plays a crucial role in hierarchical clustering to generate a tree with well-balanced and meaningful clusters. To implement this, we have two additional parameters in step 2: *trial_allowance* and *unbalanced*.

The parameter *trial_allowance* is the number of times that the program will try to split a node into two meaningful clusters. In each trial, the program will check if one of the two generated leaf nodes is an outlier set. If the outlier set is detected, the program will delete the items in the outlier set from the node being split and continue to the next trial. If all the trials are finished and the program still cannot find two meaningful clusters for this node, all the deleted items are "recycled" and placed into this node again, and this node will be labeled as a "permanent leaf node" that cannot be picked in step 1 in later iterations.

The parameter *unbalanced* is a threshold parameter to determine whether two generated leaf nodes are unbalanced. Suppose two potential leaf nodes L and R are generated from the selected node and L has fewer items than R. Let us denote the number of items in a node N as :math:`\left | N \right |`. L and R are called *unbalanced* if 

:math:`\left | L \right | < unbalanced * ( \left | L \right | + \left | R \right | )`

Note that if L and R are unbalanced, the potential node L with fewer items is not necessarily regarded as an outlier set. Please see the referenced paper for more details `[3] <http://smallk.github.io/publications/>`_.

Internally, NMF-RANK2 is applied to each leaf node to compute the score in step 3. The computed result matrices W and H in step 3 are cached so that we can avoid duplicate work in step 2 in later iterations.

The score for each leaf node is based on a modified version of the NDCG (*Normalized Discounted Cumulative Gain*) measure, a common measure in the information retrieval community. A leaf node is associated with a "topic vector", and we can define "top terms" based on the topic vector. A leaf node will receive a high score if its top terms are a good combination of the top terms of its two potential children; otherwise it receives a low score.

The hierclust application generates two output files. One file contains the assignments of documents to clusters. This file contains one integer for each document (column) of the original matrix. The integers are the cluster labels for that cluster that the document was assigned to.  If the document could not be assigned to a cluster, a -1 will be entered into the file, indicating that the document is an outlier.

The other output file contains information for each node in the factorization binary tree.  The items in this file are:

	1.  ``id``: a unique id for this node
	2.  ``level``: the level in the tree at which this node appears; the root is at level 0, the children of the root are at level 1, etc.
	3.  ``label``: the cluster label for this node (meaningful only for leaf nodes)
	4.  ``parent_id``: the unique id of the parent of this node (the root node has parent_id == 0)
	5.  ``parent_label``: the cluster label of the parent of this node
	6.  ``left_child``: a Boolean value indicating whether this node is the left or right child of its parent 
	7.  ``left_child_label``: the cluster label of the left child of this node (leaf nodes have -1 for this value)
	8.  ``right_child_label``: the cluster label of the right child of this node (leaf nodes have -1 for this value)
	9.  ``doc_count``: the number of documents that this node represents
	10. ``top_terms``: the highest probability dictionary terms for this node

The node id values and the left or right child indicators can be used to unambiguously reconstruct the factorization tree. 

Command Line Options
====================

Running the hierclust application with no command line parameters will cause the application to display all params that it supports.  These are:

.. code-block:: none

	Usage: hierclust/bin/hierclust
        --matrixfile <filename>     Filename of the matrix to be factored.
                                    Either CSV format for dense or MatrixMarket format for sparse.
        --dictfile <filename>       The name of the dictionary file.
        --clusters <integer>        The number of clusters to generate.
        [--initdir  (empty)]        Directory of initializers for all Rank2 factorizations.
                                    If unspecified, random init will be used. 
        [--tol  0.0001]             Tolerance value for each factorization. 
        [--outdir  (empty)]         Output directory.  If unspecified, results will be 
                                    written to the current directory.
        [--miniter  5]              Minimum number of iterations to perform.
        [--maxiter  5000]           Maximum number of  iterations to perform. 
        [--maxterms  5]             Number of terms per node. 
        [--maxthreads    8]         Upper limit to thread count. 
        [--unbalanced  0.1]         Threshold for determining leaf node imbalance. 
        [--trial_allowance  3]      Number of split attempts. 
        [--flat  0]                 Whether to generate a flat clustering result. 
                                        1 == yes, 0 == no
        [--verbose  1]              Whether to print updates to the screen.
                                        1 == yes, 0 == no
        [--format  XML]             Format of the output file containing the tree.
                                        XML: XML format
                                        JSON: JavaScript Object Notation
        [--treefile  tree_N.ext]    Name of the output file containing the tree.
                                    N is the number of clusters for this run.
                                    The string 'ext' depends on the desired format.
                                    This filename is relative to the outdir.
        [--assignfile assignments_N.csv]  Name of the file containing final assignments.
                                          N is the number of clusters for this run.
                                          This filename is relative to the outdir.

The ``--matrixfile``, ``--dictfile``, and ``--clusters`` options are required; all others are optional and have the default values indicated.  The meanings of the various options are as follows:

	1.  ``--matrixfile``: Filename of the matrix to be factored.  CSV files are supported for dense matrices and MTX files for sparse matrices.
	2.  ``--dictfile``: absolute or relative path to the dictionary file
	3.  ``--clusters``: the number of leaf nodes (clusters) to generate
	4.  ``--initdir``:  Initializer matrices for W and H are loaded from the initdir directory. The matrices are assumed to have the names ``Winit_1.csv``, ``Hinit_1.csv``, ``Winit_2.csv``, ``Hinit_2.csv``, etc. It is up to the user to ensure that enough matrices are present in this dir to run the HierNMF2 code to completion.  The number of matrices used is non-deterministic, so trial-and-error may be required to find a lower bound on the matrix count. This feature is used for testing (such as comparisons with Matlab), in which each factorization problem has to proceed from a known initializer. The W initializer matrices must be of shape m x 2, and the H initializer matrices must be of shape 2 x n.
	5.  ``--tol``: tolerance value for each internal NMF-RANK2 factorization; the stopping criterion is the ratio of projected gradient method
	6.  ``--outdir``: path to the folder into which to write the output files; if omitted results will be written to the current directory
	7.  ``--miniter``: minimum number of iterations to perform before checking progress on each NMF-RANK2 factorization
	8.  ``--maxiter``: the maximum number of iterations to perform on each NMF-RANK2 factorization
	9. ``--maxterms``: the number of dictionary keywords to include in each node
	10. ``--maxthreads``: the maximum number of threads to use; the default is to use as many threads as the hardware can support (your number may differ from that shown) 
	11. ``--unbalanced``: threshold value for declaring leaf node imbalance (see explanation above)
	12. ``--trial_allowance``: maximum number of split attempts for any node (see explanation above)
	13. ``--flat``: whether to generate a flat clustering result in addition to the hierarchical clustering result
	14. ``--verbose``: whether to display updates to the screen as the iterations progress
	15. ``--format``: file format to use for the clustering results
	16. ``--treefile``: name of the output file for the factorization tree; uses the format specified by the format parameter
	17. ``--assignfile``: name of the output file for the cluster assignments

Sample Runs
===========

The smallk distribution has available a ``smallk_data`` repository on github with a matrix file ``reuters.mtx`` and an associated dictionary file ``reuters_dictionary.txt``.  These files are derived from the popular Reuters data set used in machine learning experiments.  

As above, it is assumed that the `smallk_data <https://github.com/smallk/smallk_data>`_ repository was cloned into ``data`` and that the commands can be run as below or from ``/usr/local/bin``.

Suppose we want to perform hierarchical clustering on this data set and generate 10 leaf nodes.  We would do that as follows, assuming that we are in the top-level smallk folder after building the code:

.. code-block:: none

   hierclust/bin/hierclust --matrixfile data/reuters.mtx  --dictfile data/reuters_dictionary.txt --clusters 10

This will generate two result files in the current directory: ``tree_10.xml`` and ``assignments_10.csv``.

If we want to instead generate 10 clusters, each with 8 terms, using JSON output format, we would use this command line:

.. code-block:: none

   hierclust/bin/hierclust --matrixfile data/reuters.mtx  --dictfile data/reuters_dictionary.txt --clusters 10 --maxterms 8 --format JSON

Two files will be generated: ``tree_10.json`` and ``assignments_10.csv``.  The json file will have 8 keywords per node, whereas the ``tree_10.xml`` file will have only 5.

To generate a flat clustering result (in addition to the hierarchical clustering result), use this command line:

.. code-block:: none

   hierclust/bin/hierclust --matrixfile data/reuters.mtx  --dictfile data/reuters_dictionary.txt --clusters 10 --maxterms 8 --format JSON --flat 1

Two additional files will be generated this time (along with ``tree_10.json`` and ``assignments_10.csv``): ``clusters_10.json``, which contains the flat clustering results, and ``assignments_flat_10.csv``, which contains the flat clustering assignments.

*********
Flatclust
*********

Overview
========

The flatclust command line application factors the input matrix using either NMF-HALS or NMF-BPP and generates a flat clustering result.  A flatclust run generating k clusters will generally run more slowly than a hierclust run, of the same number of clusters, with the --flat option enabled.  The reason for this is that the hierclust application uses the NMF-RANK2 algorithm and always generates factor matrices with two rows or columns.  The runtime of NMF scales superlinearly with k in this case, and thus runs fastest for the smallest k value. 

The flatclust application generates two output files.  The first file contains the assignments of documents to clusters and is interpreted identically to that of the hierclust application, with the exception that there are no outliers generated by flatclust.

The second file contains the node information.  This file is much simpler than that of the hierclust application since there is no factorization tree.  The items for each node in this file are:

    1.  ``id``: the unique id of this node
    2.  ``doc_count``: the number of documents assigned to this node
    3.  ``top_terms``: the highest probability dictionary terms assigned to this node

Command Line Options
====================

Running the flatclust application with no command line parameters will cause the application to display all params that it supports.  These are:

.. code-block:: none

	Usage: flatclust
        --matrixfile <filename>      Filename of the matrix to be factored.
                                     Either CSV format for dense or MatrixMarket format for sparse.
        --dictfile <filename>        The name of the dictionary file.
        --clusters <integer>         The number of clusters to generate.
        [--algorithm  BPP]           The NMF algorithm to use: 
                                         HALS:  hierarchical alternating least squares
                                         RANK2: rank2 with optimal active set selection
                                                (for two clusters only)
                                         BPP:   block principal pivoting
        [--infile_W  (empty)]        Dense matrix to initialize W, CSV file.
                                     The matrix has m rows and 'clusters' columns.
                                     If unspecified, W will be randomly initialized.
        [--infile_H  (empty)]        Dense matrix to initialize H, CSV file. 
                                     The matrix has 'clusters' rows and n columns.
                                     If unspecified, H will be randomly initialized. 
        [--tol  0.0001]              Tolerance value for the progress metric. 
        [--outdir  (empty)]          Output directory.  If unspecified, results will be 
                                     written to the current directory.
        [--miniter  5]               Minimum number of iterations to perform.
        [--maxiter  5000]            Maximum number of  iterations to perform. 
        [--maxterms  5]              Number of terms per node. 
        [--maxthreads   8]           Upper limit to thread count. 
        [--verbose  1]               Whether to print updates to the screen.
                                         1 == yes, 0 == no
        [--format  XML]              Format of the output file containing the tree.
                                         XML: XML format
                                         JSON: JavaScript Object Notation
        [--clustfile clusters_N.ext] Name of the output XML file containing the tree.
                                     N is the number of clusters for this run.
                                     The string 'ext' depends on the desired format.
                                     This filename is relative to the outdir.
        [--assignfile assignments_N.csv]  Name of the file containing final assignments.
                                          N is the number of clusters for this run.
                                          This filename is relative to the outdir.

The ``--matrixfile``, ``--dictfile``, and ``--clusters`` options are required; all others are optional and have the default values indicated.  The meanings of the various options are as follows:

	1.  ``--matrixfile``: Filename of the matrix to be factored.  CSV files are supported for dense matrices and MTX (matrix market) files for sparse matrices.
	2.  ``--dictfile``: absolute or relative path to the dictionary file
	3.  ``--clusters``: the number of clusters to generate (equivalent to the NMF ``k`` value)
	4.  ``--algorithm``: the factorization algorithm to use 
	5.  ``--infile_W``: CSV file containing the m x ``clusters`` initial values for matrix W; if omitted, W is randomly initialized
	6.  ``--infile_H``:  CSV file containing the ``clusters`` x n initial values for matrix H; if omitted, H is randomly initialized
	7.  ``--tol``: tolerance value for the factorization; the stopping criterion is the ratio of projected gradient method
	8.  ``--outdir``: path to the folder into which to write the output files; if omitted results will be written to the current directory
	9.  ``--miniter``: minimum number of iterations to perform before checking progress 
	10. ``--maxiter``: the maximum number of iterations to perform 
	11. ``--maxterms``: the number of dictionary keywords to include in each node
	12. ``--maxthreads``: the maximum number of threads to use; the default is to use as many threads as the hardware can support (your number may differ from that shown) 
	13. ``--verbose``: whether to display updates to the screen as the iterations progress
	14. ``--format``: file format to use for the clustering results
	15. ``--clustfile``: name of the output file for the nodes; uses the format specified by the format parameter
	16. ``--assignfile``: name of the output file for the cluster assignments

Sample Runs
===========

The smallk distribution has available a ``smallk_data`` repository on github with a matrix file ``reuters.mtx`` and an associated dictionary file ``reuters_dictionary.txt``.  These files are derived from the popular Reuters data set used in machine learning experiments.  

As above, it is assumed that the `smallk_data <https://github.com/smallk/smallk_data>`_ repository was cloned into ``data`` and that the commands can be run as below or from ``/usr/local/bin``.

Suppose we want to perform flat clustering on this data set and generate 10 clusters.  We would do that as follows, assuming that we are in the top-level smallk folder after building the code:

.. code-block:: none

   flatclust/bin/flatclust --matrixfile data/reuters.mtx  --dictfile data/reuters_dictionary.txt --clusters 10

This will generate two result files in the current directory: ``clusters_10.xml`` and ``assignments_10.csv``.

If we want to instead generate 10 clusters, each with 8 terms, using JSON output format, we would use this command line:

.. code-block:: none

   flatclust/bin/flatclust --matrixfile data/reuters.mtx  --dictfile data/reuters_dictionary.txt --clusters 10 --maxterms 8 --format JSON

Two files will be generated: ``clusters_10.json`` and ``assignments_10.csv``.  The json file will have 8 keywords per node, whereas the ``clusters_10.xml`` file will have only 5.

