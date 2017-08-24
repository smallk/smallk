###########
Quick Start
###########

.. contents:: :local:

************
Introduction
************

This document describes how to use the SmallK library to perform nonnegative matrix factorization (NMF), hierarchical clustering, and flat clustering. It is assumed that the library has been installed properly, that all tests have passed, and that the user has created the SMALLK_INSTALL_DIR environment variable as described in the documentation. 
SmallK provides a very simple interface to NMF and clustering algorithms. Examples of how to use this interface are described in this document. The SmallK distribution also provides a suite of command-line tools for NMF and clustering, suitable for advanced users.

*****************
C++ Project Setup
*****************

The SmallK distribution includes an `examples` folder containing two files: smallk_examples.cpp and a Makefile. To build the example CPP file, open a terminal window, cd to the smallk/examples folder, and run the command `make`.

If the SmallK library has been installed properly and the `smallk_data <https://github.com/smallk/smallk_data>`_ repository has been cloned at the same directory level as the SmallK library repository, the project should build and the binary file bin/example will be created.  To run the example, run this command from the smallk/examples folder::

	./bin/example ../../smallk_data

Results will appear for the following algorthms::

	Running NMF-BPP using k=32
	Running NMF-HALS using k=16
	Running NMF-RANK2 with W and H initializers
	Repeating the previous run with tol = 1.0e-5
	Running HierNMF2 with 5 clusters, JSON format
	Running HierNMF2 with 10 clusters, 12 terms, XML format
	Running HierNmf2 with 18 clusters, 8 terms, with flat
		
The output files will be written to the directory where the binary ‘example’ is run. In the above, the outputs will be written to the <SmallK dir>/examples.

To experiment with the SmallK library, make a backup copy of smallk_examples.cpp as follows::

	cp smallk_examples.cpp smallk_examples.cpp.bak

The file smallk_examples.cpp can now be used for experimentation. The original file can be restored from the backup at the user’s discretion.

Delete lines 61-255 from smallk_examples.cpp (everything between the opening and closing braces of the ‘try’ block). New code will be added between these braces in the steps below.

All of the examples described in this document use a matrix derived from Reuters articles. This matrix will be referred to as the ‘Reuters’ matrix. It is a sparse matrix with 12411 rows and 7984 columns.

The SmallK documentation contains complete descriptions of all SmallK functions mentioned in this guide.

*************
Load a Matrix
*************

Suppose you want to perform NMF or clustering on a matrix. The first action to take is to load the matrix into SmallK using the *LoadMatrix* function. This function accepts either dense matrices in CSV format or sparse matrices in MatrixMarket format.  Since we want to perform NMF and clustering on the Reuters matrix, we need to supply the path to the Reuters matrix file (reuters.mtx) as an argument to LoadMatrix.  This path has already been setup in the code; the appropriate string variable is called `filepath_matrix`.  Enter the following line after the opening brace of the try block after line 61::

	smallk::LoadMatrix(filepath_matrix);

Save the file and run the following commands, which should complete without error::

	make clean
	make

Once a matrix is loaded into SmallK it remains loaded until it is replaced with a new call to LoadMatrix. Thus, SmallK makes it easy to experiment with different factorization or clustering parameters, without having to reload a matrix each time.

********************************
Perform NMF on the Loaded Matrix
********************************

Having loaded the Reuters matrix, we can now run different NMF algorithms and factor the matrix in various ways. The SmallK code factors the loaded matrix (denoted by A) as A ~ W*H, where A is mxn, W is mxk, and H is kxn.  The NMF is a low-rank approximation where the value of k, the rank, is an input parameter to the factorization routines, and is generally much smaller than either m or n. Matrix A can be either sparse or dense; matrices W and H are always dense.

NMF-BPP
=======

Let’s use the default NMF-BPP algorithm to factor the 12411 x 7984 Reuters matrix into W and H with a k value of 32.  Add the following lines to the code::

	MsgBox(“Running NMF-BPP using k=32”);
	smallk::Nmf(32);

Build the code as described above; then run it with this command::

	./bin/example ../smallk_data

The MsgBox function prints the string supplied as argument to the screen; this function is purely for annotating the output.  The Nmf function performs the factorization and generates two output files, w.csv and h.csv, which contain the matrix factors.  The files are written to the current directory.  SmallK can write these files to a specified output directory via the SetOutputDir function, but we will use the current directory for the examples in this guide.

NMF-HALS
========

Now suppose we want to repeat the factorization, this time using the NMF-HALS algorithm with a k value of 16.  Since the BPP algorithm is the default, we need to explicitly specify the algorithm as an argument to the Nmf function.  Add these lines to the code::

	MsgBox(“Running NMF-HALS using k=16”)
	smallk::Nmf(16, smallk::Algorithm::HALS);

Build and run the code again; you should observe that the code now performs two separate factorizations.

NMF Initialization
==================

The SmallK library provides the capability to explicitly initialize the W and H factors.  For the previous two examples, these matrices were randomly initialized, since no initializers were provided in the call to the Nmf function. The data directory contains initializer matrices for the W and H factors of the Reuters matrix, assuming that k has a value of 2. To illustrate the use of initializers, we will use the RANK2 algorithm to factor the Reuters matrix again, using a k-value of 2, but with explicit initializers.  Add these lines to the code::

	MsgBox("Running NMF-RANK2 with W and H initializers");
	smallk::Nmf(2, smallk::Algorithm::RANK2, filepath_w, filepath_h);

Build and run the code again, and observe that the code performs three separate factorizations.

The string arguments `filepath_w` and `filepath_h` are configured to point to the W and H initializer matrices in the data directory. Note how these are supplied as the third and fourth arguments to Nmf. For general matrix initializers, the W initializer must be a fully-dense matrix, in CSV format, with dimensions mxk, and the H initializer must be a fully-dense matrix, in CSV format, with dimensions kxn.

The main purpose of using initializer matrices is to generate deterministic output, such as for testing, benchmarking, and performance studies. You will notice that if you run the code repeatedly, the first two factorizations, which use random initializers, generate results that vary slightly from run to run. The third factorization, which uses initializers, always generates the same output on successive runs.

Typically the use of initializers is not required.

***********************
Hierarchical Clustering
***********************

Now let’s perform hierarchical clustering on the Reuters matrix. To do this, we must first load the dictionary (or vocabulary) file associated with the Reuters data (a file called `reuters_dictionary.txt`).  A string variable containing the full path to this file is provided in the ‘filepath_dict’ variable.  Add the following line to the code to load the Reuters dictionary::

	smallk::LoadDictionary(filepath_dict);

As with the matrix file, the dictionary file remains loaded until it is replaced by another call to LoadDictionary.

With the matrix file and the dictionary file both loaded, we can perform hierarchical clustering on the Reuters data. For the first attempt we will generate a factorization tree containing five clusters.  The number of clusters is specified as an argument to the clustering function.   Add these lines to the code::

	MsgBox("Running HierNMF2 with 5 clusters, JSON format");
	smallk::HierNmf2(5);

Build and run the code.

The hierarchical clustering function is called `HierNmf2`. In the call above it will generate five clusters and generate two output files. One file will be called `assignments_5.csv`, a CSV file containing the cluster labels.  The first entry in the file is the label for the first column (document) of the matrix; the second entry is the label for the second column, etc.  Any entries that contain -1 are outliers; these represent the documents that were not assigned to any cluster.

The other output file will be called `tree_5.json`, a JSON file containing the cluster information. This file contains sufficient information to unambiguously reconstruct the factorization tree.  If you open the file and examine the contents you can see the top terms assigned to each node.  Leaf nodes have -1 for their left and right child indices.  From an examination of the keywords at the leaf nodes, it is evident that this collection of Reuters documents is concerned with financial topics.

***************
Flat Clustering
***************

For the final example, let’s generate a flat clustering result in addition to the hierarchical clustering result. We will also increase the number of terms per node to 8 and the number of clusters to 18.  Add the following lines to the code::

	MsgBox("Running HierNmf2 with 18 clusters, 8 terms, with flat");
	smallk::SetMaxTerms(8);
	smallk::HierNmf2WithFlat(18);

Build and run the code.

The call to SetMaxTerms increases the number of top terms per node. The next line runs the hierarchical clustering algorithm and also generates a flat clustering result.  This time, four output files are generated.  They are:

1. ‘assignments_18.csv’: assignments from hierarchical clustering
2. ‘assignments_flat_18.csv’: assignments from flat clustering
3. ‘tree_18.json’, the hierarchical factorization tree
4. ‘clusters_18.json’, the flat clustering results
 
These examples demonstrate how easy it is to use SmallK for NMF and clustering. There are additional functions in the SmallK interface, described in the documentation, installation section, which allows users to set various parameters that affect the NMF-based algorithms of SmallK.  The default values for all such parameters are very reasonable, and most users will likely not ever need to change these parameters.

The smallk_examples.cpp file and the associated makefile can be used as a starting point for your own NMF and clustering projects.

**********
Disclaimer
**********

This software is a work in progress.  It will be updated throughout the course of the XDATA program with additional algorithms and examples.  The distributed NMF factorization routine uses sequential algorithms, but it replaces the matrices and matrix operations with distributed versions.  The GA Tech research group is working on proper distributed NMF algorithms, and when such algorithms are available they will be added to the library.  Thus the performance of the distributed code should be viewed as being the baseline for our future distributed NMF implementations.

************
Contact Info
************

For comments, questions, bug reports, suggestions, etc., contact:

.. include:: contact.inc

