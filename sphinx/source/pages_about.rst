#####
About
#####

.. toctree::
   :maxdepth: 2

SmallK is a high performance software package for low rank matrix approximation via the nonnegative matrix factorization (NMF). Algorithms for NMF compute the low rank factors of a matrix producing two nonnegative matrices whose product approximates the original matrix. The role of NMF in data analytics has been as significant as the singular value decomposition (SVD). However, due to nonnegativity constraints, NMF has far superior interpretability of its results for many practical problems such as image processing, chemometrics, bioinformatics, topic modeling for text analytics and many more. Our approach to solving the NMF nonconvex optimization problem has proven convergence properties and is one of the most efficient methods developed to date.

********************
Distributed Versions
********************

Recently open sourced: MPI-FAUN! Both MPI and OPENMP implementations for MU, HALS and ANLS/BPP based NMF algorithms are now available. The implementations can run off the shelf or can be easily integrated into other source code. These are very highly tuned NMF algorithms to work on super computers. We have tested this software in NERSC as well OLCF cluster. The openmp implementation is tested on many different linux variants with intel processors. The library works well for both sparse and dense matrices.

Please visit `MPI-FAUN text <https://github.com/ramkikannan/nmflibrary>`_ for more information and source code.

**************************************************************
Ground truth data for graph clustering and community detection
**************************************************************

Community discovery is an important task for revealing structures in large networks. The massive size of contemporary social networks poses a tremendous challenge to the scalability of traditional graph clustering algorithms and the evaluation of discovered communities. 

Please visit `dblp ground truth data <https://github.com/smallk/smallk_data/tree/master/dblp_ground_truth>`_ to obtain the data.

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
`barry.drake@gtri.gatech.edu <mailto:barry.drake@gtri.gatech.edu>`_

Stephen Lee-Urban 
Research Scientist
Information and Communications Laboratory (ICL)
Information and Cyber Sciences Directorate (ICSD)
Georgia Tech Research Institute (GTRI) 
75 5TH St. NW STE 900 
ATLANTA, GA 30308-1018
`stephen.lee-urban@gtri.gatech.edu <stephen.lee-urban@gtri.gatech.edu>`_
