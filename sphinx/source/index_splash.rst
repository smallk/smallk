**********************
Introduction to SmallK
**********************

SmallK is a high performance software package for low rank matrix approximation via the nonnegative matrix factorization (NMF). NMF is a low-rank approximation where a matrix is approximated by a product of two nonnegative factors. The role of NMF in data analytics has been as significant as the singular value decomposition (SVD). However, due to nonnegativity constraints, NMF has far superior interpretability of its results for many practical problems such as image processing, chemometrics, bioinformatics, topic modeling for text analytics and many more. Our approach to solving the NMF nonconvex optimization
problem has proven convergence properties and is one of the most efficient methods developed to date.

********************
Distributed Versions
********************

Recently open sourced: MPI-FAUN! Both MPI and OPENMP implementations for MU, HALS and ANLS/BPP based NMF algorithms are now available. The implementations can run off the shelf or can be easily integrated into other source code. These are very highly tuned NMF algorithms to work on super computers. We have tested this software in NERSC as well OLCF cluster. The openmp implementation is tested on many different linux variants with intel processors. The library works well for both sparse and dense matrices.

Please visit `MPI-FAUN <https://github.com/ramkikannan/nmflibrary>`_ for more information and source code.

**************************************************************
Ground truth data for graph clustering and community detection
**************************************************************

Community discovery is an important task for revealing structures in large networks. The massive size of contemporary social networks poses a tremendous challenge to the scalability of traditional graph clustering algorithms and the evaluation of discovered communities. 

Our methodology uses a divide-and-conquer strategy to discover hierarchical community structure, non-overlapping within each level. Our algorithm is based on the highly efficient Rank-2 Symmetric Nonnegative Matrix Factorization. We solve several implementation challenges to boost its efficiency on modern CPU architectures, specifically for very sparse adjacency matrices that represent a wide range of social networks. Empirical results have shown that our algorithm has competitive overall efficiency, and that the non-overlapping communities found by our algorithm recover the ground-truth communities better than state-of-the-art algorithms for overlapping
community detection.

Please visit `dblp ground truth data <https://github.com/smallk/smallk_data/tree/master/dblp_ground_truth>`_ to obtain the data.

.. include:: acknowledgements.inc

******************************
Copyright and Software License
******************************

SmallK is under copyright by the Georgia Institute of Technology, 2017. All source code is released under the `Apache 2.0 <http://www.apache.org/licenses/LICENSE-2.0>`_ license.

