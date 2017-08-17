Notes for smallk release 2017/07/21:

      1. All code compiled with gcc 7.1.0.
      2. Vagrant installation upgraded.
      3. Docker installation available.
      4. OSX Sierra SIP (system integrity protection) issue resolved.
      
      dblp: computer science bibliography ground truth data for graph analytics

We provide new data sets of the DBLP computer science bibliography network with richer metadata and verifiable ground-truth knowledge, which can foster future research in community finding and interpretation of communities in large networks.

There are six files in total:

dblp15_graph.mtx The adjacency matrix of the graph
dblp15_graph_weighted.mtx Weighted adjacency matrix, the weight means how many times two authors have collaborated
dblp15_ground_truth.mtx ground truth matrix,  where the (i,j) entry equaling 1 means that author i published in venue j
dblp15_ground_truth_split.mtx split ground truth matrix, where the original ground truth communities are split into connected components
dblp15_authors.txt list of author names, as appeared in the dblp.xml file, the order of which is consistent with all the matrices
dblp15_venues.txt list of venue keys, as described in the paper, the order of which is consistent with the matrix in dblp15_ground_truth.mtx
Community discovery is an important task for revealing structures in large networks. The massive size of contemporary social networks poses a tremendous challenge to the scalability of traditional graph clustering algorithms and the evaluation of discovered communities. Our methodology uses a divide-and-conquer strategy to discover hierarchical community structure, non-overlapping within each level. Our algorithm is based on the highly efficient Rank-2 Symmetric Nonnegative Matrix Factorization. We solve several implementation challenges to boost its efficiency on modern CPU architectures, specifically for very sparse adjacency matrices that represent a wide range of social networks. Empirical results have shown that our algorithm has competitive overall efficiency, and that the non-overlapping communities found by our algorithm recover the ground-truth communities better than state-of-the-art algorithms for overlapping community detection. These results are part of an upcoming publication cited below.

Rundong Du, Da Kuang, Barry Drake and Haesun Park, Georgia Institute of Technology, "Hierarchical Community Detection via Rank-2 Symmetric Nonnegative Matrix Factorization", submitted 2017.
