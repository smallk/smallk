This data set contains 13 subsets of the US Patent data set obtained from http://www.patentsview.org in August 2016, corresponding to 13 Cooperative Patent Classification (CPC) classes: A22 A42 B06 B09 B68 C06 C13 C14 C40 D02 D10 F22 Y04.

In each sub-folder, there are the following files:

* term_frequency.mtx     The term-document matrix of the patent claims, where the (i,j) entry is the frequency of term i appeared in the claims of patent j;
* graph_directed.mtx     The adjacency matrix of the citation network, where the (i,j) entry equals one if and only if patent i cites patent j;
* dictionary.txt     List of terms, where the order of terms is consitent with the matrices;
* documents.txt     List of patent ids, where the order of patents is consitent with the matrices;
* ground_truth.mtx     Ground truth matrix, where the (i,j) entry equaling 
1 means that patent i belongs to group j
* S.mtx     Normalized graph matrix, as used in [1]
* X.mtx     Normalized term-document matrix, as used in [1]

1. Rundong Du, Barry Drake and Haesun Park. Hybrid Clustering based on Content and Connection Structure using Joint Nonnegative Matrix Factorization. Journal of Global Optimization, to appear.
