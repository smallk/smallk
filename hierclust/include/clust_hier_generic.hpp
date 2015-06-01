// Copyright 2014 Georgia Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <map>
#include <list>
#include <thread>
#include <sstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include "clust.hpp"
#include "dense_matrix.hpp"
#include "clust_hier_util.hpp"
#include "matrix_generator.hpp"
#include "nmf_solve_generic.hpp"
#include "setdiff.hpp"
#include "tree.hpp"
#include "terms.hpp"
#include "random.hpp"
#include "file_loader.hpp"

//#include <iomanip>
//#include "timer.hpp"

// TODO:
//     -- done -- 1. Much better to initialize matrices W and H in ActualSplit, since only
//                   a small submatrix is needed.
//     -- done -- 2. Add code for top term generation.
//     -- done -- 3. Add code to build tree on the fly.
//     -- done -- 4. Add code for generation of output files.
//     -- done -- 5. Remove the topic_vector argument; get from tree.  Use leaf node
//                   topic vectors to initialize flat clustering.
//     -- done -- 6.  Add outliers to hierclust assignment file.
//     7.  Move W_buffer and H_buffer into the tree class
//     8.  Initialize W and H in separate threads; use two random number generators
//     9.  What happens if a node's document count goes to 0?
//    10.  What to do if no further factorization possible?  Compute top terms for
//         existing nodes?  Compute assignments so far?

namespace
{
    template <typename T,
              template <typename> class MatrixType>
    struct ClustData
    {
        static MatrixType<T> Asubset;
    };

    // these are here to avoid repeated allocation

    // definition of the static member Asubset
    template <typename T, 
              template <typename> class MatrixType>
    MatrixType<T> ClustData<T, MatrixType>::Asubset;

    std::vector<unsigned int> old_to_new_rows, new_to_old_rows;
}

//-----------------------------------------------------------------------------
template <typename T,  
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
bool ClustHier(const MatrixType<T>& A,     // m x n
               Solver<T, MatrixType>& solver,
               ProgressEst<T, MatrixType>* progress_est,
               const ClustOptions& clust_opts,
               Tree<T>& tree,
               ClustStats& clust_stats,
               Random& rng)
{
    T min_priority, max_priority, priority_one;
    
    // This code assumes that the input matrix A has no duplicate columns
    // and no zero rows.

    unsigned int m = A.Height();
    unsigned int n = A.Width();

    // a nonempty initdir means to load initializer matrices for W and H
    bool load_initializers = !clust_opts.initdir.empty();
    
    // Each cluster is represented by a leaf node in the binary tree.  A binary
    // tree with 'num_clusters' leaf nodes contains 2*(num_clusters - 1) nodes
    // (excluding the root node), since each split generates two nodes, and a
    // total of 'num_clusters' - 1 splits are performed.

    // The root node is not represented in any of the data structures below;
    // only the non-root nodes matter.
    
    unsigned int num_clusters = clust_opts.num_clusters;
    if (num_clusters <= 1)
        throw std::runtime_error("HierNMF2: number of clusters must be >= 2");

    unsigned int node_count = 2*(num_clusters - 1);    

    tree.Init(num_clusters, node_count, m, n);
    
    //-------------------------------------------------------------------------
    //
    //                      Factor the root node.
    //
    //-------------------------------------------------------------------------
    
    NmfStats nmf_stats;
    DenseMatrix<T> W0(m, 2), H0(2, n);

    // make at most three attempts to factor the root
    bool success = false;
    for (unsigned int q=0; q<3; ++q)
    {
        // initialize W and H
        if (load_initializers)
            LoadInitializers(clust_opts.initdir, m, n, W0, H0);
        else
            RandomInitializers(rng, W0, H0);

        // run Rank2 NMF on matrix A
        success = NmfSolve(A, W0, H0, solver, progress_est,
                           clust_opts.nmf_opts, nmf_stats);

        if (success)
        {
            clust_stats.nmf_count += 1;
            if (nmf_stats.iteration_count == clust_opts.nmf_opts.max_iter)
                clust_stats.max_count += 1;

            break;
        }

        // If here the factorization failed, probably because a singular
        // matrix was encountered.  Re-initialize W and H and try again.
        std::cout << "\nRoot node factorization failed, "
                  << "retrying with new initializers..." << std::endl;
    }

    if (!success)
        throw std::runtime_error("HierNMF2: root node factorization failed after three attempts");

    unsigned int split_index = 0;

    // W and H point to the most recently computed solutions for W and H.
    DenseMatrix<T> *W = &W0, *H = &H0;
    std::vector<DenseMatrix<T> > W_buffer(node_count);
    std::vector<DenseMatrix<T> > H_buffer(node_count);

    //-------------------------------------------------------------------------
    //
    //                        Main loop
    //
    //-------------------------------------------------------------------------
    for (unsigned int i=0; i<(num_clusters - 1); ++i)
    {
        if (0 == i)
        {
            min_priority = INFINITY;
            tree.SplitRoot(W, H);
        }
        else
        {
            // find min/max priority leaf nodes and index of max priority node
            tree.MinMaxLeafPriorities(min_priority, max_priority, split_index);
            if (max_priority < T(0))
            {
                std::cout << "\nHierNMF2: no further factorization possible.\n"
                          << std::endl;
                break;
            }
            
            W = &(W_buffer[split_index]);
            H = &(H_buffer[split_index]);
            tree.Split(split_index, W, H);
        }

        unsigned int index0 = tree.LeftChildIndex();
        unsigned int index1 = tree.RightChildIndex();
        assert(index0 < node_count);
        assert(index1 < node_count);
        
        priority_one = TrialSplit(A,
                                  tree.LeftChildDocs(),
                                  min_priority,
                                  tree.LeftChildTopicVector(),
                                  W_buffer[index0],
                                  H_buffer[index0],
                                  clust_opts,
                                  solver,
                                  progress_est,
                                  rng,
                                  clust_stats);

        tree.SetNodePriority(index0, priority_one);
        
        priority_one = TrialSplit(A,
                                  tree.RightChildDocs(),
                                  min_priority,
                                  tree.RightChildTopicVector(),
                                  W_buffer[index1],
                                  H_buffer[index1],
                                  clust_opts,
                                  solver,
                                  progress_est,
                                  rng,
                                  clust_stats);

        tree.SetNodePriority(index1, priority_one);
        
        // print iteration number to the screen
        if (clust_opts.verbose)
        {
            std::cout <<"[" << (i+1) << "] ";
            std::cout.flush();
        }
    }

    tree.ComputeTopTerms(clust_opts.maxterms);
    tree.ComputeAssignments();

    // print debug info if using init matrices
    //if (load_initializers && clust_opts.verbose)
    //    tree.Print();
    
    std::cout << std::endl;
    return true;
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
T TrialSplit(const MatrixType<T>& A,
             std::vector<unsigned int>& subset,          // [in, out]
             const T min_priority,                       // [in]
             DenseMatrix<T>& W_parent,                   // [in]
             DenseMatrix<T>& W,
             DenseMatrix<T>& H,
             const ClustOptions& clust_opts,             // [in]
             Solver<T, MatrixType>& solver,              // [in]
             ProgressEst<T, MatrixType>* progress_est,   // [in]
             Random& rng,                                // [in]
             ClustStats& clust_stats)                    // [out]
{
    std::vector<int> counts(2);

    // save a backup in case two clusters can't be found
    std::vector<unsigned int> subset_backup(subset), subset_small, cluster_subset;

    int trial = 0;
    T priority_one = T(-2);
    DenseMatrix<T> Wtemp, Htemp;
    
    while (trial < clust_opts.trial_allowance)
    {
        // attempt to split the current node
        priority_one = ActualSplit(A, subset, W_parent, W, H,
                                   clust_opts, solver, progress_est,
                                   rng, cluster_subset, clust_stats);

        assert(cluster_subset.size() == subset.size());
        
        if (priority_one < T(0))
            break;

        // We expect two clusters from the splitting process, each labeled by
        // 0 or 1.  Ensure that two clusters have actually been generated, and
        // determine the size of each.
        counts[0] = counts[1] = 0;
        for (unsigned int q=0; q<cluster_subset.size(); ++q)
        {
            unsigned int label = cluster_subset[q];
            if ( (0 == label) || (1 == label))
                counts[label] += 1;
            else
                throw std::runtime_error("TrialSplit: found an invalid subcluster");
        }

        // If a cluster is too small, compute its priority and compare with the
        // existing minimum priority score.  If this cluster has an even lower
        // priority score, subtract this cluster from the original 'subset' cluster
        // and keep trying for an improved factorization.

        int smallest_size = std::min(counts[0], counts[1]);
        if (smallest_size < clust_opts.unbalanced * cluster_subset.size())
        {
            unsigned int label_of_smallest = ( (smallest_size == counts[0]) ? 0 : 1);

            // extract all cols from 'subset' at indices matching the smallest label
            subset_small.clear();
            for (unsigned int q=0; q<cluster_subset.size(); ++q)
            {
                if (label_of_smallest == cluster_subset[q])
                    subset_small.push_back(subset[q]);
            }

            // Compute the priority of this small subset, which requires a
            // factorization and the creation of two potential children.
            // The topic vector for this factorization is the column of W
            // (computed via the previous factorization) at the index of
            // the smallest label.

            std::vector<unsigned int> cluster_subset_small;
            DenseMatrix<T> W_col, W_buffer_one_small, H_buffer_one_small;

            // create a 'view' of the appropriate column of W and factor again
            //View(W_col, *W_buffer_one, 0, label_of_smallest, m, 1);
            View(W_col, W, 0, label_of_smallest, W.Height(), 1);
            T priority_one_small = ActualSplit(A, subset_small, W_col,
                                               Wtemp,
                                               Htemp,
                                               clust_opts, solver,
                                               progress_est, rng,
                                               cluster_subset_small,
                                               clust_stats);

            assert(cluster_subset_small.size() == subset_small.size());

            // Now we have the priority score of the cluster 'subset_small'.
            // If this cluster has the smallest priority score so far, remove
            // the members of 'subset_small' from 'subset' and keep trying for
            // an improved factorization.
            if (priority_one_small < min_priority)
            {
                trial += 1;
                if (trial < clust_opts.trial_allowance)
                {
                    std::cout << "dropping " << subset_small.size()
                              << " items ..." << std::endl;
                    subset = SetDiff(subset, subset_small);
                }
            }
            else
            {
                break;
            }
        }
        else
        {
            break;
        }
    }

    assert(W.Height() > 0); assert(W.Width() > 0);
    assert(H.Height() > 0); assert(H.Width() > 0);
    
    assert(trial <= clust_opts.trial_allowance);
    if (trial == clust_opts.trial_allowance) 
    {
        // exhausted all split attempts, so this node becomes a permanent leaf
        if (clust_opts.verbose)
        {
            std::cout << "recycling " << subset_small.size() << " items ..." << std::endl;
        }

        subset = subset_backup;
        MakeZeros(W);
        H.Resize(2, subset.size());
        MakeZeros(H);
        priority_one = T(-2);
    }

    return priority_one;
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
T ActualSplit(const MatrixType<T>& A,                     // [in]
              std::vector<unsigned int>& subset,          // [in]
              DenseMatrix<T>& W_parent,                   // [in]
              DenseMatrix<T>& W,                          // [out] resized by this function
              DenseMatrix<T>& H,                          // [out] resized by this function
              const ClustOptions& clust_opts,             // [in]
              Solver<T, MatrixType>& solver,              // [in]
              ProgressEst<T, MatrixType>* progress_est,   // [in]
              Random& rng,                                // [in]
              std::vector<unsigned int>& cluster_subset,  // [out]
              ClustStats& clust_stats)                    // [out]
{
    unsigned int m = A.Height();
    bool load_initializers = !clust_opts.initdir.empty();
    
    DenseMatrix<T> Wtmp, Htmp;
    
    // W and H are assumed to contain initializer values when this function
    // is called.  Upon return, W and H are overwritten with the computed
    // factors.
    
    if (subset.size() <= 3) 
    { 
        // quit when the number of items in this subset is too small
        cluster_subset.resize(subset.size());
        for (unsigned int q=0; q<cluster_subset.size(); ++q)
            cluster_subset[q] = 1u;

        W.Resize(m, 2);
        H.Resize(2, subset.size());
        MakeZeros(W);
        MakeZeros(H);
        return T(-1);
    } 
    else 
    {
        // Create a new sparse matrix 'Asubset' consisting only of the cols
        // of matrix A contained in 'subset'.  This will be the matrix to
        // be factored.  Extract the most compact submatrix possible, since
        // extraction of subcolumns may result in a submatrix with a
        // different distribution of row indices, and hence a different height.
        A.SubMatrixColsCompact(ClustData<T, MatrixType>::Asubset, 
                               subset, old_to_new_rows, new_to_old_rows);

        unsigned int new_height = ClustData<T, MatrixType>::Asubset.Height();

        DenseMatrix<T> Wsubset(new_height, 2);
        DenseMatrix<T> Hsubset(2, subset.size());
                
        // perform the rank2 factorization; make at most three attempts
        NmfStats nmf_stats;
        bool success = false;
        for (unsigned int q=0; q<3; ++q)
        {
            // Create new initializer matrices for W and H by extracting the
            // appropriate subsets from the original initializers.  These
            // extractions only take a miniscule fraction of the runtime, so
            // no need to parallelize.  The NMF solver dominates the runtime.
            
            if (load_initializers)
            {
                LoadInitializers(clust_opts.initdir, A.Height(), A.Width(), Wtmp, Htmp);
                ExtractSubmatrices(subset, new_to_old_rows, Wsubset, Hsubset, Wtmp, Htmp);
            }
            else
            {
                RandomInitializers(rng, Wsubset, Hsubset);
            }
            
            success = NmfSolve(ClustData<T, MatrixType>::Asubset,
                               Wsubset, Hsubset, solver, progress_est,
                               clust_opts.nmf_opts, nmf_stats);
            if (success)
            {
                clust_stats.nmf_count += 1;
                if (clust_opts.nmf_opts.max_iter == nmf_stats.iteration_count)
                    clust_stats.max_count += 1;

                break;
            }

            // If here the factorization failed, probably because a singular
            // matrix was encountered.  Re-initialize W and H and try again.
            std::cout << "\nNode factorization failed, "
                      << "retrying with new initializers..." << std::endl;
        }

        // if failure after three attempts, something is really wrong
        if (!success)
            throw std::runtime_error("HierNMF2: node factorization failed after three attempts.");

        // Generate the subclusters (hopefully two of them); subcluster
        // membership is determied by the row index of the max value of
        // H(:, c) for all cols c.
        
        bool has_0 = false, has_1 = false;
        
        cluster_subset.clear();
        for (int c=0; c<Hsubset.Width(); ++c)
        {
            if (Hsubset.Get(0, c) > Hsubset.Get(1, c))
            {
                // element in row 0 is greater
                cluster_subset.push_back(0u);
                has_0 = true;
            }
            else
            {
                // element in row 1 is greater
                cluster_subset.push_back(1u);
                has_1 = true;
            }
        }

        // update the topic vectors for the affected rows
        W.Resize(m, 2);
        MakeZeros(W);
        for (unsigned int i=0; i<new_height; ++i)
        {
            W.Set(new_to_old_rows[i], 0, Wsubset.Get(i,0));
            W.Set(new_to_old_rows[i], 1, Wsubset.Get(i,1));
        }

        // set H_result equal to Hsubset; note that the width changes with the clusters
        H.Resize(2, subset.size());
        Copy(Hsubset, H);
        
        // compute the priority score if indeed two potential clusters are generated
        T priority = T(-1);
        if (has_0 && has_1)
            priority = compute_priority(W_parent, W);

        return priority;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ExtractSubmatrices(std::vector<unsigned int>& subset,
                        std::vector<unsigned int> new_to_old_rows,
                        DenseMatrix<T>& Wsub,
                        DenseMatrix<T>& Hsub,
                        DenseMatrix<T>& Wsource,
                        DenseMatrix<T>& Hsource)
{
    unsigned int new_height = Wsub.Height();
     
    // extract the initializer submatrix for W
    for (unsigned int r=0; r<new_height; ++r)
    {
        Wsub.Set(r, 0, Wsource.Get(new_to_old_rows[r], 0));
        Wsub.Set(r, 1, Wsource.Get(new_to_old_rows[r], 1));
    }
    
    // extract the initializer submatrix for H
    for (unsigned int c=0; c<subset.size(); ++c) 
    {
        // subset[c] is a column index
        Hsub.Set(0, c, Hsource.Get(0, subset[c]));
        Hsub.Set(1, c, Hsource.Get(1, subset[c]));
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void RandomInitializers(Random& rng,
                        DenseMatrix<T>& W,
                        DenseMatrix<T>& H)
{
//    Timer timer;
//    timer.Start();
    
    RandomMatrix(W.Buffer(), W.LDim(), W.Height(), W.Width(), rng, T(0.5), T(0.5));
    RandomMatrix(H.Buffer(), H.LDim(), H.Height(), H.Width(), rng, T(0.5), T(0.5));

//    timer.Stop();
//    std::cout << "\nInit time for W (" << std::setw(7) << W.Height()
//              << " x " << W.Width()
//              << ") and H (" << H.Height() << " x "
//              << std::setw(7) << H.Width() << "):\t"
//              << timer.ReportMicroseconds() << " us." << std::endl;
}

//-----------------------------------------------------------------------------
template <typename T>
void LoadInitializers(const std::string INIT_DIR,
                      const unsigned int m,  
                      const unsigned int n,
                      DenseMatrix<T>& W,     // W is mx2
                      DenseMatrix<T>& H)     // H is 2xn
{
    // On each call to this function, initializer matrices for W and H are
    // loaded from the INIT_DIR hardcoded below.  The matrices are assumed to
    // have the names Winit_1.csv, Hinit_1.csv, Winit_2.csv, Hinit_2.csv, etc.
    // It is up to the user to ensure that enough matrices are present in this
    // dir to run the HierNMF2 code to completion.  The number of matrices used
    // is non-deterministic, so trial-and-error may be required to find a lower
    // bound on the matrix count.

    // This function is used for testing (such as comparisons with Matlab), in
    // which each factorization problem has to proceed from a known
    // initializer.

    static int counter = 1;  // file counter; initial factorization uses matrix 1
    
    W.Resize(m, 2);
    T* buf = W.Buffer();
    std::ostringstream filename_w;
    filename_w << INIT_DIR << "Winit_" << counter << ".csv";
    if (!LoadDelimitedFile(buf, m, 2, filename_w.str()))
    {
        std::string msg = "Load failed for file " + filename_w.str();
        throw std::runtime_error(msg);
    }
    
    H.Resize(2, n);
    buf = H.Buffer();
    std::ostringstream filename_h;
    filename_h << INIT_DIR << "Hinit_" << counter << ".csv";
    if (!LoadDelimitedFile(buf, 2, n, filename_h.str()))
    {
        std::string msg = "Load failed for file " + filename_h.str();
        throw std::runtime_error(msg);
    }

    // one more set of initializer matrices has now been loaded
    ++counter;

    //std::cout << "counter: " << counter << std::endl;
    
    //std::cout << "\t\tLoaded " << filename_w.str() << std::endl;
    //std::cout << "\t\tNorm: " << Norm(Wtmp) << std::endl;
    //std::cout << "\t\t" << Wtmp.Get(0, 0) << "  " << Wtmp.Get(0, 1) << std::endl;
    //std::cout << "\t\t" << Wtmp.Get(1, 0) << "  " << Wtmp.Get(1, 1) << std::endl;
    
    //std::cout << "\t\tLoaded " << filename_h.str() << std::endl;
    //std::cout << "\t\tNorm: " << Norm(Htmp) << std::endl;
    //std::cout << "\t\t" << Htmp.Get(0, 0) << "  " << Htmp.Get(0, 1) << std::endl;
    //std::cout << "\t\t" << Htmp.Get(1, 0) << "  " << Htmp.Get(1, 1) << std::endl;
}
