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
#include <cassert>
#include <iostream>
#include <algorithm>
#include "clust.hpp"
#include "dense_matrix.hpp"
#include "clust_hier_util.hpp"
#include "matrix_generator.hpp"
#include "nmf_solve_generic.hpp"
#include "clust_hier_params.hpp"
#include "setdiff.hpp"
#include "tree.hpp"
#include "terms.hpp"
#include "random.hpp"

// In the comments of the code, we frequently mention
// the "index" and "linear index" of a leaf node.
// The index of a leaf node has the same meaning as
// the "index" in the documentation;
// The linear index of a leaf node is the index of that node
// in the sorted array of indices of all current leaf nodes.

// Generates a binary tree.  For a given cluster label, the left children
// are labeled 'label + 0', and the right children are labeled 'label + 0.5'.

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
               std::vector<DenseMatrix<T> >& matrices_W,
               std::vector<DenseMatrix<T> >& matrices_H,
               std::vector<int>& assignments,
               std::vector<DenseMatrix<T> >& topic_vectors,
               Solver<T, MatrixType>& solver,
               ProgressEst<T, MatrixType>* progress_est,
               const ClustOptions& clust_opts,
               Tree& tree,
               ClustStats& clust_stats)
{
    int m = A.Height();
    int n = A.Width();

    // num_clusters must be signed
    int num_clusters = clust_opts.num_clusters;
    if (num_clusters <= 0)
        throw std::runtime_error("ClustHier: invalid number of clusters");

    unsigned int num_initializers = matrices_W.size();

    if (matrices_H.size() != num_initializers)
        throw std::runtime_error("ClustHier: invalid number of H initializers");

    if (num_initializers != 2u*num_clusters)
        throw std::runtime_error("ClustHier: cluster-initializer mismatch");

    int max_terms = clust_opts.maxterms;
    std::vector<int> sort_indices(m);
    std::vector<int> term_indices(max_terms);
    
    T unbalanced        = clust_opts.unbalanced;
    int trial_allowance = clust_opts.trial_allowance;
    
    // Intialize some structs to be filled in during the process of building the tree

    // the level of each leaf node; e.g. the root node has level 0
    std::vector<int> levels;

    // need 'num_clusters' topic vectors, each of which is mx1
    if (topic_vectors.size() < static_cast<unsigned int>(num_clusters))
        topic_vectors.resize(num_clusters);

    for (int i=0; i<num_clusters; ++i)
        topic_vectors[i].Resize(m, 1);

    // current priority scores of all leaf nodes
    std::vector<T> priority;

    // cluster labels of all items (may have fractional part)
    std::vector<T> labels(n);

    // rank2 NMF results (mx2 W matrices only) for current leaf nodes
    std::list<DenseMatrix<T> > W_buffer;

    // number the nodes as they are created
    // TBD - move this to the tree class
    int node_index = 0;

    // Initialize the parameter struct to be passed between function calls   TBD - get rid of this
    ClustHierParams<T, MatrixType, Solver> params(clust_opts.nmf_opts,
                                                  trial_allowance, unbalanced,
                                                  labels, priority, W_buffer, 
                                                  solver);
    
    // Track the columns/rows of Winit/Hinit used in rank-2 NMF calls
    unsigned int init_used = 0;
    
    // Run the first rank-2 NMF, applied to the root node

    NmfStats nmf_stats;
    bool success = NmfSolve(A, matrices_W[0], matrices_H[0], solver, 
                            progress_est, clust_opts.nmf_opts, nmf_stats);
    assert(success);  // TBD - fix this

    clust_stats.nmf_count += 1;
    if (nmf_stats.iteration_count == clust_opts.nmf_opts.max_iter)
        clust_stats.max_count += 1;

    // Assign one of the two cluster labels to each document.  Assignments
    // are made by examining the two entries in each column of H.  Element
    // H(p, q) indicates the degree to which document q belongs to cluster p.
    // The larger of the two entries determines the right child.
    for (int i = 0; i < n; ++i) 
    {
        if (matrices_H[0].Get(0, i) < matrices_H[0].Get(1, i))
            labels[i] = static_cast<T>(0.5); // right child
        else
            labels[i] = static_cast<T>(0.0); // left child
    }

    levels.push_back(0);
    priority.push_back(std::numeric_limits<T>::infinity());
    W_buffer.push_back(DenseMatrix<T>(matrices_W[0]));

    auto it_wb = W_buffer.begin();
    
    // root node
    TreeNode node;
    node.label = 0;
    node.index = node_index++;
    node.is_left_child = false;
    node.doc_count = -1;
    node.parent_label = -1;
    node.parent_index = -1;
    node.left_child_label = -1;
    node.right_child_label = -1;
    tree.InsertNode(0, node);
    
    // iterate num_clusters-1 times
    for (int i=0; i < (num_clusters-1); ++i)
    {
        // round 'labels' to get the cluster assignments at this snapshot;
        // ignore distinction between left and right children
        for (int j=0; j<n; ++j)
            assignments[j] = (int)floor(labels[j]);

        std::map<int, int> unique_label_map;
        for (int j=0; j<n; ++j)
        {
            // outliers have a negative label
            if (assignments[j] >= 0)
            {
                int label = assignments[j];
                auto it = unique_label_map.find(label);
                if (unique_label_map.end() == it)
                {
                    // new label, so add to map with count of 1
                    unique_label_map.insert(std::make_pair(label, 1));
                }
                else
                {
                    // increment the count at this label
                    it->second += 1;
                }
            }
        }
        
        // The node to be split is the leaf node with the max priority score.
        // Store its linear index in 'split_index'.
        int split_index = 0;
        T max_priority = priority[0];
        for (unsigned int j=1; j < priority.size(); ++j) 
        {
            if (max_priority < priority[j]) 
            {
                max_priority = priority[j];
                split_index = j;
            }
        }
        
        // Quit if every current leaf node is a permanent leaf node
        if (max_priority < 0)
            break;

        // Determine the label of the leaf node to split (as indicated by 
        // 'label').  The two new leaf nodes will be added as children of the
        // current node with index 'label'.  The left child will have index 
        // 'label'.  The right child will have index 'label+split', where 
        // 'split' is the difference between the indices of the two children.
        int q=0, label=-1;
        for (auto it=unique_label_map.cbegin(); it != unique_label_map.cend(); ++it)
        {
            if (q == split_index)
            {
                label = it->first;
                break;
            }
            ++q;
        }

        assert(-1 != label);
        if (-1 == label)
            throw std::runtime_error("ClustHierSparse: invalid split node label");

        int split = 1;
        for (int j=0; j < levels[split_index]; ++j)
            split *= 2;

        // existing right children of the split node have label+0.5
        T rc_label = static_cast<T>(label) + T(0.5);
        T new_rc_label = label + split;
        int new_rc_label_int = label + split;

        // Change the cluster labels in the right children from label+0.5 to label+split
        int rc_count = 0;
        for (int j=0; j<n; ++j) 
        {
            if (labels[j] == rc_label)  // label + 0.5
            {
                assignments[j] = new_rc_label_int; // label + split
                labels[j] = new_rc_label;          // label + split
                ++rc_count;
            }
        }
        
        // Insert the new label for the right children in the unique label
        // map, along with the number of new right children.  Decrement the
        // entry for the left children by the number of new right children.
        unique_label_map.insert(std::make_pair(new_rc_label_int, rc_count));
        unique_label_map[label] -= rc_count;

        // Determine the linear indices of 'label' and 'label+split' and store
        // them in 'label_pos' and 'label_split_pos' respectively.
        //unique_label_set.insert(new_rc_label_int);
        int set_index = 0, label_pos = -1, label_split_pos = -1;
        for (auto it = unique_label_map.begin(); it != unique_label_map.end(); ++it)
        {
            if (it->first == label)
                label_pos = set_index;

            if (it->first == new_rc_label_int) 
            {
                label_split_pos = set_index;
                set_index++;
                break;
            }
            set_index++;
        }

        assert(-1 != label_pos);
        assert(-1 != label_split_pos);
        assert(label_pos < num_clusters);
        assert(label_split_pos < num_clusters);
        
        // We need to move some elements of the structs to the right by 1 position.
        // 'num_move' is the number of elements to be moved.
        int num_move = unique_label_map.size() - set_index;
        
        // Augment the structs by one element
        std::vector<T> Pr_j;
        levels.push_back(0);
        priority.push_back(0);

        // Move 'num_move' elements to the right by 1 position to make room 
        // for the new leaf node which has index 'label+split'. Also, write 
        // the topic vectors for the left and right children.
        int levels_size   = static_cast<int>(levels.size());
        int priority_size = static_cast<int>(priority.size());
        assert(levels_size >= 1);
        assert( (levels_size-num_move-1) >= 0);
        assert(priority_size >= 1);
        assert( (priority_size-num_move-1) >= 0);
        for (int j = levels_size-1; j >= levels_size-num_move; --j)
            levels[j] = levels[j-1];
        for (int j = priority_size-1; j >= priority_size-num_move; --j)
            priority[j] = priority[j-1];

        DenseMatrix<T> tv_copy_from, tv_copy_to;
        if (i > 0) 
        {
            // copy selected topic vector cols to the right by one position
            assert( (i+1) < num_clusters);
            assert( (i+2-num_move) >= 0);
            for (int j = i+1; j >= i+2-num_move; --j) 
                Copy(topic_vectors[j-1], topic_vectors[j]);
        }

        // advance the W buffer iterator to the correct position
        it_wb = W_buffer.begin();
        for (int q=0; q<split_index; ++q)
            ++it_wb;

        // copy W_buffer[split_index](:,0) to topic_vectors(:, label_pos)
        //elem::View(tv_copy_from, W_buffer[split_index], 0, 0, m, 1);
        LockedView(tv_copy_from, *it_wb, 0, 0, m, 1);
        Copy(tv_copy_from, topic_vectors[label_pos]);

        // copy W_buffer[split_index](:,1) to topic_vectors(:, label_split_pos)
        //elem::View(tv_copy_from, W_buffer[split_index], 0, 1, m, 1);
        LockedView(tv_copy_from, *it_wb, 0, 1, m, 1);
        Copy(tv_copy_from, topic_vectors[label_split_pos]);

        int w_buffer_size = static_cast<int>(W_buffer.size());
        assert(w_buffer_size >= 1);
        assert(w_buffer_size - num_move - 1 >= 0);
        // for (int j = w_buffer_size-1; j >= w_buffer_size-num_move; --j)
        //     elem::Copy(W_buffer[j-1], W_buffer[j]);
        it_wb = W_buffer.end(); 
        --it_wb;  // point to last element
        for (int j=w_buffer_size-1; j >= (w_buffer_size - num_move); --j)
            --it_wb;
        auto itsave = it_wb;
        W_buffer.insert(++it_wb, *itsave);  // insert a copy of *it_save
        
        // increment the level at the label linear index as well as at the 
        // split linear index; this indicates that the split index node has 
        // been split into two nodes
        levels[label_pos]++;
        levels[label_split_pos] = levels[label_pos];

        // get the index of the node to be split, which is at the previous level
        int parent_level = levels[label_pos]-1, parent_index;
        if (!tree.NodeIndex(parent_level, label, parent_index))
            throw std::runtime_error("ClustHierSparse: split node missing from tree");

        // update the parent node's child labels
        if (!tree.UpdateChildLabels(parent_level, label, label, new_rc_label_int))
            throw std::runtime_error("ClustHierSparse: missing parent node");

        // add a left child node at this level
        node.label = label;
        node.index = node_index++;
        node.is_left_child = true;
        node.doc_count = unique_label_map[label];
        node.parent_label = label;
        node.parent_index = parent_index;
        node.right_child_label = -1;
        node.left_child_label = -1;
        tree.InsertNode(levels[label_pos], node);

        // compute top terms for the left child
        TopTerms(max_terms, topic_vectors[label_pos], sort_indices, term_indices);
        if (!tree.Update(levels[label_pos], label, -1, -1, term_indices))
            throw std::runtime_error("ClustHierSparse: invalid left child node");

        // add a right child node at the this level
        node.label = new_rc_label_int;
        node.index = node_index++;
        node.is_left_child = false;
        node.doc_count = unique_label_map[new_rc_label_int];
        node.parent_label = label;
        node.parent_index = parent_index;
        node.left_child_label = -1;
        node.right_child_label = -1;
        tree.InsertNode(levels[label_pos], node);

        // compute top terms for the right child
        TopTerms(max_terms, topic_vectors[label_split_pos], sort_indices, term_indices);
        if (!tree.Update(levels[label_split_pos], new_rc_label_int, -1, -1, term_indices))
            throw std::runtime_error("ClustHierSparse: invalid right child node");

        // Some variables used in computing rank-2 NMF for the two new leaf nodes
        std::vector<int> subset;  // The subset of indices of the items in the new leaf node

        // Compute rank-2 NMF for the new leaf node on the left
        for (int j=0; j < n; ++j) 
        {
            if (labels[j] == label)
                subset.push_back(j);
        }

        init_used++;
        assert(init_used < num_initializers);

        // W_parent == topic_vectors[label_pos]; topic vector of new leaf node
        TrialSplit(A, subset, label, label_pos, topic_vectors[label_pos],
                   init_used, progress_est, matrices_W[init_used], matrices_H[init_used], 
                   params, clust_stats, clust_opts);

        subset.clear();

        // Compute rank-2 NMF for the new leaf node on the right
        for (int j = 0; j < n; ++j) 
        {
            if (labels[j] == label+split)
                subset.push_back(j);
        }

        init_used++;
        assert(init_used < num_initializers);

        // W_parent == topic_vectors[label_split_pos]; topic vector of new leaf node
        TrialSplit(A, subset, label+split, label_split_pos, topic_vectors[label_split_pos],
                   init_used, progress_est, matrices_W[init_used], matrices_H[init_used], 
                   params, clust_stats, clust_opts);

        // print the iteration number to the screen
        if (clust_opts.verbose)
            std::cout << "[" << (i+1) << "] ";
    }

    std::cout << std::endl;
    return true;
}

//-----------------------------------------------------------------------------
// This function attempts to precompute rank-2 NMF to split a node
// into two meaningful clusters and cache the results.
// It reads and writes the structs in 'params'.
template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
void TrialSplit(const MatrixType<T>& A, 
                std::vector<int>& subset,
                int new_label,
                int new_label_pos,
                DenseMatrix<T>& W_parent, 
                int init_used,
                ProgressEst<T, MatrixType>* progress_est,
                DenseMatrix<T>& Winit,
                DenseMatrix<T>& Hinit,
                ClustHierParams<T, MatrixType, Solver>& params,
                ClustStats& clust_stats,
                const ClustOptions& clust_opts) 
{
    int m = A.Height();
    const T OUTLIER_LABEL = ClustHierParams<T, MatrixType, Solver>::CLUST_OUTLIER_LABEL;
    const T LEAF_PRIORITY = ClustHierParams<T, MatrixType, Solver>::CLUST_LEAF_PRIORITY;

    // save a backup in case two clusters can't be found
    std::vector<int> subset_backup(subset);

    // cluster labels of the current node
    std::vector<T> cluster_subset;

    // topic vectors of the current node
    DenseMatrix<T> W_result(m, 2), W_result_small(m, 2);

    // priority score of the current node
    T subset_priority = T(0);

    std::map<T, int> unique_cluster_map;

    // 'params.trial_allowance' is the max number of split attempts
    // for the current node
    int trial = 0;
    while (trial < params.trial_allowance)
    {
        cluster_subset.clear();

        // Try splitting the current node by rank-2 NMF
        subset_priority = ActualSplit(A, subset, new_label, W_parent, 
                                      cluster_subset, W_result, progress_est, 
                                      Winit, Hinit, params, clust_stats);
        if (subset_priority < 0)
            break;

        // count the number of clusters and their sizes
        unique_cluster_map.clear();
        for (unsigned int i=0; i != cluster_subset.size(); ++i)
        {
            T label = cluster_subset[i];
            auto it = unique_cluster_map.find(label);
            if (it != unique_cluster_map.end())
            {
                // label has already been seen, so increment occurrence count
                it->second += 1;
            }
            else
            {
                // new label, so insert into map with a count of 1
                unique_cluster_map.insert(std::make_pair(label, 1));
            }
        }

        if (2 != unique_cluster_map.size())
            throw std::runtime_error("TrialSplit: invalid number of unique subclusters");
        
        // get the mininum size for the two clusters, the label, and the index

        // it1 accesses the first cluster
        auto it1 = unique_cluster_map.begin();

        // it2 accesses the second cluster
        auto it2 = it1;
        ++it2;

        int size_of_smaller = it1->second;
        T label_of_smaller = it1->first;
        int index_of_smaller = 0;
        if (it2->second < size_of_smaller)
        {
            size_of_smaller = it2->second;
            label_of_smaller = it2->first;
            index_of_smaller = 1;
        }

        // If one of the two potential clusters is too small, compute the 
        // priority score of that small cluster via rank2 NMF on that cluster.
        if (size_of_smaller < params.unbalanced * cluster_subset.size())
        {
            std::vector<int> subset_small;
            std::vector<T> cluster_subset_small;
            DenseMatrix<T> W_result_view;

            for (unsigned int i=0; i < cluster_subset.size(); ++i) 
            {
                if (cluster_subset[i] == label_of_smaller)
                    subset_small.push_back(subset[i]);
            }

            View(W_result_view, W_result, 0, index_of_smaller, m, 1);
            T temp_priority = ActualSplit(A, subset_small, new_label, 
                                          W_result_view, cluster_subset_small, 
                                          W_result_small, progress_est, 
                                          Winit, Hinit, params, clust_stats);

            // If the priority score of the small cluster is smaller than any 
            // positive priority score at present, delete the small cluster 
            // from 'subset' and try splitting the reduced subset on the 
            // next iteration.
            T min_priority = std::numeric_limits<T>::infinity();
            for (unsigned int i=0; i < params.priority.size(); ++i) 
            {
                T priority_i = params.priority[i];
                if ( (priority_i > 0) && (priority_i < min_priority))
                    min_priority = priority_i;
            }

            if (temp_priority < min_priority)
            {
                trial++;
                if (trial < params.trial_allowance) 
                {
                    if (clust_opts.verbose)
                    {
                        std::cout << "dropping " << subset_small.size() 
                                  << " items ..." << std::endl;
                    }

                    // make outliers from all entries in subset_small
                    for (unsigned int i=0; i != subset_small.size(); ++i)
                        params.labels[subset_small[i]] = OUTLIER_LABEL;

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

    auto it = params.W_buffer.begin();
    for (int q=0; q<new_label_pos; ++q)
        ++it;

    // store the results in the structs in 'params'
    if (trial == params.trial_allowance) 
    { 
        // exhausted all split attempts, so make the current node
        // a permanent leaf node
        int num_outliers = 0;
        for (unsigned int i=0; i < subset_backup.size(); ++i) 
        {
            if (params.labels[subset_backup[i]] == OUTLIER_LABEL)
                ++num_outliers;
            params.labels[subset_backup[i]] = new_label;
        }

        if (clust_opts.verbose)
        {
            std::cout << "recycling " << num_outliers 
                      << " items ..." << std::endl;
        }

        params.priority[new_label_pos] = LEAF_PRIORITY;

        //elem::MakeZeros(params.W_buffer[new_label_pos]);
        MakeZeros(*it);
    } 
    else 
    { 
        // Update the structs 
        for (unsigned int i=0; i<subset.size(); i++)
            params.labels[subset[i]] = cluster_subset[i];
        
        params.priority[new_label_pos] = subset_priority;
        //elem::Copy(W_result, params.W_buffer[new_label_pos]);
        Copy(W_result, *it);
    }
}

//-----------------------------------------------------------------------------
// This function does actual work of calling the rank-2 NMF routine.
// It *only* reads the structs in params.
template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
T ActualSplit(const MatrixType<T>& A,
              std::vector<int>& subset,
              int new_label,
              DenseMatrix<T>& W_parent,
              std::vector<T>& cluster_subset,
              DenseMatrix<T>& W_result,
              ProgressEst<T, MatrixType>* progress_est,
              DenseMatrix<T>& Winit,
              DenseMatrix<T>& Hinit,
              ClustHierParams<T, MatrixType, Solver>& params,
              ClustStats& clust_stats)
{
    if (subset.size() <= 3) 
    { 
        // quit when the number of items in this subset is too small
        cluster_subset.insert(cluster_subset.begin(), subset.size(), new_label);
        MakeZeros(W_result);
        return -1;
    } 
    else 
    {
        // convert column indices from int to unsigned int via explicit copy
        // TBD - make this unnecessary
        std::vector<unsigned int> col_indices(subset.size());
        for (unsigned int i = 0; i < subset.size(); ++i)
            col_indices[i] = subset[i];

        // Create a new sparse matrix 'Asubset' consisting only of the cols
        // of matrix A contained in 'col_indices'.  This will be the matrix to
        // be factored.  Extract the most compact submatrix possible, since
        // extraction of subcolumns may result in a submatrix with a
        // different distribution of row indices, and hence a different height.
        A.SubMatrixColsCompact(ClustData<T, MatrixType>::Asubset, 
                               col_indices, old_to_new_rows, new_to_old_rows);

        unsigned int new_height = ClustData<T, MatrixType>::Asubset.Height();

        // Create new initializer matrices for W and H by extracting the
        // appropriate subsets from the original initializers.  These
        // extractions only take a miniscule fraction of the runtime, so
        // no need to parallelize.  The NMF solver dominates the runtime.

        // extract the initializer submatrix for H
        DenseMatrix<T> Hsubset(2, subset.size());
        for (unsigned int j=0; j<subset.size(); ++j) 
        {
            // subset[j] is a column index
            Hsubset.Set(0, j, Hinit.Get(0, subset[j]));
            Hsubset.Set(1, j, Hinit.Get(1, subset[j]));
        }

        // extract the initializer submatrix for W
        DenseMatrix<T> Wsubset(new_height, 2);
        for (unsigned int i=0; i<new_height; ++i)
        {
            Wsubset.Set(i, 0, Winit.Get(new_to_old_rows[i], 0));
            Wsubset.Set(i, 1, Winit.Get(new_to_old_rows[i], 1));
        }

        // perform the rank2 factorization; make at most three attempts
        NmfStats nmf_stats;
        bool success = false;
        for (unsigned int q=0; q<3; ++q)
        {
            success = NmfSolve(ClustData<T, MatrixType>::Asubset, 
                               Wsubset, Hsubset, params.solver,
                               progress_est, params.nmf_opts, nmf_stats);
            if (success)
            {
                clust_stats.nmf_count += 1;
                if (params.nmf_opts.max_iter == nmf_stats.iteration_count)
                    clust_stats.max_count += 1;

                break;
            }

            // If here the factorization failed, probably because a singular
            // matrix was encountered.  Re-initialize W and H and try again.
            
            Random rng;
            rng.SeedFromTime();

            T* buf_w = Wsubset.Buffer();
            T* buf_h = Hsubset.Buffer();
            RandomMatrix(buf_w, Wsubset.LDim(), 
                         Wsubset.Height(), Wsubset.Width(), rng, T(0.5), T(0.5));
            RandomMatrix(buf_h, Hsubset.LDim(), 
                         Hsubset.Height(), Hsubset.Width(), rng, T(0.5), T(0.5));
        }

        // if failure after three attempts, something is really wrong
        if (!success)
            throw std::runtime_error("NMF solver failed after three attempts.");
        
        // assign cluster labels based on the values of H
        T new_label_T = static_cast<T>(new_label);
        for (int i=0; i != Hsubset.Width(); ++i) 
        {
            if (Hsubset.Get(0,i) < Hsubset.Get(1,i))
                cluster_subset.push_back(new_label_T + T(0.5));  // right child
            else
                cluster_subset.push_back(new_label_T);        // left child
        }

        // update the topic vectors for the affected rows
        MakeZeros(W_result);
        for (unsigned int i=0; i<new_height; ++i)
        {
            W_result.Set(new_to_old_rows[i], 0, Wsubset.Get(i,0));
            W_result.Set(new_to_old_rows[i], 1, Wsubset.Get(i,1));
        }

        // compute the priority score if indeed two potential clusters are generated
        std::set<T> unique_cluster_set;
        for (unsigned int i=0; i != cluster_subset.size(); ++i)
            unique_cluster_set.insert(cluster_subset[i]);

        if (unique_cluster_set.size() >= 1) 
            return compute_priority(W_parent, W_result);
        else 
        {
            return ClustHierParams<T, MatrixType, Solver>::CLUST_LEAF_PRIORITY;
        }
    }
}
