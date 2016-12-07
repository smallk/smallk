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

// This file defines some auxiliary functions for
// hierarchical clustering based on rank-2 NMF.
// Should be included in:
//   clust_hier.hpp
//   clust_hier_sparse.hpp

#pragma once

#include <vector>
#include <algorithm>
#include "dense_matrix.hpp"

//-----------------------------------------------------------------------------
// Auxiliary function
template <typename T>
std::vector<int> ordered(std::vector<T> const& values) 
{
    std::vector<int> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<int>(0));
    std::sort(begin(indices), end(indices), 
              [&](int a, int b)
              {
                  return values[a] < values[b] || (values[a]==values[b] && a<b);
              });

    return indices;
}

//-----------------------------------------------------------------------------
// Auxiliary function
template <typename T>
std::vector<int> desc_ordered(std::vector<T> const& values) 
{
    std::vector<int> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<int>(0));
    std::sort(begin(indices), end(indices), 
              [&](int a, int b)
              {
                  return values[a] > values[b] || (values[a]==values[b] && a<b);
              });

    return indices;
}

//-----------------------------------------------------------------------------
// Compute the modified NDCG value
template <typename T>
T NDCG_part(std::vector<int>& ground, 
            std::vector<int>& test, 
            std::vector<T>&weight, 
            std::vector<T>& weight_part) 
{
    std::vector<int> seq_idx;
    seq_idx = ordered(ground);
    std::vector<T> temp_weight_part((int)weight_part.size());
    for (unsigned int i = 0; i < weight_part.size(); ++i) {
        temp_weight_part[i] = weight_part[seq_idx[i]];
    }

    int n = test.size();
    std::vector<T> uncum_score(n), cum_score(n);
    for (int i = 0; i < n; ++i) {
        uncum_score[i] = temp_weight_part[test[i]];
        if (i > 0) {
            uncum_score[i] = uncum_score[i] / log2(i+1);
            cum_score[i] = cum_score[i-1] + uncum_score[i];
        } else {
            cum_score[i] = uncum_score[i];
        }
    }
    std::vector<T> ideal_score = weight;
    std::sort(ideal_score.begin(), ideal_score.end(), std::greater<T>());
    
    std::vector<T> cum_ideal_score(n);
    for (int i = 0; i < n; ++i) {
        if (i > 0) {
            ideal_score[i] = ideal_score[i] / log2(i+1);
            cum_ideal_score[i] = cum_ideal_score[i-1] + ideal_score[i];
        } else {
            cum_ideal_score[i] = ideal_score[i];
        }
    }
    
    return cum_score.back() / cum_ideal_score.back();
}

//-----------------------------------------------------------------------------
// Compute the priority score for a node (represented by W_parent)
// based on the rank-2 NMF result on it (represented by W_child)
template <typename T>
T compute_priority(DenseMatrix<T>& W_parent, // m by 1
                   DenseMatrix<T>& W_child)  // m by 2
{
    int n = W_parent.Height();
    int n_child = W_child.Height();
    std::vector<T> val_parent(n), val_child1(n_child), val_child2(n_child);
    std::vector<int> idx_parent, idx_child1, idx_child2;
    
    int n_part = 0;
    for (int i = 0; i < n; ++i) {
        val_parent[i] = W_parent.Get(i,0);
        if (W_parent.Get(i,0) != 0)
            n_part++;
    }
    for (int i = 0; i < n_child; ++i) {
        val_child1[i] = W_child.Get(i,0);
        val_child2[i] = W_child.Get(i,1);
    }
    idx_parent = desc_ordered(val_parent);
    idx_child1 = desc_ordered(val_child1);
    idx_child2 = desc_ordered(val_child2);
    
    if (n_part <= 1) {
        return -3;
    } else {
        std::vector<T> weight(n);
        for (int i = n; i >= 1; --i) {
            weight[n-i] = log(i);
        }
        int first_zero = -1;
        for (unsigned int i = 0; i < idx_parent.size(); ++i) {
            if (val_parent[idx_parent[i]] == 0) {
                first_zero = i;
                break;
            }
        }
        if (first_zero > -1) {
            for (int i = first_zero; i < n; ++i) {
                weight[i] = 1;
            }
        }
        std::vector<T> weight_part(n);
        for (int i = n_part; i >= 1; --i) {
            weight_part[n_part-i] = log(i);
        }
        for (int i = n_part; i < n; ++i) {
            weight_part[i] = 0;
        }
        std::vector<int> idx1, idx2;
        idx1 = ordered(idx_child1);
        idx2 = ordered(idx_child2);
        std::vector<int> max_pos(idx1.size());
        for (unsigned int i = 0; i < idx1.size(); ++i) {
            max_pos[i] = (idx1[i] < idx2[i] ? idx2[i] : idx1[i]);
        }
        std::vector<T> discount(idx_parent.size());
        for (unsigned int i = 0; i < idx_parent.size(); ++i) {
            discount[i] = log(n - max_pos[idx_parent[i]]);
            if (discount[i] == 0)
                discount[i] = log(2);
            
            weight[i] = weight[i] / discount[i];
            weight_part[i] = weight_part[i] / discount[i];
        }
        
        return NDCG_part(idx_parent, idx_child1, weight, weight_part) * 
               NDCG_part(idx_parent, idx_child2, weight, weight_part);
    }
}

