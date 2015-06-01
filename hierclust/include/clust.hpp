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

#include <vector>
#include "nmf.hpp"
#include "tree.hpp"
#include "random.hpp"
#include "sparse_matrix.hpp"

typedef double R;

struct ClustStats
{
    ClustStats() : nmf_count(0),  max_count(0) {}

    // number of factorizations performed
    int nmf_count;

    // number of factorizations that reached iter limit
    int max_count;
};

struct ClustOptions
{
    NmfOptions nmf_opts;

    int maxterms;
    R unbalanced;
    int trial_allowance;
    int num_clusters;
    bool verbose;
    bool flat;
    std::string initdir;
};

bool IsValid(const ClustOptions& opts, bool validate_matrix = true);

Result Clust(const ClustOptions& options,
             R* buf_A, int ldim_A,
             R* buf_w, R* buf_h,
             Tree<R>& tree,
             ClustStats& stats,
             Random& rng);

Result ClustSparse(const ClustOptions& options,
                   const SparseMatrix<R>& A,
                   R* buf_w, R* buf_h,
                   Tree<R>& tree,
                   ClustStats& stats,
                   Random& rng);

