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

// Given the topic vectors that result from hierarchical clustering, 
// initialize a W matrix with these vectors and solve the NNLS problem 
// for H, thus generating a flat clustering result.

#include <iostream>
#include <stdexcept>
#include "nnls.hpp"
#include "clust.hpp"
#include "random.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "matrix_generator.hpp"

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType>
bool ClustFlat(const ClustOptions& opts,
               Tree<T>& tree,
               Random& rng,
               const MatrixType<T>& A,
               T* buf_w,
               T* buf_h)
{
    const unsigned int MAX_ATTEMPTS = 3;
    
    int m = opts.nmf_opts.height;
    int n = opts.nmf_opts.width;
    int k = opts.num_clusters;
    
    DenseMatrix<T> W(m, k, buf_w, m);
    
    // load matrix W with the topic vectors generated from hierclust
    if (!tree.FlatclustInitW(W, m, k))
        return false;

    bool ok = false;
    for (unsigned int q=0; q<MAX_ATTEMPTS; ++q)
    {
        // randomly initialize matrix H
        RandomMatrix(buf_h, k, k, n, rng, T(0.5), T(0.5));
        DenseMatrix<T> H(k, n, buf_h, k);
        
        ok = NnlsHals(A, W, H, 
                      opts.nmf_opts.tol, 
                      opts.verbose,
                      opts.nmf_opts.max_iter);
        if (ok)
            break;
    }

    if (!ok)
    {
        std::cout << "Flatclust NNLS solver failed after "
                  << MAX_ATTEMPTS << " attempts." << std::endl;
    }
    
    return ok;
}
