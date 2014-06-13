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
               const std::vector<DenseMatrix<T> >& topic_vectors,
               const MatrixType<T>& A,
               T* buf_w,
               T* buf_h)
{
    int m = opts.nmf_opts.height;
    int n = opts.nmf_opts.width;
    int k = opts.num_clusters;
    
    DenseMatrix<T> W(m, k, buf_w, m);
    
    // load matrix W with the topic vectors
    DenseMatrix<T> col_w;
    for (int c=0; c<k; ++c)
    {
        // create a view of the cth column of W
        View(col_w, W, 0, c, m, 1);
        
        // copy the next topic vector into place
        Copy(topic_vectors[c], col_w);
    }
    
    // randomly initialize matrix H
    Random rng;
    rng.SeedFromTime();
    RandomMatrix(buf_h, k, k, n, rng, T(0.5), T(0.5));
    DenseMatrix<T> H(k, n, buf_h, k);

    bool ok = NnlsHals(A, W, H, 
                       opts.nmf_opts.tol, 
                       opts.verbose,
                       opts.nmf_opts.max_iter);
    return ok;
}
