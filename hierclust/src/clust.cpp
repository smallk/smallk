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

#include <cassert>
#include <iostream>
#include <stdexcept>
#include "nmf.hpp"
#include "nnls.hpp"
#include "tree.hpp"
#include "clust.hpp"
#include "utils.hpp"
#include "size_check.hpp"
#include "thread_utils.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "nmf_solver_rank2.hpp"
#include "clust_flat_generic.hpp"
#include "clust_hier_generic.hpp"
#include "nmf_progress_estimation.hpp"

using std::cout;
using std::cerr;
using std::endl;

typedef double R;

//-----------------------------------------------------------------------------
Result RunClust(const ClustOptions& opts,
                DenseMatrix<R>& A,
                R* buf_w, R* buf_h,
                Tree<R>& tree,
                ClustStats& stats,
                Random& rng)
{
    Solver_Generic_Rank2<R, DenseMatrix> solver;

    ProgEstGeneric<R, DenseMatrix>* progress_estimator = 
        ProgEstGeneric<R, DenseMatrix>::Create(opts.nmf_opts.algorithm, 
                                               opts.nmf_opts.prog_est_algorithm);

    if (nullptr == progress_estimator)
        throw std::runtime_error("invalid progress estimation algorithm");

    bool ok = ClustHier<R>(A, solver, progress_estimator, opts, tree, stats, rng);

    if (!ok)
        return Result::FAILURE;

    if (opts.flat)
    {
        ok = ClustFlat<R>(opts, tree, rng, A, buf_w, buf_h);
        if (!ok)
        {
            cerr << "Flat clustering failed." << endl;
            return Result::FLATCLUST_FAILURE;
        }
    }

    return Result::OK;
}

//-----------------------------------------------------------------------------
Result RunClust(const ClustOptions& opts,
                const SparseMatrix<R>& A,
                R* buf_w, R* buf_h,
                Tree<R>& tree,
                ClustStats& stats,
                Random& rng)
{
    Solver_Generic_Rank2<R, SparseMatrix> solver;

    ProgEstGeneric<R, SparseMatrix>* progress_estimator = 
        ProgEstGeneric<R, SparseMatrix>::Create(opts.nmf_opts.algorithm, 
                                                opts.nmf_opts.prog_est_algorithm);

    if (nullptr == progress_estimator)
        throw std::runtime_error("invalid progress estimation algorithm");

    bool ok = ClustHier<R>(A, solver, progress_estimator, opts, tree, stats, rng);
    if (!ok)
        return Result::FAILURE;

    if (opts.flat)
    {
        ok = ClustFlat<R>(opts, tree, rng, A, buf_w, buf_h);
        if (!ok)
        {
            cerr << "Flat clustering failed." << endl;
            return Result::FLATCLUST_FAILURE;
        }
    }
    
    return Result::OK;
}

//-----------------------------------------------------------------------------
Result Clust(const ClustOptions& options,
             R *buf_a, const int ldim_a,
             R* buf_w, 
             R* buf_h,
             Tree<R>& tree,
             ClustStats& stats,
             Random& rng)
{
    if (!EL::Initialized())
    {
        cerr << "clustlib error: nmf_initialize() must be called prior to "
        << "any clustering routine\n" << endl;
        return Result::NOTINITIALIZED;
    }
   
    // check the params in case the user did not
    if (!IsValid(options))
        return Result::BAD_PARAM;
    
    int m = options.nmf_opts.height;
    int n = options.nmf_opts.width;

    if (ldim_a < m)
        throw std::logic_error("invalid leading dimension for input matrix");
    
    // check W matrix size
    uint64_t required_size = static_cast<uint64_t>(m);
    required_size *= 2u;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "W matrix size too large" << endl;
        return Result::SIZE_TOO_LARGE;
    }

    // check H matrix
    required_size = static_cast<uint64_t>(n);
    required_size *= 2u;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "H matrix size too large" << endl;
        return Result::SIZE_TOO_LARGE;
    }
    
    SetMaxThreadCount(options.nmf_opts.max_threads);

    // create the dense input matrix
    DenseMatrix<R> A(m, n, buf_a, ldim_a);

    return RunClust(options, A, buf_w, buf_h, tree, stats, rng);
}

//-----------------------------------------------------------------------------
Result ClustSparse(const ClustOptions& options,
                   const SparseMatrix<R>& A,
                   R* buf_w, 
                   R* buf_h,
                   Tree<R>& tree,
                   ClustStats& stats,
                   Random& rng)
{
    if (!EL::Initialized())
    {
        cerr << "clustlib error: nmf_initialize() must be called prior to "
        << "any clustering routine\n" << endl;
        return Result::NOTINITIALIZED;
    }
    
    // check the params in case the user did not
    if (!IsValid(options))
        return Result::BAD_PARAM;
    
    int m = options.nmf_opts.height;
    int n = options.nmf_opts.width;

    // check W matrix size
    uint64_t required_size = static_cast<uint64_t>(m);
    required_size *= 2u;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "W matrix size too large" << endl;
        return Result::SIZE_TOO_LARGE;
    }

    // check H matrix
    required_size = static_cast<uint64_t>(n);
    required_size *= 2u;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "H matrix size too large" << endl;
        return Result::SIZE_TOO_LARGE;
    }

    SetMaxThreadCount(options.nmf_opts.max_threads);

    return RunClust(options, A, buf_w, buf_h, tree, stats, rng);
}
