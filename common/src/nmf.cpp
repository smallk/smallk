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
#include "utils.hpp"
#include "size_check.hpp"
#include "print_norms.hpp"
#include "thread_utils.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "openmp_pragma.hpp"
#include "nmf_solve_generic.hpp"
#include "nmf_progress_estimation.hpp"

using std::cout;
using std::cerr;
using std::endl;

typedef double R;

//-----------------------------------------------------------------------------
void NmfInitialize(int argc, char * argv[]) 
{
    elem::Initialize(argc, argv);
}

//-----------------------------------------------------------------------------
NmfResult NmfIsInitialized()
{
    return elem::Initialized() ? NmfResult::INITIALIZED : 
                                 NmfResult::NOTINITIALIZED;
}

//-----------------------------------------------------------------------------
void NmfFinalize()
{
    elem::Finalize();
}

//-----------------------------------------------------------------------------
NmfResult RunNmf(const NmfOptions& opts,
                 const DenseMatrix<R>& A,
                 DenseMatrix<R>& W,
                 DenseMatrix<R>& H,
                 NmfStats& stats)
{
    bool success = false;

    ProgEstGeneric<R, DenseMatrix>* progress_estimator = 
        ProgEstGeneric<R, DenseMatrix>::Create(opts.algorithm, opts.prog_est_algorithm);

    if (nullptr == progress_estimator)
        throw std::runtime_error("invalid progress estimation algorithm");

    if (NmfAlgorithm::MU == opts.algorithm)
    {
        // multiplicative updating
        Solver_Generic_MU<R, DenseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else if (NmfAlgorithm::HALS == opts.algorithm)
    {
        // HALS
        Solver_Generic_HALS_Da<R, DenseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else if (NmfAlgorithm::RANK2 == opts.algorithm)
    {
        // Rank 2
        assert(2 == opts.k);
        if (2 != opts.k)
        {
            delete progress_estimator;
            throw std::runtime_error("rank2 algorithm requires k == 2");
        }

        Solver_Generic_Rank2<R, DenseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else if (NmfAlgorithm::BPP == opts.algorithm)
    {
        // BPP
        Solver_Generic_BPP<R, DenseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else
    {
        delete progress_estimator;
        throw std::runtime_error("unknown NMF algorithm");
    }

    if (success && opts.verbose)
        PrintNorms(A, W, H);

    delete progress_estimator;
    return success ? NmfResult::OK : NmfResult::FAILURE;
}

//-----------------------------------------------------------------------------
NmfResult RunNmf(const NmfOptions& opts,
                 const SparseMatrix<R>& A,
                 DenseMatrix<R>& W,
                 DenseMatrix<R>& H,
                 NmfStats& stats)
{
    bool success = false;

    ProgEstGeneric<R, SparseMatrix>* progress_estimator = 
        ProgEstGeneric<R, SparseMatrix>::Create(opts.algorithm, opts.prog_est_algorithm);

    if (nullptr == progress_estimator)
        throw std::runtime_error("invalid progress estimation algorithm");

    if (NmfAlgorithm::MU == opts.algorithm)
    {
        // multiplicative updating
        Solver_Generic_MU<R, SparseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else if (NmfAlgorithm::HALS == opts.algorithm)
    {
        // HALS
        Solver_Generic_HALS_Da<R, SparseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else if (NmfAlgorithm::RANK2 == opts.algorithm)
    {
        // Rank 2
        assert(2 == opts.k);
        if (2 != opts.k)
        {
            delete progress_estimator;
            throw std::runtime_error("rank2 algorithm requires k == 2");
        }

        Solver_Generic_Rank2<R, SparseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else if (NmfAlgorithm::BPP == opts.algorithm)
    {
        // BPP
        Solver_Generic_BPP<R, SparseMatrix> solver;
        success = NmfSolve<R>(A, W, H, solver, progress_estimator, opts, stats);
    }
    else
    {
        delete progress_estimator;
        throw std::runtime_error("unknown NMF algorithm");
    }

    if (success && opts.verbose)
        PrintNorms(A, W, H);

    delete progress_estimator;
    return success ? NmfResult::OK : NmfResult::FAILURE;
}

//-----------------------------------------------------------------------------
NmfResult Nmf(const NmfOptions& options,
              double *buf_a, const int ldim_a,
              double *buf_w, const int ldim_w,
              double *buf_h, const int ldim_h,
              NmfStats& stats)
{
    if (!elem::Initialized())
    {
        cerr << "nmflib error: nmf_initialize() must be called prior to "
             << "any factorization routine\n" << endl;
        return NmfResult::NOTINITIALIZED;
    }

    // check the params in case the user did not
    if (!IsValid(options))
        return NmfResult::BAD_PARAM;

    int m = options.height;
    int n = options.width;
    int k = options.k;

    // check W matrix size
    uint64_t required_size = static_cast<uint64_t>(m);
    required_size *= k;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "W matrix size too large" << endl;
        return NmfResult::SIZE_TOO_LARGE;
    }

    // check H matrix size
    required_size = static_cast<uint64_t>(n);
    required_size *= k;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "H matrix size too large" << endl;
        return NmfResult::SIZE_TOO_LARGE;
    }

    // ensure that leading dims are large enough
    assert(ldim_w >= m);
    if (ldim_w < m)
        throw std::logic_error("nmflib error: leading dimension of W return buffer too small");

    assert(ldim_h >= k);
    if (ldim_h < k)
        throw std::logic_error("nmflib error: leading dimension of H return buffer too small");

    SetMaxThreadCount(options.max_threads);

    // create Elemental matrices using the provided buffers
    DenseMatrix<R> A(m, n, buf_a, ldim_a);
    DenseMatrix<R> W(m, k, buf_w, ldim_w);
    DenseMatrix<R> H(k, n, buf_h, ldim_h);

    return RunNmf(options, A, W, H, stats);
}

//-----------------------------------------------------------------------------
NmfResult NmfSparse(const NmfOptions& options,
                    const unsigned int height,
                    const unsigned int width,
                    const unsigned int nz,
                    const unsigned int* col_offsets,
                    const unsigned int* row_indices,
                    const double* data,              
                    double *buf_w, const int ldim_w,
                    double *buf_h, const int ldim_h,
                    NmfStats& stats)
{

    if (!elem::Initialized())
    {
        cerr << "nmflib error: nmf_initialize() must be called prior to "
             << "any factorization routine\n" << endl;
        return NmfResult::NOTINITIALIZED;
    }

    // check the params in case the user did not
    if (!IsValid(options))
        return NmfResult::BAD_PARAM;

    int m = options.height;
    int n = options.width;
    int k = options.k;

    // check W matrix size
    uint64_t required_size = static_cast<uint64_t>(m);
    required_size *= k;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "W matrix size too large" << endl;
        return NmfResult::SIZE_TOO_LARGE;
    }

    // check H matrix size
    required_size = static_cast<uint64_t>(n);
    required_size *= k;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "H matrix size too large" << endl;
        return NmfResult::SIZE_TOO_LARGE;
    }

    // ensure that leading dims are large enough
    assert(ldim_w >= m);
    if (ldim_w < m)
        throw std::logic_error("nmflib error: leading dimension of W return buffer too small");

    assert(ldim_h >= k);
    if (ldim_h < k)
        throw std::logic_error("nmflib error: leading dimension of H return buffer too small");

    SetMaxThreadCount(options.max_threads);

    SparseMatrix<double> A(height, width, nz, col_offsets, row_indices, data);

    // create Elemental matrices using the provided buffers
    DenseMatrix<R> W(m, k, buf_w, ldim_w);
    DenseMatrix<R> H(k, n, buf_h, ldim_h);

    return RunNmf(options, A, W, H, stats);
}
