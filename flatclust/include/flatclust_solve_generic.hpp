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

#include <iostream>
#include <cassert>
#include <thread>
#include "timer.hpp"
#include "normalize.hpp"
#include "dense_matrix.hpp"
#include "nmf_solver_mu.hpp"
#include "sliding_window.hpp"
#include "nmf_solver_hals.hpp"
#include "nmf_solver_rank2.hpp"
#include "nmf_progress_estimation.hpp"

namespace FLATCLUST
{
    const size_t WINDOW_SIZE = 5;
    const double WINDOW_EPS  = 1.0e-8;
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
bool FlatClustSolve(const MatrixType<T>& A, // m x n, sparse or dense
                    DenseMatrix<T>& W,      // m x k
                    DenseMatrix<T>& H,      // k x n
                    Solver<T, MatrixType>& solver,
                    ProgressEst<T, MatrixType>* progress_est,
                    const NmfOptions& opts,
                    NmfStats& stats)
{
    int count;
    Timer timer;
    uint64_t total_us = 0;

    int m = A.Height();
    int n = A.Width();
    int k = W.Width();
  
    int min_iter    = opts.min_iter;
    int max_iter    = opts.max_iter;
    T tolerance     = opts.tol;
    int tolcount    = opts.tolcount;
    bool verbose    = opts.verbose;
    bool normalize  = opts.normalize;
    int max_threads = opts.max_threads;
  
    // kxk matrices W'W and HH'
    DenseMatrix<T> WtW(k, k), HHt(k, k);

    // kxn matrices W'A, (W'W)*H
    DenseMatrix<T> WtA(k, n);
    DenseMatrix<T> WtWH(k, n);
    DenseMatrix<T> gradH(k, n);

    // mxk matrices AH', WHH'
    DenseMatrix<T> AHt(m, k);
    DenseMatrix<T> WHHt(m, k);
    DenseMatrix<T> gradW(m, k);

    // kx1 matrix of scale factors
    DenseMatrix<T> Norms(k, 1);

    // sliding window detects level-off behavior in the projected gradient
    SlidingWindow<T> window(FLATCLUST::WINDOW_SIZE);

    timer.Start();

    // initialize the solver
    solver.Init(A, W, H, WtA, WtW, max_threads);

    // initialize the progress estimator
    progress_est->Init(A, W, H, gradW, gradH, max_threads);

    auto prev_precision = std::cout.precision(15);

    bool success = false;
    int success_count = 0;
    for (count=0; count<max_iter; ++count)
    {
        // run the solver for the current iteration
        if (!solver(A, W, H, WtW, WtA, WtWH, HHt, AHt, WHHt, Norms, max_threads))
        {
            std::cerr << "\tNMF solver failure on iteration " 
                      << count+1 << std::endl;
            break;
        }
      
        // perform 'min_iter' iterations before progress estimation
        if (count < min_iter)
        {
            // compute the PG after the initial iteration for PG_RATIO
            if (0 == count)
                progress_est->Update(count, W, H, WtW, HHt, WtA, AHt, Norms, gradW, gradH);

            if (verbose)
                std::cout << count+1 << ":\tprogress metric: \t(min_iter)" 
                          << std::endl;

            continue;
        }

        // check progress
        T progress_metric = progress_est->Update(count, W, H, WtW, HHt, WtA, AHt, Norms, gradW, gradH);
        if (verbose)
            ReportProgress(count+1, progress_metric);

        // update the window
        window.Update(progress_metric);
        
        if (progress_metric <= tolerance)
        {
            ++success_count;
            if (success_count >= tolcount)
            {
                success = true;
                if (verbose)
                {
                    std::cout << "\nSolution converged after " << count+1 
                              << " iterations." << std::endl << std::endl;
                }

                break;
            }
        }
        else
        {
            // metric exceeded tolerance, so reset success count
            success_count = 0;

            // check window for level-off behavior
            if (window.AllBounded(FLATCLUST::WINDOW_EPS))
            {
                //ComputeAssignments(assignments, H.LockedBuffer(), H.LDim(), k, n);
                //TopTerms(opts.clust_opts.maxterms, W.LockedBuffer(), W.LDim(), m, k, term_indices);
                // uint64_t h1 = SpookyHash(assignments);
                // uint64_t h2 = SpookyHasn(term_indices);
                // store in a struct and store in a new sliding window
                // declare success when all hashes in the window are equal

                success = true;
                if (verbose)
                {
                    std::cout << "\nProgress metric values are stationary "
                              << "after " << count+1 << " iterations." 
                              << std::endl << std::endl;
                }

                break;
            }
        }
    }

    // optional normalization occurs only after iterations have completed
    if (normalize)
        NormalizeAndScale(W, H);

    timer.Stop();
    total_us += timer.ReportMicroseconds();

    if (verbose)
    {
        if (!success)
        {
            if (max_iter == count)
            {
                std::cout << "\nReached iteration limit." 
                          << std::endl << std::endl;
            }
            else
            {
                // numerical error
                std::cout << "\nSolver failure." << std::endl;
            }
        }
    }

    std::cout.precision(prev_precision);

    // reaching the iteration limit is not an error - the user may
    // want to run with a limited number of iterations, for instance
    if (!success && (max_iter == count))
        success = true;

    stats.elapsed_us = total_us;
    stats.iteration_count = count;

    return success;
}

