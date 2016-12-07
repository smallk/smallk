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
#include "nmf_solver_bpp.hpp"
#include "nmf_solver_hals.hpp"
#include "nmf_solver_rank2.hpp"
#include "nmf_progress_estimation.hpp"

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType,
          template <typename, template <typename> class Matrix> class Solver,
          template <typename, template <typename> class Matrix> class ProgressEst>
bool NmfSolve(const MatrixType<T>& A, // m x n, sparse or dense
              DenseMatrix<T>& W,      // m x k
              DenseMatrix<T>& H,      // k x n
              Solver<T, MatrixType>& solver,
              ProgressEst<T, MatrixType>* progress_est,
              const NmfOptions& opts,
              NmfStats& stats)
{
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
  
    DenseMatrix<T> gradH(k, n);
    DenseMatrix<T> gradW(m, k);

    timer.Start();

    // do one-time init prior to beginning iterations
    solver.Init(A, W, H);
    progress_est->Init(A, W, H);

    bool success = false;
    int iter=0, success_count = 0;
    for (iter=0; iter<max_iter; ++iter)
    {
        // run the solver for the current iteration
        if (!solver(A, W, H, gradW, gradH))
        {
            timer.Stop();
            total_us += timer.ReportMicroseconds();
            stats.elapsed_us = total_us;
            stats.iteration_count = iter;
            std::cerr << "\tNMF solver failure on iteration " 
                      << iter+1 << std::endl;
            return false;
        }

        if (iter < min_iter)
        {
            // compute norm of projected gradient for initial iteration
            if (0 == iter)
                progress_est->Update(iter, W, H, gradW, gradH);

            if (verbose)
            {
                std::cout << iter+1 << ":\tprogress metric: \t(min_iter)" 
                          << std::endl;
            }

            // perform 'min_iter' iterations before progress estimation
            continue;
        }

        // check progress
        T progress_metric = progress_est->Update(iter, W, H, gradW, gradH);
        if (verbose)
            ReportProgress(iter+1, progress_metric);

        if (progress_metric <= tolerance)
        {
            ++success_count;
            if (success_count >= tolcount)
            {
                success = true;
                if (verbose)
                {
                    std::cout << "\nSolution converged after " << iter+1 
                              << " iterations." << std::endl << std::endl;
                }

                break;
            }
        }
        else
        {
            // metric exceeded tolerance, so reset success count
            success_count = 0;
        }

    }

    // optional normalization occurs only after iterations have completed
    if (normalize)
        NormalizeAndScale(W, H);

    timer.Stop();
    total_us += timer.ReportMicroseconds();

    // reaching the iteration limit is not an error - the user may
    // want to run with a limited number of iterations, for instance
    if (!success && (max_iter == iter))
        success = true;

    stats.elapsed_us = total_us;
    stats.iteration_count = iter;
    return success;
}

