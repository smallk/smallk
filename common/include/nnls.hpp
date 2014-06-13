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

// Solve the nonnegative least squares problem ||WH - A||^2 for H.

#include <iostream>
#include <stdexcept>
#include "normalize.hpp"
#include "bit_matrix.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "bit_matrix_ops.hpp"
#include "nmf_solver_bpp.hpp"
#include "nmf_solver_hals.hpp"
#include "projected_gradient.hpp"
#include "nmf_progress_estimation.hpp"

void UpdatePassiveSet(BitMatrix& passive_set,
                      const int PBAR,
                      const unsigned int n,
                      const BitMatrix& not_opt_cols,
                      const BitMatrix& nonopt_set,
                      const BitMatrix& infeas_set,
                      const std::vector<int>& not_good,
                      std::vector<int>& P, 
                      std::vector<int>& Ninf);

//-----------------------------------------------------------------------------
template <typename T>
void BppUpdateSets(BitMatrix& nonopt_set,
                   BitMatrix& infeas_set,
                   const BitMatrix& not_opt_mask,
                   const DenseMatrix<T>& X,
                   const DenseMatrix<T>& Y,
                   const BitMatrix& passive_set)
{
    // This function performs the equivalent of these operations:
    //
    //     nonopt_set = not_opt_mask & (Y < T(0)) & ~passive_set;
    //     infeas_set = not_opt_mask & (X < T(0)) & passive_set;

    const unsigned int height = not_opt_mask.Height();
    const unsigned int width  = not_opt_mask.Width();

    if ( (static_cast<unsigned int>(X.Height()) != height) || 
         (static_cast<unsigned int>(Y.Height()) != height) || 
         (passive_set.Height() != height))
        throw std::logic_error("BppUpdateSets: height mismatch");

    if ( (static_cast<unsigned int>(X.Width()) != width) || 
         (static_cast<unsigned int>(Y.Width()) != width) || 
         (passive_set.Width() != width))
        throw std::logic_error("BppUpdateSets: width mismatch");

    nonopt_set.Resize(height, width);
    infeas_set.Resize(height, width);

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;
    const unsigned int MASK   = nonopt_set.Mask();
    assert(infeas_set.Mask() == MASK);

    unsigned int* buf_r = nonopt_set.Buffer();
    const unsigned int ldim_r = nonopt_set.LDim();

    unsigned int* buf_i = infeas_set.Buffer();
    const unsigned int ldim_i = infeas_set.LDim();

    const unsigned int* buf_m = not_opt_mask.LockedBuffer();
    const unsigned int ldim_m = not_opt_mask.LDim();

    const T* buf_x = X.LockedBuffer();
    const unsigned int ldim_x = X.LDim();

    const T* buf_y = Y.LockedBuffer();
    const unsigned int ldim_y = Y.LDim();

    const unsigned int* buf_p = passive_set.LockedBuffer();
    const unsigned int ldim_p = passive_set.LDim();

    const unsigned int full_wds = height / BITS;
    const unsigned int extra    = height - BITS*full_wds;

    assert(ldim_r >= ldim_m);
    assert(ldim_r >= ldim_p);

    OPENMP_PRAGMA(omp parallel for)
    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int offset_r = c*ldim_r;
        unsigned int offset_i = c*ldim_i;
        unsigned int offset_m = c*ldim_m;
        unsigned int offset_x = c*ldim_x;
        unsigned int offset_y = c*ldim_y;
        unsigned int offset_p = c*ldim_p;

        unsigned int r_wd = 0, r=0;
        for (; r_wd<full_wds; ++r_wd)
        {
            unsigned int x_wd = 0, y_wd = 0;
            for (unsigned int q=0; q<BITS; ++q, ++r)
            {
                if (buf_x[offset_x + r] < T(0))
                    x_wd |= (1 << q);
                if (buf_y[offset_y + r] < T(0))
                    y_wd |= (1 << q);
            }
            
            buf_r[offset_r + r_wd] = buf_m[offset_m + r_wd] & y_wd & ~buf_p[offset_p + r_wd];
            buf_i[offset_i + r_wd] = buf_m[offset_m + r_wd] & x_wd & buf_p[offset_p + r_wd];
        }

        if (extra > 0)
        {
            unsigned int x_wd = 0, y_wd = 0;
            for (unsigned int q=0; q<extra; ++q, ++r)
            {
                if (buf_x[offset_x + r] < T(0))
                    x_wd |= (1 << q);
                if (buf_y[offset_y + r] < T(0))
                    y_wd |= (1 << q);
            }

           buf_r[offset_r + r_wd] = MASK & buf_m[offset_m + r_wd] & y_wd & ~buf_p[offset_p + r_wd];
           buf_i[offset_i + r_wd] = MASK & buf_m[offset_m + r_wd] & x_wd & buf_p[offset_p + r_wd];
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
bool NnlsBlockpivot(const DenseMatrix<T>& LHS,
                    const DenseMatrix<T>& RHS,
                    DenseMatrix<T>& X, // input as xinit
                    DenseMatrix<T>& Y) // gradX
{
    // Solve (LHS)*X = RHS for X by block principal pivoting. Matrix LHS
    // is assumed to be symmetric positive definite.

    const int PBAR = 3;
    const unsigned int WIDTH  = RHS.Width();
    const unsigned int HEIGHT = RHS.Height();
    const unsigned int MAX_ITER = HEIGHT*5;

    BitMatrix passive_set = (X > T(0));

    std::vector<unsigned int> tmp_indices(WIDTH);
    for (unsigned int i=0; i<WIDTH; ++i)
        tmp_indices[i] = i;

    MakeZeros(X);
    if (!BppSolveNormalEqNoGroup(tmp_indices, passive_set, LHS, RHS, X))
        return false;

    // Y = LHS * X - RHS
    Gemm(NORMAL, NORMAL, T(1), LHS, X, T(0), Y);
    Axpy( T(-1), RHS, Y);

    std::vector<int> P(WIDTH, PBAR), Ninf(WIDTH, HEIGHT+1);

    BitMatrix nonopt_set = (Y < T(0)) & ~passive_set;
    BitMatrix infeas_set = (X < T(0)) & passive_set;
    std::vector<int> col_sums(WIDTH);
    std::vector<int> not_good(WIDTH);

    nonopt_set.SumColumns(not_good);
    infeas_set.SumColumns(col_sums);
    not_good += col_sums;
    BitMatrix not_opt_cols = (not_good > 0);
    BitMatrix not_opt_mask;

    std::vector<unsigned int> non_opt_col_indices(WIDTH);
    not_opt_cols.Find(non_opt_col_indices);

    DenseMatrix<double> RHSsub(HEIGHT, WIDTH);
    DenseMatrix<double> Xsub(HEIGHT, WIDTH);
    DenseMatrix<double> Ysub(HEIGHT, WIDTH);

    unsigned int iter = 0;
    while (!non_opt_col_indices.empty())
    {
        // exit if not getting anywhere
        if (iter >= MAX_ITER)
            return false;

        UpdatePassiveSet(passive_set, PBAR, HEIGHT,
                         not_opt_cols, nonopt_set, infeas_set,
                         not_good, P, Ninf);

        // equivalent of repmat(NotOptCols, HEIGHT, 1)
        not_opt_mask = MatrixFromColumnMask(not_opt_cols, HEIGHT);

        // Setup for the normal equation solver by extracting submatrices
        // from RHS and X.  The normal equation solver will extract further
        // subproblems from RHSsub and Xsub and write all updated values
        // back into RHSsub and Xsub.
        RHS.SubmatrixFromCols(RHSsub, non_opt_col_indices);
        X.SubmatrixFromCols(Xsub, non_opt_col_indices);

        if (!BppSolveNormalEqNoGroup(non_opt_col_indices, passive_set, LHS, RHSsub, Xsub))
            return false;

        ZeroizeSmallValues(Xsub, 1.0e-12);

        // compute Ysub = LHS * Xsub - RHSsub
        Ysub.Resize(RHSsub.Height(), RHSsub.Width());
        Gemm(NORMAL, NORMAL, T(1), LHS, Xsub, T(0), Ysub);
        Axpy( T(-1), RHSsub, Ysub);

        // update Y and X using the new values in Ysub and Xsub
        OverwriteCols(Y, Ysub, non_opt_col_indices, non_opt_col_indices.size());
        OverwriteCols(X, Xsub, non_opt_col_indices, non_opt_col_indices.size());

        ZeroizeSmallValues(X, 1.0e-12);
        ZeroizeSmallValues(Y, 1.0e-12);

        // Check optimality - BppUpdateSets does the equivalent of the next two lines.
        // nonopt_set = not_opt_mask & (Y < T(0)) & ~passive_set;
        // infeas_set = not_opt_mask & (X < T(0)) & passive_set;
        BppUpdateSets(nonopt_set, infeas_set, not_opt_mask, X, Y, passive_set);

        nonopt_set.SumColumns(not_good);
        infeas_set.SumColumns(col_sums);
        not_good += col_sums;
        not_opt_cols = (not_good > 0);
        not_opt_cols.Find(non_opt_col_indices);

        ++iter;
    }

    return true;
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType>
bool NnlsHals(const MatrixType<T>& A,
              DenseMatrix<T>& W, 
              DenseMatrix<T>& H,
              const T tol,
              const bool verbose,
              const unsigned int max_iter)
{
    unsigned int n = A.Width();
    unsigned int k = W.Width();

    if (static_cast<unsigned int>(W.Height()) != static_cast<unsigned int>(A.Height()))
        throw std::logic_error("NnlsHals: W and A must have identical height");
    if (static_cast<unsigned int>(H.Width()) != static_cast<unsigned int>(A.Width()))
        throw std::logic_error("NnlsHals: H and A must have identical width");
    if (H.Height() != W.Width())
        throw std::logic_error("NnlsHals: non-conformant W and H");

    DenseMatrix<T> WtW(k, k), WtA(k, n), WtWH_r(1, n), gradH(k, n);

    if (verbose)
        std::cout << "\nRunning NNLS solver..." << std::endl;
    
    // compute W'W and W'A for the normal equations
    Gemm(TRANSPOSE, NORMAL, T(1.0), W, W, T(0.0), WtW);
    Gemm(TRANSPOSE, NORMAL, T(1.0), W, A, T(0.0), WtA);

    bool success = false;

    T pg0 = T(0), pg;
    for (unsigned int i=0; i<max_iter; ++i)
    {
        // compute the new matrix H
        UpdateH_Hals(H, WtWH_r, WtW, WtA);
        
        // compute gradH = WtW*H - WtA
        Gemm(NORMAL, NORMAL, T(1.0), WtW, H, T(0.0), gradH);
        Axpy( T(-1.0), WtA, gradH);

        // compute progress metric
        if (0 == i)
        {
            pg0 = ProjectedGradientNorm(gradH, H);
            if (verbose)
                ReportProgress(i+1, T(1.0));
            continue;
        }
        else
        {
            pg = ProjectedGradientNorm(gradH, H);
        }

        if (verbose)
            ReportProgress(i+1, pg/pg0);
        
        // check progress vs. desired tolerance
        if (pg < tol * pg0)
        {
            success = true;
            NormalizeAndScale<T>(W, H);
            break;
        }
    }

    if (!success)
        std::cerr << "NNLS solver reached iteration limit." << std::endl;
    
    return success;
}
