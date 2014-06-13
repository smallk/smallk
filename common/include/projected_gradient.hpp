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

#include <cmath>
#include <cassert>
#include <stdexcept>
#include "enums.hpp"
#include "sparse_matrix.hpp"
#include "dense_matrix.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void GradientCommon(const DenseMatrix<T>& W, // m x k
                    const DenseMatrix<T>& H, // k x n
                    DenseMatrix<T>& gradW,   // m x k
                    DenseMatrix<T>& gradH)   // k x n
{
    int k = W.Width();

    // compute kxk matrix H*H' = HHt
    DenseMatrix<T> HHt(k, k);
    Gemm(NORMAL,
         TRANSPOSE,
         T(1.0), H, H, T(0.0), HHt);

    // compute kxk W'*W = WtW
    DenseMatrix<T> WtW(k, k);
    Gemm(TRANSPOSE,
         NORMAL,
         T(1.0), W, W, T(0.0), WtW);

    // compute mxk matrix W*HHt and store in gradW
    Gemm(NORMAL,
         NORMAL,
         T(1.0), W, HHt, T(0.0), gradW);

    // compute kxn matrix WtW*H and store in gradH
    Gemm(NORMAL, 
         NORMAL,
         T(1.0), WtW, H, T(0.0), gradH);
}

//-----------------------------------------------------------------------------
template <typename T>
T GradientNorm(const DenseMatrix<T>& gradW,
               const DenseMatrix<T>& gradH)
{
    T normW = T(0.0), normH = T(0.0);

    int height = gradW.Height();
    int width  = gradW.Width();
    for (int c=0; c < width; ++c)
    {
        for (int r=0; r < height; ++r)
        {
            T elt = gradW.Get(r, c);
            normW += elt*elt;
        }
    }

    height = gradH.Height();
    width  = gradH.Width();
    for (int c=0; c < width; ++c)
    {
        for (int r=0; r < height; ++r)
        {
            T elt = gradH.Get(r, c);
            normH += elt*elt;
        }
    }

    T gradient = sqrt(normW + normH);
    if (std::isnan(gradient))
        throw std::runtime_error("GradientNorm: NaN");

    return gradient;
}

//-----------------------------------------------------------------------------
template <typename T>
T ProjectedGradientNorm(const DenseMatrix<T>& gradM, const DenseMatrix<T>& M)
{
    // Computes the norm of the PG for a single matrix.

    const T ZERO(0);
    T normM = ZERO;

    int height = gradM.Height();
    int width  = gradM.Width();

    for (int c=0; c<width; ++c)
    {
        for (int r=0; r<height; ++r)
        {
            T gradM_rc = gradM.Get(r, c);
            if ( (gradM_rc < ZERO) || (M.Get(r, c) > ZERO))
            {
                normM += (gradM_rc * gradM_rc);
            }
        }
    }

    T projgrad = sqrt(normM);
    if (std::isnan(projgrad))
        throw std::runtime_error("ProjectedGradientNorm: NaN");

    return projgrad;
}

//-----------------------------------------------------------------------------
template <typename T>
T ProjectedGradientNorm(const DenseMatrix<T>& gradW,
                        const DenseMatrix<T>& gradH,
                        const DenseMatrix<T>& W,
                        const DenseMatrix<T>& H)
{
    // Compute the norm of the PG for W and H.  Note that the return value
    // is the square root of the sum, not the sum of the square roots.

    const T ZERO(0);
    T normW = ZERO, normH = ZERO;

    int height = gradW.Height();
    int width  = gradW.Width();

    for (int c=0; c<width; ++c)
    {
        for (int r=0; r<height; ++r)
        {
            T gradW_rc = gradW.Get(r, c);
            if ( (gradW_rc < ZERO) || (W.Get(r, c) > ZERO))
            {
                normW += (gradW_rc * gradW_rc);
            }
        }
    }
    
    height = gradH.Height();
    width  = gradH.Width();

    for (int c=0; c<width; ++c)
    {
        for (int r=0; r<height; ++r)
        {
            T gradH_rc = gradH.Get(r, c);
            if ( (gradH_rc < ZERO) || (H.Get(r, c) > ZERO))
            {
                normH += (gradH_rc * gradH_rc);
            }
        }
    }

    T projgrad = sqrt(normW + normH);
    if (std::isnan(projgrad))
        throw std::runtime_error("ProjectedGradientNorm: NaN");

    return projgrad;
}

//-----------------------------------------------------------------------------
template <typename T>
T Gradient(const DenseMatrix<T>& A, // m x n
           const DenseMatrix<T>& W, // m x k
           const DenseMatrix<T>& H, // k x n
           DenseMatrix<T>& gradW,   // m x k
           DenseMatrix<T>& gradH,   // k x n
           const unsigned int max_threads = 0)
{
    int m = A.Height();
    int n = A.Width();
    int k = W.Width();

    GradientCommon<T>(W, H, gradW, gradH);

    // compute mxk matrix A*H' and store in AHt
    DenseMatrix<T> AHt(m, k);
    Gemm(NORMAL,
         TRANSPOSE,
         T(1.0), A, H, T(0.0), AHt);

    // compute mxk matrix gradW = (-1.0)*AH' + WHHt  (-1.0)*AHt + gradW
    Axpy( T(-1.0), AHt, gradW);

    // compute kxn matrix W'*A and store in WtA
    DenseMatrix<T> WtA(k, n);
    Gemm(TRANSPOSE,
         NORMAL,
         T(1.0), W, A, T(0.0), WtA);

    // compute kxn matrix gradH = WtW*H - W'*A = (-1.0)*WtA + gradH
    Axpy( T(-1.0), WtA, gradH);

    return GradientNorm(gradW, gradH);
}

//-----------------------------------------------------------------------------
template <typename T>
T Gradient(const SparseMatrix<T>& A, // m x n
           const DenseMatrix<T>& W,  // m x k
           const DenseMatrix<T>& H,  // k x n
           DenseMatrix<T>& gradW,    // m x k
           DenseMatrix<T>& gradH,    // k x n
           const unsigned int max_threads)
{
    // gradW = W*HHt - A*H'
    // gradH = WtW*H - W'*A

    int m = A.Height();
    int n = A.Width();
    int k = W.Width();

    GradientCommon<T>(W, H, gradW, gradH);

    // compute mxk matrix A*H' and store in AHt
    DenseMatrix<T> AHt(m, k);
    Gemm(NORMAL,
         TRANSPOSE,
         T(1.0), A, H, T(0.0), AHt, max_threads);

    // compute mxk matrix gradW = (-1.0)*AH' + WHHt  (-1.0)*AHt + gradW
    Axpy( T(-1.0), AHt, gradW);

    // compute kxn matrix W'*A and store in WtA
    DenseMatrix<T> WtA(k, n);
    Gemm(TRANSPOSE,
         NORMAL,
         T(1.0), W, A, T(0.0), WtA, max_threads);

    // compute kxn matrix gradH = WtW*H - W'*A = (-1.0)*WtA + gradH
    Axpy( T(-1.0), WtA, gradH);

    return GradientValue(gradW, gradH);
}

//-----------------------------------------------------------------------------
template <typename T>
T ProjectedGradient(const DenseMatrix<T>& W, 
                    const DenseMatrix<T>& H, 
                    const DenseMatrix<T>& WtW, 
                    const DenseMatrix<T>& HHt,
                    const DenseMatrix<T>& WtA, 
                    const DenseMatrix<T>& AHt,
                    DenseMatrix<T>& gradW,
                    DenseMatrix<T>& gradH)
{
    //
    // compute gradW = W*HHt - AHt
    //

    // compute the mxk matrix W*HHt and store in gradW
    Gemm(NORMAL,
         NORMAL,
         T(1.0), W, HHt, T(0.0), gradW);

    // compute gradW = (WHHt - AHt) = (-1.0)*AHt + WHHt = (-1.0)*AHt + gradW
    Axpy( T(-1.0), AHt, gradW);

    //
    // compute gradH = WtW*H - WtA
    //

    // compute the kxn matrix WtW*H and store in gradH
    Gemm(NORMAL,
         NORMAL,
         T(1.0), WtW, H, T(0.0), gradH);

    // compute gradH = (WtWH - WtA) = (-1.0)*WtA + WtWH = (-1.0)*WtA + gradH
    Axpy( T(-1.0), WtA, gradH);

    return ProjectedGradientNorm(gradW, gradH, W, H);
}

// //-----------------------------------------------------------------------------
// template <typename T>
// T ProjectedGradientRank2(const DenseMatrix<T>& W, 
//                          const DenseMatrix<T>& H, 
//                          const DenseMatrix<T>& WtW, 
//                                DenseMatrix<T>& HHt,
//                          const DenseMatrix<T>& WtA, 
//                                DenseMatrix<T>& AHt,
//                          const DenseMatrix<T>& ScaleFactors,
//                          DenseMatrix<T>& gradW,
//                          DenseMatrix<T>& gradH)
// {
//     // The Rank2 solver does not update HH' or AH' after solving for H, 
//     // since it does so prior to solving for W.  The projected gradient,
//     // however, requires full updates for all matrices.  Thus the updated
//     // HH' and AH' matrices are computed prior to the PG calculation.

//     T s0 = ScaleFactors.Get(0, 0);
//     T s1 = ScaleFactors.Get(1, 0);
    
//     // update HH' for the rescaled H
//     T elt00 = HHt.Get(0,0);
//     T elt01 = HHt.Get(0,1);
//     T elt11 = HHt.Get(1,1);
//     HHt.Set(0, 0, elt00*s0*s0);
//     HHt.Set(0, 1, elt01*s0*s1);
//     HHt.Set(1, 0, elt01*s0*s1);  // since HH' is symmetric
//     HHt.Set(1, 1, elt11*s1*s1);
    
//     // update the mxk matrix AHt =  A * H' by scaling each column
//     DiagonalScale(RIGHT, NORMAL, ScaleFactors, AHt);
    
//     return ProjectedGradient(W, H, WtW, HHt, WtA, AHt, gradW, gradH);
// }
