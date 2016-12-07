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

#include "dense_matrix.hpp"
#include "openmp_pragma.hpp"

namespace SolverMU
{
    const double EPSILON = 1.0e-13;
}

//-----------------------------------------------------------------------------
template <typename T>
void Update_H_MU(DenseMatrix<T>& H,
                 DenseMatrix<T>& WtA,
                 DenseMatrix<T>& WtWH)
{
    int height = H.Height();
    int width  = H.Width();

    OPENMP_PRAGMA(omp parallel for)
    for (int c=0; c<width; ++c)
    {
        for (int r=0; r<height; ++r)
        {
            T H_rc    = H.Get(r, c);
            T WtA_rc  = WtA.Get(r, c);
            T WtWH_rc = WtWH.Get(r, c);

            H_rc *= (WtA_rc / (WtWH_rc + SolverMU::EPSILON));
            H.Set(r, c, H_rc);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void Update_W_MU(DenseMatrix<T>& W,
                 DenseMatrix<T>& AHt,
                 DenseMatrix<T>& WHHt)
{
    int height = W.Height();
    int width  = W.Width();

    OPENMP_PRAGMA(omp parallel for)
    for (int c=0; c<width; ++c)
    {
        for (int r=0; r<height; ++r)
        {
            T W_rc    = W.Get(r, c);
            T AHt_rc  = AHt.Get(r, c);
            T WHHt_rc = WHHt.Get(r, c);

            W_rc *= (AHt_rc / (WHHt_rc + SolverMU::EPSILON));
            W.Set(r, c, W_rc);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class Matrix>
class Solver_Generic_MU
{
    // This class is a function object that solves for H and W using the 
    // multiplicative update rules:
    //
    //     H = H .* (W'A) ./ (W'WH + epsilon)
    //
    //     W = W .* (AH') ./ (WHH' + epsilon)
    //
    // Matlab 'dot' notation is used in this expression, indicating
    // elementwise matrix operations.  The ' character means transpose.
    //
    // A is mxn
    // W is mxk
    // H is kxn
public:

    // ------------------------------------------------------------------------
    //
    //         Init function - call prior to start of iterations
    //
    // ------------------------------------------------------------------------
    void Init (const Matrix<T>& A,
               DenseMatrix<T>& W,
               DenseMatrix<T>& H)
    {
        HHt.Resize(H.Height(), H.Height());
        WtW.Resize(W.Width(), W.Width());
        WtA.Resize(W.Width(), A.Width());
        AHt.Resize(A.Height(), H.Height());
        WtWH.Resize(W.Width(), H.Width());
        WHHt.Resize(W.Height(), H.Height());

        // compute the kxn matrix WtA =  W' * A
        Gemm(TRANSPOSE, NORMAL, T(1.0), W, A, T(0.0), WtA);
        
        // compute the kxk matrix WtW = W' * W
        Gemm(TRANSPOSE, NORMAL, T(1.0), W, W, T(0.0), WtW);
    }

    // ------------------------------------------------------------------------
    //
    //         Update function - call on each iteration
    //
    // ------------------------------------------------------------------------
    bool operator() (const Matrix<T>& A,
                     DenseMatrix<T>& W,
                     DenseMatrix<T>& H,
                     DenseMatrix<T>& gradW,
                     DenseMatrix<T>& gradH)
    {
        // solve for H

        // compute the kxn matrix WtWH = (W'W) * H = WtW * H
        Gemm(NORMAL, NORMAL, T(1.0), WtW, H, T(0.0), WtWH);

        Update_H_MU(H, WtA, WtWH);

        // solve for W

        // compute the kxk matrix HHt = H * H'
        Gemm(NORMAL, TRANSPOSE, T(1.0), H, H, T(0.0), HHt);
        
        // compute the mxk matrix AHt =  A * H'
        Gemm(NORMAL, TRANSPOSE, T(1.0), A, H, T(0.0), AHt);
        
        // compute the mxk matrix WHHt = W * HHt
        Gemm(NORMAL,  NORMAL, T(1.0), W, HHt, T(0.0), WHHt);

        Update_W_MU(W, AHt, WHHt);

        // compute gradW and gradH

        // compute the kxn matrix WtA =  W' * A
        Gemm(TRANSPOSE, NORMAL, T(1.0), W, A, T(0.0), WtA);
        
        // compute the kxk matrix WtW = W' * W
        Gemm(TRANSPOSE, NORMAL, T(1.0), W, W, T(0.0), WtW);

        // compute gradW = W*HHt - AHt
        Gemm(NORMAL, NORMAL, T(1.0), W, HHt, T(0.0), gradW);
        Axpy( T(-1.0), AHt, gradW);        
        
        // compute gradH = WtW*H - WtA
        Gemm(NORMAL, NORMAL, T(1.0), WtW, H, T(0.0), gradH);
        Axpy( T(-1.0), WtA, gradH);        

        return true;
    }

private:

    DenseMatrix<T> WtW, WtA, WtWH, HHt, AHt, WHHt;
};
