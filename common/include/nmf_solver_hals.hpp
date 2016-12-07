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

#include <limits>
#include <algorithm>
#include <iostream>
#include "dense_matrix.hpp"
#include "openmp_pragma.hpp"
#include "normalize.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void UpdateH_Hals(DenseMatrix<T>& H,
                  DenseMatrix<T>& WtWH_r,
                  DenseMatrix<T>& WtW,       // kxk matrix W' * W
                  const DenseMatrix<T>& WtA) // kxN matrix W' * A
{
    // this is a view, not an allocated matrix
    DenseMatrix<T> WtW_r;

    // for each row r of H
    for (int r=0; r<H.Height(); ++r)
    {
        // create a view of row r of matrix W'W
        View(WtW_r, WtW, r, 0, 1, WtW.Width());

        // compute the 1xn vector WtWH_r = W'W_r * H == H' * (W'W_r)'
        // Elemental interprets the vector on the rhs as a column vector
        // independently of how it is stored, so this becomes
        // WtWH_r = H' * (W'W_r)
        Gemv(TRANSPOSE, T(1.0), H, WtW_r, T(0.0), WtWH_r);

        // diagonal element of WtW
        T WtW_rr = WtW.Get(r, r);

        for (int c=0; c<H.Width(); ++c)
        {
            T h    = H.Get(r, c);
            T wta  = WtA.Get(r, c);
            T wtwh = WtWH_r.Get(0, c);

            h = h + (wta - wtwh) / WtW_rr;
            if (std::isnan(h) || (h < T(0)))
                h = T(0.0);

            H.Set(r, c, h);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void UpdateW_Hals(DenseMatrix<T>& W,
                  DenseMatrix<T>& WHHt_c,
                  DenseMatrix<T>& HHt,        // kxk matrix H * H'
                  const DenseMatrix<T>& AHt)  // mxk matrix A * H'
{
    // these are views
    DenseMatrix<T> HHt_c, W_c;

    int zero_col_count = 0;
    for (int c=0; c<W.Width(); ++c)
    {
        // create a view of column c of matrix HHt
        View(HHt_c, HHt, 0, c, HHt.Height(), 1);

        // compute the vector WHHt_c = W * HHt_c
        Gemv(NORMAL, T(1.0), W, HHt_c, T(0.0), WHHt_c);

        // diagonal element of HHt
        T HHt_cc = HHt.Get(c, c);

        int num_zeros = 0;
        for (int r=0; r<W.Height(); ++r)
        {
            T w    = W.Get(r, c);
            T aht  = AHt.Get(r, c);
            T whht = WHHt_c.Get(r, 0);

            w = w + (aht - whht) / HHt_cc;
            if (std::isnan(w) || (w < T(0.0)))
            {
                w = T(0.0);
                ++num_zeros;
            }

            W.Set(r, c, w);
        }

        if (W.Height() == num_zeros)
        {
            // column contains all zeros
            ++zero_col_count;
            for (int r=0; r<W.Height(); ++r)
                W.Set(r, c, std::numeric_limits<T>::epsilon());
        }

        // view of column c of matrix W
        View(W_c, W, 0, c, W.Height(), 1);
        T norm = Norm(W_c, FROBENIUS_NORM);
        T inv_norm = T(1.0) / norm;
        Scal(inv_norm, W_c);
    }
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class Matrix >
class Solver_Generic_HALS_Da
{
    // This class is a function object that solves for H and W using the 
    // HALS update rules (expressed Matlab notation):
    //
    //     H(r,:) = max( H(r,:) + W'A(r,:) - W'W(r,:)*H, epsilon)
    //
    //     W(:,c) = max( W(:,c) * HH'(c,c) + AH'(:,c) - W*HH'(:,c), epsilon)
    //     normalize the columns of W
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
    void Init(const Matrix<T>& A,
              DenseMatrix<T>& W,
              DenseMatrix<T>& H)
    {
        WtWH_r.Resize(1, H.Width());
        WHHt_c.Resize(W.Height(), 1);

        HHt.Resize(H.Height(), H.Height());
        WtW.Resize(W.Width(),  W.Width());
        AHt.Resize(A.Height(), H.Height());
        WtA.Resize(W.Width(),  A.Width());

        // compute the kxk matrix HHt = H * H'
        Gemm(NORMAL, TRANSPOSE, T(1.0), H, H, T(0.0), HHt);
        
        // compute the mxk matrix AHt =  A * H'
        Gemm(NORMAL, TRANSPOSE, T(1.0), A, H, T(0.0), AHt);
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
        // copute the new W matrix
        UpdateW_Hals(W, WHHt_c, HHt, AHt);

        // compute the kxk matrix WtW = W' * W
        Gemm(TRANSPOSE, NORMAL, T(1.0), W, W, T(0.0), WtW);
        
        // compute the kxn matrix WtA =  W' * A
        Gemm(TRANSPOSE, NORMAL, T(1.0), W, A, T(0.0), WtA);

        // compute the new H matrix
        UpdateH_Hals(H, WtWH_r, WtW, WtA);

        // compute gradH = WtW * H - WtA
        Gemm(NORMAL, NORMAL, T(1.0), WtW, H, T(0.0), gradH);
        Axpy( T(-1.0), WtA, gradH);

        // compute the kxk matrix HHt = H * H'
        Gemm(NORMAL, TRANSPOSE, T(1.0), H, H, T(0.0), HHt);

        // compute the mxk matrix AHt =  A * H'
        Gemm(NORMAL, TRANSPOSE, T(1.0), A, H, T(0.0), AHt);

        // compute gradW = W * HHt - AHt
        Gemm(NORMAL, NORMAL, T(1.0), W, HHt, T(0.0), gradW);
        Axpy( T(-1.0), AHt, gradW);        

        return true;
    }

private:

    DenseMatrix<T> WtWH_r;  // rth row of kxn matrix (W'W)H
    DenseMatrix<T> WHHt_c;  // cth col of mxk matrix W(HH')

    DenseMatrix<T> HHt, AHt, WtW, WtA;
};

