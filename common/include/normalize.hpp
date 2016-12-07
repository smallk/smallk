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
#include <cstdlib>
#include <cassert>
#include "dense_matrix.hpp"
#include "openmp_pragma.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void NormalizeColumns(DenseMatrix<T>& W, DenseMatrix<T>& Norms)
{
    // Normalize the columns of matrix W to unit L2 norm and return the norms.
    // For matrices stored in column-major order this code should be very fast.

    const T eps = std::numeric_limits<T>::epsilon();
    DenseMatrix<T> w_col;

    const int m = W.Height();
    const int k = W.Width();
    assert(Norms.Height() == k);

    for (int c=0; c<k; ++c)
    {
        // create a view of the c'th column
        View(w_col, W, 0, c, m, 1);
        
        // compute the 2-norm of this col vector
        T two_norm = Nrm2(w_col);
        if (std::abs(two_norm) < eps)
            throw std::runtime_error("Normalize: column norm < machine epsilon");

        // scale this column by the reciprocal of the 2-norm
        T inv_two_norm = T(1.0) / two_norm;
        Scal(inv_two_norm, w_col);

        Norms.Set(c, 0, two_norm);
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ScaleRowsRank2(DenseMatrix<T>& H, DenseMatrix<T>& scale_factors)
{
    // Specialization of ScaleRows for rank2.

    const int k = H.Height();
    const int n = H.Width();
    assert(2 == k);
    assert(scale_factors.Height() >= k);

    const T f0 = scale_factors.Get(0, 0);
    const T f1 = scale_factors.Get(1, 0);

    T* buf_h = H.Buffer();
    int ldim_h = H.LDim();

    OPENMP_PRAGMA(omp parallel for)
    for (int c=0; c<n; ++c)
    {
        unsigned int col_offset = c*ldim_h;

        T h0 = buf_h[col_offset + 0];
        T h1 = buf_h[col_offset + 1];
        
        h0 *= f0;
        h1 *= f1;
        
        buf_h[col_offset + 0] = h0;
        buf_h[col_offset + 1] = h1;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ScaleRows(DenseMatrix<T>& H, DenseMatrix<T>& scale_factors)
{
    // Multiply each row of H by the corresponding scale factor.  This could
    // be slow if H.Height() << H.Width(), as it usually is for NMF, since
    // Elemental stores matrices in column-major order.

    const int k = H.Height();
    const int n = H.Width();
    assert(scale_factors.Height() >= k);

    if (2 == k)
        ScaleRowsRank2(H, scale_factors);
    else
    {
        DenseMatrix<T> h_row;    
        for (int r=0; r<k; ++r)
        {
            T f = scale_factors.Get(r, 0);
            
            // create a view of the rth row of H and scale
            View(h_row, H, r, 0, 1, n);
            Scal(f, h_row);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void NormalizeAndScale(DenseMatrix<T>& W, 
                       DenseMatrix<T>& H)
 {
    // Compute a diagonal matrix D that simultaneously normalizes the columns
    // of W and scales the rows of H so that the product W*H is unchanged.
    //
    // In other words, compute matrix D so that:
    //
    //    W*H == W * D * D' * H
    //
    //    L2 norm of each col of W is 1
    //
    // The matrix D does not need to be explicitly formed.

    if (W.Width() != H.Height())
        throw std::logic_error("NormalizeAndScale: non-conformant W and H");

    const int k = W.Width();
    DenseMatrix<T> col_norms(k, 1);

    NormalizeColumns(W, col_norms);
    ScaleRows(H, col_norms);
}

//-----------------------------------------------------------------------------
template <typename T>
void NormalizeAndScale(DenseMatrix<T>& W,
                       DenseMatrix<T>& H,
                       DenseMatrix<T>& col_norms) // [out] col norms of W
{
    // Same as NormalizeAndScale, but returns the column norms of W.

    if (W.Width() != H.Height())
        throw std::logic_error("NormalizeAndScale: non-conformant W and H");

    const int k = W.Width();

    if (col_norms.Height() < k)
        throw std::logic_error("NormalizeAndScale: norm array too small");

    NormalizeColumns(W, col_norms);
    ScaleRows(H, col_norms);

}
