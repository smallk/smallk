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

//-----------------------------------------------------------------------------
template <typename T>
DenseMatrix<T> MakeDense(const SparseMatrix<T>& A)
{
    // Return a fully dense matrix equivalent to SparseMatrix A.

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();
    const unsigned int height = A.Height();
    const unsigned int width  = A.Width();

    DenseMatrix<T> D(height, width);
    MakeZeros(D);

    for (unsigned int c=0; c != width; ++c)
    {
        unsigned int start = cols_a[c];
        unsigned int end   = cols_a[c+1];
        for (unsigned int offset=start; offset != end; ++offset)
        {
            unsigned int row = rows_a[offset];
            D.Set(row, c, data_a[offset]);
        }
    }

    return D;
}

//-----------------------------------------------------------------------------
template <typename T>
void Add(DenseMatrix<T>& D, const T& alpha, const SparseMatrix<T>& A)
{
    // matrix D is dense, matrix A is sparse
    // computes D = D + alpha*A

    if ( (unsigned int)D.Width() != A.Width())
        throw std::logic_error("Add: widths of dense and sparse matrix do not match");
    if ( (unsigned int)D.Height() != A.Height())
        throw std::logic_error("Add: heights of dense and sparse matrix do not match");

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();
    const unsigned int width_a = A.Width();

    for (unsigned int c=0; c != width_a; ++c)
    {
        unsigned int start = cols_a[c];
        unsigned int   end = cols_a[c+1];
        for (unsigned int offset=start; offset != end; ++offset)
        {
            unsigned int row = rows_a[offset];
            T A_rc = data_a[offset];

            T D_rc = D.Get(row, c);
            D_rc += alpha*A_rc;
            D.Set(row, c, D_rc);
        }
    }
}
