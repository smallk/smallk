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
#include <iostream>
#include <stdexcept>
#include "dense_matrix.hpp"
#include "dense_matrix_ops.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void CholeskyIndexed(const std::vector<unsigned int>& ri,
                     const unsigned int num_rows,
                     const DenseMatrix<T>& WtW,
                     DenseMatrix<T>& U)
{
    for (unsigned int r=0; r<num_rows; ++r)
    {
        unsigned int row_index = ri[r];
        
        // compute diagonal element
        T Urr = WtW.Get(row_index, row_index);
        for (unsigned int q=0; q<r; ++q)
        {
            T Uqr = U.Get(q, r);
            Urr -= Uqr*Uqr;
        }
        
        if (Urr <= T(0))
            throw std::logic_error("CholeskyIndexed: Matrix is not SPD.");

        Urr = std::sqrt(Urr);
        U.Set(r, r, Urr);

        // solve for remaining elements in row r
        T diag = T(1)/Urr;
        for (unsigned int c=r+1; c<num_rows; ++c)
        {
            T Urc = WtW.Get(row_index, ri[c]);
            for (unsigned int q=0; q<r; ++q)
                Urc -= U.Get(q, r)*U.Get(q, c);

            U.Set(r, c, Urc * diag);
        }
    }    
}

//-----------------------------------------------------------------------------
template <typename T>
void CholeskySolveIndexed(const std::vector<unsigned int>& ri, // indices
                          const unsigned int num_rows, // use this many indices
                          const DenseMatrix<T>& WtW, // cross-product matrix
                          const DenseMatrix<T>& R,   // subset of the RHS vector
                          DenseMatrix<T>& X)         // subset of the solution vector
{
    // Solve WtW(ri, ri) * X = R.  Keep data in place and use
    // indices to find elements on which to operate.

    // Step 1: Cholesk factorization of WtW(ri, ri) = U'U
    DenseMatrix<T> U(num_rows, num_rows);
    CholeskyIndexed(ri, num_rows, WtW, U);
    
    // Step 2: Solve U'UX(ri) = R(ri).  Let UX(ri) = Y, solve U'Y = R(ri) 
    // for Y.  This is a lower-triangular system, so solve by forward 
    // substitution.

    std::vector<T> Y(num_rows);
    for (unsigned int r=0; r<num_rows; ++r)
    {
        T Urr = U.Get(r, r);

        T sum = T(0);
        for (unsigned int c=0; c<r; ++c)
            sum += U.Get(c, r)*Y[c];
        
        Y[r] = (R.Get(r, 0) - sum) / Urr;
    }

    // Step3: Solve UX(ri) = Y by backsubstitution.

    for (unsigned int r=num_rows-1; r >= 0; --r)
    {
        T Urr = U.Get(r, r);

        T sum = T(0);
        for (unsigned int q=r+1; q<num_rows; ++q)
            sum += U.Get(r, q)*X.Get(q, 0);

        X.Set(r, 0, (Y[r] - sum) / Urr);

        // prevent wrap-around at 0
        if (0 == r)
            break;
    }
}
                          
