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

#include <vector>
#include "random.hpp"

//-----------------------------------------------------------------------------
template <typename T>
bool RandomMatrix(std::vector<T>& buf,
                  const unsigned int height, // number of data rows, <= ldim
                  const unsigned int width,  // number of data cols
                  Random& rng,
                  const T rng_center = T(0.5), 
                  const T rng_radius = T(0.5))
{
    // Generate a height*width matrix of random numbers.  Store data in
    // buf in column-major order.

    unsigned int required_size = height*width;
    if (buf.size() < required_size)
        buf.resize(required_size);
    
    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int col_offset = c*height;
        for (unsigned int r=0; r<height; ++r)
        {
            buf[r + col_offset] = 
                static_cast<T>(rng.RandomDouble(rng_center, rng_radius));
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
bool RandomMatrix(T* buf,
                  const unsigned int ldim,   // leading dimension of buffer
                  const unsigned int height, // number of data rows, <= ldim
                  const unsigned int width,  // number of data cols
                  Random& rng,
                  const T rng_center = T(0.5), 
                  const T rng_radius = T(0.5))
{
    // Generate a height*width matrix of random numbers.  Store data in
    // buf in column-major order.  Assumes the buffer is large enough.

    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int col_offset = c*ldim;
        for (unsigned int r=0; r<height; ++r)
        {
            buf[r + col_offset] = 
                static_cast<T>(rng.RandomDouble(rng_center, rng_radius));
        }
    }

    return true;
}

// //-----------------------------------------------------------------------------
// template <typename T>
// bool RandomMatrix(elem::Matrix<T>& M, 
//                   Random& rng,
//                   const T rng_center = T(0.5), 
//                   const T rng_radius = T(0.5))
// {
//     // Fill Elemental matrix M with random entries.

//     int height    = M.Height();
//     int width     = M.Width();
    
//     for (int c=0; c<width; ++c)
//     {
//         for (int r=0; r<height; ++r)
//         {
//             M.Set(r, c, static_cast<T>(rng.RandomDouble(rng_center, rng_radius)));
//         }
//     }

//     return true;
// }

// //-----------------------------------------------------------------------------
// template <typename T>
// void RandomSparseMatrix(SparseMatrix<T>& A, 
//                         const unsigned int nonzeros_per_column,
//                         Random& rng,
//                         const unsigned int height,
//                         const unsigned int width)
// {
//     // Generate a random sparse matrix with 'nonzeros_per_column' nonzero
//     // entries in each column.  The entries in each column are in random
//     // rows.

//     std::vector<unsigned int> rows(height, 0);
//     std::vector<unsigned int> cols(width, 0);
//     std::vector<unsigned int> nz(width, 0);

//     unsigned int nz_per_col = nonzeros_per_column;
//     if (0 == nz_per_col)
//         nz_per_col = 1;
//     if (nz_per_col > height)
//         nz_per_col = height;

//     // aux vector of row indices
//     for (unsigned int k=0; k != height; ++k)
//         rows[k] = k;

//     unsigned int nnz = nz_per_col * width;

//     A.Clear();
//     A.Reserve(height, width, nnz);

//     // load num random entries of A
//     A.BeginLoad();
    
//     for (unsigned int c=0; c != width; ++c)
//     {
//         // generate row indices in random order
//         std::random_shuffle(rows.begin(), rows.begin() + height);
        
//         for (unsigned int k=0; k != nz_per_col; ++k)
//             A.Load(rows[k], c, rng.RandomDouble(0.5, 0.5));
//     }
    
//     A.EndLoad();
// }

