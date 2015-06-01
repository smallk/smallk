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

#include <cstring>
#include <cassert>
#include <stdexcept>
#include <vector>
#include "enums.hpp"

#if ELEM_VER >= 85
#include "El.hpp"
namespace EL = El;
#else
#include "elemental.hpp"
namespace EL = elem;
#endif

#include "openmp_pragma.hpp"

//-----------------------------------------------------------------------------
template <typename T>
class DenseMatrix
{
    // This class is a thin wrapper for a subset of the methods of
    // elem::Matrix<T, Int>.  A few simple methods have also been added.
    // The purpose of this class is to simplify generic NMF solvers.  The 
    // 'DenseMatrix' wrapper hides the Elemental 'Int' template param, which 
    // we don't need.  It also provides a mechanism to hide any changes in 
    // the Elemental API.

public:

    DenseMatrix() 
        : M() {}
    DenseMatrix(int height, int width) 
        : M(height, width) {}
    DenseMatrix(int height, int width, T* buffer, int ldim) 
        : M(height, width, buffer, ldim) {}

    int LDim()   const {return M.LDim();}
    int Width()  const {return M.Width();}
    int Height() const {return M.Height();}
    
    // given zero-based row and column indices, return the zero-based linear index
    int Sub2Ind(const int row_index, const int col_index) const;

    // given a zero-based linear index, return the zero-based row and column indices
    void Ind2Sub(const int linear_index, int& row_index, int& col_index) const;

    T* Buffer()                     {return M.Buffer();}
    const T* LockedBuffer() const   {return M.LockedBuffer();}

    bool Locked() const             {return M.Locked();}

    T Get(int i, int j) const       {return M.Get(i, j);}
    void Set(int i, int j, T alpha) {M.Set(i, j, alpha);}

    void Attach(int height, int width, T* buffer, int ldim);

    DenseMatrix<T>& operator=(const T val);
    const DenseMatrix<T>& operator=(const DenseMatrix<T>& rhs);

    // extract entire columns
    void SubmatrixFromCols(DenseMatrix<T>& result, 
                           const std::vector<unsigned int>& col_indices) const;

    // extract entire rows
    void SubmatrixFromRows(DenseMatrix<T>& result, 
                           const std::vector<unsigned int>& row_indices) const;

    // assign val to the elements at the given linear indices
    void LinearIndexedAssign(const std::vector<int>& linear_indices, const T val);

    void Empty() {M.Empty();}
    void Resize(int height, int width);

    void SubMatrixColsCompact(DenseMatrix<T>& result,
                              const std::vector<unsigned int>& col_indices,
                              std::vector<unsigned int>& old_to_new_rows,
                              std::vector<unsigned int>& new_to_old_rows) const;

    // extract submatrix at the intersection of the row and column indices
    void Submatrix(DenseMatrix<T>& result,
                   const std::vector<unsigned int>& row_indices,
                   const std::vector<unsigned int>& col_indices,
                   const unsigned int num_rows,
                   const unsigned int num_cols) const;

    EL::Matrix<T> M;
};
