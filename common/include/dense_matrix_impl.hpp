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

#include <algorithm>

//-----------------------------------------------------------------------------
template <typename T>
void 
DenseMatrix<T>::Resize(int height, int width)
{
    // Name changed from 'ResizeTo' to just 'Resize' at Elemental release 0.83.

#if ELEM_VER >= 83
    M.Resize(height, width);
#else
    M.ResizeTo(height, width);
#endif
}

//-----------------------------------------------------------------------------
template <typename T>
int 
DenseMatrix<T>::Sub2Ind(const int row_index, const int col_index) const
{
    return col_index*M.LDim() + row_index;
}

//-----------------------------------------------------------------------------
template <typename T>
void 
DenseMatrix<T>::Ind2Sub(const int linear_index,  int& row_index, int& col_index) const
{
    col_index = linear_index / M.LDim();
    row_index = linear_index - col_index * M.LDim();
}

//-----------------------------------------------------------------------------
template <typename T>
void 
DenseMatrix<T>::Attach(int height, int width, T* buffer, int ldim) 
{
    M.Attach(height, width, buffer, ldim);
}

//-----------------------------------------------------------------------------
template <typename T>
DenseMatrix<T>& 
DenseMatrix<T>::operator=(const T val)
{
    T* buf     = M.Buffer();
    int ldim   = M.LDim();
    int height = M.Height();
    int width  = M.Width();

    for (int c=0; c<width; ++c)
    {
        int offset = c*ldim;
        for (int r=0; r<height; ++r)
            buf[offset + r] = val;
    }

    return *this;
}

//-----------------------------------------------------------------------------
template <typename T>
const DenseMatrix<T>& 
DenseMatrix<T>::operator=(const DenseMatrix<T>& rhs) 
{
    if (this != &rhs)
        M = rhs.M;
    
    return *this;
}

//-----------------------------------------------------------------------------
template <typename T>
void 
DenseMatrix<T>::LinearIndexedAssign(const std::vector<int>& linear_indices, 
                                    const T val)
{
    T* buf = M.Buffer();
    const int max_index = M.LDim() * (M.Width()-1) + (M.Height()-1);

    for (unsigned int i=0; i != linear_indices.size(); ++i)
    {
        int index = linear_indices[i];
        assert(index >= 0);
        assert(index < max_index);

        buf[index] = val;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void
DenseMatrix<T>::SubmatrixFromCols(DenseMatrix<T>& result,
                                  const std::vector<unsigned int>& col_indices) const
{
    const unsigned int new_width = col_indices.size();
    if (0 == new_width)
        throw std::logic_error("DenseMatrix::operator(): empty index set");

    const T* buf_m = M.LockedBuffer();
    const unsigned int ldim_m = M.LDim();
    const unsigned int height = M.Height();
    const unsigned int width  = M.Width();

    result.Resize(height, new_width);
    T* buf_r = result.Buffer();
    const unsigned int ldim_r = result.LDim();

    OPENMP_PRAGMA(omp parallel for)
    for (unsigned int c=0; c<new_width; ++c)
    {
        unsigned int source_col = col_indices[c];
        assert(source_col < width);
        
        unsigned int dest_col_offset   = c*ldim_r;
        unsigned int source_col_offset = source_col*ldim_m;

        memcpy(&buf_r[dest_col_offset], 
               &buf_m[source_col_offset], 
               height*sizeof(T));
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void
DenseMatrix<T>::SubmatrixFromRows(DenseMatrix<T>& result,
                                  const std::vector<unsigned int>& row_indices) const
{
    const unsigned int new_height = row_indices.size();
    if (0 == new_height)
        throw std::logic_error("DenseMatrix::operator(): empty row index set");

    const unsigned int ldim   = M.LDim();
    const unsigned int height = M.Height();
    const unsigned int width  = M.Width();
    const T* buf_m = M.LockedBuffer();

    result.Resize(new_height, width);
    T* buf_r = result.Buffer();
    const unsigned int ldim_r = result.LDim();

    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int source_col_offset = c*ldim; 
        unsigned int dest_col_offset   = c*ldim_r;
        for (unsigned int r=0; r<new_height; ++r)
        {
            unsigned int row_index = row_indices[r];
            assert(row_index < height);
            buf_r[dest_col_offset + r] =
                buf_m[source_col_offset + row_index];
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void 
DenseMatrix<T>::Submatrix(DenseMatrix<T>& result,
                          const std::vector<unsigned int>& row_indices,
                          const std::vector<unsigned int>& col_indices,
                          const unsigned int num_rows,
                          const unsigned int num_cols) const
{
    const unsigned int new_width = num_cols;
    const unsigned int new_height = num_rows;

    if (0 == new_width)
        throw std::logic_error("DenseMatrix::Submatrix: empty column index set");
    if (0 == new_height)
        throw std::logic_error("DenseMatrix::Submatrix: empty row index set");
    if (static_cast<unsigned int>(M.Height()) < new_height)
        throw std::logic_error("DenseMatrix::Submatrix: submatrix height too large");
    if (static_cast<unsigned int>(M.Width()) < new_width)
        throw std::logic_error("DenseMatrix::Submatrix: submatrix width too large");

    const unsigned int ldim   = M.LDim();
    const unsigned int height = M.Height();
    const unsigned int width  = M.Width();
    const T* buf_m = M.LockedBuffer();

    result.Resize(new_height, new_width);
    T* buf_r = result.Buffer();
    const unsigned int ldim_r = result.LDim();

    for (unsigned int c=0; c<new_width; ++c)
    {
        unsigned int source_col = col_indices[c];
        assert(source_col < width);
        
        unsigned int dest_col_offset   = c*ldim_r;
        unsigned int source_col_offset = source_col*ldim; 
        for (unsigned int r=0; r<new_height; ++r)
        {
            unsigned int row_index = row_indices[r];
            assert(row_index < height);
            buf_r[dest_col_offset + r] = buf_m[source_col_offset + row_index];
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void 
DenseMatrix<T>::SubMatrixColsCompact(DenseMatrix<T>& result,
                                     const std::vector<unsigned int>& col_indices,
                                     std::vector<unsigned int>& old_to_new_rows,
                                     std::vector<unsigned int>& new_to_old_rows) const
{
    int height = this->M.Height();
    int width = this->M.Width();

    unsigned int uh = static_cast<unsigned int>(height);
    unsigned int uw = static_cast<unsigned int>(width);

    // Extract entire columns from the source matrix to form the dest matrix.
    // Unlike the sparse case, complete columns are extracted, so the height
    // of the submatrix is always the same as the source matrix, and the row
    // mappings are one-to-one.
    
    unsigned int new_width = col_indices.size();
    if (0u == new_width)
        throw std::logic_error("DenseMatrix::SubMatrixColsCompact: empty column set");
    
    // check the column indices for validty
    for (auto it=col_indices.begin(); it != col_indices.end(); ++it)
    {
        // index of next source column
        unsigned int c = *it;
        if (c >= uw)
            throw std::logic_error("DenseMatrix::SubMatrixColsCompact: column index out of range");
    }

    // allocate memory in the result
    result.Resize(height, new_width);

    // construct the one-to-one row mappings
    if (old_to_new_rows.size() < uh)
        old_to_new_rows.resize(height);
    if (new_to_old_rows.size() < uh)
        new_to_old_rows.resize(height);

    for (unsigned int r=0; r<uh; ++r)
        old_to_new_rows[r] = r;
    for (unsigned int r=0; r<uh; ++r)
        new_to_old_rows[r] = r;

    EL::Matrix<T> source_col, dest_col;

    // copy the columns from the source matrix to the dest matrix
    int c_dest = 0;
    for (auto it=col_indices.begin(); it != col_indices.end(); ++it, ++c_dest)
    {
        // index of the next source column
        int c = *it;

        EL::LockedView(source_col, this->M,   0, c,      height, 1);
        EL::View(dest_col,   result.M,  0, c_dest, height, 1);
        Copy(source_col, dest_col);
    }
}
