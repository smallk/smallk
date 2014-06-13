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

#include <iostream>
#include "dense_matrix.hpp"
#include "bit_matrix_ops.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
bool IsEmpty(const BitMatrix& A)
{
    return (0 == A.SizeInBits());
}

//-----------------------------------------------------------------------------
bool All(const BitMatrix& A)
{
    // return true if all bits are set, false otherwise

    const unsigned int height = A.Height();
    const unsigned int width  = A.Width();
    std::vector<unsigned int> col_sums(width);

    // Sum the columns of A; if any sum is not equal to the height
    // that column must contain at least one zero bit.

    A.SumColumns(col_sums);

    bool ok = true;
    for (unsigned int c=0; c<width; ++c)
    {
        if (height != col_sums[c])
        {
            ok = false;
            break;
        }
    }

    return ok;
}

//-----------------------------------------------------------------------------
bool AllCols(const BitMatrix& A, const std::vector<unsigned int>& col_indices)
{
    // return true if all bits in the specified cols are set

    const unsigned int height = A.Height();
    const unsigned int num_cols = col_indices.size();

    for (unsigned int c=0; c<num_cols; ++c)
    {
        unsigned int col_index = col_indices[c];
        if (height != A.ColumnSum(col_index))
            return false;
    }

    return true;
}

//-----------------------------------------------------------------------------
bool AllRows(const BitMatrix& A, const std::vector<unsigned int>& row_indices)
{
    // return true if all bits in the specified rows are set

    const unsigned int width = A.Height();
    const unsigned int num_rows = row_indices.size();

    for (unsigned int r=0; r<num_rows; ++r)
    {
        unsigned int row_index = row_indices[r];
        if (width != A.RowSum(row_index))
            return false;
    }

    return true;
}

//-----------------------------------------------------------------------------
BitMatrix ColumnwiseAND(const BitMatrix& B, const BitMatrix& mask)
{
    const unsigned int height = B.Height();
    const unsigned int width  = B.Width();

    if (mask.Height() != width)
        throw std::logic_error("ColumnwiseAND: mask height must match matrix width");
    if (mask.Width() != 1)
        throw std::logic_error("ColumnwiseAND: mask must be a column vector");
    
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    // For each bit in mask ('width' of them), compute the bitwise AND with
    // an entire column of B.  Return the resulting BitMatrix.

    const unsigned int* buf_m = mask.LockedBuffer();
    const unsigned int* buf_b = B.LockedBuffer();
    const unsigned int ldim_b = B.LDim();

    BitMatrix result(height, width);
    unsigned int* buf_r = result.Buffer();
    unsigned int ldim_r = result.LDim();

    const unsigned int full_wds = height / BITS;
    const unsigned int extra    = height - full_wds*BITS;

    for (unsigned int c=0; c<width; ++c)
    {
        // get bit c of mask and set the mask value
        unsigned int full_mask_wds   = c / BITS;
        unsigned int extra_mask_bits = c - full_mask_wds*BITS;
        unsigned int mask_val        = buf_m[full_mask_wds] & (1 << extra_mask_bits);
        mask_val = (mask_val > 0 ? 0xFFFFFFFF : 0x00);

        unsigned int r_wd = 0u;
        unsigned int col_offset_r = c*ldim_r;
        unsigned int col_offset_b = c*ldim_b;

        for (; r_wd < full_wds; ++r_wd)
            buf_r[col_offset_r + r_wd] = buf_b[col_offset_b + r_wd] & mask_val;

        if (extra > 0)
            buf_r[col_offset_r + r_wd] = buf_b[col_offset_b + r_wd] & mask_val;
    }

    return result;
}

//-----------------------------------------------------------------------------
BitMatrix MatrixFromColumnMask(const BitMatrix& mask, const unsigned int height)
{
    // Perform the equivalent of repmat(mask, height, 1); mask.Height() is the
    // WIDTH of the new matrix (one bit per column).  Each bit in the mask is 
    // the value to replicate down each column.

    const unsigned int width = mask.Height();

    if (mask.Width() != 1)
        throw std::logic_error("MatrixFromColumnMask: mask must be a column vector");
    
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    const unsigned int* buf_m = mask.LockedBuffer();

    BitMatrix result(height, width);
    unsigned int* buf_r = result.Buffer();
    unsigned int ldim_r = result.LDim();

    const unsigned int full_wds = height / BITS;
    const unsigned int extra    = height - full_wds*BITS;

    for (unsigned int c=0; c<width; ++c)
    {
        // get bit c of mask and set the mask value
        unsigned int full_mask_wds   = c / BITS;
        unsigned int extra_mask_bits = c - full_mask_wds*BITS;
        unsigned int mask_val        = buf_m[full_mask_wds] & (1 << extra_mask_bits);
        mask_val = (mask_val > 0 ? 0xFFFFFFFF : 0x00);

        unsigned int r_wd = 0u;
        unsigned int col_offset_r = c*ldim_r;

        for (; r_wd < full_wds; ++r_wd)
            buf_r[col_offset_r + r_wd] = mask_val;

        if (extra > 0)
            buf_r[col_offset_r + r_wd] = mask_val;
    }

    return result;
}
