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
#include "population_count.hpp"

//-----------------------------------------------------------------------------
class BitMatrix
{
    // This is a Boolean matrix in which each element is a single bit.  
    // The internal buffer is an array of unsigned ints, packed columnwise.
    // There are 'width' columns of unsigned ints in the underlying buffer.  
    // Unused bits in the final word in each col are set to zero.  This 
    // condition is enforced in all operations.

public:

    static constexpr unsigned int BITS_PER_WORD = sizeof(unsigned int) * 8;

    BitMatrix() {}
    BitMatrix(const BitMatrix& rhs);

    // height and width specify the number of bits in each dimension
    BitMatrix(const unsigned int height, const unsigned int width);

    // creates a single-column bit vector
    explicit BitMatrix(const unsigned int size) : BitMatrix(size, 1) {}

    // sets all words equal to val
    BitMatrix(const unsigned int height, 
              const unsigned int width, 
              const unsigned int val);

    // Copy 'count' words from 'buf' into the internal buffer.  The final 
    // bits in each column will be zeroed according to the internal MASK_.  
    // The value of 'count' must be less than or equal to SizeInWords().
    void Load(const unsigned int* buf, const unsigned int count);

    // need move c'tor and move assignment op - TBD

    unsigned int Mask() const            {return MASK_;}
    unsigned int Width() const           {return width_;}
    unsigned int Height() const          {return height_;}
    unsigned int LDim() const            {return ldim_wds_;}
    unsigned int SizeInBits() const      {return height_ * width_;}
    unsigned int SizeInWords() const     {return words_.size();}

          unsigned int* Buffer()             {return &words_[0];}
    const unsigned int* LockedBuffer() const {return &words_[0];}

    // set all words equal to 'val'; final word in each column is masked
    void Fill(const unsigned int val);

    // resize the BitMatrix; reallocates only if existing buffer too small
    void Resize(const unsigned int new_height, const unsigned int new_width);

    static unsigned int HeightInWords(const unsigned int height)
    {return (height + BITS_PER_WORD-1)/BITS_PER_WORD;} 

    // return a new BitMatrix with all bits flipped
    BitMatrix  operator~() const;

    BitMatrix& operator&=(const BitMatrix& rhs);
    BitMatrix& operator|=(const BitMatrix& rhs);
    BitMatrix& operator^=(const BitMatrix& rhs);
    BitMatrix& operator=(const BitMatrix& rhs);

    // set all bits to either 1 or 0
    BitMatrix& operator=(const unsigned int val);

    // sum the bits in all rows; result will be resized if necessary
    template <typename IntType>
    void SumRows(std::vector<IntType>& row_sums) const;

    // sum the bits in all columns; result will be resized if necessary
    template <typename IntType>
    void SumColumns(std::vector<IntType>& col_sums) const;

    // hash the columns of the bit matrix
    void GroupIdenticalColumns(std::vector<unsigned int>& indices,
                               std::vector<unsigned int>& groups) const;

    void GroupIdenticalColumns(std::vector<unsigned int>& indices,
                               std::vector<unsigned int>& groups,
                               const std::vector<unsigned int>& col_indices) const;

    // sum the bits in the specified column and return the result
    unsigned int ColumnSum(const unsigned int col_index) const;

    // sum the bits in the specified row and return the result
    unsigned int RowSum(const unsigned int row_index) const;

    // return the linear indices of all nonzero bits
    template <typename IntType>
    void Find(std::vector<IntType>& linear_indices) const;

    // assign 0 or 1 to the bits at the specified linear indices
    void LinearIndexedAssign(const std::vector<unsigned int>& linear_indices,
                             const unsigned int val);

    // return a vector of row indices for all nonzero bits in column c
    void RowIndices(const unsigned int c,
                    std::vector<unsigned int>& row_indices,
                    unsigned int& num_rows) const;

    // return a vector of col indices for all nonzero bits in row r
    void ColIndices(const unsigned int r,
                    std::vector<unsigned int>& col_indices,
                    unsigned int& num_cols) const;

    // return the maximum row index of all the nonzero bits in column c
    unsigned int MaxRowIndex(const unsigned int c) const;

    // extract a submatrix from the given column indices
    void SubmatrixFromCols(BitMatrix& result, const BitMatrix& mask) const;
    void SubmatrixFromCols(BitMatrix& result, const unsigned int col_index) const;
    void SubmatrixFromCols(BitMatrix& result, const std::vector<unsigned int>& col_indices) const;

    // set or clear bits wherever mask == 1
    void SetBits(const BitMatrix& mask);
    void ClearBits(const BitMatrix& mask);
    void ToggleBit(const unsigned int row, const unsigned int col);

    void Print() const;

private:
    
    // dimensions 
    unsigned int height_;    // height of the buffer in bits
    unsigned int width_;     // width in bits
    unsigned int ldim_wds_;  // height of the buffer in words

    unsigned int full_wds_;  // number of complete words in each col
    unsigned int MASK_;      // mask for final word in each col

    std::vector<unsigned int> words_;

    void SetMask();
    void Assign(const BitMatrix& rhs);
};

//-----------------------------------------------------------------------------
template <typename IntType>
void 
BitMatrix::SumRows(std::vector<IntType>& row_sums) const
{
    if (row_sums.size() < height_)
        row_sums.resize(height_);

    for (unsigned int r=0; r<height_; ++r)
        row_sums[r] = (IntType)0;

    // for each word in each column, scan the bits of this word and
    // accumulate in the appropriate result element

    unsigned int r_wd=0, offset;
    for (; r_wd < full_wds_; ++r_wd)
    {
        // offset in the row_sums array
        offset = r_wd * BITS_PER_WORD;
        for (unsigned int c=0; c<width_; ++c)
        {
            unsigned int wd = words_[c*ldim_wds_ + r_wd];
            for (unsigned int q=0; q<BITS_PER_WORD; ++q)
            {
                if (wd & (1 << q))
                    row_sums[offset + q] += 1;
            }
        }
    }
    
    offset = r_wd*BITS_PER_WORD;
    unsigned int extra = height_ - full_wds_*BITS_PER_WORD;

    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int wd = words_[c*ldim_wds_ + r_wd];
        for (unsigned int q=0; q<extra; ++q)
        {
            if (wd & (1 << q))
                row_sums[offset + q] += 1;
        }
    }
}

//-----------------------------------------------------------------------------
template <typename IntType>
void 
BitMatrix::SumColumns(std::vector<IntType>& col_sums) const
{
    if (col_sums.size() < width_)
        col_sums.resize(width_);

    for (unsigned int c=0; c<width_; ++c)
    {
        // offset to column c
        unsigned int col_offset = c*ldim_wds_, r_wd = 0u;

        IntType sum=0;
        for (; r_wd<full_wds_; ++r_wd)
            sum += (IntType) PopulationCount(words_[col_offset + r_wd]);

        if (MASK_ > 0)
            sum += (IntType) PopulationCount(MASK_ & words_[col_offset + r_wd]);

        col_sums[c] = sum;
    }
}

//-----------------------------------------------------------------------------
template <typename IntType>
void
BitMatrix::Find(std::vector<IntType>& linear_indices) const
{
    // mimics Matlab 'find': returns LINEAR indices for nonzero elements
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    linear_indices.clear();
    unsigned int extra = height_ - BITS*full_wds_;

    for (unsigned int c=0; c != width_; ++c)
    {
        unsigned int wd, r_wd = 0u, r=0u, col_offset = c*ldim_wds_;

        // scan all bits in each full word in column c
        for (; r_wd < full_wds_; ++r_wd)
        {
            wd = words_[col_offset + r_wd];
            for (unsigned int q=0u; q<BITS; ++q, ++r)
            {
                if (wd & (1 << q))
                    linear_indices.push_back( IntType(c*height_ + r) );
            }
        }

        // scan extra bits, if any, in the final word of column c
        if (extra > 0)
        {
            wd = words_[col_offset + r_wd];
            for (unsigned int q=0u; q<extra; ++q, ++r)
            {
                if (wd & (1 << q))
                    linear_indices.push_back( IntType(c*height_ + r) );
            }
        }
    }
}
