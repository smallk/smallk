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

#include <algorithm>
#include <iostream>
#include <cassert>
#include <limits>
#include "utils.hpp"
#include "spooky_v2.hpp"
#include "bit_matrix.hpp"
#include "bit_matrix_ops.hpp"

using std::cout;
using std::endl;

//-----------------------------------------------------------------------------
void BitMatrix::SetMask()
{
    // The mask contains 1 bits for all valid rows in the final word.

    assert(height_ > 0);
    full_wds_ = height_ / BITS_PER_WORD;
    unsigned int extra = height_ - (full_wds_ * BITS_PER_WORD);
    assert(extra >= 0);
    assert(extra < BITS_PER_WORD);

    // bit mask for the final word in each column; set as many bits as
    // the value of extra
    //MASK_ = std::numeric_limits<unsigned int>::max() >> (BITS_PER_WORD - extra);
    MASK_ = 0;
    for (unsigned int q=0; q<extra; ++q)
        MASK_ |= (1 << q);
}

//-----------------------------------------------------------------------------
BitMatrix::BitMatrix(const unsigned int height, const unsigned int width)
{
    // dimensions in bits
    height_ = height;
    width_ = width;
    ldim_wds_ = HeightInWords(height_);
    assert(ldim_wds_ >= 1);
    SetMask();
    words_.resize(ldim_wds_ * width_);
}

//-----------------------------------------------------------------------------
BitMatrix::BitMatrix(const unsigned int height, 
                     const unsigned int width, 
                     const unsigned int val)
    : BitMatrix(height, width)
{
    std::fill(words_.begin(), words_.end(), val);
}

//-----------------------------------------------------------------------------
void BitMatrix::Load(const unsigned int* buf, const unsigned int count)
{
    if (count > words_.size())
        throw std::logic_error("BitMatrix::Load: count too large");

    unsigned int buf_index = 0;
    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int r_wd = 0;
        unsigned int col_offset = c*ldim_wds_;

        // copy complete words
        for (; r_wd < full_wds_; ++r_wd)
        {
            words_[col_offset + r_wd] = buf[buf_index++];
            if (count == buf_index)
                return;
        }

        // copy final word and mask off unused bits
        if (MASK_ > 0)
            words_[col_offset + r_wd] = MASK_ & buf[buf_index++];
    }
}

//-----------------------------------------------------------------------------
void BitMatrix::Resize(const unsigned int new_height, 
                       const unsigned int new_width)
{
    height_ = new_height;
    width_  = new_width;
    ldim_wds_ = HeightInWords(height_);

    unsigned int required_size_wds = ldim_wds_ * new_width;
    if (words_.size() < required_size_wds)
        words_.resize(required_size_wds);

    SetMask();
}

//-----------------------------------------------------------------------------
void BitMatrix::Assign(const BitMatrix& rhs)
{
    height_   = rhs.height_;
    width_    = rhs.width_;
    ldim_wds_ = rhs.ldim_wds_;
    full_wds_ = rhs.full_wds_; 
    MASK_     = rhs.MASK_;
    words_    = rhs.words_;
}

//-----------------------------------------------------------------------------
BitMatrix BitMatrix::operator~() const
{
    BitMatrix result(height_, width_);
    unsigned int* buf = result.Buffer();

    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int wd_offset = c*ldim_wds_;

        unsigned int r_wd=0;
        for (; r_wd<full_wds_; ++r_wd)
            buf[wd_offset + r_wd] = ~words_[wd_offset + r_wd];

        // final word
        if (MASK_ > 0)
            buf[wd_offset + r_wd] = MASK_ & ~words_[wd_offset + r_wd];
    }

    return result;
}

//-----------------------------------------------------------------------------
BitMatrix& BitMatrix::operator&=(const BitMatrix& rhs)
{
    if ( (rhs.Height() != height_) && (rhs.Width() != width_))
        throw std::logic_error("BitMatrix::operator&=: non-conformant matrices");

    const unsigned int ldim_rhs = rhs.LDim();
    const unsigned int* buf_rhs = rhs.LockedBuffer();

    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int r_wd = 0;
        unsigned int col_offset = c*ldim_wds_;
        unsigned int col_offset_rhs = c*ldim_rhs;

        for (; r_wd < full_wds_; ++r_wd)
            words_[col_offset + r_wd] &= buf_rhs[col_offset_rhs + r_wd];

        if (MASK_ > 0)
            words_[col_offset + r_wd] &= (MASK_ & buf_rhs[col_offset_rhs + r_wd]);
    }

    return *this;
} 

//-----------------------------------------------------------------------------
BitMatrix& BitMatrix::operator|=(const BitMatrix& rhs)
{
    if ( (rhs.Height() != height_) && (rhs.Width() != width_))
        throw std::logic_error("BitMatrix::operator&=: non-conformant matrices");

    const unsigned int ldim_rhs = rhs.LDim();
    const unsigned int* buf_rhs = rhs.LockedBuffer();

    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int r_wd = 0;
        unsigned int col_offset = c*ldim_wds_;
        unsigned int col_offset_rhs = c*ldim_rhs;

        for (; r_wd < full_wds_; ++r_wd)
            words_[col_offset + r_wd] |= buf_rhs[col_offset_rhs + r_wd];

        if (MASK_ > 0)
            words_[col_offset + r_wd] |= (MASK_ & buf_rhs[col_offset_rhs + r_wd]);
    }

    return *this;
}

//-----------------------------------------------------------------------------
BitMatrix& BitMatrix::operator^=(const BitMatrix& rhs)
{
    if ( (rhs.Height() != height_) && (rhs.Width() != width_))
        throw std::logic_error("BitMatrix::operator&=: non-conformant matrices");

    const unsigned int ldim_rhs = rhs.LDim();
    const unsigned int* buf_rhs = rhs.LockedBuffer();

    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int r_wd = 0;
        unsigned int col_offset = c*ldim_wds_;
        unsigned int col_offset_rhs = c*ldim_rhs;

        for (; r_wd < full_wds_; ++r_wd)
            words_[col_offset + r_wd] ^= buf_rhs[col_offset_rhs + r_wd];

        if (MASK_ > 0)
            words_[col_offset + r_wd] ^= (MASK_ & buf_rhs[col_offset_rhs + r_wd]);
    }

    return *this;
} 

//-----------------------------------------------------------------------------
BitMatrix operator&(const BitMatrix& A, const BitMatrix& B)
{
    BitMatrix result(A);
    return (result &= B);
}

//-----------------------------------------------------------------------------
BitMatrix operator|(const BitMatrix& A, const BitMatrix& B)
{
    BitMatrix result(A);
    return (result |= B);
}

//-----------------------------------------------------------------------------
BitMatrix operator^(const BitMatrix& A, const BitMatrix& B)
{
    BitMatrix result(A);
    return (result ^= B);
}

//-----------------------------------------------------------------------------
BitMatrix& BitMatrix::operator=(const BitMatrix& rhs)
{
    if (this != &rhs)
        Assign(rhs);

    return *this;
}

//-----------------------------------------------------------------------------
BitMatrix::BitMatrix(const BitMatrix& rhs)
{
    Assign(rhs);
}

//-----------------------------------------------------------------------------
BitMatrix& BitMatrix::operator=(const unsigned int val)
{
    unsigned int fill_val = 
        val > 0 ? std::numeric_limits<unsigned int>::max() : 0;
    std::fill(words_.begin(), words_.end(), fill_val);
    return *this;
}

//-----------------------------------------------------------------------------
void BitMatrix::LinearIndexedAssign(const std::vector<unsigned int>& linear_indices,
                                    const unsigned int val)
{
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    for (unsigned int i=0; i != linear_indices.size(); ++i)
    {
        unsigned int linear_index = linear_indices[i];
        unsigned int c = linear_index / height_;
        assert(c < width_);

        unsigned int r = linear_index - c*height_;
        assert(r < height_);

        // find the word containing bit (r, c)
        unsigned int r_wd = r / BITS;
        unsigned int extra = r - BITS*r_wd;
        assert(r_wd <= ldim_wds_);

        unsigned int wd = words_[c*ldim_wds_ + r_wd];
        unsigned int bit_mask = (1 << extra);

        if (val > 0)
        {
            // turn on the bit
            wd |= bit_mask;
        }
        else
        {
            // turn off the bit
            bit_mask = ~bit_mask;       
            wd &= bit_mask;
        }
        
        words_[c*ldim_wds_ + r_wd] = wd;
    }
}

//-----------------------------------------------------------------------------
unsigned int BitMatrix::ColumnSum(const unsigned int c) const
{
    unsigned int col_offset = c*ldim_wds_;
    unsigned int sum=0, r_wd=0;
    for (; r_wd<full_wds_; ++r_wd)
        sum += PopulationCount(words_[col_offset + r_wd]);

    if (MASK_ > 0)
        sum += PopulationCount(MASK_ & words_[col_offset + r_wd]);

    return sum;
}

//-----------------------------------------------------------------------------
unsigned int BitMatrix::RowSum(const unsigned int r) const
{
    // find the word and bit offset for this row
    unsigned int r_wd = r / BITS_PER_WORD;
    unsigned int q = r - BITS_PER_WORD*r_wd;
    unsigned int MASK = (1 << q);

    unsigned int sum=0;
    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int col_offset = c*ldim_wds_;
        if (words_[col_offset + r_wd] & MASK)
            ++sum;
    }

    return sum;
}

//-----------------------------------------------------------------------------
void BitMatrix::SubmatrixFromCols(BitMatrix& result, const unsigned int c) const
{
    if (c >= width_)
        throw std::logic_error("BitMatrix::SubmatrixFromCols: invalid column index");

    result.Resize(height_, 1);
    assert(result.LDim() >= ldim_wds_);
    memcpy(&result.words_[0], &words_[c*ldim_wds_], height_ * sizeof(unsigned int));
}

//-----------------------------------------------------------------------------
void BitMatrix::SubmatrixFromCols(BitMatrix& result, 
                                  const std::vector<unsigned int>& col_indices) const
{
    const unsigned int new_width = col_indices.size();

    result.Resize(height_, new_width);
    unsigned int* buf_r = result.Buffer();
    const unsigned int ldim_r = result.LDim();
    assert(ldim_r >= ldim_wds_);

    unsigned int dest_indx = 0;
    for (unsigned int c=0; c<new_width; ++c)
    {
        unsigned int src_indx = col_indices[c];
        assert(src_indx < width_);
        assert(dest_indx < new_width);
        memcpy(&buf_r[dest_indx*ldim_r], 
               &words_[src_indx*ldim_wds_], 
               height_ * sizeof(unsigned int));
        ++dest_indx;
    }

    assert(new_width == dest_indx);
}

//-----------------------------------------------------------------------------
void BitMatrix::SubmatrixFromCols(BitMatrix& result, 
                                  const BitMatrix& mask) const
{
    // The mask must be a column vector with one bit per column.
    if ( (mask.Height() != width_) || (mask.Width() != 1))
        throw std::logic_error("BitMatrix::SubmatrixFromCols: invalid mask dimensions");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    unsigned int new_width = mask.ColumnSum(0);

    result.Resize(height_, new_width);
          unsigned int* buf_r = result.Buffer();
    const unsigned int ldim_r = result.LDim();

    assert(ldim_r >= ldim_wds_);

    const unsigned int* buf_m = mask.LockedBuffer();

    // scan the bits of mask; the index of each nonzero bit is a column to copy
    
    const unsigned int full_wds = mask.Height() / BITS;
    const unsigned int extra    = mask.Height() - full_wds*BITS;

    unsigned int r_wd = 0, wd=0, dest_col=0, c=0;
    for (; r_wd < full_wds; ++r_wd)
    {
        wd = buf_m[r_wd];
        for (unsigned int q=0; q<BITS; ++q, ++c)
        {
            if (wd & (1 << q))
            {
                assert(c < width_);
                assert(dest_col < new_width);
                memcpy(&buf_r[dest_col*ldim_r], 
                       &words_[c*ldim_wds_], 
                       height_ * sizeof(unsigned int));
                ++dest_col;
            }
        }
    }

    if (extra > 0)
    {
        wd = buf_m[r_wd];
        for (unsigned int q=0; q<extra; ++q, ++c)
        {
            if (wd & (1 << q))
            {
                assert(c < width_);
                assert(dest_col < new_width);
                memcpy(&buf_r[dest_col*ldim_r],
                       &words_[c*ldim_wds_],
                       height_ * sizeof(unsigned int));
                ++dest_col;
            }
        }
    }
    assert(new_width == dest_col);
}

//-----------------------------------------------------------------------------
unsigned int BitMatrix::MaxRowIndex(const unsigned int c) const
{
    if (c >= width_)
        throw std::logic_error("BitMatrix::MaxRowIndex: invalid column index");

    const unsigned int col_offset = c*ldim_wds_;
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;
    int r_wd_start = static_cast<int>(ldim_wds_ - 1);

    // scan backwards from the final word in column c
    if (MASK_ > 0)
    {
        unsigned int extra = height_ - BITS*full_wds_;
        unsigned int wd    = MASK_ & words_[col_offset + r_wd_start];
        for (int q=extra-1; q>=0; --q)
        {
            if (wd & (1 << q))
                return (full_wds_*BITS + q);
        }

        --r_wd_start;
    }

    for (int r_wd = r_wd_start; r_wd >= 0; --r_wd)
    {
        unsigned int wd = words_[col_offset + r_wd];
        for (int q=BITS-1; q>=0; --q)
        {
            if (wd & (1 << q))
            {
                if (r_wd > 0)
                    return (r_wd-1)*BITS + q;
                else
                    return q;
            }
        }
    }

    return 0;
}

//-----------------------------------------------------------------------------
void BitMatrix::ColIndices(const unsigned int r,
                           std::vector<unsigned int>& col_indices,
                           unsigned int& num_cols) const
{
    if (r >= height_)
        throw std::logic_error("BitMatrix::ColIndices: invalid row index");

    if (col_indices.size() < width_)
        col_indices.resize(width_);

    // find the word and bit that contains this row
    unsigned int r_wd = r/BITS_PER_WORD;
    unsigned int q = r - BITS_PER_WORD*r_wd;
    unsigned int MASK = (1 << q);

    num_cols = 0;
    for (unsigned int c=0; c<width_; ++c)
    {
        unsigned int col_offset = c*ldim_wds_;
        if (words_[col_offset + r_wd] & MASK)
            col_indices[num_cols++] = c;
    }

    // if (0 == num_cols)
    // {
    //     cout << "*** BitMatrix::RowIndices: found a zero row. ***" << endl;
    //     cout << "height_: " << height_ << ", width_: " << width_ << endl;
    //     cout << "r_wd: " << r_wd << endl;
    //     cout << "ldim_wds_: " << ldim_wds_ << endl;
    //     cout << "q: " << q << endl;
    //     cout << "MASK: " << std::hex << MASK_ << std::dec << endl;
    //     assert(0 != num_cols);        
    // }
}

//-----------------------------------------------------------------------------
void BitMatrix::RowIndices(const unsigned int c,
                           std::vector<unsigned int>& row_indices,
                           unsigned int& num_rows) const
{
    if (c >= width_)
        throw std::logic_error("BitMatrix::RowIndices: invalid column index");

    if (row_indices.size() < height_)
        row_indices.resize(height_);

    num_rows = 0;
    unsigned int offset = c*ldim_wds_, r_wd = 0, r=0;
    for (; r_wd < full_wds_; ++r_wd)
    {
        unsigned int wd = words_[offset + r_wd];
        for (unsigned int q=0; q<BITS_PER_WORD; ++q, ++r)
        {
            if (wd & (1 << q))
            {
                assert(r < height_);
                row_indices[num_rows++] = r;
            }
        }
    }

    if (MASK_ > 0)
    {
        unsigned int extra = height_ - BITS_PER_WORD*full_wds_;
        unsigned int wd = MASK_ & words_[offset + r_wd];
        for (unsigned int q=0; q<extra; ++q, ++r)
        {
            if (wd & (1 << q))
            {
                assert(r < height_);
                row_indices[num_rows++] = r;
            }
        }
    }

    // if (0 == num_rows)
    // {
    //     cout << "*** BitMatrix::RowIndices: found a zero column. ***" << endl;
    //     cout << "height_: " << height_ << ", width_: " << width_ << endl;
    //     cout << "MASK_: " << std::hex << MASK_ << std::dec << endl;
    //     cout << "full_wds_: " << full_wds_ << endl;
    //     cout << "extra: " << height_ - BITS_PER_WORD*full_wds_ << endl;
    //     cout << "wd: " << (MASK_ & words_[offset + r_wd]) << endl;
    //     cout << "column sum for col " << c << ": " << ColumnSum(c) << endl;

    //     unsigned int r_wd = 0;
    //     for (; r_wd<full_wds_; ++r_wd)
    //         cout << "words_[" << r_wd << "]: " << std::hex << words_[offset + r_wd] << std::dec << endl;
    //     if (MASK_ > 0)
    //     {
    //         unsigned int wd = MASK_ & words_[offset + r_wd];
    //         cout << "words_[" << r_wd << "]: " << std::hex << wd << std::dec << endl;
    //     }

    //     assert(0 != num_rows);
    // }
}

//-----------------------------------------------------------------------------
void BitMatrix::SetBits(const BitMatrix& mask)
{
    if ( (mask.Height() != height_) || (mask.Width() != width_))
        throw std::logic_error("BitMatrix::SetBits: mask size mismatch");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    unsigned int extra = height_ - BITS*full_wds_;
    for (unsigned int c=0; c != width_; ++c)
    {
        unsigned int wd, mask_wd, r_wd = 0u, r=0u, col_offset = c*ldim_wds_;

        // set the bits in each full word in column c
        for (; r_wd < full_wds_; ++r_wd)
        {
            wd      = words_[col_offset + r_wd];
            mask_wd = mask.words_[col_offset + r_wd];
            for (unsigned int q=0u; q<BITS; ++q, ++r)
            {
                if (mask_wd & (1 << q))
                    wd |= (1 << q);
            }

            words_[col_offset + r_wd] = wd;
        }

        // scan extra bits, if any, in the final word of column c
        if (extra > 0)
        {
            wd      = words_[col_offset + r_wd];
            mask_wd = mask.words_[col_offset + r_wd];
            for (unsigned int q=0u; q<extra; ++q, ++r)
            {
                if (mask_wd & (1 << q))
                    wd |= (1 << q);
            }

            words_[col_offset + r_wd] = wd;
        }
    }    
}

//-----------------------------------------------------------------------------
void BitMatrix::ClearBits(const BitMatrix& mask)
{
    if ( (mask.Height() != height_) || (mask.Width() != width_))
        throw std::logic_error("BitMatrix::ClearBits: mask size mismatch");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    unsigned int extra = height_ - BITS*full_wds_;
    for (unsigned int c=0; c != width_; ++c)
    {
        unsigned int wd, mask_wd, r_wd = 0u, r=0u, col_offset = c*ldim_wds_;

        // set the bits in each full word in column c
        for (; r_wd < full_wds_; ++r_wd)
        {
            wd      = words_[col_offset + r_wd];
            mask_wd = mask.words_[col_offset + r_wd];
            for (unsigned int q=0u; q<BITS; ++q, ++r)
            {
                if (mask_wd & (1 << q))
                    wd &= ~(1 << q);
            }

            words_[col_offset + r_wd] = wd;
        }

        // scan extra bits, if any, in the final word of column c
        if (extra > 0)
        {
            wd      = words_[col_offset + r_wd];
            mask_wd = mask.words_[col_offset + r_wd];
            for (unsigned int q=0u; q<extra; ++q, ++r)
            {
                if (mask_wd & (1 << q))
                    wd &= ~(1 << q);
            }

            words_[col_offset + r_wd] = wd;
        }
    }
}

//-----------------------------------------------------------------------------
void BitMatrix::Fill(const unsigned int val)
{
    // fill the buffer with the desired value
    std::fill(words_.begin(), words_.end(), val);

    // mask off the final word in each column
    for (unsigned int c=0; c<width_; ++c)
    {
        // final word in column c is at word index (ldim_wds_ - 1)
        words_[c*ldim_wds_ + (ldim_wds_ - 1)] &= MASK_;
    }
}

//-----------------------------------------------------------------------------
void BitMatrix::ToggleBit(const unsigned int row, const unsigned int col)
{
    if ( (row >= height_) || (col >= width_))
        throw std::logic_error("BitMatrix::ToggleBit: invalid indices");

    // get the column, the word, and the offset for the desired bit
    unsigned int col_offset = col * ldim_wds_;
    unsigned int r_wd = row / BITS_PER_WORD;
    unsigned int q = row - r_wd*BITS_PER_WORD;

    // isolate the bit
    unsigned int MASK = (1 << q);
    unsigned int wd = words_[col_offset + r_wd];
    unsigned int bitval = wd & MASK;

    if (bitval > 0)
    {
        // toggle bit off
        wd &= ~MASK;
    }
    else
    {
        // toggle bit on
        wd |= MASK;
    }

    words_[col_offset + r_wd] = wd;
}

//-----------------------------------------------------------------------------
void BitMatrix::Print() const
{
    for (unsigned int r=0; r<height_; ++r)
    {
        unsigned int r_wd = r/BITS_PER_WORD;
        unsigned int start_row = r_wd * BITS_PER_WORD;

        cout << "[" << r << "]: ";
        for (unsigned int c=0; c<width_; ++c)
        {
            unsigned int wd = words_[c*ldim_wds_ + r_wd];
            unsigned int val = wd & (1 << (r-start_row));
            cout << (val ? 1 : 0) << ", ";
        }

        cout << endl;
    }
}

//-----------------------------------------------------------------------------
void FindGroups(const unsigned int num_cols,
                const std::vector<uint64_t>& hashes,
                const std::vector<unsigned int>& indices,
                std::vector<unsigned int>& groups)
{
    unsigned int c1=0, c2;
    uint64_t hashval = hashes[indices[c1]];
    while (c1 < (num_cols-1))
    {
        c2 = c1 + 1;
        while (c2 < num_cols)
        {
            if (hashes[indices[c2]] == hashval)
                ++c2;
            else
            {
                // found a different hash at c2
                hashval = hashes[indices[c2]];
                break;
            }
        }
        
        // finished with this group, so save the size
        unsigned int equal_count = c2 - c1;
        assert(0 != equal_count);
        groups.push_back(equal_count);
        
        // Check for collisions? - TBD
        
        c1 = c2;
        if (c1 == (num_cols - 1))
        {
            // only a single unique column left
            groups.push_back(1);
            break;
        }
    }
}

//-----------------------------------------------------------------------------
void BitMatrix::GroupIdenticalColumns(std::vector<unsigned int>& indices,
                                      std::vector<unsigned int>& groups,
                                      const std::vector<unsigned int>& col_indices) const
{
    // This method finds all groups of identical columns in the BitMatrix.
    // The returned 'groups' array contains the number of columns in each
    // group.  The returned 'indices' array contains the column indices for
    // each group.  To illustrate, let
    //
    //        groups  = {1, 2, 1, 3}
    //        indices = {6, 0, 1, 4, 2, 3, 5}
    //
    // In this example, the elements of the 'groups' array sum to seven, so 
    // this example matrix has seven columns. The groups can be read from
    // the indices array as follows:
    //
    //        group 0: {column 6}
    //        group 1: {column 0, column 1}
    //        group 2: {column 4}
    //        group 3: {column 2, column 3, column 5}
    //
    // Thus columns 0 and 1 are identical, and columns 2, 3, and 5 are also.
    // Columns 4 and 6 are unique.

    const unsigned int num_cols = col_indices.size();

    if (indices.size() < num_cols)
        indices.resize(num_cols);

    groups.clear();

    // special case of a single-column matrix
    if (1 == num_cols)
    {
        indices[0] = 0;
        groups.push_back(1);
        return;
    }

    SpookyHash spooky;
    const uint64_t SPOOKY_HASH_SEED = 42u;
    std::vector<uint64_t> hashes(num_cols);

    // hash the columns with the 64-bit variant of SpookyHash
    for (unsigned int c=0; c != num_cols; ++c)
    {
        unsigned int col_index = col_indices[c];
        uint64_t hashval = spooky.Hash64(&words_[col_index*ldim_wds_],
                                         ldim_wds_ * sizeof(unsigned int),
                                         SPOOKY_HASH_SEED);

        hashes[c] = hashval;
    }

    // find groups of identical hashes by sorting the hash values

    // first setup the sort indices (number from 0 for application to submatrices)
    for (unsigned int c=0; c<num_cols; ++c)
        indices[c] = c;

    // Rearrange the indices by comparing hash values.  IMPORTANT: this gives
    // an 'alloc_dealloc_mismatch' error in AddressSanitizer; disable this 
    // error with this: export ASAN_OPTIONS=alloc_dealloc_mismatch=0 if this
    // is the only such error.
    std::stable_sort(indices.begin(), indices.begin() + num_cols,
                     [&hashes](const unsigned int i1, const unsigned int i2)
                     {
                         return hashes[i1] < hashes[i2];
                     });

    // walk the (conceptually sorted) hash array and compute group sizes
    FindGroups(num_cols, hashes, indices, groups);

    // cout << "GroupIdenticalColumns: found " << groups.size() << " hash groups." << endl;
    // cout << "MASK_: " << std::hex << MASK_ << std::dec << endl;

    // unsigned int col_count = 0;
    // for (unsigned int q=0; q<groups.size(); ++q)
    //     col_count += groups[q];
    // cout << "GroupIdenticalColumns: column count: " << col_count << endl;
    // assert(col_count == width_);

    // unsigned int count = 0;
    // for (unsigned int g=0; g<groups.size(); ++g)
    // {
    //     unsigned int group_size = groups[g];
    //     cout << "\tgroup[" << g << "] size: " << group_size << endl;
    //     for (unsigned int c=count; c<(count + group_size); ++c)
    //         cout << "hashes[" << c << "] (col " << indices[c] << "):\t" << hashes[indices[c]] << endl;
    //     cout << endl;
    //     count += group_size;
    // }
}

//-----------------------------------------------------------------------------
void BitMatrix::GroupIdenticalColumns(std::vector<unsigned int>& indices,
                                      std::vector<unsigned int>& groups) const
{
    // This method finds all groups of identical columns in the BitMatrix.
    // The returned 'groups' array contains the number of columns in each
    // group.  The returned 'indices' array contains the column indices for
    // each group.  To illustrate, let
    //
    //        groups  = {1, 2, 1, 3}
    //        indices = {6, 0, 1, 4, 2, 3, 5}
    //
    // In this example, the elements of the 'groups' array sum to seven, so 
    // this example matrix has seven columns. The groups can be read from
    // the indices array as follows:
    //
    //        group 0: {column 6}
    //        group 1: {column 0, column 1}
    //        group 2: {column 4}
    //        group 3: {column 2, column 3, column 5}
    //
    // Thus columns 0 and 1 are identical, and columns 2, 3, and 5 are also.
    // Columns 4 and 6 are unique.

    if (indices.size() < width_)
        indices.resize(width_);

    groups.clear();

    // special case of a single-column matrix
    if (1 == width_)
    {
        indices[0] = 0;
        groups.push_back(1);
        return;
    }

    SpookyHash spooky;
    const uint64_t SPOOKY_HASH_SEED = 42u;
    std::vector<uint64_t> hashes(width_);

    // hash the columns with the 64-bit variant of SpookyHash
    for (unsigned int c=0; c != width_; ++c)
    {
        uint64_t hashval = spooky.Hash64(&words_[c*ldim_wds_],
                                         ldim_wds_ * sizeof(unsigned int),
                                         SPOOKY_HASH_SEED);

        hashes[c] = hashval;
    }

    // find groups of identical hashes by sorting the hash values

    // first setup the sort indices
    for (unsigned int c=0; c<width_; ++c)
        indices[c] = c;

    // Rearrange the indices by comparing hash values.  IMPORTANT: this gives
    // an 'alloc_dealloc_mismatch' error in AddressSanitizer; disable this 
    // error with this: export ASAN_OPTIONS=alloc_dealloc_mismatch=0 if this
    // is the only such error.
    std::stable_sort(indices.begin(), indices.begin() + width_,
                     [&hashes](const unsigned int i1, const unsigned int i2)
                     {
                         return hashes[i1] < hashes[i2];
                     });

    // walk the (conceptually sorted) hash array and compute group sizes
    FindGroups(width_, hashes, indices, groups);
}
