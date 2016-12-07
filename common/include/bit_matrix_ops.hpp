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
#include <stdexcept>
#include <functional>
#include "bit_matrix.hpp"
#include "dense_matrix.hpp"
#include "openmp_pragma.hpp"

//-----------------------------------------------------------------------------
template <typename T, typename BinaryPredicate>
BitMatrix Apply(const DenseMatrix<T>& M, const T val, BinaryPredicate f)
{
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    unsigned int height = M.Height();
    unsigned int width  = M.Width();
    BitMatrix result(height, width);
    
    unsigned int* buf_r       = result.Buffer();
    const unsigned int ldim_r = result.LDim();
    const unsigned int MASK   = result.Mask();

    const T* buf_m = M.LockedBuffer();
    const unsigned int ldim_m = M.LDim();

    const unsigned int full_wds = height / BITS;
    const unsigned int extra    = height - BITS*full_wds;
    
    OPENMP_PRAGMA(omp parallel for)
    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int offset_r = c*ldim_r;
        unsigned int offset_m = c*ldim_m;

        unsigned int r_wd=0, r=0;
        for (; r_wd < full_wds; ++r_wd)
        {
            unsigned int wd = 0u;
            for (unsigned int q=0; q<BITS; ++q, ++r)
            {
                if ( f(buf_m[offset_m + r], val))
                    wd |= (1 << q);
            }

            buf_r[offset_r + r_wd] = wd;
        }

        if (extra > 0)
        {
            unsigned int wd = 0u;
            for (unsigned int q=0; q<extra; ++q, ++r)
            {
                if ( f(buf_m[offset_m + r], val))
                    wd |= (1 << q);
            }

            buf_r[offset_r + r_wd] = MASK & wd;
        }
    }

    return result;
}

//-----------------------------------------------------------------------------
template <typename T, typename BinaryPredicate>
BitMatrix Apply(const std::vector<T>& V, const T val, BinaryPredicate f)
{
    const unsigned int BITS = BitMatrix::BITS_PER_WORD;
    const unsigned int size = V.size();

    BitMatrix result(size);
    unsigned int* buf_r = result.Buffer();
    
    unsigned int full_wds = size / BITS;
    unsigned int extra = size - BITS * full_wds;

    unsigned int r_wd = 0, wd, r=0;
    for (; r_wd < full_wds; ++r_wd)
    {
        wd = 0;
        for (unsigned int q=0; q<BITS; ++q, ++r)
        {
            if ( f(V[r], val))
                wd |= (1 << q);
        }

        buf_r[r_wd] = wd;
    }

    if (extra > 0)
    {
        wd = 0;
        for (unsigned int q=0; q<extra; ++q, ++r)
        {
            if ( f(V[r], val))
                wd |= (1 << q);
        }

        buf_r[r_wd] = wd;
    }

    return result;   
}

//-----------------------------------------------------------------------------
template <typename T>
DenseMatrix<T>& CompoundAND(DenseMatrix<T>& D, const BitMatrix& B)
{
    // Set elements of D to zero wherever the BitMatrix contains a zero bit.

    const unsigned int height = D.Height();
    const unsigned int width  = D.Width();

    if (B.Height() != height)
         throw std::logic_error("operator&=: matrix and mask have different heights");
    if (B.Width() != width)
        throw std::logic_error("operator&= matrix and mask have different widths");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    const unsigned int* buf_b  = B.LockedBuffer();
    const unsigned int  ldim_b = B.LDim();
    const unsigned int full_wds = height / BITS;
    const unsigned int extra    = height - full_wds*BITS;

    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int offset_b = c*ldim_b;

        unsigned int r_wd = 0, r=0;
        for (; r_wd < full_wds; ++r_wd)
        {
            unsigned int wd = buf_b[offset_b + r_wd];
            for (unsigned int q=0; q<BITS; ++q, ++r)
            {
                if (0 == (wd & (1 << q)))
                    D.Set(r, c, T(0.0));
            }
        }

        if (extra > 0)
        {
            unsigned int wd = buf_b[offset_b + r_wd];
            for (unsigned int q=0; q<extra; ++q, ++r)
            {
                if (0 == (wd & (1 << q)))
                    D.Set(r, c, T(0.0));
            }
        }
    }

    return D;
}

//-----------------------------------------------------------------------------
template <typename T, typename BinaryOp>
BitMatrix VectorCompare(const std::vector<T>& A, 
                        const std::vector<T>& B,
                        BinaryOp f)
{
    const int size = A.size();
    if (B.size() != A.size())
        throw std::logic_error("VectorCompare: vectors must have same size");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;

    BitMatrix result(size);
    unsigned int* buf_r = result.Buffer();

    unsigned int full_wds = size / BITS;
    unsigned int extra = size - BITS * full_wds;

    unsigned int r_wd = 0, wd, r=0;
    for (; r_wd < full_wds; ++r_wd)
    {
        wd = 0;
        for (unsigned int q=0; q<BITS; ++q, ++r)
        {
            if ( f(A[r], B[r]))
                wd |= (1 << q);
        }

        buf_r[r_wd] = wd;
    }

    if (extra > 0)
    {
        wd = 0;
        for (unsigned int q=0; q<extra; ++q, ++r)
        {
            if ( f(A[r], B[r]))
                wd |= (1 << q);
        }

        buf_r[r_wd] = wd;
    }

    return result;
}

//-----------------------------------------------------------------------------
template <typename T>
void IndexedUpdate(std::vector<T>& A, const BitMatrix& B, const T val)
{
    // A(B) += val; update indexed elements of A with val

    const unsigned int size = A.size();
    if (0u == size)
        return;

    if ( (B.Height() != size) || (B.Width() != 1))
        throw std::logic_error("IndexedAssign: BitMatrix must be a column vector");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;
    const unsigned int* buf_b = B.LockedBuffer();

    const unsigned int full_wds = size / BITS;
    const unsigned int extra    = size - full_wds * BITS;

    unsigned int wd, r_wd = 0, r=0;
    for (; r_wd < full_wds; ++r_wd)
    {
        wd = buf_b[r_wd];
        for (unsigned int q=0; q<BITS; ++q, ++r)
        {
            if (wd & (1 << q))
                A[r] += val;
        }
    }

    if (extra > 0)
    {
        wd = buf_b[r_wd];
        for (unsigned int q=0; q<extra; ++q, ++r)
        {
            if (wd & (1 << q))
                A[r] += val;
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void IndexedAssign(std::vector<T>& A, const BitMatrix& B, const T val)
{
    // assign val to elements of A using the mask B

    const unsigned int size = A.size();
    if (0u == size)
        return;

    if ( (B.Height() != size) || (B.Width() != 1))
        throw std::logic_error("IndexedAssign: BitMatrix must be a column vector");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;
    const unsigned int* buf_b = B.LockedBuffer();

    const unsigned int full_wds = size / BITS;
    const unsigned int extra    = size - full_wds * BITS;

    unsigned int wd, r_wd = 0, r=0;
    for (; r_wd < full_wds; ++r_wd)
    {
        wd = buf_b[r_wd];
        for (unsigned int q=0; q<BITS; ++q, ++r)
        {
            if (wd & (1 << q))
                A[r] = val;
        }
    }

    if (extra > 0)
    {
        wd = buf_b[r_wd];
        for (unsigned int q=0; q<extra; ++q, ++r)
        {
            if (wd & (1 << q))
                A[r] = val;
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void IndexedAssign(std::vector<T>& A, 
                   const BitMatrix& B, 
                   const std::vector<T>& C)
{
    const unsigned int size = A.size();
    if (0u == size)
        return;
    
    if (C.size() != size)
        throw std::logic_error("IndexedAssign: vectors must have same size");
    if ( (B.Height() != size) || (B.Width() != 1))
        throw std::logic_error("IndexedAssign: BitMatrix must be a column vector");

    const unsigned int BITS = BitMatrix::BITS_PER_WORD;
    const unsigned int* buf_b = B.LockedBuffer();

    const unsigned int full_wds = size / BITS;
    const unsigned int extra    = size - full_wds * BITS;

    unsigned int wd, r_wd = 0, r=0;
    for (; r_wd < full_wds; ++r_wd)
    {
        wd = buf_b[r_wd];
        for (unsigned int q=0; q<BITS; ++q, ++r)
        {
            if (wd & (1 << q))
                A[r] = C[r];
        }
    }

    if (extra > 0)
    {
        wd = buf_b[r_wd];
        for (unsigned int q=0; q<extra; ++q, ++r)
        {
            if (wd & (1 << q))
                A[r] = C[r];
        }
    }
}

// returns true if A contains no data
bool IsEmpty(const BitMatrix& A);

// returns true if all bits in A or in the specified cols of A are set to 1
bool All(const BitMatrix& A);
bool AllCols(const BitMatrix& A, const std::vector<unsigned int>& col_indices);
bool AllRows(const BitMatrix& A, const std::vector<unsigned int>& row_indices);

// equivalent to this Matlab operation: B & repmat(mask, B.Height(), 1)
BitMatrix ColumnwiseAND(const BitMatrix& B, const BitMatrix& mask);

// generate a new matrix of dimension height x mask.Height(); each bit of
// mask is replicated down each column
BitMatrix MatrixFromColumnMask(const BitMatrix& mask, const unsigned int height);

// bitwise operations with BitMatrix
BitMatrix operator&(const BitMatrix& A, const BitMatrix& B);
BitMatrix operator|(const BitMatrix& A, const BitMatrix& B);
BitMatrix operator^(const BitMatrix& A, const BitMatrix& B);

// set elements of DenseMatrix D to zero wherever BitMatrix B has a zero bit
template <typename T>
DenseMatrix<T>& operator&=(DenseMatrix<T>& D, const BitMatrix& B)
{return CompoundAND(D, B);}

template <typename T>
DenseMatrix<T> operator&(const DenseMatrix<T>& D, const BitMatrix& B)
{
    DenseMatrix<T> result(D);
    return CompoundAND(result, B);
}

// logical operations with DenseMatrix
template <typename T>
BitMatrix operator < (const DenseMatrix<T>& M, const T val)
{return Apply(M, val, std::less<T>());}

template <typename T>
BitMatrix operator <= (const DenseMatrix<T>& M, const T val)
{return Apply(M, val, std::less_equal<T>());}

template <typename T>
BitMatrix operator > (const DenseMatrix<T>& M, const T val)
{return Apply(M, val, std::greater<T>());}

template <typename T>
BitMatrix operator >= (const DenseMatrix<T>& M, const T val)
{return Apply(M, val, std::greater_equal<T>());}

template <typename T>
BitMatrix operator == (const DenseMatrix<T>& M, const T val)
{return Apply(M, val, std::equal_to<T>());}

template <typename T>
BitMatrix operator != (const DenseMatrix<T>& M, const T val)
{return Apply(M, val, std::not_equal_to<T>());}

// logical operations with std::vector
template <typename T>
BitMatrix operator < (const std::vector<T>& V, const T val)
{return Apply(V, val, std::less<T>());}

template <typename T>
BitMatrix operator <= (const std::vector<T>& V, const T val)
{return Apply(V, val, std::less_equal<T>());}

template <typename T>
BitMatrix operator > (const std::vector<T>& V, const T val)
{return Apply(V, val, std::greater<T>());}

template <typename T>
BitMatrix operator >= (const std::vector<T>& V, const T val)
{return Apply(V, val, std::greater_equal<T>());}

template <typename T>
BitMatrix operator == (const std::vector<T>& V, const T val)
{return Apply(V, val, std::equal_to<T>());}

template <typename T>
BitMatrix operator != (const std::vector<T>& V, const T val)
{return Apply(V, val, std::not_equal_to<T>());}

// operations on two std::vector types that yield a logical result
template <typename T>
BitMatrix operator < (const std::vector<T>& A, const std::vector<T>& B)
{return VectorCompare(A, B, std::less<T>());}

template <typename T>
BitMatrix operator >= (const std::vector<T>& A, const std::vector<T>& B)
{return VectorCompare(A, B, std::greater_equal<T>());}

