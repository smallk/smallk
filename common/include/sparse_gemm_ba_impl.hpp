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

#include <thread>
#include <mutex>
#include <vector>
#include <algorithm>
#include <cassert>

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmBA(const unsigned int start_col,
                      const unsigned int num_threads,
                      const T alpha,
                      const DenseMatrix<T>& B,
                      const SparseMatrix<T>& A,
                      const T beta,

                      // For some reason, passing a reference to matrix C
                      // causes a "double free" error from the thread library
                      // building with g++-4.7 and g++-4.8 in debug mode.  
                      // I don't know the reason for this (yet), so instead 
                      // of declaring the next parameter as 
                      //
                      //       "DenseMatrix<T>& C",
                      //
                      // which is what I would prefer to use, I'll use the 
                      // following instead...)
                      const unsigned int height_c,
                      T* data_c,
                      const unsigned int ldim_c)
{

    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    //T*                  data_c = C.Buffer();  // this line causes the problem

    const unsigned int ldim_b = B.LDim();
    //const unsigned int ldim_c = C.LDim();

    // for column j of the result (A and C have the same width)
    unsigned int j = start_col;
    while (j < width_a)
    {
        // offset to jth column of C
        unsigned int offset_c = j*ldim_c;

        // compute C(:,j) = beta*C(:,j); B and C have the same height
        if (T(0) == beta)
        {
            for (unsigned int r=0; r != height_b; ++r)
                data_c[offset_c + r] = T(0);
        }
        else
        {
            for (unsigned int r=0; r != height_b; ++r)
                data_c[offset_c + r] *= beta;
        }

        // bounds on data in the jth column of A 
        unsigned int start = cols_a[j];
        unsigned int end   = cols_a[j+1];
        for (unsigned int offset = start; offset != end; ++offset)
        {
            // the next nonzero row index in column j of A
            unsigned int row = rows_a[offset];

            // alpha*A(row, j)
            T alpha_A_rj = alpha*data_a[offset];

            // scale B(:, row) by alpha*A(row, j) and add to C(:,j)
            unsigned int offset_b = row*ldim_b;
            for (unsigned int r=0; r != height_b; ++r)
                data_c[offset_c + r] += (alpha_A_rj * data_b[offset_b + r]);
        }

        j += num_threads;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmBtA(const unsigned int start_col,
                       const unsigned int num_threads,
                       const T alpha,
                       const DenseMatrix<T>& B,
                       const SparseMatrix<T>& A,
                       const T beta,
                       const unsigned int height_c,
                       T* data_c,
                       const unsigned int ldim_c)
{

    const unsigned int width_a  = A.Width();
    const unsigned int width_b  = B.Width();

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();
    const T*            data_b = B.LockedBuffer();
    const unsigned int  ldim_b = B.LDim();

    // for column j of the result (A and C have the same width)
    unsigned int j = start_col;
    while (j < width_a)
    {
        // offset to jth column of C
        unsigned int offset_c = j*ldim_c;

        // compute C(:,j) = beta*C(:,j)
        if (T(0) == beta)
        {
            // Only need a write if beta is precisely zero.
            for (unsigned int r=0; r != height_c; ++r)
                data_c[offset_c + r] = T(0);
        }
        else
        {
            // must read, multiply, then write
            for (unsigned int r=0; r != height_c; ++r)
                data_c[offset_c + r] *= beta;
        }

        // bounds on data in the jth column of A 
        unsigned int start = cols_a[j];
        unsigned int end   = cols_a[j+1];
        for (unsigned int offset = start; offset != end; ++offset)
        {
            // the next nonzero row index in column j of A
            unsigned int row = rows_a[offset];

            // alpha*A(row, j)
            T alpha_A_rj = alpha*data_a[offset];

            // scale B'(:, row) == B(row, :) by alpha*A(row, j) and add to C(:,j)
            for (unsigned int c=0; c != width_b; ++c)
            {
                unsigned int offset_b = c*ldim_b;

                assert( (offset_c + c) < (ldim_c*A.Width()) );
                data_c[offset_c + c] += (alpha_A_rj * data_b[offset_b + row]);
            }
        }

        j += num_threads;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmRank2BA(const unsigned int j, // column index
                           const T alpha,
                           const DenseMatrix<T>& B,
                           const SparseMatrix<T>& A,
                           const T beta,
                           const unsigned int height_c,
                           T* data_c,
                           const unsigned int ldim_c)
{
    // matrices A and C each have two columns
    assert(j < 2);

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const unsigned int  ldim_b  = B.LDim();
    const T*            data_b  = B.LockedBuffer();

    // compute offset to jth column of C
    unsigned int offset_c = j*ldim_c;
    
    // bounds on data in the jth column of A
    unsigned int start = cols_a[j];
    unsigned int end   = cols_a[j+1];

    // loop over nonzeros in column j of A
    for (unsigned int offset=start; offset != end; ++offset)
    {
        // the next nonzero row index in column j of A
        unsigned int row = rows_a[offset];

        // alpha*A(row, j)
        T alpha_A_rj = alpha*data_a[offset];

        // scale B(:,row) by alpha*A(row, j) and add to C(:, j)
        unsigned int offset_b = row*ldim_b;
        for (unsigned int r=0; r != height_c; ++r)
            data_c[offset_c + r] += (alpha_A_rj * data_b[offset_b + r]);
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmRank2BtA(const unsigned int j, // column index
                            const T alpha,
                            const DenseMatrix<T>& B,
                            const SparseMatrix<T>& A,
                            const T beta,
                            const unsigned int height_c,
                            T* data_c,
                            const unsigned int ldim_c)
{
    assert(j < 2);
    assert(2 == A.Width());

    const unsigned int width_b = B.Width();

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    const unsigned int  ldim_b = B.LDim();

    // compute offset to jth column of C
    unsigned int offset_c = j*ldim_c;
    
    // bounds on data in the jth column of A 
    unsigned int start = cols_a[j];
    unsigned int end   = cols_a[j+1];
    for (unsigned int offset = start; offset != end; ++offset)
    {
        // the next nonzero row index in column j of A
        unsigned int row = rows_a[offset];

        // alpha*A(row, j)
        T alpha_A_rj = alpha*data_a[offset];

        // scale B'(:, row) == B(row, :) by alpha*A(row, j) and add to C(:,j)
        for (unsigned int c=0; c != width_b; ++c)
        {
            unsigned int offset_b = c*ldim_b;
            data_c[offset_c + c] += (alpha_A_rj * data_b[offset_b + row]);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmRank2BA(const T alpha, 
                 const DenseMatrix<T>& B, 
                 const SparseMatrix<T>& A, 
                 const T beta,  
                 DenseMatrix<T>& C,
                 const unsigned int num_threads)
{
    assert(2 == num_threads);
    assert(2 == C.Width());
    assert(2 == A.Width());

    // use sequential code to handle tiny matrices
    if (1 == A.Width())
        return GemmSequentialBA(alpha, B, A, beta, C);

    // Scale C by beta as a dense operation.  Timing experiments show that
    // this is generally faster for a two-column matrix than having each 
    // thread traverse the columns individually, as is done in
    // gemm_rank2_thread_func.
    Scal(beta, C);

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmRank2BA<T>, 
                                      k,
                                      alpha, std::cref(B), std::cref(A),
                                      beta,
                                      C.Height(),
                                      C.Buffer(),
                                      C.LDim()));
    }

    // wait for threads to finish
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmRank2BtA(const T alpha, 
                  const DenseMatrix<T>& B, 
                  const SparseMatrix<T>& A, 
                  const T beta,  
                  DenseMatrix<T>& C,
                  const unsigned int num_threads)
{
    assert(2 == num_threads);
    assert(2 == C.Width());
    assert(2 == A.Width());

    // use sequential code to handle tiny matrices
    if (1 == A.Width())
        return GemmSequentialBtA(alpha, B, A, beta, C);

    // Scale C by beta as a dense operation.  Timing experiments show that
    // this is generally faster for a two-column matrix than having each 
    // thread traverse the columns individually, as is done in
    // gemm_rank2_thread_func.
    Scal(beta, C);

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmRank2BtA<T>, 
                                      k,
                                      alpha, std::cref(B), std::cref(A),
                                      beta,
                                      C.Height(),
                                      C.Buffer(),
                                      C.LDim()));
    }

    // wait for threads to finish
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmSequentialBA(const T alpha, 
                      const DenseMatrix<T>& B, 
                      const SparseMatrix<T>& A, 
                      const T beta, 
                      DenseMatrix<T>& C)
{
    // compute C = alpha*B*A + beta*C, A is sparse, B and C are dense

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (height_a != width_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (width_a != width_c)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (height_b != height_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    T*                  data_c = C.Buffer();

    const unsigned int ldim_b = B.LDim();
    const unsigned int ldim_c = C.LDim();

    for (unsigned int j=0; j<width_c; ++j)
    {
        // offset to jth column of C
        unsigned int offset_c = j*ldim_c;

        // compute C(:,j) = beta*C(:,j)
        if (T(0) == beta)
        {
            for (unsigned int r=0; r != height_c; ++r)
                data_c[offset_c + r] = T(0);
        }
        else
        {
            for (unsigned int r=0; r != height_c; ++r)
                data_c[offset_c + r] *= beta;
        }

        // bounds on data in the jth column of A 
        unsigned int start = cols_a[j];
        unsigned int end   = cols_a[j+1];
        for (unsigned int offset = start; offset != end; ++offset)
        {
            // the next nonzero row index in column j of A
            unsigned int row = rows_a[offset];

            // alpha*A(row, j)
            T alpha_A_rj = alpha*data_a[offset];

            // scale B(:, row) by alpha*A(row, j) and add to C(:,j)
            unsigned int offset_b = row*ldim_b;
            for (unsigned int r=0; r != height_b; ++r)
                data_c[offset_c + r] += (alpha_A_rj * data_b[offset_b + r]);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmSequentialBtA(const T alpha, 
                       const DenseMatrix<T>& B, 
                       const SparseMatrix<T>& A, 
                       const T beta, 
                       DenseMatrix<T>& C)
{
    // compute C = alpha*B'*A + beta*C, A is sparse, B and C are dense

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (height_a != height_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (width_a != width_c)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (width_b != height_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    T*                  data_c = C.Buffer();

    const unsigned int ldim_b = B.LDim();
    const unsigned int ldim_c = C.LDim();

    for (unsigned int j=0; j<width_c; ++j)
    {
        // offset to jth column of C
        unsigned int offset_c = j*ldim_c;

        // compute C(:,j) = beta*C(:,j)
        if (T(0) == beta)
        {
            for (unsigned int r=0; r != height_c; ++r)
                data_c[offset_c + r] = T(0);
        }
        else
        {
            for (unsigned int r=0; r != height_c; ++r)
                data_c[offset_c + r] *= beta;
        }

        // bounds on data in the jth column of A 
        unsigned int start = cols_a[j];
        unsigned int end   = cols_a[j+1];
        for (unsigned int offset = start; offset != end; ++offset)
        {
            // the next nonzero row index in column j of A
            unsigned int row = rows_a[offset];

            // alpha*A(row, j)
            T alpha_A_rj = alpha*data_a[offset];

            // scale B'(:, row) == B(row, :) by alpha*A(row, j) and add to C(:,j)
            for (unsigned int c=0; c != width_b; ++c)
            {
                unsigned int offset_b = c*ldim_b;
                data_c[offset_c + c] += (alpha_A_rj * data_b[offset_b + row]);
            }
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmBA(const T alpha, 
            const DenseMatrix<T>& B, 
            const SparseMatrix<T>& A, 
            const T beta,  
            DenseMatrix<T>& C,
            const unsigned int max_threads)
{
    // C = alpha*B*A + beta*C

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (height_a != width_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (width_a != width_c)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (height_b != height_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");
    
    if (1 == max_threads)
        return GemmSequentialBA(alpha, B, A, beta, C);

    // handle rank 2 as a special case
    if (2 == width_c)
        return GemmRank2BA(alpha, B, A, beta, C, 2);

    // do not use more threads than columns of the result (C)
    unsigned int num_threads = max_threads;
    if (width_c < num_threads)
        num_threads = width_c;

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmBA<T>,
                                      k, num_threads,
                                      alpha, std::cref(B), std::cref(A),
                                      beta,
                                      C.Height(),
                                      C.Buffer(),
                                      C.LDim()));
    }

    // wait for threads to finish
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmBtA(const T alpha, 
             const DenseMatrix<T>& B, 
             const SparseMatrix<T>& A, 
             const T beta,  
             DenseMatrix<T>& C,
             const unsigned int max_threads)
{
    // C = alpha*B'A + beta*C

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (height_a != height_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (width_a != width_c)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (width_b != height_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");
    
    if (1 == max_threads)
        return GemmSequentialBtA(alpha, B, A, beta, C);

    // handle rank 2 as a special case
    if (2 == width_c)
        return GemmRank2BtA(alpha, B, A, beta, C, 2);

    // do not use more threads than columns of the result (C)
    unsigned int num_threads = max_threads;
    if (width_c < num_threads)
        num_threads = width_c;

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmBtA<T>,
                                      k, num_threads,
                                      alpha, std::cref(B), std::cref(A),
                                      beta,
                                      C.Height(),
                                      C.Buffer(),
                                      C.LDim()));
    }

    // wait for threads to finish
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}

