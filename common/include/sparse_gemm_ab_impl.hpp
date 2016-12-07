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

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmAB(const unsigned int start_col,
                      const unsigned int num_threads,
                      const T alpha,
                      const SparseMatrix<T>& A,
                      const DenseMatrix<T>& B,
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
    const unsigned int width_b  = B.Width();
    //const unsigned int height_c = C.Height();

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    //T*                  data_c = C.Buffer();  // this line causes the problem

    const unsigned int ldim_b = B.LDim();
    //const unsigned int ldim_c = C.LDim();

    unsigned int j = start_col;
    while (j < width_b)
    {
        // compute offset to jth column of B and C
        unsigned int offset_b = j*ldim_b;
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

        // C(:,j) += alpha*A*B(:,j)
        
        // loop over columns of A
        for (unsigned int c=0; c != width_a; ++c)
        {
            unsigned int start = cols_a[c];
            unsigned int end   = cols_a[c+1];

            // alpha*B(c,j); this value multiplies the entire column c of A
            T alpha_b_cj = alpha*data_b[offset_b + c];

            for (unsigned int offset=start; offset != end; ++offset)
            {
                // the next nonzero row in column j of A
                unsigned int row = rows_a[offset];
                data_c[offset_c + row] += alpha_b_cj*data_a[offset];
            }
        }

        j += num_threads;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmABt(const unsigned int start_col,
                       const unsigned int num_threads,
                       const T alpha,
                       const SparseMatrix<T>& A,
                       const DenseMatrix<T>& B,
                       const T beta,
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
    const unsigned int  ldim_b = B.LDim();

    unsigned int j = start_col;
    while (j < height_b)
    {
        // compute offset to jth column of C
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

        // C(:,j) += alpha*A*B'(:,j), or
        // C(:,j) += alpha*A*B(j, :)
        
        // loop over columns of A
        for (unsigned int c=0; c != width_a; ++c)
        {
            unsigned int start = cols_a[c];
            unsigned int end   = cols_a[c+1];

            // alpha*B'(c,j) == alpha*B(j, c)
            T alpha_b_cj = alpha*data_b[c*ldim_b + j];

            for (unsigned int offset=start; offset != end; ++offset)
            {
                // the next nonzero row in column j of A
                unsigned int row = rows_a[offset];
                data_c[offset_c + row] += alpha_b_cj*data_a[offset];
            }
        }

        j += num_threads;
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmRank2AB(const unsigned int j, // column index
                           const unsigned int num_threads,
                           const T alpha,
                           const SparseMatrix<T>& A,
                           const DenseMatrix<T>& B,
                           const T beta,
                           const unsigned int height_c,
                           T* data_c,
                           const unsigned int ldim_c)
{

    const unsigned int width_a  = A.Width();
    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    const unsigned int  ldim_b = B.LDim();

    // compute offset to jth column of B and C
    unsigned int offset_b = j*ldim_b;
    unsigned int offset_c = j*ldim_c;
    
    // C(:,j) += alpha*A*B(:,j)
    
    // loop over columns of A
    for (unsigned int c=0; c != width_a; ++c)
    {
        unsigned int start = cols_a[c];
        unsigned int end   = cols_a[c+1];
        
        // alpha*B(c,j); this value multiplies the entire column c of A
        T alpha_b_cj = alpha*data_b[offset_b + c];
        
        for (unsigned int offset=start; offset != end; ++offset)
        {
            // the next nonzero row in column j of A
            unsigned int row = rows_a[offset];
            data_c[offset_c + row] += alpha_b_cj*data_a[offset];
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ThreadFuncGemmRank2ABt(const unsigned int j, // column index
                            const unsigned int num_threads,
                            const T alpha,
                            const SparseMatrix<T>& A,
                            const DenseMatrix<T>& B,
                            const T beta,
                            const unsigned int height_c,
                            T* data_c,
                            const unsigned int ldim_c)
{

    const unsigned int width_a  = A.Width();
    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    const unsigned int  ldim_b = B.LDim();

    // compute offset to jth column of C
    unsigned int offset_c = j*ldim_c;
    
    // C(:,j) += alpha*A*B'(:,j), or
    // C(:,j) += alpha*A*B(j,:)
    
    // loop over columns of A
    for (unsigned int c=0; c != width_a; ++c)
    {
        unsigned int start = cols_a[c];
        unsigned int end   = cols_a[c+1];
        
        // alpha*B(c,j) == alpha*B(j, c)
        T alpha_b_cj = alpha*data_b[c*ldim_b + j];
        
        for (unsigned int offset=start; offset != end; ++offset)
        {
            // the next nonzero row in column j of A
            unsigned int row = rows_a[offset];
            data_c[offset_c + row] += alpha_b_cj*data_a[offset];
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmRank2AB(const T alpha, 
                 const SparseMatrix<T>& A, 
                 const DenseMatrix<T>& B, 
                 const T beta,  
                 DenseMatrix<T>& C,
                 const unsigned int num_threads)
{
    const unsigned int height_a = A.Height();

    // use sequential code to handle tiny matrices
    if (height_a < num_threads)
        return GemmSequentialAB(alpha, A, B, beta, C);

    // Scale C by beta as a dense operation.  Timing experiments show that
    // this is generally faster for a two-column matrix than having each 
    // thread traverse the columns individually, as is done in
    // gemm_rank2_thread_func.
    Scal(beta, C);

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmRank2AB<T>, 
                                      k,
                                      num_threads,
                                      alpha, std::cref(A), std::cref(B),
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
void GemmRank2ABt(const T alpha, 
                  const SparseMatrix<T>& A, 
                  const DenseMatrix<T>& B, 
                  const T beta,  
                  DenseMatrix<T>& C,
                  const unsigned int num_threads)
{
    const unsigned int height_a = A.Height();

    // use sequential code to handle tiny matrices
    if (height_a < num_threads)
        return GemmSequentialABt(alpha, A, B, beta, C);

    // Scale C by beta as a dense operation.  Timing experiments show that
    // this is generally faster for a two-column matrix than having each 
    // thread traverse the columns individually, as is done in
    // gemm_rank2_thread_func.
    Scal(beta, C);

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmRank2ABt<T>, 
                                      k,
                                      num_threads,
                                      alpha, std::cref(A), std::cref(B),
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
void GemmSequentialAB(const T alpha, 
                      const SparseMatrix<T>& A, 
                      const DenseMatrix<T>& B, 
                      const T beta, 
                      DenseMatrix<T>& C)
{
    // compute C = alpha*A*B + beta*C, A is sparse, B and C are dense

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (width_a != height_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (height_c != height_a)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (width_b != width_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");

    const unsigned int* cols_a = A.LockedColBuffer();
    const unsigned int* rows_a = A.LockedRowBuffer();
    const T*            data_a = A.LockedDataBuffer();

    const T*            data_b = B.LockedBuffer();
    T*                  data_c = C.Buffer();

    const unsigned int ldim_b = B.LDim();
    const unsigned int ldim_c = C.LDim();

    for (unsigned int j=0; j<width_b; ++j)
    {
        // offset to jth column of B and C
        unsigned int offset_b = j*ldim_b;
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

        // C(:,j) += alpha*A*B(:,j)

        // loop over columns of A
        for (unsigned int c=0; c != width_a; ++c)
        {
            unsigned int start = cols_a[c];
            unsigned int end   = cols_a[c+1];

            // alpha*B(c, j)
            T alpha_b_cj = alpha*data_b[offset_b + c];

            for (unsigned int offset = start; offset != end; ++offset)
            {
                // the next nonzero row in column j of A
                unsigned int row = rows_a[offset];
                data_c[offset_c + row] += alpha_b_cj*data_a[offset];
            }
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmSequentialABt(const T alpha, 
                       const SparseMatrix<T>& A, 
                       const DenseMatrix<T>& B, 
                       const T beta, 
                       DenseMatrix<T>& C)
{
    // compute C = alpha*A*B' + beta*C, A is sparse, B and C are dense

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (width_a != width_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (height_c != height_a)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (height_b != width_c)
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

        // C(:,j) += alpha*A*B'(:,j), or
        // C(:,j) += alpha*A*B(j, :)

        // loop over columns of A
        for (unsigned int c=0; c != width_a; ++c)
        {
            unsigned int start = cols_a[c];
            unsigned int end   = cols_a[c+1];

            // alpha*B'(c, j) == alpha*B(j, c)
            T alpha_b_cj = alpha*data_b[c*ldim_b + j];

            for (unsigned int offset = start; offset != end; ++offset)
            {
                // the next nonzero row in column j of A
                unsigned int row = rows_a[offset];
                data_c[offset_c + row] += alpha_b_cj*data_a[offset];
            }
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void GemmAB(const T alpha, 
            const SparseMatrix<T>& A, 
            const DenseMatrix<T>& B, 
            const T beta,  
            DenseMatrix<T>& C,
            const unsigned int max_threads)
{
    // C = alpha*A*B + beta*C

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (width_a != height_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (height_c != height_a)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (width_b != width_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");
    
    if (1 == max_threads)
        return GemmSequentialAB(alpha, A, B, beta, C);

    // handle rank 2 as a special case
    if (2 == width_b)
        return GemmRank2AB(alpha, A, B, beta, C, 2);

    // do not use more threads than columns of B
    unsigned int num_threads = max_threads;
    if (width_b < num_threads)
        num_threads = width_b;

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmAB<T>,
                                      k, num_threads,
                                      alpha, std::cref(A), std::cref(B),
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
void GemmABt(const T alpha, 
             const SparseMatrix<T>& A, 
             const DenseMatrix<T>& B, 
             const T beta,  
             DenseMatrix<T>& C,
             const unsigned int max_threads)
{
    // C = alpha*A*B' + beta*C

    const unsigned int height_a = A.Height();
    const unsigned int width_a  = A.Width();
    const unsigned int height_b = B.Height();
    const unsigned int width_b  = B.Width();
    const unsigned int height_c = C.Height();
    const unsigned int width_c  = C.Width();

    if (width_a != width_b)
        throw std::logic_error("Gemm: matrices A and B are non-conformant.");
    if (height_c != height_a)
        throw std::logic_error("Gemm: matrices A and C are non-conformant.");
    if (height_b != width_c)
        throw std::logic_error("Gemm: matrices B and C are non-conformant.");
    
    if (1 == max_threads)
        return GemmSequentialABt(alpha, A, B, beta, C);

    // handle rank 2 as a special case
    if (2 == height_b)
        return GemmRank2ABt(alpha, A, B, beta, C, 2);

    // do not use more threads than heightof B
    unsigned int num_threads = max_threads;
    if (height_b < num_threads)
        num_threads = height_b;

    std::vector<std::thread> threads;
    for (unsigned int k=0; k<num_threads; ++k)
    {
        threads.push_back(std::thread(ThreadFuncGemmABt<T>,
                                      k, num_threads,
                                      alpha, std::cref(A), std::cref(B),
                                      beta,
                                      C.Height(),
                                      C.Buffer(),
                                      C.LDim()));
    }

    // wait for threads to finish
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}

