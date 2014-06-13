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

#include "enums.hpp"
#include "thread_utils.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "sparse_gemm_ab_impl.hpp"
#include "sparse_gemm_ba_impl.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void Gemm(const Orientation orientationA,
          const Orientation orientationB,
          const T alpha, 
          const SparseMatrix<T>& A, 
          const DenseMatrix<T>& B, 
          const T beta,  
          DenseMatrix<T>& C)
{
    // Compute C = alpha*A*op(B) + beta*C, A is sparse, B and C are dense.
    
    // Matrix A can only have NORMAL orientation.
    // Matrix B can have either NORMAL or TRANSPOSE orientation.

    if (NORMAL != orientationA)
        throw std::logic_error("Gemm: TRANSPOSE orientation for sparse matrices is unsupported");

    unsigned int max_threads = GetMaxThreadCount();

    if (NORMAL == orientationB)
        GemmAB(alpha, A, B, beta, C, max_threads);
    else
        GemmABt(alpha, A, B, beta, C, max_threads);
}

//-----------------------------------------------------------------------------
template <typename T>
void Gemm(const Orientation orientationB,
          const Orientation orientationA,
          const T alpha, 
          const DenseMatrix<T>& B, 
          const SparseMatrix<T>& A, 
          const T beta,  
          DenseMatrix<T>& C)
{
    // Compute C = alpha*op(B)*A + beta*C, A is sparse, B and C are dense.

    // Matrix A can only have NORMAL orientation.
    // Matrix B can have either NORMAL or TRANSPOSE orientation.

    if (NORMAL != orientationA)
        throw std::logic_error("Gemm: TRANSPOSE orientation for sparse matrices is unsupported");

    unsigned int max_threads = GetMaxThreadCount();

    if (NORMAL == orientationB)
        GemmBA(alpha, B, A, beta, C, max_threads);
    else
        GemmBtA(alpha, B, A, beta, C, max_threads);
}

