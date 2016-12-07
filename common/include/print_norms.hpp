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

#include <iostream>
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void PrintNormsWH(const DenseMatrix<T>& W,
                  const DenseMatrix<T>& H)
{
    // compute norms for W
    T max_norm = Norm(W, MAX_NORM);
    T one_norm = Norm(W, ONE_NORM);
    T inf_norm = Norm(W, INFINITY_NORM);
    T f_norm   = Norm(W, FROBENIUS_NORM);

    std::cout << "Norms for matrix W: " << std::endl;
    std::cout << "|| W ||_max = " << max_norm << "\n"
              << "|| W ||_1   = " << one_norm << "\n"
              << "|| W ||_oo  = " << inf_norm << "\n"
              << "|| W ||_F   = " << f_norm << "\n" << std::endl;

    // compute norms for H
    max_norm = Norm(H, MAX_NORM);
    one_norm = Norm(H, ONE_NORM);
    inf_norm = Norm(H, INFINITY_NORM);
    f_norm   = Norm(H, FROBENIUS_NORM);

    std::cout << "Norms for matrix H: " << std::endl;
    std::cout << "|| H ||_max = " << max_norm << "\n"
              << "|| H ||_1   = " << one_norm << "\n"
              << "|| H ||_oo  = " << inf_norm << "\n"
              << "|| H ||_F   = " << f_norm << "\n" << std::endl;    
}
//-----------------------------------------------------------------------------
template <typename T>
void PrintNorms(const DenseMatrix<T>& A, 
                const DenseMatrix<T>& W,
                const DenseMatrix<T>& H)                
{
    PrintNormsWH(W, H);

    // compute the difference matrix between A and W*H
    DenseMatrix<T> Diff(A);
    Gemm(NORMAL,
         NORMAL,
         -1.0, W, H, 1.0, Diff);

    T max_norm = Norm(Diff, MAX_NORM);
    T one_norm = Norm(Diff, ONE_NORM);
    T inf_norm = Norm(Diff, INFINITY_NORM);
    T   f_norm = Norm(Diff, FROBENIUS_NORM);

    std::cout << "Norms for the residual: " << std::endl;
    std::cout << "||A - W H||_max = " << max_norm << "\n"
              << "||A - W H||_1   = " << one_norm << "\n"
              << "||A - W H||_oo  = " << inf_norm << "\n"
              << "||A - W H||_F   = " << f_norm << "\n" << std::endl;
}

//-----------------------------------------------------------------------------
template <typename T>
void PrintNorms(const SparseMatrix<T>& A, 
                const DenseMatrix<T>& W,
                const DenseMatrix<T>& H)                
{
    PrintNormsWH(W, H);

    // Computation of the residual could be a time-consuming or impossible 
    // computation if the sparse matrix is large enough, so skip this.

    // // compute the product W*H
    // DenseMatrix<T> C(W.Height(), H.Width());
    // elem::Gemm(elem::NORMAL, elem::NORMAL, 1.0, W, H, 0.0, C);

    // // subtract sparse matrix A from WH
    // Add(C, -1.0, A);

    // T max_norm = elem::Norm(C, elem::MAX_NORM);
    // T one_norm = elem::Norm(C, elem::ONE_NORM);
    // T inf_norm = elem::Norm(C, elem::INFINITY_NORM);
    // T   f_norm = elem::Norm(C, elem::FROBENIUS_NORM);

    // std::cout << "Norms for the residual: " << std::endl;
    // std::cout << "||A - W H||_max = " << max_norm << "\n"
    //           << "||A - W H||_1   = " << one_norm << "\n"
    //           << "||A - W H||_oo  = " << inf_norm << "\n"
    //           << "||A - W H||_F   = " << f_norm << "\n" << std::endl;
}

