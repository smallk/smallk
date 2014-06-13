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

// Various routines for solving (W'W)*X = W'A for H.  The matrix on the
// left is assumed to be symmetric positive-definite (SPD).

#include <vector>
#include <iostream>
#include <stdexcept>
#include "dense_matrix.hpp"

//-----------------------------------------------------------------------------
template <typename T>
bool SolveNormalEq(DenseMatrix<T>& LHS, // kxk
                   DenseMatrix<T>& X)   // kxn
{
    // Solve LHS * X = RHS by Cholesky factorization and two triangular
    // solves.  X is assumed to be initialized with the desired RHS.
    
    bool success = true;

    try
    {
        HPDSolve(UPPER, NORMAL, LHS, X);
    }
    catch (elem::NonHPSDMatrixException& e)
    {
        std::cerr << "Cholesky factorization failure - ";
        std::cerr << "matrix was not symmetric positive-definite." << std::endl;
        success = false;
    }
    catch (std::logic_error& e)
    {
        std::cerr << "Cholesky factorization failure - ";
        std::cerr << "matrix was not symmetric positive-definite." << std::endl;
        success = false;
    }

    return success;
    
}

//-----------------------------------------------------------------------------
template <typename T>
bool SolveNormalEq(const DenseMatrix<T>& LHS, // kxk
                   const DenseMatrix<T>& RHS, // kxn
                   DenseMatrix<T>& X)         // kxn
{
    // Solve LHS * X = RHS for X, where LHS is assumed SPD.

    if (LHS.Width() != LHS.Height())
        throw std::logic_error("SolveNormalEq: expected square matrix on LHS");
    if ( (X.Height() != LHS.Height()) || (X.Width() != RHS.Width()))
        throw std::logic_error("SolveNormalEq: non-conformant matrix X");
    
    // copy the input, since the solver overwrites it
    DenseMatrix<T> M(LHS);
    X = RHS;

    return SolveNormalEq(M, X);
}

//-----------------------------------------------------------------------------
template <typename T>
bool SolveNormalEqLeft(const DenseMatrix<T>& LHS, // kxk
                       const DenseMatrix<T>& RHS, // mxk
                       DenseMatrix<T>& X)         // mxk
{
    // Solve X * LHS = RHS for X, where LHS is symmetric positive-definite.

    if (LHS.Width() != LHS.Height())
        throw std::logic_error("SolveNormalEqLeft: expected square matrix on LHS");
    if ( (X.Width() != LHS.Height()) || (X.Height() != RHS.Height()))
        throw std::logic_error("SolveNormalEqLeft: non-conformant matrix X");
    
    // // copy the input, since the solver overwrites it
    DenseMatrix<T> U(LHS);

    // Compute the Cholesky factor LHS = U'U
    bool success = true;
    try
    {
        Cholesky(UPPER, U);
    }
    catch (elem::NonHPSDMatrixException& e)
    {
        std::cerr << "Cholesky factorization failure - ";
        std::cerr << "matrix was not symmetric positive-definite." << std::endl;
        success = false;
    }
    catch (std::logic_error& e)
    {
        std::cerr << "Cholesky factorization failure - ";
        std::cerr << "matrix was not symmetric positive-definite." << std::endl;
        success = false;
    }

    if (!success)
        return false;

    // Solve X (U'U) = RHS as follows:
    //
    // Group the two left matrices:
    // 
    //     (XU') U = RHS
    //
    // Let Y = XU', so the equation becomes:
    //
    //     YU = RHS
    //
    // Do a triangular solve: Y = (RHS) * inverse(U)

    X = RHS;
    Trsm(RIGHT, UPPER, NORMAL, NON_UNIT, T(1), U, X);

    // The solution Y is stored in X.  
    // 
    // Now solve XU' = Y, so X = Y * inverse(U')
    Trsm(RIGHT, UPPER, TRANSPOSE, NON_UNIT, T(1), U, X);

    // check
    // DenseMatrix<T> temp(m, k);
    // Gemm(NORMAL, NORMAL, T(1), X, LHS, T(0), temp);
    // Axpy(-1.0, RHS, temp);
    // double norm = Norm(temp, FROBENIUS_NORM);
    // std::cout << "\tSolveNormalEqLeft: norm = " << norm << std::endl;

    return true;
}
