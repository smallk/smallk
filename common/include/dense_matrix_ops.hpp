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
#include <iomanip>
#include <cstdlib>
#include "bit_matrix.hpp"
#include "openmp_pragma.hpp"

//-----------------------------------------------------------------------------
template <typename T>
void Cholesky(const UpperOrLower uplo, DenseMatrix<T>& A)
{
    EL::UpperOrLower elem_uplo = 
        (UPPER == uplo ? EL::UPPER : EL::LOWER);
    EL::Cholesky(elem_uplo, A.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void SolveAfterCholesky(const UpperOrLower uplo,
                        const Orientation orientation,
                        const DenseMatrix<T>& A,
                        DenseMatrix<T>& B)
{
    EL::UpperOrLower elem_uplo = 
        (UPPER == uplo ? EL::UPPER : EL::LOWER);
    EL::Orientation elem_orientation =
        (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);

    EL::cholesky::SolveAfter(elem_uplo, elem_orientation, A.M, B.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void SolveAfterLU(const Orientation orientation,
                  const DenseMatrix<T>& A,
                  const DenseMatrix<int>& P,
                  DenseMatrix<T>& B)
{
    EL::Orientation elem_orientation =
        (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);

    EL::lu::SolveAfter(elem_orientation, A.M, P.M, B.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void LU(DenseMatrix<T>& A, DenseMatrix<int>& P)
{EL::LU(A.M, P.M);}

//-----------------------------------------------------------------------------
template <typename T>
void GaussianElimination(DenseMatrix<T>& A, DenseMatrix<T>& B)
{EL::GaussianElimination(A.M, B.M);}

//-----------------------------------------------------------------------------
template <typename T>
void Pseudoinverse(DenseMatrix<T>& A)
{EL::Pseudoinverse(A.M);}

//-----------------------------------------------------------------------------
template <typename T>
void Trsm(const LeftOrRight side, 
          const UpperOrLower uplo,
          const Orientation orientation,
          const UnitOrNonUnit diag,
          const T alpha,
          const DenseMatrix<T>& A,
          DenseMatrix<T>& B)
{
    EL::LeftOrRight elem_side = 
        (LEFT == side ? EL::LEFT : EL::RIGHT);
    EL::UpperOrLower elem_uplo =
        (UPPER == uplo ? EL::UPPER : EL::LOWER);
    EL::Orientation elem_orient =
        (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);
    EL::UnitOrNonUnit elem_diag =
        (UNIT == diag ? EL::UNIT : EL::NON_UNIT);

    EL::Trsm(elem_side, elem_uplo, elem_orient, elem_diag, alpha, A.M, B.M);
}

//-----------------------------------------------------------------------------
template <typename T>
bool HasNaNs(const DenseMatrix<T>& M)
{
    // returns true if any matrix element is NaN

    int height = M.Height();
    int width  = M.Width();
    for (int c=0; c<width; ++c)
    {
        for (int r=0; r<height; ++r)
        {
            T elt = M.Get(r, c);
            if (std::isnan(elt))
                return true;
        }
    }

    return false;
}

//-----------------------------------------------------------------------------
template <typename T>
T Norm(const DenseMatrix<T>& A, 
       const NormType norm_type = NormType::FROBENIUS_NORM)
{
    switch(norm_type)
    {
    case MAX_NORM:
        return EL::Norm(A.M, EL::MAX_NORM);
        break;
    case ONE_NORM:
        return EL::Norm(A.M, EL::ONE_NORM);
        break;
    case INFINITY_NORM:
        return EL::Norm(A.M, EL::INFINITY_NORM);
        break;
    case FROBENIUS_NORM:
        return EL::Norm(A.M, EL::FROBENIUS_NORM);
        break;
    default:
        throw std::logic_error("Norm: unknown norm type");
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void Axpy(const T alpha, const DenseMatrix<T>& X, DenseMatrix<T>& Y)
{EL::Axpy(alpha, X.M, Y.M);}

//-----------------------------------------------------------------------------
template <typename T>
void View(DenseMatrix<T>& A, DenseMatrix<T>& B,
          const int i, const int j, const int height, const int width)
{EL::View(A.M, B.M, i, j, height, width);}

//-----------------------------------------------------------------------------
template <typename T>
void LockedView(DenseMatrix<T>& A, const DenseMatrix<T>& B,
                const int i, const int j, const int height, const int width)
{EL::LockedView(A.M, B.M, i, j, height, width);}

//-----------------------------------------------------------------------------
template <typename T>
T Nrm2(const DenseMatrix<T>& X)
{return EL::Nrm2(X.M);}

//-----------------------------------------------------------------------------
template <typename T>
void Scal(const T alpha, DenseMatrix<T>& X)
{
#if ELEM_VER >= 84
    EL::Scale(alpha, X.M);
#else
    EL::Scal(alpha, X.M);
#endif
}

//-----------------------------------------------------------------------------
template <typename T>
void MakeZeros(DenseMatrix<T>& X)
{
#if ELEM_VER >= 85
    EL::Zeros(X.M, X.M.Height(), X.M.Width());
#else
    EL::MakeZeros(X.M);
#endif
}

//-----------------------------------------------------------------------------
template <typename T>
void MakeUniform(DenseMatrix<T>& X)
{EL::MakeUniform(X.M);}

//-----------------------------------------------------------------------------
template <typename T>
T Dot(const DenseMatrix<T>& x, const DenseMatrix<T>& y)
{return EL::Dot(x.M, y.M);}

//-----------------------------------------------------------------------------
template <typename T>
void Copy(const DenseMatrix<T>& X, DenseMatrix<T>& Y)
{EL::Copy(X.M, Y.M);}

//-----------------------------------------------------------------------------
template <typename T>
void Transpose(const DenseMatrix<T>& X, DenseMatrix<T>& Y)
{EL::Transpose(X.M, Y.M);}

//-----------------------------------------------------------------------------
template <typename T>
void Print(const DenseMatrix<T>& X)
{
    //EL::Print(X.M);

    const unsigned int height = X.M.Height();
    const unsigned int width  = X.M.Width();

    std::cout << std::fixed;
    for (unsigned int r=0; r<height; ++r)
    {
        for (unsigned int c=0; c<width-1; ++c)
            std::cout << X.M.Get(r, c) << " ";
        std::cout << X.M.Get(r, width-1) << std::endl;
    }
    std::cout.unsetf(std::ios_base::floatfield);
}

//-----------------------------------------------------------------------------
template <typename T>
void HPDSolve(const UpperOrLower uplo, 
              const Orientation orientation,
              DenseMatrix<T>& A, 
              DenseMatrix<T>& B)
{
    EL::UpperOrLower elem_uplo = 
        (UPPER == uplo ? EL::UPPER : EL::LOWER);
    EL::Orientation elem_orientation =
        (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);

    EL::HPDSolve(elem_uplo, elem_orientation, A.M, B.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void LeastSquares(const Orientation orientation,
                  DenseMatrix<T>& A,
                  const DenseMatrix<T>& B,
                  DenseMatrix<T>& X)
{
    EL::Orientation elem_orientation =
        (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);

    EL::LeastSquares(elem_orientation, A.M, B.M, X.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void Gemm(const Orientation orientA,
          const Orientation orientB,
          const T alpha, 
          const DenseMatrix<T>& A,
          const DenseMatrix<T>& B,
          const T beta, 
          DenseMatrix<T>& C,
          const unsigned int max_threads = 0)
{
    EL::Orientation orientation_A = 
        (NORMAL == orientA ? EL::NORMAL : EL::TRANSPOSE);
    EL::Orientation orientation_B =
        (NORMAL == orientB ? EL::NORMAL : EL::TRANSPOSE);

    EL::Gemm(orientation_A, orientation_B, alpha, A.M, B.M, beta, C.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void DiagonalScale(const LeftOrRight side,
                   const Orientation orientation,
                   const DenseMatrix<T>& d,
                   DenseMatrix<T>& X)
{
    EL::LeftOrRight lr = (LEFT == side ? EL::LEFT : EL::RIGHT);
    EL::Orientation orient = (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);

    EL::DiagonalScale(lr, orient, d.M, X.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void Gemv(const Orientation orientation, 
          const T alpha, 
          const DenseMatrix<T>& A,
          const DenseMatrix<T>& x, 
          const T beta, 
          DenseMatrix<T>& y)
{
    EL::Orientation orient = (NORMAL == orientation ? EL::NORMAL : EL::TRANSPOSE);
    EL::Gemv(orient, alpha, A.M, x.M, beta, y.M);
}

//-----------------------------------------------------------------------------
template <typename T>
void OverwriteCols(DenseMatrix<T>& A, 
                   const DenseMatrix<T>& B, 
                   const std::vector<unsigned int>& col_indices,
                   const unsigned int num_cols)
{
    const unsigned int height = A.Height();

    // Overwrite columns of A with B.
    if (B.Height() != A.Height())
        throw std::logic_error("OverwriteCols: height mismatch");
    if (num_cols > static_cast<unsigned int>(A.Width()))
        throw std::logic_error("OverwriteCols: col indices out of bounds");

    T* buf_a = A.Buffer();
    const unsigned int ldim_a = A.LDim();

    const T* buf_b = B.LockedBuffer();
    const unsigned int ldim_b = B.LDim();

    for (unsigned int c=0; c<num_cols; ++c)
    {
        unsigned int col_a = col_indices[c];
        unsigned int offset_a = col_a*ldim_a;
        unsigned int offset_b = c*ldim_b;
            
        memcpy(&buf_a[offset_a], &buf_b[offset_b], height * sizeof(T));
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void Overwrite(DenseMatrix<T>& A, 
               const DenseMatrix<T>& B, 
               const std::vector<unsigned int>& row_indices,
               const std::vector<unsigned int>& col_indices,
               const unsigned int num_rows,
               const unsigned int num_cols)
{
    // Overwrite entries in A with entries in B.  The row and column
    // index arrays contain the destination indices to be overwritten.
    // Matrix B has size num_rows x num_cols.

    if (num_rows > static_cast<unsigned int>(A.Height()))
        throw std::logic_error("Overwrite: row indices out of bounds");
    if (num_cols > static_cast<unsigned int>(A.Width()))
        throw std::logic_error("Overwrite: col indices out of bounds");

    const T* buf_b = B.LockedBuffer();
    const unsigned int ldim_b = B.LDim();
    
    T* buf_a = A.Buffer();
    const unsigned int ldim_a = A.LDim();

    for (unsigned int c=0; c<num_cols; ++c)
    {
        unsigned int col_a = col_indices[c];
        unsigned int offset_a = ldim_a * col_a;
        unsigned int offset_b = ldim_b * c;

        for (unsigned int r=0; r<num_rows; ++r)
        {
            unsigned int row_a = row_indices[r];
            //T val = B.Get(r, c);
            //A.Set(row_a, col_a, val);
            buf_a[offset_a + row_a] = buf_b[offset_b + r];
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void ZeroizeSmallValues(DenseMatrix<T>& A, const T tol)
{
    // set Aij to zero if |Aij| < tol

    T* buf = A.Buffer();
    const unsigned int ldim = A.LDim();

    const unsigned int height = A.Height();
    const unsigned int width  = A.Width();

    OPENMP_PRAGMA(omp parallel for)
    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int col_offset = c*ldim;
        for (unsigned int r=0; r<height; ++r)
        {
            T val = buf[col_offset + r];
            if (std::abs(val) < tol)
            {
                buf[col_offset + r] = T(0);
            }
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void MakeDiagonallyDominant(DenseMatrix<T>& M)
{
    // Make the diagonal element larger than the row sum, to ensure that
    // the matrix is nonsingular.  All entries in the matrix are nonnegative, 
    // so no absolute values are needed.

    for (int r=0; r<M.Height(); ++r)
    {
        T row_sum = 0.0;
        for (int c=0; c<M.Width(); ++c)
            row_sum += M.Get(r, c);
        M.Set(r, r, row_sum + T(1));
    }
}
