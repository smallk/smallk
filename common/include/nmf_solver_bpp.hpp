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

#include <cassert>
#include <stdexcept>
#include "nnls.hpp"
#include "normal_eq.hpp"
#include "bit_matrix.hpp"
#include "bit_matrix_ops.hpp"
#include "dense_matrix.hpp"
#include "vector_utils.hpp"
#include "openmp_pragma.hpp"

//-----------------------------------------------------------------------------
template <typename T>
bool BppSolveNormalEq(const std::vector<unsigned int>& col_indices, 
                      const BitMatrix& passive_set, // full passive set
                      const DenseMatrix<T>& WtW,    // kxk, full cross-product matrix
                      const DenseMatrix<T>& WtAsub, // kxn, WtA(:, col_indices)
                      DenseMatrix<T>& Xsub)         // kxn, X(:, col_indices)
{
    // Solve W'W * Xsub = WtAsub for Xsub.  The 'col_indices' array is the set 
    // of non_opt_cols from the NnlsBlockpivot function.  The 'passive_set' 
    // matrix is the full passive set.  Matrix WtW is the full cross-product 
    // matrix, assumed to be symmetric positive definite. 
    // Matrix WtAsub is equal to WtA(:, col_indices).
    // Matrix Xsub is equal to X(:, col_indices).

    const unsigned int k = WtW.Height();
    const unsigned int n = WtAsub.Width();

    assert(col_indices.size() == n);

    if (WtW.Width() != WtW.Height())
        throw std::logic_error("BppSolveNormalEq: LHS matrix is not square");
    if (WtAsub.Height() != WtW.Height())
        throw std::logic_error("BppSolveNormalEq: LHS and RHS matrices are non-conformant");
    if ( (Xsub.Height() != WtW.Height()) || (Xsub.Width() != WtAsub.Width()))
        throw std::logic_error("BppSolveNormalEq: non-conformant solution matrix");

    if (IsEmpty(passive_set) || AllCols(passive_set, col_indices))
    {
        // solve WtW * Xsub = WtAsub using full matrices
        return SolveNormalEq(WtW, WtAsub, Xsub);
    }

    MakeZeros(Xsub);
    unsigned int num_rows = 0;
    std::vector<unsigned int> indices(n), groups, ri(k), ci(n);
    
    DenseMatrix<T> Lsub(k, k);
    DenseMatrix<T> Rsub(k, n);
    
    if (1 == n)
    {
        // The passive set has a single column, so there is a single
        // 'grouped' column (column 0), with a group size of 1.
        ci[0] = 0;
        
        // The call to RowIndices requires column indices from the
        // full passive set, which here is col_indices[0].
        passive_set.RowIndices(col_indices[0], ri, num_rows);
        
        // Solve the normal equation subproblem; will throw if non-HPD.
        WtW.Submatrix(Lsub, ri, ri, num_rows, num_rows);
        WtAsub.Submatrix(Rsub, ri, ci, num_rows, 1);
        
        if (!SolveNormalEq(Lsub, Rsub))
            return false;

        Overwrite(Xsub, Rsub, ri, ci, num_rows, 1);
    }
    else
    {
        // Find all columns having identical bit patterns in the passive
        // set matrix.  The normal equation solver (in the loop below) 
        // will combine the identical columns into a single RHS matrix 
        // and solve all of these columns en masse.  The 'indices' array
        // contains the column indices sorted by hash value, so groups of
        // identical columns appear together.  The 'groups' array contains
        // the number of identical columns in each group.  Only passive 
        // set columns in the 'col_indices' array are considered.
        
        // IMPORTANT: the indices array is used to group columns of 
        // matrix WtAsub, so these indices span the range 
        // [0, WtAsub.Width()).  These are NOT the indices in the original
        // matrix WtA, which likely exceed the width of WtAsub.
        passive_set.GroupIdenticalColumns(indices, groups, col_indices);
        
        // Solve the normal equations for each set of grouped columns.
        unsigned int count=0;
        for (unsigned int g=0; g<groups.size(); ++g)
        {
            // number of columns in this group
            unsigned int group_size = groups[g];
            assert(group_size > 0);
            
            // extract column indices for this group
            for (unsigned int c=count; c < (count + group_size); ++c)
                ci[c-count] = indices[c];
            
            // Row indices are extracted from the full passive set,
            // so the original indices are required. Observe that the
            // value q = indices[count] is a submatrix index, spanning the
            // range [0, col_indices.size()-1].  Thus the original index
            // is col_indices[q].
            passive_set.RowIndices(col_indices[indices[count]], ri, num_rows);
            if (num_rows > 0)
            {
                // Solve the normal equation subproblem; will throw if non-HPD.
                WtW.Submatrix(Lsub, ri, ri, num_rows, num_rows);
                WtAsub.Submatrix(Rsub, ri, ci, num_rows, group_size);

                if (!SolveNormalEq(Lsub, Rsub))
                    return false;

                Overwrite(Xsub, Rsub, ri, ci, num_rows, group_size);
            }
            
            // have processed this many columns so far
            count += group_size;
        }
        
        // this condition is true if all cols were processed
        assert(count == n);
    }

    return true;
}
                      
//-----------------------------------------------------------------------------
template <typename T>
bool BppSolveNormalEqNoGroup(const std::vector<unsigned int>& col_indices,
                             const BitMatrix& passive_set, // full passive set
                             const DenseMatrix<T>& WtW,    // kxk, full cross-product matrix
                             const DenseMatrix<T>& WtAsub, // kxn, WtA(:, col_indices)
                             DenseMatrix<T>& Xsub)         // kxn, X(:, col_indices)
{
    // Solve W'W * Xsub = WtAsub for Xsub.  The 'col_indices' array is the set
    // of non_opt_cols from the NnlsBlockpivot function.  The 'passive_set'
    // matrix is the full passive set.  Matrix WtW is the full cross-product
    // matrix, assumed to be symmetric positive definite.
    // Matrix WtAsub is equal to WtA(:, col_indices).
    // Matrix Xsub is equal to X(:, col_indices).
    // This function does not use column grouping, but instead solves all
    // problems in parallel.

    const unsigned int k = WtW.Height();
    const unsigned int n = WtAsub.Width();

    assert(col_indices.size() == n);

    if (WtW.Width() != WtW.Height())
        throw std::logic_error("BppSolveNormalEqNoGroup: LHS matrix is not square");
    if (WtAsub.Height() != WtW.Height())
        throw std::logic_error("BppSolveNormalEqNoGroup: LHS and RHS matrices are non-conformant");
    if ( (Xsub.Height() != WtW.Height()) || (Xsub.Width() != WtAsub.Width()))
        throw std::logic_error("BppSolveNormalEqNoGroup: non-conformant solution matrix");
    
    if (IsEmpty(passive_set) || AllCols(passive_set, col_indices))
    {
        // solve WtW * Xsub = WtAsub using full matrices
        return SolveNormalEq(WtW, WtAsub, Xsub);
    }

    MakeZeros(Xsub);    
    bool success = true;

    OPENMP_PRAGMA(omp parallel for)
    for (unsigned int c=0; c<n; ++c)
    {
        unsigned int num_rows;
        std::vector<unsigned int> ri(k), ci(1);
        DenseMatrix<T> Lsub(k, k), Rsub(k, 1), x_view, rhs_view;
        
        // col index in the original passive set
        unsigned int col_index = col_indices[c];
        
        // create views of the cth column of Xsub and WtAsub
        View      (x_view,     Xsub, 0, c, k, 1);
        LockedView(rhs_view, WtAsub, 0, c, k, 1);
        
        // row indices are from the 'col_index' column of the passive set
        passive_set.RowIndices(col_index, ri, num_rows);
        if (num_rows > 0)
        {
            // the cth column of WtAsub is the 0th column of the view
            ci[0] = 0;
            
            // solve the normal equation subproblem
            WtW.Submatrix(Lsub, ri, ri, num_rows, num_rows);
            rhs_view.Submatrix(Rsub, ri, ci, num_rows, 1);
            
            if (!SolveNormalEq(Lsub, Rsub))
            {
                // note the failure and keep going until all threads finish
                success = false;
                continue;
            }
            
            Overwrite(x_view, Rsub, ri, ci, num_rows, 1);
        }
    }

    return success;
}

//-----------------------------------------------------------------------------
template <typename T>
bool BppSolveNormalEqLeftNoGroup(const std::vector<unsigned int>& row_indices,
                                 const BitMatrix& passive_set, // full passive set
                                 const DenseMatrix<T>& HHt,    // kxk, full cross-product matrix
                                 const DenseMatrix<T>& AHtsub, // mxk, AHt(row_indices, :)
                                 DenseMatrix<T>& Xsub)         // mxk, X(row_indices, :)
{
    // Solve Xsub * HHt = AHtsub for Xsub.  The 'col_indices' array is the set
    // of non_opt_cols from the NnlsBlockpivot function.  The 'passive_set'
    // matrix is the full passive set.  Matrix HHt is the full cross-product
    // matrix, assumed to be symmetric positive definite.
    // Matrix AHtsub is equal to AHt(row_indices, :).
    // Matrix Xsub is equal to X(row_indices, :).
    // This function does not use column grouping, but instead solves all
    // problems in parallel.

    const unsigned int k = HHt.Height();
    const unsigned int m = AHtsub.Height();

    assert(row_indices.size() == m);

    if (HHt.Width() != HHt.Height())
        throw std::logic_error("BppSolveNormalEqLeftNoGroup: LHS matrix is not square");
    if (AHtsub.Width() != HHt.Width())
        throw std::logic_error("BppSolveNormalEqLeftNoGroup: non-conformant RHS matrix");
    if ( (Xsub.Height() != AHtsub.Height()) || (Xsub.Width() != HHt.Height()))
        throw std::logic_error("BppSolveNormalEqLeftNoGroup: non-conformant solution matrix");
    
    if (IsEmpty(passive_set) || AllRows(passive_set, row_indices))
    {
        // solve Xsub * HHt = AHtsub using full matrices
        return SolveNormalEqLeft(HHt, AHtsub, Xsub);
    }

    MakeZeros(Xsub);
    bool success = true;
    
    OPENMP_PRAGMA(omp parallel for)
    for (unsigned int r=0; r<m; ++r)
    {
        unsigned int num_cols;
        std::vector<unsigned int> ri(1), ci(k);
        DenseMatrix<T> Lsub(k, k), Rsub(1, k), x_tmp(1, k), x_view, rhs_view;
        
        // row index in the original passive set
        unsigned int row_index = row_indices[r];
        
        // create views of the rth row of Xsub and AHtsub
        View      (x_view,     Xsub, r, 0, 1, k);
        LockedView(rhs_view, AHtsub, r, 0, 1, k);
        
        // col indices are from the 'row_index' row of the passive set
        passive_set.ColIndices(row_index, ci, num_cols);
        if (num_cols > 0)
        {
            x_tmp.Resize(1, num_cols);
            
            // the rth row of AHtsub is te 0th row of the view
            ri[0] = 0;
            HHt.Submatrix(Lsub, ci, ci, num_cols, num_cols);
            rhs_view.Submatrix(Rsub, ri, ci, 1, num_cols);

            if (!SolveNormalEqLeft(Lsub, Rsub, x_tmp))
            {
                // note the failure and keep going until all threads finish
                success = false;
                continue;
            }

            Overwrite(x_view, x_tmp, ri, ci, 1, num_cols);
        }
    }

    return success;
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class Matrix >
class Solver_Generic_BPP
{
public:

    // ------------------------------------------------------------------------
    //
    //         Init function - call prior to start of iterations
    //
    // ------------------------------------------------------------------------
    void Init(const Matrix<T>& A,
              DenseMatrix<T>& W,
              DenseMatrix<T>& H)
    {
        unsigned int m = A.Height();
        unsigned int n = A.Width();
        unsigned int k = W.Width();

        // This solver requires the transpose of A, so compute once and save.
        Transpose(A, At);

        // set sizes of auxiliary matrices
        Wt.Resize(k, m);
        gradWt.Resize(k, m);
        WtA.Resize(k, n);
        HAt.Resize(k, m);
        WtW.Resize(k, k);
        HHt.Resize(k, k);

        // compute WtW and WtA, to setup for iterations
        Gemm(TRANSPOSE, NORMAL, T(1), W, W, T(0), WtW);
        Gemm(TRANSPOSE, NORMAL, T(1), W, A, T(0), WtA);

        // load Wt with the initial values of W
        Transpose(W, Wt);
    }

    // ------------------------------------------------------------------------
    //
    //         Update function - call on each iteration
    //
    // ------------------------------------------------------------------------
    bool operator() (const Matrix<T>& A,
                     DenseMatrix<T>& W,
                     DenseMatrix<T>& H,
                     DenseMatrix<T>& gradW,
                     DenseMatrix<T>& gradH)
    {
        // solve (WtW)*H = (WtA) for H and gradH
        if (!NnlsBlockpivot(WtW, WtA, H, gradH))
            return false;
        
        // solve (HHt)*W' = HA' for W' and gradW'

        // compute HHt and HAt
        Gemm(NORMAL, TRANSPOSE, T(1), H, H,  T(0), HHt);
        Gemm(NORMAL, NORMAL,    T(1), H, At, T(0), HAt);

        // solve for W' and gradW'
        if (!NnlsBlockpivot(HHt, HAt, Wt, gradWt))
            return false;

        // transpose back
        Transpose(Wt, W);
        Transpose(gradWt, gradW);

        // compute gradH = (W'W)*H - W'A  using the updated W

        // WtW and WtA
        Gemm(TRANSPOSE, NORMAL, T(1), W, W, T(0), WtW);
        Gemm(TRANSPOSE, NORMAL, T(1), W, A, T(0), WtA);

        // gradH
        Gemm(NORMAL, NORMAL, 1.0, WtW, H, 0.0, gradH);
        Axpy(-1.0, WtA, gradH);

        return true;
    }

private:

    Matrix<T> At;
    DenseMatrix<T> Wt, gradWt, WtA, HAt, WtW, HHt;
};
