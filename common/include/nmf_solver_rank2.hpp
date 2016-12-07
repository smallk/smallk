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

#include <cmath>
#include <limits>
#include <stdexcept>
#include "dense_matrix.hpp"
#include "normalize.hpp"

//-----------------------------------------------------------------------------
template <typename T>
bool SystemSolveH(DenseMatrix<T>& X, // H    2 x n
                  DenseMatrix<T>& A, // WtW  2 x 2
                  DenseMatrix<T>& B) // WtA  2 x n
{
    // This function solves the system A*X = B for X.  The columns of X and B
    // are treated as separate subproblems.  Each subproblem can be expressed
    // as follows, for column i of X and B (the ' means transpose):
    //
    //              A * [x0  x1]' = [b0  b1]'
    //
    // The solution proceeds by transforming matrix A to upper triangular form
    // via a fast Givens rotation.  The resulting triangular system for the 
    // two unknowns in each column of X is then solved by backsubstitution.
    //
    // The fast Givens method relies on the fact that matrix A can be 
    // expressed as the product of a diagonal matrix D and another matrix Y,
    // and that this form is preserved under the rotation.  If J is the 
    // Givens rotation matrix, then:
    //
    //           AX    = B
    //          JAX    = JB        J is the Givens rotation matrix
    //         (JA)X   = JB        letting A = D1Y1
    //        (JD1Y1)X = JB
    //        (D2Y2)X  = JB
    //
    // The matrix D1 is the identity matrix, and the matrix Y1 is the
    // original matrix A.  The matrix J has the form [c -s; c s] in Matlab
    // notation.  NOTE: with this convention, the tangent has a minus sign:
    //
    //            t = -A(1, 0) / A(0, 0)
    //
    // Matrix A has the form [a b; c d], and after the rotation it is
    // transformed to [a2 b2; 0 d2].
    //
    // There are two forms for the D2 and Y2 matrices, depending on whether
    // cosines or sines are 'factored out'.  If |A00| >= |A01|, then the
    // upper left element is factored out and the 'cosine' formulation is
    // used.  If |A00| < |A01|, then the upper right element is factored
    // out and the 'sine' formulation is used.
    //
    // A general reference for this code is the paper 'Fast Plane Rotations
    // with Dynamic Scaling', by A. Anda and H. Park, SIAM Journal on Matrix
    // Analysis and Applications, Vol 15, no. 1, pp. 162-174, Jan. 1994.

    int n = B.Width();
    T abs_A00 = std::abs(A.Get(0, 0));
    T abs_A01 = std::abs(A.Get(0, 1));
    const T epsilon = std::numeric_limits<T>::epsilon();

    if ( (abs_A00 < epsilon) && (abs_A01 < epsilon))
    {
        std::cerr << "SystemSolveH: singular matrix" << std::endl;
        return false;
    }

    T a2, b2, d2, e2, f2, inv_a2, inv_d2, x_1;
    if (abs_A00 >= abs_A01)
    {   
        // use 'cosine' formulation; t is the tangent
        T t = -A.Get(1,0) / A.Get(0,0);
        a2 = A.Get(0,0) - t*A.Get(1,0);
        b2 = A.Get(0,1) - t*A.Get(1,1);
        d2 = A.Get(1,1) + t*A.Get(0,1);

        // precompute 1/a2 and 1/d2 to avoid repeated division
        inv_a2 = T(1.0) / a2;
        inv_d2 = T(1.0) / d2;

        // a2 is guaranteed to be positive
        if (std::abs(d2/a2) < epsilon)
            return false;

        // solve the upper triangular systems by backsubstitution
        for (int i=0; i<n; ++i)
        {
            e2 = B.Get(0,i) - t*B.Get(1,i);
            f2 = B.Get(1,i) + t*B.Get(0,i);
            x_1 = f2 * inv_d2;
            X.Set(1, i, x_1);
            X.Set(0, i, (e2 - b2*x_1)*inv_a2);
        }
    }
    else
    {
        // use 'sine' formulation; ct is the cotangent
        T ct = -A.Get(0,0) / A.Get(1,0);
        a2 = -A.Get(1,0) + ct*A.Get(0,0);
        b2 = -A.Get(1,1) + ct*A.Get(0,1);
        d2 =  A.Get(0,1) + ct*A.Get(1,1);

        // precompute 1/a2 and 1/d2 to avoid repeated division
        inv_a2 = T(1.0) / a2;
        inv_d2 = T(1.0) / d2;

        // a2 is guaranteed to be positive
        if (std::abs(d2/a2) < epsilon)
            return false;

        // solve the upper triangular systems by backsubstitution
        for (int i=0; i<n; ++i)
        {
            e2 = -B.Get(1,i) + ct*B.Get(0,i);
            f2 =  B.Get(0,i) + ct*B.Get(1,i);
            x_1 = f2 * inv_d2;
            X.Set(1, i, x_1);
            X.Set(0, i, (e2 - b2*x_1) * inv_a2);
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
bool SystemSolveW(DenseMatrix<T>& X,  // m x 2
                  DenseMatrix<T>& A,  // 2 x 2
                  DenseMatrix<T>& B)  // m x 2
{
    // Solve XA = B; use the code from the previous solver, but transpose
    // everything and change t to -t.

    const T eps = std::numeric_limits<double>::epsilon();
    int m = B.Height();

    if (std::abs(A.Get(0,0)) < eps && std::abs(A.Get(0,1)) < eps) 
    {
        std::cerr << "SystemSolveW: singular matrix" << std::endl;
        return false;
    }

    T a2, b2, d2, e2, f2, x_1;
    T inv_a2, inv_d2;

    if (std::abs(A.Get(0,0)) >= std::abs(A.Get(0, 1))) 
    {
        // use 'cosine' formulation; t is the tangent
        T t = A.Get(0, 1) / A.Get(0,0);
        a2 = A.Get(0,0) + t*A.Get(0,1);
        b2 = A.Get(1,0) + t*A.Get(1,1);
        d2 = A.Get(1,1) - t*A.Get(1,0);

        // precompute 1/a2 and 1/d2 to avoid repeated division
        inv_a2 = T(1.0) / a2;
        inv_d2 = T(1.0) / d2;

        // a2 is guaranteed to be positive
        if (std::abs(d2/a2) < eps)
            return false;

        // solve the upper triangular systems by backsubstitution
        for (int i=0; i<m; ++i)
        {
            e2 = B.Get(i,0) + t*B.Get(i,1);
            f2 = B.Get(i,1) - t*B.Get(i,0);
            x_1 = f2 * inv_d2;
            X.Set(i, 1, x_1);
            X.Set(i, 0, (e2 - b2*x_1)*inv_a2);
        }
    } 
    else 
    {
        // use 'sine' formulation; ct is the cotangent
        T ct = A.Get(0,0) / A.Get(0,1);
        a2 = -A.Get(0,1) - ct*A.Get(0,0);
        b2 = -A.Get(1,1) - ct*A.Get(1,0);
        d2 =  A.Get(1,0) - ct*A.Get(1,1);

        // precompute 1/a2 and 1/d2 to avoid repeated division
        inv_a2 = T(1.0) / a2;
        inv_d2 = T(1.0) / d2;

        // a2 is guaranteed to be positive
        if (std::abs(d2/a2) < eps)
            return false;

        // solve the upper triangular systems by backsubstitution
        for (int i=0; i<m; ++i)
        {
            e2 = -B.Get(i,1) - ct*B.Get(i,0);
            f2 =  B.Get(i,0) - ct*B.Get(i,1);
            x_1 = f2 * inv_d2;
            X.Set(i, 1,  x_1);
            X.Set(i, 0, (e2 - b2*x_1) * inv_a2);
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
template <typename T>
void OptimalActiveSetH(DenseMatrix<T>& H,         // 2 x n
                       const DenseMatrix<T>& WtW, // 2 x 2
                       const DenseMatrix<T>& WtA) // 2 x n
{
    // Remove negative entries in each of the n columns of matrix H, but
    // do so in a manner that minimizes the overall objective function. 
    // Each column can be considered in isolation.
    //
    // The problem for each column i of H, h = H(:,i),  is as follows:
    //
    //     min_{h>=0} |Wh - a|^2
    //
    // This is minimized when h* = w'a / (w'w).  Expressed in terms of the
    // individual elements, this becomes:
    //
    //    h*[0] = w[0]'a / (w[0]'w[0])
    //    h*[1] = w[1]'a / (w[1]'w[1])

    // If both elements of h* are nonnegative, this is the optimal solution.  
    // If any elements are negative, the Rank2 algorithm is used to adjust
    // their values.
    
    int width = H.Width();
    
    const T wtw_00 = WtW.Get(0,0);
    const T wtw_11 = WtW.Get(1,1);
    const T inv_wtw_00 = T(1.0) / wtw_00;
    const T inv_wtw_11 = T(1.0) / wtw_11;
    const T sqrt_wtw_00 = sqrt(wtw_00);
    const T sqrt_wtw_11 = sqrt(wtw_11);

    for (int i=0; i<width; ++i)
    {
        T v1 = WtA.Get(0,i) * inv_wtw_00;
        T v2 = WtA.Get(1,i) * inv_wtw_11;
        T vv1 = v1 * sqrt_wtw_00;
        T vv2 = v2 * sqrt_wtw_11;

        if (vv1 >= vv2)
            v2 = T(0);
        else
            v1 = T(0);

        if ( (H.Get(0,i) <= T(0)) || (H.Get(1,i) <= T(0)))
        {
            H.Set(0, i, v1);
            H.Set(1, i, v2);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T>
void OptimalActiveSetW(DenseMatrix<T>& W,         // m x 2
                       const DenseMatrix<T>& HHt, // 2 x 2
                       const DenseMatrix<T>& AHt) // m x 2
{
    // Remove negative entries in each of the m rows of matrix W, but
    // do so in a manner that minimizes the overall objective function. 
    // Each row can be considered in isolation.
    //
    // The problem for each row i of W, w = W(i,:), is as follows:
    //
    //     min_{w>=0} |wH - a|^2
    //
    // This is minimized when w'* = Ha'/(hh'), or w* = aH'/(hh').  
    // Expressed in terms of the individual elements, this becomes:
    //
    //    w*[0] = ah'[0] / (h[0]h'[0])
    //    w*[1] = ah'[1] / (h[1]h'[1])

    // If both elements of w* are nonnegative, this is the optimal solution.  
    // If any elements are negative, the Rank2 algorithm is used to adjust
    // their values.
    
    int height = W.Height();
    
    const T hht_00 = HHt.Get(0,0);
    const T hht_11 = HHt.Get(1,1);
    const T inv_hht_00 = T(1.0) / hht_00;
    const T inv_hht_11 = T(1.0) / hht_11;
    const T sqrt_hht_00 = sqrt(hht_00);
    const T sqrt_hht_11 = sqrt(hht_11);

    for (int i=0; i<height; ++i)
    {
        T v1 = AHt.Get(i, 0) * inv_hht_00;
        T v2 = AHt.Get(i, 1) * inv_hht_11;
        T vv1 = v1 * sqrt_hht_00;
        T vv2 = v2 * sqrt_hht_11;

        if (vv1 >= vv2)
            v2 = T(0);
        else
            v1 = T(0);

        if ( (W.Get(i, 0) <= T(0)) || (W.Get(i, 1) <= T(0)))
        {
            W.Set(i, 0, v1);
            W.Set(i, 1, v2);
        }
    }
}

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class Matrix>
class Solver_Generic_Rank2
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
        WtW.Resize(W.Width(), W.Width());
        WtA.Resize(W.Width(), A.Width());
        HHt.Resize(H.Height(), H.Height());
        AHt.Resize(A.Height(), H.Height());
        ScaleFactors.Resize(H.Height(), 1);

        // compute the kxk matrix WtW = W' * W
        Gemm(TRANSPOSE, NORMAL, T(1), W, W, T(0), WtW);

        // compute the kxn matrix WtA =  W' * A
        Gemm(TRANSPOSE, NORMAL, T(1), W, A, T(0), WtA);
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
        //
        // Solve this problem for H:  min_{H>=0} |WH - A|^2
        //
        //         W is mx2, H is 2xn, A is mxn.
        //
        // Step 1: Solve the unconstrained least squares problem for H.
        //
        //         Objective function phi(H) = (WH - A)'(WH-A)
        //
        // The objective function is minimized for H satisfying this equation:
        //  
        //         W'W * H = W'A
        // 
        if (!SystemSolveH(H, WtW, WtA))
            return false;

        // Step 2: Adjust the solution of the unconstrained problem by 
        //         searching for the optimal active set.  
        //
        // The solver may return an H matrix that contains negative elements.
        // These negative elements need to be removed.  Remove them in a way
        // that minimizes the overall objective function.
        //
        OptimalActiveSetH(H, WtW, WtA);

        // Do the same for W:

        // compute the kxk matrix HHt = H * H'
        Gemm(NORMAL, TRANSPOSE, T(1), H, H, T(0), HHt);
        
        // compute the mxk matrix AHt =  A * H'
        Gemm(NORMAL, TRANSPOSE, T(1), A, H, T(0), AHt);

        //
        // Solve this problem for W:  min_{W>=0} |H'W' - A'|^2
        //
        //         W is mx2, H is 2xn, A is mxn
        //
        // Step 1: Solve the unconstrained least squares problem for W.
        //
        //         Objective function phi(W) = (H'W' - A')'(H'W' - A')
        //
        // The objective function is minimized for W' satisfying this equation:
        //
        //         HH' * W' = H'A, or the equivalent after a transpose
        //         W * (HH') = AH'
        //      
        // Solve this system for the unknown matrix W.
        //
        if (!SystemSolveW(W, HHt, AHt))
            return false;

        // Step 2: Adjust the solution of the unconstrained problem by 
        //         searching for the optimal active set.  
        //
        // The solver may return a W matrix that contains negative elements.
        // These negative elements need to be removed.  Remove them in a way
        // that minimizes the overall objective function.
        //
        OptimalActiveSetW(W, HHt, AHt);

        NormalizeAndScale(W, H, ScaleFactors);

        // At this point both W and H are fully updated.  Now compute
        // gradW = W*HHt - AHt and gradH = WtW*H - WtA.
        
        // First compute HHt and AHt.  HHt is a 2x2 matrix and can be updated
        // quickly by using the scale factors returned by the call to 
        // NormalizeAndScale.  A call to Gemm is NOT necessary.

        T s0 = ScaleFactors.Get(0, 0);
        T s1 = ScaleFactors.Get(1, 0);
        
        // update HH' for the rescaled H
        T elt00 = HHt.Get(0,0);
        T elt01 = HHt.Get(0,1);
        T elt11 = HHt.Get(1,1);
        HHt.Set(0, 0, elt00*s0*s0);
        HHt.Set(0, 1, elt01*s0*s1);
        HHt.Set(1, 0, elt01*s0*s1);  // since HH' is symmetric
        HHt.Set(1, 1, elt11*s1*s1);
        
        // update the mxk matrix AHt =  A * H' by scaling each column
        DiagonalScale(RIGHT, NORMAL, ScaleFactors, AHt);
        
        // Now gradW can be computed, since HHt and AHt are both current.
        Gemm(NORMAL, NORMAL, T(1), W, HHt, T(0), gradW);
        Axpy( T(-1), AHt, gradW);

        // Compute gradH by first updating WtW and WtA.
        Gemm(TRANSPOSE, NORMAL, T(1), W, W, T(0), WtW);
        Gemm(TRANSPOSE, NORMAL, T(1), W, A, T(0), WtA);
        Gemm(NORMAL, NORMAL, T(1), WtW, H, T(0), gradH);
        Axpy( T(-1), WtA, gradH);

        return true;
    }

private:

    DenseMatrix<T> ScaleFactors;
    DenseMatrix<T> WtW, WtA, HHt, AHt;
};
