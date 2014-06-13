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

#include <stdexcept>
#include "nmf.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "projected_gradient.hpp"

template <typename T,
          template <typename> class MatrixType>
class ProgEstGeneric;

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType>
class ProgEstGenericDeltaW : public ProgEstGeneric<T, MatrixType>
{
    // This class estimates progress by observing the relative
    // change in the Frobenius norm of W on successive iterations.
public:

    T Init(const MatrixType<T>& A,
           const DenseMatrix<T>& W,
           const DenseMatrix<T>& H)
    {
        Wprev_.Resize(W.Height(), W.Width());
        MakeZeros(Wprev_);
        return this->Compute(W);
    }

    T Update(const unsigned int iter,
             const DenseMatrix<T>& W, 
             const DenseMatrix<T>& H, 
             const DenseMatrix<T>& gradW,
             const DenseMatrix<T>& gradH)
    {
        return Compute(W);
    }

private:
    
    DenseMatrix<T> Wprev_;

    T Compute(const DenseMatrix<T>& W)
    {
        // compute diff of previous W and current W (Wprev = -W + Wprev)
        Axpy( T(-1.0), W, Wprev_);
        double norm_diff   = Norm(Wprev_, FROBENIUS_NORM);
        double norm_cur    = Norm(W,      FROBENIUS_NORM);
        double fnorm_ratio = norm_diff / norm_cur;

        // W now becomes Wprev
        Wprev_ = W;   
        return fnorm_ratio;
    }
};

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType>
class ProgEstGenericPgRatio : public ProgEstGeneric<T, MatrixType>
{
    // Method of projected gradient ratios.
public:

    T Init(const MatrixType<T>& A,
           const DenseMatrix<T>& W,
           const DenseMatrix<T>& H)
    {
        return T(1.0);
    }

    T Update(const unsigned int iter,
             const DenseMatrix<T>& W, 
             const DenseMatrix<T>& H, 
             const DenseMatrix<T>& gradW,
             const DenseMatrix<T>& gradH)
    {
        if (0 == iter)
        {
            // save value after first iteration
            pg0_ = ProjectedGradientNorm<T>(gradW, gradH, W, H);
            return T(1.0);
        }
        else
        {
            T proj_gradient = ProjectedGradientNorm<T>(gradW, gradH, W, H);
            return proj_gradient / pg0_;
        }
    }

private:

    T pg0_;
};

//-----------------------------------------------------------------------------
template <typename T,
          template <typename> class MatrixType>
class ProgEstGeneric
{
public:
    virtual ~ProgEstGeneric() {}

    // factory method
    static ProgEstGeneric<T, MatrixType>* Create(const NmfAlgorithm algorithm,
                                                 const NmfProgressAlgorithm prog_method)
    {
        ProgEstGeneric<T, MatrixType>* result = nullptr;

        if (NmfProgressAlgorithm::DELTA_FNORM == prog_method)
        {
            result = new ProgEstGenericDeltaW<T, MatrixType>;
        }
        else if (NmfProgressAlgorithm::PG_RATIO == prog_method)
        {
            result = new ProgEstGenericPgRatio<T, MatrixType>;
        }
        else
            throw std::logic_error("Unknown progress estimation algorithm.");
        
        return result;
    }

    virtual T Init(const MatrixType<T>& A,
                   const DenseMatrix<T>& W,
                   const DenseMatrix<T>& H) {return T(1.0);}

    // compute an updated metric value
    virtual T Update(const unsigned int iter,
                     const DenseMatrix<T>& W, 
                     const DenseMatrix<T>& H, 
                     const DenseMatrix<T>& gradW,
                     const DenseMatrix<T>& gradH) {return T(1.0);}
};
