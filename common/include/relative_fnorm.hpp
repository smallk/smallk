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
#include "elemental.hpp"
#include "sparse_matrix.hpp"

//-----------------------------------------------------------------------------
template <typename T>
class RelativeFnormW
{
public:
    RelativeFnormW() {}

    T Init(const elem::Matrix<T>& W)
    {
        Wprev_.ResizeTo(W.Height(), W.Width());
        elem::MakeZeros(Wprev_);
        return this->operator()(W);
    }

    T operator () (const elem::Matrix<T>& W)
    {        
        // compute diff of previous W and current W (Wprev = -W + Wprev)
        elem::Axpy(-1.0, W, Wprev_);
        double norm_diff   = elem::Norm(Wprev_, elem::FROBENIUS_NORM);
        double norm_cur    = elem::Norm(W,      elem::FROBENIUS_NORM);
        double fnorm_ratio = norm_diff / norm_cur;
        
        // W now becomes Wprev
        Wprev_ = W;   
        return fnorm_ratio;
    }

private:

    elem::Matrix<T> Wprev_;
};

