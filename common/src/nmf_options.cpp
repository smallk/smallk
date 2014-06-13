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

#include "nmf.hpp"
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
bool IsValid(const NmfOptions& opts, bool validate_matrix)
{
    if (opts.k <= 0)
    {
        cerr << "nmflib error: k-value must be a positive integer" << endl;
        return false;
    }

    if (validate_matrix)
    {
        // m, n, and k are required to be nonzero

        // m, n, k must all be > 0
        if (opts.height <= 0)
        {
            cerr << "nmflib error: matrix height must be a positive integer" << endl;
            return false;
        }
        if (opts.width <= 0)
        {
            cerr << "nmflib error: matrix width must be a positive integer" << endl;
            return false;
        }
        
        // k <= n
        if (opts.k > opts.width)
        {
            cerr << "nmflib error: k value cannot exceed the number of columns" << endl;
            return false;
        }
    }

    // 0 < tolerance < 1.0
    if ( (opts.tol <= 0.0) || (opts.tol >= 1.0))
    {
        cerr << "nmflib error: tolerance must be in the interval (0.0, 1.0)" << endl;
        return false;
    }

    // min iterations > 0
    if (opts.min_iter <= 0)
    {
        cerr << "nmflib error: miniter must be a positive integer" << endl;
        return false;
    }

    // max iterations > 0
    if (opts.max_iter <= 0)
    {
        cerr << "nmflib error: maxiter must be a positive integer" << endl;
        return false;
    }

    // tolcount > 0
    if (opts.tolcount <= 0)
    {
        cerr << "nmflib error: tolcount must be a positive integer" << endl;
        return false;
    }

    // an algorithm must be specified
    if ( (NmfAlgorithm::MU    != opts.algorithm) &&
         (NmfAlgorithm::HALS  != opts.algorithm) &&
         (NmfAlgorithm::RANK2 != opts.algorithm) &&
         (NmfAlgorithm::BPP   != opts.algorithm))
    {
        cerr << "nmflib error: unknown NMF algorithm specified" << endl;
        return false;
    }

    // if using RANK2 algorithm, must have k==2
    if (NmfAlgorithm::RANK2 == opts.algorithm)
    {
        if (2 != opts.k)
        {
            cerr << "nmflib error: RANK2 algorithm requires k == 2" << endl;
            return false;
        }
    }

    if ( //(NmfProgressAlgorithm::PG != opts.prog_est_algorithm) &&
         (NmfProgressAlgorithm::PG_RATIO != opts.prog_est_algorithm) &&
         (NmfProgressAlgorithm::DELTA_FNORM != opts.prog_est_algorithm))
    {
        cerr << "nmflib error: unknown stopping criterion specified" << endl;
        return false;
    }

    return true;
}
