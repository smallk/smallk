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
#include "clust.hpp"
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
bool IsValid(const ClustOptions& opts, bool validate_matrix)
{
    if (validate_matrix)
    {
        // m, n, k must all be > 0
        if (opts.nmf_opts.height <= 0)
        {
            cerr << "error: matrix height must be a positive integer" << endl;
            return false;
        }
        if (opts.nmf_opts.width <= 0)
        {
            cerr << "error: matrix width must be a positive integer" << endl;
            return false;
        }
        if (opts.nmf_opts.k <= 0)
        {
            cerr << "error: cluster count must be a positive integer" << endl;
            return false;
        }
        
        // k <= n
        if (opts.nmf_opts.k > opts.nmf_opts.width)
        {
            cerr << "error: k value cannot exceed the matrix width" << endl;
            return false;
        }
    }

    // number of clusters must be >= 2 (num_clusters == 1 isn't meaningful...)
    if (opts.num_clusters <= 1)
    {
        cerr << "error: number of clusters must be >= 2" << endl;
        return false;
    }

    // 0 < tolerance < 1.0
    if ( (opts.nmf_opts.tol <= 0.0) || (opts.nmf_opts.tol >= 1.0))
    {
        cerr << "error: tolerance must be in the interval (0.0, 1.0)" << endl;
        return false;
    }

    // min iterations > 0
    if (opts.nmf_opts.min_iter <= 0)
    {
        cerr << "error: miniter must be a positive integer" << endl;
        return false;
    }

    // max iterations > 0
    if (opts.nmf_opts.max_iter <= 0)
    {
        cerr << "error: maxiter must be a positive integer" << endl;
        return false;
    }

    // maxterms > 0
    if (opts.maxterms <= 0)
    {
        cerr << "error: maxterms must be a positive integer" << endl;
        return false;
    }

    // trial_allowance for hierarchical clustering should not be negative
    if (opts.trial_allowance < 0)
    {
        cerr << "error: trial_allowance for hierarchical clustering is negative" << endl;
        return false;
    }

    // the unbalanced parameter for hierarchical clustering should be in [0, 1) interval
    if ( (opts.unbalanced < 0.0) || (opts.unbalanced >= 1.0))
    {
        cerr << "error: the unbalanced value should be in the interval [0, 1)" << endl;
        return false;
    }

    if ( (NmfProgressAlgorithm::PG_RATIO != opts.nmf_opts.prog_est_algorithm) &&
         (NmfProgressAlgorithm::DELTA_FNORM != opts.nmf_opts.prog_est_algorithm))
    {
        cerr << "error: unknown stopping criterion " << endl;
        return false;
    }

    return true;
}
