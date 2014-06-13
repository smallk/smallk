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

// generate matrices of various types

#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include <vector>
#include "utils.hpp"
#include "random.hpp"
#include "matrix_types.hpp"
#include "command_line.hpp"
#include "sparse_matrix_decl.hpp"
#include "sparse_matrix_impl.hpp"
#include "sparse_matrix_ops.hpp"
#include "sparse_matrix_io.hpp"
#include "delimited_file.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    Random rng;
    rng.SeedFromTime();

    // print usage info if no command line args were specified
    std::string prog_name(argv[0]);
    if (1 == argc)
    {
        print_usage(prog_name);
        return 0;
    }
    
    CommandLineOptions opts;
    parse_command_line(argc, argv, opts);
    if (!validate_opts(opts))
        return -1;

    unsigned int m = opts.height;
    unsigned int n = opts.width;
    double center = opts.rng_center;
    double radius = opts.rng_radius;

    SparseMatrix<double> S;
    std::vector<double> A(m*n);

    bool is_sparse = false;
    if (CaseInsensitiveEquals(TYPE_UNIFORM, opts.type))
    {
        for (unsigned int c=0; c != n; ++c)
        {
            unsigned int col_offset = c*m;
            for (unsigned int r=0; r != m; ++r)
                A[col_offset + r] = rng.RandomDouble(center, radius);
        }
    }
    else if (CaseInsensitiveEquals(TYPE_DENSEDIAG, opts.type))
    {
        for (unsigned int c=0; c != n; ++c)
        {
            unsigned int col_offset = c*m;
            for (unsigned int r=0; r != m; ++r)
                A[col_offset + r] = 0.0;
        }

        for (unsigned int c=0; c != n; ++c)
            A[c*m + c] = rng.RandomDouble(center, radius);
    }
    else if (CaseInsensitiveEquals(TYPE_SPARSEDIAG, opts.type))
    {
        is_sparse = true;
        S.Reserve(n, n, n);

        S.BeginLoad();
        for (unsigned int c=0; c != n; ++c)
            S.Load(c, c, rng.RandomDouble(center, radius));
        S.EndLoad();
    }
    else if (CaseInsensitiveEquals(TYPE_ONES, opts.type))
    {
        for (unsigned int c=0; c != n; ++c)
        {
            unsigned int col_offset = c*m;
            for (unsigned int r=0; r != m; ++r)
                A[col_offset + r] = 1.0;
        }
    }
    else if (CaseInsensitiveEquals(TYPE_ZEROS, opts.type))
    {
        for (unsigned int c=0; c != n; ++c)
        {
            unsigned int col_offset = c*m;
            for (unsigned int r=0; r != m; ++r)
                A[col_offset + r] = 0.0;
        }
    }
    else if (CaseInsensitiveEquals(TYPE_SPARSE, opts.type))
    {
        is_sparse = true;
        RandomSparseMatrix(rng, S, opts.nz_per_col, m, m, n, n);
    }


    if (is_sparse)
    {
        if (!WriteMatrixMarketFile(opts.filename, S, opts.precision))
            cerr << "matrixgen error - sparse matrix file write failed" << endl;
     }
    else
    {
        // write matrix to file in CSV format
        if (!WriteDelimitedFile(&A[0], m, m, n, opts.filename, opts.precision))
             cerr << "matrixgen error - file write failed" << endl;
     }

    return 0;
}
