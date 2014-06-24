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

#include <string>
#include <vector>
#include <iostream>
#include "nmf.hpp"
#include "timer.hpp"
#include "utils.hpp"
#include "random.hpp"
#include "constants.hpp"
#include "size_check.hpp"
#include "file_loader.hpp"
#include "command_line.hpp"
#include "sparse_matrix.hpp"
#include "delimited_file.hpp"
#include "matrix_generator.hpp"

using std::cout;
using std::cerr;
using std::endl;

// datatype of all matrices
typedef double R;

static const R RNG_CENTER = R(0.5);
static const R RNG_RADIUS = R(0.5);

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    Random rng;
    rng.SeedFromTime();

    CommandLineOptions opts;
    if (!ParseCommandLine(argc, argv, opts))
    {
        if (opts.show_help)
        {
            ShowHelp(argv[0]);
            return 0;
        }
        else
        {
            // command line error
            return -1;
        }
    }

    // validate command line options
    if (!IsValid(opts))
        return -1;

    NmfInitialize(argc, argv);

    bool ok = true;
    unsigned int m, n, k, nnz, ldim_a, ldim_w, ldim_h;

    //-------------------------------------------------------------------------
    //
    //                 load matrix A, the matrix to be factored
    //
    //-------------------------------------------------------------------------
    if (opts.nmf_opts.verbose)
        cout << "Loading matrix..." << endl;

    SparseMatrix<R> A;
    std::vector<R> buf_a;

    if (IsSparse(opts.infile_A))
    {
        if (!LoadSparseMatrix(opts.infile_A, A, m, n, nnz))
        {
            cerr << "\nload failed for file " << opts.infile_A << endl;
            NmfFinalize();
            return -1;
        }
    }
    else if (IsDense(opts.infile_A))
    {
        ok = LoadDenseMatrix(opts.infile_A, buf_a, m, n);
        if (!ok || (buf_a.size() < m*n))
        {
            cerr << "\nload failed for file " << opts.infile_A << endl;
            NmfFinalize();
            return -1;
        }
    }
    else
    {
        cerr << "\nunsupported file type: " << opts.infile_A << endl;
        NmfFinalize();
        return -1;
    }

    opts.nmf_opts.height = m;
    opts.nmf_opts.width = n;
    k = opts.nmf_opts.k;
    
    // leading dimensions for data buffers
    ldim_a = m;
    ldim_w = m;
    ldim_h = k;

    // Elemental's offset computations are performed with the same datatype
    // as the index type, which for this code is a signed 32-bit integer.
    // Check the required sizes of the W and H matrices and see if they fit.
    
    // check W matrix size
    uint64_t required_size = static_cast<uint64_t>(m);
    required_size *= k;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "W matrix size too large" << endl;
        NmfFinalize();
        return -1;
    }

    // check H matrix size
    required_size = static_cast<uint64_t>(n);
    required_size *= k;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "H matrix size too large" << endl;
        NmfFinalize();
        return -1;
    }

    //-------------------------------------------------------------------------
    //
    //                 initialize factor matrices W and H
    //
    //-------------------------------------------------------------------------

    unsigned int height_w = m, width_w = k;
    unsigned int height_h = k, width_h = n;

    // allocate space for W and H matrices
    std::vector<R> buf_w(m*k); // W is m x k
    std::vector<R> buf_h(k*n); // H is k x n

    if (opts.nmf_opts.verbose)
        cout << "Initializing matrix W..." << endl;

    // load W, the leftmost factor, always dense
    if (opts.infile_W.empty())
        ok = RandomMatrix(buf_w, m, k, rng, RNG_CENTER, RNG_RADIUS);
    else
        ok = LoadDelimitedFile(buf_w, height_w, width_w, opts.infile_W);
    if (!ok)
    {
        cerr << "\nload failed for file " << opts.infile_W << endl;
        NmfFinalize();
        return -1;
    }

    // check dimensions of W if loaded from file
    if ( (height_w != m) || (width_w != k))
    {
        cerr << "\tdimensions of matrix W are " << height_w
             << " x " << width_w << endl;
        cerr << "\texpected " << m << " x " << k << endl;
        NmfFinalize();
        return -1;
    }

    if (opts.nmf_opts.verbose)
        cout << "Initializing matrix H..." << endl;

    // load H, the rightmost factor, always dense
    if (opts.infile_H.empty())
        ok = RandomMatrix(buf_h, k, n, rng, RNG_CENTER, RNG_RADIUS);
    else
        ok = LoadDelimitedFile(buf_h, height_h, width_h, opts.infile_H);
    if (!ok)
    {
        cerr << "\nload failed for file " << opts.infile_H << endl;
        NmfFinalize();
        return -1;
    }

    // check dimensions of H if loaded from file
    if ( (height_h != k) || (width_h != n))
    {
        cerr << "\tdimensions of matrix H are " << height_h
             << " x " << width_h << endl;
        cerr << "\texpected " << k << " x " << n << endl;
        NmfFinalize();
        return -1;
    }

    // now that all matrices are initialized, print a summary of all options
    if (opts.nmf_opts.verbose)
        PrintOpts(opts);

    //-------------------------------------------------------------------------
    //
    //                 run the selected NMF algorithm
    //
    //-------------------------------------------------------------------------

    NmfStats stats;
    Result result = Result::OK;

    if (A.Size() > 0)
    {
        result = NmfSparse(opts.nmf_opts, 
                           A.Height(), A.Width(), A.Size(),
                           A.LockedColBuffer(), 
                           A.LockedRowBuffer(), 
                           A.LockedDataBuffer(),
                           &buf_w[0], ldim_w, 
                           &buf_h[0], ldim_h,
                           stats);
    }
    else
    {
        result = Nmf(opts.nmf_opts, 
                     &buf_a[0], ldim_a, 
                     &buf_w[0], ldim_w, 
                     &buf_h[0], ldim_h, 
                     stats);
    }

    std::cout << "Elapsed wall clock time: ";
    std::cout << ElapsedTime(stats.elapsed_us) << std::endl;
    std::cout << std::endl;

    if (Result::OK != result)
    {
        cerr << "\nNMF solver failure." << endl;
    }
    else
    {
        if (opts.nmf_opts.verbose)
            cout << "Writing output files..." << endl;

        WriteDelimitedFile(&buf_w[0], ldim_w, m, k, opts.outfile_W, opts.output_precision);
        WriteDelimitedFile(&buf_h[0], ldim_h, k, n, opts.outfile_H, opts.output_precision);
    }

    NmfFinalize();
    return 0;
}
