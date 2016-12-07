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
#include <limits>
#include <vector>
#include <sstream>
#include <iostream>
#include "nmf.hpp"
#include "tree.hpp"
#include "timer.hpp"
#include "clust.hpp"
#include "utils.hpp"
#include "terms.hpp"
#include "random.hpp"
#include "size_check.hpp"
#include "assignments.hpp"
#include "file_loader.hpp"
#include "command_line.hpp"
#include "sparse_matrix.hpp"
#include "delimited_file.hpp"
#include "matrix_generator.hpp"
#include "flat_clust_output.hpp"
#include "hierclust_writer_factory.hpp"

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
    Timer timer;

    Random rng;
    rng.SeedFromTime();
    //rng.SeedFromInt(78);

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
    unsigned int m, n, nnz, ldim_a, num_clusters;

    //-------------------------------------------------------------------------
    //
    //                 load the dictionary file
    //
    //-------------------------------------------------------------------------
    
    if (opts.clust_opts.verbose)
        cout << "loading dictionary..." << endl;

    std::vector<std::string> dictionary;
    if (!LoadStringsFromFile(opts.dictfile, dictionary))
    {
        cerr << "\ncould not load dictionary file " << opts.dictfile << endl;
        NmfFinalize();
        return -1;
    }

    //-------------------------------------------------------------------------
    //
    //                 load matrix A, the data matrix
    //
    //-------------------------------------------------------------------------
    if (opts.clust_opts.verbose)
        cout << "loading matrix..." << endl;

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

    num_clusters = opts.clust_opts.num_clusters;

    // HierNMF2 requires W and H initializer matrices of dimensions
    //
    //    W: m x 2
    //    H: 2 x n
    //
    // Elemental uses a default signed 32-bit integer for its index type,
    // and it computes offsets into the data buffer with this same type.
    // The following lines check to see that the W and H matrices fit
    // within these limits.
    
    // check W matrix (2*m elements required)
    uint64_t required_size = static_cast<uint64_t>(m);
    required_size *= 2u;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "W matrix size too large" << endl;
        NmfFinalize();
        return -1;
    }

    // check H matrix (2*n elements required)
    required_size = static_cast<uint64_t>(n);
    required_size *= 2u;
    if (!FitsWithin<int>(required_size))
    {
        cerr << "H matrix size too large" << endl;
        NmfFinalize();
        return -1;
    }

    opts.clust_opts.nmf_opts.height = m;
    opts.clust_opts.nmf_opts.width  = n;
    opts.clust_opts.nmf_opts.k      = 2;

    // leading dimensions for dense matrix A data buffer
    ldim_a = m;
    
    // print a summary of all options
    if (opts.clust_opts.verbose)
        PrintOpts(opts);

    //-------------------------------------------------------------------------
    //
    //                 run the selected clustering algorithm
    //
    //-------------------------------------------------------------------------

    // W and H buffer for flat clustering
    std::vector<R> buf_w(m*num_clusters);
    std::vector<R> buf_h(num_clusters*n);

    Tree<R> tree;
    ClustStats stats;
    std::vector<float> probabilities;
    std::vector<unsigned int> assignments_flat;
    std::vector<int> term_indices(opts.clust_opts.maxterms * num_clusters);
    Result result = Result::OK;

    timer.Start();

    if (A.Size() > 0)
    {
        result = ClustSparse(opts.clust_opts, A,
                             &buf_w[0], &buf_h[0], tree, stats, rng);
    }
    else
    {
        result = Clust(opts.clust_opts, &buf_a[0], ldim_a,
                       &buf_w[0], &buf_h[0], tree, stats, rng);
    }
    
    if (opts.clust_opts.flat)
    {
        // compute flat clustering assignments and top terms
        unsigned int k = num_clusters;
        ComputeFuzzyAssignments(probabilities, &buf_h[0], k, k, n);
        ComputeAssignments(assignments_flat, &buf_h[0], k, k, n);
        TopTerms(opts.clust_opts.maxterms, &buf_w[0], m, m, k, term_indices);        
    }
 
    timer.Stop();
    double elapsed = timer.ReportMilliseconds();

    cout << "\nElapsed wall clock time: ";
    if (elapsed < 1000.0)
        cout << elapsed << " ms." << endl;
    else
        cout << elapsed*0.001 << " s." << endl;

    int num_converged = stats.nmf_count - stats.max_count;
    cout << num_converged << "/" << stats.nmf_count << " factorizations"
         << " converged." << endl << endl;

    //-------------------------------------------------------------------------
    //
    //                 write results
    //
    //-------------------------------------------------------------------------

    if (Result::FAILURE == result)
    {
        cerr << "\nHierarchical clustering fatal error." << endl;
    }
    else
    {    
        if (opts.clust_opts.verbose)
            cout << "Writing output files..." << endl;

        if (!tree.WriteAssignments(opts.assignfile))
            cerr << "\terror writing assignments file" << endl;

        IHierclustWriter* writer = CreateHierclustWriter(opts.format);
        if (!tree.WriteTree(writer, opts.treefile, dictionary))
            cerr << "\terror writing factorization file" << endl;

        if (opts.clust_opts.flat && Result::FLATCLUST_FAILURE != result)
        {
            FlatClustWriteResults(opts.outdir, assignments_flat, probabilities,
                                  dictionary, term_indices, opts.format,
                                  opts.clust_opts.maxterms, n,
                                  opts.clust_opts.num_clusters);
        }

        delete writer;
    }

    NmfFinalize();
    return 0;
}
