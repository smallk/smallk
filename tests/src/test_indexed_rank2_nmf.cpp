#include <iostream>
#include <thread>
#include "nmf.hpp"
#include "elemental.hpp"
#include "sparse_matrix.hpp"
#include "matrix_generator.hpp"
#include "delimited_file.hpp"
#include "index_set.hpp"
#include "random.hpp"
#include "timer.hpp"
#include "file_loader.hpp"
#include "nmf_solve_sparse.hpp"
#include "nmf_solve_indexed.hpp"
#include "nmf_progress_estimation.hpp"

using std::cout;
using std::cerr;
using std::endl;

namespace IndexedRank2NmfTest
{
    const unsigned int NUM_RUNS = 16;
    const double MAX_ACCEPTABLE_FNORM = 1.0e-10;
    const unsigned int MAX_DIM = 1024;
}

//-----------------------------------------------------------------------------
bool TestIndexedRank2Nmf()
{    
    Random rng;
    Timer timer;
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (0 == max_threads)
        max_threads = 2;

    SparseMatrix<double> A, Asubset;
    elem::Matrix<double> H, W, Hinit, Winit, W1, W2, H1, H2;
    elem::Matrix<double> Hresult1, Hresult2, Wresult1, Wresult2;
    std::vector<unsigned int> old_to_new_rows, new_to_old_rows;
    std::vector<unsigned int> old_to_new_rows2, new_to_old_rows2;

    double min_sparsity = 1.0, max_sparsity = 0.0;
    double max_norm = 0.0, elapsed_subset = 0.0, elapsed_indexed = 0.0;

    cout << "Running indexed rank2 Nmf test..." << endl;

    ProgressEstimator<double> progress_est1(NmfProgressAlgorithm::PG);
    ProgressEstimator<double> progress_est2(NmfProgressAlgorithm::PG);
    Solver_Sparse_Rank2<double> solver;
    Solver_Sparse_Indexed_Rank2<double> solver_indexed;

    NmfStats nmf_stats;
    NmfOptions nmf_opts;
    nmf_opts.tol = 0.0001;
    nmf_opts.algorithm = NmfAlgorithm::RANK2;
    nmf_opts.prog_est_algorithm = NmfProgressAlgorithm::PG;
    nmf_opts.k = 2;
    nmf_opts.min_iter = 5;
    nmf_opts.max_iter = 5000;
    nmf_opts.tolcount = 1;
    nmf_opts.max_threads = max_threads;
    nmf_opts.verbose = false;
    nmf_opts.normalize = true;

    // unsigned int m, n, nnz;
    // if (!LoadSparseMatrix(std::string("/Users/richardboyd/repos/bitbucket/xdata2/data/reuters.mtx"),
    //                       A, m, n, nnz))
    // {
    //     cerr << "\nmatrix load failed" << endl;
    //     return false;
    // }

    for (unsigned int i=0; i<IndexedRank2NmfTest::NUM_RUNS; ++i)
    {
        unsigned int height_a = rng.RandomRangeInt(512, IndexedRank2NmfTest::MAX_DIM);
        unsigned int width_a  = rng.RandomRangeInt(512, IndexedRank2NmfTest::MAX_DIM);

        // sparse_fraction ranges from (0.1 - 0.09, 0.1 + 0.09)
        double sparse_fraction = rng.RandomDouble(0.1, 0.09);
        unsigned int nonzeros_per_col = sparse_fraction * height_a;
        if (0 == nonzeros_per_col)
            nonzeros_per_col = 1;        
        if (sparse_fraction < min_sparsity)
            min_sparsity = sparse_fraction;
        if (sparse_fraction > max_sparsity)
            max_sparsity = sparse_fraction;

        RandomSparseMatrix(A, nonzeros_per_col, rng, height_a, width_a);

        // create random H and W
        Hinit.ResizeTo(2, width_a);
        RandomMatrix(Hinit, rng);

        Winit.ResizeTo(height_a, 2);
        RandomMatrix(Winit, rng);

        // resize the result matrices
        Hresult1.ResizeTo(2, width_a);
        Hresult2.ResizeTo(2, width_a);
        Wresult1.ResizeTo(height_a, 2);
        Wresult2.ResizeTo(height_a, 2);
        
        // generate a random set of column indices of A
        size_t index_set_size = 1 + rng.RandomRangeInt(1, width_a/2);
        std::vector<unsigned int> col_indices(index_set_size);
        RandomIndexSet(rng, width_a, index_set_size, col_indices);

        // compute the height of the most compact submatrix of A
        std::vector<unsigned int> old_to_new_rows, new_to_old_rows;
        A.SubMatrixRowsCompact(col_indices, old_to_new_rows, new_to_old_rows);
        unsigned int new_height = new_to_old_rows.size();

        // resize the submatrices to matxh
        H1.ResizeTo(2, index_set_size);
        H2.ResizeTo(2, index_set_size);
        W1.ResizeTo(new_height, 2);
        W2.ResizeTo(new_height, 2);

        // load the subset initializer matrices
        for (unsigned int c=0; c != index_set_size; ++c)
        {
            H1.Set(0, c, Hinit.Get(0, col_indices[c]));
            H1.Set(1, c, Hinit.Get(1, col_indices[c]));
            H2.Set(0, c, Hinit.Get(0, col_indices[c]));
            H2.Set(1, c, Hinit.Get(1, col_indices[c]));
        }

        for (unsigned int r=0; r != new_height; ++r)
        {
            W1.Set(r, 0, Winit.Get(new_to_old_rows[r], 0));
            W1.Set(r, 1, Winit.Get(new_to_old_rows[r], 1));
            W2.Set(r, 0, Winit.Get(new_to_old_rows[r], 0));
            W2.Set(r, 1, Winit.Get(new_to_old_rows[r], 1));
        }

        timer.Start();

        // extract a compact subset of A
        A.SubMatrixColsCompact(Asubset, col_indices, old_to_new_rows, new_to_old_rows);
        assert(Asubset.Height() == new_height);
        
        // run Rank2 NMF on the compact subproblem
        nmf_opts.height = Asubset.Height();
        nmf_opts.width  = Asubset.Width();
        bool ok = NmfSolve(Asubset, W1, H1, solver, progress_est1, nmf_opts, nmf_stats);
        if (!ok)
        {
            cerr << "[" << (i+1) << "]: subset NMF solver failure" << endl;
            continue;
        }

        timer.Stop();
        elapsed_subset += timer.ReportMicroseconds();

        // update results
        elem::MakeZeros(Hresult1);
        for (int c=0; c != H1.Width(); ++c)
        {
            Hresult1.Set(0, col_indices[c], H1.Get(0, c));
            Hresult1.Set(1, col_indices[c], H1.Get(1, c));
        }

        // update W
        elem::MakeZeros(Wresult1);
        for (int r=0; r != W1.Height(); ++r)
        {
            Wresult1.Set(new_to_old_rows[r], 0, W1.Get(r, 0));
            Wresult1.Set(new_to_old_rows[r], 1, W1.Get(r, 1));
        }

        // solve the indexed NMF problem

        timer.Start();

        nmf_opts.height = new_height;
        nmf_opts.width  = col_indices.size();
        ok = NmfSolveIndexed(A, W2, H2, solver_indexed, progress_est2, nmf_opts,
                             Asubset.Height(), old_to_new_rows, col_indices);
        if (!ok)
        {
            cerr << "[" << (i+1) << "]: indexed NMF solver failure" << endl;
            continue;
        }

        timer.Stop();
        elapsed_indexed += timer.ReportMicroseconds();

        // update results
        elem::MakeZeros(Hresult2);
        for (int c=0; c != H2.Width(); ++c)
        {
            Hresult2.Set(0, col_indices[c], H2.Get(0, c));
            Hresult2.Set(1, col_indices[c], H2.Get(1, c));
        }

        // update W
        elem::MakeZeros(Wresult2);
        for (int r=0; r != W2.Height(); ++r)
        {
            Wresult2.Set(new_to_old_rows[r], 0, W2.Get(r, 0));
            Wresult2.Set(new_to_old_rows[r], 1, W2.Get(r, 1));
        }

        // check residuals

        // H2 = -H1 + H2
        elem::Axpy(-1.0, Hresult1, Hresult2);
        double fnorm_h = elem::Norm(Hresult2, elem::FROBENIUS_NORM);

        // W2 = -W1 + W2
        elem::Axpy(-1.0, Wresult1, Wresult2);
        double fnorm_w = elem::Norm(Wresult2, elem::FROBENIUS_NORM);

        cout << "\t[" << (i+1) << "/" << IndexedRank2NmfTest::NUM_RUNS <<
            "] fnorm of W residual: " << fnorm_w << ", fnorm of H residual: "
             << fnorm_h << endl;

        max_norm = std::max(max_norm, fnorm_h);
        max_norm = std::max(max_norm, fnorm_w);
        if (max_norm > IndexedRank2NmfTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            break;
        }
    }

    // compute average runtimes per loop
    double t_copying = elapsed_subset * 0.001 / IndexedRank2NmfTest::NUM_RUNS;
    double t_indexed = elapsed_indexed * 0.001 / IndexedRank2NmfTest::NUM_RUNS;    

    cout << endl;
    cout << "\t**** Results for Indexed Rank2 NMF Test ****" << endl;
    cout << endl;
    cout << "\t\t" << IndexedRank2NmfTest::NUM_RUNS << " runs " << endl;
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;
    cout << "\t\tElapsed time with data copying:\t" << t_copying << " ms." << endl;
    cout << "\t\tElapsed time with indexing:\t" << t_indexed << " ms." << endl;
    cout << endl;
    cout << "\t********************************************" << endl;
    cout << endl;

    return (max_norm < IndexedRank2NmfTest::MAX_ACCEPTABLE_FNORM);
}
