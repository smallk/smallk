#include <iostream>
#include <thread>
#include "elemental.hpp"
#include "sparse_matrix.hpp"
#include "matrix_generator.hpp"
#include "delimited_file.hpp"
#include "index_set.hpp"
#include "random.hpp"
#include "timer.hpp"

using std::cout;
using std::cerr;
using std::endl;

// For the NMF formulation A~WH, we need the products AH' and W'A.

bool SimpleTest(Random& rng);
bool TestSparseGemmIndexedABt_Rank2 (Random& rng);
bool TestSparseGemmIndexedBtA_Rank2 (Random& rng);

namespace GemmIndexedTest
{
    const unsigned int NUM_RUNS = 64;
    const double MAX_ACCEPTABLE_FNORM = 1.0e-10;
    const unsigned int MAX_DIM = 8192;
}

//-----------------------------------------------------------------------------
bool TestSparseGemmIndexed()
{
    Random rng;
    rng.SeedFromTime();

    bool result0 = SimpleTest(rng);

    //bool result1 = TestSparseGemmIndexedAB(rng);
    bool result2 = TestSparseGemmIndexedABt_Rank2(rng);
    //bool result3 = TestSparseGemmIndexedBA(rng);
    bool result4 = TestSparseGemmIndexedBtA_Rank2(rng);

    return (result0 && result2 && result4);
}

//-----------------------------------------------------------------------------
bool SimpleTest(Random& rng)
{
    cout << "Running SimpleTest...";
    cout.flush();

    // A = 1 0 0 0  H = 1 2 3 4  == W'
    //     0 2 0 0      5 6 7 8
    //     0 0 3 0
    //     1 2 3 4

    SparseMatrix<double> A(4, 4, 7);
    A.BeginLoad();
    A.Load(0, 0, 1.0);
    A.Load(3, 0, 1.0);
    A.Load(1, 1, 2.0);
    A.Load(3, 1, 2.0);
    A.Load(2, 2, 3.0);
    A.Load(3, 2, 3.0);
    A.Load(3, 3, 4.0);
    A.EndLoad();

    elem::Matrix<double> H(2, 4);
    H.Set(0, 0, 1.0); H.Set(0, 1, 2.0); H.Set(0, 2, 3.0); H.Set(0, 3, 4.0);
    H.Set(1, 0, 5.0); H.Set(1, 1, 6.0); H.Set(1, 2, 7.0); H.Set(1, 3, 8.0);

    elem::Matrix<double> W(4, 2);
    W.Set(0, 0, 1.0); W.Set(0, 1, 5.0);
    W.Set(1, 0, 2.0); W.Set(1, 1, 6.0);
    W.Set(2, 0, 3.0); W.Set(2, 1, 7.0);
    W.Set(3, 0, 4.0); W.Set(3, 1, 8.0);

    elem::Matrix<double> C1(2, 2), C2(2, 2);

    std::vector<unsigned int> col_indices = {1, 3};
    std::vector<unsigned int> old_to_new_rows = {0xFFFFFFFF, 0, 0xFFFFFFFF, 1};
    std::vector<unsigned int> new_to_old_rows = {1, 3};
    unsigned int max_threads = 1;

    elem::Matrix<double> Hsubset(2, 2);
    Hsubset.Set(0, 0, H.Get(0, col_indices[0]));
    Hsubset.Set(0, 1, H.Get(0, col_indices[1]));
    Hsubset.Set(1, 0, H.Get(1, col_indices[0]));
    Hsubset.Set(1, 1, H.Get(1, col_indices[1]));

    //cout << "\nHsubset: " << endl;
    //elem::Print(Hsubset);

    elem::Matrix<double> Wsubset(2, 2);
    Wsubset.Set(0, 0, W.Get(new_to_old_rows[0], 0));
    Wsubset.Set(0, 1, W.Get(new_to_old_rows[0], 1));
    Wsubset.Set(1, 0, W.Get(new_to_old_rows[1], 0));
    Wsubset.Set(1, 1, W.Get(new_to_old_rows[1], 1));

    //cout << "\nWsubset: " << endl;
    //elem::Print(Wsubset);

    // compute the product C1 = AHsubset' =  4  12
    //                                      20  44
    GemmIndexed(NORMAL, TRANSPOSE,
                1.0, A, Hsubset, 0.0, C1, max_threads, 
                new_to_old_rows.size(), old_to_new_rows, col_indices);

    //cout << "\nC1: " << endl;
    //elem::Print(C1);

    bool ok1 = (4.0 == C1.Get(0, 0))  && (12.0 == C1.Get(0, 1)) &&
               (20.0 == C1.Get(1, 0)) && (44.0 == C1.Get(1, 1));

    // compute the product C2 = Wsubset'A =  5 12
    //                                      13 24
    GemmIndexed(TRANSPOSE, NORMAL,
                1.0, Wsubset, A, 0.0, C2, max_threads, 
                new_to_old_rows.size(), old_to_new_rows, col_indices);

    //cout << "\nC2: " << endl;
    //elem::Print(C2);

    bool ok2 = (12.0 == C2.Get(0, 0)) && (16.0 == C2.Get(0, 1)) &&
               (28.0 == C2.Get(1, 0)) && (32.0 == C2.Get(1, 1));

    cout << ( (ok1 && ok2) ? "passed." : "failed.") << endl;
    return (ok1 && ok2);
}

//-----------------------------------------------------------------------------
bool TestSparseGemmIndexedABt_Rank2(Random& rng)
{
    // C = alpha*A*B' + beta*C
    // 
    // Matrix A is sparse and is mxn, but use s cols from col index set S.
    // Matrix B has 2 rows; matrix C has 2 cols.

    Timer timer;
    unsigned int max_threads = 2;

    SparseMatrix<double> A, Asubset;
    elem::Matrix<double> B, Bsubset, C1, C2, C1subset, C2subset;

    double min_sparsity = 1.0, max_sparsity = 0.0;
    double max_norm = 0.0, elapsed_subset = 0.0, elapsed_indexed = 0.0;

    cout << "Running indexed GEMM rank2 AB' test..." << endl;
    
    for (unsigned int i=0; i<GemmIndexedTest::NUM_RUNS; ++i)
    {
        unsigned int height_a = rng.RandomRangeInt(512, GemmIndexedTest::MAX_DIM);
        unsigned int width_a  = rng.RandomRangeInt(512, GemmIndexedTest::MAX_DIM);

        // sparse_fraction ranges from (0.1 - 0.09, 0.1 + 0.09)
        double sparse_fraction = rng.RandomDouble(0.1, 0.09);
        unsigned int nonzeros_per_col = sparse_fraction * height_a;
        if (0 == nonzeros_per_col)
            nonzeros_per_col = 1;        
        if (sparse_fraction < min_sparsity)
            min_sparsity = sparse_fraction;
        if (sparse_fraction > max_sparsity)
            max_sparsity = sparse_fraction;

        // generate a random alpha and beta on [-1, 1]
        double alpha = rng.RandomDouble(0.0, 1.0);
        double beta  = rng.RandomDouble(0.0, 1.0);

        RandomSparseMatrix(A, nonzeros_per_col, rng, height_a, width_a);

        // set the dimensions of B and randomly initialize
        B.ResizeTo(2, width_a);
        RandomMatrix(B, rng);
        
        // set the dimensions of C1 and C2 and randomly initialize
        C1.ResizeTo(height_a, 2);
        RandomMatrix(C1, rng);
        C2 = C1;

        // generate a random set of column indices of A
        size_t index_set_size = rng.RandomRangeInt(1, width_a/2);
        std::vector<unsigned int> col_indices(index_set_size);
        RandomIndexSet(rng, width_a, index_set_size, col_indices);

        // compute the height of the most compact submatrix of A
        std::vector<unsigned int> old_to_new_rows, new_to_old_rows;
        A.SubMatrixRowsCompact(col_indices, old_to_new_rows, new_to_old_rows);
        unsigned int new_height = new_to_old_rows.size();
        
        // resize Bsubset and Csubset to match
        Bsubset.ResizeTo(2, index_set_size);
        C1subset.ResizeTo(new_height, 2);
        C2subset = C1subset;

        // extract Bsubset from B (cols of B from index set)
        for (unsigned int c=0; c<index_set_size; ++c)
        {
            for (int r=0; r<Bsubset.Height(); ++r)
            {
                double val = B.Get(r, col_indices[c]);
                Bsubset.Set(r, c, val);
            }
        }

        old_to_new_rows.clear();
        new_to_old_rows.clear();

        timer.Start();

        // extract a compact submatrix of A using the columns in col_indices
        A.SubMatrixColsCompact(Asubset, col_indices, old_to_new_rows, new_to_old_rows);
        assert(Asubset.Height() == new_height);

        // perform sparse Gemm on the subset
        Gemm(NORMAL, TRANSPOSE,
             alpha, Asubset, Bsubset, beta, C1subset, max_threads);

        timer.Stop();
        elapsed_subset += timer.ReportMicroseconds();

        timer.Start();

        // perform indexed Gemm on the original matrix
        GemmIndexed(NORMAL, TRANSPOSE,
                    alpha, A, Bsubset, beta, C2subset, max_threads, 
                    new_height, old_to_new_rows, col_indices);

        timer.Stop();
        elapsed_indexed += timer.ReportMicroseconds();

        // check residuals
        
        // C2 = -C1 + C2
        elem::Axpy(-1.0, C1subset, C2subset);
        double fnorm = elem::Norm(C2subset, elem::FROBENIUS_NORM);
        cout << "\t[" << (i+1) << "/" << GemmIndexedTest::NUM_RUNS <<
            "]: fnorm of residual: " << fnorm << endl;
        if (fnorm > max_norm)
            max_norm = fnorm;
        
        if (max_norm > GemmIndexedTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            //WriteDelimitedFile(D.LockedBuffer(), D.LDim(), D.Height(), D.Width(), "Derror.csv", 6);
            //WriteDelimitedFile(B.LockedBuffer(), B.LDim(), B.Height(), B.Width(), "Berror.csv", 6);
            break;
        }
    }

    // compute average runtimes per loop
    double t_copying = elapsed_subset * 0.001 / GemmIndexedTest::NUM_RUNS;
    double t_indexed = elapsed_indexed * 0.001 / GemmIndexedTest::NUM_RUNS;

    cout << endl;
    cout << "\t**** Results for Rank2 Indexed Sparse Gemm AB' Test ****" << endl;
    cout << endl;
    cout << "\t\t" << GemmIndexedTest::NUM_RUNS << " runs " << endl;
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << "\t\tElapsed time with data copying:\t" << t_copying << " ms." << endl;
    cout << "\t\tElapsed time with indexing:\t" << t_indexed << " ms." << endl;
    cout << endl;
    cout << "\t********************************************************" << endl;
    cout << endl;

    return (max_norm < GemmIndexedTest::MAX_ACCEPTABLE_FNORM);
}

//-----------------------------------------------------------------------------
bool TestSparseGemmIndexedBtA_Rank2(Random& rng)
{
    // C = alpha*B'*A + beta*C
    // 
    // Matrix A is sparse and is mxn, but use s cols from col index set S.
    // Matrix B has 2 cols; matrix C has 2 rows.

    Timer timer;
    unsigned int hw_max_threads = std::thread::hardware_concurrency();

    SparseMatrix<double> A, Asubset;
    elem::Matrix<double> B, Bsubset, C1, C2, C1subset, C2subset, D;

    // this variant of indexed gemm requires beta == 0
    double beta  = 0.0;
    double min_sparsity = 1.0, max_sparsity = 0.0;
    double max_norm = 0.0, elapsed_subset = 0.0, elapsed_indexed = 0.0;
    
    cout << "Running indexed GEMM rank2 B'A test..." << endl;

    for (unsigned int i=0; i<GemmIndexedTest::NUM_RUNS; ++i)
    {
        unsigned int height_a = rng.RandomRangeInt(512, GemmIndexedTest::MAX_DIM);
        unsigned int width_a  = rng.RandomRangeInt(512, GemmIndexedTest::MAX_DIM);

        unsigned int max_threads = 1 + (rng.RandomInt() % hw_max_threads);

        // sparse_fraction ranges from (0.1 - 0.09, 0.1 + 0.09)
        double sparse_fraction = rng.RandomDouble(0.1, 0.09);
        unsigned int nonzeros_per_col = sparse_fraction * height_a;
        if (0 == nonzeros_per_col)
            nonzeros_per_col = 1;        
        if (sparse_fraction < min_sparsity)
            min_sparsity = sparse_fraction;
        if (sparse_fraction > max_sparsity)
            max_sparsity = sparse_fraction;

        // generate a random alpha on [-1, 1]
        double alpha = rng.RandomDouble(0.0, 1.0);

        RandomSparseMatrix(A, nonzeros_per_col, rng, height_a, width_a);

        // set the dimensions of B and randomly initialize
        B.ResizeTo(height_a, 2);
        RandomMatrix(B, rng);
                
        // set the dimensions of C1 and C2 and randomly initialize
        C1.ResizeTo(2, width_a);
        RandomMatrix(C1, rng);
        C2 = C1;

        // generate a random set of column indices of A
        size_t index_set_size = rng.RandomRangeInt(1, width_a/2);

        // force rank2 some of the time
        if (0 == (i % 3))
            index_set_size = 2;

        std::vector<unsigned int> col_indices(index_set_size);
        RandomIndexSet(rng, width_a, index_set_size, col_indices);
        
        // compute the height of the most compact submatrix of A
        std::vector<unsigned int> old_to_new_rows, new_to_old_rows;
        A.SubMatrixRowsCompact(col_indices, old_to_new_rows, new_to_old_rows);
        unsigned int new_height = new_to_old_rows.size();

        // resize Bsubset and Csubset to match
        Bsubset.ResizeTo(new_height, 2);
        C1subset.ResizeTo(2, index_set_size);
        C2subset = C1subset;

        // extract B subset from B
        for (int c=0; c != Bsubset.Width(); ++c)
        {
            for (int r=0; r != Bsubset.Height(); ++r)
            {
                double val = B.Get(new_to_old_rows[r], c);
                Bsubset.Set(r, c, val);
            }
        }
        
        old_to_new_rows.clear();
        new_to_old_rows.clear();

        timer.Start();

        // extract a submatrix of A using the columns in col_indices
        A.SubMatrixColsCompact(Asubset, col_indices, old_to_new_rows, new_to_old_rows);
        assert(Asubset.Height() == new_height);
               
        // perform sparse Gemm on the subset
        Gemm(TRANSPOSE, NORMAL,
             alpha, Bsubset, Asubset, beta, C1subset, max_threads);

        timer.Stop();
        elapsed_subset += timer.ReportMicroseconds();

        timer.Start();

        // perform indexed Gemm 
        GemmIndexed(TRANSPOSE, NORMAL,
                    alpha, Bsubset, A, beta, C2subset, max_threads, 
                    new_height, old_to_new_rows, col_indices);

        timer.Stop();
        elapsed_indexed += timer.ReportMicroseconds();

        // check residuals
        
        // C2 = -C1 + C2
        elem::Axpy(-1.0, C1subset, C2subset);
        double fnorm = elem::Norm(C2subset, elem::FROBENIUS_NORM);
        cout << "\t[" << (i+1) << "/" << GemmIndexedTest::NUM_RUNS <<
            "]: fnorm of residual: " << fnorm << endl;
        if (fnorm > max_norm)
            max_norm = fnorm;
        
        if (max_norm > GemmIndexedTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            //WriteDelimitedFile(D.LockedBuffer(), D.LDim(), D.Height(), D.Width(), "Derror.csv", 6);
            //WriteDelimitedFile(B.LockedBuffer(), B.LDim(), B.Height(), B.Width(), "Berror.csv", 6);
            break;
        }
    }

    // compute average runtimes per loop
    double t_copying = elapsed_subset * 0.001 / GemmIndexedTest::NUM_RUNS;
    double t_indexed = elapsed_indexed * 0.001 / GemmIndexedTest::NUM_RUNS;

    cout << endl;
    cout << "\t**** Results for Rank2 Indexed Sparse Gemm B'A Test ****" << endl;
    cout << endl;
    cout << "\t\t" << GemmIndexedTest::NUM_RUNS << " runs " << endl;
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << "\t\tElapsed time with data copying:\t" << t_copying << " ms." << endl;
    cout << "\t\tElapsed time with indexing:\t" << t_indexed << " ms." << endl;
    cout << endl;
    cout << "\t********************************************************" << endl;
    cout << endl;

    return (max_norm < GemmIndexedTest::MAX_ACCEPTABLE_FNORM);
}
