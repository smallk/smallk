#include <thread>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include "tests.hpp"
#include "timer.hpp"
#include "random.hpp"
#include "normal_eq.hpp"
#include "bit_matrix.hpp"
#include "file_loader.hpp"
#include "dense_matrix.hpp"
#include "thread_utils.hpp"
#include "vector_utils.hpp"
#include "sparse_matrix.hpp"
#include "bit_matrix_ops.hpp"
#include "nmf_solver_bpp.hpp"
#include "matrix_generator.hpp"
#include "fisher_yates_shuffle.hpp"
#include "delimited_file.hpp"

// Test various operations needed for the NMF-BPP implementation.

using std::cout;
using std::cerr; 
using std::endl;
using std::setw;

bool TestBitMatrix(Random& rng);
bool SolveNormalEqComb(Random& rng);
bool TestNormalEqSolver(Random& rng);
bool TestNormalEqSolverLeft(Random& rng);
bool TestNnlsBlockpivot(Random& rng, const std::string& data_dir);

void RandomPassiveSet(BitMatrix& passive_set,
                      Random& rng, 
                      const unsigned int height,
                      const unsigned int width);
    
namespace BppTest
{
    const unsigned int NUM_RUNS = 32;
    const double MAX_ACCEPTABLE_FNORM = 1.0e-10;
}

//-----------------------------------------------------------------------------
bool TestBpp(const std::string& data_dir)
{
    Random rng;
    rng.SeedFromTime();
    //rng.SeedFromInt(42);

    bool result1 = TestNnlsBlockpivot(rng, data_dir);
    bool result2 = TestNormalEqSolver(rng);
    bool result3 = TestNormalEqSolverLeft(rng);

    return (result1 && result2 && result3);
}

//-----------------------------------------------------------------------------
bool TestNnlsBlockpivot(Random& rng, const std::string& data_dir)
{
    // Run the NNLS-BPP solver and compare results with Matlab.

    unsigned int kvals[] = {16, 33};
    std::string matrices[] = {"reuters.mtx", "test/reduced_matrix_20news.mtx"};
    std::string h_inits[] = {"h_init_k_16.csv", "h_init_k_33.csv"};
    std::string w_inits[] = {"w_init_k_16.csv", "w_init_k_33.csv"};
    std::string matlab_h[] = {"nnls_bpp_reuters_k_16_matlab_H_result.csv",
                            "nnls_bpp_20news_k_33_matlab_H_result.csv"};
    std::string matlab_gradH[] = {"nnls_bpp_reuters_k_16_matlab_gradH_result.csv",
                                "nnls_bpp_20news_k_33_matlab_gradH_result.csv"};
    Timer timer;
    SparseMatrix<double> A;
    double max_norm = 0.0;
    SetMaxThreadCount(std::thread::hardware_concurrency());

    for (unsigned int i=0; i<2; ++i)
    {
        cout << "Running NNLS-BPP Matlab comparison test on file: data/" 
             << matrices[i] << ", "  << "with k = " << kvals[i] << "." << endl;

        unsigned int k = kvals[i];
        unsigned int height_a = 0, width_a = 0, nnz_a = 0;

        std::string input_file = data_dir + matrices[i];
        if (!LoadSparseMatrix(input_file, A, height_a, width_a, nnz_a))
        {
            cerr << "could not load file" << input_file << endl;
            return false;
        }

        unsigned int m = height_a;
        unsigned int n = width_a;
        std::vector<double> buf_w(m*k), buf_h(k*n), buf_a(m*n);
        DenseMatrix<double> W((int)m, (int)k), H((int)k, (int)n), Winit, Hinit;
        DenseMatrix<double> WtW((int)k, (int)k), WtA((int)k, (int)n), gradH((int)k, (int)n);

        std::string test_dir = data_dir + std::string("test/");
        if ( (!LoadDelimitedFile(buf_h, k, n, test_dir + h_inits[i])) ||
             (!LoadDelimitedFile(buf_w, m, k, test_dir + w_inits[i])))
        {
            cerr << "could not load initializer matrices for W and H" << endl;
            return false;
        }
    
        Winit.Attach((int)m, (int)k, &buf_w[0], (int)m);
        Hinit.Attach((int)k, (int)n, &buf_h[0], (int)k);
        
        double elapsed_ms = 0.0;        
        timer.Start();
        
        for (unsigned int q=0; q<BppTest::NUM_RUNS; ++q)
        {
            W = Winit;
            H = Hinit;
            
            // compute WtW and WtA
            Gemm(TRANSPOSE, NORMAL, 1.0, W, W, 0.0, WtW);
            Gemm(TRANSPOSE, NORMAL, 1.0, W, A, 0.0, WtA);

            // solve the nonnegative least squares problem using BPP
            NnlsBlockpivot(WtW, WtA, H, gradH);
        }
        
        timer.Stop();
        elapsed_ms = timer.ReportMilliseconds();
        cout << "Avg elapsed time: " << elapsed_ms/BppTest::NUM_RUNS << " ms." << endl;
        cout << endl;
        
        // load the Matlab results
        std::vector<double> matlab_h_buf(k*n), matlab_gradH_buf(k*n);
        if (!LoadDelimitedFile(matlab_h_buf, k, n, test_dir + matlab_h[i]))
        {
            cerr << "could not load Matlab H result" << endl;
            return false;
        }
        if (!LoadDelimitedFile(matlab_gradH_buf, k, n, test_dir + matlab_gradH[i]))
        {
            cerr << "could not load Matlab gradH result" << endl;
            return false;
        }
        
        DenseMatrix<double> H_matlab    ( (int)k, (int)n, &matlab_h_buf[0], (int)k);
        DenseMatrix<double> gradH_matlab( (int)k, (int)n, &matlab_gradH_buf[0], (int)k);
        
        // H_matlab = H_matlab - H
        Axpy(-1.0, H, H_matlab);
        
        // gradH_matlab = gradH_matlab - gradH
        Axpy(-1.0, gradH, gradH_matlab);
        
        double norm = Norm(H_matlab, FROBENIUS_NORM);
        max_norm = std::max(norm, max_norm);
        norm = Norm(gradH_matlab, FROBENIUS_NORM);
        max_norm = std::max(norm, max_norm);
    }
    
    cout << endl;
    cout << "\t**** Results for NNLS BPP solver tests ****" << endl;
    cout << endl;
    cout << "\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*******************************************" << endl;
    cout << endl;

    return true;
}

//-----------------------------------------------------------------------------
bool TestNormalEqSolver(Random& rng)
{
    // This test generates random normal equations to be solved by the BPP
    // solver.  Diagonal dominance is enforced on the cross-product matrix to
    // ensure that it is nonsingular.  A random passive set is also generated,
    // and it is designed to always have at least a single '1' entry in any
    // column.  The normal equation system is solved twice, one time with
    // column grouping and one time without.  The solutions are compared and
    // the residual norms printed out. 
    
    BitMatrix passive_set;
    SparseMatrix<double> A;
    DenseMatrix<double> W, WtW, WtA, X, X2;
    SetMaxThreadCount(std::thread::hardware_concurrency());
    double max_norm = 0.0, min_occupancy = 1.0, max_occupancy = 0.0;

    std::vector<unsigned int> col_indices;

    cout << "Running BPP normal equation solver test..." << endl;

    for (unsigned int i=0; i<BppTest::NUM_RUNS; ++i)
    {
        unsigned int m = rng.RandomRangeInt(16,  768);
        unsigned int n = rng.RandomRangeInt(64, 1024);

        // confine k to [4, 256], but also no larger than min(m, n)
        unsigned int s = std::min(m, n);
        unsigned int k = rng.RandomRangeInt( 4,  std::min(s, 256u));

        // force n == 1 occasionally
        if (0 == (i % 5))
            n = 1;

        col_indices.resize(n);
        for (unsigned int q=0; q<n; ++q)
            col_indices[q] = q;

        // generate random W
        W.Resize(m, k);
        RandomMatrix(W.Buffer(), W.LDim(), W.Height(), W.Width(), rng);

        // generate random sparse A with occupancy (0.1-0.09, 0.1+0.09)
        double occupancy = rng.RandomDouble(0.1, 0.09);
        unsigned int nz_per_col = occupancy * m;
        if (0 == nz_per_col)
            nz_per_col = 1;
        min_occupancy = std::min(min_occupancy, occupancy);
        max_occupancy = std::max(max_occupancy, occupancy);
        RandomSparseMatrix(rng, A, nz_per_col, m, m, n, n);

        // compute WtW and WtA
        WtW.Resize(k, k);
        WtA.Resize(k, n);
        Gemm(TRANSPOSE, NORMAL, 1.0, W, W, 0.0, WtW);
        Gemm(TRANSPOSE, NORMAL, 1.0, W, A, 0.0, WtA);

        // make WtW diagonally dominant to guarantee nonsingular
        MakeDiagonallyDominant(WtW);

        // generate a random passive set with random numbers of identical cols
        passive_set.Resize(k, n);
        RandomPassiveSet(passive_set, rng, k, n);

        // X will contain the column-grouped solution, X2 will contain the
        // solution with each column solved independently
        X.Resize(k, n);
        X2.Resize(k, n);
        
        // group columns and solve the linear systems
        BppSolveNormalEq(col_indices, passive_set, WtW, WtA, X);

        // solve again, but this time with no column grouping
        BppSolveNormalEqNoGroup(col_indices, passive_set, WtW, WtA, X2);

        // compute difference matrix X = X - X2
        Axpy(-1.0, X2, X);

        // check norm
        double norm = Norm(X, FROBENIUS_NORM);
        max_norm = std::max(norm, max_norm);
        if (norm > 1.0e-6)
        {
            // error - the norm of the difference matrix should be tiny
            WriteDelimitedFile(WtW.LockedBuffer(), WtW.LDim(), 
                               WtW.Height(), WtW.Width(), "LHS_error.csv", 6);
            WriteDelimitedFile(WtA.LockedBuffer(), WtA.LDim(), 
                               WtA.Height(), WtA.Width(), "RHS_error.csv", 6);
            cout << "**** ERROR EXIT ****" << endl;
            break;
        }

        cout << "[" << setw(4) << i << "/" << setw(4) << BppTest::NUM_RUNS 
             <<"] m: " << setw(6) << m << ", n: " << setw(6) << n << ", k: " 
             << setw(3) << k << ", norm of residual: " << norm << endl;
    }

    cout << endl;
    cout << "\t**** Results for BPP Normal Eq. Solver Test ****" << endl;
    cout << endl;
    cout << "\t\t" << BppTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_occupancy << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_occupancy << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*************************************************" << endl;
    cout << endl;

    return (max_norm < BppTest::MAX_ACCEPTABLE_FNORM);
}

//-----------------------------------------------------------------------------
bool TestNormalEqSolverLeft(Random& rng)
{
    // This test generates random normal equations to be solved by the BPP
    // solvers.  Diagonal dominance is enforced on the cross-product matrix to
    // ensure that it is nonsingular.  A random passive set is also generated,
    // and it is designed to always have at least a single '1' entry in any
    // column.  The normal equation system is solved twice, one time using
    // the 'left-oriented' solver, and one time using the standard solver.
    // The solutions are compared and the residual norms printed out. 
    
    BitMatrix passive_set_1, passive_set_2;
    SparseMatrix<double> A;
    DenseMatrix<double> H, HHt, AHt, AHtt, X, Xt, temp, PSG, PSGt;
    SetMaxThreadCount(std::thread::hardware_concurrency());
    double max_norm = 0.0, min_occupancy = 1.0, max_occupancy = 0.0;

    std::vector<unsigned int> row_indices, col_indices;

    for (unsigned int i=0; i<BppTest::NUM_RUNS; ++i)
    {
        unsigned int m = rng.RandomRangeInt(16,  768);
        unsigned int n = rng.RandomRangeInt(64, 1024);

        // confine k to [4, 256], but also no larger than min(m, n)
        unsigned int s = std::min(m, n);
        unsigned int k = rng.RandomRangeInt( 4,  std::min(s, 256u));

        // force m == 1 occasionally (single row)
        if (0 == (i % 5))
            m = 1;

        row_indices.resize(m);
        for (unsigned int q=0; q<m; ++q)
            row_indices[q] = q;

        col_indices.resize(m);
        for (unsigned int q=0; q<m; ++q)
            col_indices[q] = q;

        // generate random H
        H.Resize(k, n);
        RandomMatrix(H.Buffer(), H.LDim(), H.Height(), H.Width(), rng);

        // generate random sparse A with occupancy (0.1-0.09, 0.1+0.09)
        double occupancy = rng.RandomDouble(0.1, 0.09);
        unsigned int nz_per_col = occupancy * m;
        if (0 == nz_per_col)
            nz_per_col = 1;
        min_occupancy = std::min(min_occupancy, occupancy);
        max_occupancy = std::max(max_occupancy, occupancy);
        RandomSparseMatrix(rng, A, nz_per_col, m, m, n, n);

        // compute HHt AHt, and (AHt)'
        HHt.Resize(k, k);
        AHt.Resize(m, k);
        AHtt.Resize(k, m);
        Gemm(NORMAL, TRANSPOSE, 1.0, H, H, 0.0, HHt);
        Gemm(NORMAL, TRANSPOSE, 1.0, A, H, 0.0, AHt);
        Transpose(AHt, AHtt);

        // make HHt diagonally dominant to guarantee nonsingular
        MakeDiagonallyDominant(HHt);

        // X will contain the column-grouped solution, X2 will contain the
        // solution with each column solved independently
        X.Resize(m, k);
        Xt.Resize(k, m);

        // load X and X'
        RandomMatrix(X.Buffer(), X.LDim(), X.Height(), X.Width(), rng);
        Transpose(X, Xt);

        // generate a passive set and the corresponding col and row indices
        PSG.Resize(m, k);
        RandomMatrix(PSG.Buffer(), PSG.LDim(), PSG.Height(), PSG.Width(), rng, 0.4, 0.45);
        passive_set_1.Resize(m, k);
        passive_set_1 = (PSG > 0.0);

        PSGt.Resize(k, m);
        Transpose(PSG, PSGt);
        passive_set_2.Resize(k, m);
        passive_set_2 = (PSGt > 0.0);

        // solve X * HHt = AHt   [mxk][kxk] = [mxk]
        //SolveNormalEqLeft(HHt, AHt, X);
        if (!BppSolveNormalEqLeftNoGroup(row_indices, passive_set_1, HHt, AHt, X))
        {
            cerr << "\tRandom matrix was rank-deficient; skipping..." << endl;
            continue;
        }

        // solve HHt * X' = (AHt)' with the standard solver  [kxk][kxm] = [kxm]
        //SolveNormalEq(HHt, AHtt, Xt);
        if (!BppSolveNormalEqNoGroup(col_indices, passive_set_2, HHt, AHtt, Xt))
        {
            cerr << "\tRandom matrix was rank-deficient; skipping..." << endl;
            continue;
        }

        // check the residual norm
        temp.Resize(m, k);
        Transpose(Xt, temp);
        Axpy(-1.0, X, temp);
        double norm = Norm(temp, FROBENIUS_NORM);
        max_norm = std::max(norm, max_norm);

        cout << "[" << setw(4) << i << "/" << setw(4) << BppTest::NUM_RUNS 
             <<"] m: " << setw(6) << m << ", n: " << setw(6) << n << ", k: " 
             << setw(3) << k << ", norm of residual: " << norm << endl;
    }
    
    cout << endl;
    cout << "\t****** Results for TestNormalEqSolverLeft *******" << endl;
    cout << endl;
    cout << "\t\t" << BppTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_occupancy << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_occupancy << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*************************************************" << endl;
    cout << endl;

    return (max_norm < BppTest::MAX_ACCEPTABLE_FNORM);
}

//-----------------------------------------------------------------------------
void RandomPassiveSet(BitMatrix& passive_set, 
                      Random& rng,
                      const unsigned int height,
                      const unsigned int width)
{
    passive_set.Resize(height, width);

    unsigned int* buf_ps = passive_set.Buffer();
    const unsigned int ldim = passive_set.LDim();
    const unsigned int full_wds = height / BitMatrix::BITS_PER_WORD;
    const unsigned int extra = height - BitMatrix::BITS_PER_WORD*full_wds;
    const unsigned int MASK  = passive_set.Mask();

    // initial nonzero bit pattern to fill columns with
    unsigned int pattern = 0x01;

    // randomly shuffle the column indices of the passive set
    std::vector<unsigned int> col_indices(width);
    for (unsigned int q=0; q<width; ++q)
        col_indices[q] = q;
    
    if (width > 1)
        FisherYatesShuffle(col_indices.begin(), col_indices.end(), rng);
    
    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int col_index  = col_indices[c];
        unsigned int col_offset = col_index * ldim;
        
        // Fill this column with the bit pattern, ensuring that the
        // final bits in each column (outside the MASK for the passive
        // set matrix) are zero.
        unsigned int r_wd=0;
        for (; r_wd < full_wds; ++r_wd)
            buf_ps[col_offset + r_wd] = pattern;
        if (extra > 0)
        {
            unsigned int wd = MASK & pattern;
            assert(wd > 0);
            buf_ps[col_offset + r_wd] = wd;
        }
        
        // shift the pattern to the left one bit and OR in a new bit
        pattern <<= 1;
        pattern |= 1;
        
        // start over if all ones
        if (0xFFFFFFFF == pattern)
            pattern = 0x01;            
    }

    // None of the columns should contain all zero bits.
    std::vector<unsigned int> col_sums(width);
    passive_set.SumColumns(col_sums);
    for (unsigned int c=0; c<width; ++c)
    {
        if (0 == col_sums[c])
        {
            unsigned int col_offset = c*ldim;
            cout << "RandomPassiveSet: column " << c << " is a zero column." << endl;
            cout << "Pattern: " << std::hex << pattern << std::dec << endl;
            cout << "Height: " << height << ", width: " << width << endl;
            cout << "full_wds: " << full_wds << endl;
            cout << "extra: " << extra << endl;
            cout << "Ldim: " << ldim << ", extra: " << extra << endl;
            cout << "MASK: " << std::hex << MASK << std::dec << endl;
            cout << "col_offset: " << col_offset << endl;
            unsigned int r_wd = 0;
            for (; r_wd<full_wds; ++r_wd)
                cout << "words_[" << r_wd << "]: " << std::hex 
                     << buf_ps[col_offset + r_wd] << std::dec << endl;
            if (MASK > 0)
            {
                unsigned int wd = MASK & buf_ps[col_offset + r_wd];
                cout << "words_[" << r_wd << "]: " << std::hex << wd << std::dec << endl;
            }

            assert(0 != col_sums[c]);
        }
    }
}
