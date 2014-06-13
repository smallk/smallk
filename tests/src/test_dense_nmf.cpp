#include <thread>
#include <limits>
#include <cassert>
#include <iostream>
#include "nmf.hpp"
#include "tests.hpp"
#include "random.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "matrix_generator.hpp"
#include "nmf_solver_hals.hpp"
#include "projected_gradient.hpp"
#include "file_loader.hpp"
#include "delimited_file.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

bool HalsTestDense(Random& rng);
bool ComparisonTest(Random& rng);
void UpdateW_Local(DenseMatrix<double>& W,
                   DenseMatrix<double>& WHHt_c,
                   DenseMatrix<double>& HHt,
                   const DenseMatrix<double>& AHt);

namespace NmfTest
{
    const unsigned int NUM_RUNS = 64;
    const double MAX_ACCEPTABLE_FNORM = 1.0e-8;
}

//-----------------------------------------------------------------------------
bool TestDenseNmf()
{
    Random rng;
    rng.SeedFromTime();

    //bool result1 = HalsTestDense(rng); return result1;
    bool result2 = ComparisonTest(rng);

    return result2;
}

//-----------------------------------------------------------------------------
bool HalsTestDense(Random& rng)
{
/*
    unsigned int m = 4;
    unsigned int n = 6;
    unsigned int k = 4;
    unsigned int max_threads = 4;

    DenseMatrix<double> WtWH_r, WHHt_c;
    WtWH_r.Resize(1, n);
    WHHt_c.Resize(m, 1);

    DenseMatrix<double> W_c;

    DenseMatrix<double> A(m, n), W(m, k), H(k, n), HHt(k, k), AHt(m, k);
    DenseMatrix<double> WtW(k, k), WtA(k, n), gradW(m, k), gradH(k, n);

    for (unsigned int c=0; c<n; ++c)
        for (unsigned int r=0; r<m; ++r)
            A.Set(r, c, c*m + r);

    for (unsigned int c=0; c<k; ++c)
        for (unsigned int r=0; r<m; ++r)
            W.Set(r, c, c*m + r);

    for (unsigned int c=0; c<n; ++c)
        for (unsigned int r=0; r<k; ++r)
            H.Set(r, c, c*m + r);
*/

    const unsigned int max_threads = 8;
    const unsigned int MAXITER = 5000;
    const double tol = 1.0e-4;
    std::string data_dir("/Users/richardboyd/Downloads/nmf_package/");

    bool ok;
    std::vector<double> buf_a;
    unsigned int m, n, ldim_a, ldim_w, ldim_h;

    // dense test
    unsigned int k = 16;
    std::string name_A("rnd_256_256.csv");
    ok = LoadDenseMatrix(data_dir + name_A, buf_a, m, n); assert(ok);
    ldim_a = m;
    DenseMatrix<double> A(m, n, &buf_a[0], ldim_a);
    cout << "m: " << m << ", n: " << n << ", k: " << k << endl;
    cout << "size of buf_a: " << buf_a.size() << endl;
    std::string name_W("w_init_256_16.csv");
    std::string name_H("h_init_16_256.csv");

    // // sparse test
    // unsigned int k=128;
    // SparseMatrix<double> A;
    // std::string name_A("reuters.mtx");
    // ok = LoadSparseMatrix(data_dir + name_A, A, m, n, nnz);
    // ldim_a = m;
    // cout << "m: " << m << ", n: " << n << ", k: " << k << endl;
    // std::string name_W("reuters_w_init_128.csv");
    // std::string name_H("reuters_h_init_128.csv");
    
    std::vector<double> buf_w(m*k), buf_h(k*n);
    ldim_w = m;
    ldim_h = k;

    ok = LoadDelimitedFile(buf_w, m, k, data_dir + name_W); assert(ok);
    ok = LoadDelimitedFile(buf_h, k, n, data_dir + name_H); assert(ok);
    
    DenseMatrix<double> W(m, k, &buf_w[0], ldim_w);
    DenseMatrix<double> H(k, n, &buf_h[0], ldim_h);
    DenseMatrix<double> HHt(k, k), WtW(k, k), AHt(m, k), WtA(k,n);
    DenseMatrix<double> gradW(m, k), gradH(k, n);
    DenseMatrix<double> WtWH_r(1, n), WHHt_c(m, 1);
    
    auto prev_precision = cout.precision(15);

    cout << "dimensions of A: " << A.Height() << " x " << A.Width() << endl;
    cout << "dimensions of W: " << W.Height() << " x " << W.Width() << endl;
    cout << "dimensions of H: " << H.Height() << " x " << H.Width() << endl;
    cout << "norm of W: " << Norm(W, FROBENIUS_NORM) << endl;
    cout << "norm of H: " << Norm(H, FROBENIUS_NORM) << endl;

    // compute HHt and AHt
    Gemm(NORMAL, TRANSPOSE, 1.0, H, H, 0.0, HHt);
    Gemm(NORMAL, TRANSPOSE, 1.0, A, H, 0.0, AHt, max_threads);

    std::cout << "norm of HHt: " << Norm(HHt, FROBENIUS_NORM) << std::endl;
    std::cout << "norm of AHt: " << Norm(AHt, FROBENIUS_NORM) << std::endl;

    double initgrad = 0.0, projnorm;

    for (unsigned int iter=1; iter<=MAXITER; ++iter)
    {
        cout << "iteration: " << iter << endl;
        
        UpdateW_Hals(W, WHHt_c, HHt, AHt);
        std::cout << "norm of W: " << Norm(W, FROBENIUS_NORM) << std::endl;
        
        // compute WtW and WtA
        Gemm(TRANSPOSE, NORMAL, 1.0, W, W, 0.0, WtW);
        Gemm(TRANSPOSE, NORMAL, 1.0, W, A, 0.0, WtA, max_threads);
        
        std::cout << "norm of WtW: " << Norm(WtW, FROBENIUS_NORM) << std::endl;
        std::cout << "norm of WtA: " << Norm(WtA, FROBENIUS_NORM) << std::endl;
        
        UpdateH_Hals(H, WtWH_r, WtW, WtA);
        
        std::cout << "norm of H: " << Norm(H, FROBENIUS_NORM) << std::endl;
        
        // compute gradH = WtW * H - WtA
        Gemm(NORMAL, NORMAL, 1.0, WtW, H, 0.0, gradH);  // gradH = WtW * H
        Axpy( -1.0, WtA, gradH);
        
        std::cout << "norm of gradH: " << Norm(gradH, FROBENIUS_NORM) << std::endl;    
        
        // compute the kxk matrix HHt = H * H'
        Gemm(NORMAL, TRANSPOSE, 1.0, H, H, 0.0, HHt);
        
        std::cout << "norm of HHt: " << Norm(HHt, FROBENIUS_NORM) << std::endl;
        
        // compute the mxk matrix AHt =  A * H'
        Gemm(NORMAL, TRANSPOSE, 1.0, A, H, 0.0, AHt, max_threads);
        
        std::cout << "norm of AHt: " << Norm(AHt, FROBENIUS_NORM) << std::endl;
        
        // compute gradW = W * HHt - AHt
        Gemm(NORMAL, NORMAL, 1.0, W, HHt, 0.0, gradW);  // gradW = W * HHt
        Axpy( -1.0, AHt, gradW);        
        
        std::cout << "norm of gradW: " << Norm(gradW, FROBENIUS_NORM) << std::endl;


        if (1 == iter)
        {
            initgrad = ProjectedGradientNorm(gradW, gradH, W, H);
            cout << "initgrad value: " << initgrad << endl;
            continue;
        }
        else
        {
            projnorm = ProjectedGradientNorm(gradW, gradH, W, H);
        }

        double progress_metric = projnorm / initgrad;
        cout << "\tprogress metric: " << progress_metric << endl;
        
        if (projnorm < tol*initgrad)
        {
            cout << "converged on iteration " << iter << endl;
            break;
        }

    }

    cout.precision(prev_precision);
    return true;
}

//-----------------------------------------------------------------------------
bool ComparisonTest(Random& rng)
{
    SparseMatrix<double> A;
    DenseMatrix<double> D, W, W2, H, H2, diff, Wsave, Hsave;
    double max_norm = 0.0, min_sparsity = 1.0, max_sparsity = 0.0;

    // Generate a random sparse matrix and its fully-dense equivalent, solve
    // both NMF problems, and compare results.  Some of the randomly-generated
    // sparse matrices may be rank-deficient, in which case the factorization
    // is likely to fail.  For algorithms such as BPP, which requires the
    // input matrix to have full rank, such failures do not indicate a problem
    // with the NMF code.  The BPP solver will write a message to cerr 
    // indicating rank-deficiency.  The NMF code has a problem if the 
    // input matrix has full rank and the comparison test fails.

    for (unsigned int i=0; i<NmfTest::NUM_RUNS; ++i)
    {
        unsigned int height    = rng.RandomRangeInt(64, 512);
        unsigned int width     = rng.RandomRangeInt(64, 512);
        unsigned int k         = rng.RandomRangeInt(16, 64);

        // sparse_fraction ranges from (0.1 - 0.09, 0.1 + 0.09)
        double sparse_fraction = rng.RandomDouble(0.1, 0.09);
        unsigned int nonzeros_per_col = sparse_fraction * height;
        if (0 == nonzeros_per_col)
            nonzeros_per_col = 1;
        
        if (sparse_fraction < min_sparsity)
            min_sparsity = sparse_fraction;
        if (sparse_fraction > max_sparsity)
            max_sparsity = sparse_fraction;

        // generate a random sparse matrix and its fully dense equivalent
        RandomSparseMatrix(rng, A, nonzeros_per_col, height, height, width, width);
        D = MakeDense(A);
        
        // randomly initialize W and H
        W.Resize(height, k);
        W2.Resize(height, k);
        Wsave.Resize(height, k);
        H.Resize(k, width);
        H2.Resize(k, width);
        Hsave.Resize(k, width);

        RandomMatrix(W.Buffer(), W.LDim(), W.Height(), W.Width(), rng);
        RandomMatrix(H.Buffer(), H.LDim(), H.Height(), H.Width(), rng);

        // save the initial W and H
        Wsave = W;
        Hsave = H;
        
        // copy for the dense run
        W2 = W;
        H2 = H;
        
        NmfOptions opts;
        opts.tol = 0.005;

        // Switch between the MU, RANK2, and BPP algorithms.  HALS is 
        // EXTREMELY sensitive to the epsilon value; the sparse and dense 
        // code need a HALS epsilon of 1.0e-8 to show agreement to 1.0e-6.
        // We are using 1.0e-13 for other reasons, so skip HALS in this test.

        int q = i % 3;
        if (0 == q)
            opts.algorithm = NmfAlgorithm::MU;
        //else if (1 == q)
        //opts.algorithm = NmfAlgorithm::HALS;
        else if (1 == q)
        {
            opts.algorithm = NmfAlgorithm::RANK2;
            k = 2;
        }
        else
            opts.algorithm = NmfAlgorithm::BPP;

        opts.prog_est_algorithm = NmfProgressAlgorithm::DELTA_FNORM;

        opts.height      = height;
        opts.width       = width;
        opts.k           = (NmfAlgorithm::RANK2 == opts.algorithm ? 2 : k);
        opts.min_iter    = 5;
        opts.max_iter    = 10000;
        opts.tolcount    = 1;

        int hw_threads = std::thread::hardware_concurrency();
        if (0 == hw_threads)
            opts.max_threads = 2;
        else
            opts.max_threads = hw_threads;
        opts.verbose     = false;
        opts.normalize   = true;
        
        NmfStats stats;
        NmfResult result_sparse, result_dense;

        // run both sparse and dense NMF
        result_sparse = NmfSparse(opts, 
                                  A.Height(), A.Width(), A.Size(),
                                  A.LockedColBuffer(), 
                                  A.LockedRowBuffer(), 
                                  A.LockedDataBuffer(),
                                  W.Buffer(), W.LDim(),
                                  H.Buffer(), H.LDim(),
                                  stats);
        if (NmfResult::OK != result_sparse)
        {
            cerr << "\tSkipping this rank-deficient problem..." << endl;
            continue;
        }

        result_dense = Nmf(opts, D.Buffer(), D.LDim(), 
                           W2.Buffer(), W2.LDim(), H2.Buffer(), H2.LDim(), stats);
        if (NmfResult::OK != result_dense)
            cerr << "TestDenseNmf: dense nmf computation failed" << endl;
        
        std::string algorithm;
        if (NmfAlgorithm::MU == opts.algorithm)
            algorithm = "MU";
        else if (NmfAlgorithm::RANK2 == opts.algorithm)
            algorithm = "Rank2";
        else if (NmfAlgorithm::BPP == opts.algorithm)
            algorithm = "BPP";
        else if (NmfAlgorithm::HALS == opts.algorithm)
            algorithm = "HALS";

        // check residuals for W and H
        
        // W2 = -W + W2
        Axpy(-1.0, W, W2);
        double fnorm_w = Norm(W2, FROBENIUS_NORM);
        if (fnorm_w > max_norm)
            max_norm = fnorm_w;
        
        // H2 = -H + H2
        Axpy(-1.0, H, H2);
        double fnorm_h = Norm(H2, FROBENIUS_NORM);
        if (fnorm_h > max_norm)
            max_norm = fnorm_h;

        cout << "[" << setw(4) << (i+1) << "/" << setw(4) << NmfTest::NUM_RUNS 
             <<"] m: " << setw(6) << height << ", n: " << setw(6) << width 
             << ", k: " << setw(3) << k << ", algorithm: " << setw(5) << algorithm 
             << ", max norm of residuals: " << std::max(fnorm_w, fnorm_h) << endl;

        if (max_norm > NmfTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            WriteDelimitedFile(Wsave.LockedBuffer(), Wsave.LDim(), 
                               Wsave.Height(), Wsave.Width(), "Werror.csv", 6);
            WriteDelimitedFile(Hsave.LockedBuffer(), Hsave.LDim(), 
                               Hsave.Height(), Hsave.Width(), "Herror.csv", 6);
            WriteDelimitedFile(D.LockedBuffer(), D.LDim(), 
                               D.Height(), D.Width(), "Derror.csv", 6);
            break;
        }
    }

    cout << endl;
    cout << "\t************** Results for Dense NMF Test *************" << endl;
    cout << endl;
    cout << "\t\t" << NmfTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_sparsity << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_sparsity << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*******************************************************" << endl;
    cout << endl;

    return (max_norm < NmfTest::MAX_ACCEPTABLE_FNORM);
};
//-----------------------------------------------------------------------------
void UpdateW_Local(DenseMatrix<double>& W,
                   DenseMatrix<double>& WHHt_c,
                   DenseMatrix<double>& HHt,        // kxk matrix H * H'
                   const DenseMatrix<double>& AHt)  // mxk matrix A * H'
{
    // these are views
    DenseMatrix<double> HHt_c, W_c;

    cout << "\t\tUpdateLocal: norm of W: " << Norm(W, FROBENIUS_NORM) << endl;
    cout << "\t\tUpdateLocal: norm of HHt: " << Norm(HHt, FROBENIUS_NORM) << endl;
    cout << "\t\tUpdateLocal: norm of AHt: " << Norm(AHt, FROBENIUS_NORM) << endl;

    int zero_col_count = 0;
    for (int c=0; c<W.Width(); ++c)
    {
        // create a view of column c of matrix HHt
        View(HHt_c, HHt, 0, c, HHt.Height(), 1);

        // compute the vector WHHt_c = W * HHt_c
        Gemv(NORMAL, 1.0, W, HHt_c, 0.0, WHHt_c);

        // diagonal element of HHt
        double HHt_cc = HHt.Get(c, c);

        cout << "\t\tUpdateLocal: diagonal element: " << HHt_cc << endl;

        int num_zeros = 0;
        for (int r=0; r<W.Height(); ++r)
        {
            double w    = W.Get(r, c);
            double aht  = AHt.Get(r, c);
            double whht = WHHt_c.Get(r, 0);

            w = w + (aht - whht) / HHt_cc;
            if (std::isnan(w) || (w < 0.0))
            {
                w = 0.0;
                ++num_zeros;
            }

            W.Set(r, c, w);
        }

        if (W.Height() == num_zeros)
        {
            // column contains all zeros
            ++zero_col_count;
            for (int r=0; r<W.Height(); ++r)
                W.Set(r, c, std::numeric_limits<double>::epsilon());            
        }

        // view of column c of matrix W
        View(W_c, W, 0, c, W.Height(), 1);
        double norm = Norm(W_c, FROBENIUS_NORM);
        double inv_norm = 1.0 / norm;
        Scal(inv_norm, W_c);
    }

    std::cout << "Found " << zero_col_count << " zero columns." << std::endl;
}
