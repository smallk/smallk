#include <thread>
#include <iostream>
#include "random.hpp"
#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"
#include "matrix_generator.hpp"
#include "delimited_file.hpp"
#include "thread_utils.hpp"

using std::cout;
using std::cerr;
using std::endl;

bool TestSparseGemmAB  (Random& rng);
bool TestSparseGemmABt (Random& rng);
bool TestSparseGemmBA  (Random& rng);
bool TestSparseGemmBtA (Random& rng);

namespace GemmTest
{
const unsigned int NUM_RUNS = 256;
const double MAX_ACCEPTABLE_FNORM = 1.0e-10;
}

//-----------------------------------------------------------------------------
bool TestSparseGemm()
{
    Random rng;
    rng.SeedFromTime();

    bool result1 = TestSparseGemmAB(rng);
    bool result2 = TestSparseGemmABt(rng);
    bool result3 = TestSparseGemmBA(rng);
    bool result4 = TestSparseGemmBtA(rng);

    return (result1 && result2 && result3 && result4);
}

//-----------------------------------------------------------------------------
bool TestSparseGemmAB(Random& rng)
{
    // C = alpha*A*B + beta*C

    SparseMatrix<double> A;
    DenseMatrix<double> D, B, C1, C2, diff;

    double alpha, beta;
    double max_norm = 0.0, min_sparsity = 1.0, max_sparsity = 0.0;

    unsigned int hw_max_threads = std::thread::hardware_concurrency();
    if (0 == hw_max_threads)
        hw_max_threads = 2;

    for (unsigned int i=0; i<GemmTest::NUM_RUNS; ++i)
    {
        unsigned int max_threads = 1 + (rng.RandomInt() % hw_max_threads);
        SetMaxThreadCount(max_threads);

        unsigned int height_a = rng.RandomRangeInt(64, 512);
        unsigned int width_a  = rng.RandomRangeInt(64, 512);
        unsigned int width_b  = rng.RandomRangeInt(16, 768);

        // need to force rank2 occasionally (B.Width() == C.Width() == 2)
        if (0 == (i % 3))
            width_b = 2;

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
        alpha = rng.RandomDouble(0.0, 1.0);
        beta  = rng.RandomDouble(0.0, 1.0);

        // generate a random sparse matrix and its fully dense equivalent
        RandomSparseMatrix(rng, A, nonzeros_per_col, height_a, height_a, width_a, width_a);
        D = MakeDense(A);
        
        // set the dimensions of B and randomly initialize
        B.Resize(width_a, width_b);
        RandomMatrix(B.Buffer(), B.LDim(), B.Height(), B.Width(), rng);
        
        // set the dimensions of C1 and C2 and randomly initialize
        C1.Resize(height_a, width_b);
        RandomMatrix(C1.Buffer(), C1.LDim(), C1.Height(), C1.Width(), rng);
        C2 = C1;

        // sparse Gemm
        Gemm(NORMAL, NORMAL, alpha, A, B, beta, C1);
        
        // dense Gemm
        Gemm(NORMAL, NORMAL, alpha, D, B, beta, C2);

        // check residuals
        
        // C2 = -C1 + C2
        Axpy(-1.0, C1, C2);
        double fnorm = Norm(C2, FROBENIUS_NORM);
        //cout << "[" << i << "]: fnorm of C residual: " << fnorm << endl;
        if (fnorm > max_norm)
            max_norm = fnorm;
        
        if (max_norm > GemmTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            //WriteDelimitedFile(D.LockedBuffer(), D.LDim(), D.Height(), D.Width(), "Derror.csv", 6);
            //WriteDelimitedFile(B.LockedBuffer(), B.LDim(), B.Height(), B.Width(), "Berror.csv", 6);
            break;
        }
    }

    cout << endl;
    cout << "\t******** Results for Sparse Gemm AB Test ********" << endl;
    cout << endl;
    cout << "\t\t" << GemmTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_sparsity << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_sparsity << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*************************************************" << endl;
    cout << endl;

    return (max_norm < GemmTest::MAX_ACCEPTABLE_FNORM);
}

//-----------------------------------------------------------------------------
bool TestSparseGemmABt(Random& rng)
{
    // C = alpha*A*B' + beta*C

    SparseMatrix<double> A;
    DenseMatrix<double> D, B, C1, C2, diff;

    double alpha, beta;
    double max_norm = 0.0, min_sparsity = 1.0, max_sparsity = 0.0;

    unsigned int hw_max_threads = std::thread::hardware_concurrency();
    if (0 == hw_max_threads)
        hw_max_threads = 2;

    for (unsigned int i=0; i<GemmTest::NUM_RUNS; ++i)
    {
        unsigned int max_threads = 1 + (rng.RandomInt() % hw_max_threads);
        SetMaxThreadCount(max_threads);

        unsigned int height_a = rng.RandomRangeInt(64, 512);
        unsigned int width_a  = rng.RandomRangeInt(64, 512);
        unsigned int width_bt  = rng.RandomRangeInt(16, 768);

        // need to force rank2 occasionally (B.Width() == C.Width() == 2)
        if (0 == (i % 3))
            width_bt = 2;

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
        alpha = rng.RandomDouble(0.0, 1.0);
        beta  = rng.RandomDouble(0.0, 1.0);

        // generate a random sparse matrix and its fully dense equivalent
        RandomSparseMatrix(rng, A, nonzeros_per_col, height_a, height_a, width_a, width_a);
        D = MakeDense(A);
        
        // set the dimensions of B and randomly initialize
        B.Resize(width_bt, width_a);
        RandomMatrix(B.Buffer(), B.LDim(), B.Height(), B.Width(), rng);
        
        // set the dimensions of C1 and C2 and randomly initialize
        C1.Resize(height_a, width_bt);
        RandomMatrix(C1.Buffer(), C1.LDim(), C1.Height(), C1.Width(), rng);
        C2 = C1;

        // sparse Gemm
        Gemm(NORMAL, TRANSPOSE, alpha, A, B, beta, C1);
        
        // dense Gemm
        Gemm(NORMAL, TRANSPOSE, alpha, D, B, beta, C2);

        // check residuals
        
        // C2 = -C1 + C2
        Axpy(-1.0, C1, C2);
        double fnorm = Norm(C2, FROBENIUS_NORM);
        //cout << "[" << i << "]: fnorm of C residual: " << fnorm << endl;
        if (fnorm > max_norm)
            max_norm = fnorm;
        
        if (max_norm > GemmTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            //WriteDelimitedFile(D.LockedBuffer(), D.LDim(), D.Height(), D.Width(), "Derror.csv", 6);
            //WriteDelimitedFile(B.LockedBuffer(), B.LDim(), B.Height(), B.Width(), "Berror.csv", 6);
            break;
        }
    }

    cout << endl;
    cout << "\t******** Results for Sparse Gemm AB' Test ********" << endl;
    cout << endl;
    cout << "\t\t" << GemmTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_sparsity << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_sparsity << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t**************************************************" << endl;
    cout << endl;

    return (max_norm < GemmTest::MAX_ACCEPTABLE_FNORM);
}

//-----------------------------------------------------------------------------
bool TestSparseGemmBA(Random& rng)
{
    // C = alpha*B*A + beta*C

    SparseMatrix<double> A;
    DenseMatrix<double> D, B, C1, C2, diff;

    double alpha, beta;
    double max_norm = 0.0, min_sparsity = 1.0, max_sparsity = 0.0;

    unsigned int hw_max_threads = std::thread::hardware_concurrency();
    if (0 == hw_max_threads)
        hw_max_threads = 2;

    for (unsigned int i=0; i<GemmTest::NUM_RUNS; ++i)
    {
        unsigned int max_threads = 1 + (rng.RandomInt() % hw_max_threads);
        SetMaxThreadCount(max_threads);

        unsigned int height_a = rng.RandomRangeInt(64, 512);
        unsigned int width_a  = rng.RandomRangeInt(64, 512);
        unsigned int height_b = rng.RandomRangeInt(16, 768);

        // generate a random alpha and beta on [-1, 1]
        alpha = rng.RandomDouble(0.0, 1.0);
        beta  = rng.RandomDouble(0.0, 1.0);

        // need to force rank2 occasionally (C.Width() == A.Width() == 2)
        if (0 == (i % 3))
            width_a = 2;

        // sparse_fraction ranges from (0.1 - 0.09, 0.1 + 0.09)
        double sparse_fraction = rng.RandomDouble(0.1, 0.09);
        unsigned int nonzeros_per_col = sparse_fraction * height_a;
        if (0 == nonzeros_per_col)
            nonzeros_per_col = 1;
        
        if (sparse_fraction < min_sparsity)
            min_sparsity = sparse_fraction;
        if (sparse_fraction > max_sparsity)
            max_sparsity = sparse_fraction;

        // generate a random sparse matrix and its fully dense equivalent
        RandomSparseMatrix(rng, A, nonzeros_per_col, height_a, height_a, width_a, width_a);
        D = MakeDense(A);
        
        // set the dimensions of B and randomly initialize
        B.Resize(height_b, height_a);
        RandomMatrix(B.Buffer(), B.LDim(), B.Height(), B.Width(), rng);

        // set the dimensions of C1 and C2 and randomly initialize
        C1.Resize(height_b, width_a);
        RandomMatrix(C1.Buffer(), C1.LDim(), C1.Height(), C1.Width(), rng);
        C2 = C1;

        // sparse Gemm
        Gemm(NORMAL, NORMAL, alpha, B, A, beta, C1);
        
        // dense Gemm
        Gemm(NORMAL, NORMAL, alpha, B, D, beta, C2);

        // check residuals
        
        // C2 = -C1 + C2
        Axpy(-1.0, C1, C2);
        double fnorm = Norm(C2, FROBENIUS_NORM);
        //cout << "[" << i << "]: fnorm of C residual: " << fnorm << endl;
        if (fnorm > max_norm)
            max_norm = fnorm;
        
        if (max_norm > GemmTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            //WriteDelimitedFile(D.LockedBuffer(), D.LDim(), D.Height(), D.Width(), "Derror.csv", 6);
            //WriteDelimitedFile(B.LockedBuffer(), B.LDim(), B.Height(), B.Width(), "Berror.csv", 6);
            break;
        }
    }

    cout << endl;
    cout << "\t******** Results for Sparse Gemm BA Test ********" << endl;
    cout << endl;
    cout << "\t\t" << GemmTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_sparsity << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_sparsity << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*************************************************" << endl;
    cout << endl;

    return (max_norm < GemmTest::MAX_ACCEPTABLE_FNORM);
}

//-----------------------------------------------------------------------------
bool TestSparseGemmBtA(Random& rng)
{
    // C = alpha*B'*A + beta*C

    SparseMatrix<double> A;
    DenseMatrix<double> D, B, C1, C2, diff;

    double alpha, beta;
    double max_norm = 0.0, min_sparsity = 1.0, max_sparsity = 0.0;

    unsigned int hw_max_threads = std::thread::hardware_concurrency();
    if (0 == hw_max_threads)
        hw_max_threads = 2;

    for (unsigned int i=0; i<GemmTest::NUM_RUNS; ++i)
    {
        unsigned int max_threads = 1 + (rng.RandomInt() % hw_max_threads);
        SetMaxThreadCount(max_threads);

        unsigned int height_a = rng.RandomRangeInt(64, 512);
        unsigned int width_a  = rng.RandomRangeInt(64, 512);
        unsigned int height_bt = rng.RandomRangeInt(16, 768);

        // generate a random alpha and beta on [-1, 1]
        alpha = rng.RandomDouble(0.0, 1.0);
        beta  = rng.RandomDouble(0.0, 1.0);

        // need to force rank2 occasionally (C.Width() == A.Width() == 2)
        if (0 == (i % 3))
            width_a = 2;

        // sparse_fraction ranges from (0.1 - 0.09, 0.1 + 0.09)
        double sparse_fraction = rng.RandomDouble(0.1, 0.09);
        unsigned int nonzeros_per_col = sparse_fraction * height_a;
        if (0 == nonzeros_per_col)
            nonzeros_per_col = 1;
        
        if (sparse_fraction < min_sparsity)
            min_sparsity = sparse_fraction;
        if (sparse_fraction > max_sparsity)
            max_sparsity = sparse_fraction;

        // generate a random sparse matrix and its fully dense equivalent
        RandomSparseMatrix(rng, A, nonzeros_per_col, height_a, height_a, width_a, width_a);
        D = MakeDense(A);
        
        // set the dimensions of B and randomly initialize
        B.Resize(height_a, height_bt);
        RandomMatrix(B.Buffer(), B.LDim(), B.Height(), B.Width(), rng);

        // set the dimensions of C1 and C2 and randomly initialize
        C1.Resize(height_bt, width_a);
        RandomMatrix(C1.Buffer(), C1.LDim(), C1.Height(), C1.Width(), rng);
        C2 = C1;

        // sparse Gemm
        Gemm(TRANSPOSE, NORMAL, alpha, B, A, beta, C1);
        
        // dense Gemm
        Gemm(TRANSPOSE, NORMAL, alpha, B, D, beta, C2);

        // check residuals
        
        // C2 = -C1 + C2
        Axpy(-1.0, C1, C2);
        double fnorm = Norm(C2, FROBENIUS_NORM);
        //cout << "[" << i << "]: fnorm of C residual: " << fnorm << endl;
        if (fnorm > max_norm)
            max_norm = fnorm;
        
        if (max_norm > GemmTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            //WriteDelimitedFile(D.LockedBuffer(), D.LDim(), D.Height(), D.Width(), "Derror.csv", 6);
            //WriteDelimitedFile(B.LockedBuffer(), B.LDim(), B.Height(), B.Width(), "Berror.csv", 6);
            break;
        }
    }

    cout << endl;
    cout << "\t******** Results for Sparse Gemm B'A Test ********" << endl;
    cout << endl;
    cout << "\t\t" << GemmTest::NUM_RUNS << " runs " << endl;
    auto prec = cout.precision();
    cout.precision(4);
    cout << "\t\tMin sparse percentage: " << 100.0*min_sparsity << endl;
    cout << "\t\tMax sparse percentage: " << 100.0*max_sparsity << endl;
    cout.precision(prec);
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t**************************************************" << endl;
    cout << endl;

    return (max_norm < GemmTest::MAX_ACCEPTABLE_FNORM);
}
