#include <iostream>
#include "tests.hpp"
#include "random.hpp"
#include "dense_matrix.hpp"
#include "matrix_generator.hpp"
#include "nmf_solver_rank2.hpp"

// Test the solver for the system AX = B, where A is 2x2, X is 2xn, and
// B is 2xn. This is needed in the Rank2 NMF code.

using std::cout;
using std::cerr; 
using std::endl;

bool TestSystemSolveH(Random& rng);
bool TestSystemSolveW(Random& rng);

namespace Rank2SystemSolveTest
{
    const unsigned int NUM_RUNS = 256;
    const double MAX_ACCEPTABLE_FNORM = 1.0e-10;
}

//-----------------------------------------------------------------------------
bool TestRank2SystemSolve()
{
    Random rng;
    rng.SeedFromTime();

    bool result1 = TestSystemSolveH(rng);
    bool result2 = TestSystemSolveW(rng);

    return (result1 & result2);
}

//-----------------------------------------------------------------------------
bool TestSystemSolveH(Random& rng)
{    
    double norm, max_norm = -1.0;
    DenseMatrix<double> A(2, 2), X, B, D;
    
    for (unsigned int i=0; i<Rank2SystemSolveTest::NUM_RUNS; ++i)
    {
        // width of matrices X and B
        unsigned int n = rng.RandomRangeInt(16, 16384);

        X.Resize(2, n);
        B.Resize(2, n);
        D.Resize(2, n);

        // use random right-hand sides
        RandomMatrix(B.Buffer(), B.LDim(), B.Height(), B.Width(), rng);

        // generate a random A, but make it diagonally dominant
        RandomMatrix(A.Buffer(), A.LDim(), A.Height(), A.Width(), rng);
        A.Set(0, 0, 1.0 + A.Get(0, 0));
        A.Set(1, 1, 1.0 + A.Get(1, 1));
        //A.Set(0, 0, 2.0*A.Get(1, 0));    // force cosine formulation
        //A.Set(1, 0, 2.0*A.Get(0, 0));    // force sine formulation

        if (!SystemSolveH(X, A, B))
        {
            cerr << "TestRank2SystemSolve: SystemSolveH failure "
                 << "on iteration " << i << endl;
            continue;
        }

        // compute the residual and check the norm

        // D = A*X
        Gemm(NORMAL, NORMAL, 1.0, A, X, 0.0, D);

        // D = -1.0*B + D
        Axpy( -1.0, B, D);

        norm = Norm(D, FROBENIUS_NORM);
        if (norm > max_norm)
            max_norm = norm;

        if (max_norm > Rank2SystemSolveTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            cerr << "[" << i << "] residual fnorm: " << max_norm << endl;
            break;
        }
    }

    cout << endl;
    cout << "\t******** Results for Rank2 System Solve H Test ********" << endl;
    cout << endl;
    cout << "\t\t" << Rank2SystemSolveTest::NUM_RUNS << " runs " << endl;
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*******************************************************" << endl;
    cout << endl;

    return (max_norm < Rank2SystemSolveTest::MAX_ACCEPTABLE_FNORM);    
}

//-----------------------------------------------------------------------------
bool TestSystemSolveW(Random& rng)
{    
    double norm, max_norm = -1.0;
    DenseMatrix<double> A(2, 2), X, B, D;
    DenseMatrix<double> At(2,2), Xt, Bt, Dt;
    
    for (unsigned int i=0; i<Rank2SystemSolveTest::NUM_RUNS; ++i)
    {
        // width of matrices X and B
        unsigned int m = rng.RandomRangeInt(16, 16384);

        X.Resize(m, 2);
        B.Resize(m, 2);
        D.Resize(m, 2);

        // use random right-hand sides
        RandomMatrix(B.Buffer(), B.LDim(), B.Height(), B.Width(), rng);

        // generate a random A, but make it diagonally dominant
        RandomMatrix(A.Buffer(), A.LDim(), A.Height(), A.Width(), rng);
        A.Set(0, 0, 1.0 + A.Get(0, 0));
        A.Set(1, 1, 1.0 + A.Get(1, 1));
        //A.Set(0, 0, 2.0*A.Get(0, 1));  // force cosine formulation
        //A.Set(0, 1, 2.0*A.Get(0, 0));  // force sine formulation

        // works
        // elem::Transpose(A, At);
        // elem::Transpose(X, Xt);
        // elem::Transpose(B, Bt);
        // elem::Transpose(D, Dt);
        // SystemSolveH(Xt, At, Bt);
        // Gemm(NORMAL, NORMAL, 1.0, At, Xt, 0.0, Dt);
        // Axpy( -1.0, Bt, Dt);
        // norm = Norm(Dt, FROBENIUS_NORM);

        if (!SystemSolveW(X, A, B))
        {
            cerr << "TestRank2SystemSolve: SystemSolveW failure "
                 << "on iteration " << i << endl;
            continue;
        }

        // compute the residual and check the norm

        // D = X*A
        Gemm(NORMAL, NORMAL, 1.0, X, A, 0.0, D);

        // D = -1.0*B + D
        Axpy( -1.0, B, D);
        norm = Norm(D, FROBENIUS_NORM);

        if (norm > max_norm)
            max_norm = norm;

        if (max_norm > Rank2SystemSolveTest::MAX_ACCEPTABLE_FNORM)
        {
            cerr << "*** ERROR ***" << endl;
            cerr << "[" << i << "] residual fnorm: " << max_norm << endl;
            break;
        }
    }

    cout << endl;
    cout << "\t******** Results for Rank2 System Solve W Test ********" << endl;
    cout << endl;
    cout << "\t\t" << Rank2SystemSolveTest::NUM_RUNS << " runs " << endl;
    cout << "\t\tMax residual Frobenius norm: " << max_norm << endl;    
    cout << endl;
    cout << "\t*******************************************************" << endl;
    cout << endl;

    return (max_norm < Rank2SystemSolveTest::MAX_ACCEPTABLE_FNORM);    
}
