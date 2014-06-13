// This is an OpenMP test that can be built independently of Elemental
// and OpenBLAS.  Build the test as follows:
//
// /path/to-g++-4.8 -std=c++11 -o omp_test omp_test.cpp ../../common/src/system_{mac, posix}.cpp -I ../../common/include -I ../include -O3 -fopenmp -lgomp
//
// Run the test as follows:
//
// ./omp_test
//

#include <omp.h>
#include <cmath>
#include <thread>
#include <limits>
#include <vector>
#include <cassert>
#include <iostream>
#include "timer.hpp"
#include "tests.hpp"
#include "random.hpp"

#if defined _OPENMP
#include <omp.h>
#define OPENMP_PRAGMA(x) _Pragma(#x)
#else
#define OPENMP_PRAGMA(x)
#endif

using std::cout;
using std::cerr;
using std::endl;

bool DotProductTest(Random& rng);

namespace OpenMPTest
{
    const unsigned int NUM_RUNS = 16;
}

//-----------------------------------------------------------------------------
int main()
{
    Random rng;
    rng.SeedFromTime();

    bool result1 = DotProductTest(rng);
    return result1;
}

//-----------------------------------------------------------------------------
bool DotProductTest(Random& rng)
{
    // Fill two large arrays with random values on (-1, 1) and compute the
    // dot product.  Use OpenMP parallelization and measure the runtime.

    Timer timer;
    double s_time = 0.0, p_time = 0.0, max_residual_norm = 0.0;

    cout << "\nRunning OpenMP test..." << endl;

    // Count threads by parallel reduction, since omp_get_num_threads() 
    // does not work outside of an OpenMP parallel region.
    int num_threads = 0;
    OPENMP_PRAGMA(omp parallel reduction(+:num_threads))
    {
        num_threads += 1;
    }

    cout << "\tOpenMP is using " << num_threads
         << (1 == num_threads ? " thread." : " threads.") << endl;

    int hw_threads = std::thread::hardware_concurrency();
    cout << "\tThe HW supports " << hw_threads << " threads." << endl;

    if (num_threads < hw_threads)
    {
        omp_set_dynamic(0);
        omp_set_num_threads(hw_threads);
    }

    // count the OpenMP threads again
    num_threads = 0;
    OPENMP_PRAGMA(omp parallel reduction(+:num_threads))
    {
        num_threads += 1;
    }

    cout << "\tOpenMP is now using " << num_threads
         << (1 == num_threads ? " thread." : " threads.") << endl;
    cout << endl;

    const unsigned int N = 256*1024*1024;
    std::vector<double> A(N), B(N), C1(N), C2(N);

    cout << "\tInitializing test arrays..." << endl;

    // fill the arrays with random values
    for (unsigned int i=0; i<N; ++i)
    {
        A[i] = rng.RandomRangeDouble(-1.0, 1.0);
        B[i] = rng.RandomRangeDouble(-1.0, 1.0);
    }

    for (unsigned int i=0; i<OpenMPTest::NUM_RUNS; ++i)
    {
        timer.Start();
        for (unsigned int r=0; r<N; ++r)
        {
            C1[i] = A[i] * B[i];
        }

        timer.Stop();
        s_time += timer.ReportMilliseconds();
    
        timer.Start();
        OPENMP_PRAGMA(omp parallel)
        {
            OPENMP_PRAGMA(omp for nowait)
            for (unsigned int r=0; r<N; ++r)
            {
                C2[i] = A[i] * B[i];
            }
        }
        timer.Stop();
        p_time += timer.ReportMilliseconds();

        double residual_norm = 0.0;
        for (unsigned int r=0; r<N; ++r)
        {
            double val = C2[i] - C1[i];
            residual_norm += val*val;
        }
        residual_norm = sqrt(residual_norm);

        if (residual_norm > max_residual_norm)
            max_residual_norm = residual_norm;

        cout << "\tCompleted run " << i+1 << "/" << OpenMPTest::NUM_RUNS 
             << "." << endl;
    }

    double speedup = s_time / p_time;

    cout << endl;
    cout << "\t*********** Results for OpenMP Test ************" << endl;
    cout << endl;
    cout << "\t\t" << OpenMPTest::NUM_RUNS << " runs " << endl;
    cout << endl;
    cout << "\t\tSequential runtime: " << s_time << " ms." << endl;    
    cout << "\t\t  Parallel runtime: " << p_time << " ms." << endl;
    cout << endl;
    cout << "\t\t           Speedup: " << speedup << "." << endl;
    cout << endl;
    cout << "\t\tMax residual norm: " << max_residual_norm << "." << endl;
    cout << endl;
    cout << "\t*************************************************" << endl;
    cout << endl;

    return (speedup > 1.0);
}
