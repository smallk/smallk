#include <thread>
#include <vector>
#include <iostream>
#include "tests.hpp"
#include "timer.hpp"
#include "random.hpp"
#include "dense_matrix.hpp"
#include "matrix_generator.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

//-----------------------------------------------------------------------------
bool TestParallelInit()
{
    // timing experiments on parallel initialization of random W and H matrices
    
    Random rng;
    rng.SeedFromTime();

    Timer timer;

    unsigned int max_threads = std::thread::hardware_concurrency();    
    SetMaxThreadCount(max_threads);
    cout << "\n\nRandom Initialization Test: running with "
         << GetMaxThreadCount() << " threads." << endl;
    
    // W matrix is MAX_DIM x 2 and H matrix is 2 x MAX_DIM
    const unsigned int MAX_DIM = 512 * 1024 * 1024;
    
    // use a single buffer for W and H
    std::vector<double> buf(2 * MAX_DIM);

    // try various matrix sizes and compare sequential vs. parallel init time
    cout << "\n    Matrix size\t\t\tWseq\t\tWpar\t\tHseq\t\tHpar\t  Speedup" << endl;
    for (unsigned int s=1024; s<=MAX_DIM; s *= 2)
    {
        DenseMatrix<double> W(s, 2, &buf[0], s);

        // initialize W sequentially
        timer.Start();
        RandomMatrixSequential(W.Buffer(), W.LDim(), W.Height(), W.Width(), rng, 0.5, 0.5);
        timer.Stop();
        uint64_t w_seq = timer.ReportMicroseconds();
        
        DenseMatrix<double> H(2, s, &buf[0], 2);

        // initialize H sequentially
        timer.Start();
        RandomMatrixSequential(H.Buffer(), H.LDim(), H.Height(), H.Width(), rng, 0.5, 0.5);
        timer.Stop();
        uint64_t h_seq = timer.ReportMicroseconds();

        std::vector<std::thread> threads;
        
        // assumes s > num_threads
        
        // initialize W in parallel        
        timer.Start();

        for (unsigned int k=0; k<max_threads; ++k)
        {            
            threads.push_back(std::thread(ThreadFuncRank2W<double>,
                                          rng.RandomInt(),
                                          max_threads, k,
                                          W.Buffer(),
                                          W.LDim(),
                                          W.Height(),
                                          0.5, 0.5));
        }

        // wait for threads to finish
        std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
        
        timer.Stop();
        uint64_t w_par = timer.ReportMicroseconds();

        // // check W
        // for (unsigned int c=0; c<2; ++c)
        // {
        //     unsigned int offset = c*W.LDim();
        //     for (unsigned int r=0; r<W.Height(); ++r)
        //     {
        //         double val = W.Get(r, c);
        //         if (val != offset + r)
        //         {
        //             cout << "W has errors." << endl;
        //             return false;
        //         }
        //     }
        // }
        
        // initialize H in parallel
        threads.clear();
        
        timer.Start();

        for (unsigned int k=0; k<max_threads; ++k)
        {
            threads.push_back(std::thread(ThreadFuncRank2H<double>,
                                          rng.RandomInt(),
                                          max_threads, k,
                                          H.Buffer(),
                                          H.LDim(),                                          
                                          H.Width(),
                                          0.5, 0.5));
        }

        // wait for threads to finish
        std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
        
        timer.Stop();
        uint64_t h_par = timer.ReportMicroseconds();

        // // check H
        // for (unsigned int c=0; c<H.Width(); ++c)
        // {
        //     unsigned int offset = c*H.LDim();
        //     for (unsigned int r=0; r<2; ++r)
        //     {
        //         double val = H.Get(r, c);
        //         if (val != offset + r)
        //         {
        //             cout << "H has errors." << endl;
        //             return false;
        //         }
        //     }
        // }

        float w_seq_ms = w_seq * 0.001f;
        float w_par_ms = w_par * 0.001f;
        float h_seq_ms = h_seq * 0.001f;
        float h_par_ms = h_par * 0.001f;

        // save state
        std::ios::fmtflags flags(cout.flags());

        cout << std::setprecision(3);
        cout << "[" << setw(14) << (2*s) << "]\t"
             << setw(10) << std::fixed << std::right << w_seq_ms << " ms" << "\t"
             << setw(10) << std::fixed << std::right << w_par_ms << " ms" << "\t"
             << setw(10) << std::fixed << std::right << h_seq_ms << " ms" << "\t"
             << setw(10) << std::fixed << std::right << h_par_ms << " ms";
        
        // cout << "[" << setw(14) << (2*s) << "]\t"
        //      << setw(9) << (w_seq > 10000 ? w_seq/1000 : w_seq)
        //      << (w_seq > 10000 ? " ms" : " us") << "\t"
        //      << setw(9) << (w_par > 10000 ? w_par/1000 : w_par)
        //      << (w_par > 10000 ? " ms" : " us") << "\t"
        //      << setw(9) << (h_seq > 10000 ? h_seq/1000 : h_seq)
        //      << (h_seq > 10000 ? " ms" : " us") << "\t"
        //      << setw(9) << (h_par > 10000 ? h_par/1000 : h_par)
        //      << (h_par > 10000 ? " ms" : " us");

        double speedup = (w_seq + h_seq) / static_cast<double>(w_par + h_par);
        cout << "\t" << setw(9) << std::setprecision(4) << speedup << endl;

        // restore state
        cout.flags(flags);
    }
    
    return true;
}

// //-----------------------------------------------------------------------------
// void ThreadFuncInitW(const unsigned int rng_seed,
//                      const unsigned int num_threads,
//                      const unsigned int index,
//                      double* buf,
//                      const unsigned int m,
//                      const unsigned int ldim)                     
// {
//     // the W matrix is mx2; each thread initializes 1/num_threads rows in
//     // each column to avoid wrap-around issues

//     Random rng;
//     rng.SeedFromInt(rng_seed);
    
//     unsigned int num_elts = m / num_threads;
//     unsigned int start = index * num_elts;
//     unsigned int end = (index+1) * num_elts;

//     // column 0
//     for (unsigned int r=start; r<end; ++r)
//         buf[r] = rng.RandomDouble(0.5, 0.5);

//     // column 1
//     unsigned int offset = ldim;
//     for (unsigned int r=start; r<end; ++r)
//         buf[offset + r] = rng.RandomDouble(0.5, 0.5);
// }

// //-----------------------------------------------------------------------------
// void ThreadFuncInitH(const unsigned int rng_seed,
//                      const unsigned int num_threads,
//                      const unsigned int index,
//                      double* buf,
//                      const unsigned int n,
//                      const unsigned int ldim)
// {
//     // the H matrix is 2xn; each thread initializes 1/num_threads of
//     // the columns

//     Random rng;
//     rng.SeedFromInt(rng_seed);

//     unsigned int num_cols = n / num_threads;
//     unsigned int start = index * num_cols;
//     unsigned int end = (index+1) * num_cols;

//     for (unsigned int c=start; c<end; ++c)
//     {
//         unsigned int offset = c*ldim;

//         // row 0
//         buf[offset + 0] = rng.RandomDouble(0.5, 0.5);//offset + 0;

//         // row 1
//         buf[offset + 1] = rng.RandomDouble(0.5, 0.5);//offset + 1;
//     }
// }
