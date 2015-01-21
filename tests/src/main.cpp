#include <iostream>
#include <string>
#include "nmf.hpp"
#include "tests.hpp"
#include "utils.hpp"

using std::cout;
using std::cerr;
using std::endl;

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    bool all_ok = true;

    if (1 == argc)
    {
        cerr << "Usage: " << argv[0] << "  <path_to_data_dir>" << endl;
        return -1;
    }

    std::string data_dir = EnsureTrailingPathSep(std::string(argv[1]));

    NmfInitialize(argc, argv);

    if (!TestOpenMP())
    {
        all_ok = false;
        cerr << "OpenMP test failed" << endl;
    }

    if (!TestBpp(data_dir))
    {
        all_ok = false;
        cerr << "BPP operations test failed" << endl;
    }

    // if (!TestDenseNmf())
    // {
    //     all_ok = false;
    //     cerr << "dense nmf test failed " << endl;
    // }

    if (!TestRank2SystemSolve())
    {
        all_ok = false;
        cerr << "Rank2 system solver test failed " << endl;
    }

    if (!TestSparseGemm())
    {
        all_ok = false;
        cerr << "sparse gemm test failed " << endl;
    }

    NmfFinalize();

    if (all_ok)
        cout << "all tests passed" << endl;

    return 0;
}

    // if (!TestSparseGemmIndexed())
    // {
    //     all_ok = false;
    //     cerr << "sparse gemm indexed test failed " << endl;
    // }

    // if (!TestIndexedRank2Nmf())
    // {
    //     all_ok = false;
    //     cerr << "indexed rank2 nmf test failed" << endl;
    // }
