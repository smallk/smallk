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
#include <cassert>
#include <iostream>
#include <stdexcept>
#include "smallk.hpp"

using std::cout;
using std::cerr;
using std::endl;

const std::string FILENAME_W("nmf_init_w.csv");
const std::string FILENAME_H("nmf_init_h.csv");
const std::string FILENAME_MATRIX("reuters.mtx");
const std::string FILENAME_DICT("reuters_dictionary.txt");

static const char PATH_SEP = '/';

//-----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (1 == argc)
    {
        cerr << "usage: " << argv[0] << " <path_to_data_dir> " << endl;
        return -1;
    }

    // ensure that data_dir  ends in a trailing '/' character
    std::string data_dir(argv[1]);
    if (data_dir.size() >= 1)
    {
        auto sz = data_dir.size();
        if (PATH_SEP != data_dir[sz-1])
            data_dir += PATH_SEP;
    }
  
    std::string filepath_matrix = data_dir + FILENAME_MATRIX;
    std::string filepath_w      = data_dir + FILENAME_W;
    std::string filepath_h      = data_dir + FILENAME_H;
    std::string filepath_dict   = data_dir + FILENAME_DICT;

    smallk::Initialize(argc, argv);
    assert(smallk::IsInitialized());

    cout << "Smallk major version: " << smallk::GetMajorVersion() << endl;
    cout << "Smallk minor version: " << smallk::GetMinorVersion() << endl;
    cout << "Smallk patch level:   " << smallk::GetPatchLevel() << endl;
    cout << "Smallk version string: " << smallk::GetVersionString() << endl;

    try
    {
        ///////////////////////////////////////////////////////////////////////
        //
        //                set params and verify values
        //
        ///////////////////////////////////////////////////////////////////////

        smallk::SetOutputPrecision(8);
        if (8 != smallk::GetOutputPrecision())
            throw std::runtime_error("SetOutputPrecision failed");

        smallk::SetNmfTolerance(1.0e-5);
        if (1.0e-5 != smallk::GetNmfTolerance())
            throw std::runtime_error("SetNmfTolerance failed");

        smallk::SetMaxIter(2500);
        if (2500 != smallk::GetMaxIter())
            throw std::runtime_error("SetMaxIter failed");

        smallk::SetMinIter(10);
        if (10 != smallk::GetMinIter())
            throw std::runtime_error("GetMaxIter failed");
        
        smallk::SetMaxThreads(3);
        if (3 != smallk::GetMaxThreads())
            throw std::runtime_error("SetMaxThreads failed");

        smallk::SetMaxTerms(6);
        if (6 != smallk::GetMaxTerms())
            throw std::runtime_error("SetMaxTerms failed");

        smallk::SetOutputFormat(smallk::OutputFormat::XML);
        if (smallk::OutputFormat::XML != smallk::GetOutputFormat())
            throw std::runtime_error("SetOutputFormat failed");

        smallk::SetHierNmf2Tolerance(0.01);
        if (0.01 != smallk::GetHierNmf2Tolerance())
            throw std::runtime_error("SetHierNmf2Tolerance");

        ///////////////////////////////////////////////////////////////////////
        //
        // Reset params to defaults and load the Reuters matrix.
        //
        ///////////////////////////////////////////////////////////////////////

        smallk::Reset();

        smallk::LoadMatrix(filepath_matrix);
        assert(smallk::IsMatrixLoaded());
        
        ///////////////////////////////////////////////////////////////////////
        //
        // Factor the Reuters matrix using BPP and k == 8. Use 
        // initializer matrices for W and H.
        //
        ///////////////////////////////////////////////////////////////////////
        
        cout << "\nRunning NMF-BPP...\n" << endl;

        smallk::SetMinIter(1);
        smallk::SetOutputPrecision(6);
        smallk::Nmf(8, smallk::Algorithm::BPP, filepath_w, filepath_h);

        ///////////////////////////////////////////////////////////////////////
        //
        // Run a hierarchical clustering problem and generate 5 clusters.
        // Use XML format for the clustering result file.
        //
        ///////////////////////////////////////////////////////////////////////

        cout << "\nRunning HierNmf2...\n" << endl;
        
        smallk::SetOutputFormat(smallk::OutputFormat::XML);
        smallk::LoadDictionary(filepath_dict);
        smallk::HierNmf2(5);
    }
    catch (std::exception& e)
    {
        cerr << e.what() << endl;
    }
    
    smallk::Finalize();
    return 0;
}

//-----------------------------------------------------------------------------
void MinAndMax(double& minval, 
               double& maxval,
               const double* buf, 
               const unsigned int ldim, 
               const unsigned int height, 
               const unsigned int width)
{
    maxval = 0.0;
    minval = std::numeric_limits<double>::max();

    // the array is stored in column-major order
    for (unsigned int c=0; c<width; ++c)
    {
        unsigned int col_offset = c*ldim;
        for (unsigned int r=0; r<height; ++r)
        {
            double val = buf[col_offset + r];
            if (val < minval)
                minval = val;
            if (val > maxval)
                maxval = val;
        }
    }
}
