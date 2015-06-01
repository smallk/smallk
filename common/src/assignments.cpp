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

#include <iomanip>
#include "assignments.hpp"

static const char COMMA = ',';

//-----------------------------------------------------------------------------
bool WriteAssignmentsFile(const std::vector<unsigned int>& labels, 
                          const std::string& file_path)
{
    // Each entry in 'labels' represents a document, and the value in the 
    // column is the cluster id to which that document was assigned.

    std::ofstream outfile(file_path);
    if (!outfile)
        return false;

    unsigned int len = labels.size();
    if (len > 0)
        outfile << labels[0];
    for (unsigned int i = 1; i < len; ++i)
        outfile << COMMA << labels[i];
    outfile << std::endl;

    outfile.close();
    return true;
}

//-----------------------------------------------------------------------------
bool WriteFuzzyAssignmentsFile(const std::vector<float>& probabilities,
                               const unsigned int k, // num clusters
                               const unsigned int n, // num docs
                               const std::string& file_path)
{
    // use three digits of precison, to limit file size
    const unsigned int P = 3;
    
    std::ofstream outfile(file_path);
    if (!outfile)
        return false;

    for (unsigned int c=0; c<n; ++c)
    {
        unsigned int offset = c*k;
        outfile << std::scientific << std::setprecision(P)
                << probabilities[offset + 0];

        for (unsigned int r=1; r<k; ++r)
            outfile << COMMA << std::scientific << std::setprecision(P)
                    << probabilities[offset + r];

        outfile << std::endl;
    }
    
    outfile.close();
    return true;
}
