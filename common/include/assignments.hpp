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

#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

bool WriteAssignmentsFile(const std::vector<int>& labels, 
                          const std::string& file_path);

//-----------------------------------------------------------------------------
template <typename T>
void ComputeAssignments(std::vector<int>& assignments, 
                        const T* buf_h,            // H.Buffer()
                        const unsigned int ldim_h, // H.LDim()
                        const unsigned int k,      // H.Height()
                        const unsigned int n)      // H.Width()
{
    // This function determines cluster assignments using the H matrix.

    if (k > n)
        throw std::logic_error("ComputeAssignments: dimensions of matrix H are invalid");

    assignments.clear();
    if (assignments.size() < n)
        assignments.resize(n);
    
    // Examine each column of H and find the maximum element.  The index of
    // the row containing the max element is the cluster assignment.  
    // Elements of H are assumed to be stored in column-major order.
    
    for (unsigned int c=0; c<n; ++c)
    {
        // offset to column c in H.Buffer()
        unsigned int offset_c = c*ldim_h;

        T max_row = 0;
        T max_elt = buf_h[offset_c + 0];
        for (unsigned int r=1; r<k; ++r)
        {
            T elt = buf_h[offset_c + r];
            if (elt > max_elt)
            {
                max_elt = elt;
                max_row = r;
            }
        }

        assignments[c] = max_row;
    }
}
