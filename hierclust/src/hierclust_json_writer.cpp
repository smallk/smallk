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

#include <sstream>
#include <iostream>
#include "hierclust_json_writer.hpp"

using std::cout;
using std::endl;

// indentation
const std::string S4("    ");
const std::string S8  = S4 + S4;
const std::string S12 = S4 + S8;
const std::string S16 = S4 + S12;

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteHeader(std::ofstream& outfile, 
                                     const int doc_count)
{
    std::ostringstream count_str;
    count_str << doc_count;

    outfile << "{" << endl;
    outfile << S4 << "\"doc_count\": " << count_str.str() << "," << endl;
    outfile << S4 << "\"nodes\": [" << endl;
    nodes_written_ = 0u;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteNodeBegin(std::ofstream& outfile, 
                                        const int node_id)
{
    if (nodes_written_ > 0u)
    {
        // write the trailing comma for the previous element
        outfile << "," << endl;
    }
    outfile << S8 << "{" << endl;
    outfile << S12 << "\"id\": " << node_id << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteLevel(std::ofstream& outfile, 
                                    const int level)
{
    outfile << S12 << "\"level\": " << level << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteLabel(std::ofstream& outfile, 
                                    const int label)
{
    outfile << S12 << "\"label\": " << label << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteParentId(std::ofstream& outfile, 
                                       const int parent_id)
{
    outfile << S12 << "\"parent_id\": " << parent_id << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteParentLabel(std::ofstream& outfile, 
                                          const int parent_label)
{
    outfile << S12 << "\"parent_label\": " << parent_label << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteLeftChild(std::ofstream& outfile, 
                                        const bool is_left_child,
                                        const int lc_label)
{
    outfile << S12 << "\"left_child\": " 
            << (is_left_child ? "true" : "false") << "," << endl;
    
    outfile << S12 << "\"left_child_label\": " << lc_label << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteRightChild(std::ofstream& outfile, 
                                         const int rc_label)
{
    outfile << S12 << "\"right_child_label\": " << rc_label << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteDocCount(std::ofstream& outfile, 
                                       const int count)
{
    outfile << S12 << "\"doc_count\": " << count << "," << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteTopTerms(std::ofstream& outfile, 
                                       const std::vector<int>& term_indices,
                                       const std::vector<std::string>& dictionary)
{
    unsigned int num_terms = term_indices.size();
    if (0u == num_terms)
        return;

    outfile << S12 << "\"top_terms\": [" << endl;
    if (num_terms > 1)
    {
        for (size_t q=0; q != (num_terms-1); ++q)
        {
            int index = term_indices[q];
            outfile << S16 << "\"" << dictionary[index] << "\"," << endl;
        }
    }

    // final term
    int index = term_indices[num_terms - 1];
    outfile << S16 << "\"" << dictionary[index] << "\"" << endl;
    outfile << S12 << "]" << endl;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteNodeEnd(std::ofstream& outfile)
{
    outfile << S8 << "}";
    ++nodes_written_;
}

//-----------------------------------------------------------------------------
void HierclustJsonWriter::WriteFooter(std::ofstream& outfile)
{
    outfile << endl;
    outfile << S4 << "]" << endl;
    outfile << "}" << endl;
}

